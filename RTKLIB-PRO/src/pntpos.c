/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*           2020/11/30 1.7  support NavIC/IRNSS in pntpos()
*                           no support IONOOPT_LEX option in ioncorr()
*                           improve handling of TGD correction for each system
*                           use E1-E5b for Galileo dual-freq iono-correction
*                           use API sat2freq() to get carrier frequency
*                           add output of velocity estimation error in estvel()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))
#define MAX(x,y)    ((x)>=(y)?(x):(y))

#define QZSDT /* enable GPS-QZS time offset estimation */
#ifdef QZSDT
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
# define MAX_GDOP   30          /* max gdop for valid solution  */

/* 计算单个卫星伪距观测值误差方差，用以随机模型定权------------------------------------
 args   :const prcopt_t *opt   I   处理过程选项
         const ssat_t   *ssat  I   卫星状态
         const obsd_t   *obs   I   观测量数据
               double   el     I   卫星高度角 (rad)
               int      sys    I   卫星系统 (SYS_???)
               int      f2     I   IFLC第二频率索引（-1表示非IFLC）
 return  : 测量误差方差 (m^2)
------------------------------------*/
static double varerr(const prcopt_t *opt, const ssat_t *ssat, const obsd_t *obs, double el, int sys, int f2)
{
    double fact=1.0,varr,snr_rover;
    double iflc_factor=1.0; /* IFLC模式方差放大因子 */

    switch (sys) { // 确定各导航系统误差放大因子
        case SYS_GPS: fact *= EFACT_GPS; break;
        case SYS_GLO: fact *= EFACT_GLO; break;
        case SYS_SBS: fact *= EFACT_SBS; break;
        case SYS_CMP: fact *= EFACT_CMP; break;
        case SYS_QZS: fact *= EFACT_QZS; break;
        case SYS_IRN: fact *= EFACT_IRN; break;
        default:      fact *= EFACT_GPS; break;
    }
    if (el<MIN_EL) el=MIN_EL;
    /* var = R^2*(a^2 + (b^2/sin(el) + c^2*(10^(0.1*(snr_max-snr_rover)))) + (d*rcv_std)^2) */
    varr=SQR(opt->err[1])+SQR(opt->err[2])/sin(el);
    if (opt->err[6]>0.0) {  /* if snr term not zero */
        snr_rover=(ssat)?ssat->snr_rover[0]:opt->err[5]; // 默认使用第一频点的观测SNR
        varr+=SQR(opt->err[6])*pow(10,0.1*MAX(opt->err[5]-snr_rover,0));
    }
    varr*=SQR(opt->eratio[0]);
    if (opt->err[7]>0.0) {
        varr+=SQR(opt->err[7]*obs->Pstd[0]);
    }
    /* 根据IFLC模式和退化情况设置不同方差放大因子 */
    if (f2>=0) {
        int has_f2 = (obs->P[f2]!=0.0); // 第二频率是否存在
        if (has_f2) {
            /* 标准IFLC：L1-L5组合方差放大约5倍，L1-L2组合方差放大约9倍 */
            iflc_factor = (f2==2) ? SQR(2.2) : SQR(3.0);
        } 
    }
    varr*=iflc_factor;
    return SQR(fact)*varr;
}
/* 计算指定卫星的群延迟参数TGD (s→m) ---------------------------------------------
args   :       int        sat     I   卫星编号
         const nav_t     *nav     I   导航数据
		       int        type    I   TGD类型 (0:TGD_L1L2, 1:TGD_L2L5, 2:TGD_B1I_B1Cp, 3:TGD_B2I_B2bI, 4:ISC_B1Cd)
（需要调试！是否会改正L1与L5之间TGD参数（目前似乎看只改正L1-L2））
return : TGD m
---------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    
	if (sys == SYS_GLO) { // GLONASS卫星单独计算
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* test SNR mask 信噪比阈值测试 -------------------------------------------------------------
args:    const obsd_t   *obs     I   观测量数据
         const double   *azel    I   卫星方位角和高度角 {az,el} (rad)
		 const prcopt_t *opt     I   处理过程选项
         int            f2       I   IFLC第二频率索引（-1表示非IFLC）

return: 1:ok, 0: 拒绝卫星
------------------------------------------------------------- */
static int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt, int f2)
{
    if (testsnr(0, 0, azel[1], obs->SNR[0], &opt->snrmask)) { //第一频点SNR门限检查，若不满足则直接拒绝卫星
        return 0;
    }
    if (opt->ionoopt==IONOOPT_IFLC && f2>=0) {
        /* 只在第二频点观测值存在时才检查其SNR门限，避免L5缺失导致卫星被误拒 */
        if (obs->SNR[f2]>0.0 &&
            testsnr(0,f2,azel[1],obs->SNR[f2],&opt->snrmask)) return 0;
    }
    return 1;
}
/* 无电离层或“伪无电离”，伪距差分码偏差校正DCB、IFLC、TGD----------------------------------------
args:  const obsd_t      *obs     I   观测量数据
       const nav_t       *nav     I   导航数据
       const prcopt_t    *opt     I   处理过程选项
             int          f2      I   IFLC第二频率索引（-1表示非IFLC）
	         double      *var     O   伪距差分码偏差DCB校正的方差 (m^2)

return: 伪距差分码偏差DCB校正后的伪距值 (m)
------------------------------------------------------------- */
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt, int f2,
                     double *var)
{
    double P1,P2=0.0,gamma,b1,b2;
	int sat,sys,bias_ix,use_iflc;

	sat=obs->sat; // 卫星编号
	sys=satsys(sat,NULL); // 卫星系统
	P1=obs->P[0];
	if (f2>=0) P2=obs->P[f2];
	*var=0.0;

	if (P1==0.0 && (f2<0 || P2==0.0)) return 0.0; /* L1和L5均不存在则跳过 */

	use_iflc = (opt->ionoopt==IONOOPT_IFLC && f2>=0 && P1!=0.0 && P2!=0.0); // 确认是否真的使用IFLC

	/* DCB校正：L1 DCB校正 */
	if (P1!=0.0) {
		bias_ix=code2bias_ix(sys,obs->code[0]);  /* L1 DCB校正 */
		if (bias_ix>0) { /* 若为0代表为基准，无需校正 */
			P1+=nav->cbias[sat-1][0][bias_ix-1];
		}
	}

    /* DCB校正：L2/L5 DCB校正 */
    if (f2>=0&&P2!=0.0) {
        bias_ix=code2bias_ix(sys,obs->code[f2]);
        if (bias_ix>0) { 
            P2+=nav->cbias[sat-1][f2][bias_ix-1];
        }
    }

    /* 构建IFLC，TGD校正 */
    if (use_iflc) {
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 or L1-L5 */
            gamma=f2==1?SQR(FREQL1/FREQL2):SQR(FREQL1/FREQL5);
            if (f2!=1) { /* L1-L5组合偏离了星历钟差的L1-L2基准，需要修正 */
                b1=gettgd(sat,nav,0);
                b2=gettgd(sat,nav,1);
                return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
            }
            /* GPS/QZS L1-L2无需TGD校正（星历钟差基准一致） */
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 or G1-G3 */
            gamma=f2==1?SQR(FREQ1_GLO/FREQ2_GLO):SQR(FREQ1_GLO/FREQ3_GLO);
            b1=0.0; /* GLONASS星历钟差基准即为G1，因此G1无需群延迟修正 */
            b2=(f2==1)?gettgd(sat,nav,0):0.0;  /* 对G2修正，G3通常没有TGD参数，因此不校正 */
            return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b or E1-E5a */
            gamma=f2==1?SQR(FREQL1/FREQE5b):SQR(FREQL1/FREQL5);
            if (f2==1&&getseleph(SYS_GAL)) { /* F/NAV 类型导航文件以E5a为基准 */
                P2-=gettgd(sat,nav,0)-gettgd(sat,nav,1); /* 当L1-L2组合，TGD:E5a→E5b */
            }
            /* Galileo E1-E5a组合无需TGD校正（两者均以E1为基准） */
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_CMP) { /* 应该有bug：B1-B2/ B1-B2a */
			gamma = SQR(((obs->code[0] == CODE_L2I) ? FREQ1_CMP : FREQL1) / FREQ2_CMP); // B1I使用FREQ1_CMP，B1C使用FREQL1; 第二频率默认使用FREQ2_CMP（B2I/B2b）
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            b2=gettgd(sat,nav,1); /* TGD_B2I/B2bI (m) */
            return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
        }
        else if (sys==SYS_IRN) { /* L5-S */
            gamma=SQR(FREQL5/FREQs);
            return (P2-gamma*P1)/(1.0-gamma);
        }
    }

    /* 单频处理：IFLC模式缺失L5时退化，或仅有L5/L1频段 */
    else {
        *var = SQR(ERR_CBIAS); // 单频的未建模误差(Code Bias等)，因此方差设置为一个较大的值以降低其权重

        /* 标准单频TGD校正（L5） */
        if (P2!=0.0 && P1==0.0) {
            /* 仅L5存在时：退化单频L5消电离层处理（智能手机仅支持L5频段情况） */
            if (sys==SYS_GPS||sys==SYS_QZS) {
                /* GPS/QZS 星历钟差基于 L1/L2 无电离层组合。
                   L5 需同时修正 TGD (L1-L2 IF基准) 和 ISC_L5 (L1-L5偏差) */
                double tgd = gettgd(sat, nav, 0); /* TGD:L5→L2 */
                double isc = gettgd(sat, nav, 1); /* TGD:L2→L1 */
                return P2 - tgd - isc;
            }
            else if (sys == SYS_GAL) {
                /* Galileo 星历钟差基于 E1-E5a/E5b 无电离层组合，单频 E5a/E5b 需乘 gamma 补偿 BGD */
                if (getseleph(SYS_GAL)) { /* F/NAV 以 E1-E5a 为基准 */
                    gamma = SQR(FREQL1 / FREQL5);
                    b1 = gettgd(sat, nav, 0); /* BGD(E1,E5a) */
                }
                else {                  /* I/NAV 以 E1-E5b 为基准 */
                    gamma = SQR(FREQL1 / FREQE5b);
                    b1 = gettgd(sat, nav, 1); /* BGD(E1,E5b) */
                }
                return P2 - gamma * b1;
            }
            else if (sys==SYS_CMP) { /* B2a */
                b1=gettgd(sat,nav,3); /* TGD_B2I_B2bI */
                return P2-b1;
            }
        }

        /* 标准单频TGD校正（L1） */
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys == SYS_GLO) { /* G1 */
            return P1; /* GLONASS G1: 星历钟差已经对准G1，无需群延迟改正 */
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* F/NAV基于E5a */
            else                    b1=gettgd(sat,nav,1); /* I/NAV基于E5b */
            return P1-b1;
        }
        else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            return P1-b1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
    }
    return P1;
}
/* 电离层修正 ------------------------------------------------------
* args   : gtime_t time     I   时间
*          nav_t  *nav      I   导航数据
*          int    sat       I   卫星编号
*          double *pos      I   接收机位置 {lat,lon,h} (rad|m)
*          double *azel     I   方位角/高度角 {az,el} (rad)
*          int    ionoopt   I   电离层校正选项 (IONOOPT_???)
*          double *ion      O   电离层延迟 (L1) (m)
*          double *var      O   电离层延迟 (L1) 方差 (m^2)
* return : 状态 (1:成功, 0:错误)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    int err=0;

    char tstr[40];
    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time2str(time,tstr,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* SBAS ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {    
        if (sbsioncorr(time,nav,pos,azel,ion,var)) return 1;
        err=1;
    }
    /* IONEX TEC model */
    if (ionoopt==IONOOPT_TEC) {
        if (iontec(time,nav,pos,azel,1,ion,var)) return 1;
        err=1;
    }
    /* QZSS broadcast ionosphere model（Klobuchar） */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* GPS broadcast ionosphere model（Klobuchar） */
    if (ionoopt==IONOOPT_BRDC||err==1) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* no correction */
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* 对流层校正 -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   时间
*          nav_t  *nav      I   卫星数据
*          double *pos      I   接收机位置 {lat,lon,h} (rad|m)
*          double *azel     I   方位角/高度角 {az,el} (rad)
*          int    tropopt   I   对流层校正项 (TROPOPT_???)
*          double *trp      O   对流层延迟 (m)
*          double *var      O   对流层延迟方差 (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    char tstr[40];
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time2str(time,tstr,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* Saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* SBAS (MOPS) troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* 残差计算、设计矩阵构建 -----------------------------------------------------
    int      iter      I   迭代次数，在estpos()里迭代调用，第i次迭代就传i
    obsd_t   *obs      I   观测量数据
    int      n         I   观测量数据的数量
    double   *rs       I   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
    double   *dts      I   卫星钟差，长度为2*n， {bias,drift} (s|s/s)
    double   *vare     I   卫星位置和钟差的协方差 (m^2)
    int      *svh      I   卫星健康标志 (-1:correction not available)
    nav_t    *nav      I   导航数据
    double   *x        I   本次迭代开始之前的定位值,7*1,前3个是本次迭代开始之前的定位值，第4个是钟差，后三个分别是gps系统与glonass、galileo、bds系统的钟差。
    prcopt_t *opt      I   处理过程选项
    ssat_t   *ssat     I   卫星状态
    double   *v        O   定位方程的右端部分，伪距残差
    double   *H        O   定位方程中的几何矩阵
    double   *var      O   参与定位的伪距残差的方差
    double   *azel     O   对于当前定位值，所有观测卫星的 {方位角、高度角} (2*n)
    int      *vsat     O   所有观测卫星在当前定位时是否有效 (1*n)
    double   *resp     O   所有观测卫星的伪距残差，(P-(r+c*dtr-c*dts+I+T)) (1*n)
    int      *ns       O   参与定位的卫星的个数
return : nv
----------------------------------------------------- */
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   const ssat_t *ssat, double *v, double *H, double *var,
                   double *azel, int *vsat, double *resp, int *ns)
{
    gtime_t time;
    double r,freq,dion=0.0,dtrp=0.0,vmeas,vion=0.0,vtrp=0.0,rr[3],pos[3],dtr,e[3],P;
	int i,j,nv=0,sat,sys,mask[NX-3]={0},f2; // nv表示有效观测数（计数器）；mask数组用于标记不同导航系统的时间偏差是否被处理过

    //将之前得到的定位解信息赋值给rr和dtr数组
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];

    ecef2pos(rr,pos); // rr{x,y,z}->pos{lat,lon,h}
    trace(3,"rescode: rr=%.3f %.3f %.3f\n",rr[0], rr[1], rr[2]);

    //遍历当前历元所有OBS[]
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0; // 初始化，因为在前后两次定位结果中，每颗卫星的上述信息都会发生变化。
        time=obs[i].time; // 时间
        sat=obs[i].sat; // 卫星编号
        if (!(sys=satsys(sat,NULL))) continue; //1.调用satsys()函数，验证卫星编号是否合理，查找其所属的导航系统
        f2=(opt->ionoopt==IONOOPT_IFLC)?seliflc(opt->nf,sys):-1; // 确定IFLC第二频率索引，在循环开始处确定以供全程使用
        
        // 剔除重复的观测数据
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            char tstr[40];
            trace(2,"duplicated obs data %s sat=%d\n",time2str(time,tstr,3),sat);
            i++;
            continue;
        }
        /* 2.排除事先指定的卫星 */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* 3-4.卫星与接收机之间几何距离和高度角mask */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
        
        if (iter>0) { // 第一次迭代时候缺少定位数据，强行经过下面修正会产生巨大误差
            /* 5.test SNR mask */
            if (!snrmask(obs+i,azel+i*2,opt,f2)) continue;

			/* 6.电离层校正 */
			/* 若为IFLC模式但缺失任一频率，退化为BRDC模型进行单频补偿 */
			int use_iflc = (opt->ionoopt==IONOOPT_IFLC && f2>=0 && obs[i].P[0]!=0.0 && obs[i].P[f2]!=0.0);
			int ionoopt_cur = (opt->ionoopt==IONOOPT_IFLC && !use_iflc) ? IONOOPT_BRDC : opt->ionoopt;

			if (ionoopt_cur != IONOOPT_IFLC) { // 单频电离层校正，或IFLC模式退化为单频处理
				if (!ionocorr(time,nav,sat,pos,azel+i*2,ionoopt_cur,&dion,&vion)) {
					continue;
				}
				/* 单频以实际使用的频率为基准 */
				uint8_t curr_code = obs[i].code[0];
				if (obs[i].P[0] == 0.0 && f2>=0 && obs[i].P[f2] != 0.0) {
					curr_code = obs[i].code[f2];
				}
				if ((freq=sat2freq(sat,curr_code,nav))==0.0) continue;
				dion*=SQR(FREQL1/freq);
				vion*=SQR(SQR(FREQL1/freq));
			} else {
				dion=0.0; vion=0.0; // IFLC
			}

            /* 7.对流层校正 */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
            }
        }
        /* 8.计算校正伪距（DCB、IFLC、TGD） */
        if ((P=prange(obs+i,nav,opt,f2,&vmeas))==0.0) continue;

        /* 9.计算此时伪距残差，累加测距误差 */
        // 计算伪距残差(P-(r+c*dtr-c*dts+I+T)),程序中dtr单位为m
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        trace(4,"sat=%d: v=%.3f P=%.3f r=%.3f dtr=%.6f dts=%.6f dion=%.3f dtrp=%.3f\n",
            sat,v[nv],P,r,dtr,dts[i*2],dion,dtrp);
       
        // 10.设计矩阵
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0); // 前3列是坐标改正数（赋卫星与接收机之间的单位矢量），第四列是GPS钟差（赋1），后面是不同导航系统之间的时间偏差（赋初值0）
        }
        // 处理不同导航系统之间的时间偏差（伪距残差减去与基准GPS的接收机钟差），修改矩阵 H
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
#ifdef QZSDT
        else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
#endif
        else mask[0]=1; // GPS

        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        // 11.计算用户测距总方差，用以随机模型定权
        var[nv]=vare[i]+vmeas+vion+vtrp;
        if (ssat)
            var[nv++]+=varerr(opt,&ssat[i],&obs[i],azel[1+i*2],sys,f2);
        else
            var[nv++]+=varerr(opt,NULL,&obs[i],azel[1+i*2],sys,f2);
        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }
    // 为防止不满秩，检查是否有缺失的导航系统，若有进行下述处理，以防止矩阵H秩亏
    for (i=0;i<NX-3;i++) {
		if (mask[i]) continue; // 表示该导航系统的时间偏差已经被处理过了（有该系统），跳过
        v[nv]=0.0; // 否则就（新增一个观测方程）伪距残差设置为0，设计矩阵设置为1，方差设置为0.01以补满秩
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0; // 新增的方程代表的意思就是接收机的这个导航系统钟差为0
        var[nv++]=0.01;
    }
	return nv; // 返回有效观测数
}
/* 对定位结果进行卡方检验和GDOP检验---------------------------------------------------------
    const double   *azel     方位角、高度角
    const int      *vsat     观测卫星在当前定位时是否有效 (1*n)
          int      n         观测值个数
    const prcopt_t *opt      处理选项
    const double   *v        定位方程的右端部分，伪距残差
          int      nv        观测值数
          int      nx        待估计参数数
          char     *msg      错误消息
---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    /* 对残差卡方检验 */
    vv=dot(v,v,nv);  // chisqr:卡方值表
    if (nv>nx&&vv>chisqr[nv-nx-1]) {  //观测值数大于待估计参数数  nv-nx-1:多余观测数
        sprintf(msg,"Warning: large chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        /* return 0; */ /* 阈值对所有用例都过于严格，报告错误但继续前进 */
    }
    /* GDOP检验 */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>MAX_GDOP) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}
/* 用伪距估算接收机位置及钟差 ------------------------------------------------
    obsd_t   *obs      I   观测量数据
    int      n         I   观测量数据的数量
    double   *rs       I   卫星位置和速度，长度为6*n，{x,y,z,vx,vy,vz}(ecef)(m,m/s)
    double   *dts      I   卫星钟差，长度为2*n，{dt,dt_sap}(s)
    double   *vare     I   卫星位置和钟差的协方差 (m^2)
    int      *svh      I   卫星健康标志 (-1:correction not available)
    nav_t    *nav      I   导航数据
    prcopt_t *opt      I   处理过程选项
    ssat_t   *ssat     I   卫星状态
    sol_t    *sol      IO  结果
    double   *azel     IO  方位角和俯仰角 (rad)
    int      *vsat     IO  卫星在定位时是否有效
    double   *resp     IO  定位后伪距残差 (P-(r+c*dtr-c*dts+I+T))
    char     *msg      O   错误消息
retrun : 0:解无效，1:单点定位解，2:SBAS单点定位解
------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, const ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+NX-3,1); H=mat(NX,n+NX-3); var=mat(n+NX-3,1); // 分配内存，长度为n+NX-3，是因为在rescode函数中，如果参与定位的卫星系统个数小于待估计系统钟差参数个数
                                                          // （各个卫星系统的接收机钟差都是待估参数，但不一定都有数据），则需要补满秩，用于补满秩的方程个数为NX-3
    
    for (i=0;i<3;i++) x[i]=sol->rr[i]; // 初始化接收机位置（将上一历元的位置作为初值，若初次则赋为0）

    // 开始迭代定位计算
    for (i=0;i<MAXITR;i++) {

        /*  1.调用 rescode 函数，计算当前迭代的伪距残差 v、几何矩阵 H、
            伪距残差的方差 var、所有观测卫星的方位角和仰角 azel、定位时有效性 vsat、
            定位后伪距残差 resp、参与定位的卫星个数 ns 和方程个数 nv */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,vsat,resp,
                   &ns);
        
        if (nv<NX) { // 确定方程组中方程的个数要大于未知数的个数
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        /* 以伪距残差的标准差的倒数作为权重，对H和v分别左乘权重对角阵，得到加权之后的H和v */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* 2.调用lsq(最小二乘估计)函数,得到当前x的修改量dx和定位误差协方差矩阵中的权系数阵Q */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) { // 更新估计参数
            x[j]+=dx[j];
        }
        // 如果求得的待估参数变化量小于截断因子(目前是1E-4)，则将x[j]作为最终的定位结果，
        // 对 sol 的相应参数赋值,之后再调用 valsol 函数确认当前解是否符合要求,参考 RTKLIB Manual P162
        // 否则，进行下一次循环。
        if (norm(dx,NX)<1E-4) {
            sol->type=0;
			// dtr：接收机钟差(秒)，x[i]单位是m，因此需要除以光速CLIGHT将其转换为秒
			sol->time = timeadd(obs[0].time, -x[3] / CLIGHT); // 接收机时间 = 观测时间 - 接收机钟差
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
#ifdef QZSDT
            sol->dtr[5]=x[8]/CLIGHT; /* QZS-GPS time offset (s) */
#endif
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;
            
            /* 3.对定位结果进行卡方检验和GDOP检验 */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    //如果超过了规定的循环次数，则输出发散信息后，return 0
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
    free(v); free(H); free(var);
    return 0;
}
/* RAIM FDE (failure detection and exclusion) -------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, const ssat_t *ssat, 
                    sol_t *sol, double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[40],name[8],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time2str(obs[0].time,tstr,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclusion */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,ssat,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        sol_e.eventime = sol->eventime;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
#ifdef TRACE
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
#endif
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* range rate residuals ------------------------------------------------------*/
static int resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double err, double *v,
                  double *H)
{
    double freq,rate,pos[3],E[9],a[3],e[3],vs[3],cosel,sig;
    int i,j,nv=0;
    
    trace(3,"resdop  : n=%d\n",n);
    
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        freq=sat2freq(obs[i].sat,obs[i].code[0],nav);
        
        if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,E,a,e);
        
        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }
        /* range rate with earth rotation correction */
        rate=dot3(vs,e)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                     rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* Std of range rate error (m/s) */
        sig=(err<=0.0)?1.0:err*CLIGHT/freq;
        
        /* range rate residual (m/s) */
        v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
        
        /* design matrix */
        for (j=0;j<4;j++) {
            H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
        }
        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
    double err=opt->err[4]; /* Doppler error (Hz) */
    int i,j,nv;
    
    v=mat(n,1); H=mat(4,n);
    
    for (i=0;i<MAXITR;i++) {
        
        /* range rate residuals (m/s) */
        if ((nv=resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H))<4) {
            break;
        }
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            trace(3,"estvel : vx=%.3f vy=%.3f vz=%.3f, n=%d\n",x[0],x[1],x[2],n);
            matcpy(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Q[0];  /* xx */
            sol->qv[1]=(float)Q[5];  /* yy */
            sol->qv[2]=(float)Q[10]; /* zz */
            sol->qv[3]=(float)Q[1];  /* xy */
            sol->qv[4]=(float)Q[6];  /* yz */
            sol->qv[5]=(float)Q[2];  /* zx */
            break;
        }
    }
    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data OBS观测数据
*          int    n         I   number of observation data OBS数量
*          nav_t  *nav      I   navigation data NAV导航电文数据
*          prcopt_t *opt    I   processing options 处理过程选项
*          sol_t  *sol      IO  solution 结果
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output) 方位角和高度角
*          ssat_t *ssat     IO  satellite status              (NULL: no output) 卫星状态
*          char   *msg      O   error message for error exit 错误信息
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    
    char tstr[40];
    trace(3,"pntpos  : tobs=%s n=%d\n",time2str(obs[0].time,tstr,3),n);
    
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {     //检验观测值数是否大于0
        strcpy(msg,"no observation data");
        return 0;
    }
    sol->time=obs[0].time; //sol->time赋值第一个观测值的时间
    msg[0]='\0';
    sol->eventime = obs[0].eventime;
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].snr_rover[0]=0;
            ssat[i].snr_base[0]=0;
        }
        for (i=0;i<n;i++)
            ssat[obs[i].sat-1].snr_rover[0]=obs[i].SNR[0];
    }
    
    if (opt_.mode!=PMODE_SINGLE) { //如果处理选项不是SPP
        opt_.ionoopt=IONOOPT_BRDC; //电离层矫正选Klobuchar广播星历模型
        opt_.tropopt=TROPOPT_SAAS; //对流层矫正采用Saastmoinen模型
    }
    /* 1.计算卫星位置、速度和钟差 */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* 2.用伪距估算接收机位置及钟差 */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    
    /* 3.接收机自主完好性监测与粗差剔除 */
    if (!stat&&n>=6&&opt->posopt[4]) {      //estpos中valsol检验失败，即位置估计失败，
                                            //会调用RAIM接收机自主完好性监测重新估计，
                                            //前提是卫星数>6、对应参数解算设置opt->posopt[4]=1（似乎现在配置文件没）
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }
    /* 4.用多普勒估算接收机速度 */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];   //存入方位角和高度角
    }
    if (ssat) {     //赋值卫星状态结构体ssat
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
        }
        for (i=0;i<n;i++) {
			ssat[obs[i].sat-1].azel[0] = azel_[i*2]; // 方位角
			ssat[obs[i].sat-1].azel[1] = azel_[1+i*2]; // 高度角
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1; // 可用状态
            ssat[obs[i].sat-1].resp[0]=resp[i]; // 伪距残差
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp); // 后置数据清理
    return stat;
}
