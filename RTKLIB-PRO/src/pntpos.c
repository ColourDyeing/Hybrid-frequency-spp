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
/* NX: 统一11个，mask[NX-3]=mask[8]满足所有星座模式 */
/* IFLC模式参数: x[3]=GPS_BRDC(统一基准), x[4]=GLO钟差, x[5]=GAL钟差, x[6]=CMP钟差, x[7]=IRN钟差, x[8]=QZS钟差, x[9]=GPS_IFLC钟差, x[10]=GAL_IFLC钟差 */
/* 非IFLC模式参数: x[3]=GPS, x[4]=GLO, x[5]=GAL, x[6]=CMP, x[7]=IRN, x[8]=QZS */
#define NX          11          /* 统一最大待估参数数量 */
#define NMASK       9           /* mask数组大小，IFLC需要9个slot存放GLO/GAL/CMP/IRN/QZS/GPS_IFLC/GAL_IFLC+GPS_BRDC基准+1个留余 */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
# define MAX_GDOP   30          /* max gdop for valid solution  */

/* 计算伪距测量误差方差，用以随机模型定权 ------------------------------------
 args   :const prcopt_t *opt      I   处理过程选项
         const ssat_t   *ssat     I   卫星状态
         const obsd_t   *obs      I   观测数据
               double   el        I   卫星高度角 (rad)
               int      sys       I   卫星系统 (SYS_???)
			   int      use_iflc  I   是否使用IFLC组合
               int      f2        I   IFLC第二频率索引，-1表示非IFLC模式

 return  : 测量误差方差 (m^2)
------------------------------------*/
static double varerr(const prcopt_t *opt, const ssat_t *ssat, const obsd_t *obs, double el, int sys, int use_iflc, int f2)
{
    double fact=1.0,varr,snr_rover;

    switch (sys) {
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
        snr_rover=obs->SNR[0]!=0?obs->SNR[0]:opt->err[5]; // 修复bug：直接使用obs中的snr而非ssat中的snr
        varr+=SQR(opt->err[6])*pow(10,0.1*MAX(opt->err[5]-snr_rover,0));
    }

    /* 根据是否使用IFLC模式选择对应的eratio */
    int er_ix=use_iflc?f2:0;
    varr*=SQR(opt->eratio[er_ix]);

    /* 接收机标准差项 */
    if (opt->err[7]>0.0) {
        varr+=SQR(opt->err[7]*obs->Pstd[0]);
    }
    /* 乘以IFLC组合模式对应的噪声放大因子 */
    if (use_iflc) {
        varr*=(f2==2)?SQR(2.2):SQR(3.0);
    }
    return SQR(fact)*varr;
}
/* 获取TGD群延迟参数 (s -> m) ---------------------------------------------
args   :       int        sat     I   卫星号
         const nav_t     *nav     I   导航电文
		       int        type    I   TGD类型   GPS/QZS:tgd[0]=TGD
                                                GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b
                                                CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp
                                                tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd,tgd[5]=ISC_B2ad
return : TGD m
---------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);

	if (sys == SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        /* 先查找BIA文件提供的TGD改正数 */
        double tgd_bia=nav->tgd_bia[sat-1][type];
        if (tgd_bia!=0.0) return tgd_bia*CLIGHT; /* BIA-derived TGD: s -> m */
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* 测试信噪比SNR是否满足遮罩要求 -------------------------------------------------------------
args:    const obsd_t   *obs     I   观测数据
         const double   *azel    I   卫星方位角/高度角 {az,el} (rad)
		 const prcopt_t *opt     I   处理过程选项
		 int            use_iflc I   是否使用IFLC组合
         int            f2       I   IFLC第二频率索引

return: 1:通过, 0:未通过
------------------------------------------------------------- */
static int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt, int use_iflc, int f2)
{
    if (testsnr(0,0,azel[1],obs->SNR[0],&opt->snrmask)) {
        return 0;
    }
    if (use_iflc) {
        if (testsnr(0,f2,azel[1],obs->SNR[f2],&opt->snrmask)) return 0;
    }
    return 1;
}
/* 计算伪距并应用DCB和IFLC组合的TGD修正----------------------------------------
args:  const obsd_t      *obs     I   观测数据
       const nav_t       *nav     I   导航电文
       const prcopt_t    *opt     I   处理过程选项
			 int         use_iflc I   是否使用IFLC组合
             int          f2      I   IFLC第二频率索引
	         double      *var     O   输出伪距方差，包含DCB修正项方差 (m^2)

return: 伪距改正后值，已包含DCB和TGD改正 (m)
------------------------------------------------------------- */
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt, int use_iflc, int f2,
                     double *var)
{
    double P1,P2=0.0,gamma,b1,b2;
	int sat,sys,bias_ix;

	sat=obs->sat;
	sys=satsys(sat,NULL);
	P1=obs->P[0];
	if (use_iflc) P2=obs->P[f2];
	*var=0.0;

	if (P1==0.0) return 0.0; /* L1频点无有效数据 */

    /* SPP单点定位需要的所有DCB修正：P1-C1、P2-C2等DCB码偏差改正 */
	/* L1 DCB修正 */
	bias_ix=code2bias_ix(sys,obs->code[0]);  /* L1 DCB索引 */
	if (bias_ix>0) { /* 索引大于0表示存在有效的DCB改正 */
		P1+=nav->cbias[sat-1][0][bias_ix-1];
        trace(3,"prange: sat=%3d sys=%d code=%d bias_ix=%d dcb=%7.3f
",sat,sys,obs->code[0],bias_ix,nav->cbias[sat-1][0][bias_ix-1]);
	}
    /* L2/L5 DCB修正 */
    if (use_iflc) {
        int f2_bias_ix=code2bias_ix(sys,obs->code[f2]);
        int f2_freq;
        if (sys==SYS_GPS||sys==SYS_QZS) {
            f2_freq=f2==1?1:2; /* f2==1(L2)->freq1, f2!=1(L5)->freq2 */
        } else if (sys==SYS_GAL) {
            /* GAL: freq=0(E1), freq=1(E5a), freq=2(E5b), freq=3(E6) */
            /* RTKLIB约定f2==1代表E5b, f2!=1代表E5a */
            f2_freq=f2==1?2:1;
        } else if (sys==SYS_CMP) {
            /* CMP: freq=0(B1), freq=1(B2), freq=2(B3), freq=3(B2a) */
            f2_freq=f2==1?1:3; /* 北斗B2->freq1, B2a->freq3 */
        } else {
            f2_freq=1; /* 默认 */
        }
        if (f2_bias_ix>0) {
            P2+=nav->cbias[sat-1][f2_freq][f2_bias_ix-1];
            trace(3,"prange: sat=%3d sys=%d code=%d freq=%d bias_ix=%d dcb=%7.3f
",sat,sys,obs->code[f2],f2_freq,f2_bias_ix,nav->cbias[sat-1][f2_freq][f2_bias_ix-1]);
        }
    }

    /* 计算IFLC组合的TGD改正 */
    if (use_iflc) {
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 or L1-L5 */
            if (f2==1) { /* L1-L2 */
                return (P2-SQR(FREQL1/FREQL2)*P1)/(1.0-SQR(FREQL1/FREQL2));
            }
            gamma=SQR(FREQL1/FREQL5);
            b1=gettgd(sat,nav,0); /* TGD_L1L2 (m) */
            b2=gettgd(sat,nav,5); /* TGD_L1L5 (m) */
            trace(3,"prange: sat=%2d L1-L5 gamma=%.6f b1=%.3f b2=%.3f TGD_corr=%.3f
",
                  sat,gamma,b1,b2,(b2-gamma*b1));
            if (b1!=0.0||b2!=0.0) {
                return (P2-gamma*P1+(b2-gamma*b1))/(1.0-gamma);
            }
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 or G1-G3 */
            gamma=f2==1?SQR(FREQ1_GLO/FREQ2_GLO):SQR(FREQ1_GLO/FREQ3_GLO);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b or E1-E5a */
            if (f2==1) { /* E1-E5b (I/NAV) */
                /* I/NAV播发星历中E1/E5b的BGD_E1E5b即等于TGD */
                gamma=SQR(FREQL1/FREQE5b);
                b1=gettgd(sat,nav,1); /* BGD_E1E5b */
                if (b1!=0.0) {
                    return (P2-gamma*P1+gamma*b1)/(1.0-gamma);
                }
                return (P2-gamma*P1)/(1.0-gamma);
            } else { /* E1-E5a (F/NAV): E1-E5a is already ionosphere-free (f1?f5a), no TGD needed */
                gamma=SQR(FREQL1/FREQL5);
                return (P2-gamma*P1)/(1.0-gamma);
            }
        }
        else if (sys==SYS_CMP) { /* 修复BDS的bug：B1-B2/B1-B2a组合 */
			gamma = SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQL1)/FREQ2_CMP); /* B1I对应FREQ1_CMP，B1C对应FREQL1；频率比公式中FREQ2_CMP代表B2I/B2b的频率 */
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

    /* 非IFLC模式的单频BRDC消电离层组合 */
    else {
        *var = SQR(ERR_CBIAS); /* 代码偏差方差(Code Bias：P1-P2、P1-C1等码偏差)会导致测距误差，对应分配一定的方差 */

        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys == SYS_GLO) { /* G1 */
            return P1; /* GLONASS G1: GLONASS卫星不播发TGD参数，G1频点无需TGD改正 */
            //gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            //b1=gettgd(sat,nav,0); /* -dtaun (m) */
            //return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* F/NAV播发时取E5a */
            else                    b1=gettgd(sat,nav,1); /* I/NAV播发时取E5b */
            return P1-b1;
        }
        else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            return P1-b1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            //b1=gettgd(sat,nav,0); /* TGD (m) */
            //return P1-b1;
            gamma=SQR(FREQs/FREQL5);
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-gamma*b1;
        }
    }
    return P1;
}
/* 电离层延迟改正 ------------------------------------------------------
* args   : gtime_t time     I   观测时间
*          nav_t  *nav      I   导航电文
*          int    sat       I   卫星号
*          double *pos      I   接收机位置 {lat,lon,h} (rad|m)
*          double *azel     I   方位角/高度角 {az,el} (rad)
*          int    ionoopt   I   电离层改正模式 (IONOOPT_???)
*          double *ion      O   电离层延迟 (L1) (m)
*          double *var      O   电离层延迟方差 (L1) (m^2)
* return : 状态 (1:成功, 0:失败)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    int err=0;

    char tstr[40];
    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f
",
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
    /* QZSS广播电离层模型（Klobuchar模型） */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* GPS广播电离层模型（Klobuchar模型） */
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
/* 对流层延迟改正 -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   观测时间
*          nav_t  *nav      I   导航电文
*          double *pos      I   接收机位置 {lat,lon,h} (rad|m)
*          double *azel     I   方位角/高度角 {az,el} (rad)
*          int    tropopt   I   对流层改正模式 (TROPOPT_???)
*          double *trp      O   对流层延迟 (m)
*          double *var      O   对流层延迟方差 (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    char tstr[40];
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f
",
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
/* 计算伪距残差 -----------------------------------------------------
    int      iter      I   迭代次数，estpos()中第i次迭代时的第i次
    obsd_t   *obs      I   观测数据
    int      n         I   观测数据数量
    double   *rs       I   卫星位置速度，6*n数组{x,y,z,vx,vy,vz}(ecef)(m,m/s)
    double   *dts      I   卫星钟差，2*n数组 {bias,drift} (s|s/s)
    double   *vare     I   卫星位置误差方差 (m^2)
    int      *svh      I   卫星健康状态标识 (-1:correction not available)
    nav_t    *nav      I   导航电文
    double   *x        I   状态向量参数，11*1，前3个为位置xyz，第4个为接收机钟差，后续为各系统间钟差（gps/glonass/galileo/bds等系统的接收机钟差基准）
    prcopt_t *opt      I   处理过程选项
    ssat_t   *ssat     I   卫星状态
    double   *v        O   伪距残差向量（观测值减去计算值）
    double   *H        O   设计矩阵（偏导数矩阵）
    double   *var      O   残差方差
    double   *azel     O   卫星方位/高度角 {方位角,高度角} (2*n)
    int      *vsat     O   卫星有效标识：1表示卫星参与了定位解算 (1*n)
    double   *resp     O   伪距残差(P-(r+c*dtr-c*dts+I+T)) (1*n)
    int      *ns       O   有效观测数量
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
	int i,j,nv=0,sat,sys,mask[NMASK]={0},use_iflc,f2; /* nv为有效观测数，mask数组用于标识各系统是否有卫星参与定位，用于约束法方程防止秩亏 */

    /* 从状态向量x中提取接收机位置rr和钟差dtr */
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];

    ecef2pos(rr,pos); // rr{x,y,z}->pos{lat,lon,h}
    trace(3,"iter=%d, rescode: x=%.3f y=%.3f z=%.3f
", iter, rr[0], rr[1], rr[2]);

    /* 遍历所有OBS[]观测数据 */
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0; /* 初始化：先清零所有输出，后续对有效卫星重新赋值，避免残留值影响 */
        time=obs[i].time; /* 观测时刻 */
        sat=obs[i].sat; /* 卫星编号 */
        if (!(sys=satsys(sat,NULL))) continue; /* 1.调用satsys()获取卫星系统类型，若失败（sat=0）则跳过该卫星 */
        f2=seliflc(opt->nf,sys);
        use_iflc=(opt->ionoopt==IONOOPT_IFLC&&obs[i].P[0]!=0.0&&obs[i].P[f2]!=0.0);


        /* 检测重复观测数据（同一历元同一卫星的重复记录） */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            char tstr[40];
            trace(2,"duplicated obs data %s sat=%d
",time2str(time,tstr,3),sat);
            i++;
            continue;
        }
        /* 2.剔除健康状况异常的卫星（satellite health check） */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;

        /* 3-4.几何距离计算和高度角mask筛选 */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;

        if (iter>0) { /* 仅在iter>0时（非第一次迭代）执行精度相关的筛选和改正 */
            /* 5.信噪比SNR mask筛选 */
            if (!snrmask(obs+i,azel+i*2,opt,use_iflc,f2)) continue;

			/* 6.电离层延迟改正 */
			if (use_iflc) {
				/* IFLC组合本身就是消电离层组合，无需再施加电离层改正 */
				dion=0.0;
				vion=0.0;
			} else {
				/* 对于非IFLC卫星使用BRDC/SBAS/TEC模型计算电离层延迟，IFLC卫星也降级为BRDC以便统一处理 */
				int iono_opt=(opt->ionoopt==IONOOPT_IFLC)?IONOOPT_BRDC:opt->ionoopt;
				if (!ionocorr(time,nav,sat,pos,azel+i*2,iono_opt,&dion,&vion)) {
					continue;
				}
            if ((freq=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
            /* Convert from FREQL1 to freq */
            dion*=SQR(FREQL1/freq);
            vion*=SQR(SQR(FREQL1/freq));
			}

            /* 7.对流层延迟改正 */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
        }
            }
        /* 8.计算伪距，并应用DCB和IFLC组合的TGD改正 */
        if ((P=prange(obs+i,nav,opt,use_iflc,f2,&vmeas))==0.0) continue;

        /* 9.构建残差方程
        残差 = 伪距 - (几何距离 + 接收机钟差影响 - 卫星钟差影响 + 电离层延迟 + 对流层延迟)，其中dtr需要乘以光速c换算成m */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        trace(3,"sat=%2d: v=%.3f P=%.3f r=%.3f dtr_gps=%.6f dts=%.6f dion=%.3f dtrp=%.3f\n",
            sat,v[nv],P,r,dtr,dts[i*2],dion,dtrp);

        /* 10.构建设计矩阵 */
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0); /* 前3列对应xyz的偏导数，第4列对应接收机钟差偏导数为1.0 */
        }
        /* 构建设计矩阵的钟差列和约束mask
        IFLC模式: x[3]=GPS_BRDC, x[4]=GLO, x[5]=GAL, x[6]=CMP, x[7]=IRN, x[8]=QZS, x[9]=GPS_IFLC, x[10]=GAL_IFLC
        BRDC模式: x[3]=GPS, x[4]=GLO, x[5]=GAL, x[6]=CMP, x[7]=IRN, x[8]=QZS
        mask数组用于标识各系统是否有卫星参与，以约束法方程防止秩亏 */
        if (opt->ionoopt==IONOOPT_IFLC) {
            if (sys==SYS_GPS) {
                if (use_iflc) {v[nv]-=x[9]; H[9+nv*NX]=1.0; mask[6]=1;}
                else mask[0]=1;}
            else if (sys==SYS_GAL) {
                if (use_iflc) {v[nv]-=x[10]; H[10+nv*NX]=1.0; mask[7]=1;}
                else {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}}
            else if (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
            else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
            else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
            else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
        } else {
            if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
            else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
            else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
            else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
            else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
            else mask[0]=1;
        }

        vsat[i]=1; resp[i]=v[nv]; (*ns)++;

        /* 11.累计残差方差（卫星位置误差+伪距测量误差+电离层方差+对流层方差） */
        var[nv]=vare[i]+vmeas+vion+vtrp;
        if (ssat)
            var[nv++]+=varerr(opt,&ssat[i],&obs[i],azel[1+i*2],sys,use_iflc,f2);
        else
            var[nv++]+=varerr(opt,NULL,&obs[i],azel[1+i*2],sys,use_iflc,f2);
        trace(3,"        azel=%5.1f %4.1f res=%7.3f sig=%5.3f use_iflc=%d
",
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]),use_iflc);
    }
    /* 对未参与定位的系统添加虚拟观测约束，防止法方程秩亏 */
    for (i=0;i<NMASK;i++) {
		if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
	return nv; /* 返回有效观测数 */
}
/* 检验解算结果的有效性，计算GDOP --------------------------------------------------
    const double   *azel     卫星方位/高度角
    const int      *vsat     卫星有效标识：标识哪些卫星参与了定位解算 (1*n)
          int      n         观测数据总数
    const prcopt_t *opt      处理过程选项
    const double   *v        伪距残差向量
          int      nv        有效观测数
          int      nx        状态向量维度
          char     *msg      错误信息
---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;

    trace(3,"valsol  : n=%d nv=%d
",n,nv);

    /* 卡方检验 */
    vv=dot(v,v,nv);  // chisqr:卡方分布阈值
    if (nv>nx&&vv>chisqr[nv-nx-1]) {  /* 若残差卡方值超过阈值说明解算质量差  nv-nx-1:自由度 */
        sprintf(msg,"Warning: large chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        /* return 0; */ /* 暂不启用此严格检查，防止误删好解 */
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
/* 估计接收机位置 ------------------------------------------------
    obsd_t   *obs      I   观测数据
    int      n         I   观测数据数量
    double   *rs       I   卫星位置速度，6*n数组{x,y,z,vx,vy,vz}(ecef)(m,m/s)
    double   *dts      I   卫星钟差，2*n数组{dt,dt_sap}(s)
    double   *vare     I   卫星位置误差方差 (m^2)
    int      *svh      I   卫星健康状态标识 (-1:correction not available)
    nav_t    *nav      I   导航电文
    prcopt_t *opt      I   处理过程选项
    ssat_t   *ssat     I   卫星状态
    sol_t    *sol      IO  解算结果
    double   *azel     IO  方位/高度角 (rad)
    int      *vsat     IO  卫星有效标识
    double   *resp     IO  伪距残差 (P-(r+c*dtr-c*dts+I+T))
    char     *msg      O   错误信息
retrun : 0:失败  1:成功  2:SBAS解算
------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, const ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;

    trace(3,"estpos  : n=%d
",n);

    v=mat(n+NX-3,1); H=mat(NX,n+NX-3); var=mat(n+NX-3,1); /* 动态分配：n + NX - 3是因为rescode中会额外添加NX-3个虚拟约束观测，这部分空间也需分配
                                                           实际使用最多n个有效观测 + NX-3个约束 = n+NX-3个观测方程 */

    for (i=0;i<3;i++) x[i]=sol->rr[i]; /* 用上次的解作为本次迭代的初始值，若无可用初值则自动为0（第一次定位时） */

    /* 迭代求解 */
    for (i=0;i<MAXITR;i++) {

        /*  1.调用 rescode 计算伪距残差和设计矩阵，输出：v（残差向量）、H（设计矩阵）、var（方差）、azel（方位/高度角）、vsat（有效卫星）、resp（残差）、ns（有效卫星数）、nv（有效观测数） */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,vsat,resp,
                   &ns);

        if (nv<NX) { /* 有效观测数少于待估参数数，法方程秩亏无法求解 */
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        /* 对残差和设计矩阵进行方差归一化，即残差除以sigma，设计矩阵每列除以sigma */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* 2.调用lsq(最小二乘)求解，得到dx（状态增量）和Q（协方差矩阵） */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) { // 更新状态向量
            x[j]+=dx[j];
        }
        /* 当状态增量dx的模小于阈值(1E-4)时，认为收敛，更新sol解算结果，然后调用valsol做质量检核，具体参见 RTKLIB Manual P162
         收敛条件：dx的模小于1E-4米 */
        if (norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
#ifdef QZSDT
            sol->dtr[5]=x[8]/CLIGHT; /* QZS-GPS time offset (s) */
#endif
            if (opt->ionoopt==IONOOPT_IFLC) {
                sol->dtr[6]=x[9]/CLIGHT;  /* GPS_IFLC-GPS clock bias (s) */
                sol->dtr[7]=x[10]/CLIGHT; /* GAL_IFLC-GPS clock bias (s) */
            }
            trace(3,"estpos  : x=%.3f y=%.3f z=%.3f dtr(m):gps_brdc=%.6f glo=%.6f gal=%.6f cmp=%.6f irn=%.6f qzs=%.6f gps_iflc=%.6f gal_iflc=%.6f
",
                  x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10]);

            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;

            /* 3.用残差进行卡方检验和GDOP检核，通过后设置解状态 */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    /* 达到最大迭代次数仍未收敛，返回失败 */
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

    trace(3,"raim_fde: %s n=%2d
",time2str(obs[0].time,tstr,0),n);

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
            trace(3,"raim_fde: exsat=%2d (%s)
",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d
",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);

        trace(3,"raim_fde: exsat=%2d rms=%8.3f
",obs[i].sat,rms_e);

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
        trace(2,"%s: %s excluded by raim
",tstr+11,name);
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

    trace(3,"resdop  : n=%d
",n);

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
            trace(3,"estvel : vx=%.3f vy=%.3f vz=%.3f, n=%d
",x[0],x[1],x[2],n);
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
* args   : obsd_t *obs      I   观测数据 OBS数组指针
*          int    n         I   观测数据数量 OBS数组长度
*          nav_t  *nav      I   导航电文 NAV数组指针
*          prcopt_t *opt    I   处理过程选项
*          sol_t  *sol      IO  解算结果
*          double *azel     IO  方位/高度角 (rad) (NULL: 不输出)
*          ssat_t *ssat     IO  卫星状态信息              (NULL: 不输出)
*          char   *msg      O   错误信息输出
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
    trace(3,"pntpos  : tobs=%s n=%d
",time2str(obs[0].time,tstr,3),n);

    sol->stat=SOLQ_NONE;

    if (n<=0) {     /* 观测数据数量为0 */
        strcpy(msg,"no observation data");
        return 0;
    }
    sol->time=obs[0].time; /* sol->time先用观测时刻初始化，后续被估计的接收机钟差会修正它 */
    msg[0]=' ';
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

    if (opt_.mode!=PMODE_SINGLE) { /* 非单点定位模式（如PPP）时 */
        opt_.ionoopt=IONOOPT_BRDC; /* 强制使用Klobuchar广播模型进行电离层改正 */
        opt_.tropopt=TROPOPT_SAAS; /* 强制使用Saastamoinen模型进行对流层改正 */
    }
    /* 1.计算卫星位置、速度、钟差 */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);

    /* 2.估计接收机位置和钟差 */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);

    /* 3.若位置解算失败且满足条件，尝试RAIM/FDE（故障检测与排除） */
    if (!stat&&n>=6&&opt->posopt[4]) {      /* estpos和valsol都失败时启用
                                            RAIM算法逐个剔除卫星重新定位检测粗差
                                            要求观测数>6且opt->posopt[4]=1才启用RAIM算法 */
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }
    /* 4.估计接收机速度（使用多普勒观测值） */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];   /* 将方位/高度角结果复制到输出数组 */
    }
    if (ssat) {     /* 更新卫星状态信息到ssat */
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
        }
        for (i=0;i<n;i++) {
			ssat[obs[i].sat-1].azel[0] = azel_[i*2]; /* 方位角 */
			ssat[obs[i].sat-1].azel[1] = azel_[1+i*2]; /* 高度角 */
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1; /* 有效卫星 */
            ssat[obs[i].sat-1].resp[0]=resp[i]; /* 残差 */
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp); /* 释放动态分配的内存 */
    return stat;
}
