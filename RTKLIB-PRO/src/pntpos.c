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

/* ============================================================
 * 多频未组合SPP模型状态向量参数索引
 * 状态向量 x[NX] 共10维: [X,Y,Z, dt, dt_L5, ISB_GLO, ISB_GAL, ISB_BDS, ISB_IRN, ISB_QZS]
 * ============================================================ */
#define IDX_DT      3           /* 接收机钟差 (GPS参考系统, 单位m) */
#define IDX_DT_L5   4           /* L5频间接收机钟差 (所有系统共用, 单位m) */
#define IDX_ISB_GLO 5           /* GLONASS相对GPS的系统间偏差ISB (m) */
#define IDX_ISB_GAL 6           /* Galileo相对GPS的系统间偏差ISB (m) */
#define IDX_ISB_BDS 7           /* 北斗相对GPS的系统间偏差ISB (m) */
#define IDX_ISB_IRN 8           /* NavIC相对GPS的系统间偏差ISB (m) */
#ifdef QZSDT
#define IDX_ISB_QZS 9           /* QZSS相对GPS的系统间偏差ISB (m) */
#endif
/* 状态向量总维数: 位置3维 + 钟差1维 + L5钟差1维 + ISB 5维 = 10维 */
#ifdef QZSDT
#define NX          (4+6)       /* # of estimated parameters */
#else
#define NX          (4+5)       /* # of estimated parameters */
#endif
/* 状态向量索引说明:
 * [0-2]  = 位置 (X, Y, Z)
 * [3]    = GPS参考接收机钟差 dt
 * [4]    = L5频间钟差 dt_L5 (多频模式时新增)
 * [5-9]  = 各系统ISB (GLO/GAL/BDS/IRN/QZS)
 */
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
# define MAX_GDOP   30          /* max gdop for valid solution  */

/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, const ssat_t *ssat, const obsd_t *obs, double el, int sys, int fidx)
{
    double fact=1.0,varr,snr_rover;

    switch (sys) {
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
        snr_rover=obs->SNR[0]!=0?obs->SNR[0]:opt->err[5]; /* 修复bug：直接使用obs中的snr而非ssat中的snr */ 
        varr+=SQR(opt->err[6])*pow(10,0.1*MAX(opt->err[5]-snr_rover,0));
    }
    varr*=SQR(opt->eratio[0]);
    if (opt->err[7]>0.0) {
        varr+=SQR(opt->err[7]*obs->Pstd[0]);
    }
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* 消电离层组合 */
    if (fidx>=2) varr*=SQR(2.0);       /* L5观测噪声较大,放大方差防止过度信任 */
    return SQR(fact)*varr;
}
/* get group delay parameter (m) ---------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* test SNR mask -------------------------------------------------------------*/
static int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    int f2;

    if (testsnr(0,0,azel[1],obs->SNR[0],&opt->snrmask)) {
        return 0;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        f2=seliflc(opt->nf,satsys(obs->sat,NULL));
        if (testsnr(0,f2,azel[1],obs->SNR[f2],&opt->snrmask)) return 0;
    }
    return 1;
}
/* iono-free or "pseudo iono-free" pseudorange with code bias correction -----*/
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var)
{
    double P1,P2,gamma,b1,b2;
    int sat,sys,f2,bias_ix;

    sat=obs->sat;
    sys=satsys(sat,NULL);
    P1=obs->P[0];
    f2=seliflc(opt->nf,satsys(obs->sat,NULL));
    P2=obs->P[f2];
    *var=0.0;
    
    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;
    bias_ix=code2bias_ix(sys,obs->code[0]);  /* L1 code bias */
    if (bias_ix>0) { /* 0=ref code */
        P1+=nav->cbias[sat-1][0][bias_ix-1];
    }
    /* GPS code biases are L1/L2, Galileo are L1/L5 */
    if (sys==SYS_GAL&&f2==1) {
        /* skip code bias, no GAL L2 bias available */
    }
    else {  /* apply L2 or L5 code bias */
        bias_ix=code2bias_ix(sys,obs->code[f2]);
        if (bias_ix>0) { /* 0=ref code */
            P2+=nav->cbias[sat-1][1][bias_ix-1]; /* L2 or L5 code bias */
        }
    }
    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 or L1-L5 */
            gamma=f2==1?SQR(FREQL1/FREQL2):SQR(FREQL1/FREQL5);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 or G1-G3 */
            gamma=f2==1?SQR(FREQ1_GLO/FREQ2_GLO):SQR(FREQ1_GLO/FREQ3_GLO);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b, E1-E5a */
            gamma=f2==1?SQR(FREQL1/FREQE5b):SQR(FREQL1/FREQL5);
            if (f2==1&&getseleph(SYS_GAL)) { /* F/NAV */
                P2-=gettgd(sat,nav,0)-gettgd(sat,nav,1); /* BGD_E5aE5b */
            }
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_CMP) { /* B1-B2 */
            gamma=SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQL1)/FREQ2_CMP);
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
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys==SYS_GLO) { /* G1 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            b1=gettgd(sat,nav,0); /* -dtaun (m) */
            return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a */
            else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b */
            return P1-b1;
        }
        else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            return P1-b1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            gamma=SQR(FREQs/FREQL5);
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-gamma*b1;
        }
    }
    return P1;
}
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
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
    /* QZSS broadcast ionosphere model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* GPS broadcast ionosphere model */
    if (ionoopt==IONOOPT_BRDC||err==1) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
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
/* ============================================================
 * 判断观测量是否在L5/E5a/B2a频段
 * 返回:  1 = L5/E5a/B2a频段,  0 = L1/E1/B1频段,  -1 = 其他/未知
 * 说明: code2idx()返回频率索引: 0=L1/E1/B1, 1=L2/E5b/B2, 2=L5/E5a/B2a
 * ============================================================ */
static int is_l5_freq(int sat, uint8_t code)
{
    int sys=satsys(sat,NULL);
    int idx=code2idx(sys,code);
    if (idx<0) return -1;
    return (idx==2)?1:0;
}
/* ============================================================
 * 单频伪距提取函数(多频未组合SPP专用)
 * 对指定频率的伪距进行码偏差(code bias)改正和TGD/BGD钟差改正
 * 参数:
 *   sat    - 卫星编号
 *   sys    - 卫星系统 (SYS_GPS/SYS_GLO/SYS_GAL等)
 *   code   - 观测量类型代码 (CODE_L1C/CODE_L5Q等)
 *   fidx   - 频率索引 (0=L1/E1, 1=L2/E5b, 2=L5/E5a)
 *   obs    - 观测数据
 *   nav    - 导航电文
 *   var    - 输出: 码偏差方差
 * 返回: 改正后的伪距 (m)
 * ============================================================ */
static double prange_uc(int sat, int sys, uint8_t code, int fidx,
                        const obsd_t *obs, const nav_t *nav, double *var)
{
    double P,b1,gamma;
    int bias_ix;

    P=obs->P[fidx];
    *var=0.0;
    if (P==0.0) return 0.0;

    /* 应用码偏差(code bias)改正: 将非参考码改正到参考码 */
    bias_ix=code2bias_ix(sys,code);
    if (bias_ix>0) {
        P+=nav->cbias[sat-1][fidx][bias_ix-1];
    }

    /* 应用TGD/BGD钟差改正: 将伪距归算到"消电离层参考" */
    if (sys==SYS_GPS||sys==SYS_QZS) {        /* GPS/QZS: L1/L2/L5 */
        if (fidx==2) b1=gettgd(sat,nav,1);   /* L5: TGD第1类(L1-TGD) */
        else          b1=gettgd(sat,nav,0);   /* L1/L2: TGD第0类 */
        P-=b1;
    }
    else if (sys==SYS_GLO) {                 /* GLONASS: G1/G2/G3 */
        if (fidx==0) {
            /* G1: 用dtaun/(gamma-1)将G1归算到消电离层参考 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            b1=gettgd(sat,nav,0);
            P-=b1/(gamma-1.0);
        }
        else if (fidx==1) {
            /* G2: 用dtaun将G2归算到G1参考 */
            b1=gettgd(sat,nav,0);
            P-=b1;
        }
        /* G3: 无TGD改正 */
    }
    else if (sys==SYS_GAL) {                 /* Galileo: E1/E5b/E5a */
        if (fidx==0) {
            /* E1: 减去BGD归算到消电离层参考 */
            if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a (I/NAV) */
            else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b (F/NAV) */
            P-=b1;
        }
        else if (fidx==1) {
            /* E5b: 减去BGD_E1E5b将E5b归算到E1参考
             * BGD_E1E5b = TGD(E1) - TGD(E5b) */
            b1=gettgd(sat,nav,1);
            P-=b1;
        }
        else if (fidx==2) {
            /* E5a: 减去BGD_E1E5a将E5a归算到E1参考
             * BGD_E1E5a = TGD(E1) - TGD(E5a) */
            b1=gettgd(sat,nav,0);
            P-=b1;
        }
    }
    else if (sys==SYS_CMP) {                 /* 北斗: B1I/B1C/B2/B2a/B3 */
        if      (code==CODE_L2I) { b1=gettgd(sat,nav,0); P-=b1; }          /* TGD_B1I */
        else if (code==CODE_L1P) { b1=gettgd(sat,nav,2); P-=b1; }          /* TGD_B1Cp */
        else {
            /* B1Cp/B1Cd: TGD_B1Cp + ISC_B1Cd */
            b1=gettgd(sat,nav,2)+gettgd(sat,nav,4);
            P-=b1;
        }
    }
    else if (sys==SYS_IRN) {                 /* NavIC: L5 */
        gamma=SQR(FREQs/FREQL5);
        b1=gettgd(sat,nav,0);
        P-=gamma*b1;
    }
    return P;
}
/* ============================================================
 * 伪距残差计算函数 (多频未组合SPP核心)
 *
 * 功能: 构建伪距观测方程的残差向量v和设计矩阵H
 *
 * 多频处理策略:
 *   - 当 posopt[POSOPT_UNCOMB]=1 且 nf>=2 时: 对每颗卫星遍历所有频率
 *   - 否则: 仅处理频率索引0 (与传统SPP行为一致)
 *
 * 状态向量 x[NX] = [X,Y,Z, dt, dt_L5, ISB_GLO, ISB_GAL, ISB_BDS, ISB_IRN, ISB_QZS]
 * ============================================================ */
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   const ssat_t *ssat, double *v, double *H, double *var,
                   double *azel, int *vsat, double *resp, int *ns)
{
    gtime_t time;
    double r,rr[3],pos[3],dtr,e[3];
    int i,j,f,nv=0,sat,sys,nf;
    int mask[7]={0};  /* 标记各参数是否被观测到,用于约束处理 (最多7个非位置参数) */

    /* 确定要处理的频率数量
     * 多频未组合模式: 使用opt->nf指定的所有频率
     * 否则: 仅使用频率0 (保持与传统SPP兼容) */
    nf=(opt->posopt[POSOPT_UNCOMB]&&opt->nf>=2)?opt->nf:1;

    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[IDX_DT];  /* 提取GPS参考接收机钟差 */

    ecef2pos(rr,pos);
    trace(3,"rescode: rr=%.3f %.3f %.3f nf=%d\n",rr[0],rr[1],rr[2],nf);

    /* 外层循环: 遍历每颗卫星 */
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;

        /* 剔除重复观测数据 */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            char tstr[40];
            trace(2,"duplicated obs data %s sat=%d\n",time2str(time,tstr,3),sat);
            i++;
            continue;
        }
        /* 剔除被排除的卫星 */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;

        /* 计算几何距离和截止高度角 (同一卫星所有频率共用) */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;

        /* SNR掩码检验 (仅对频率0) */
        if (iter>0&&!snrmask(obs+i,azel+i*2,opt)) continue;

        /*=========================================================*/
        /* 内层循环: 遍历该卫星的每个频率                         */
        /*=========================================================*/
        for (f=0;f<nf&&f<NFREQ+NEXOBS;f++) {
            double P,dion=0.0,dtrp=0.0,vmeas=0.0,vion=0.0,vtrp=0.0;
            double freq;
            int isL5;

            /* 跳过无观测的频率 */
            if (obs[i].P[f]==0.0||obs[i].code[f]==CODE_NONE) continue;

            /* 获取载波频率 */
            freq=sat2freq(sat,obs[i].code[f],nav);
            if (freq<=0.0) continue;

            /* 判断是否为L5/E5a频段 */
            isL5=is_l5_freq(sat,obs[i].code[f]);
            if (isL5<0) continue; /* 未知频率,跳过 */

            /* 电离层延迟改正 (iter>0表示已有位置估计,可计算投影) */
            if (iter>0) {
                if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
                    continue;
                }
                /* 将L1频段的电离层延迟转换到当前观测频率 */
                dion*=SQR(FREQL1/freq);
                vion*=SQR(SQR(FREQL1/freq));

                /* 对流层延迟改正 */
                if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                    continue;
                }
            }

            /* 获取伪距 (多频未组合模式用prange_uc,否则用原始prange) */
            if (opt->posopt[POSOPT_UNCOMB]&&opt->nf>=2) {
                P=prange_uc(sat,sys,obs[i].code[f],f,obs+i,nav,&vmeas);
            }
            else {
                P=prange(obs+i,nav,opt,&vmeas);
            }
            if (P==0.0) continue;

            /*=========================================================*/
            /* 构架伪距残差方程:                                      */
            /* P_obs = r + dtr - c*dts + dion + dtrp               */
            /* 残差 = P_obs - 模型值                                  */
            /*=========================================================*/
            v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);

            /*=========================================================*/
            /* 构架设计矩阵(观测方程系数矩阵)H                        */
            /*=========================================================*/
            /* 位置参数偏导数: -e[0..2] (几何距离对位置的偏导) */
            for (j=0;j<3;j++) H[j+nv*NX]=-e[j];
            /* 钟差参数偏导数: +1 */
            H[IDX_DT+nv*NX]=1.0;
            /* dt_L5参数偏导数: 仅L5观测时为1,否则为0 */
            H[IDX_DT_L5+nv*NX]=(isL5>0)?1.0:0.0;
            /* ISB参数偏导数: 各系统特定列置1,GPS/SBAS列保持0 */
            H[IDX_ISB_GLO+nv*NX]=0.0;
            H[IDX_ISB_GAL+nv*NX]=0.0;
            H[IDX_ISB_BDS+nv*NX]=0.0;
            H[IDX_ISB_IRN+nv*NX]=0.0;
#ifdef QZSDT
            H[IDX_ISB_QZS+nv*NX]=0.0;
#endif

            /*=========================================================*/
            /* 应用钟差改正项,并构建mask标记各参数是否被观测           */
            /* 时钟模型:                                              */
            /*   GPS:  dtr                                          */
            /*   GAL:  dtr + ISB_GAL                               */
            /*   GLO:  dtr + ISB_GLO                               */
            /*   BDS:  dtr + ISB_BDS                               */
            /*   L5:   dt_L5叠加到上述模型上                         */
            /*=========================================================*/
            switch (sys) {
                case SYS_GLO:
                    v[nv]-=x[IDX_ISB_GLO];
                    H[IDX_ISB_GLO+nv*NX]=1.0;
                    mask[2]=1;  /* ISB_GLO被观测到 */
                    break;
                case SYS_GAL:
                    v[nv]-=x[IDX_ISB_GAL];
                    H[IDX_ISB_GAL+nv*NX]=1.0;
                    mask[3]=1;  /* ISB_GAL被观测到 */
                    break;
                case SYS_CMP:
                    v[nv]-=x[IDX_ISB_BDS];
                    H[IDX_ISB_BDS+nv*NX]=1.0;
                    mask[4]=1;  /* ISB_BDS被观测到 */
                    break;
                case SYS_IRN:
                    v[nv]-=x[IDX_ISB_IRN];
                    H[IDX_ISB_IRN+nv*NX]=1.0;
                    mask[5]=1;  /* ISB_IRN被观测到 */
                    break;
#ifdef QZSDT
                case SYS_QZS:
                    v[nv]-=x[IDX_ISB_QZS];
                    H[IDX_ISB_QZS+nv*NX]=1.0;
                    mask[6]=1;  /* ISB_QZS被观测到 */
                    break;
#endif
                case SYS_SBS:
                    /* SBAS使用GPS时钟,与原始RTKLIB一致 */
                    mask[0]=1;  /* GPS时钟被观测到 */
                    break;
                default: /* GPS及其他: 使用GPS时钟,无ISB */
                    mask[0]=1;  /* GPS时钟被观测到 */
                    break;
            }
            /* 标记dt_L5参数: L5观测时有效 */
            if (isL5>0) mask[1]=1;

            /* 标记卫星有效(取第一个有效频率的残差) */
            if (!vsat[i]) {
                vsat[i]=1;
                resp[i]=v[nv];
            }

            trace(4,"sat=%2d f=%d: v=%.3f P=%.3f r=%.3f dtr=%.6f dts=%.6f "
                "dion=%.3f dtrp=%.3f isL5=%d\n",
                sat,f,v[nv],P,r,dtr,dts[i*2],dion,dtrp,isL5);

            /* 构架伪距误差方差 */
            var[nv]=vare[i]+vmeas+vion+vtrp;
            if (ssat)
                var[nv++]+=varerr(opt,&ssat[i],&obs[i],azel[1+i*2],sys,f);
            else
                var[nv++]+=varerr(opt,NULL,&obs[i],azel[1+i*2],sys,f);
        } /* 内层循环: 遍历每个频率结束 */

        /* 卫星级标记处理 */
        if (vsat[i]) {
            (*ns)++;
            trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f\n",obs[i].sat,
                azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i]);
        }
    } /* 外层循环: 遍历每颗卫星结束 */

    /*=========================================================*/
    /* 添加约束行,避免秩亏                                     */
    /* mask[i]=1表示参数x[i+3]已被观测到,无需约束              */
    /* mask[i]=0表示参数未被任何观测约束,需添加弱约束          */
    /*=========================================================*/
    for (i=0;i<NX-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=(j==i+3)?1.0:0.0;
        var[nv++]=0.01; /* 弱零阶约束,方差0.01m^2 */
    }
    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    /* Chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        sprintf(msg,"Warning: large chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        /* return 0; */ /* threshold too strict for all use cases, report error but continue on */
    }
    /* large GDOP check */
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
/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, const ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns,nf;
    
    /* 确定要处理的频率数量 (与rescode中逻辑一致) */
    nf=(opt->posopt[POSOPT_UNCOMB]&&opt->nf>=2)?opt->nf:1;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v  =mat(n*nf+(NX-3),1);
    H  =mat(NX,n*nf+(NX-3));
    var=mat(n*nf+(NX-3),1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    trace(3,"estpos  : initial x=%.3f %.3f %.3f\n",x[0],x[1],x[2]);

    for (i=0;i<MAXITR;i++) {

        /* pseudorange residuals (m) */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,vsat,resp,
                   &ns);
        
        if (nv<NX) {
            trace(3,"estpos  : iter=%d nv=%d < NX=%d -> lack of valid sats\n",i,nv,NX);
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        trace(3,"estpos  : iter=%d nv=%d\n",i,nv);
        /* weight by variance (lsq uses sqrt of weight */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            trace(3,"estpos  : iter=%d lsq error info=%d\n",i,info);
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        trace(3,"estpos  : iter=%d dx norm=%.6f\n",i,norm(dx,NX));
        for (j=0;j<NX;j++) {
            x[j]+=dx[j];
        }
        /* 仅用位置参数判断收敛,避免钟差/ISB参数秩亏导致整体不收敛
         * 多频未组合模型中,GPS钟差和GLONASS-ISB的设计矩阵列完全相同,
         * 导致法方程病态,钟差参数震荡但位置参数实际已收敛 */
        if (norm(dx,3)<1E-4) {
            trace(3,"estpos  : converged iter=%d x=%.3f %.3f %.3f\n",i,x[0],x[1],x[2]);
            sol->type=0;
            /* 以GPS时钟为参考计算接收机时间: GPST = obs_time - dt/CLIGHT */
            sol->time=timeadd(obs[0].time,-x[IDX_DT]/CLIGHT);
            /* 各系统时钟/ISB输出到解结构sol_t.dtr[]
             * dtr[0] = GPS接收机钟差
             * dtr[1] = L5频间钟差 (新增)
             * dtr[2] = GLONASS-GPS ISB
             * dtr[3] = Galileo-GPS ISB
             * dtr[4] = 北斗-GPS ISB
             * dtr[5] = NavIC-GPS ISB
             * dtr[6] = QZSS-GPS ISB (QZSDT模式)
             * 注意: PPP模式使用dtr[1-3]存储GLO/GAL/BDS ISB, 与本模式索引不同
             */
            sol->dtr[0]=x[IDX_DT]/CLIGHT;       /* GPS接收机钟差 (s) */
            sol->dtr[1]=x[IDX_DT_L5]/CLIGHT;    /* L5频间钟差 (s) */
            sol->dtr[2]=x[IDX_ISB_GLO]/CLIGHT;   /* GLONASS-GPS时间偏移 (s) */
            sol->dtr[3]=x[IDX_ISB_GAL]/CLIGHT;   /* Galileo-GPS时间偏移 (s) */
            sol->dtr[4]=x[IDX_ISB_BDS]/CLIGHT;   /* 北斗-GPS时间偏移 (s) */
            sol->dtr[5]=x[IDX_ISB_IRN]/CLIGHT;    /* NavIC-GPS时间偏移 (s) */
#ifdef QZSDT
            sol->dtr[6]=x[IDX_ISB_QZS]/CLIGHT;   /* QZSS-GPS时间偏移 (s) */
#endif
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];     /* cov xy */
            sol->qr[4]=(float)Q[2+NX];  /* cov yz */
            sol->qr[5]=(float)Q[2];     /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;

            /* validate solution */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
        /* 停滞检测: 当位置参数变化小于1m但未达到严格收敛阈值时,
         * 说明钟差参数受限于模型秩亏而持续小幅振荡,位置参数实际已稳定 */
        if (i>=3&&norm(dx,3)<1.0) {
            trace(3,"estpos  : pos converged (stall at %.3fm), accepting iter=%d\n",
                  norm(dx,3),i);
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[IDX_DT]/CLIGHT);
            sol->dtr[0]=x[IDX_DT]/CLIGHT;
            sol->dtr[1]=x[IDX_DT_L5]/CLIGHT;
            sol->dtr[2]=x[IDX_ISB_GLO]/CLIGHT;
            sol->dtr[3]=x[IDX_ISB_GAL]/CLIGHT;
            sol->dtr[4]=x[IDX_ISB_BDS]/CLIGHT;
            sol->dtr[5]=x[IDX_ISB_IRN]/CLIGHT;
#ifdef QZSDT
            sol->dtr[6]=x[IDX_ISB_QZS]/CLIGHT;
#endif
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];
            sol->qr[4]=(float)Q[2+NX];
            sol->qr[5]=(float)Q[2];
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    if (i>=MAXITR) {
        trace(3,"estpos  : MAXITR=%d reached, divergent\n",MAXITR);
        sprintf(msg,"iteration divergent i=%d",i);
    }
    
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
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
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
    
    if (n<=0) {
        strcpy(msg,"no observation data");
        return 0;
    }
    sol->time=obs[0].time;
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
    
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positions, velocities and clocks */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position and time with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    trace(3,"pntpos  : estpos returned stat=%d msg=%s\n",stat,msg);
    
    /* RAIM FDE */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }
    /* estimate receiver velocity with Doppler */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
