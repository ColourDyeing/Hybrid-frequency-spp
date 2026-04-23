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
/* NX: ???11??mask[NX-3]=mask[8]?????????? */
/* IFLC??????: x[3]=GPS_BRDC(?????), x[4]=GLO???, x[5]=GAL???, x[6]=CMP???, x[7]=IRN???, x[8]=QZS???, x[9]=GPS_IFLC???, x[10]=GAL_IFLC??? */
/* ??IFLC??????: x[3]=GPS, x[4]=GLO, x[5]=GAL, x[6]=CMP, x[7]=IRN, x[8]=QZS */
#define NX          11          /* ???????????????? */
#define NMASK       9           /* mask????????IFLC???9??slot??GLO/GAL/CMP/IRN/QZS/GPS_IFLC/GAL_IFLC+GPS_BRDC???+1?????? */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
# define MAX_GDOP   30          /* max gdop for valid solution  */

/* ?????????????????????????????????------------------------------------
 args   :const prcopt_t *opt      I   ???????????
         const ssat_t   *ssat     I   ??????
         const obsd_t   *obs      I   ?????????
               double   el        I   ??????? (rad)
               int      sys       I   ?????? (SYS_???)
			   int      use_iflc  I   ??????IFLC???
               int      f2        I   IFLC????????????-1?????IFLC??

 return  : ????????? (m^2)
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
        snr_rover=obs->SNR[0]!=0?obs->SNR[0]:opt->err[5]; // ????????bug???????????????obs????snr
        varr+=SQR(opt->err[6])*pow(10,0.1*MAX(opt->err[5]-snr_rover,0));
    }

    /* ??????????????eratio */
    int er_ix=use_iflc?f2:0;
    varr*=SQR(opt->eratio[er_ix]);

    /* ?????????? */
    if (opt->err[7]>0.0) { 
        varr+=SQR(opt->err[7]*obs->Pstd[0]);
    }
    /* ????IFLC??????????????????????????? */
    if (use_iflc) {
        varr*=(f2==2)?SQR(2.2):SQR(3.0);
    }
    return SQR(fact)*varr;
}
/* ???????????????????TGD (s??m) ---------------------------------------------
args   :       int        sat     I   ??????
         const nav_t     *nav     I   ????????
		       int        type    I   TGD????   GPS/QZS:tgd[0]=TGD 
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
        /* ????BIA?????????????? */
        double tgd_bia=nav->tgd_bia[sat-1][type];
        if (tgd_bia!=0.0) return tgd_bia*CLIGHT; /* BIA-derived TGD: s -> m */
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* test SNR mask ???????????? -------------------------------------------------------------
args:    const obsd_t   *obs     I   ?????????
         const double   *azel    I   ???????????? {az,el} (rad)
		 const prcopt_t *opt     I   ???????????
		 int            use_iflc I   ??????IFLC???
         int            f2       I   IFLC??????????

return: 1:ok, 0: ???????
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
/* ????????????????????????????DCB??IFLC??TGD----------------------------------------
args:  const obsd_t      *obs     I   ?????????
       const nav_t       *nav     I   ????????
       const prcopt_t    *opt     I   ???????????
			 int         use_iflc I   ??????IFLC???
             int          f2      I   IFLC??????????
	         double      *var     O   ??????????DCB???????? (m^2)

return: ??????????DCB?????????? (m)
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

	if (P1==0.0) return 0.0; /* L1???????????? */

    /* SPP???????????????????????DCB?????????DCB?????????????? */
	/* L1 DCB??? */
	bias_ix=code2bias_ix(sys,obs->code[0]);  /* L1 DCB??? */
	if (bias_ix>0) { /* ???0????????????????? */
		P1+=nav->cbias[sat-1][0][bias_ix-1];
        trace(3,"prange: sat=%3d sys=%d code=%d bias_ix=%d dcb=%7.3f\n",sat,sys,obs->code[0],bias_ix,nav->cbias[sat-1][0][bias_ix-1]);
	}
    /* L2/L5 DCB??? */
    if (use_iflc) {
        int f2_bias_ix=code2bias_ix(sys,obs->code[f2]);
        int f2_freq;
        if (sys==SYS_GPS||sys==SYS_QZS) {
            f2_freq=f2==1?1:2; /* f2==1(L2)->freq1, f2!=1(L5)->freq2 */
        } else if (sys==SYS_GAL) {
            /* GAL: freq=0(E1), freq=1(E5a), freq=2(E5b), freq=3(E6) */
            /* RTKLIB??f2==1???E5b, f2!=1???E5a */
            f2_freq=f2==1?2:1;
        } else if (sys==SYS_CMP) {
            /* CMP: freq=0(B1), freq=1(B2), freq=2(B3), freq=3(B2a) */
            f2_freq=f2==1?1:3; /* ????B2->freq1, B2a->freq3 */
        } else {
            f2_freq=1; /* ??? */
        }
        if (f2_bias_ix>0) {
            P2+=nav->cbias[sat-1][f2_freq][f2_bias_ix-1];
            trace(3,"prange: sat=%3d sys=%d code=%d freq=%d bias_ix=%d dcb=%7.3f\n",sat,sys,obs->code[f2],f2_freq,f2_bias_ix,nav->cbias[sat-1][f2_freq][f2_bias_ix-1]);
        }
    }

    /* ????IFLC??TGD??? */
    if (use_iflc) {
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 or L1-L5 */
            if (f2==1) { /* L1-L2 */
                return (P2-SQR(FREQL1/FREQL2)*P1)/(1.0-SQR(FREQL1/FREQL2));
            }
            gamma=SQR(FREQL1/FREQL5);
            b1=gettgd(sat,nav,0); /* TGD_L1L2 (m) */
            b2=gettgd(sat,nav,5); /* TGD_L1L5 (m) */
            trace(3,"prange: sat=%2d L1-L5 gamma=%.6f b1=%.3f b2=%.3f TGD_corr=%.3f\n",
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
                /* I/NAV??????????E1/E5b??BGD_E1E5b = TGD */
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
        else if (sys==SYS_CMP) { /* ????????bug??B1-B2/ B1-B2a */
			gamma = SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQL1)/FREQ2_CMP); /* B1I???FREQ1_CMP??B1C???FREQL1; ????????????FREQ2_CMP??B2I/B2b?? */
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

    /* ?????????BRDC??IFLC??? */
    else {
        *var = SQR(ERR_CBIAS); /* ????????????(Code Bias??)??????????????????????????????? */

        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys == SYS_GLO) { /* G1 */
            return P1; /* GLONASS G1: ?????????????G1????????????? */
            //gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            //b1=gettgd(sat,nav,0); /* -dtaun (m) */
            //return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* F/NAV????E5a */
            else                    b1=gettgd(sat,nav,1); /* I/NAV????E5b */
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
/* ????????? ------------------------------------------------------
* args   : gtime_t time     I   ???
*          nav_t  *nav      I   ????????
*          int    sat       I   ??????
*          double *pos      I   ???????? {lat,lon,h} (rad|m)
*          double *azel     I   ?????/???? {az,el} (rad)
*          int    ionoopt   I   ??????????? (IONOOPT_???)
*          double *ion      O   ???????? (L1) (m)
*          double *var      O   ???????? (L1) ???? (m^2)
* return : ?? (1:???, 0:????)
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
    /* QZSS broadcast ionosphere model??Klobuchar?? */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* GPS broadcast ionosphere model??Klobuchar?? */
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
/* ????????? -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   ???
*          nav_t  *nav      I   ????????
*          double *pos      I   ???????? {lat,lon,h} (rad|m)
*          double *azel     I   ?????/???? {az,el} (rad)
*          int    tropopt   I   ??????????? (TROPOPT_???)
*          double *trp      O   ????????? (m)
*          double *var      O   ???????????? (m^2)
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
/* ?????????????? -----------------------------------------------------
    int      iter      I   ????????????estpos()????????????i????????i
    obsd_t   *obs      I   ?????????
    int      n         I   ??????????????
    double   *rs       I   ?????????????????6*n??{x,y,z,vx,vy,vz}(ecef)(m,m/s)
    double   *dts      I   ????????????2*n?? {bias,drift} (s|s/s)
    double   *vare     I   ????????????????? (m^2)
    int      *svh      I   ?????????? (-1:correction not available)
    nav_t    *nav      I   ????????
    double   *x        I   ?????????????????,7*1,?3????????????????????????4??????????????????gps????glonass??galileo??bds???????
    prcopt_t *opt      I   ???????????
    ssat_t   *ssat     I   ??????
    double   *v        O   ????????????????????
    double   *H        O   ????????????????
    double   *var      O   ????????????????
    double   *azel     O   ??????????????????????? {??????????} (2*n)
    int      *vsat     O   ???????????????????????? (1*n)
    double   *resp     O   ????????????????(P-(r+c*dtr-c*dts+I+T)) (1*n)
    int      *ns       O   ???????????????
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
	int i,j,nv=0,sat,sys,mask[NMASK]={0},use_iflc,f2; /* nv???????????????????????mask??????????????????????????????????? */

    /* ?????????????????????rr??dtr???? */
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];

    ecef2pos(rr,pos); // rr{x,y,z}->pos{lat,lon,h}
    trace(3,"iter=%d, rescode: x=%.3f y=%.3f z=%.3f\n", iter, rr[0], rr[1], rr[2]);

    /* ??????????????OBS[] */
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0; /* ???????????????????????????????????????????????????? */
        time=obs[i].time; /* ??? */
        sat=obs[i].sat; /* ?????? */
        if (!(sys=satsys(sat,NULL))) continue; /* 1.????satsys()???????????????????????????????????????? */
        f2=seliflc(opt->nf,sys);
        use_iflc=(opt->ionoopt==IONOOPT_IFLC&&obs[i].P[0]!=0.0&&obs[i].P[f2]!=0.0);

        
        /* ?????????????? */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            char tstr[40];
            trace(2,"duplicated obs data %s sat=%d\n",time2str(time,tstr,3),sat);
            i++;
            continue;
        }
        /* 2.???????????????? */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* 3-4.??????????????????????mask */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
        
        if (iter>0) { /* ?????????????????????????????????????????????? */
            /* 5.test SNR mask */
            if (!snrmask(obs+i,azel+i*2,opt,use_iflc,f2)) continue;

			/* 6.???????? */
			if (use_iflc) {
				/* IFLC???????????????????????????????????????????????? */
				dion=0.0;
				vion=0.0;
			} else {
				/* ??????????????????BRDC/SBAS/TEC????????IFLC???????????? */
				int iono_opt=(opt->ionoopt==IONOOPT_IFLC)?IONOOPT_BRDC:opt->ionoopt;
				if (!ionocorr(time,nav,sat,pos,azel+i*2,iono_opt,&dion,&vion)) {
					continue;
				}
            if ((freq=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
            /* Convert from FREQL1 to freq */
            dion*=SQR(FREQL1/freq);
            vion*=SQR(SQR(FREQL1/freq));
			}

            /* 7.????????? */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
        }
            }
        /* 8.??????????DCB??IFLC??TGD?? */
        if ((P=prange(obs+i,nav,opt,use_iflc,f2,&vmeas))==0.0) continue;

        /* 9.??????????????????? 
        /* ?????????(P - (r + c * dtr - c * dts + I + T)), ??????dtr????m */ 
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        trace(3,"sat=%2d: v=%.3f P=%.3f r=%.3f dtr_gps=%.6f dts=%.6f dion=%.3f dtrp=%.3f\n",
            sat,v[nv],P,r,dtr,dts[i*2],dion,dtrp);
       
        /* 10.?????? */ 
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0); /* ?3?????????????????????????????????????? */
        }
        /* ????????????????????? 
        IFLC??????: x[3]=GPS_BRDC, x[4]=GLO, x[5]=GAL, x[6]=CMP, x[7]=IRN, x[8]=QZS, x[9]=GPS_IFLC, x[10]=GAL_IFLC
        BRDC??????: x[3]=GPS, x[4]=GLO, x[5]=GAL, x[6]=CMP, x[7]=IRN, x[8]=QZS
        mask?????????GPS_BRDC/GLO/GAL/CMP/IRN/QZS/GPS_IFLC/GAL_IFLC+1???? */
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
        
        /* 11.??????????????????????????? */
        var[nv]=vare[i]+vmeas+vion+vtrp;
        if (ssat)
            var[nv++]+=varerr(opt,&ssat[i],&obs[i],azel[1+i*2],sys,use_iflc,f2);
        else
            var[nv++]+=varerr(opt,NULL,&obs[i],azel[1+i*2],sys,use_iflc,f2);
        trace(3,"        azel=%5.1f %4.1f res=%7.3f sig=%5.3f use_iflc=%d\n",
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]),use_iflc);
    }
    /* ???????????????????????????? */
    for (i=0;i<NMASK;i++) {
		if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
	return nv; /* ???????????? */ 
}
/* ???????????????????GDOP????---------------------------------------------------------
    const double   *azel     ??????????
    const int      *vsat     ????????????????????? (1*n)
          int      n         ????????
    const prcopt_t *opt      ???????
    const double   *v        ????????????????????
          int      nv        ??????
          int      nx        ???????????
          char     *msg      ???????
---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    /* ?????????? */
    vv=dot(v,v,nv);  // chisqr:???????
    if (nv>nx&&vv>chisqr[nv-nx-1]) {  /* ????????????????????  nv-nx-1:???????? */
        sprintf(msg,"Warning: large chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        /* return 0; */ /* ???????????????????????????????????? */
    }
    /* GDOP???? */
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
/* ??????????????????? ------------------------------------------------
    obsd_t   *obs      I   ?????????
    int      n         I   ??????????????
    double   *rs       I   ?????????????????6*n??{x,y,z,vx,vy,vz}(ecef)(m,m/s)
    double   *dts      I   ????????????2*n??{dt,dt_sap}(s)
    double   *vare     I   ????????????????? (m^2)
    int      *svh      I   ?????????? (-1:correction not available)
    nav_t    *nav      I   ????????
    prcopt_t *opt      I   ???????????
    ssat_t   *ssat     I   ??????
    sol_t    *sol      IO  ???
    double   *azel     IO  ??????????? (rad)
    int      *vsat     IO  ???????????????
    double   *resp     IO  ?????????? (P-(r+c*dtr-c*dts+I+T))
    char     *msg      O   ???????
retrun : 0:???????1:???????2:SBAS???????
------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, const ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+NX-3,1); H=mat(NX,n+NX-3); var=mat(n+NX-3,1); /* ????????????n + NX - 3?????????rescode???????????????????????????????????????????????
                                                           ???????????????????????????????????????????????????????????????????????????????NX-3 */
    
    for (i=0;i<3;i++) x[i]=sol->rr[i]; /* ????????????????????????????????????????????0?? */

    /* ?????????????? */
    for (i=0;i<MAXITR;i++) {

        /*  1.???? rescode ?????????????????????? v????????? H??
            ????????? var?????????????????????? azel??????????? vsat??
            ?????????? resp???????????????? ns ???????? nv */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,vsat,resp,
                   &ns);
        
        if (nv<NX) { /* ???????????????????????????????? */
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        /* ?????????????????????????H??v???????????????????????H??v */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* 2.????lsq(??????????)????,??????x???????dx????????????????????????Q */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) { // ??????????
            x[j]+=dx[j];
        }
        /* ???????????????????????????(????1E-4)????x[j]???????????????
         ?? sol ????????????,???????? valsol ????????????????????,??? RTKLIB Manual P162
         ??????????????????*/
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
            trace(3,"estpos  : x=%.3f y=%.3f z=%.3f dtr(m):gps_brdc=%.6f glo=%.6f gal=%.6f cmp=%.6f irn=%.6f qzs=%.6f gps_iflc=%.6f gal_iflc=%.6f\n",
                  x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10]);

            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;
            
            /* 3.???????????????????GDOP???? */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    /* ??????????????????????????????????return 0 */
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
* args   : obsd_t *obs      I   observation data OBS???????
*          int    n         I   number of observation data OBS????
*          nav_t  *nav      I   navigation data NAV????????????
*          prcopt_t *opt    I   processing options ???????????
*          sol_t  *sol      IO  solution ???
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output) ?????????
*          ssat_t *ssat     IO  satellite status              (NULL: no output) ??????
*          char   *msg      O   error message for error exit ???????
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
    
    if (n<=0) {     /* ???????????????0 */
        strcpy(msg,"no observation data");
        return 0;
    }
    sol->time=obs[0].time; /* sol->time????????????????? */
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
    
    if (opt_.mode!=PMODE_SINGLE) { /* ????????????SPP */
        opt_.ionoopt=IONOOPT_BRDC; /* ?????????Klobuchar????????? */
        opt_.tropopt=TROPOPT_SAAS; /* ?????????????Saastmoinen??? */
    }
    /* 1.??????????????????? */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* 2.??????????????????? */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    
    /* 3.??????????????????????? */
    if (!stat&&n>=6&&opt->posopt[4]) {      /* estpos??valsol????????????????????
                                            ?????RAIM????????????????????????
                                            ???????????>6?????????????????opt->posopt[4]=1??????????????????? */
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }
    /* 4.????????????????? */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];   /* ??????????? */
    }
    if (ssat) {     /* ?????????????ssat */
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
        }
        for (i=0;i<n;i++) {
			ssat[obs[i].sat-1].azel[0] = azel_[i*2]; /* ????? */
			ssat[obs[i].sat-1].azel[1] = azel_[1+i*2]; /* ???? */
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1; /* ?????? */
            ssat[obs[i].sat-1].resp[0]=resp[i]; /* ????? */
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp); /* ???????????? */
    return stat;
}
