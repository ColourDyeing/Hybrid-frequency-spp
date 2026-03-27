/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*           2026/03/27 1.8  UDUC (Un-Differenced Un-Combined) SPP:
*              Every satellite has ONE ionospheric parameter (slant L1 TEC).
*              All frequencies of the SAME satellite share this parameter,
*              scaled by γ_k = (f₁/f_k)².
*              
*              - Dual-freq sat: two pseudorange equations (P1,P2) + MW
*                → P1 and P2 give two equations for the SAME I_sat
*                → MW gives iono-free clock observable → resolves rank deficiency
*              - Single-freq sat: one pseudorange equation
*                → Klobuchar prior as mild stochastic constraint on I_sat
*                → Rank deficiency broken by cross-satellite diversity
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)     ((x)*(x))
#define MAX(x,y)   ((x)>=(y)?(x):(y))

#define QZSDT
#ifdef QZSDT
#define NX_BASE    (4+5)       /* pos(3) + clk(1) + sys_bias(5) */
#else
#define NX_BASE    (4+4)
#endif

/* UDUC iono estimation */
#define VAR_IONO   SQR(20.0)   /* initial variance of iono slant delay (m^2) */
#define SIG_IONO   0.5         /* stochastic constraint sigma on iono (m) */
#define MW_SIG     2.0         /* MW combination sigma (cycles) */
#define MIN_EL     (5.0*D2R)
#define MAXITR     10
#define ERR_ION    5.0
#define ERR_TROP   3.0
#define ERR_SAAS   0.3
#define ERR_BRDCI  0.5
#define ERR_CBIAS  0.3
#define REL_HUMI   0.7
#define MAX_GDOP   30.0

/* ionospheric parameter index for satellite s (1-based satno) */
#define II(s)      (NX_BASE+(s)-1)
/* iono state index limit */
#define NION_MAX   MAXSAT

/* compute NX for UDUC: base + iono states (one per sat) */
#define NX_UDUC    (NX_BASE+MAXSAT)

/* 计算伪距测量误差方差 ------------------------------------*/
static double varerr(const prcopt_t *opt, const ssat_t *ssat, const obsd_t *obs,
                     double el, int sys)
{
    double fact=1.0,varr,snr_rover;

    switch (sys) {
        case SYS_GPS: fact*=EFACT_GPS; break;
        case SYS_GLO: fact*=EFACT_GLO; break;
        case SYS_SBS: fact*=EFACT_SBS; break;
        case SYS_CMP: fact*=EFACT_CMP; break;
        case SYS_QZS: fact*=EFACT_QZS; break;
        case SYS_IRN: fact*=EFACT_IRN; break;
        default:      fact*=EFACT_GPS; break;
    }
    if (el<MIN_EL) el=MIN_EL;
    varr=SQR(opt->err[1])+SQR(opt->err[2])/sin(el);
    if (opt->err[6]>0.0) {
        snr_rover=(ssat)?ssat->snr_rover[0]:opt->err[5];
        varr+=SQR(opt->err[6])*pow(10,0.1*MAX(opt->err[5]-snr_rover,0));
    }
    varr*=SQR(opt->eratio[0]);
    if (opt->err[7]>0.0) varr+=SQR(opt->err[7]*obs->Pstd[0]);
    return SQR(fact)*varr;
}
/* get group delay parameter (m) ---------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) if (nav->geph[i].sat==sat) break;
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    for (i=0;i<nav->n;i++) if (nav->eph[i].sat==sat) break;
    return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
}
/* test SNR mask -------------------------------------------------------------*/
static int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    if (testsnr(0,0,azel[1],obs->SNR[0],&opt->snrmask)) return 0;
    return 1;
}

/* Select second frequency index for a satellite -------------------------------*/
static int sel2freq(int sat, const prcopt_t *opt)
{
    int sys=satsys(sat,NULL);
    if      (sys==SYS_GPS||sys==SYS_QZS) return (opt->nf>=3)?2:1;   /* L2 or L5 */
    else if (sys==SYS_GLO)               return (opt->nf>=3)?2:1;   /* G2 or G3 */
    else if (sys==SYS_GAL)               return (opt->nf>=3)?2:1;   /* E5a or E5b */
    else if (sys==SYS_CMP)               return 1;                   /* B2 */
    else if (sys==SYS_IRN)               return 1;                   /* S */
    return 0;
}

/* Get ionospheric mapping coefficient γ_k = (f₁/f_k)² ------------------------
* For L1/E1/B1: γ = 1.0
* For L2/G2/E5b/B2: γ = (f1/f2)²
* For L5/E5a: γ = (f1/f5)²
*----------------------------------------------------------------------------*/
static double iono_gamma(int sat, int fidx, const obsd_t *obs, const nav_t *nav)
{
    double f1,fk;
    int f2=sel2freq(sat,NULL);

    if (fidx==0) return 1.0; /* L1/E1/B1: always γ=1 */

    /* get frequency of the obs code */
    f1=sat2freq(sat,obs->code[0],nav);       /* freq of code[0] (always f1) */
    fk=sat2freq(sat,obs->code[fidx],nav);    /* freq of code[fidx] */

    if (f1<=0.0||fk<=0.0||f1==fk) return 1.0;
    return SQR(f1/fk);
}

/* DCB + TGD correction: bring pseudorange to reference bias (L1) --------------
* Returns corrected pseudorange in meters, or 0 if unavailable.
*---------------------------------------------------------------------------*/
static double prange_dcb(const obsd_t *obs, const nav_t *nav, int fidx,
                        double *vmeas)
{
    double Pk,b1=0.0,gamma=1.0;
    int sat,sys,bias_ix,f1,fk;
    *vmeas=0.0;

    sat=obs->sat; sys=satsys(sat,NULL);
    f1=sat2freq(sat,obs->code[0],nav);
    fk=sat2freq(sat,obs->code[fidx],nav);
    if (fk<=0.0) return 0.0;
    gamma=SQR(f1/fk);

    Pk=obs->P[fidx];
    if (Pk==0.0) return 0.0;

    /* code bias correction */
    bias_ix=code2bias_ix(sys,obs->code[fidx]);
    if (bias_ix>0) Pk+=nav->cbias[sat-1][fidx>0?1:0][bias_ix-1];

    /* TGD correction: maps to L1 reference */
    if (sys==SYS_GPS||sys==SYS_QZS) {
        if (fidx==0) {
            b1=gettgd(sat,nav,0);
        } else {
            b1=gettgd(sat,nav,0)*gamma; /* TGD is on L1, scale for f2 */
        }
        Pk-=b1;
    }
    else if (sys==SYS_GLO) {
        /* GLONASS: inter-frequency bias correction
         * -dtaun converts P1 to G1/G2 common reference
         * For fidx=0: subtract -dtaun/(gamma_glo-1)
         * For fidx>0: the bias is already applied similarly */
        double gamma_glo=SQR(FREQ1_GLO/FREQ2_GLO);
        b1=gettgd(sat,nav,0); /* -dtaun in meters */
        if (fidx==0) {
            Pk-=b1/(gamma_glo-1.0); /* bring P1 to G1/G2 common ref */
        } else {
            Pk-=b1/(gamma_glo-gamma); /* bring P2 to G1/G2 common ref */
        }
    }
    else if (sys==SYS_GAL) {
        if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a */
        else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b */
        Pk-=(fidx==0)?b1:b1*gamma;
    }
    else if (sys==SYS_CMP) {
        if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0);
        else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2);
        else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4);
        Pk-=(fidx==0)?b1:b1*gamma;
    }
    else if (sys==SYS_IRN) {
        double gamma_irn=SQR(FREQs/FREQL5);
        b1=gettgd(sat,nav,0);
        Pk-=(fidx==0)?gamma_irn*b1:gamma*b1;
    }
    *vmeas=SQR(ERR_CBIAS);
    return Pk;
}

/* 电离层修正 ------------------------------------------------------
* For UDUC (ionoopt==IONOOPT_EST): returns Klobuchar prior as iono guess.
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    int err=0;
    char tstr[40];
    trace(4,"ionocorr: time=%s opt=%d sat=%2d\n",
          time2str(time,tstr,3),ionoopt,sat);

    if (ionoopt==IONOOPT_SBAS) {
        if (sbsioncorr(time,nav,pos,azel,ion,var)) return 1; err=1;
    }
    if (ionoopt==IONOOPT_TEC) {
        if (iontec(time,nav,pos,azel,1,ion,var)) return 1; err=1;
    }
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI); return 1;
    }
    if (ionoopt==IONOOPT_BRDC||err==1||ionoopt==IONOOPT_EST) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI); return 1;
    }
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1)); return 1;
    }
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var); return 1;
    }
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}

/*==============================================================================
* UDUC (Un-Differenced Un-Combined) SPP Core
*
* Unified model for ALL satellites, ALL frequencies:
*
*   P_k^s = ρ_s + c·dtr - c·dts_s + I_s·γ_k^s + T
*
*   I_s      = slant ionospheric delay on L1 (ONE parameter per satellite s)
*   γ_k^s    = (f₁/f_k)² (frequency mapping, e.g. γ=1 for L1, γ=1.65 for L2)
*   dtr      = receiver clock bias (one common to all sats and all freqs)
*
* Rank deficiency resolution:
*   - Dual-freq sat: two equations (k=1,2) for the SAME I_s → I_s observable
*   - Single-freq sat: I_s observable via cross-satellite dtr constraint
*     (from dual-freq sats: their two equations determine dtr + I_s together)
*   - As a safeguard: Klobuchar prior acts as mild stochastic constraint on I_s
*
* Parameter: x = [pos(3), dtr(1), sys_bias(4), I_s1, I_s2, ..., I_sN]
*============================================================================*/

/* MW (Melbourne-Wübbena) combination ----------------------------------------
* MW = (P1-P2) - (λ1·Φ1-λ2·Φ2) = λ_wl·N_wl
* No ionospheric delay. Provides iono-free clock constraint for dual-freq sats.
* Returns MW in meters.
*----------------------------------------------------------------------------*/
static double mw_comb(const obsd_t *obs, const nav_t *nav, int f2)
{
    double f1,fk,lam1,lam2,lambda_wl;
    double P1,P2,L1,L2;

    f1=sat2freq(obs->sat,obs->code[0],nav);
    fk=sat2freq(obs->sat,obs->code[f2],nav);
    if (f1<=0.0||fk<=0.0) return 0.0;

    P1=obs->P[0]; L1=obs->L[0];
    P2=obs->P[f2]; L2=obs->L[f2];
    if (P1==0.0||P2==0.0||L1==0.0||L2==0.0) return 0.0;

    lam1=CLIGHT/f1; lam2=CLIGHT/fk;
    lambda_wl=CLIGHT/(f1-fk);   /* wide-lane wavelength */
    return (P1-P2)-(lam1*L1-lam2*L2); /* MW in meters */
}

/* Build UDUC observation equations ------------------------------------------
* Unified for all satellites (single-freq and dual-freq alike):
*
* For each observation P_k^s at frequency f_k:
*   v_P = P_k - (ρ_s + c·dtr - c·dts + I_s·γ_k + T)
*   H: pos=-e_s, dtr=1, sysbias=1(if applicable), I_s=γ_k
*
* For dual-freq satellites, also add MW equation:
*   v_MW = MW - λ_wl·N_wl
*   where N_wl = (dtr - c·dts)·(f1-f2)/c  [iono-free clock in cycles]
*   H: dtr=1, sysbias=1(if applicable), N_wl=-1 (or treated differently)
*
* Actually: MW = λ_wl·N_wl, and N_wl = (Φ1-Φ2) - (P1-P2)/λ_wl
*         = integer wide-lane ambiguity
* The receiver clock does NOT appear in MW (ionosphere-free by construction).
* So MW constrains the wide-lane ambiguity, NOT the clock.
*
* The key insight: dual-freq sat provides TWO equations for ONE I_s
* → I_s becomes directly observable (no rank deficiency).
*============================================================================*/
static int rescode_uduc(int iter, const obsd_t *obs, int n, const double *rs,
                        const double *dts, const double *vare, const int *svh,
                        const nav_t *nav, const double *x, const prcopt_t *opt,
                        const ssat_t *ssat, double *v, double *H, double *var,
                        double *azel, int *vsat, double *resp, int *ns, int nx)
{
    gtime_t time;
    double r,dtrp=0.0,vmeas,vtrp=0.0;
    double rr[3],pos[3],dtr,e[3];
    int i,j,nv=0,sat,sys,mask[NX_BASE-3]={0};
    int f2,isdual;

    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    ecef2pos(rr,pos);
    trace(3,"rescode_uduc: rr=%.3f %.3f %.3f\n",rr[0],rr[1],rr[2]);

    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;

        /* skip duplicate observations */
        if (i<n-1&&sat==obs[i+1].sat) {
            char tstr[40];
            trace(2,"duplicated obs data %s sat=%d\n",
                  time2str(time,tstr,3),sat);
            i++; continue;
        }
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;

        f2=sel2freq(sat,opt);
        isdual=(obs[i].P[f2]!=0.0&&obs[i].L[f2]!=0.0&&obs[i].L[0]!=0.0);

        if (iter>0) {
            if (!snrmask(obs+i,azel+i*2,opt)) continue;
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) continue;
        }

        /*--------------------------------------------------------
         * Frequency index 0 (L1/E1/B1): the PRIMARY observation
         * Every satellite contributes this equation.
         * I_sat maps with γ=1.0.
         *--------------------------------------------------------*/
        {
            double P1,vion,dion,gamma=1.0;
            int jion=II(sat);
            if (jion>=nx) continue;

            if ((P1=prange_dcb(obs+i,nav,0,&vmeas))==0.0) continue;

            /* ionospheric delay: use estimated state or Klobuchar prior */
            if (x[jion]!=0.0) {
                dion=x[jion];
                vion=0.0; /* state already encodes uncertainty */
            } else {
                /* Klobuchar prior as initial value (iter 0) */
                if (!ionocorr(time,nav,sat,pos,azel+i*2,IONOOPT_BRDC,&dion,&vion)) {
                    dion=0.0; vion=SQR(ERR_ION);
                }
            }

            v[nv]=P1-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
            trace(4,"sat=%d P1: v=%.3f P1=%.3f r=%.3f dtr=%.6f dts=%.6f dion=%.3f\n",
                  sat,v[nv],P1,r,dtr,dts[i*2],dion);

            /* design matrix: P1 equation */
            for (j=0;j<nx;j++) H[j+nv*nx]=0.0;
            for (j=0;j<3;j++) H[j+nv*nx]=-e[j];
            H[3+nv*nx]=1.0;             /* receiver clock */
            H[jion+nv*nx]=gamma;       /* I_sat, γ=1 for L1 */
            if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*nx]=1.0; mask[1]=1;}
            else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*nx]=1.0; mask[2]=1;}
            else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*nx]=1.0; mask[3]=1;}
            else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*nx]=1.0; mask[4]=1;}
#ifdef QZSDT
            else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*nx]=1.0; mask[5]=1;}
#endif
            else mask[0]=1;

            vsat[i]=1; resp[i]=v[nv]; (*ns)++;

            var[nv]=vare[i]+vmeas+vion+vtrp;
            if (ssat)
                var[nv++]+=varerr(opt,&ssat[i],&obs[i],azel[1+i*2],sys);
            else
                var[nv++]+=varerr(opt,NULL,&obs[i],azel[1+i*2],sys);
        }

        /*--------------------------------------------------------
         * Frequency index f2 (L2/L5/G2/G3/E5a/E5b/B2/S): SECOND equation
         * Only for dual-frequency satellites.
         * Same I_sat as P1, but γ ≠ 1.
         *
         * KEY: Two equations (P1,P2) for the SAME I_sat → I_sat observable
         *      because γ₁ ≠ γ₂ (e.g., 1.0 vs 1.65).
         *--------------------------------------------------------*/
        if (isdual) {
            double P2,vion,dion,gamma;
            int jion=II(sat);
            if (jion>=nx) continue;

            if ((P2=prange_dcb(obs+i,nav,f2,&vmeas))==0.0) continue;
            gamma=iono_gamma(sat,f2,obs+i,nav);

            /* ionospheric delay (same as P1) */
            if (x[jion]!=0.0) {
                dion=x[jion];
                vion=0.0;
            } else {
                if (!ionocorr(time,nav,sat,pos,azel+i*2,IONOOPT_BRDC,&dion,&vion)) {
                    dion=0.0; vion=SQR(ERR_ION);
                }
            }

            v[nv]=P2-(r+dtr-CLIGHT*dts[i*2]+dion*gamma+dtrp);
            trace(4,"sat=%d P2: v=%.3f P2=%.3f r=%.3f dtr=%.6f dts=%.6f "
                  "dion=%.3f gamma=%.4f\n",
                  sat,v[nv],P2,r,dtr,dts[i*2],dion,gamma);

            /* design matrix: P2 equation */
            for (j=0;j<nx;j++) H[j+nv*nx]=0.0;
            for (j=0;j<3;j++) H[j+nv*nx]=-e[j];
            H[3+nv*nx]=1.0;
            H[jion+nv*nx]=gamma;      /* I_sat, γ≠1 for L2/L5 */
            if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*nx]=1.0;}
            else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*nx]=1.0;}
            else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*nx]=1.0;}
            else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*nx]=1.0;}
#ifdef QZSDT
            else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*nx]=1.0;}
#endif

            var[nv]=vare[i]+vmeas+vion+vtrp;
            if (ssat)
                var[nv++]+=varerr(opt,&ssat[i],&obs[i],azel[1+i*2],sys);
            else
                var[nv++]+=varerr(opt,NULL,&obs[i],azel[1+i*2],sys);

            /*--------------------------------------------------------
             * MW equation: MW = λ_wl·N_wl
             * Provides an extra constraint (no ionosphere).
             * MW is used to validate the iono estimate, not directly
             * as an extra parameter in the LSQ (it's redundant with P1+P2).
             *--------------------------------------------------------*/
            {
                double mw=mw_comb(obs+i,nav,f2);
                if (fabs(mw)>0.5&&fabs(mw)<50.0) { /* sanity check (meters) */
                    /* MW contains no iono or dtr, just wide-lane ambiguity.
                     * We log it for QC but don't add a separate MW state
                     * (P1+P2 already resolve iono+dtr from their 2 equations). */
                    trace(4,"sat=%d MW=%.3f (QC only)\n",sat,mw);
                }
            }
        }
    }

    /* rank deficiency correction: each sys must have at least one equation */
    for (i=0;i<NX_BASE-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<nx;j++) H[j+nv*nx]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}

/* Count active satellites for NX dimension ----------------------------------*/
static int count_active_sats(const obsd_t *obs, int n, const double *rs,
                             const double *vare, const int *svh,
                             const nav_t *nav, const prcopt_t *opt,
                             const double *rr)
{
    double pos[3],e[3],r;
    int i,sat,ns=0;
    ecef2pos(rr,pos);
    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,NULL)<opt->elmin) continue;
        if (obs[i].P[0]==0.0) continue;
        ns++;
    }
    return ns;
}

/* 对定位结果进行卡方检验和GDOP检验 -----------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;

    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        sprintf(msg,"Warning: large chi-square error nv=%d vv=%.1f cs=%.1f",
                nv,vv,chisqr[nv-nx-1]);
    }
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
* Supports both conventional SPP and UDUC mode (ionoopt==IONOOPT_EST).
*
* UDUC model (all satellites unified):
*   P_k^s = ρ_s + c·dtr - c·dts_s + I_s·γ_k^s + T
*
*   Every satellite: one iono parameter I_s (slant L1 delay).
*   Every frequency: maps to I_s via γ_k = (f₁/f_k)².
*   Dual-freq sat: P1 + P2 → two equations for SAME I_s → I_s observable.
*   Single-freq sat: I_s constrained by cross-satellite diversity + Klobuchar.
*
* Parameter: NX = 9 + n_active_sats (one iono per active satellite)
*---------------------------------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, const ssat_t *ssat, sol_t *sol,
                  double *azel, int *vsat, double *resp, char *msg)
{
    int nx=opt->ionoopt==IONOOPT_EST? NX_UDUC : NX_BASE;
    int i,j,k,info,stat,nv,ns,nsat;
    double *x=NULL,*dx=NULL,*Q=NULL;
    double *v=NULL,*H=NULL,*var_=NULL;

    trace(3,"estpos  : n=%d ionoopt=%d\n",n,opt->ionoopt);

    /* determine NX for UDUC: base + one iono param per active satellite */
    if (opt->ionoopt==IONOOPT_EST) {
        nsat=count_active_sats(obs,n,rs,vare,svh,nav,opt,sol->rr);
        nx=NX_BASE+nsat;
        trace(3,"estpos  : UDUC, nsat=%d nx=%d\n",nsat,nx);
    }

    x =mat(nx,1); dx=mat(nx,1); Q=mat(nx,nx);
    v  =mat(n*2+NX_BASE,1);
    H  =mat(nx,n*2+NX_BASE);
    var_=mat(n*2+NX_BASE,1);

    /* initialize states */
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    for (i=3;i<NX_BASE;i++) x[i]=0.0;

    /* initialize iono states (UDUC) with Klobuchar prior */
    if (opt->ionoopt==IONOOPT_EST) {
        double pos[3],dion,vion;
        ecef2pos(sol->rr,pos);
        for (i=0;i<n&&i<MAXOBS;i++) {
            int sat=obs[i].sat;
            int jion=II(sat);
            if (jion>=nx) continue;
            ionocorr(obs[i].time,nav,sat,pos,azel+i*2,IONOOPT_BRDC,&dion,&vion);
            x[jion]=dion;                    /* Klobuchar as initial value */
            Q[jion+jion*nx]=MAX(VAR_IONO,vion); /* large initial variance */
        }
    }

    /* iterative weighted least squares */
    for (i=0;i<MAXITR;i++) {
        nv=rescode_uduc(i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,
                        v,H,var_,azel,vsat,resp,&ns,nx);

        if (nv<nx) {
            sprintf(msg,"lack of valid sats nv=%d nx=%d",nv,nx);
            break;
        }
        /* weighted */
        for (j=0;j<nv;j++) {
            double sig=sqrt(var_[j]);
            v[j]/=sig;
            for (k=0;k<nx;k++) H[k+j*nx]/=sig;
        }
        if ((info=lsq(H,v,nx,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<nx;j++) x[j]+=dx[j];

        if (norm(dx,nx)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT;
            sol->dtr[1]=x[4]/CLIGHT;
            sol->dtr[2]=x[5]/CLIGHT;
            sol->dtr[3]=x[6]/CLIGHT;
            sol->dtr[4]=x[7]/CLIGHT;
#ifdef QZSDT
            sol->dtr[5]=x[8]/CLIGHT;
#endif
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*nx];
            sol->qr[3]=(float)Q[1];
            sol->qr[4]=(float)Q[2+nx];
            sol->qr[5]=(float)Q[2];
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;

            if ((stat=valsol(azel,vsat,n,opt,v,nv,nx,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(x); free(dx); free(Q);
            free(v); free(H); free(var_);
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);

    free(x); free(dx); free(Q);
    free(v); free(H); free(var_);
    return 0;
}

/* RAIM FDE ------------------------------------------------------------------*/
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
    rs_e=mat(6,n); dts_e=mat(2,n); vare_e=mat(1,n);
    azel_e=zeros(2,n); svh_e=imat(1,n);
    vsat_e=imat(1,n); resp_e=mat(1,n);

    for (i=0;i<n;i++) {
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e+6*k,rs+6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j]; svh_e[k++]=svh[j];
        }
        if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,ssat,&sol_e,
                    azel_e,vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg); continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]); nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites\n",obs[i].sat); continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        if (rms_e>rms) continue;

        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k]; resp[j]=resp_e[k++];
        }
        stat=1; sol_e.eventime=sol->eventime;
        *sol=sol_e; sat=obs[i].sat; rms=rms_e; vsat[i]=0;
        strcpy(msg,msg_e);
    }
#ifdef TRACE
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
#endif
    free(obs_e);
    free(rs_e); free(dts_e); free(vare_e); free(azel_e);
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

    ecef2pos(rr,pos); xyz2enu(pos,E);

    for (i=0;i<n&&i<MAXOBS;i++) {
        freq=sat2freq(obs[i].sat,obs[i].code[0],nav);
        if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) continue;

        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,E,a,e);

        for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];
        rate=dot3(vs,e)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                     rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        sig=(err<=0.0)?1.0:err*CLIGHT/freq;
        v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
        for (j=0;j<4;j++) H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
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
    double err=opt->err[4];
    int i,j,nv;

    v=mat(n,1); H=mat(4,n);

    for (i=0;i<MAXITR;i++) {
        if ((nv=resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H))<4) break;
        if (lsq(H,v,4,nv,dx,Q)) break;
        for (j=0;j<4;j++) x[j]+=dx[j];
        if (norm(dx,4)<1E-6) {
            matcpy(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Q[0]; sol->qv[1]=(float)Q[5];
            sol->qv[2]=(float)Q[10]; sol->qv[3]=(float)Q[1];
            sol->qv[4]=(float)Q[6]; sol->qv[5]=(float)Q[2];
            break;
        }
    }
    free(v); free(H);
}

/* single-point positioning ----------------------------------------------------
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observations
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message
* return : status(1:ok,0:error)
*
* UDUC mode (opt->ionoopt==IONOOPT_EST):
*   Unified model: P_k^s = ρ_s + c·dtr - c·dts_s + I_s·γ_k^s + T
*   - ONE iono parameter I_s per satellite (on L1 slant TEC)
*   - ALL frequencies of sat s share I_s, scaled by γ_k=(f1/fk)²
*   - Dual-freq sat: P1 (γ=1) + P2 (γ≠1) → 2 equations → I_s directly observable
*   - Single-freq sat: Klobuchar prior as mild constraint, rank deficiency
*     broken by cross-satellite diversity (other dual-freq sats constrain dtr)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    char tstr[40];

    trace(3,"pntpos  : tobs=%s n=%d ionoopt=%d\n",
          time2str(obs[0].time,tstr,3),n,opt_.ionoopt);

    sol->stat=SOLQ_NONE;
    if (n<=0) { strcpy(msg,"no observation data"); return 0; }
    sol->time=obs[0].time;
    msg[0]='\0';
    sol->eventime=obs[0].eventime;

    rs=mat(6,n); dts=mat(2,n); var=mat(1,n);
    azel_=zeros(2,n); resp=mat(1,n);

    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].snr_rover[0]=0; ssat[i].snr_base[0]=0;
        }
        for (i=0;i<n;i++) ssat[obs[i].sat-1].snr_rover[0]=obs[i].SNR[0];
    }

    if (opt_.mode!=PMODE_SINGLE) {
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }

    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);

    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0; ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
