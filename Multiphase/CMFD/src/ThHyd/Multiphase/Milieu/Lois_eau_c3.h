/////////////////////////////////////////////////////////////////////////////
//
// File      : Lois_Eau_C3.h
// Directory : $MATHYS_ROOT/src/ThHyd/Multiphase/Milieu/Lois_Eau_C3.h
//
/////////////////////////////////////////////////////////////////////////////
#include <arch.h>
#define HAVE_LIBC3 0
/* prototypes des fonctions Fortran */
#ifdef HAVE_LIBC3
extern "C" {
#ifdef F77_Majuscule
  int F77DECLARE(FTSATP)(const int *nsca, const int* ienc, const double *pv,
                         double *tspv, double *dtspvv, double *d2tsdpvv,
                         double *hlsv, double *dhlsvv,
                         double *hvsv, double *dhvsvv,
                         double *cplsv, double *dclsvv,
                         double *cpvsv, double *dcvsvv,
                         double *rlsv, double *drlsvv,
                         double *rvsv, double *drvsvv);

  int F77DECLARE(FPSATT)(const int *ienc, const double *tsp0,
                         double *p, double *dtsp1, double *d2tsp1,
                         double *hlsp, double *dhlsp1, double *hvsp, double *dhvsp1,
                         double *cplsp, double *dclsp1, double *cpvsp, double *dcvsp1,
                         double *rlsp, double *drlsp1, double *rvsp, double *drvsp1);

  int F77DECLARE(FHLIQ)(const int *nsca, const int* ienc,
                        const double *p, const double*hl,
                        const double *tspv, const double *dtspvv,
                        const double *hlsv, const double *dhlsvv,
                        const double *cplsv, const double *dclsvv,
                        const double *rlsv, const double *drlsvv,
                        double *tl,double *dtl1,double *dtl2,
                        double *cpl,double *dcpl1,double *dcpl2,
                        double *rl,double *drl1,double *drl2,
                        double *el,double *del1,double *del2,
                        double *betal,double *dbetal1,double *dbetal2);

  int F77DECLARE(FHVAPP)(const int* nsca, const int* itermin, const int* ienc,
                         const double *p, const double*hv, const double* tgini,
                         const double *tspv, const double *dtspvv,
                         const double *hvsv, const double *dhvsvv,
                         const double *cpvsv, const double *dcvsvv,
                         const double *rvsv, const double *drvsvv,
                         double *tv,double *dtv1,double *dtv2,
                         double *cpv,double *dcpv1,double *dcpv2,
                         double *rv,double *drv1,double *drv2,
                         double *ev,double *dev1,double *dev2,
                         double *Tgk, double *hvmhs, int* ill, int* ivstat);

  int F77DECLARE(FHLIQA)(const int *nsca,
                         const double *p, const double*hl,
                         const double *tl, const double *dtl1, const double *dtl2,
                         double *tlal,double *dtlal1,double *dtlal2,
                         double *tmul,double *dtmul1,double *dtmul2);

  int F77DECLARE(FHVAPA)(const int *nsca, const int* ienc,
                         const double *pv, const double* Tg,
                         const double *tspv, const double *dtspvv,
                         double *tlav,double *dlavpv,double *dlavtg,
                         double *tmuv,double *dmuvpv,double *dmuvtg,
                         double* sigma, double* dsigpv);

  int F77DECLARE(FTLIQ)(const int *nsca, const double *p, const double* Tl,
                        double *hl,double *dhldp,double *dhldtl,
                        double *cpl,double *dcpldp,double *dcpldtl,
                        double *rl,double *drldp,double *drldtl);

  int F77DECLARE(FTVAP)(const int *nsca, const int *ienc, int *ier, int *itest,
                        const double *pv, const double *tg,
                        const double *tspv, const double *dtspvv, const double *hvsv, const double *dhvsvv,
                        const double *vapa, const double *vapb, const double *vapc,
                        const double *vapdb, const double *vapdc,
                        double *hv, double *dhv1,
                        double *cpv, double *dcpvpv, double *dcpvtg
                        double *rv, double *drv1, double *drv3, double *hvmhvs);

#else
  int F77DECLARE(ftsatp)(const int *nsca, const int* ienc, const double *pv,
                         double *tspv, double *dtspvv, double *d2tsdpvv,
                         double *hlsv, double *dhlsvv,
                         double *hvsv, double *dhvsvv,
                         double *cplsv, double *dclsvv,
                         double *cpvsv, double *dcvsvv,
                         double *rlsv, double *drlsvv,
                         double *rvsv, double *drvsvv);

  int F77DECLARE(fpsatt)(const int *ienc, const double *tsp0,
                         double *p, double *dtsp1, double *d2tsp1,
                         double *hlsp, double *dhlsp1, double *hvsp, double *dhvsp1,
                         double *cplsp, double *dclsp1, double *cpvsp, double *dcvsp1,
                         double *rlsp, double *drlsp1, double *rvsp, double *drvsp1);

  int F77DECLARE(fhliq)(const int *nsca, const int* ienc,
                        const double *p, const double*hl,
                        const double *tspv, const double *dtspvv,
                        const double *hlsv, const double *dhlsvv,
                        const double *cplsv, const double *dclsvv,
                        const double *rlsv, const double *drlsvv,
                        double *tl,double *dtl1,double *dtl2,
                        double *cpl,double *dcpl1,double *dcpl2,
                        double *rl,double *drl1,double *drl2,
                        double *el,double *del1,double *del2,
                        double *betal,double *dbetal1,double *dbetal2);

  int F77DECLARE(fhvapp)(const int* nsca, const int* itermin, const int* ienc,
                         const double *p, const double*hv, const double* tgini,
                         const double *tspv, const double *dtspvv,
                         const double *hvsv, const double *dhvsvv,
                         const double *cpvsv, const double *dcvsvv,
                         const double *rvsv, const double *drvsvv,
                         double *tv,double *dtv1,double *dtv2,
                         double *cpv,double *dcpv1,double *dcpv2,
                         double *rv,double *drv1,double *drv2,
                         double *ev,double *dev1,double *dev2,
                         double *Tgk, double *hvmhs, int* ill, int* ivstat);

  int F77DECLARE(fhliqa)(const int *nsca,
                         const double *p, const double*hl,
                         const double *tl, const double *dtl1, const double *dtl2,
                         double *tlal,double *dtlal1,double *dtlal2,
                         double *tmul,double *dtmul1,double *dtmul2);

  int F77DECLARE(fhvapa)(const int *nsca, const int* ienc,
                         const double *pv, const double* Tg,
                         const double *tspv, const double *dtspvv,
                         double *tlav,double *dlavpv,double *dlavtg,
                         double *tmuv,double *dmuvpv,double *dmuvtg,
                         double* sigma, double* dsigpv);

  int F77DECLARE(ftliq)(const int *nsca, const double *p, const double* Tl,
                        double *hl,double *dhldp,double *dhldtl,
                        double *cpl,double *dcpldp,double *dcpldtl,
                        double *rl,double *drldp,double *drldtl);

  int F77DECLARE(ftvap)(const int *nsca, const int *ienc, int *ier, int *itest,
                        const double *pv, const double *tg,
                        const double *tspv, const double *dtspvv, const double *hvsv, const double *dhvsvv,
                        const double *vapa, const double *vapb, const double *vapc,
                        const double *vapdb, const double *vapdc,
                        double *hv, double *dhv1,
                        double *cpv, double *dcpvpv, double *dcpvtg,
                        double *rv, double *drv1, double *drv3, double *hvmhvs);

#endif

#ifndef F77_Majuscule
  inline int F77NAME(FTSATP)(const int *nsca, const int* ienc,const double *pv,
                             double *tspv, double *dtspvv, double *d2tsdpvv,
                             double *hlsv, double *dhlsvv,
                             double *hvsv, double *dhvsvv,
                             double *cplsv, double *dclsvv,
                             double *cpvsv, double *dcvsvv,
                             double *rlsv, double *drlsvv,
                             double *rvsv, double *drvsvv)
  {
    return F77NAME(ftsatp)(nsca,ienc,pv,tspv,dtspvv,d2tsdpvv,hlsv,dhlsvv,hvsv,dhvsvv,
                           cplsv,dclsvv,cpvsv,dcvsvv,rlsv,drlsvv,rvsv,drvsvv);
  }

  inline int F77NAME(FPSATT)(const int *ienc, const double *tsp0,
                             double *p, double *dtsp1, double *d2tsp1,
                             double *hlsp, double *dhlsp1, double *hvsp, double *dhvsp1,
                             double *cplsp, double *dclsp1, double *cpvsp, double *dcvsp1,
                             double *rlsp, double *drlsp1, double *rvsp, double *drvsp1)
  {
    return F77NAME(fpsatt)(ienc, tsp0, p, dtsp1, d2tsp1, hlsp, dhlsp1, hvsp, dhvsp1, cplsp, dclsp1, cpvsp, dcvsp1, rlsp, drlsp1, rvsp, drvsp1);
  }

  inline int F77NAME(FHLIQ)(const int *nsca, const int* ienc,
                            const double *p, const double*hl,
                            const double *tspv, const double *dtspvv,
                            const double *hlsv, const double *dhlsvv,
                            const double *cplsv, const double *dclsvv,
                            const double *rlsv, const double *drlsvv,
                            double *tl,double *dtl1,double *dtl2,
                            double *cpl,double *dcpl1,double *dcpl2,
                            double *rl,double *drl1,double *drl2,
                            double *el,double *del1,double *del2,
                            double *betal,double *dbetal1,double *dbetal2)

  {
    return F77NAME(fhliq)(nsca,ienc,
                          p,hl,
                          tspv,dtspvv,
                          hlsv,dhlsvv,
                          cplsv,dclsvv,
                          rlsv,drlsvv,
                          tl,dtl1,dtl2,
                          cpl,dcpl1,dcpl2,
                          rl,drl1,drl2,
                          el,del1,del2,
                          betal,dbetal1,dbetal2);
  }

  inline int F77NAME(FHVAPP)(const int* nsca, const int* itermin, const int* ienc,
                             const double *p, const double*hv, const double* tgini,
                             const double *tspv, const double *dtspvv,
                             const double *hvsv, const double *dhvsvv,
                             const double *cpvsv, const double *dcvsvv,
                             const double *rvsv, const double *drvsvv,
                             double *tv,double *dtv1,double *dtv2,
                             double *cpv,double *dcpv1,double *dcpv2,
                             double *rv,double *drv1,double *drv2,
                             double *ev,double *dev1,double *dev2,
                             double *Tgk, double *hvmhs, int* ill, int* ivstat)
  {
    return F77NAME(fhvapp)(nsca,itermin,ienc,
                           p,hv,tgini,
                           tspv,dtspvv,
                           hvsv,dhvsvv,
                           cpvsv,dcvsvv,
                           rvsv,drvsvv,
                           tv,dtv1,dtv2,
                           cpv,dcpv1,dcpv2,
                           rv,drv1,drv2,
                           ev,dev1,dev2,
                           Tgk,hvmhs,ill,ivstat);
  }

  inline int F77NAME(FHLIQA)(const int *nsca,
                             const double *p, const double*hl,
                             const double *tl, const double *dtl1, const double *dtl2,
                             double *tlal,double *dtlal1,double *dtlal2,
                             double *tmul,double *dtmul1,double *dtmul2)
  {
    return F77NAME(fhliqa)(nsca,
                           p,hl,
                           tl,dtl1,dtl2,
                           tlal,dtlal1,dtlal2,
                           tmul,dtmul1,dtmul2);

  }

  inline int F77NAME(FHVAPA)(const int *nsca, const int* ienc,
                             const double *pv, const double* Tg,
                             const double *tspv, const double *dtspvv,
                             double *tlav,double *dlavpv,double *dlavtg,
                             double *tmuv,double *dmuvpv,double *dmuvtg,
                             double* sigma, double* dsigpv)
  {
    return F77NAME(fhvapa)(nsca,ienc,
                           pv,Tg,
                           tspv,dtspvv,
                           tlav,dlavpv,dlavtg,
                           tmuv,dmuvpv,dmuvtg,
                           sigma,dsigpv);

  }

  inline int F77NAME(FTLIQ)(const int *nsca, const double *p, const double* Tl,
                            double *hl,double *dhldp,double *dhldtl,
                            double *cpl,double *dcpldp,double *dcpldtl,
                            double *rl,double *drldp,double *drldtl)

  {

    return F77NAME(ftliq) (nsca,p,Tl,
                           hl,dhldp,dhldtl,
                           cpl,dcpldp,dcpldtl,
                           rl,drldp,drldtl);


  }

  inline int F77NAME(FTVAP)(const int *nsca, const int *ienc, int *ier, int *itest,
                            const double *pv, const double *tg,
                            const double *tspv, const double *dtspvv, const double *hvsv, const double *dhvsvv,
                            const double *vapa, const double *vapb, const double *vapc,
                            const double *vapdb, const double *vapdc,
                            double *hv, double *dhv1,
                            double *cpv, double *dcpvpv, double *dcpvtg,
                            double *rv, double *drv1, double *drv3, double *hvmhvs)
  {
    return F77NAME(ftvap)(nsca, ienc, ier, itest, pv, tg, tspv, dtspvv, hvsv, dhvsvv, vapa, vapb, vapc, vapdb, vapdc,
                          hv, dhv1, cpv, dcpvpv, dcpvtg, rv, drv1, drv3, hvmhvs);
  }
#endif
}
#endif
