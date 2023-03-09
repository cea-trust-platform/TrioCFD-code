/////////////////////////////////////////////////////////////////////////////
//
// File      : Lois_R12_c1.h
// Directory : $MATHYS_ROOT/src/ThHyd/Multiphase/Milieu/Lois_R12_c1.h
//
/////////////////////////////////////////////////////////////////////////////
#include <arch.h>
#include <Lois_eau_c3.h>
#define HAVE_LIBC3 0
/* prototypes des fonctions Fortran */
#ifdef HAVE_LIBC3
extern "C" {
#ifdef F77_Majuscule
// FPSATR12 Sat temperature from pressure
  int F77DECLARE(FPSATR12)(const int *nsca, const double *p,
                           double *tsp, double *dtsp1,
                           int* ill, int* ivstat, int* ierrth);

// FHSATR12 Sat enthalpy from sat temperature, densities and pressure
  int F77DECLARE(FHSATR12)(const int *nsca, const double *tsp, const double *rvsat, const double *rlsat, const double *p,
                           double *hvsp, double *hlsp,
                           const double *dtsp1, const double *drvsat1, const double *drlsat1,
                           double *dhvsp1, double *dhlsp1);

// FROVLR12 Densities from gas&liquid temperature, pressure
  int F77DECLARE(FROVLR12)(const int *nsca, const double *p, const double *tg, const double *tl,
                           double *rg, double *rl,
                           const double *dtl1,
                           double *drl1, double *drg1,
                           int* ill, int* ivstat, int* ierrth);

// FSIGMAR12 Get surface tension from sat temperature
  int F77DECLARE(FSIGMAR12)(const int *nsca, const double* tsp, const double* dtsp1,
                            double* sigma, double* dsig1);

// FCONLR12 Liquid phase conductivity from liquid temperature
  int F77DECLARE(FCONLR12)(const int *nsca, const double *tl, const double *dtl1, const double *dtl2,
                           double *lambl, double *dlambl1, double *dlambl2);

// FCONVR12 Gas phase conductivity from liquid temperature
  int F77DECLARE(FCONVR12)(const int *nsca, const double *tg, const double *dtg1, const double *dtg3,
                           double *lambv, double *dlambv1, double *dlambv3);

// FCONLR12 Liquid phase viscosity from liquid temperature
  int F77DECLARE(FMULR12)(const int *nsca, const double *tl, const double *dtl1, const double *dtl2,
                          double *tmul, double *dtmul1, double *dtmul2);

// FCONVR12 Gas phase viscosity from liquid temperature
  int F77DECLARE(FMUVR12)(const int *nsca, const double *p, const double *tg, const double *dtg1, const double *dtg3,
                          double *tmug, double *dtmug1, double *dtmug3);

// FCPLR12 Liquid phase Cp & numerous derivatives from tl
  int F77DECLARE(FCPLR12)(const int *nsca, const double *tl,
                          double * cpl, double* dtl1, double* dtl2,
                          double* drl1, double* drl2, double* dcpl1, double*dcpl2,
                          int* ill, int* ivstat, int* ierrth);

// FCPLR12 Gas phase Cp & derivatives from tg and rhog
  int F77DECLARE(FCPVR12)(const int *nsca, const double *tg, const double *dtg1, const double *dtg3,
                          const double* rg, const double* drg1, const double* drg3,
                          double* cpv, double* dcpv1, double* dcpv3);

// FCPLR12 Gas phase enthalpy from pressure and gas temperature
  int F77DECLARE(FPTHGR12)(const int *nsca, const double *p, const double *tg,
                           double* hg,
                           int* ill, int* ivstat, int* ierrth);

// FCPLR12 Liquid phase enthalpy from liquid temperature
  int F77DECLARE(FPTHLR12)(const int *nsca, const double *tl,
                           double* hl,
                           int* ill, int* ivstat, int* ierrth);

// FCPLR12 Gas phase temperature and derivatives from enthalpy and pressure
  int F77DECLARE(FTGR12)(const int *nsca, const double *p, const double *hg, const double *tsp,  const double *rvsat,
                         double* tg, double* dtg1, double* dtg3,
                         double* rg, double* drg1, double* drg3,
                         int* ill, int* ivstat, int* ierrth);

// FCPLR12 Liquid phase temperature from enthalpy and pressure
  int F77DECLARE(FTLR12)(const int *nsca, const double *p, const double *tsp, const double *rvsat,   const double *hvsat,  const double *hl,
                         double* tl, double* rl,
                         int* ill, int* ivstat, int* ierrth);

#else
// FPSATR12 Sat temperature from pressure
  int F77DECLARE(fpsatr12)(const int *nsca, const double *p,
                           double *tsp, double *dtsp1,
                           int* ill, int* ivstat, int* ierrth);

// FHSATR12 Sat enthalpy from sat temperature, densities and pressure
  int F77DECLARE(fhsatr12)(const int *nsca, const double *tsp, const double *rvsat, const double *rlsat, const double *p,
                           double *hvsp, double *hlsp,
                           const double *dtsp1, const double *drvsat1, const double *drlsat1,
                           double *dhvsp1, double *dhlsp1);

// FROVLR12 Densities from gas&liquid temperature, pressure
  int F77DECLARE(frovlr12)(const int *nsca, const double *p, const double *tg, const double *tl,
                           double *rg, double *rl,
                           const double *dtl1,
                           double *drl1, double *drg1,
                           int* ill, int* ivstat, int* ierrth);

// FSIGMAR12 Get surface tension from sat temperature
  int F77DECLARE(fsigmar12)(const int *nsca, const double* tsp, const double* dtsp1,
                            double* sigma, double* dsig1);

// FCONLR12 Liquid phase conductivity from liquid temperature
  int F77DECLARE(fconlr12)(const int *nsca, const double *tl, const double *dtl1, const double *dtl2,
                           double *lambl, double *dlambl1, double *dlambl2);

// FCONVR12 Gas phase conductivity from liquid temperature
  int F77DECLARE(fconvr12)(const int *nsca, const double *tg, const double *dtg1, const double *dtg3,
                           double *lambv, double *dlambv1, double *dlambv3);

// FCONLR12 Liquid phase viscosity from liquid temperature
  int F77DECLARE(fmulr12)(const int *nsca, const double *tl, const double *dtl1, const double *dtl2,
                          double *tmul, double *dtmul1, double *dtmul2);

// FCONVR12 Gas phase viscosity from liquid temperature
  int F77DECLARE(fmuvr12)(const int *nsca, const double *p, const double *tg, const double *dtg1, const double *dtg3,
                          double *tmug, double *dtmug1, double *dtmug3);

// FCPLR12 Liquid phase Cp & numerous derivatives from tl
  int F77DECLARE(fcplr12)(const int *nsca, const double *tl,
                          double * cpl, double* dtl1, double* dtl2,
                          double* drl1, double* drl2, double* dcpl1, double*dcpl2,
                          int* ill, int* ivstat, int* ierrth);

// FCPLR12 Gas phase Cp & derivatives from tg and rhog
  int F77DECLARE(fcpvr12)(const int *nsca, const double *tg, const double *dtg1, const double *dtg3,
                          const double* rg, const double* drg1, const double* drg3,
                          double* cpv, double* dcpv1, double* dcpv3);

// FCPLR12 Gas phase enthalpy from pressure and gas temperature
  int F77DECLARE(fpthgr12)(const int *nsca, const double *p, const double *tg,
                           double* hg,
                           int* ill, int* ivstat, int* ierrth);

// FCPLR12 Liquid phase enthalpy from liquid temperature
  int F77DECLARE(fpthlr12)(const int *nsca, const double *tl,
                           double* hl,
                           int* ill, int* ivstat, int* ierrth);

// FCPLR12 Gas phase temperature and derivatives from enthalpy and pressure
  int F77DECLARE(ftgr12)(const int *nsca, const double *p, const double *hg, const double *tsp,  const double *rvsat,
                         double* tg, double* dtg1, double* dtg3,
                         double* rg, double* drg1, double* drg3,
                         int* ill, int* ivstat, int* ierrth);

// FCPLR12 Liquid phase temperature from enthalpy and pressure
  int F77DECLARE(ftlr12)(const int *nsca, const double *p, const double *tsp, const double *rvsat, const double *hvsat,  const double *hl,
                         double* tl, double* rl,
                         int* ill, int* ivstat, int* ierrth);


#endif

#ifndef F77_Majuscule
  inline int F77NAME(FPSATR12)(const int *nsca, const double *p,
                               double *tsp, double *dtsp1,
                               int* ill, int* ivstat, int* ierrth)
  {
    return F77NAME(fpsatr12)(nsca, p,tsp, dtsp1,ill, ivstat, ierrth);
  }

  inline int F77NAME(FHSATR12)(const int *nsca, const double *tsp, const double *rvsat, const double *rlsat, const double *p,
                               double *hvsp, double *hlsp,
                               const double *dtsp1, const double *drvsat1, const double *drlsat1,
                               double *dhvsp1, double *dhlsp1)
  {
    return F77NAME(fhsatr12)(nsca, tsp, rvsat, rlsat, p,hvsp, hlsp,dtsp1, drvsat1, drlsat1,dhvsp1, dhlsp1);
  }

  inline int F77NAME(FROVLR12)(const int *nsca, const double *p, const double *tg, const double *tl,
                               double *rg, double *rl,
                               const double *dtl1,
                               double *drl1, double *drg1,
                               int* ill, int* ivstat, int* ierrth)
  {
    return F77NAME(frovlr12)(nsca, p, tg, tl,rg, rl,dtl1,drl1, drg1,ill, ivstat, ierrth);
  }

  inline int F77NAME(FSIGMAR12)(const int *nsca, const double* tsp, const double* dtsp1,
                                double* sigma, double* dsig1)
  {
    return F77NAME(fsigmar12)( nsca, tsp, dtsp1,sigma, dsig1);
  }

  inline int F77NAME(FCONLR12)(const int *nsca, const double *tl, const double *dtl1, const double *dtl2,
                               double *lambl, double *dlambl1, double *dlambl2)
  {
    return F77NAME(fconlr12)(nsca,  tl,  dtl1,  dtl2,lambl, dlambl1, dlambl2);
  }

  inline int F77NAME(FCONVR12)(const int *nsca, const double *tg, const double *dtg1, const double *dtg3,
                               double *lambv, double *dlambv1, double *dlambv3)
  {
    return F77NAME(fconvr12)( nsca,  tg,  dtg1,  dtg3,lambv, dlambv1, dlambv3);
  }

  inline int F77NAME(FMULR12)(const int *nsca, const double *tl, const double *dtl1, const double *dtl2,
                              double *tmul, double *dtmul1, double *dtmul2)
  {
    return F77NAME(fmulr12)( nsca,  tl,  dtl1,  dtl2,tmul, dtmul1, dtmul2);
  }

  inline int F77NAME(FMUVR12)(const int *nsca, const double *p, const double *tg, const double *dtg1, const double *dtg3,
                              double *tmug, double *dtmug1, double *dtmug3)
  {
    return F77NAME(fmuvr12)( nsca, p, tg,  dtg1,  dtg3,tmug, dtmug1, dtmug3);
  }

  inline int F77NAME(FCPLR12)(const int *nsca, const double *tl,
                              double * cpl, double* dtl1, double* dtl2,
                              double* drl1, double* drl2, double* dcpl1, double*dcpl2,
                              int* ill, int* ivstat, int* ierrth)
  {
    return F77NAME(fcplr12)( nsca,  tl, cpl, dtl1, dtl2,drl1, drl2, dcpl1, dcpl2,ill,  ivstat,  ierrth);
  }

  inline int F77NAME(FCPVR12)(const int *nsca, const double *tg, const double *dtg1, const double *dtg3,
                              const double* rg, const double* drg1, const double* drg3,
                              double* cpv, double* dcpv1, double* dcpv3)
  {
    return F77NAME(fcpvr12)( nsca,  tg,  dtg1,  dtg3, rg,  drg1,  drg3,cpv,  dcpv1,  dcpv3);
  }

  inline int F77NAME(FPTHGR12)(const int *nsca, const double *p, const double *tg,
                               double* hg,
                               int* ill, int* ivstat, int* ierrth)
  {
    return F77NAME(fpthgr12)( nsca,  p,  tg,  hg,  ill,  ivstat,  ierrth);
  }

  inline int F77NAME(FPTHLR12)(const int *nsca, const double *tl,
                               double* hl,
                               int* ill, int* ivstat, int* ierrth)
  {
    return F77NAME(fpthlr12)( nsca,  tl, hl, ill, ivstat, ierrth);
  }

  inline int F77NAME(FTGR12)(const int *nsca, const double *p, const double *hg, const double *tsp,  const double *rvsat,
                             double* tg, double* dtg1, double* dtg3,
                             double* rg, double* drg1, double* drg3,
                             int* ill, int* ivstat, int* ierrth)
  {
    return F77NAME(ftgr12)(nsca,  p,  hg,  tsp,   rvsat, tg, dtg1, dtg3,rg, drg1,  drg3,ill, ivstat, ierrth);
  }

  inline int F77NAME(FTLR12)(const int *nsca, const double *p, const double *tsp, const double *rvsat, const double *hvsat,  const double *hl,
                             double* tl, double* rl,
                             int* ill, int* ivstat, int* ierrth)
  {
    return F77NAME(ftlr12)( nsca, p, tsp, rvsat, hvsat, hl,tl,  rl,  ill,  ivstat,  ierrth);
  }
#endif
}
#endif
