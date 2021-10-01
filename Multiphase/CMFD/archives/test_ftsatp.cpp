
#include <Objet_U.h>
extern "C"
{
#ifdef F77_Majuscule
  entier F77DECLARE(FTSATP)(const entier *nsca, const entier* ienc, const double *pv,
			   double *tspv, double *dtspvv, double *d2tsdpvv,
                            double *hlsv, double *dhlsvv,
			   double *hvsv, double *dhvsvv,
                           double *cplsv, double *dclsvv,
                           double *cpvsv, double *dcvsvv,
                           double *rlsv, double *drlsvv,
                           double *rvsv, double *drvsvv);
#else
  entier F77DECLARE(ftsatp)(const entier *nsca, const entier* ienc, const double *pv,
			   double *tspv, double *dtspvv, double *d2tsdpvv,
			   double *hlsv, double *dhlsvv,
			   double *hvsv, double *dhvsvv,
                           double *cplsv, double *dclsvv,
                           double *cpvsv, double *dcvsvv,
                           double *rlsv, double *drlsvv,
                           double *rvsv, double *drvsvv);
#endif
}

#ifndef F77_Majuscule
inline entier F77NAME(FTSATP)(const entier *nsca, const entier* ienc,const double *pv,
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
#endif


int main()
{
  
  int nsca=1;
  int ienc=0;
  double pv=1e5,tspv,dtspvv,d2tsdpvv,hlsv,dhlsvv,hvsv,dhvsv;
  double cplsv,dclsvv,cpvsv,dcvsvv,rlsv,drlsvv,rvsv,drvsvv;
  F77NAME(FTSATP)(&nsca,&ienc,&pv,&tspv,&dtspvv,&d2tsdpvv,&hlsv,&dhlsvv,&hvsv,&dhvsv,
                  &cplsv,&dclsvv,&cpvsv,&dcvsvv,&rlsv,&drlsvv,&rvsv,&drvsvv);
  
  return 0;
}
