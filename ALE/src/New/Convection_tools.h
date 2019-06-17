#include <ArrOfInt.h>
#include <Double.h>

//////////////////////////////////////////////////////////////
//   Fonctions limiteurs de MUSCL
////////////////////////////////////////////////////////////////


double minmod(double grad1, double grad2)
{
  double gradlim=0.;
  if(grad1*grad2>0.) (dabs(grad1)<dabs(grad2)) ? gradlim=grad1 : gradlim=grad2 ;
  return gradlim;
}

double vanleer(double grad1, double grad2)
{
  double gradlim=0.;
  if(grad1*grad2>0.) gradlim=2.*grad1*grad2/(grad1+grad2) ;
  return gradlim;
}

double vanalbada(double grad1, double grad2)
{
  double gradlim=0.;
  if(grad1*grad2>0.) gradlim=grad1*grad2*(grad1+grad2)/(grad1*grad1+grad2*grad2) ;
  return gradlim;
}


double chakravarthy(double grad1, double grad2)
{
  /*
    Cerr << " limiteur chakavarthy non preconise (non symetrique) " << finl;
    exit();
    return 0;
  */
  double gradlim=0.;
  if ((grad1*grad2)>0)
    {
      gradlim=dmin(grad1/grad2,1.8); // 1<<beta<<2
      gradlim=dmax(gradlim,0.);
      gradlim*=grad2;
    }
  return gradlim;
}

double superbee(double grad1, double grad2)
{
  /*
    Cerr << " limiteur superbee non preconise (source d'instabilites) " << finl;
    exit();
    return 0;
  */
  double gradlim=0.;
  if ((grad1*grad2)>0)
    {
      double gradlim1,gradlim2;
      gradlim1=dmin(2*(grad1/grad2),1);
      gradlim2=dmin(grad1/grad2,2);
      gradlim=dmax(gradlim1,gradlim2);
      gradlim=dmax(gradlim,0.);
      gradlim*=grad2;
    }
  return gradlim;
}
