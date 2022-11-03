/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        FTd_tools.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#include <FTd_tools.h>
#include <math.h>
#include <TRUSTTabFT.h>

/*! @brief Cette fonction calcule le max de la norme d'un ensemble de vecteurs
 *
 * @param (vecteurs) vecteurs a mesurer
 * @return (double) le max du deplacement
 */
double FTd_calculer_max_norme_vecteurs(const DoubleTab& vecteurs)
{
  int som,k;
  double norme2max = 0., norme2;
  int nbsom = vecteurs.dimension(0);
  for (som=0 ; som<nbsom ; som++)
    {
      norme2 = vecteurs(som,0) * vecteurs(som,0);
      for (k=1 ; k<Objet_U::dimension ; k++)
        {
          norme2 += vecteurs(som,k) * vecteurs(som,k);
        }
      if (norme2>norme2max)
        {
          norme2max = norme2;
        }
    }

  return sqrt(norme2max);
}

// The area of a triangle is half of the cross-product!
double FTd_calculer_aire_triangle(const FTd_vecteur2& coord_som0,
                                  const FTd_vecteur2& coord_som1,
                                  const FTd_vecteur2& coord_som2)
{
  double x0 = coord_som0[0];
  double y0 = coord_som0[1];

  double x1 = coord_som1[0];
  double y1 = coord_som1[1];

  double x2 = coord_som2[0];
  double y2 = coord_som2[1];

  double aire_tr = 0.5*((x1-x0) * (y2-y0) - (y1-y0) * (x2-x0));

  return aire_tr;
}
double FTd_calculer_volume_tetraedre(const FTd_vecteur3& coord_som0,
                                     const FTd_vecteur3& coord_som1,
                                     const FTd_vecteur3& coord_som2,
                                     const FTd_vecteur3& coord_som3)
{
  double x0 = coord_som0[0];
  double y0 = coord_som0[1];
  double z0 = coord_som0[2];

  double x1 = coord_som1[0];
  double y1 = coord_som1[1];
  double z1 = coord_som1[2];

  double x2 = coord_som2[0];
  double y2 = coord_som2[1];
  double z2 = coord_som2[2];

  double x3 = coord_som3[0];
  double y3 = coord_som3[1];
  double z3 = coord_som3[2];

  double vol = (  (x1-x0)*( (y2-y0)*(z3-z0)-(y3-y0)*(z2-z0) )
                  - (x2-x0)*( (y1-y0)*(z3-z0)-(y3-y0)*(z1-z0) )
                  + (x3-x0)*( (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0) ) )/6;

  return vol;
}

//en numerotation globale, cette fonction renvoie
// -1 si som0 < som1
//  0 si som0 = som1
//  1 si som0 > som1
//En numerotation globale,
// som0=som1 si
//  - peO==pe1 && som0==som1
// som0<som1
//  - pe0<pe1 ou (pe0==pe1 && som0<som1)
// som0>som1
//  - pe0>pe1 ou (pe0==pe1 && som0>som1)
int FTd_compare_sommets_global(int pe0, int numOwner0, int pe1, int numOwner1)
{
  if (pe0==pe1)
    {
      if (numOwner0==numOwner1)
        {
          return 0;
        }
      else if (numOwner0<numOwner1)
        {
          return -1;
        }
    }
  else if (pe0<pe1)
    {
      return -1;
    }

  return 1;
}
