/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Loi_paroi_adaptative.cpp
// Directory:   $TRUST_ROOT/Multiphase/CMFD/src/Turbulence/Lois_paroi
//
//////////////////////////////////////////////////////////////////////////////

#include <Loi_paroi_adaptative.h>
#include <Navier_Stokes_std.h>
#include <Correlation_base.h>
#include <QDM_Multiphase.h>
#include <TRUSTTab_parts.h>
#include <Cond_lim_base.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Motcle.h>
#include <Param.h>
#include <math.h>
#include <Nom.h>
#include <Champ_Face_base.h>

Implemente_instanciable(Loi_paroi_adaptative, "Loi_paroi_adaptative", Loi_paroi_log);

Sortie& Loi_paroi_adaptative::printOn(Sortie& os) const
{
  return os;
}

Entree& Loi_paroi_adaptative::readOn(Entree& is)
{
  return Loi_paroi_log::readOn(is);
}

double Loi_paroi_adaptative::u_plus_de_y_plus(double y_p)  // Blended Reichardt model
{
  double reichardt = std::log(1+0.4*y_p)/von_karman_;
  reichardt += 7.8;
  reichardt += -7.8*std::exp(-y_p/11);
  reichardt += -7.8*y_p/11*std::exp(-y_p/3);

  double log_law = std::log(y_p+limiteur_y_p)/von_karman_ + 5.1;

  double blending = std::tanh( y_p/27*y_p/27*y_p/27*y_p/27);

  return (1-blending)*reichardt + blending*log_law;
}

double Loi_paroi_adaptative::deriv_u_plus_de_y_plus(double y_p)
{
  double reichardt = std::log(1+0.4*y_p)/von_karman_ + 7.8 -7.8*std::exp(-y_p/11) -7.8*y_p/11*std::exp(-y_p/3);
  double log_law = std::log(y_p+limiteur_y_p)/von_karman_ + 5.1;
  double blending = std::tanh( y_p/27*y_p/27*y_p/27*y_p/27);

  double d_reichardt = 0.4/(1+0.4*y_p)*1/von_karman_;
  d_reichardt += 7.8/11*std::exp(-y_p/11);
  d_reichardt += -7.8/11*std::exp(-y_p/3) + 7.8*y_p/33*std::exp(-y_p/3) ;
  double d_log_law = 1/(y_p+limiteur_y_p)*1./von_karman_;
  double d_blending = 4/27*(y_p/27*y_p/27*y_p/27)*(1-blending*blending);

  return (1-blending)*d_reichardt - reichardt*d_blending + blending*d_log_law + d_blending*log_law ;
}

