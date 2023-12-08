/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Paroi_frottante_loi.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/28
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_frottante_loi.h>

#include <Viscosite_turbulente_k_omega.h>
#include <Op_Dift_Multiphase_VDF_Face.h>
#include <Viscosite_turbulente_k_tau.h>
#include <Op_Diff_PolyMAC_P0_base.h>
#include <Loi_paroi_adaptative.h>
#include <Frontiere_dis_base.h>
#include <Navier_Stokes_std.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>
#include <Frontiere.h>
#include <Motcle.h>

#include <math.h>

Implemente_instanciable(Paroi_frottante_loi,"Paroi_frottante_loi", Paroi_frottante_simple);
// XD Paroi_frottante_loi condlim_base Paroi_frottante_loi 1 Adaptive wall-law boundary condition for velocity

Sortie& Paroi_frottante_loi::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Paroi_frottante_loi::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("y_p_prod_k", &y_p_prod_k_);
  param.ajouter("fac_prod_k", &fac_prod_k_);
  param.ajouter("y_p_prod_k_grand", &y_p_prod_k_grand_);
  param.ajouter("fac_prod_k_grand", &fac_prod_k_grand_);
  param.lire_avec_accolades_depuis(s);

  Paroi_frottante_simple::readOn(s);

  return s;
}

void Paroi_frottante_loi::completer()
{
  if (!ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).is_turb()) Process::exit(que_suis_je() + " : diffusion operator must be turbulent !");

  if (is_externe_)
    {
      if (fac_prod_k_<-1.e7) fac_prod_k_ =  0.;
      if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
      if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ =  0.;
      if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 120.;
    }
  else
    {
      if sub_type(Viscosite_turbulente_k_tau, (*ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).correlation_viscosite_turbulente()).valeur())
        {
          if (fac_prod_k_<-1.e7) fac_prod_k_ = 1.2;
          if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
          if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ = .2;
          if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 150.;
        }
      else if sub_type(Viscosite_turbulente_k_omega, (*ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).correlation_viscosite_turbulente()).valeur())
        {
          if (fac_prod_k_<-1.e7) fac_prod_k_ = 1.0;
          if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
          if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ = .6;
          if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 120.;
        }
      else
        {
          if (fac_prod_k_<-1.e7) fac_prod_k_ =  0.;
          if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
          if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ =  0.;
          if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 120.;
        }
    }
}

double Paroi_frottante_loi::fac_coeff_grad(double y_p) const
{
  return 1 + fac_prod_k_ * std::tanh( (y_p/y_p_prod_k_)*(y_p/y_p_prod_k_) ) - fac_prod_k_grand_ * std::tanh( (y_p/y_p_prod_k_grand_)*(y_p/y_p_prod_k_grand_)*(y_p/y_p_prod_k_grand_) );
}