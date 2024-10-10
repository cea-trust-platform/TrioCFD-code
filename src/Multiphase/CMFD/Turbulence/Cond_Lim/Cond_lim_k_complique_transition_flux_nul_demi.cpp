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
// File:        Cond_lim_k_complique_transition_flux_nul_demi.cpp
// Directory:   $TRUST_ROOT/src/ThSol
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_lim_k_complique_transition_flux_nul_demi.h>
#include <Energie_cinetique_turbulente.h>
#include <Op_Dift_Multiphase_VDF_Elem.h>
#include <Viscosite_turbulente_base.h>
#include <Transport_turbulent_base.h>
#include <Op_Diff_PolyMAC_P0_base.h>
#include <Loi_paroi_adaptative.h>
#include <Convection_diffusion_turbulence_multiphase.h>
#include <Champ_Face_base.h>
#include <Probleme_base.h>
#include <Domaine_VF.h>

Implemente_instanciable(Cond_lim_k_complique_transition_flux_nul_demi,"Cond_lim_k_complique_transition_flux_nul_demi",Echange_global_impose_turbulent);
// XD Cond_lim_k_complique_transition_flux_nul_demi condlim_base Cond_lim_k_complique_transition_flux_nul_demi 0 Adaptive wall law boundary condition for turbulent kinetic energy

Sortie& Cond_lim_k_complique_transition_flux_nul_demi::printOn(Sortie& s ) const {return Echange_global_impose_turbulent::printOn(s);}

Entree& Cond_lim_k_complique_transition_flux_nul_demi::readOn(Entree& s ) {return Echange_global_impose_turbulent::readOn(s);}

void Cond_lim_k_complique_transition_flux_nul_demi::completer()
{
  if (!sub_type(Energie_cinetique_turbulente, domaine_Cl_dis().equation())) Process::exit("Cond_lim_k_simple : equation must be k !");
  if (domaine_Cl_dis().equation().inconnue().valeurs().line_size() != 1)  Process::exit("Cond_lim_k_simple : Only one phase for turbulent wall law is coded for now");
}

void Cond_lim_k_complique_transition_flux_nul_demi::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_Cl_dis().equation().domaine_dis());
  const DoubleTab&       yp = corr_loi_paroi.get_tab("y_plus"), &u_tau = corr_loi_paroi.get_tab("u_tau");
  const DoubleTab& nu_visc = ref_cast(Convection_diffusion_turbulence_multiphase, domaine_Cl_dis().equation()).diffusivite_pour_pas_de_temps().passe(),
                   &mu_visc = ref_cast(Convection_diffusion_turbulence_multiphase, domaine_Cl_dis().equation()).diffusivite_pour_transport().passe();

  // On va chercher le mu turbulent de polymac et celui de vdf et on prend le bon dans la suite
  const DoubleTab* mu_poly = domaine.que_suis_je().debute_par("Domaine_PolyMAC") ? &ref_cast(Op_Diff_PolyMAC_P0_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu() : nullptr,
                   *mu_vdf = domaine.que_suis_je().debute_par("Domaine_VDF") ? &ref_cast(Op_Dift_Multiphase_VDF_Elem, domaine_Cl_dis().equation().operateur(0).l_op_base()).get_diffusivite_turbulente() : nullptr;
  assert((mu_poly) || (mu_vdf));

  int nf = la_frontiere_dis->frontiere().nb_faces(), f1 = la_frontiere_dis->frontiere().num_premiere_face();
  const IntTab& f_e = domaine.face_voisins();

  int n = 0 ; // Carrying phase is 0 for turbulent flows

  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      int e_domaine = (f_e(f_domaine,0)>=0) ? f_e(f_domaine,0) : f_e(f_domaine,1) ; // Make orientation vdf-proof
      double y_loc = f_e(f_domaine,0)>=0 ? domaine.dist_face_elem0(f_domaine,e_domaine) : domaine.dist_face_elem1(f_domaine,e_domaine) ;
      double mu_tot_loc = (mu_poly) ? (*mu_poly)(e_domaine,n) : (mu_vdf) ? (*mu_vdf)(e_domaine,n) + mu_visc(e_domaine,n) : -1;

      h_(f, 0) = 2.*mu_tot_loc/y_loc * ( 1 - std::tanh(  std::pow( yp(f_domaine, 0)/50.,3)  ) );
      h_grad_(f, 0) = 2./y_loc * ( 1 - std::tanh(  std::pow( yp(f_domaine, 0)/50.,3)  ) );
      T_(f, 0) = calc_k(y_loc/2., u_tau(f_domaine, 0), nu_visc(e_domaine, 0));
    }

  h_.echange_espace_virtuel();
  h_grad_.echange_espace_virtuel();
  T_.echange_espace_virtuel();
}

double Cond_lim_k_complique_transition_flux_nul_demi::calc_k(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;

  double  f1 = (y_p-1) *(y_p-1)*(y_p-1)/30 ;

  double  f2 = 1/std::sqrt(beta_k_)-.08* std::pow(     std::abs(4.6 - std::log(y_p))    , 3 ) ;

  double  b1 = std::tanh(  std::pow(y_p/4.,10)  ) ;

  double  b2 = std::tanh(  std::pow(y_p/2500.,1.4)  ) ;

  double  f3 = 1/std::sqrt(beta_k_)*.25 ;

  return u_tau*u_tau*(std::max((1-b1)*f1 + b1*f2, 0.)*(1-b2)+b2*f3) ;
}

