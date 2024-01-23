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
// File:        Cond_lim_k_simple_flux_nul.cpp
// Directory:   $TRUST_ROOT/src/ThSol
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_lim_k_simple_flux_nul.h>

#include <Energie_cinetique_turbulente.h>
#include <Op_Dift_Multiphase_VDF_Elem.h>
#include <Op_Diff_PolyMAC_P0_base.h>
#include <Loi_paroi_adaptative.h>
#include <Domaine_VF.h>

Implemente_instanciable(Cond_lim_k_simple_flux_nul,"Cond_lim_k_simple_flux_nul",Echange_global_impose_turbulent);
// XD Cond_lim_k_simple_flux_nul condlim_base Cond_lim_k_simple_flux_nul 0 Adaptive wall law boundary condition for turbulent kinetic energy


Sortie& Cond_lim_k_simple_flux_nul::printOn(Sortie& s ) const {return Echange_global_impose_turbulent::printOn(s);}

Entree& Cond_lim_k_simple_flux_nul::readOn(Entree& s ) {return Echange_global_impose_turbulent::readOn(s);}

void Cond_lim_k_simple_flux_nul::completer()
{
  if (!sub_type(Energie_cinetique_turbulente, domaine_Cl_dis().equation())) Process::exit("Cond_lim_k_simple : equation must be k !");
  if (domaine_Cl_dis().equation().inconnue().valeurs().line_size() != 1)  Process::exit("Cond_lim_k_simple : Only one phase for turbulent wall law is coded for now");
}

int Cond_lim_k_simple_flux_nul::initialiser(double temps)
{
  int init = Echange_global_impose_turbulent::initialiser(temps);

  int nf = la_frontiere_dis->frontiere().nb_faces();
  for (int f =0 ; f < nf ; f++) T_(f, 0) = 0 ; // K is 0 on the wall

  return init;
}

void Cond_lim_k_simple_flux_nul::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_Cl_dis().equation().domaine_dis().valeur());
  const DoubleTab&   u_tau = corr_loi_paroi.get_tab("u_tau");
  const DoubleTab&  visc_c = ref_cast(Convection_diffusion_turbulence_multiphase, domaine_Cl_dis().equation()).diffusivite_pour_pas_de_temps().passe(),
                    &mu_visc  = ref_cast(Convection_diffusion_turbulence_multiphase, domaine_Cl_dis().equation()).diffusivite_pour_transport().passe();

  const int cvisc = visc_c.dimension(0) == 1, cmu = mu_visc.dimension(0) == 1;
  // On va chercher le mu turbulent de polymac et celui de vdf et on prend le bon dans la suite
  const DoubleTab* mu_poly = domaine.que_suis_je().debute_par("Domaine_PolyMAC") ? &ref_cast(Op_Diff_PolyMAC_P0_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu() : nullptr,
                   *mu_vdf = domaine.que_suis_je().debute_par("Domaine_VDF") ? &ref_cast(Op_Dift_Multiphase_VDF_Elem, domaine_Cl_dis().equation().operateur(0).l_op_base()).get_diffusivite_turbulente() : nullptr;
  assert((mu_poly) || (mu_vdf));

  int nf = la_frontiere_dis->frontiere().nb_faces(), f1 = la_frontiere_dis->frontiere().num_premiere_face();
  const IntTab& f_e = domaine.face_voisins();

  int n=0; // Turbulence only in first phase for now

  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      int e_domaine = (f_e(f_domaine,0)>=0) ? f_e(f_domaine,0) : f_e(f_domaine,1) ; // Make orientation vdf-proof
      double mu_tot_loc = (mu_poly) ? (*mu_poly)(e_domaine,n) : (mu_vdf) ? (*mu_vdf)(e_domaine,n) + mu_visc(!cmu * e_domaine,n) : -1;
      double y_loc = f_e(f_domaine,0)>=0 ? domaine.dist_face_elem0(f_domaine,e_domaine) : domaine.dist_face_elem1(f_domaine,e_domaine) ;

      h_(f, 0) = mu_tot_loc / y_loc * (1-std::tanh( std::pow(y_loc*u_tau(f_domaine, 0)/visc_c(!cvisc * e_domaine, 0)/10.,2))); // Coeff d'echange de mu/y ; /20 avant modif
      h_grad_(f, 0) = 1. / y_loc * (1-std::tanh( std::pow(y_loc*u_tau(f_domaine, 0)/visc_c(!cvisc * e_domaine, 0)/10.,2)));
    }

  h_.echange_espace_virtuel();
  h_grad_.echange_espace_virtuel();
}
