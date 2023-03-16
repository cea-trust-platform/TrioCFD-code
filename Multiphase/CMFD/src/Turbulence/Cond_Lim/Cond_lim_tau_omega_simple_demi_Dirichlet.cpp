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
// File:        Cond_lim_tau_omega_simple_demi_Dirichlet.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/28
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_lim_tau_omega_simple_demi_Dirichlet.h>

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Op_Diff_PolyMAC_base.h>
#include <Loi_paroi_adaptative.h>
#include <Champ_Face_base.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Motcle.h>

#include <math.h>

Implemente_instanciable(Cond_lim_tau_omega_simple_demi_Dirichlet,"Cond_lim_tau_simple_demi_Dirichlet|Cond_lim_omega_simple_demi_Dirichlet",Dirichlet_loi_paroi);

Sortie& Cond_lim_tau_omega_simple_demi_Dirichlet::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Cond_lim_tau_omega_simple_demi_Dirichlet::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.ajouter("beta_k", &beta_k);
  param.ajouter("von_karman", &von_karman_);
  param.ajouter("facteur_paroi", &facteur_paroi_);
  param.lire_avec_accolades_depuis(s);

  le_champ_front.typer("Champ_front_vide");

  return s;
}

void Cond_lim_tau_omega_simple_demi_Dirichlet::completer()
{
  if (sub_type(Echelle_temporelle_turbulente, domaine_Cl_dis().equation())) is_tau_ = 1;
  else if (sub_type(Taux_dissipation_turbulent, domaine_Cl_dis().equation())) is_tau_ = 0;
  else Process::exit(que_suis_je() + " : equation must be tau/omega !");

  int N = domaine_Cl_dis().equation().inconnue().valeurs().line_size();
  if (N > 1)  Process::exit(que_suis_je() + " : Only one phase for turbulent wall law is coded for now");
}

void Cond_lim_tau_omega_simple_demi_Dirichlet::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_Cl_dis().equation().domaine_dis().valeur());
  const DoubleTab& y = corr_loi_paroi.get_tab("y");
  const DoubleTab& nu_visc = ref_cast(Navier_Stokes_std, domaine_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().passe(),
                   &vit = domaine_Cl_dis().equation().probleme().get_champ("vitesse").passe();


  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), nf_tot = domaine.nb_faces_tot(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face(), D = dimension;
  const IntTab& f_e = domaine.face_voisins();
  const DoubleTab& n_f = domaine.face_normales();
  const DoubleVect& fs = domaine.face_surfaces();

  int Nv = vit.line_size();

  DoubleTab pvit_elem(0, Nv * dimension);
  if (nf_tot == vit.dimension_tot(0))
    {
      const Champ_Face_base& ch = ref_cast(Champ_Face_base, domaine_Cl_dis().equation().probleme().equation(0).inconnue().valeur());
      domaine.domaine().creer_tableau_elements(pvit_elem);
      ch.get_elem_vector_field(pvit_elem, true);
    }

  int n = 0 ; // Carrying phase is 0 for turbulent flows

  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      int e_domaine = (f_e(f_domaine,0)>=0) ? f_e(f_domaine,0) : f_e(f_domaine,1) ; // Make orientation vdf-proof

      double u_orth = 0 ;
      DoubleTrav u_parallel(D);
      if (nf_tot == vit.dimension_tot(0)) // <=> VDF
        {
          for (int d = 0; d <D ; d++) u_orth -= pvit_elem(e_domaine, Nv*d+n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          for (int d = 0 ; d < D ; d++) u_parallel(d) = pvit_elem(e_domaine, Nv*d+n) - u_orth*(-n_f(f_domaine,d))/fs(f_domaine) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
        }
      else // <=> PolyMAC
        {
          for (int d = 0; d <D ; d++) u_orth -= vit(nf_tot + e_domaine * D+d, n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          for (int d = 0 ; d < D ; d++) u_parallel(d) = vit(nf_tot + e_domaine * D + d, n) - u_orth*(-n_f(f_domaine,d))/fs(f_domaine) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
        }
      double norm_u_parallel = std::sqrt(domaine.dot(&u_parallel(0), &u_parallel(0)));

      if (is_tau_ == 1)
        {
          double u_tau_demi = corr_loi_paroi.calc_u_tau_loc(norm_u_parallel, nu_visc(e_domaine, 0), y(f_domaine, 0)/2.);
          d_(f, n) = .5*facteur_paroi_*calc_tau(y(f_domaine, n)/2., u_tau_demi, nu_visc(e_domaine, n));
        }
      if (is_tau_ == 0)
        {
          double u_tau_demi = corr_loi_paroi.calc_u_tau_loc(norm_u_parallel, nu_visc(e_domaine, 0), y(f_domaine, 0)/2.);
          d_(f, n) = 2. * facteur_paroi_*calc_omega(y(f_domaine, n)/2., u_tau_demi, nu_visc(e_domaine, n));
        }
      d_.echange_espace_virtuel();
    }

}

double Cond_lim_tau_omega_simple_demi_Dirichlet::calc_tau(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double w_vis = 6 * visc / (beta_omega * y * y);
  double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
  double w_1 = w_vis + w_log ;
  double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
  double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

  return 1 / (blending * w_1 + (1-blending) * w_2) ;

}

double Cond_lim_tau_omega_simple_demi_Dirichlet::calc_omega(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double w_vis = 6 * visc / (beta_omega * y * y);
  double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
  double w_1 = w_vis + w_log ;
  double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
  double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

  return blending * w_1 + (1-blending) * w_2 ;
}

