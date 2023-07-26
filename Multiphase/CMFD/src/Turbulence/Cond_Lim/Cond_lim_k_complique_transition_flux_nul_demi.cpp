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
#include <Loi_paroi_adaptative.h>
#include <Frontiere_dis_base.h>
#include <Pb_Multiphase.h>
#include <Domaine_Poly_base.h>
#include <Op_Diff_PolyMAC_base.h>
#include <Op_Diff_PolyMAC_P0_base.h>

Implemente_instanciable(Cond_lim_k_complique_transition_flux_nul_demi,"Cond_lim_k_complique_transition_flux_nul_demi",Echange_global_impose);


Sortie& Cond_lim_k_complique_transition_flux_nul_demi::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Cond_lim_k_complique_transition_flux_nul_demi::readOn(Entree& s )
{
  h_imp_.typer("Champ_front_vide");
  le_champ_front.typer("Champ_front_vide");

  return s;
}

void Cond_lim_k_complique_transition_flux_nul_demi::completer()
{
  if (!sub_type(Energie_cinetique_turbulente, domaine_Cl_dis().equation())) Process::exit("Cond_lim_k_simple : equation must be k !");
  if (domaine_Cl_dis().equation().inconnue().valeurs().line_size() != 1)  Process::exit("Cond_lim_k_simple : Only one phase for turbulent wall law is coded for now");
}

void Cond_lim_k_complique_transition_flux_nul_demi::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis->frontiere().nb_faces(), f1 = la_frontiere_dis->frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++)
    for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

int Cond_lim_k_complique_transition_flux_nul_demi::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Turbulence="Turbulence";

  if (dom_app==Turbulence)
    return 1;
  else err_pas_compatible(eqn);
  return 0;
}

int Cond_lim_k_complique_transition_flux_nul_demi::initialiser(double temps)
{
  h_.resize(0,domaine_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(h_);

  h_grad_.resize(0,domaine_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(h_grad_);

  K_.resize(0,domaine_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(K_);

  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");

  int nf = la_frontiere_dis->frontiere().nb_faces();
  for (int f =0 ; f < nf ; f++) K_(f, 0) = 0 ; // K is 0 on the wall

  return 1;
}

double Cond_lim_k_complique_transition_flux_nul_demi::T_ext(int i) const
{
  return K_(i,0);
}

double Cond_lim_k_complique_transition_flux_nul_demi::T_ext(int i, int j) const
{
  return K_(i,j);
}

double Cond_lim_k_complique_transition_flux_nul_demi::h_imp(int i) const
{
  return h_(i,0);
}

double Cond_lim_k_complique_transition_flux_nul_demi::h_imp(int i, int j) const
{
  return h_(i,j);
}

double Cond_lim_k_complique_transition_flux_nul_demi::h_imp_grad(int i) const
{
  return h_grad_(i,0);
}

double Cond_lim_k_complique_transition_flux_nul_demi::h_imp_grad(int i, int j) const
{
  return h_grad_(i,j);
}

void Cond_lim_k_complique_transition_flux_nul_demi::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

void Cond_lim_k_complique_transition_flux_nul_demi::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Domaine_Poly_base& domaine = ref_cast(Domaine_Poly_base, domaine_Cl_dis().equation().domaine_dis().valeur());
  const DoubleTab&       y = corr_loi_paroi.get_tab("y");
  const DoubleTab&      mu = sub_type(Op_Diff_PolyMAC_base, domaine_Cl_dis().equation().operateur(0).l_op_base()) ? ref_cast(Op_Diff_PolyMAC_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu() :
                             ref_cast(Op_Diff_PolyMAC_P0_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu(),
                             &nu_visc = ref_cast(Navier_Stokes_std, domaine_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs(),
                              &vit = domaine_Cl_dis().equation().probleme().get_champ("vitesse").valeurs();

  int nf = la_frontiere_dis->frontiere().nb_faces(), f1 = la_frontiere_dis->frontiere().num_premiere_face();
  int D = dimension, nb_faces_tot = domaine.nb_faces_tot() ;
  const DoubleTab& n_f = domaine.face_normales();
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();

  if (mu.nb_dim() >= 3) Process::exit("Cond_lim_k_simple : transport of k must be SGDH !");

  int n = 0 ; // Carrying phase is 0 for turbulent flows

  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      int e_domaine = f_e(f_domaine,0);

      double u_orth = 0 ;
      for (int d = 0; d <D ; d++) u_orth -= vit(nb_faces_tot + e_domaine * D+d, n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -

      DoubleTrav u_parallel(D);
      for (int d = 0 ; d < D ; d++) u_parallel(d) = vit(nb_faces_tot + e_domaine * D + d, n) - u_orth*(-n_f(f_domaine,d))/fs(f_domaine) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
      double norm_u_parallel = std::sqrt(domaine.dot(&u_parallel(0), &u_parallel(0)));

      double u_tau_demi = corr_loi_paroi.calc_u_tau_loc(norm_u_parallel, nu_visc(e_domaine, 0), y(f_domaine, 0)/2.);
      double y_p = y(f_domaine, 0) * norm_u_parallel / nu_visc(e_domaine, 0);

      h_(f, 0) = 2.*mu(e_domaine, 0)/y(f_domaine, 0)   * ( 1 - std::tanh(  std::pow( y_p/600.,3)  ) );
      h_grad_(f, 0) = 2./y(f_domaine, 0)            * ( 1 - std::tanh(  std::pow( y_p/600.,3)  ) );
      K_(f, 0) = calc_k(y(f_domaine, 0)/2., u_tau_demi, nu_visc(e_domaine, 0));

    }

  h_.echange_espace_virtuel();
  h_grad_.echange_espace_virtuel();
  K_.echange_espace_virtuel();
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

