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
// File:        Loi_paroi_faible_tau.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/28
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_lim_tau_omega_simple_demi.h>
#include <Motcle.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Convection_Diffusion_Concentration.h>
#include <Loi_paroi_adaptative.h>
#include <Frontiere_dis_base.h>
#include <Frontiere.h>
#include <Pb_Multiphase.h>
#include <Navier_Stokes_std.h>
#include <Zone_VF.h>
#include <Operateur_Diff_base.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Op_Diff_PolyMAC_base.h>
#include <Op_Diff_PolyMAC_P0_base.h>
#include <Op_Diff_Tau_PolyMAC_P0_Elem.h>
#include <TRUSTTrav.h>
#include <Zone_Poly_base.h>

#include <math.h>

Implemente_instanciable(Cond_lim_tau_omega_simple_demi,"Cond_lim_tau_simple_demi|Cond_lim_omega_simple_demi",Echange_global_impose);

Sortie& Cond_lim_tau_omega_simple_demi::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Cond_lim_tau_omega_simple_demi::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.ajouter("beta_k", &beta_k);
  param.ajouter("von_karman", &von_karman_);
  param.lire_avec_accolades_depuis(s);

  le_champ_front.typer("Champ_front_vide");

  return s;
}

void Cond_lim_tau_omega_simple_demi::completer()
{
  if (sub_type(Echelle_temporelle_turbulente, zone_Cl_dis().equation())) is_tau_ = 1;
  else if (sub_type(Taux_dissipation_turbulent, zone_Cl_dis().equation())) is_tau_ = 0;
  else Process::exit("Neumann_loi_paroi_faible_tau_omega : equation must be tau/omega !");
}

void Cond_lim_tau_omega_simple_demi::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++) for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

int Cond_lim_tau_omega_simple_demi::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Turbulence="Turbulence";

  if (dom_app==Turbulence)
    return 1;
  else err_pas_compatible(eqn);
  return 0;
}

double Cond_lim_tau_omega_simple_demi::T_ext(int i) const
{
  return d_(i,0);
}

double Cond_lim_tau_omega_simple_demi::T_ext(int i, int j) const
{
  return d_(i,j);
}

double Cond_lim_tau_omega_simple_demi::h_imp(int i) const
{
  return h_(i,0);
}

double Cond_lim_tau_omega_simple_demi::h_imp(int i, int j) const
{
  return h_(i,j);
}

void Cond_lim_tau_omega_simple_demi::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

int Cond_lim_tau_omega_simple_demi::initialiser(double temps)
{
  h_.resize(0,zone_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(h_);

  d_.resize(0,zone_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(d_);

  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, zone_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");

  return 1;
}

void Cond_lim_tau_omega_simple_demi::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Zone_Poly_base& zone = ref_cast(Zone_Poly_base, zone_Cl_dis().equation().zone_dis().valeur());
  const DoubleTab&       y = corr_loi_paroi.get_tab("y");
  const DoubleTab&      mu = sub_type(Op_Diff_PolyMAC_base, zone_Cl_dis().equation().operateur(0).l_op_base()) ? ref_cast(Op_Diff_PolyMAC_base, zone_Cl_dis().equation().operateur(0).l_op_base()).nu() :
                             ref_cast(Op_Diff_PolyMAC_P0_base, zone_Cl_dis().equation().operateur(0).l_op_base()).nu() ,
                             &nu_visc = ref_cast(Navier_Stokes_std, zone_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs(),
                              &vit = zone_Cl_dis().equation().probleme().get_champ("vitesse").valeurs();

  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = zone_Cl_dis().equation().inconnue().valeurs().line_size(), D = dimension ;
  int nb_faces_tot = zone.nb_faces_tot();
  const DoubleTab& n_f = zone.face_normales();
  const DoubleVect& fs = zone.face_surfaces();
  const IntTab& f_e = zone.face_voisins();

  if (mu.nb_dim() >= 3) Process::exit("Neumann_loi_paroi_faible_tau_omega : transport of tau/omega must be SGDH !");
  if (N > 1)  Process::exit("Neumann_loi_paroi_faible_tau : Only one phase for turbulent wall law is coded for now");

  int n = 0 ; // Carrying phase is 0 for turbulent flows

  if (is_tau_ == 1)
    {
      for (int f =0 ; f < nf ; f++)
        {
          int f_zone = f + f1; // number of the face in the zone
          int e_zone = f_e(f_zone,0);

          double u_orth = 0 ;
          for (int d = 0; d <D ; d++) u_orth -= vit(nb_faces_tot + e_zone * D+d, n)*n_f(f_zone,d)/fs(f_zone); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -

          DoubleTrav u_parallel(D);
          for (int d = 0 ; d < D ; d++) u_parallel(d) = vit(nb_faces_tot + e_zone * D + d, n) - u_orth*(-n_f(f_zone,d))/fs(f_zone) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          double norm_u_parallel = std::sqrt(zone.dot(&u_parallel(0), &u_parallel(0)));

          double u_tau_demi = corr_loi_paroi.calc_u_tau_loc(norm_u_parallel, nu_visc(e_zone, 0), y(f_zone, 0)/2.);
          h_(f, 0) = 2.*mu(e_zone, 0)/y(f_zone, 0) ;
          d_(f, 0) = calc_tau(y(f_zone, 0)/2., u_tau_demi, nu_visc(e_zone, 0));
        }
    }
  if (is_tau_ == 0)
    {
      for (int f =0 ; f < nf ; f++)
        {
          int f_zone = f + f1; // number of the face in the zone
          int e_zone = f_e(f_zone,0);

          double u_orth = 0 ;
          for (int d = 0; d <D ; d++) u_orth -= vit(nb_faces_tot + e_zone * D+d, n)*n_f(f_zone,d)/fs(f_zone); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -

          DoubleTrav u_parallel(D);
          for (int d = 0 ; d < D ; d++) u_parallel(d) = vit(nb_faces_tot + e_zone * D + d, n) - u_orth*(-n_f(f_zone,d))/fs(f_zone) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          double norm_u_parallel = std::sqrt(zone.dot(&u_parallel(0), &u_parallel(0)));

          double u_tau_demi = corr_loi_paroi.calc_u_tau_loc(norm_u_parallel, nu_visc(e_zone, 0), y(f_zone, 0)/2.);
          h_(f, 0) = 2.*mu(e_zone, 0)/y(f_zone, 0) ;
          d_(f, 0) = calc_omega(y(f_zone, 0)/2., u_tau_demi, nu_visc(e_zone, 0));
        }
    }

  h_.echange_espace_virtuel();
  d_.echange_espace_virtuel();
}

double Cond_lim_tau_omega_simple_demi::calc_tau(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double w_vis = 6 * visc / (beta_omega * y * y);
  double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
  double w_1 = w_vis + w_log ;
  double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
  double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

  return 1 / (blending * w_1 + (1-blending) * w_2) ;

}

double Cond_lim_tau_omega_simple_demi::calc_omega(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double w_vis = 6 * visc / (beta_omega * y * y);
  double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
  double w_1 = w_vis + w_log ;
  double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
  double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

  return blending * w_1 + (1-blending) * w_2 ;
}