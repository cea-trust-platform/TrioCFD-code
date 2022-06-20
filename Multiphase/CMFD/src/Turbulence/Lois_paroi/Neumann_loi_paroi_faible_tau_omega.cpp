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

#include <Neumann_loi_paroi_faible_tau_omega.h>
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

#include <math.h>

Implemente_instanciable(Neumann_loi_paroi_faible_tau_omega,"Neumann_loi_paroi_faible_tau|Neumann_loi_paroi_faible_omega",Neumann_loi_paroi);
// XD Neumann_loi_paroi_faible_omega condlim_base Neumann_loi_paroi_faible_omega 1 not_set
// XD Neumann_loi_paroi_faible_tau condlim_base Neumann_loi_paroi_faible_tau 1 not_set

Sortie& Neumann_loi_paroi_faible_tau_omega::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Neumann_loi_paroi_faible_tau_omega::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.ajouter("beta_k", &beta_k);
  param.ajouter("von_karman", &von_karman_);
  param.lire_avec_accolades_depuis(s);

  le_champ_front.typer("Champ_front_vide");

  return s;
}

void Neumann_loi_paroi_faible_tau_omega::completer()
{
  if (sub_type(Echelle_temporelle_turbulente, zone_Cl_dis().equation())) is_tau_ = 1;
  else if (sub_type(Taux_dissipation_turbulent, zone_Cl_dis().equation())) is_tau_ = 0;
  else Process::exit("Neumann_loi_paroi_faible_tau_omega : equation must be tau/omega !");
}

void Neumann_loi_paroi_faible_tau_omega::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++)
    for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

int Neumann_loi_paroi_faible_tau_omega::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Turbulence="Turbulence";

  if (dom_app==Turbulence)
    return 1;
  else err_pas_compatible(eqn);
  return 0;
}

double Neumann_loi_paroi_faible_tau_omega::flux_impose(int i) const
{
  return valeurs_flux_(i,0);
}

double Neumann_loi_paroi_faible_tau_omega::flux_impose(int i,int j) const
{
  return valeurs_flux_(i,j);
}

int Neumann_loi_paroi_faible_tau_omega::initialiser(double temps)
{
  valeurs_flux_.resize(0,zone_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(valeurs_flux_);
  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, zone_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");
  return 1;
}

void Neumann_loi_paroi_faible_tau_omega::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

void Neumann_loi_paroi_faible_tau_omega::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Zone_VF& zone = ref_cast(Zone_VF, zone_Cl_dis().equation().zone_dis().valeur());
  const DoubleTab&   u_tau = corr_loi_paroi.get_tab("u_tau");
  const DoubleTab&       y = corr_loi_paroi.get_tab("y");
  const DoubleTab& visc_c  = ref_cast(Navier_Stokes_std, zone_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs();
  const DoubleTab&      mu = sub_type(Op_Diff_PolyMAC_base, zone_Cl_dis().equation().operateur(0).l_op_base()) ? ref_cast(Op_Diff_PolyMAC_base, zone_Cl_dis().equation().operateur(0).l_op_base()).nu() :
                             ref_cast(Op_Diff_PolyMAC_P0_base, zone_Cl_dis().equation().operateur(0).l_op_base()).nu() ;

  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = zone_Cl_dis().equation().inconnue().valeurs().line_size() ;
  const IntTab& f_e = zone.face_voisins();

  if (mu.nb_dim() >= 3) Process::exit("Neumann_loi_paroi_faible_tau_omega : transport of tau/omega must be SGDH !");
  if (N > 1)  Process::exit("Neumann_loi_paroi_faible_tau : Only one phase for turbulent wall law is coded for now");

  if (is_tau_ == 1)
    {
      if (sub_type(Op_Diff_Tau_PolyMAC_P0_Elem, zone_Cl_dis().equation().operateur(0).l_op_base()))
        {
          double limiter = ref_cast(Op_Diff_Tau_PolyMAC_P0_Elem, zone_Cl_dis().equation().operateur(0).l_op_base()).limiter_tau();
          const DoubleTab& tau = zone_Cl_dis().equation().inconnue().passe();
          for (int f =0 ; f < nf ; f++)
            {
              int f_zone = f + f1; // number of the face in the zone
              int e_zone = f_e(f_zone,0);
              valeurs_flux_(f, 0) = - mu(e_zone, 0) * dy_tau(y(f_zone, 0), u_tau(f_zone, 0), visc_c(e_zone, 0)) / ((tau(e_zone,0) > limiter) ? limiter/(tau(e_zone,0)*tau(e_zone,0)) : 1/limiter); // flux de Neumann = -mu * dy_tau car flux selon - grad ; besoin de multiplier par tau**2 à cause de la forme partiucliere de la diffusion
            }
        }
      else
        {
          for (int f =0 ; f < nf ; f++)
            {
              int f_zone = f + f1; // number of the face in the zone
              int e_zone = f_e(f_zone,0);
              valeurs_flux_(f, 0) = - mu(e_zone, 0) * dy_tau(y(f_zone, 0), u_tau(f_zone, 0), visc_c(e_zone, 0)) ; // flux de Neumann = -mu * dy_tau car flux selon - grad ; besoin de multiplier par tau**2 à cause de la forme partiucliere de la diffusion
            }
        }
    }
  if (is_tau_ == 0)
    {
      for (int f =0 ; f < nf ; f++)
        {
          int f_zone = f + f1; // number of the face in the zone
          int e_zone = f_e(f_zone,0);
          valeurs_flux_(f, 0) = - mu(e_zone, 0) * dy_omega(y(f_zone, 0), u_tau(f_zone, 0), visc_c(e_zone, 0)); // flux de Neumann = -mu * dy_omega car flux selon - grad
        }
    }

  valeurs_flux_.echange_espace_virtuel();
}

double Neumann_loi_paroi_faible_tau_omega::dy_tau(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double w_vis = 6 * visc / (beta_omega * y * y);
  double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
  double w_1 = w_vis + w_log ;
  double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
  double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

  double d_w_vis = - 12 * visc / (beta_omega * y * y * y);
  double d_w_log = - u_tau / ( std::sqrt(beta_k) * von_karman_ * y * y);
  double d_w_1 = d_w_vis + d_w_log ;
  double d_w_2 = ( d_w_vis * std::pow(w_vis, 0.2) + d_w_log * std::pow(w_log, 0.2) )*std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2-1 );
  double d_blending = 4*u_tau/(visc*10)*(y_p/10*y_p/10*y_p/10)*(1-blending*blending);

  return - ( d_blending * w_1 + blending * d_w_1 - d_blending * w_2 + (1-blending) * d_w_2 ) / ((blending * w_1 + (1-blending) * w_2) * (blending * w_1 + (1-blending) * w_2)) ;

}

double Neumann_loi_paroi_faible_tau_omega::dy_omega(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double w_vis = 6 * visc / (beta_omega * y * y);
  double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
  double w_1 = w_vis + w_log ;
  double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
  double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

  double d_w_vis = - 12 * visc / (beta_omega * y * y * y);
  double d_w_log = - u_tau / ( std::sqrt(beta_k) * von_karman_ * y * y);
  double d_w_1 = d_w_vis + d_w_log ;
  double d_w_2 = ( d_w_vis * std::pow(w_vis, 0.2) + d_w_log * std::pow(w_log, 0.2) )*std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2-1 );
  double d_blending = 4*u_tau/(visc*10)*(y_p/10*y_p/10*y_p/10)*(1-blending*blending);

  return d_blending * w_1 + blending * d_w_1 - d_blending * w_2 + (1-blending) * d_w_2 ;
}
