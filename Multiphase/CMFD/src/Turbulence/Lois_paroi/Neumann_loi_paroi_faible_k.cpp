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

#include <Neumann_loi_paroi_faible_k.h>
#include <Motcle.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Convection_Diffusion_Concentration.h>
#include <Loi_paroi_adaptative.h>
#include <Frontiere_dis_base.h>
#include <Frontiere.h>
#include <Pb_Multiphase.h>
#include <Navier_Stokes_std.h>
#include <Op_Diff_CoviMAC_base.h>
#include <Zone_CoviMAC.h>
#include <Energie_cinetique_turbulente.h>

#include <math.h>

Implemente_instanciable(Neumann_loi_paroi_faible_k,"Neumann_loi_paroi_faible_k",Neumann_loi_paroi);

Sortie& Neumann_loi_paroi_faible_k::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Neumann_loi_paroi_faible_k::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.ajouter("beta_k", &beta_k);
  param.ajouter("von_karman", &von_karman_);
  param.lire_avec_accolades_depuis(s);

  le_champ_front.typer("Champ_front_vide");

  return s;
}

void Neumann_loi_paroi_faible_k::completer()
{
  if (!sub_type(Energie_cinetique_turbulente, zone_Cl_dis().equation())) Process::exit("Neumann_loi_paroi_faible_k : equation must be k !");
}

void Neumann_loi_paroi_faible_k::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++) for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

int Neumann_loi_paroi_faible_k::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Turbulence="Turbulence";

  if (dom_app==Turbulence)
    return 1;
  else err_pas_compatible(eqn);
  return 0;
}

double Neumann_loi_paroi_faible_k::flux_impose(int i) const
{
  return valeurs_flux_(i,0);
}

double Neumann_loi_paroi_faible_k::flux_impose(int i,int j) const
{
  return valeurs_flux_(i,j);
}

int Neumann_loi_paroi_faible_k::initialiser(double temps)
{
  valeurs_flux_.resize(0,zone_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(valeurs_flux_);
  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, zone_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");
  return 1;
}

void Neumann_loi_paroi_faible_k::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

void Neumann_loi_paroi_faible_k::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, zone_Cl_dis().equation().zone_dis().valeur());
  const DoubleTab&   u_tau = corr_loi_paroi.get_tab("u_tau");
  const DoubleTab&       y = corr_loi_paroi.get_tab("y");
  const DoubleTab&  visc_c = ref_cast(Navier_Stokes_std, zone_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs();
  const DoubleTab&      mu = ref_cast(Op_Diff_CoviMAC_base, zone_Cl_dis().equation().operateur(0).l_op_base()).nu();

  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = zone_Cl_dis().equation().inconnue().valeurs().line_size();
  const IntTab& f_e = zone.face_voisins();

  if (mu.nb_dim() >= 3) Process::exit("Neumann_loi_paroi_faible_k : transport of k must be SGDH !");

  for (int f =0 ; f < nf ; f++)
    {
      int f_zone = f + f1; // number of the face in the zone
      int e_zone = f_e(f_zone,0);
      valeurs_flux_(f, 0) = mu(e_zone, 0) * u_tau(f_zone, 0)*u_tau(f_zone, 0)*u_tau(f_zone, 0)/visc_c(e_zone, 0) *calc_dyplus_kplus(y(f_zone, 0)*u_tau(f_zone, 0)/visc_c(e_zone, 0)); // flux de Neumann = -(-mu * dy_k) car 1er - avec orientation face de bord et 2e car flux selon - grad
    }
  for (int n =1 ; n < N ; n++) for (int f =0 ; f < nf ; f++)
      {
        Process::exit("Neumann_loi_paroi_faible_k : Only one phase for turbulent wall law is coded for now");
      }

  valeurs_flux_.echange_espace_virtuel();
}

double Neumann_loi_paroi_faible_k::calc_dyplus_kplus(double y_p)
{
  double arg = y_p/10.;
  double phi = std::tanh(arg*arg*arg*arg);
  double w_vis_p = 6./(beta_omega*y_p*y_p);
  double w_log_p = 1/(std::sqrt(beta_k)*von_karman_*y_p);
  double w_1_p = w_vis_p+w_log_p;
  double w_2_p = std::pow( std::pow(w_vis_p,1.2)+std::pow(w_log_p,1.2),1./1.2);
  double w = phi * w_1_p + (1-phi)*w_2_p;

  double d_w_vis = - 12. / (beta_omega * y_p*y_p*y_p);
  double d_w_log = - 1 / ( std::sqrt(beta_k) * von_karman_ * y_p*y_p);
  double d_w_1 = d_w_vis + d_w_log ;
  double d_w_2 = ( d_w_vis * std::pow(w_vis_p,0.2) + d_w_log * std::pow(w_log_p,0.2))*std::pow( std::pow(w_vis_p,1.2) + std::pow(w_log_p,1.2), 1./1.2-1.);
  double d_phi = 4./10.*(arg*arg*arg)*(1-phi*phi);
  double d_w = d_phi * w_1_p + phi * d_w_1 - d_phi * w_2_p + (1-phi) * d_w_2 ;

  double d_u_plus   = deriv_u_plus_de_y_plus(y_p);
  double d_2_u_plus = deriv_u_plus_de_y_plus_2(y_p);

  double kplus_ypgrand = (1/d_u_plus-1)*w;
  double d_kplus_ypgrand = (1/d_u_plus-1)*d_w + w*(-d_2_u_plus)/(d_u_plus*d_u_plus);


  double arg2 = y_p/4.;
  double phi2 = std::tanh(arg2*arg2*arg2*arg2);
  double d_phi2 = arg2*arg2*arg2*(1-phi2*phi2);

  double kplus_yppetit = 5e-2*std::pow(y_p,1.5);
  double d_kplus_yppetit = 1.5*5e-2*std::pow(y_p,0.5);

  return phi2 * d_kplus_ypgrand + d_phi2 * kplus_ypgrand + (1-phi2)*d_kplus_yppetit - d_phi2 * kplus_yppetit;
}

double Neumann_loi_paroi_faible_k::deriv_u_plus_de_y_plus(double y_p)
{
  double reichardt = std::log(1+0.4*y_p)/von_karman_;
  reichardt += 7.8;
  reichardt += -7.8*std::exp(-y_p/11);
  reichardt += -7.8*y_p/11*std::exp(-y_p/3);

  double log_law = log(y_p+limiteur_y_p)/von_karman_ + 5.1;

  double blending = tanh( y_p/27.*y_p/27.*y_p/27.*y_p/27.);

  double d_reichardt = 0.4/(1+0.4*y_p)*1/von_karman_;
  d_reichardt += 7.8/11*std::exp(-y_p/11);
  d_reichardt += -7.8/11*std::exp(-y_p/3) + 7.8*y_p/33*std::exp(-y_p/3) ;
  double d_log_law = 1/(y_p+limiteur_y_p)*1/von_karman_;
  double d_blending = 4./27.*(y_p/27.*y_p/27.*y_p/27.)*(1-blending*blending);

  return (1-blending)*d_reichardt - reichardt*d_blending + blending*d_log_law + d_blending*log_law ;
}

double Neumann_loi_paroi_faible_k::deriv_u_plus_de_y_plus_2(double y_p)
{
  double reichardt = std::log(1+0.4*y_p)/von_karman_ + 7.8 -7.8*std::exp(-y_p/11) -7.8*y_p/11*std::exp(-y_p/3);
  double log_law = std::log(y_p+limiteur_y_p)/von_karman_ + 5.1;
  double blending = std::tanh( y_p/27.*y_p/27.*y_p/27.*y_p/27.);

  double d_reichardt = 0.4/(1+0.4*y_p)*1/von_karman_;
  d_reichardt += 7.8/11*std::exp(-y_p/11);
  d_reichardt += -7.8/11*std::exp(-y_p/3) + 7.8*y_p/33.*std::exp(-y_p/3);
  double d_log_law = 1/(y_p+limiteur_y_p)*1/von_karman_;
  double d_blending = 4./27.*(y_p/27*y_p/27*y_p/27)*(1-blending*blending);

  double d_2_reichardt = -0.4*0.4/((1+0.4*y_p)*(1+0.4*y_p))*1/von_karman_;
  d_2_reichardt += -7.8/(11.*11.)*std::exp(-y_p/11);
  d_2_reichardt += 7.8/33*std::exp(-y_p/3) + 7.8/33*std::exp(-y_p/3) - 7.8*y_p/99*std::exp(-y_p/3) ;
  double d_2_log_law = -1/((y_p+limiteur_y_p)*(y_p+limiteur_y_p))*1/von_karman_;
  double d_2_blending = 3.*4./(27.*27.)*(y_p/27.*y_p/27.)*(1-blending*blending) + 4./27.*(y_p/27.*y_p/27.*y_p/27.)*(-2.*d_blending*blending);

  double rep = (1-blending)*d_2_reichardt-2*d_reichardt*d_blending-reichardt*d_2_blending;
  rep += blending*d_2_log_law + 2*d_blending*d_log_law + d_2_blending*log_law;
  return rep;
}
