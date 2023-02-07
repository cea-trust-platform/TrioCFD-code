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
// File:        Cond_lim_tau_omega_simple_dix.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/28
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_lim_tau_omega_simple_dix.h>
#include <Motcle.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Convection_Diffusion_Concentration.h>
#include <Loi_paroi_adaptative.h>
#include <Frontiere_dis_base.h>
#include <Frontiere.h>
#include <Pb_Multiphase.h>
#include <Navier_Stokes_std.h>
#include <Domaine_VF.h>
#include <Operateur_Diff_base.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Op_Diff_PolyMAC_base.h>
#include <Op_Diff_PolyMAC_P0_base.h>
#include <TRUSTTrav.h>
#include <Domaine_Poly_base.h>

#include <math.h>

Implemente_instanciable(Cond_lim_tau_omega_simple_dix,"Cond_lim_tau_simple_dix|Cond_lim_omega_simple_dix",Dirichlet_loi_paroi);

Sortie& Cond_lim_tau_omega_simple_dix::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Cond_lim_tau_omega_simple_dix::readOn(Entree& s )
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

void Cond_lim_tau_omega_simple_dix::completer()
{
  if (sub_type(Echelle_temporelle_turbulente, domaine_Cl_dis().equation())) is_tau_ = 1;
  else if (sub_type(Taux_dissipation_turbulent, domaine_Cl_dis().equation())) is_tau_ = 0;
  else Process::exit(que_suis_je() + " : equation must be tau/omega !");

  int N = domaine_Cl_dis().equation().inconnue().valeurs().line_size();
  if (N > 1)  Process::exit(que_suis_je() + " : Only one phase for turbulent wall law is coded for now");
}

void Cond_lim_tau_omega_simple_dix::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++)
    for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

int Cond_lim_tau_omega_simple_dix::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Turbulence="Turbulence";

  if (dom_app==Turbulence)
    return 1;
  else err_pas_compatible(eqn);
  return 0;
}

void Cond_lim_tau_omega_simple_dix::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

int Cond_lim_tau_omega_simple_dix::initialiser(double temps)
{
  d_.resize(0,domaine_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(d_);

  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");

  return 1;
}

void Cond_lim_tau_omega_simple_dix::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Domaine_Poly_base& domaine = ref_cast(Domaine_Poly_base, domaine_Cl_dis().equation().domaine_dis().valeur());
  const DoubleTab&   u_tau = corr_loi_paroi.get_tab("u_tau");
  const DoubleTab&       y = corr_loi_paroi.get_tab("y");
  const DoubleTab&      nu_visc = ref_cast(Navier_Stokes_std, domaine_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs();

  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  const IntTab& f_e = domaine.face_voisins();

  int n = 0 ; // Carrying phase is 0 for turbulent flows

  if (is_tau_ == 1)
    {
      for (int f =0 ; f < nf ; f++)
        {
          int f_domaine = f + f1; // number of the face in the domaine
          int e_domaine = f_e(f_domaine,0);

          d_(f, n) = facteur_paroi_*calc_tau(y(f_domaine, n), u_tau(f_domaine, n), nu_visc(e_domaine, n));
        }
    }
  if (is_tau_ == 0)
    {
      for (int f =0 ; f < nf ; f++)
        {
          int f_domaine = f + f1; // number of the face in the domaine
          int e_domaine = f_e(f_domaine,0);

          d_(f, n) = facteur_paroi_*calc_omega(y(f_domaine, n), u_tau(f_domaine, n), nu_visc(e_domaine, n));
       }
    }
  d_.echange_espace_virtuel();
}

double Cond_lim_tau_omega_simple_dix::calc_tau(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double w_vis = 6 * visc / (beta_omega * y * y);
  double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
  double w_1 = w_vis + w_log ;
  double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
  double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

  return 1 / (blending * w_1 + (1-blending) * w_2) ;

}

double Cond_lim_tau_omega_simple_dix::calc_omega(double y, double u_tau, double visc)
{
  /*  double y_p = y * u_tau / visc;
    double w_vis = 6 * visc / (beta_omega * y * y);
    double w_log = u_tau / ( std::sqrt(beta_k) * von_karman_ * y);
    double w_1 = w_vis + w_log ;
    double w_2 = std::pow( std::pow(w_vis, 1.2) + std::pow(w_log, 1.2) , 1/1.2 );
    double blending = std::tanh( y_p/10*y_p/10*y_p/10*y_p/10);

    return blending * w_1 + (1-blending) * w_2 ;*/
  return 6 * visc / (beta_omega * y * y);
}

