/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Echange_global_impose_adaptatif.cpp
// Directory:   $TRUST_ROOT/src/ThSol
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_global_impose_adaptatif.h>
#include <Motcle.h>
#include <Frontiere.h>
#include <Pb_Multiphase.h>
#include <Energie_Multiphase.h>
#include <Loi_paroi_adaptative.h>
#include <Zone_CoviMAC.h>
#include <math.h>

Implemente_instanciable(Echange_global_impose_adaptatif,"Paroi_echange_global_impose_adaptatif",Echange_global_impose);

Sortie& Echange_global_impose_adaptatif::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Echange_global_impose_adaptatif::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("Prandtl_turbulent", &Pr_turb);
  param.ajouter("Diametre_hydraulique", &Diam_hyd_);
  param.ajouter("von_karman", &von_karman_);
  param.lire_avec_accolades_depuis(s);

  Motcle motlu;
  s >> motlu;
  if (motlu != Motcle("T_ext")) Process::exit("Echange_global_impose_adaptatif : the field T_ext must be specified !");
  s >> T_ext();

  h_imp_.typer("Champ_front_vide");

  return s ;
}

int Echange_global_impose_adaptatif::initialiser(double temps)
{
  valeurs_h_.resize(0,zone_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis_.valeur().frontiere().creer_tableau_faces(valeurs_h_);
  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, zone_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");
  return 1;
}

int Echange_global_impose_adaptatif::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Thermique = "Thermique";
  if (dom_app == Thermique) return 1;
  else {err_pas_compatible(eqn); return 0;}
}

void Echange_global_impose_adaptatif::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis_.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis_.valeur().frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++) for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

void Echange_global_impose_adaptatif::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

void Echange_global_impose_adaptatif::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, zone_Cl_dis().equation().zone_dis().valeur());
  const DoubleTab& u_tau = corr_loi_paroi.get_tab("u_tau"); // y_p est numerote selon les faces de la zone
  const DoubleTab& y = corr_loi_paroi.get_tab("y"); // y_p est numerote selon les faces de la zone
  const DoubleTab& visc  = ref_cast(Navier_Stokes_std, zone_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs();
  const DoubleTab& cond  = ref_cast(Energie_Multiphase, zone_Cl_dis().equation().probleme().equation(2)).diffusivite_pour_pas_de_temps().valeurs();
  const DoubleTab& rho   = zone_Cl_dis().equation().probleme().get_champ("masse_volumique").valeurs();
  const DoubleTab& Cp    = zone_Cl_dis().equation().probleme().get_champ("capacite_calorifique").valeurs();
  int nf = la_frontiere_dis_.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis_.valeur().frontiere().num_premiere_face();
  int N = zone_Cl_dis().equation().inconnue().valeurs().line_size() ;
  const IntTab& f_e = zone.face_voisins();

  for (int f =0 ; f < nf ; f++) for (int n =1 ; n < N ; n++) // Reflechir a comment gerer les phases secondaires
      {
        int f_zone = f + f1; // number of the face in the zone
        int e_zone = f_e(f_zone,0);
        double theta_plus = calc_theta_plus(y(f_zone, 0), u_tau(f_zone, 0), visc(e_zone, 0), cond(e_zone, 0));
        valeurs_h_(f, 0) = u_tau(f_zone, 0) * rho(e_zone, n) * Cp(e_zone, n) * theta_plus;
      }
  valeurs_h_.echange_espace_virtuel();
}

double Echange_global_impose_adaptatif::calc_theta_plus(double y, double u_tau, double visc, double cond)
{
  double Prandtl = visc/cond;
  double y_plus = y * u_tau/visc;
  double beta = pow(3.85*pow(Prandtl, 1./3)-1.3, 2) + 2.12*log(Prandtl);
  double gamma = 0.01*(Prandtl*y_plus)*(Prandtl*y_plus)*(Prandtl*y_plus)*(Prandtl*y_plus)/(1+5*Prandtl*Prandtl*Prandtl*y_plus);
  double y_on_D_h = (Diam_hyd_>0)? y/Diam_hyd_:0;
  return Prandtl*y_plus*exp(-gamma) + (2.12*log((1+y_plus)*1.5*(2-y_on_D_h)/(1+2*(1-y_on_D_h)*(1-y_on_D_h)))+beta)*exp(-gamma);
}