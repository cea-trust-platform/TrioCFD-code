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

#include <Cond_lim_k_simple_transition_constante_Dirichlet.h>
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
#include <Domaine_VF.h>

#include <math.h>

Implemente_instanciable(Cond_lim_k_simple_transition_constante_Dirichlet,"Cond_lim_k_simple_transition_constante_Dirichlet",Dirichlet_loi_paroi);

Sortie& Cond_lim_k_simple_transition_constante_Dirichlet::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Cond_lim_k_simple_transition_constante_Dirichlet::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("beta_k", &beta_k_);
  param.ajouter("von_karman", &von_karman_);
  param.lire_avec_accolades_depuis(s);

  le_champ_front.typer("Champ_front_vide");

  return s;
}

void Cond_lim_k_simple_transition_constante_Dirichlet::completer()
{
}

void Cond_lim_k_simple_transition_constante_Dirichlet::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++)
    for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

int Cond_lim_k_simple_transition_constante_Dirichlet::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Turbulence="Turbulence";

  if (dom_app==Turbulence)
    return 1;
  else err_pas_compatible(eqn);
  return 0;
}

void Cond_lim_k_simple_transition_constante_Dirichlet::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

int Cond_lim_k_simple_transition_constante_Dirichlet::initialiser(double temps)
{
  d_.resize(0,domaine_Cl_dis().equation().inconnue().valeurs().line_size());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(d_);

  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");

  return 1;
}

void Cond_lim_k_simple_transition_constante_Dirichlet::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_Cl_dis().equation().domaine_dis().valeur());
  const DoubleTab&   u_tau = corr_loi_paroi.get_tab("u_tau");
  const DoubleTab&       y = corr_loi_paroi.get_tab("y");
  const DoubleTab&      mu = sub_type(Op_Diff_PolyMAC_base, domaine_Cl_dis().equation().operateur(0).l_op_base()) ? ref_cast(Op_Diff_PolyMAC_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu() :
                             ref_cast(Op_Diff_PolyMAC_P0_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu() ,
                             &nu_visc = ref_cast(Navier_Stokes_std, domaine_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs();

  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = domaine_Cl_dis().equation().inconnue().valeurs().line_size();
  const IntTab& f_e = domaine.face_voisins();

  if (mu.nb_dim() >= 3) Process::exit(que_suis_je() + " : transport of k must be SGDH !");
  if (N > 1)  Process::exit(que_suis_je() + " : Only one phase for turbulent wall law is coded for now");

  int n = 0 ; // Carrying phase is 0 for turbulent flows

  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      int e_domaine = f_e(f_domaine,0);

      d_(f, n) = calc_k(y(f_domaine, n), u_tau(f_domaine, n), nu_visc(e_domaine, n));
    }
  d_.echange_espace_virtuel();
}

double Cond_lim_k_simple_transition_constante_Dirichlet::calc_k(double y, double u_tau, double visc)
{
  double y_p = y * u_tau / visc;
  double b = std::pow(   (1+std::tanh(   1.5 * std::log(y_p)-1.5   ))/2   ,3.); // std::log(y_p)-1.2 avant modif ; std::log(y_p)-1. marchait pas mal ;

  return u_tau*u_tau*std::max(1/std::sqrt(beta_k_)*b, 0.);
}

