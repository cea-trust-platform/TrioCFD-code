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
// File:        Cond_lim_omega_dix.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/28
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_lim_omega_dix.h>

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Loi_paroi_adaptative.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Motcle.h>

#include <math.h>

Implemente_instanciable(Cond_lim_omega_dix,"Cond_lim_omega_dix",Dirichlet_loi_paroi);
// XD Cond_lim_omega_dix condlim_base Cond_lim_omega_dix 1 Adaptive wall law boundary condition for turbulent dissipation rate

Sortie& Cond_lim_omega_dix::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Cond_lim_omega_dix::readOn(Entree& s )
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

void Cond_lim_omega_dix::completer()
{
  if (sub_type(Echelle_temporelle_turbulente, domaine_Cl_dis().equation())) Process::exit(que_suis_je() + " : you cannot define such a BC for tau equation, only for omega. Use scalaire_impose_paroi Champ_front_uniforme 1 0 for tau.");
  else if (!sub_type(Taux_dissipation_turbulent, domaine_Cl_dis().equation())) Process::exit(que_suis_je() + " : equation must be tau/omega !");

  int N = domaine_Cl_dis().equation().inconnue().valeurs().line_size();
  if (N > 1)  Process::exit(que_suis_je() + " : Only one phase for turbulent wall law is coded for now");

  correlation_loi_paroi_ = domaine_Cl_dis().equation().probleme().get_correlation("Loi_paroi");
}

void Cond_lim_omega_dix::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_Cl_dis().equation().domaine_dis());
  const DoubleTab&   u_tau = corr_loi_paroi.get_tab("u_tau");
  const DoubleTab&      nu_visc = ref_cast(Convection_diffusion_turbulence_multiphase, domaine_Cl_dis().equation()).diffusivite_pour_pas_de_temps().passe();

  const int cnu = nu_visc.dimension(0) == 1;
  int nf = la_frontiere_dis->frontiere().nb_faces(), f1 = la_frontiere_dis->frontiere().num_premiere_face();
  const IntTab& f_e = domaine.face_voisins();

  int n = 0 ; // Carrying phase is 0 for turbulent flows

  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      int e_domaine = (f_e(f_domaine,0)>=0) ? f_e(f_domaine,0) : f_e(f_domaine,1) ; // Make orientation vdf-proof
      double y_loc = f_e(f_domaine,0)>=0 ? domaine.dist_face_elem0(f_domaine,e_domaine) : domaine.dist_face_elem1(f_domaine,e_domaine) ;

      d_(f, n) = facteur_paroi_*calc_omega(y_loc, u_tau(f_domaine, n), nu_visc(!cnu * e_domaine, n));
    }
  d_.echange_espace_virtuel();
}

double Cond_lim_omega_dix::calc_omega(double y, double u_tau, double visc)
{
  return 6 * visc / (beta_omega * y * y);
}

