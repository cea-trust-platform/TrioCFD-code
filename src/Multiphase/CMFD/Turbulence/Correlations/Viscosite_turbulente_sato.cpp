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
// File:        Viscosite_turbulente_sato.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Viscosite_turbulente_sato.h>
#include <Pb_Multiphase.h>
#include <Masse_ajoutee_base.h>
#include <TRUSTTab_parts.h>
#include <Probleme_base.h>
#include <Champ_base.h>
#include <Param.h>
#include <Champ_Face_PolyMAC_P0.h>

Implemente_instanciable(Viscosite_turbulente_sato, "Viscosite_turbulente_sato", Viscosite_turbulente_base);

Sortie& Viscosite_turbulente_sato::printOn(Sortie& os) const
{
  return os;
}

Entree& Viscosite_turbulente_sato::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("coef_sato", &coef_sato);
  param.lire_avec_accolades_depuis(is);

  //identification des phases
  Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr ;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;

  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  pbm->creer_champ("distance_paroi_globale"); // Besoin de distance a la paroi

  return is;
}

void Viscosite_turbulente_sato::eddy_viscosity(DoubleTab& nu_t) const
{
  nu_t = 0; // La contribution de nu_sato est ajouter via un terme source dans la qdm
}

void Viscosite_turbulente_sato::reynolds_stress(DoubleTab& R_ij) const // Renvoie <u_i'u_j'>
{

  // Récupère les champs de vitesse - alpha - d_bulles (diamètre des bulles) - grad vitesse
  const Pb_Multiphase* pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr ;
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, pb_->domaine_dis());
  const DoubleTab& tab_u = pb_->get_champ("vitesse").passe();
  const DoubleTab& d_bulles = pb_->get_champ("diametre_bulles").passe();
  const DoubleTab& alpha = pb_->get_champ("alpha").passe();
  const DoubleTab& tab_grad = pbm->get_champ("gradient_vitesse").passe();

  int D = dimension;
  int nf_tot = domaine.nb_faces_tot();
  int N = alpha.dimension(1);

  // Déclarer de nouvelles variables

  double S_ij; // taux de déformation : S = 0.5*(gradV+gradV^T)
  double nu_sato; // viscosité de sato
  // Champ de vitesse
  ConstDoubleTab_parts p_u(tab_u); //en PolyMAC_P0, tab_u contient (nf.u) aux faces, puis (u_i) aux elements
  int i_part = -1;
  for (int i = 0; i < p_u.size(); i++)
    if (p_u[i].get_md_vector() == R_ij.get_md_vector()) i_part = i; //on cherche une partie ayant le meme support
  if (i_part < 0) Process::exit("Viscosite_turbulente_sato : inconsistency between velocity and Rij!");
  const DoubleTab& u = p_u[i_part]; //le bon tableau
  // vitesse relative
  DoubleTrav u_r(R_ij.dimension(0), 1);
  double u_r_carre;

  // Tenseur Reynolds de Sato

  for( int e = 0 ; e < R_ij.dimension(0) ; e++) // elements
    for (int i = 0; i < D; i++)                 // dimension i
      for (int j = 0; j < D; j++)               // dimension j
        for (int k = 0; k < N; k++)             // phases
          if (k!=n_l)                           // iteration sur phases gazeuses
            {
              // Calcul de la norme de la vitesse relative
              u_r_carre = 0.;
              for (int d = 0; d < D; d++) u_r_carre += (u(e, d, k) - u(e, d, n_l))*(u(e, d, k) - u(e, d, n_l)); // relative speed = gas speed - liquid speed
              u_r(e, 0) = std::sqrt(u_r_carre);
              // Calcul de la viscosite de Sato
              nu_sato = coef_sato * alpha(e, k) * d_bulles(e, k) * u_r(e,0);
              // Tenseur des taux de déformations Sij = gradv_ij + gradv_ji (facteur 1/2 supprimé car il se simplifie après)
              S_ij = tab_grad(nf_tot + i + e * D , D * n_l + j) + tab_grad(nf_tot + j + e * D , D * n_l + i);
              // Tenseur de Reynolds BIA
              R_ij(e, k, i, j) = 0 ; // No BIT for gas phase
              R_ij(e, 0, i, j) -= nu_sato * S_ij; // <u_i'u_j'> = - nu_sato * S_ij
            }
}

void Viscosite_turbulente_sato::k_over_eps(DoubleTab& k_sur_eps) const
{
  k_sur_eps =  0;
}

void Viscosite_turbulente_sato::eps(DoubleTab& eps_) const
{
  eps_ =  0;
}

