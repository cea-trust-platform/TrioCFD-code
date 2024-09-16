/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Source_Transport_K_Eps_Bas_Reynolds_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_Bas_Reynolds_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <Fluide_base.h>
#include <Domaine_Cl_VEF.h>
#include <TRUSTTrav.h>
#include <Domaine_VEF.h>
#include <Debog.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_VEF_Face,"Source_Transport_K_Eps_Bas_Reynolds_VEF_P1NC",Source_Transport_VEF_Face_base);

Sortie& Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::readOn(Entree& is ) { return Source_Transport_VEF_Face_base::readOn(is); }

void Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::associer_pb(const Probleme_base& pb )
{
  Source_Transport_VEF_Face_base::associer_pb(pb);
  eqn_keps_bas_re = ref_cast(Transport_K_Eps_Bas_Reynolds,equation());
}

DoubleTab& Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter(DoubleTab& resu) const
{
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter resu 0", resu);
  const Domaine_Cl_dis& zcl_keps = eqn_keps_bas_re->domaine_Cl_dis();
  const Domaine_dis_base& domaine_dis_keps = eqn_keps_bas_re->domaine_dis();
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const Domaine_Cl_VEF& domaine_Cl_VEF = le_dom_Cl_VEF.valeur();
  const DoubleTab& K_eps_Bas_Re = eqn_keps_bas_re->inconnue()->valeurs();
  const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds, eqn_keps_bas_re->modele_turbulence());
  const DoubleTab& visco_turb = mod_turb.viscosite_turbulente()->valeurs();
  const Modele_Fonc_Bas_Reynolds_Base& mon_modele_fonc = mod_turb.associe_modele_fonction().valeur();
  const Fluide_base& fluide = ref_cast(Fluide_base, eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& vit = eq_hydraulique->inconnue()->valeurs();
  const DoubleVect& vol_ent = domaine_VEF.volumes_entrelaces();
  const int nb_faces = domaine_VEF.nb_faces();

  DoubleTrav P(nb_faces), D(vol_ent), E(vol_ent), F1(nb_faces), F2(nb_faces);

  mon_modele_fonc.Calcul_D(D, domaine_dis_keps, zcl_keps, vit, K_eps_Bas_Re, ch_visco_cin);
  D.echange_espace_virtuel();
  mon_modele_fonc.Calcul_E(E, domaine_dis_keps, zcl_keps, vit, K_eps_Bas_Re, ch_visco_cin, visco_turb);
  mon_modele_fonc.Calcul_F2(F2, D, domaine_dis_keps, K_eps_Bas_Re, ch_visco_cin);

  if (_interpolation_viscosite_turbulente != 0)
    {
      Cerr << "Error 'interpolation_viscosite_turbulente' must be equal to '0' in this case." << finl;
      Process::exit();
    }
  calculer_terme_production_K(domaine_VEF, domaine_Cl_VEF, P, K_eps_Bas_Re, vit, visco_turb, _interpolation_viscosite_turbulente, _coefficient_limiteur);

  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter P 0", P);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter D 0", D);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter E 0", E);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter F1 0", F1);
  Debog::verifier("Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter F2 0", F2);

  for (int num_face = 0; num_face < nb_faces; num_face++)
    {
      if (K_eps_Bas_Re(num_face, 0) >= 1.e-20 && K_eps_Bas_Re(num_face, 1) > 1.e-20)
        {
          resu(num_face, 0) += (P(num_face) - K_eps_Bas_Re(num_face, 1) - D(num_face)) * vol_ent(num_face);
          resu(num_face, 1) += ((C1 * F1(num_face) * P(num_face) - C2 * F2(num_face) * K_eps_Bas_Re(num_face, 1)) * K_eps_Bas_Re(num_face, 1) / K_eps_Bas_Re(num_face, 0) + E(num_face))
                               * vol_ent(num_face);
        }
      else
        {
          resu(num_face, 0) += 0.;
          resu(num_face, 1) += 0.;
        }
    }
  return resu;
}
