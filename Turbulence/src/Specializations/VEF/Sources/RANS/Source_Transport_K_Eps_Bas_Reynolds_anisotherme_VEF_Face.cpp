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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_Face.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Probleme_base.h>
#include <Zone_Cl_VEF.h>
#include <Fluide_base.h>
#include <TRUSTTrav.h>
#include <Zone_VEF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_Face,"Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_P1NC",Source_Transport_K_Eps_Bas_Reynolds_VEF_Face);

Sortie& Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_Face::printOn(Sortie& s) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_Face::readOn(Entree& s) { return s ; }

void Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::verifier_milieu_anisotherme(pb,que_suis_je());
  Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::associer_pb(pb);
  //correction pour quasi compressible : on utilise pas l'eq_thermique en fait (ce erait une Convection_Diffusion_Enthalpie_Turbulent sinon)
  if (sub_type(Convection_Diffusion_Temperature,pb.equation(1)))
    {
      eq_thermique = ref_cast(Convection_Diffusion_Temperature,pb.equation(1));
      const Fluide_base& fluide_2 = eq_thermique->fluide();
      if (fluide_2.beta_t().non_nul()) beta_t = fluide_2.beta_t();
      gravite = fluide_2.gravite();
    }
  else gravite = pb.equation(0).milieu().gravite();
}

// Elie Saikali : TODO : FIXME : a factoriser avec Source_Transport_K_Eps_Bas_Reynolds_VEF_Face::ajouter
DoubleTab& Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VEF_Face::ajouter(DoubleTab& resu) const
{
  const Zone_Cl_dis& zcl = eq_hydraulique->zone_Cl_dis();
  const Zone_Cl_dis& zcl_keps = eqn_keps_bas_re->zone_Cl_dis();
  const Zone_dis& zone_dis_keps = eqn_keps_bas_re->zone_dis();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, eq_hydraulique->zone_dis().valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF, zcl.valeur());
  const Zone_Cl_VEF& zcl_VEF_th = ref_cast(Zone_Cl_VEF, eq_thermique->zone_Cl_dis().valeur());
  const DoubleTab& K_eps_Bas_Re = eqn_keps_bas_re->inconnue().valeurs();
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleTab& visco_turb = eqn_keps_bas_re->modele_turbulence().viscosite_turbulente().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base, eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs(), &g = gravite->valeurs();
  const Champ_Don& ch_beta = beta_t.valeur();
  const Fluide_base& fluide = ref_cast(Fluide_base, eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds, eqn_keps_bas_re->modele_turbulence());
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc = mod_turb.associe_modele_fonction();
  int nb_faces = zone_VEF.nb_faces();
  const DoubleVect& vol_ent = zone_VEF.volumes_entrelaces();

  DoubleTrav P(nb_faces), G(nb_faces), G1(nb_faces), D(nb_faces), E(nb_faces), F1(nb_faces), F2(nb_faces);

  mon_modele_fonc.Calcul_D(D,zone_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin);
  mon_modele_fonc.Calcul_E(E,zone_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin,visco_turb);
  mon_modele_fonc.Calcul_F2(F2,D,zone_dis_keps,K_eps_Bas_Re,ch_visco_cin);

  if (_interpolation_viscosite_turbulente != 0)
    {
      Cerr << "Error 'interpolation_viscosite_turbulente' must be equal to '0' in this case." << finl;
      Process::exit();
    }
  calculer_terme_production_K(zone_VEF,zone_Cl_VEF,P,K_eps_Bas_Re,vit,visco_turb, _interpolation_viscosite_turbulente);

  // C'est l'objet de type zone_Cl_dis de l'equation thermique qui est utilise dans le calcul de G
  // Nous utilisons le modele de fluctuation thermique pour le calcul du terme de destruction G.
  calculer_terme_destruction_K_gen(zone_VEF,zcl_VEF_th,G,scalaire,alpha_turb,ch_beta,g,0);

  for (int num_face=0; num_face<nb_faces; num_face++)
    {
      resu(num_face,0) += (P(num_face)-K_eps_Bas_Re(num_face,1)-D(num_face))*vol_ent(num_face);
      if (K_eps_Bas_Re(num_face,0) >= 1.e-20)
        resu(num_face,1) += ((C1*F1(num_face)*P(num_face)- C2*F2(num_face)*K_eps_Bas_Re(num_face,1))*K_eps_Bas_Re(num_face,1)/K_eps_Bas_Re(num_face,0)+E(num_face))*vol_ent(num_face);

      if ( (G(num_face)>0) && (K_eps_Bas_Re(num_face,1) >= 1.e-20) )
        {
          resu(num_face,0) += G(num_face)*vol_ent(num_face);
          resu(num_face,1) += C1*F1(num_face)*G(num_face)*vol_ent(num_face)*K_eps_Bas_Re(num_face,1)/K_eps_Bas_Re(num_face,0);
        }
    }

  return resu;
}
