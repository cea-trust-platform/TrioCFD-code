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
// File      : Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_Elem.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_Elem.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <Convection_Diffusion_Concentration.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <Zone_Cl_VDF.h>
#include <Champ_Face.h>
#include <TRUSTTrav.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_Elem,"Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_P0_VDF",Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem);

Sortie& Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_Elem::readOn(Entree& s) { return s ; }

void Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem::verifier_milieu_anisotherme_concen(pb,que_suis_je());
  Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem::associer_pb(pb);
  Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem::associer_pb_anisotherme_concen(pb);
}

// TODO : FIXME : a factoriser avec Source_Transport_K_Eps_Bas_Reynolds_anisotherme_VDF_Elem::ajouter
void Source_Transport_K_Eps_Bas_Reynolds_aniso_therm_concen_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& resu, const tabs_t& semi_impl) const
{
  const Zone_Cl_dis& zcl=eq_hydraulique->zone_Cl_dis();
  const Zone_Cl_dis& zcl_keps=eqn_keps_bas_re->zone_Cl_dis();
  const Zone_dis& zone_dis_keps =eqn_keps_bas_re ->zone_dis();
  const Zone_VDF& zone_VDF = ref_cast(Zone_VDF,eq_hydraulique->zone_dis().valeur());
  const Zone_Cl_VDF& zone_Cl_VDF = ref_cast(Zone_Cl_VDF,zcl.valeur());
  const Zone_Cl_VDF& zcl_VDF_th = ref_cast(Zone_Cl_VDF,eq_thermique->zone_Cl_dis().valeur());
  const Zone_Cl_VDF& zcl_VDF_co = ref_cast(Zone_Cl_VDF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& K_eps_Bas_Re = eqn_keps_bas_re->inconnue().valeurs();
  const DoubleTab& temper = eq_thermique->inconnue().valeurs(), &concen = eq_concentration->inconnue().valeurs(), &vit = eq_hydraulique->inconnue().valeurs();
  const DoubleTab& visco_turb = eqn_keps_bas_re->modele_turbulence().viscosite_turbulente().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  const DoubleVect& g = gravite->valeurs(), &volumes = zone_VDF.volumes(), &porosite_vol = zone_VDF.porosite_elem();
  const Champ_Don& ch_beta_temper = beta_t.valeur();
  const Champ_Uniforme& ch_beta_concen = ref_cast(Champ_Uniforme, beta_c->valeur());
  const Fluide_base& fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());
  const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds,eqn_keps_bas_re->modele_turbulence());
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc = ref_cast(Modele_Fonc_Bas_Reynolds,mod_turb.associe_modele_fonction());
  const int nb_elem = zone_VDF.nb_elem(), nb_elem_tot = zone_VDF.nb_elem_tot(), nb_consti = eq_concentration->nb_constituants();

  DoubleTrav P(nb_elem_tot), G_t(nb_elem_tot), G_c(nb_elem_tot), D(nb_elem_tot), E(nb_elem_tot), F1(nb_elem_tot), F2(nb_elem_tot);
  mon_modele_fonc.Calcul_D(D,zone_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin);
  mon_modele_fonc.Calcul_E(E,zone_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin,visco_turb);
  mon_modele_fonc.Calcul_F2(F2,D,zone_dis_keps,K_eps_Bas_Re,ch_visco_cin);

  if (axi) calculer_terme_production_K_Axi(zone_VDF,vitesse,P,K_eps_Bas_Re,visco_turb);
  else calculer_terme_production_K(zone_VDF,zone_Cl_VDF,P,K_eps_Bas_Re,vit,vitesse,visco_turb);

  const DoubleTab& tab_beta_t = ch_beta_temper.valeurs();

  if (sub_type(Champ_Uniforme,ch_beta_temper.valeur())) calculer_terme_destruction_K(zone_VDF,zcl_VDF_th,G_t,temper,alpha_turb,tab_beta_t(0,0),g);
  else calculer_terme_destruction_K(zone_VDF,zcl_VDF_th,G_t,temper,alpha_turb,tab_beta_t,g);

  if (nb_consti == 1) calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G_c,concen,alpha_turb,ch_beta_concen(0,0),g);
  else
    {
      const DoubleVect& d_beta_c = ch_beta_concen.valeurs();
      calculer_terme_destruction_K(zone_VDF,zcl_VDF_co,G_c,concen,alpha_turb,d_beta_c,g,nb_consti);
    }

  for (int elem = 0; elem < nb_elem; elem++)
    {
      resu(elem,0) += (P(elem)-K_eps_Bas_Re(elem,1)-D(elem))*volumes(elem)*porosite_vol(elem);
      if (K_eps_Bas_Re(elem,0) >= 1.e-20)
        resu(elem,1) += (C1*F1(elem)*P(elem)- C2*F2(elem)*K_eps_Bas_Re(elem,1))*volumes(elem)*porosite_vol(elem) *K_eps_Bas_Re(elem,0)/K_eps_Bas_Re(elem,1) + E(elem)*volumes(elem)*porosite_vol(elem);

      if ( ((G_t(elem)+G_c(elem))>0 ) && (K_eps_Bas_Re(elem,1) >= 1.e-20) )
        {
          resu(elem,0) += (G_t(elem)+G_c(elem))*volumes(elem)*porosite_vol(elem);
          resu(elem,1) += C1*F1(elem)*(G_t(elem)+G_c(elem))*volumes(elem)*porosite_vol(elem)*K_eps_Bas_Re(elem,0)/K_eps_Bas_Re(elem,1)+E(elem)*volumes(elem)*porosite_vol(elem);
        }
    }
}
