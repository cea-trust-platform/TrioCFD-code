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
// File      : Source_Transport_VDF_Elem_base.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Source_Transport_VDF_Elem_base.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <Domaine_Cl_VDF.h>
#include <Constituant.h>
#include <Champ_Face_VDF.h>
#include <Debog.h>

Implemente_base_sans_constructeur( Source_Transport_VDF_Elem_base, "Source_Transport_VDF_Elem_base", Source_base);

Sortie& Source_Transport_VDF_Elem_base::printOn( Sortie& os ) const { return os << que_suis_je(); }
Entree& Source_Transport_VDF_Elem_base::readOn( Entree& is ) { return Source_Transport_proto::readOn_proto(is,que_suis_je()); }

void Source_Transport_VDF_Elem_base::associer_domaines(const Domaine_dis& domaine_dis, const Domaine_Cl_dis&  domaine_Cl_dis)
{
  le_dom_VDF = ref_cast(Domaine_VDF, domaine_dis.valeur());
  le_dom_Cl_VDF = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis.valeur());
}

void Source_Transport_VDF_Elem_base::associer_pb(const Probleme_base& pb) { Source_Transport_proto::associer_pb_proto(pb); }

DoubleTab& Source_Transport_VDF_Elem_base::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

DoubleTab& Source_Transport_VDF_Elem_base::ajouter_keps(DoubleTab& resu) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const Champ_Face_VDF& ch_vit = ref_cast(Champ_Face_VDF,eq_hydraulique->inconnue().valeur());

  DoubleVect P; // Ajout d'un espace virtuel au tableu P
  domaine_VDF.domaine().creer_tableau_elements(P);
  calculer_terme_production(ch_vit,visco_turb,vit,P); // voir les classes filles

  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc = get_modele_fonc_bas_reyn();
  const int is_modele_fonc=(mon_modele_fonc.non_nul());

  if (is_modele_fonc)
    {
      DoubleTab& D=ref_cast_non_const(DoubleTab,mon_modele_fonc->get_champ("D").valeurs());
      DoubleTab& E=ref_cast_non_const(DoubleTab,mon_modele_fonc->get_champ("E").valeurs());
      DoubleTab& F1=ref_cast_non_const(DoubleTab,mon_modele_fonc->get_champ("F1").valeurs());
      DoubleTab& F2=ref_cast_non_const(DoubleTab,mon_modele_fonc->get_champ("F2").valeurs());
      const Fluide_base& fluide=ref_cast(Fluide_base,eq_hydraulique->milieu());
      const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
      calcul_D_E(vit,visco_turb,ch_visco_cin,D,E); // voir les classes filles
      D.echange_espace_virtuel();
      E.echange_espace_virtuel();
      const Champ_base& ch_visco_cin_ou_dyn =ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();
      DoubleTab P_tab;
      P_tab.ref(P);
      calcul_F1_F2(ch_visco_cin_ou_dyn,P_tab,D,F1,F2); // voir les classes filles

      Debog::verifier("D",D);
      Debog::verifier("E",E);
      Debog::verifier("F1",F1);
      Debog::verifier("F2",F2);
      Debog::verifier("avt",resu);

      fill_resu_bas_rey(P,D,E,F1,F2,resu); // voir les classes filles
      Debog::verifier("ap",resu);
    }
  else fill_resu(P,resu);

  resu.echange_espace_virtuel();
  return resu;
}

DoubleTab& Source_Transport_VDF_Elem_base::ajouter_anisotherme(DoubleTab& resu) const
{
  // on ajoute directement G
  const Domaine_Cl_VDF& zcl_VDF_th = ref_cast(Domaine_Cl_VDF,eq_thermique->domaine_Cl_dis().valeur());
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& g = gravite->valeurs(), &tab_beta = beta_t->valeurs();
  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  const DoubleVect& volumes = le_dom_VDF->volumes(), &porosite_vol = le_dom_Cl_VDF->equation().milieu().porosite_elem();

  // Ajout d'un espace virtuel au tableau G
  DoubleVect G;
  le_dom_VDF->domaine().creer_tableau_elements(G);


  if (sub_type(Champ_Uniforme,beta_t->valeur())) calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_th,G,scalaire,alpha_turb,tab_beta(0,0),g);
  else calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_th,G,scalaire,alpha_turb,tab_beta,g);

  fill_resu_anisotherme(G,volumes,porosite_vol,resu); // voir les classes filles
  return resu;
}

DoubleTab& Source_Transport_VDF_Elem_base::ajouter_concen(DoubleTab& resu) const
{
  const Domaine_Cl_VDF& zcl_VDF_co = ref_cast(Domaine_Cl_VDF,eq_concentration->domaine_Cl_dis().valeur());
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& diffu_turb = le_modele_scalaire.conductivite_turbulente().valeurs();
//  const DoubleTab& diffu_turb = le_modele_scalaire.diffusivite_turbulente().valeurs(); // XXX : realisable utilise ca ???? a voir
  const Champ_Uniforme& ch_beta_concen = ref_cast(Champ_Uniforme, beta_c->valeur());
  const DoubleVect& g = gravite->valeurs(), &volumes = le_dom_VDF->volumes(), &porosite_vol = le_dom_Cl_VDF->equation().milieu().porosite_elem();
  const int nb_consti = eq_concentration->constituant().nb_constituants();

  // Ajout d'un espace virtuel au tableau G
  DoubleVect G;
  le_dom_VDF->domaine().creer_tableau_elements(G);

  if (nb_consti == 1) calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_co,G,concen,diffu_turb,ch_beta_concen(0,0),g);
  else
    {
      const DoubleVect& d_beta_c = ch_beta_concen.valeurs();
      calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_co,G,concen,diffu_turb,d_beta_c,g,nb_consti);
    }

  fill_resu_concen(G,volumes,porosite_vol,resu); // voir les classes filles
  return resu;
}

// TODO : FIXME : on peut factoriser avec les 2 methodes ajouter_anisotherme et ajouter_concen
DoubleTab& Source_Transport_VDF_Elem_base::ajouter_anisotherme_concen(DoubleTab& resu) const
{
  const Domaine_Cl_VDF& zcl_VDF_th = ref_cast(Domaine_Cl_VDF,eq_thermique->domaine_Cl_dis().valeur());
  const Domaine_Cl_VDF& zcl_VDF_co = ref_cast(Domaine_Cl_VDF,eq_concentration->domaine_Cl_dis().valeur());
  const DoubleTab& temper = eq_thermique->inconnue().valeurs(), &concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const Modele_turbulence_scal_base& le_modele_scal_co = ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());

  // XXX : Elie Saikali : vaut mieux utiliser diffusivite_turbulente au lie de faire ca ....
  // voila dans Source_Transport_Eps_Realisable_aniso_therm_concen_VDF_Elem
  // const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  DoubleTab alpha_turb(le_modele_scalaire.conductivite_turbulente().valeurs()); // on veut pas modifier la ref !
  double rhocp = eq_thermique->milieu().capacite_calorifique().valeurs()(0, 0) * eq_thermique->milieu().masse_volumique().valeurs()(0, 0);
  alpha_turb /= rhocp;

  const Champ_Don& ch_beta_temper = beta_t.valeur();
  const DoubleTab& diffu_turb = le_modele_scal_co.conductivite_turbulente().valeurs(), &tab_beta_t = ch_beta_temper.valeurs();
//  const DoubleTab& diffu_turb = le_modele_scal_co.diffusivite_turbulente().valeurs(); // XXX : realisable utilise ca ???? a voir
  const Champ_Uniforme& ch_beta_concen = ref_cast(Champ_Uniforme, beta_c->valeur());
  const DoubleVect& volumes = le_dom_VDF->volumes(), &porosite_vol = le_dom_Cl_VDF->equation().milieu().porosite_elem(), &g = gravite->valeurs();
  const int nb_consti = eq_concentration->constituant().nb_constituants();

  // Ajout d'un espace virtuel au tableaux Gt et Gc
  DoubleVect G_t, G_c;
  le_dom_VDF->domaine().creer_tableau_elements(G_t);
  le_dom_VDF->domaine().creer_tableau_elements(G_c);

  if (sub_type(Champ_Uniforme,ch_beta_temper.valeur())) calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_th,G_t,temper,alpha_turb,tab_beta_t(0,0),g);
  else calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_th,G_t,temper,alpha_turb,tab_beta_t,g);

  if (nb_consti == 1) calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_co,G_c,concen,diffu_turb,ch_beta_concen(0,0),g);
  else
    {
      const DoubleVect& d_beta_c = ch_beta_concen.valeurs();
      calculer_terme_destruction_K(le_dom_VDF.valeur(),zcl_VDF_co,G_c,concen,diffu_turb,d_beta_c,g,nb_consti);
    }

  fill_resu_anisotherme_concen(G_t,G_c,volumes,porosite_vol,resu); // voir les classes filles
  return resu;
}
