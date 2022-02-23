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
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Source_Transport_VDF_Elem_base.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Fluide_Dilatable_base.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Fluide_base.h>
#include <Zone_Cl_VDF.h>
#include <Constituant.h>
#include <Champ_Face.h>
#include <Zone_VDF.h>
#include <Param.h>
#include <Debog.h>

Implemente_base_sans_constructeur( Source_Transport_VDF_Elem_base, "Source_Transport_VDF_Elem_base", Source_base);

Sortie& Source_Transport_VDF_Elem_base::printOn( Sortie& os ) const { return os << que_suis_je(); }

Entree& Source_Transport_VDF_Elem_base::readOn( Entree& is )
{
  Param param(que_suis_je());
  param.ajouter("C1_eps", &C1);
  param.ajouter("C2_eps", &C2);
  param.lire_avec_accolades(is);
  Cerr << "C1_eps = " << C1 << finl;
  Cerr << "C2_eps = " << C2 << finl;
  return is;
}

Entree& Source_Transport_VDF_Elem_base::readOn_anisotherme(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("C1_eps", &C1);
  param.ajouter("C2_eps", &C2);
  param.ajouter("C3_eps", &C3);
  param.lire_avec_accolades(is);
  Cerr << "C1_eps = " << C1 << finl;
  Cerr << "C2_eps = " << C2 << finl;
  Cerr << "C3_eps = " << C3 << finl;
  return is ;
}

Entree& Source_Transport_VDF_Elem_base::readOn_concen(Entree& is) { return readOn_anisotherme(is); }

void Source_Transport_VDF_Elem_base::verifier_pb_keps(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Hydraulique_Turbulent,pb) && !sub_type(Pb_Thermohydraulique_Turbulent_QC,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_VDF_Elem_base::verifier_pb_keps_anisotherme(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Thermohydraulique_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_VDF_Elem_base::verifier_pb_keps_concen(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Hydraulique_Concentration_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_VDF_Elem_base::verifier_milieu_anisotherme(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(1).milieu();
  if (pb.nombre_d_equations()<2) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_VDF_Elem_base::verifier_milieu_concen(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(0).milieu(); // XXX : Attention pas eq 1 car Constituant derive pas de Fluide_base !
  if (pb.nombre_d_equations()<2) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_VDF_Elem_base::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis&  zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF,zone_Cl_dis.valeur());
}

void Source_Transport_VDF_Elem_base::associer_pb(const Probleme_base& pb)
{
  eq_hydraulique = pb.equation(0);
}

void Source_Transport_VDF_Elem_base::associer_pb_anisotherme(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(1).milieu());
  beta_t = fluide.beta_t();
  gravite = fluide.gravite();
  eq_thermique = ref_cast(Convection_Diffusion_Temperature,pb.equation(1));
}

void Source_Transport_VDF_Elem_base::associer_pb_concen(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(0).milieu()); // XXX : Attention pas eq 1 car Constituant derive pas de Fluide_base !
  if (!fluide.beta_c().non_nul())
    {
      Cerr << "You forgot to define beta_co field in the fluid. It is mandatory when using the K-Eps model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      Process::exit();
    }
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
  eq_concentration = ref_cast(Convection_Diffusion_Concentration,pb.equation(1));
}

DoubleTab& Source_Transport_VDF_Elem_base::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

DoubleTab& Source_Transport_VDF_Elem_base::ajouter_keps(DoubleTab& resu) const
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const Champ_Face& ch_vit = ref_cast(Champ_Face,eq_hydraulique->inconnue().valeur());

  DoubleVect P; // Ajout d'un espace virtuel au tableu P
  zone_VDF.zone().creer_tableau_elements(P);
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
  const Zone_Cl_VDF& zcl_VDF_th = ref_cast(Zone_Cl_VDF,eq_thermique->zone_Cl_dis().valeur());
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  const DoubleTab& g = gravite->valeurs(), &tab_beta = beta_t->valeurs();
  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  const DoubleVect& volumes = la_zone_VDF->volumes(), &porosite_vol = la_zone_VDF->porosite_elem();

  // Ajout d'un espace virtuel au tableau G
  DoubleVect G;
  la_zone_VDF->zone().creer_tableau_elements(G);


  if (sub_type(Champ_Uniforme,beta_t->valeur())) calculer_terme_destruction_K(la_zone_VDF.valeur(),zcl_VDF_th,G,scalaire,alpha_turb,tab_beta(0,0),g);
  else calculer_terme_destruction_K(la_zone_VDF.valeur(),zcl_VDF_th,G,scalaire,alpha_turb,tab_beta,g);

  fill_resu_anisotherme(G,volumes,porosite_vol,resu); // voir les classes filles
  return resu;
}

DoubleTab& Source_Transport_VDF_Elem_base::ajouter_concen(DoubleTab& resu) const
{
  const Zone_Cl_VDF& zcl_VDF_co = ref_cast(Zone_Cl_VDF,eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& diffu_turb  = le_modele_scalaire.conductivite_turbulente().valeurs();
  const Champ_Uniforme& ch_beta_concen = ref_cast(Champ_Uniforme, beta_c->valeur());
  const DoubleVect& g = gravite->valeurs(), &volumes = la_zone_VDF->volumes(), &porosite_vol = la_zone_VDF->porosite_elem();
  const int nb_consti = eq_concentration->constituant().nb_constituants();

  // Ajout d'un espace virtuel au tableau G
  DoubleVect G;
  la_zone_VDF->zone().creer_tableau_elements(G);

  if (nb_consti == 1) calculer_terme_destruction_K(la_zone_VDF.valeur(),zcl_VDF_co,G, concen,diffu_turb,ch_beta_concen(0,0),g);
  else
    {
      const DoubleVect& d_beta_c = ch_beta_concen.valeurs();
      calculer_terme_destruction_K(la_zone_VDF.valeur(),zcl_VDF_co,G, concen,diffu_turb,d_beta_c,g, nb_consti);
    }

  fill_resu_concen(G,volumes,porosite_vol,resu); // voir les classes filles
  return resu;
}
