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
// File      : Source_Transport_VEF_Face_base.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS/Base
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Source_Transport_VEF_Face_base.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Fluide_Dilatable_base.h>
#include <Champ_Uniforme.h>
#include <Constituant.h>
#include <Fluide_base.h>
#include <Zone_Cl_VEF.h>
#include <Champ_P1NC.h>
#include <Zone_VEF.h>
#include <Param.h>

Implemente_base_sans_constructeur( Source_Transport_VEF_Face_base, "Source_Transport_VEF_Face_base", Source_base );

Sortie& Source_Transport_VEF_Face_base::printOn( Sortie& os ) const { return os << que_suis_je(); }

Entree& Source_Transport_VEF_Face_base::readOn( Entree& is )
{
  Param param(que_suis_je());
  param.ajouter("C1_eps", &C1);
  param.ajouter("C2_eps", &C2);
  param.lire_avec_accolades(is);
  Cerr << "C1_eps = " << C1 << finl;
  Cerr << "C2_eps = " << C2 << finl;
  return is;
}

Entree& Source_Transport_VEF_Face_base::readOn_nothing(Entree& is)
{
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}

Entree& Source_Transport_VEF_Face_base::readOn_anisotherme(Entree& is)
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

Entree& Source_Transport_VEF_Face_base::readOn_concen(Entree& is) { return readOn_anisotherme(is); }

Entree& Source_Transport_VEF_Face_base::readOn_anisotherme_concen(Entree& is) { return readOn_anisotherme(is); }

void Source_Transport_VEF_Face_base::verifier_pb_keps(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Hydraulique_Turbulent,pb) && !sub_type(Pb_Thermohydraulique_Turbulent_QC,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_VEF_Face_base::verifier_pb_keps_anisotherme(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Thermohydraulique_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_VEF_Face_base::verifier_pb_keps_concen(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Hydraulique_Concentration_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_VEF_Face_base::verifier_pb_keps_anisotherme_concen(const Probleme_base& pb, const Nom& nom)
{
  if (!sub_type(Pb_Thermohydraulique_Concentration_Turbulent,pb)) error_keps(nom,pb.que_suis_je());
}

void Source_Transport_VEF_Face_base::verifier_milieu_anisotherme(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(1).milieu(); // eq thermique
  if (pb.nombre_d_equations()<2) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_VEF_Face_base::verifier_milieu_concen(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(0).milieu(); // XXX : Attention pas eq 1 car Constituant derive pas de Fluide_base ! donc eq hydro
  if (pb.nombre_d_equations()<2) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_VEF_Face_base::verifier_milieu_anisotherme_concen(const Probleme_base& pb, const Nom& nom)
{
  const Milieu_base& milieu = pb.equation(1).milieu(); // eq thermique
  if (pb.nombre_d_equations()<3) error_keps(nom,pb.que_suis_je());
  if (sub_type(Fluide_Dilatable_base,ref_cast(Fluide_base,milieu))) error_keps(nom,milieu.que_suis_je());
}

void Source_Transport_VEF_Face_base::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF, zone_dis.valeur());
  la_zone_Cl_VEF = ref_cast(Zone_Cl_VEF, zone_Cl_dis.valeur());
}

void Source_Transport_VEF_Face_base::associer_pb(const Probleme_base& pb)
{
  eq_hydraulique = pb.equation(0);
}

void Source_Transport_VEF_Face_base::associer_pb_anisotherme(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(1).milieu());
  beta_t = fluide.beta_t();
  gravite = fluide.gravite();
  eq_thermique = ref_cast(Convection_Diffusion_Temperature,pb.equation(1));
}

void Source_Transport_VEF_Face_base::associer_pb_concen(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(0).milieu()); // XXX : Attention pas eq 1 car Constituant derive pas de Fluide_base !
  verifier_beta_concen(fluide);
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
  eq_concentration = ref_cast(Convection_Diffusion_Concentration,pb.equation(1));
}

void Source_Transport_VEF_Face_base::associer_pb_anisotherme_concen(const Probleme_base& pb)
{
  const Fluide_base& fluide = ref_cast(Fluide_base,pb.equation(1).milieu()); // a partir de l'eq thermique
  verifier_beta_concen(fluide);
  beta_t = fluide.beta_t();
  beta_c = fluide.beta_c();
  gravite = fluide.gravite();
  eq_thermique = ref_cast(Convection_Diffusion_Temperature,pb.equation(1));
  eq_concentration = ref_cast(Convection_Diffusion_Concentration,pb.equation(2));
}

void Source_Transport_VEF_Face_base::verifier_beta_concen(const Fluide_base& fluide)
{
  if (!fluide.beta_c().non_nul())
    {
      Cerr << "You forgot to define beta_co field in the fluid. It is mandatory when using the K-Eps model (buoyancy effects)." << finl;
      Cerr << "If you don't want buoyancy effects, then specify: beta_co champ_uniforme 1 0." << finl;
      Process::exit();
    }
}

DoubleTab& Source_Transport_VEF_Face_base::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

DoubleTab& Source_Transport_VEF_Face_base::ajouter_keps(DoubleTab& resu) const
{
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF, eq_hydraulique->zone_Cl_dis().valeur());
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& volumes_entrelaces = la_zone_VEF->volumes_entrelaces();
  const DoubleTab& tab = get_cisaillement_paroi(); // voir les classes filles
  const int nb_faces_ = la_zone_VEF->nb_faces();
  DoubleTrav P(nb_faces_);
  const DoubleTab& K = get_K_pour_production(); // voir les classes filles
  calculer_terme_production_K(la_zone_VEF.valeur(), zone_Cl_VEF, P, K, vit, visco_turb);

  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc = get_modele_fonc_bas_reyn(); // voir les classes filles
  const int is_modele_fonc = (mon_modele_fonc.non_nul());
  if (is_modele_fonc)
    {
      DoubleTab& D = ref_cast_non_const(DoubleTab, mon_modele_fonc->get_champ("D").valeurs());
      DoubleTab& E = ref_cast_non_const(DoubleTab, mon_modele_fonc->get_champ("E").valeurs());
      DoubleTab& F1 = ref_cast_non_const(DoubleTab, mon_modele_fonc->get_champ("F1").valeurs());
      DoubleTab& F2 = ref_cast_non_const(DoubleTab, mon_modele_fonc->get_champ("F2").valeurs());

      const Fluide_base& fluide = ref_cast(Fluide_base, eq_hydraulique->milieu());
      const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
      const Champ_base& ch_visco_cin_ou_dyn = ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();
      calcul_tabs_bas_reyn(P, vit, visco_turb, ch_visco_cin, ch_visco_cin_ou_dyn, D, E, F1, F2); // voir les classes filles

      // Pour modele EASM
      const int nb_elem_tot = la_zone_VEF->nb_elem_tot();

      int is_Reynolds_stress_isotrope = mon_modele_fonc.Calcul_is_Reynolds_stress_isotrope();
      if (is_Reynolds_stress_isotrope == 0)
        {
          Cerr << "On utilise une diffusion turbulente non linaire dans le terme source P" << finl;
          const DoubleTab& visco_scal = ch_visco_cin->valeurs();
          DoubleTab visco_tab(nb_elem_tot);
          assert(sub_type(Champ_Uniforme,ch_visco_cin.valeur()));
          visco_tab = visco_scal(0, 0);
          DoubleTab gradient_elem(nb_elem_tot, Objet_U::dimension, Objet_U::dimension);
          gradient_elem = 0.;
          Champ_P1NC::calcul_gradient(vit, gradient_elem, zone_Cl_VEF);
          /*Paroi*/
          const Nom lp = get_type_paroi(); // voir les classes filles
          if (lp != "negligeable_VEF")
            if (mon_equation->schema_temps().nb_pas_dt() > 0)
              Champ_P1NC::calcul_duidxj_paroi(gradient_elem, visco_tab, visco_turb, tab, zone_Cl_VEF);

          gradient_elem.echange_espace_virtuel();
          DoubleTab Re(gradient_elem);

          calcul_tenseur_reyn(visco_turb, gradient_elem, Re); // voir les classes filles
          Re.echange_espace_virtuel();

          calculer_terme_production_K_EASM(la_zone_VEF.valeur(), zone_Cl_VEF, P, K, gradient_elem, visco_turb, Re);
        }  // Fin pour modele EASM

      fill_resu_bas_rey(volumes_entrelaces, P, D, E, F1, F2, resu); // voir les classes filles
    }
  else fill_resu(volumes_entrelaces, P, resu); // voir les classes filles

  return resu;
}

DoubleTab& Source_Transport_VEF_Face_base::ajouter_anisotherme(DoubleTab& resu) const
{
  // on ajoute directement G
  const Zone_Cl_VEF& zcl_VEF_th = ref_cast(Zone_Cl_VEF,eq_thermique->zone_Cl_dis().valeur());
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  DoubleTab alpha_turb(le_modele_scalaire.diffusivite_turbulente().valeurs());
  const DoubleTab& g = gravite->valeurs();
  const Champ_Don& ch_beta = beta_t.valeur();
  const DoubleVect& volumes_entrelaces = la_zone_VEF->volumes_entrelaces();
  const int nb_face = la_zone_VEF->nb_faces();
  DoubleTrav G(nb_face);

  // C'est l'objet de type zone_Cl_dis de l'equation thermique qui est utilise dans le calcul de G
  calculer_terme_destruction_K_gen(la_zone_VEF.valeur(),zcl_VEF_th,G,scalaire,alpha_turb,ch_beta,g,0);

  fill_resu_anisotherme(G,volumes_entrelaces,resu); // voir les classes filles
  return resu;
}

DoubleTab& Source_Transport_VEF_Face_base::ajouter_concen(DoubleTab& resu) const
{
  // on ajoute directement G
  const Zone_Cl_VEF& zcl_VEF_co = ref_cast(Zone_Cl_VEF, eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base, eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& lambda_turb = le_modele_scalaire.conductivite_turbulente().valeurs();
//  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs(); // XXX : realisable utilise ca ???? a voir
  const DoubleVect& g = gravite->valeurs();
  const Champ_Don& ch_beta_concen = beta_c.valeur();
  const DoubleVect& volumes_entrelaces = la_zone_VEF->volumes_entrelaces();
  const int nb_face = la_zone_VEF->nb_faces(), nb_consti = eq_concentration->constituant().nb_constituants();
  DoubleTrav G(nb_face);

  calculer_terme_destruction_K_gen(la_zone_VEF.valeur(), zcl_VEF_co, G, concen, lambda_turb, ch_beta_concen, g, nb_consti);

  fill_resu_concen(G,volumes_entrelaces,resu); // voir les classes filles
  return resu;
}

// TODO : FIXME : on peut factoriser avec les 2 methodes ajouter_anisotherme et ajouter_concen
DoubleTab& Source_Transport_VEF_Face_base::ajouter_anisotherme_concen(DoubleTab& resu) const
{
  // on ajoute directement G
  const Zone_Cl_VEF& zcl_VEF_th = ref_cast(Zone_Cl_VEF, eq_thermique->zone_Cl_dis().valeur());
  const Zone_Cl_VEF& zcl_VEF_co = ref_cast(Zone_Cl_VEF, eq_concentration->zone_Cl_dis().valeur());
  const DoubleTab& temper = eq_thermique->inconnue().valeurs(), &concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base, eq_thermique->get_modele(TURBULENCE).valeur());

  // XXX : Elie Saikali : vaut mieux utiliser diffusivite_turbulente au lie de faire ca ....
  // voila dans Source_Transport_K_Eps_Realisable_aniso_therm_concen_VEF_Face
  // const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  DoubleTab alpha_turb(le_modele_scalaire.conductivite_turbulente().valeurs());
  double rhocp = eq_thermique->milieu().capacite_calorifique().valeurs()(0, 0) * eq_thermique->milieu().masse_volumique().valeurs()(0, 0);
  alpha_turb /= rhocp;
  const DoubleVect& g = gravite->valeurs(), &volumes_entrelaces = la_zone_VEF->volumes_entrelaces();
  const Champ_Don& ch_beta_temper = beta_t.valeur(), &ch_beta_concen = beta_c.valeur();
  const int nb_face = la_zone_VEF->nb_faces(), nb_consti = eq_concentration->constituant().nb_constituants();
  DoubleTrav P(nb_face), G_t(nb_face), G_c(nb_face);

  calculer_terme_destruction_K_gen(la_zone_VEF.valeur(), zcl_VEF_th, G_t, temper, alpha_turb, ch_beta_temper, g, 0);
  calculer_terme_destruction_K_gen(la_zone_VEF.valeur(), zcl_VEF_co, G_c, concen, alpha_turb, ch_beta_concen, g, nb_consti);

  fill_resu_anisotherme_concen(G_t, G_c, volumes_entrelaces,resu); // voir les classes filles
  return resu;
}

