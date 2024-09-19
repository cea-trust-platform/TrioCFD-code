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
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Source_Transport_VEF_Face_base.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Champ_Uniforme.h>
#include <Constituant.h>
#include <Fluide_base.h>
#include <Domaine_Cl_VEF.h>
#include <Champ_P1NC.h>
#include <Domaine_VEF.h>

Implemente_base_sans_constructeur( Source_Transport_VEF_Face_base, "Source_Transport_VEF_Face_base", Source_base );

Sortie& Source_Transport_VEF_Face_base::printOn( Sortie& os ) const { return os << que_suis_je(); }
Entree& Source_Transport_VEF_Face_base::readOn( Entree& is ) { return Source_Transport_proto::readOn_proto(is,que_suis_je()); }

void Source_Transport_VEF_Face_base::associer_domaines(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_VEF = ref_cast(Domaine_VEF, domaine_dis);
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis);
}

void Source_Transport_VEF_Face_base::associer_pb(const Probleme_base& pb) { Source_Transport_proto::associer_pb_proto(pb); }

DoubleTab& Source_Transport_VEF_Face_base::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

DoubleTab& Source_Transport_VEF_Face_base::ajouter_keps(DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF, eq_hydraulique->domaine_Cl_dis());
  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();
  const DoubleTab& tab = get_cisaillement_paroi(); // voir les classes filles
  const int nb_faces_ = le_dom_VEF->nb_faces();
  DoubleTrav P(nb_faces_);
  const DoubleTab& K = get_K_pour_production(); // voir les classes filles
  calculer_terme_production_K(le_dom_VEF.valeur(), domaine_Cl_VEF, P, K, vit, visco_turb, _interpolation_viscosite_turbulente, _coefficient_limiteur);

  const OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& mon_modele_fonc = get_modele_fonc_bas_reyn(); // voir les classes filles
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
      const int nb_elem_tot = le_dom_VEF->nb_elem_tot();

      int is_Reynolds_stress_isotrope = mon_modele_fonc->Calcul_is_Reynolds_stress_isotrope();
      if (is_Reynolds_stress_isotrope == 0)
        {
          Cerr << "On utilise une diffusion turbulente non linaire dans le terme source P" << finl;
          const DoubleTab& visco_scal = ch_visco_cin->valeurs();
          DoubleTab visco_tab(nb_elem_tot);
          assert(sub_type(Champ_Uniforme,ch_visco_cin.valeur()));
          visco_tab = visco_scal(0, 0);
          DoubleTab gradient_elem(nb_elem_tot, Objet_U::dimension, Objet_U::dimension);
          gradient_elem = 0.;
          Champ_P1NC::calcul_gradient(vit, gradient_elem, domaine_Cl_VEF);
          /*Paroi*/
          const Nom lp = get_type_paroi(); // voir les classes filles
          if (lp != "negligeable_VEF")
            if (mon_equation->schema_temps().nb_pas_dt() > 0)
              Champ_P1NC::calcul_duidxj_paroi(gradient_elem, visco_tab, visco_turb, tab, domaine_Cl_VEF);

          gradient_elem.echange_espace_virtuel();
          DoubleTab Re(gradient_elem);

          calcul_tenseur_reyn(visco_turb, gradient_elem, Re); // voir les classes filles
          Re.echange_espace_virtuel();

          calculer_terme_production_K_EASM(le_dom_VEF.valeur(), domaine_Cl_VEF, P, K, gradient_elem, visco_turb, Re, _interpolation_viscosite_turbulente, _coefficient_limiteur);
        }  // Fin pour modele EASM

      fill_resu_bas_rey(volumes_entrelaces, P, D, E, F1, F2, resu); // voir les classes filles
    }
  else fill_resu(volumes_entrelaces, P, resu); // voir les classes filles

  return resu;
}

DoubleTab& Source_Transport_VEF_Face_base::ajouter_anisotherme(DoubleTab& resu) const
{
  // on ajoute directement G
  const Domaine_Cl_VEF& zcl_VEF_th = ref_cast(Domaine_Cl_VEF,eq_thermique->domaine_Cl_dis());
  const DoubleTab& scalaire = eq_thermique->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_thermique->get_modele(TURBULENCE).valeur());
  DoubleTab alpha_turb(le_modele_scalaire.diffusivite_turbulente()->valeurs());
  const DoubleTab& g = gravite->valeurs();
  const Champ_Don& ch_beta = beta_t.valeur();
  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();
  const int nb_face = le_dom_VEF->nb_faces();
  DoubleTrav G(nb_face);

  // C'est l'objet de type domaine_Cl_dis de l'equation thermique qui est utilise dans le calcul de G
  calculer_terme_destruction_K_gen(le_dom_VEF.valeur(),zcl_VEF_th,G,scalaire,alpha_turb,ch_beta,g,0);

  fill_resu_anisotherme(G,volumes_entrelaces,resu); // voir les classes filles
  return resu;
}

DoubleTab& Source_Transport_VEF_Face_base::ajouter_concen(DoubleTab& resu) const
{
  // on ajoute directement G
  const Domaine_Cl_VEF& zcl_VEF_co = ref_cast(Domaine_Cl_VEF, eq_concentration->domaine_Cl_dis());
  const DoubleTab& concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base, eq_concentration->get_modele(TURBULENCE).valeur());
  const DoubleTab& lambda_turb = le_modele_scalaire.conductivite_turbulente()->valeurs();
//  const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs(); // XXX : realisable utilise ca ???? a voir
  const DoubleVect& g = gravite->valeurs();
  const Champ_Don& ch_beta_concen = beta_c.valeur();
  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();
  const int nb_face = le_dom_VEF->nb_faces(), nb_consti = eq_concentration->constituant().nb_constituants();
  DoubleTrav G(nb_face);

  calculer_terme_destruction_K_gen(le_dom_VEF.valeur(), zcl_VEF_co, G, concen, lambda_turb, ch_beta_concen, g, nb_consti);

  fill_resu_concen(G,volumes_entrelaces,resu); // voir les classes filles
  return resu;
}

// TODO : FIXME : on peut factoriser avec les 2 methodes ajouter_anisotherme et ajouter_concen
DoubleTab& Source_Transport_VEF_Face_base::ajouter_anisotherme_concen(DoubleTab& resu) const
{
  // on ajoute directement G
  const Domaine_Cl_VEF& zcl_VEF_th = ref_cast(Domaine_Cl_VEF, eq_thermique->domaine_Cl_dis());
  const Domaine_Cl_VEF& zcl_VEF_co = ref_cast(Domaine_Cl_VEF, eq_concentration->domaine_Cl_dis());
  const DoubleTab& temper = eq_thermique->inconnue().valeurs(), &concen = eq_concentration->inconnue().valeurs();
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base, eq_thermique->get_modele(TURBULENCE).valeur());

  // XXX : Elie Saikali : vaut mieux utiliser diffusivite_turbulente au lie de faire ca ....
  // voila dans Source_Transport_K_Eps_Realisable_aniso_therm_concen_VEF_Face
  // const DoubleTab& alpha_turb = le_modele_scalaire.diffusivite_turbulente().valeurs();
  DoubleTab alpha_turb(le_modele_scalaire.conductivite_turbulente()->valeurs());
  double rhocp = eq_thermique->milieu().capacite_calorifique()->valeurs()(0, 0) * eq_thermique->milieu().masse_volumique()->valeurs()(0, 0);
  alpha_turb /= rhocp;
  const DoubleVect& g = gravite->valeurs(), &volumes_entrelaces = le_dom_VEF->volumes_entrelaces();
  const Champ_Don& ch_beta_temper = beta_t.valeur(), &ch_beta_concen = beta_c.valeur();
  const int nb_face = le_dom_VEF->nb_faces(), nb_consti = eq_concentration->constituant().nb_constituants();
  DoubleTrav G_t(nb_face), G_c(nb_face);

  calculer_terme_destruction_K_gen(le_dom_VEF.valeur(), zcl_VEF_th, G_t, temper, alpha_turb, ch_beta_temper, g, 0);
  calculer_terme_destruction_K_gen(le_dom_VEF.valeur(), zcl_VEF_co, G_c, concen, alpha_turb, ch_beta_concen, g, nb_consti);

  fill_resu_anisotherme_concen(G_t, G_c, volumes_entrelaces,resu); // voir les classes filles
  return resu;
}
