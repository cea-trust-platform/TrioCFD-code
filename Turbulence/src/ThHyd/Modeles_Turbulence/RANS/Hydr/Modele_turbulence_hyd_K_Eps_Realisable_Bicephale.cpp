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
// File:        Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Fluide_Incompressible.h>
#include <Champ_Uniforme.h>
#include <Schema_Temps.h>
#include <Debog.h>
#include <stat_counters.h>
#include <Param.h>
#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Modifier_nut_pour_fluide_dilatable.h>
#include <Modele_turbulence_scal_base.h>
#include <TRUSTTrav.h>
#include <communications.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale,"Modele_turbulence_hyd_K_Epsilon_Realisable_Bicephale",Mod_turb_hyd_RANS_Bicephale);

// XD K_Epsilon_Realisable_Bicephale mod_turb_hyd_rans K_Epsilon_Realisable_Bicephale -1 Realizable Two-headed K-Epsilon Turbulence Model

//
// printOn et readOn

Sortie& Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

/*! @brief Simple appel a Mod_turb_hyd_RANS_Bicephale::readOn(Entree&)
 *
 */
Entree& Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::readOn(Entree& is)
{
  return Mod_turb_hyd_RANS_Bicephale::readOn(is);
}

void Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::set_param(Param& param)
{
  Mod_turb_hyd_RANS_Bicephale::set_param(param);
  param.ajouter_non_std("Transport_K",(this),Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (k) transportation equation.
  param.ajouter_non_std("Transport_Epsilon",(this),Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (eps) transportation equation.
  param.ajouter_non_std("Modele_Fonc_Realisable",(this),Param::REQUIRED); // XD_ADD_P Modele_Fonc_Realisable_base This keyword is used to set the model used
  param.ajouter("PRANDTL_K",&Prandtl_K,Param::REQUIRED); // XD_ADD_P double Keyword to change the Prk value (default 1.0).
  param.ajouter("PRANDTL_EPS",&Prandtl_Eps,Param::REQUIRED); // XD_ADD_P double Keyword to change the Pre value (default 1.3)
}

int Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="Transport_K")
    {
      eqn_transp_K().transporte_K( );
      eqn_transp_K().associer_modele_turbulence(*this);
      is >> eqn_transp_K();
      return 1;
    }
  else if (mot=="Transport_Epsilon")
    {
      eqn_transp_Eps().transporte_Eps( );
      eqn_transp_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_Eps();
      return 1;
    }
  else if (mot=="Modele_Fonc_Realisable")
    {
      get_modele_fonction().associer_eqn(eqn_transp_K());
      is >> mon_modele_fonc;
      mon_modele_fonc.associer_eqn_2(eqn_transp_Eps());
      get_modele_fonction().discretiser();
      Cerr << "Realizable Two-headed K_Epsilon model type " << get_modele_fonction().que_suis_je() << finl;
      return 1;
    }
  else
    return Mod_turb_hyd_RANS_Bicephale::lire_motcle_non_standard(mot,is);
}

Champ_Fonc& Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK=eqn_transp_K().inconnue().valeur();
  Nom type=chK.que_suis_je();
  const DoubleTab& tab_K = chK.valeurs();
  const Champ_base& chEps=eqn_transp_Eps().inconnue().valeur();
  const DoubleTab& tab_Eps = chEps.valeurs();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::calculer_viscosite_turbulente K",tab_K);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::calculer_viscosite_turbulente Epsilon",tab_Eps);
  DoubleTab& visco_turb =  la_viscosite_turbulente.valeurs();

  int n = tab_K.dimension(0);

  const DoubleTab& Cmu = get_modele_fonction().get_Cmu(); // attention : il faut qu'il soit deja calcule!

  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::calculer_viscosite_turbulente Cmu",Cmu);

  // dans le cas d'une domaine nulle on doit effectuer le dimensionnement
  double non_prepare=1;
  if (visco_turb.size() == n)
    non_prepare=0.;
  non_prepare=mp_max(non_prepare);

  if (non_prepare==1)
    {
      Champ_Inc visco_turb_au_format_K_eps_Rea;
      visco_turb_au_format_K_eps_Rea.typer(type);
      Champ_Inc_base& ch_visco_turb_K_eps_Rea=visco_turb_au_format_K_eps_Rea.valeur();
      ch_visco_turb_K_eps_Rea.associer_domaine_dis_base(eqn_transp_K().domaine_dis().valeur());
      ch_visco_turb_K_eps_Rea.nommer("diffusivite_turbulente");
      ch_visco_turb_K_eps_Rea.fixer_nb_comp(1);
      ch_visco_turb_K_eps_Rea.fixer_nb_valeurs_nodales(n);
      ch_visco_turb_K_eps_Rea.fixer_unite("inconnue");
      ch_visco_turb_K_eps_Rea.changer_temps(0.);

      DoubleTab& visco_turb_K_eps_Rea =  ch_visco_turb_K_eps_Rea.valeurs();
      if(visco_turb_K_eps_Rea.size() != n)
        {
          Cerr << "visco_turb_K_eps_Rea size is " << visco_turb_K_eps_Rea.size()
               << " instead of " << n << finl;
          exit();
        }

      for (int i=0; i<n; i++)
        {
          if (tab_Eps(i,1) <= DMINFLOAT)
            visco_turb_K_eps_Rea[i] = 0;
          else
            visco_turb_K_eps_Rea[i] = Cmu(i)*tab_K(i)*tab_K(i)/tab_Eps(i);

        }

      la_viscosite_turbulente->affecter(visco_turb_au_format_K_eps_Rea.valeur());

    }
  else
    {
      for (int i=0; i<n; i++)
        {
          if (tab_Eps(i,1) <= DMINFLOAT)
            visco_turb[i] = 0;
          else
            visco_turb[i] = Cmu(i)*tab_K(i)*tab_K(i)/tab_Eps(i);
        }
    }
  la_viscosite_turbulente.changer_temps(temps);
  return la_viscosite_turbulente;
}

int Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::preparer_calcul()
{
  eqn_transp_K().preparer_calcul();
  eqn_transp_Eps().preparer_calcul();
  Mod_turb_hyd_base::preparer_calcul();
  // GF pour initialiser la loi de paroi thermique en TBLE
//   if (equation().probleme().nombre_d_equations()>1)
//     {
//       const RefObjU& modele_turbulence = equation().probleme().equation(1).get_modele(TURBULENCE);
//       if (sub_type(Modele_turbulence_scal_base,modele_turbulence.valeur()))
//         {
//           Turbulence_paroi_scal& loi_paroi_T = ref_cast_non_const(Modele_turbulence_scal_base,modele_turbulence.valeur()).loi_paroi();
//           loi_paroi_T.init_lois_paroi();
//         }
//     }
  // GF quand on demarre un calcul il est bon d'utliser la ldp
  // encore plus quand on fait une reprise !!!!!!!!
  Champ_Inc& ch_K   = K();
  Champ_Inc& ch_Eps = Eps();

  const Milieu_base& mil=equation().probleme().milieu();
  if (equation().probleme().is_dilatable()) diviser_par_rho_si_dilatable(ch_K.valeurs(),mil);
  diviser_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
  loipar.calculer_hyd_BiK(ch_K,ch_Eps);
  eqn_transp_K().controler_variable();
  eqn_transp_Eps().controler_variable();
  calculer_viscosite_turbulente(K().temps());
  limiter_viscosite_turbulente();
  // on remultiplie K et Eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K.valeurs(),mil);
      multiplier_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
  return 1;

}

void Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::mettre_a_jour(double temps)
{
  Champ_Inc& ch_K   = K();
  Champ_Inc& ch_Eps = Eps();
  Schema_Temps_base& sch = eqn_transp_K().schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base
  eqn_transp_K().domaine_Cl_dis().mettre_a_jour(temps);
  if (!eqn_transp_K().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K());
  eqn_transp_K().mettre_a_jour(temps);
  if (!eqn_transp_Eps().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_Eps());
  eqn_transp_Eps().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  const Milieu_base& mil=equation().probleme().milieu();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::mettre_a_jour la_viscosite_turbulente before",la_viscosite_turbulente.valeurs());
  // on divise K_eps par rho en QC pour revenir a K et Eps
  if (equation().probleme().is_dilatable()) diviser_par_rho_si_dilatable(ch_K.valeurs(),mil);
  diviser_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
  loipar.calculer_hyd_BiK(ch_K,ch_Eps);
  eqn_transp_K().controler_variable();
  eqn_transp_Eps().controler_variable();
  calculer_viscosite_turbulente(ch_K.temps());
  limiter_viscosite_turbulente();
  // on remultiplie K_eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K.valeurs(),mil);
      multiplier_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::mettre_a_jour apres calculer_viscosite_turbulente la_viscosite_turbulente",la_viscosite_turbulente.valeurs());
  statistiques().end_count(nut_counter_);
}

bool Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::initTimeStep(double dt)
{
  return ( eqn_transport_K_Rea.initTimeStep(dt) and eqn_transport_Eps_Rea.initTimeStep(dt) );
}

const Equation_base& Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::equation_k_eps(int i) const
{
  assert ((i==0)||(i==1));
  if (i==0)
    {
      return eqn_transport_K_Rea;
    }
  else
    {
      return eqn_transport_Eps_Rea;
    }
}

const Champ_base& Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::get_champ(const Motcle& nom) const
{
  REF(Champ_base) ref_champ;
  try
    {
      return Mod_turb_hyd_RANS_Bicephale::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }

  if (mon_modele_fonc.non_nul())
    {
      try
        {
          return mon_modele_fonc->get_champ(nom);
        }
      catch (Champs_compris_erreur)
        {
        }
    }

  throw Champs_compris_erreur();
}

void Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Mod_turb_hyd_RANS_Bicephale::get_noms_champs_postraitables(nom,opt);
  if (mon_modele_fonc.non_nul())
    mon_modele_fonc->get_noms_champs_postraitables(nom,opt);
}

void Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::verifie_loi_paroi()
{
  Nom lp=loipar.valeur().que_suis_je();
  if (lp=="negligeable_VEF" || lp=="negligeable_VDF")
    if (!associe_modele_fonction().non_nul())
      {
        Cerr<<"The turbulence model of type "<<que_suis_je()<<finl;
        Cerr<<"must not be used with a wall law of type negligeable or with a modele_function."<<finl;
        Cerr<<"Another wall law must be selected with this kind of turbulence model."<<finl;
      }
}

void Modele_turbulence_hyd_K_Eps_Realisable_Bicephale::completer()
{
  eqn_transp_K().completer();
  eqn_transp_Eps().completer();
  verifie_loi_paroi();
}
