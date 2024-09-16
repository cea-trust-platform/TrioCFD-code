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
// File:        Modele_turbulence_hyd_K_Eps_Bicephale.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <Modifier_pour_fluide_dilatable.h>
#include <Modele_turbulence_scal_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <Champ_Uniforme.h>
#include <communications.h>
#include <stat_counters.h>
#include <Probleme_base.h>
#include <Fluide_base.h>
#include <TRUSTTrav.h>
#include <Param.h>
#include <Debog.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_Bicephale, "Modele_turbulence_hyd_K_Epsilon_Bicephale", Modele_turbulence_hyd_RANS_Bicephale_base);
// XD K_Epsilon_Bicephale mod_turb_hyd_rans_bicephale K_Epsilon_Bicephale -1 Turbulence model (k-eps) en formalisation bicephale.

Sortie& Modele_turbulence_hyd_K_Eps_Bicephale::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Modele_turbulence_hyd_K_Eps_Bicephale::readOn(Entree& s) { return Modele_turbulence_hyd_RANS_Bicephale_base::readOn(s); }

void Modele_turbulence_hyd_K_Eps_Bicephale::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_Bicephale_base::set_param(param);
  param.ajouter_non_std("Transport_K", (this), Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (k) transportation equation.
  param.ajouter_non_std("Transport_Epsilon", (this), Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (eps) transportation equation.
  param.ajouter_non_std("Modele_Fonc_Bas_Reynolds", (this)); // XD_ADD_P Modele_Fonc_Realisable_base This keyword is used to set the model used
  param.ajouter("CMU", &LeCmu_); // XD_ADD_P double Keyword to modify the Cmu constant of k-eps model : Nut=Cmu*k*k/eps Default value is 0.09
}

int Modele_turbulence_hyd_K_Eps_Bicephale::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Modele_Fonc_Bas_Reynolds")
    {
      Cerr << "Lecture du modele bas reynolds associe " << finl;
      Modele_Fonc_Bas_Reynolds_Base::typer_lire_Modele_Fonc_Bas_Reynolds(mon_modele_fonc_, eqn_transp_K(), is);
      mon_modele_fonc_->associer_eqn_2(eqn_transp_Eps());
      Cerr << "mon_modele_fonc.que_suis_je() avant discretisation " << mon_modele_fonc_->que_suis_je() << finl;
      mon_modele_fonc_->discretiser();
      Cerr << "mon_modele_fonc.que_suis_je() " << mon_modele_fonc_->que_suis_je() << finl;
      mon_modele_fonc_->lire_distance_paroi();
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_Bicephale_base::lire_motcle_non_standard(mot, is);
}

/*! @brief Calcule la viscosite turbulente au temps demande.
 *
 * @param (double temps) le temps auquel il faut calculer la viscosite
 * @return (Champ_Fonc&) la viscosite turbulente au temps demande
 * @throws erreur de taille de visco_turb_K_eps
 */
Champ_Fonc& Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK = eqn_transp_K().inconnue().valeur(), &chEps = eqn_transp_Eps().inconnue().valeur();
  const Domaine_dis_base& le_dom_dis = eqn_transp_K().domaine_dis();
  const Domaine_Cl_dis& le_dom_Cl_dis = eqn_transp_K().domaine_Cl_dis();

  const Nom& type = chK.que_suis_je();
  const DoubleTab& tab_K = chK.valeurs(), &tab_Eps = chEps.valeurs();

  DoubleTab& visco_turb = la_viscosite_turbulente_->valeurs();

  DoubleTrav Cmu(tab_K.dimension_tot(0));
  Cmu = 0.;

  int n = tab_K.dimension(0);
  if (n < 0)
    {
      if (sub_type(Champ_Inc_P0_base, chK))
        n = eqn_transp_K().domaine_dis().domaine().nb_elem();
      else
        {
          Cerr << "Unsupported K field in Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente" << finl;
          Process::exit(-1);
        }
      if (sub_type(Champ_Inc_P0_base, chEps))
        n = eqn_transp_Eps().domaine_dis().domaine().nb_elem();
      else
        {
          Cerr << "Unsupported epsilon field in Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente" << finl;
          Process::exit(-1);
        }
    }

  DoubleTrav Fmu, D(tab_K.dimension_tot(0));
  D = 0;
  if (mon_modele_fonc_.non_nul())
    {
      // pour avoir nu en incompressible et mu en QC
      // et non comme on a divise K et eps par rho (si on est en QC) on veut toujours nu
      const Champ_Don ch_visco = ref_cast(Fluide_base,eqn_transp_K().milieu()).viscosite_cinematique();
      const Champ_Don& ch_visco_cin = ref_cast(Fluide_base,eqn_transp_K().milieu()).viscosite_cinematique();

      const DoubleTab& tab_visco = ch_visco_cin->valeurs();

      Fmu.resize(tab_K.dimension_tot(0));

      mon_modele_fonc_->Calcul_Fmu_BiK(Fmu, le_dom_dis, le_dom_Cl_dis, tab_K, tab_Eps, ch_visco);

      int is_Cmu_constant = mon_modele_fonc_->Calcul_is_Cmu_constant();
      if (is_Cmu_constant == 0)
        {
          Cerr << " On utilise un Cmu non constant " << finl;
          const DoubleTab& vitesse = mon_equation_->inconnue()->valeurs();
          mon_modele_fonc_->Calcul_Cmu_BiK(Cmu, le_dom_dis, le_dom_Cl_dis, vitesse, tab_K, tab_Eps, EPS_MIN_);

          /*Paroi*/
          Nom lp = eqn_transp_K().modele_turbulence().loi_paroi()->que_suis_je();
          if (lp != "negligeable_VEF")
            {
              DoubleTab visco_tab(visco_turb.dimension_tot(0));
              assert(sub_type(Champ_Uniforme,ch_visco_cin.valeur()));
              visco_tab = tab_visco(0, 0);
              const int idt = mon_equation_->schema_temps().nb_pas_dt();
              const DoubleTab& tab_paroi = loi_paroi()->Cisaillement_paroi();
              mon_modele_fonc_->Calcul_Cmu_Paroi_BiK(Cmu, le_dom_dis, le_dom_Cl_dis, visco_tab, visco_turb, tab_paroi, idt, vitesse, tab_K, tab_Eps, EPS_MIN_);
            }
        }
      else
        Cerr << " On utilise un Cmu constant " << finl;
    }

  // dans le cas d'un domaine nul on doit effectuer le dimensionnement
  double non_prepare = 1;
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente la_viscosite_turbulente before", la_viscosite_turbulente_->valeurs());
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente tab_K", tab_K);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente tab_Eps", tab_Eps);
  if (visco_turb.size() == n)
    non_prepare = 0.;
  non_prepare = mp_max(non_prepare);

  if (non_prepare == 1)
    {
      Champ_Inc visco_turb_au_format_K_eps;
      visco_turb_au_format_K_eps.typer(type);
      DoubleTab& visco_turb_K_eps = complete_viscosity_field(n, eqn_transp_K().domaine_dis(), visco_turb_au_format_K_eps);

      if (visco_turb_K_eps.size() != n)
        {
          Cerr << "visco_turb_K_eps size is " << visco_turb_K_eps.size() << " instead of " << n << finl;
          exit();
        }
      // A la fin de cette boucle, le tableau visco_turb_K_eps contient les valeurs de la viscosite turbulente
      // au centre des faces du maillage.
      fill_turbulent_viscosity_tab(n, tab_K, tab_Eps, Cmu, Fmu, D, visco_turb_K_eps);

      // On connait donc la viscosite turbulente au centre des faces de chaque element
      // On cherche maintenant a interpoler cette viscosite turbulente au centre des elements.
      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_eps.valeur());
      Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente visco_turb_au_format_K_eps", visco_turb_au_format_K_eps.valeur());
    }
  else
    fill_turbulent_viscosity_tab(n, tab_K, tab_Eps, Cmu, Fmu, D, visco_turb);

  la_viscosite_turbulente_->changer_temps(temps);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente la_viscosite_turbulente after", la_viscosite_turbulente_->valeurs());
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_K_Eps_Bicephale::fill_turbulent_viscosity_tab(const int n, const DoubleTab& tab_K, const DoubleTab& tab_Eps, const DoubleTab& Cmu, const DoubleTab& Fmu, const DoubleTab& D, DoubleTab& turbulent_viscosity)
{
  for (int i = 0; i < n; i++)
    {
      if (tab_Eps(i) <= EPS_MIN_)
        {
          turbulent_viscosity[i] = 0.;
        }
      else
        {
          if (mon_modele_fonc_.non_nul())
            {
              int is_Cmu_constant = mon_modele_fonc_->Calcul_is_Cmu_constant();
              if (is_Cmu_constant)
                turbulent_viscosity[i] = Fmu(i) * LeCmu_ * tab_K(i) * tab_K(i) / (tab_Eps(i) + D(i));
              else
                turbulent_viscosity[i] = Fmu(i) * Cmu(i) * tab_K(i) * tab_K(i) / (tab_Eps(i) + D(i));
            }
          else
            turbulent_viscosity[i] = LeCmu_ * tab_K(i) * tab_K(i) / (tab_Eps(i) + D(i));
        }
    }
}

int Modele_turbulence_hyd_K_Eps_Bicephale::preparer_calcul()
{
  eqn_transp_K().preparer_calcul();
  eqn_transp_Eps().preparer_calcul();
  Modele_turbulence_hyd_base::preparer_calcul();
  // GF pour initialiser la loi de paroi thermique en TBLE
  if (equation().probleme().nombre_d_equations() > 1)
    {
      const RefObjU& modele_turbulence = equation().probleme().equation(1).get_modele(TURBULENCE);
      if (sub_type(Modele_turbulence_scal_base, modele_turbulence.valeur()))
        {
          Turbulence_paroi_scal& loi_paroi_T = ref_cast_non_const(Modele_turbulence_scal_base,modele_turbulence.valeur()).loi_paroi();
          loi_paroi_T->init_lois_paroi();
        }
    }

  calculate_limit_viscosity<MODELE_TYPE::K_EPS_BICEPHALE>(K(), Eps(), LeCmu_);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::preparer_calcul la_viscosite_turbulente", la_viscosite_turbulente_->valeurs());
  return 1;
}

bool Modele_turbulence_hyd_K_Eps_Bicephale::initTimeStep(double dt)
{
  return (eqn_transport_K_.initTimeStep(dt) and eqn_transport_Eps_.initTimeStep(dt));
}

/*! @brief Effectue une mise a jour en temps du modele de turbulence.
 *
 * Met a jour les equations de transport de K et epsilon,
 *     calcule la viscosite turbulente
 *     au nouveau temps.
 *
 * @param (double temps) le temps de mise a jour
 */
void Modele_turbulence_hyd_K_Eps_Bicephale::mettre_a_jour(double temps)
{
  Schema_Temps_base& sch = eqn_transp_K().schema_temps(), &sch2 = eqn_transp_Eps().schema_temps();
  eqn_transp_K().domaine_Cl_dis()->mettre_a_jour(temps);
  eqn_transp_Eps().domaine_Cl_dis()->mettre_a_jour(temps);

  if (!eqn_transp_K().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K());
  eqn_transp_K().mettre_a_jour(temps);

  if (!eqn_transp_Eps().equation_non_resolue())
    sch2.faire_un_pas_de_temps_eqn_base(eqn_transp_Eps());
  eqn_transp_Eps().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::mettre_a_jour la_viscosite_turbulente before", la_viscosite_turbulente_->valeurs());
  calculate_limit_viscosity<MODELE_TYPE::K_EPS_BICEPHALE>(K(), Eps(), LeCmu_);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::mettre_a_jour la_viscosite_turbulente after", la_viscosite_turbulente_->valeurs());
  statistiques().end_count(nut_counter_);
}

const Equation_base& Modele_turbulence_hyd_K_Eps_Bicephale::equation_k_eps(int i) const
{
  assert((i == 0) || (i == 1));
  if (i == 0)
    {
      return eqn_transport_K_;
    }
  else
    {
      return eqn_transport_Eps_;
    }
}
const Champ_base& Modele_turbulence_hyd_K_Eps_Bicephale::get_champ(const Motcle& nom) const
{

  try
    {
      return Modele_turbulence_hyd_RANS_Bicephale_base::get_champ(nom);
    }
  catch (Champs_compris_erreur&)
    {
    }
  if (mon_modele_fonc_.non_nul())
    {
      try
        {
          return mon_modele_fonc_->get_champ(nom);
        }
      catch (Champs_compris_erreur&)
        {
        }
    }
  throw Champs_compris_erreur();

}
void Modele_turbulence_hyd_K_Eps_Bicephale::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_RANS_Bicephale_base::get_noms_champs_postraitables(nom, opt);
  if (mon_modele_fonc_.non_nul())
    mon_modele_fonc_->get_noms_champs_postraitables(nom, opt);

}
void Modele_turbulence_hyd_K_Eps_Bicephale::verifie_loi_paroi()
{
  Nom lp = loipar_->que_suis_je();
  if (!(lp == "negligeable_VEF" || lp == "negligeable_VDF"))
    {
      Cerr << "The turbulence model of type " << que_suis_je() << finl;
      Cerr << "must be used with a wall law of type negligeable for now." << finl;
    }
  if (!associe_modele_fonction().non_nul())
    {
      Cerr << "The turbulence model of type " << que_suis_je() << finl;
      Cerr << "must be used with a modele fonction for now." << finl;
    }
}
