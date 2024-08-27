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
// File:        Modele_turbulence_hyd_K_Eps_Bas_Reynolds.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <Schema_Temps_base.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <Debog.h>
#include <stat_counters.h>
#include <Param.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_Bas_Reynolds, "Modele_turbulence_hyd_K_Epsilon_Bas_Reynolds", Modele_turbulence_hyd_RANS_K_Eps_base);

Sortie& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::readOn(Entree& is)
{
  return Modele_turbulence_hyd_RANS_K_Eps_base::readOn(is);
}

void Modele_turbulence_hyd_K_Eps_Bas_Reynolds::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_K_Eps_base::set_param(param);
  param.ajouter_non_std("Transport_K_Epsilon_Bas_Reynolds", (this), Param::REQUIRED);
  param.ajouter_non_std("Modele_Fonc_Bas_Reynolds", (this), Param::REQUIRED);
}

int Modele_turbulence_hyd_K_Eps_Bas_Reynolds::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Transport_K_Epsilon_Bas_Reynolds")
    {
      eqn_transp_K_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_K_Eps();
      Cerr << "K_Epsilon equation type " << eqn_transp_K_Eps().que_suis_je() << finl;
      return 1;
    }
  else if (mot == "Modele_Fonc_Bas_Reynolds")
    {
      mon_modele_fonc_.associer_eqn(eqn_transp_K_Eps());
      is >> mon_modele_fonc_;
      mon_modele_fonc_->discretiser();
      Cerr << "Low Reynolds number model type " << mon_modele_fonc_->que_suis_je() << finl;
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_K_Eps_base::lire_motcle_non_standard(mot, is);
}

Champ_Fonc& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::calculer_viscosite_turbulente(double temps)
{

  const Champ_base& chK_Eps = eqn_transp_K_Eps().inconnue().valeur();
  const Domaine_dis& le_dom_dis = eqn_transp_K_Eps().domaine_dis();
  const Domaine_Cl_dis& le_dom_Cl_dis = eqn_transp_K_Eps().domaine_Cl_dis();
  Nom type = chK_Eps.que_suis_je();
  const DoubleTab& tab_K_Eps = chK_Eps.valeurs();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bas_Reynolds::calculer_viscosite_turbulente K_Eps", tab_K_Eps);
  DoubleTab& visco_turb = la_viscosite_turbulente_->valeurs();
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_transp_K_Eps().milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();
  int n = tab_K_Eps.dimension(0);
  DoubleTab Fmu(n);

  mon_modele_fonc_.Calcul_Fmu(Fmu, le_dom_dis, le_dom_Cl_dis, tab_K_Eps, ch_visco_cin);

  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bas_Reynolds::calculer_viscosite_turbulente Fmu", Fmu);

  //
  //  limiteur Durbin
  //
  //        double T_durbin, T_kolmo, T_ke;
  // dans le cas d'un domaine nul on doit effectuer le dimensionnement
  double non_prepare = 1;
  if (visco_turb.size() == n)
    non_prepare = 0.;
  non_prepare = mp_max(non_prepare);

  if (non_prepare == 1)
    {
      Champ_Inc visco_turb_au_format_K_eps_Bas_Re;
      visco_turb_au_format_K_eps_Bas_Re.typer(type);
      DoubleTab& visco_turb_K_eps_Bas_Re = complete_viscosity_field(n, eqn_transp_K_Eps().domaine_dis().valeur(), visco_turb_au_format_K_eps_Bas_Re);

      if (visco_turb_K_eps_Bas_Re.size() != n)
        {
          Cerr << "visco_turb_K_eps_Bas_Re size is " << visco_turb_K_eps_Bas_Re.size() << " instead of " << n << finl;
          exit();
        }

      fill_turbulent_viscosity_tab(n, tab_K_Eps, Fmu, visco_turb_K_eps_Bas_Re);

      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_eps_Bas_Re.valeur());
    }
  else
    fill_turbulent_viscosity_tab(n, tab_K_Eps, Fmu, visco_turb);

  la_viscosite_turbulente_->changer_temps(temps);
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_K_Eps_Bas_Reynolds::fill_turbulent_viscosity_tab(const int n, const DoubleTab& tab_K_Eps, const DoubleTab& Fmu, DoubleTab& turbulent_viscosity)
{
  for (int i = 0; i < n; i++)
    {
      if (tab_K_Eps(i, 1) <= DMINFLOAT)
        turbulent_viscosity[i] = 0.;
      else
        turbulent_viscosity[i] = LeCmu_ * Fmu(i) * tab_K_Eps(i, 0) * tab_K_Eps(i, 0) / tab_K_Eps(i, 1);
    }
}

int Modele_turbulence_hyd_K_Eps_Bas_Reynolds::preparer_calcul()
{
  return eqn_transp_K_Eps().preparer_calcul();
}

void Modele_turbulence_hyd_K_Eps_Bas_Reynolds::mettre_a_jour(double temps)
{
  Schema_Temps_base& sch = eqn_transp_K_Eps().schema_temps();
  eqn_transp_K_Eps().domaine_Cl_dis()->mettre_a_jour(temps);
  sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K_Eps());
  eqn_transp_K_Eps().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  calculate_limit_viscosity<MODELE_TYPE::K_EPS_BAS_REYNOLDS>(K_Eps(), -123. /* unused */);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bas_Reynolds::mettre_a_jour apres calculer_viscosite_turbulente la_viscosite_turbulente", la_viscosite_turbulente_->valeurs());
  statistiques().end_count(nut_counter_);
}

bool Modele_turbulence_hyd_K_Eps_Bas_Reynolds::initTimeStep(double dt)
{
  return eqn_transport_K_Eps_Bas_Re_.initTimeStep(dt);
}

const Equation_base& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::equation_k_eps(int i) const
{
  assert((i == 0));
  return eqn_transport_K_Eps_Bas_Re_;
}

const Champ_base& Modele_turbulence_hyd_K_Eps_Bas_Reynolds::get_champ(const Motcle& nom) const
{
  try
    {
      return Modele_turbulence_hyd_RANS_K_Eps_base::get_champ(nom);
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

void Modele_turbulence_hyd_K_Eps_Bas_Reynolds::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_RANS_K_Eps_base::get_noms_champs_postraitables(nom, opt);

  if (mon_modele_fonc_.non_nul())
    mon_modele_fonc_->get_noms_champs_postraitables(nom, opt);

}

void Modele_turbulence_hyd_K_Eps_Bas_Reynolds::completer()
{
  eqn_transp_K_Eps().completer();
}
