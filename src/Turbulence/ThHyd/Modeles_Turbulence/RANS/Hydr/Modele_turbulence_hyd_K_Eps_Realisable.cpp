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
// File:        Modele_turbulence_hyd_K_Eps_Realisable.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Modele_turbulence_scal_base.h>
#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Fluide_Incompressible.h>
#include <Champ_Inc_P0_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <communications.h>
#include <Probleme_base.h>
#include <stat_counters.h>
#include <TRUSTTrav.h>
#include <Param.h>
#include <Debog.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_Realisable, "Modele_turbulence_hyd_K_Epsilon_Realisable", Modele_turbulence_hyd_RANS_K_Eps_base);

// XD K_Eps_Realisable mod_turb_hyd_rans K_Epsilon_Realisable -1 Realizable K-Epsilon Turbulence Model.

Sortie& Modele_turbulence_hyd_K_Eps_Realisable::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_K_Eps_Realisable::readOn(Entree& is)
{
  return Modele_turbulence_hyd_RANS_K_Eps_base::readOn(is);
}

void Modele_turbulence_hyd_K_Eps_Realisable::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_K_Eps_base::set_param(param);
  param.ajouter_non_std("Transport_K_Epsilon_Realisable", (this), Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (k-eps) transportation equation.
  param.ajouter_non_std("Modele_Fonc_Realisable", (this), Param::REQUIRED); // XD_ADD_P Modele_Fonc_Realisable_base This keyword is used to set the model used
  param.ajouter("PRANDTL_K", &Prandtl_K_, Param::REQUIRED); // XD_ADD_P double Keyword to change the Prk value (default 1.0).
  param.ajouter("PRANDTL_EPS", &Prandtl_Eps_, Param::REQUIRED); // XD_ADD_P double Keyword to change the Pre value (default 1.3)
}

int Modele_turbulence_hyd_K_Eps_Realisable::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Transport_K_Epsilon_Realisable")
    {
      eqn_transp_K_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_K_Eps();
      Cerr << "Realizable K_Epsilon equation type " << eqn_transp_K_Eps().que_suis_je() << finl;
      return 1;
    }
  else if (mot == "Modele_Fonc_Realisable")
    {
      Modele_Fonc_Realisable_base::typer_lire_Modele_Fonc_Realisable(mon_modele_fonc_, eqn_transp_K_Eps(), is);
      get_modele_fonction()->discretiser();
      Cerr << "Realizable K_Epsilon model type " << get_modele_fonction().que_suis_je() << finl;
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_K_Eps_base::lire_motcle_non_standard(mot, is);
}

Champ_Fonc& Modele_turbulence_hyd_K_Eps_Realisable::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK_Eps = eqn_transp_K_Eps().inconnue();
  const Nom& type = chK_Eps.que_suis_je();
  const DoubleTab& tab_K_Eps = chK_Eps.valeurs();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::calculer_viscosite_turbulente K_Eps", tab_K_Eps);
  DoubleTab& visco_turb = la_viscosite_turbulente_->valeurs();

  int n = tab_K_Eps.dimension(0);

  // dans le cas d'un domaine nul on doit effectuer le dimensionnement
  double non_prepare = 1;
  if (visco_turb.size() == n)
    non_prepare = 0.;
  non_prepare = mp_max(non_prepare);

  if (non_prepare == 1)
    {
      OWN_PTR(Champ_Inc_base) visco_turb_au_format_K_eps_Rea;
      visco_turb_au_format_K_eps_Rea.typer(type);
      DoubleTab& visco_turb_K_eps_Rea = complete_viscosity_field(n, eqn_transp_K_Eps().domaine_dis(), visco_turb_au_format_K_eps_Rea);

      if (visco_turb_K_eps_Rea.size() != n)
        {
          Cerr << "visco_turb_K_eps_Rea size is " << visco_turb_K_eps_Rea.size() << " instead of " << n << finl;
          exit();
        }

      fill_turbulent_viscosity_tab(n, tab_K_Eps, visco_turb_K_eps_Rea);
      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_eps_Rea.valeur());
    }
  else
    fill_turbulent_viscosity_tab(n, tab_K_Eps, visco_turb);

  la_viscosite_turbulente_->changer_temps(temps);
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_K_Eps_Realisable::fill_turbulent_viscosity_tab(const int n, const DoubleTab& tab_K_Eps, DoubleTab& turbulent_viscosity)
{
  const DoubleTab& Cmu = get_modele_fonction()->get_Cmu(); // attention : il faut qu'il soit deja calcule!
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::calculer_viscosite_turbulente Cmu", Cmu);

  for (int i = 0; i < n; i++)
    {
      if (tab_K_Eps(i, 1) <= EPS_MIN_)
        turbulent_viscosity[i] = 0.;
      else
        turbulent_viscosity[i] = Cmu(i) * tab_K_Eps(i, 0) * tab_K_Eps(i, 0) / tab_K_Eps(i, 1);
    }
}

int Modele_turbulence_hyd_K_Eps_Realisable::preparer_calcul()
{
  eqn_transp_K_Eps().preparer_calcul();
  Modele_turbulence_hyd_base::preparer_calcul();
  calculate_limit_viscosity<MODELE_TYPE::K_EPS_REALISABLE>(K_Eps(), LeCmu_);
  return 1;
}

void Modele_turbulence_hyd_K_Eps_Realisable::mettre_a_jour(double temps)
{
  Schema_Temps_base& sch = eqn_transp_K_Eps().schema_temps();
  eqn_transp_K_Eps().domaine_Cl_dis().mettre_a_jour(temps);
  if (!eqn_transp_K_Eps().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K_Eps());
  eqn_transp_K_Eps().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::mettre_a_jour la_viscosite_turbulente before", la_viscosite_turbulente_->valeurs());
  calculate_limit_viscosity<MODELE_TYPE::K_EPS_REALISABLE>(K_Eps(), LeCmu_);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::mettre_a_jour apres calculer_viscosite_turbulente la_viscosite_turbulente", la_viscosite_turbulente_->valeurs());
  statistiques().end_count(nut_counter_);
}

bool Modele_turbulence_hyd_K_Eps_Realisable::initTimeStep(double dt)
{
  return eqn_transport_K_Eps_Rea_.initTimeStep(dt);
}

const Equation_base& Modele_turbulence_hyd_K_Eps_Realisable::equation_k_eps(int i) const
{
  assert((i == 0));
  return eqn_transport_K_Eps_Rea_;
}

const Champ_base& Modele_turbulence_hyd_K_Eps_Realisable::get_champ(const Motcle& nom) const
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

void Modele_turbulence_hyd_K_Eps_Realisable::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_RANS_K_Eps_base::get_noms_champs_postraitables(nom, opt);
  if (mon_modele_fonc_.non_nul())
    mon_modele_fonc_->get_noms_champs_postraitables(nom, opt);
}

void Modele_turbulence_hyd_K_Eps_Realisable::verifie_loi_paroi()
{
  Nom lp = loipar_->que_suis_je();
  if (lp == "negligeable_VEF" || lp == "negligeable_VDF")
    if (!associe_modele_fonction().non_nul())
      {
        Cerr << "The turbulence model of type " << que_suis_je() << finl;
        Cerr << "must not be used with a wall law of type negligeable or with a modele_function." << finl;
        Cerr << "Another wall law must be selected with this kind of turbulence model." << finl;
      }
}

void Modele_turbulence_hyd_K_Eps_Realisable::completer()
{
  eqn_transp_K_Eps().completer();
  verifie_loi_paroi();
}
