/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Modele_turbulence_hyd_K_Eps_2_Couches.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_2_Couches.h>
#include <Schema_Temps_base.h>
#include <Modifier_pour_fluide_dilatable.h>
#include <Probleme_base.h>
#include <stat_counters.h>
#include <Param.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_2_Couches, "Modele_turbulence_hyd_K_Epsilon_2_Couches", Modele_turbulence_hyd_RANS_K_Eps_base);

Sortie& Modele_turbulence_hyd_K_Eps_2_Couches::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_K_Eps_2_Couches::readOn(Entree& s)
{
  return Modele_turbulence_hyd_RANS_K_Eps_base::readOn(s);
}

void Modele_turbulence_hyd_K_Eps_2_Couches::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_K_Eps_base::set_param(param);
  param.ajouter_non_std("Transport_K_KEpsilon", (this), Param::REQUIRED);
}

int Modele_turbulence_hyd_K_Eps_2_Couches::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Transport_K_KEpsilon")
    {
      eqn_transp_K_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_K_Eps();
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_K_Eps_base::lire_motcle_non_standard(mot, is);
}

/*! @brief Calcule la viscosite turbulente au temps demande.
 *
 * @param (double temps) le temps auquel il faut calculer la viscosite
 * @return (Champ_Fonc_base&) la viscosite turbulente au temps demande
 * @throws erreur de taille de visco_turb_K_eps
 */
Champ_Fonc_base& Modele_turbulence_hyd_K_Eps_2_Couches::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK_Eps = eqn_transport_K_Eps_.inconnue();
  Nom type = chK_Eps.que_suis_je();
  const DoubleTab& tab_K_Eps = chK_Eps.valeurs();
  DoubleTab& visco_turb = la_viscosite_turbulente_->valeurs();

  // K_Eps(i,0) = K au noeud i
  // K_Eps(i,1) = Epsilon au noeud i

  int n = tab_K_Eps.dimension(0);
  if (visco_turb.size() != n)
    {
      OWN_PTR(Champ_Inc_base) visco_turb_au_format_K_eps;
      visco_turb_au_format_K_eps.typer(type);
      DoubleTab& visco_turb_K_eps = complete_viscosity_field(n, eqn_transport_K_Eps_.domaine_dis(), visco_turb_au_format_K_eps);

      if (visco_turb_K_eps.size() != n)
        {
          Cerr << "visco_turb_K_eps size is " << visco_turb_K_eps.size() << " instead of " << n << finl;
          exit();
        }

      fill_turbulent_viscosity_tab(n, tab_K_Eps, visco_turb_K_eps);
      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_eps.valeur());
    }
  else
    fill_turbulent_viscosity_tab(n, tab_K_Eps, visco_turb);

  la_viscosite_turbulente_->changer_temps(temps);
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_K_Eps_2_Couches::fill_turbulent_viscosity_tab(const int n, const DoubleTab& tab_K_Eps,  DoubleTab& turbulent_viscosity)
{
  for (int i = 0; i < n; i++)
    {
      if (tab_K_Eps(i, 1) <= EPS_MIN_)
        turbulent_viscosity[i] = 0.;
      else
        turbulent_viscosity[i] = CMU * tab_K_Eps(i, 0) * tab_K_Eps(i, 0) / tab_K_Eps(i, 1);
    }
}

int Modele_turbulence_hyd_K_Eps_2_Couches::preparer_calcul()
{
  eqn_transp_K_Eps().preparer_calcul();
  calculer_viscosite_turbulente(equation().schema_temps().temps_courant());
  Modele_turbulence_hyd_base::preparer_calcul();
  calculer_viscosite_turbulente(K_Eps().temps());
  la_viscosite_turbulente_->valeurs().echange_espace_virtuel();
  return 1;
}

bool Modele_turbulence_hyd_K_Eps_2_Couches::initTimeStep(double dt)
{
  return eqn_transport_K_Eps_.initTimeStep(dt);
}

/*! @brief Effectue une mise a jour en temps du modele de turbulence.
 *
 * Met a jour l'equation de transport K-epsilon,
 *     calcule la loi de paroi et la viscosite turbulente
 *     au nouveau temps.
 *
 * @param (double temps) le temps de mise a jour
 */
void Modele_turbulence_hyd_K_Eps_2_Couches::mettre_a_jour(double temps)
{
  Schema_Temps_base& sch = eqn_transport_K_Eps_.schema_temps();
  eqn_transport_K_Eps_.domaine_Cl_dis().mettre_a_jour(temps);
  sch.faire_un_pas_de_temps_eqn_base(eqn_transport_K_Eps_);
  eqn_transport_K_Eps_.mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  calculate_limit_viscosity<MODELE_TYPE::K_EPS_2_COUCHES>(K_Eps(), LeCmu_);
  statistiques().end_count(nut_counter_);
}

/*! @brief Simple appel a Transport_K_Eps::completer()
 *
 */
void Modele_turbulence_hyd_K_Eps_2_Couches::completer()
{
  eqn_transp_K_Eps().completer();
  verifie_loi_paroi();
}
