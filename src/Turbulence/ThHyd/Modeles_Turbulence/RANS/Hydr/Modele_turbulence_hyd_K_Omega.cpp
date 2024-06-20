/****************************************************************************
* Copyright (c) 2023, CEA
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
// File:        Modele_turbulence_hyd_K_Omega.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modifier_pour_fluide_dilatable.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Modele_turbulence_scal_base.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <communications.h>
#include <Champ_Uniforme.h>
#include <TRUSTTab_parts.h>
#include <Probleme_base.h>
#include <stat_counters.h>
#include <Schema_Temps.h>
#include <Fluide_base.h>
#include <TRUSTTrav.h>
#include <Param.h>
#include <Debog.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Omega, "Modele_turbulence_hyd_K_Omega", Modele_turbulence_hyd_RANS_K_Omega_base);
// XD k_omega mod_turb_hyd_rans k_omega -1 Turbulence model (k-omega).

Sortie& Modele_turbulence_hyd_K_Omega::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_K_Omega::readOn(Entree& s)
{
  Modele_turbulence_hyd_RANS_K_Omega_base::readOn(s);

  if (model_variant_ == "SST")
    {
      Cerr << "SST model: initialize les distances paroi" << "\n";
      // equation().probleme().equation(0).creer_champ("distance_paroi_globale");
      // equation().creer_champ("distance_paroi_globale");

      Navier_Stokes_std& moneq = ref_cast(Navier_Stokes_std, equation());
      moneq.creer_champ("distance_paroi_globale");
    }

  return s;
}

void Modele_turbulence_hyd_K_Omega::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_K_Omega_base::set_param(param);
  param.ajouter_non_std("Transport_K_Omega", (this), Param::REQUIRED); // XD_ADD_P transport_k_omega Keyword to define the (k-omega) transportation equation.
  param.ajouter("PRANDTL_K", &Prandtl_K_);
  param.ajouter("PRANDTL_Omega", &Prandtl_Omega_);
  param.ajouter("model_variant", &model_variant_, Param::OPTIONAL); // XD_ADD_P chaine Model variant for k-omega (default value STD)
}

int Modele_turbulence_hyd_K_Omega::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Transport_K_Omega")
    {
      eqn_transp_K_Omega().associer_modele_turbulence(*this);
      is >> eqn_transp_K_Omega();
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_K_Omega_base::lire_motcle_non_standard(mot, is);
}

/*! @brief Calcule la viscosite turbulente au temps demande.
 *
 * @param (double temps) le temps auquel il faut calculer la viscosite
 * @return (Champ_Fonc&) la viscosite turbulente au temps demande
 * @throws erreur de taille de visco_turb_K_omega
 */
Champ_Fonc& Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK_Omega = eqn_transp_K_Omega().inconnue().valeur();
  Nom type = chK_Omega.que_suis_je();
  // const Domaine_Cl_dis& la_domaine_Cl_dis = eqn_transp_K_Omega().domaine_Cl_dis();
  const DoubleTab& tab_K_Omega = chK_Omega.valeurs();
  DoubleTab& visco_turb = la_viscosite_turbulente_.valeurs();

  // K_Omega(i, 0) = K au noeud i
  // K_Omega(i, 1) = Omega au noeud i

  int n = tab_K_Omega.dimension(0);
  if (n < 0)
    {
      if (!sub_type(Champ_Inc_P0_base, chK_Omega))
        {
          Cerr << "Unsupported K_Omega field in Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente" << finl;
          Process::exit(-1);
        }
      n = eqn_transp_K_Omega().domaine_dis().domaine().nb_elem();
    }

  // cAlan, le 20/01/2023 : sortir cette partie et en faire une fonction à part ?
  // dans le cas d'une domaine nulle on doit effectuer le dimensionnement
  double non_prepare = 1;
  Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente la_viscosite_turbulente before", la_viscosite_turbulente_.valeurs());
  Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente tab_K_Omega", tab_K_Omega);
  if (visco_turb.size() == n)
    non_prepare = 0.;
  non_prepare = mp_max(non_prepare);

  if (non_prepare == 1)
    {
      Champ_Inc visco_turb_au_format_K_Omega;
      visco_turb_au_format_K_Omega.typer(type);
      DoubleTab& visco_turb_K_Omega = complete_viscosity_field(n, eqn_transp_K_Omega().domaine_dis().valeur(), visco_turb_au_format_K_Omega);

      if (visco_turb_K_Omega.size() != n)
        {
          Cerr << "visco_turb_K_Omega size is " << visco_turb_K_Omega.size() << " instead of " << n << finl;
          exit();
        }

      // A la fin de cette boucle, le tableau visco_turb_K_Omega contient les valeurs de la viscosite turbulente
      // au centre des faces du maillage.
      fill_turbulent_viscosity_tab(n, tab_K_Omega, visco_turb_K_Omega);

      // On connait donc la viscosite turbulente au centre des faces de chaque element
      // On cherche maintenant a interpoler cette viscosite turbulente au centre des elements.
      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_Omega.valeur());
      Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente visco_turb_au_format_K_Omega", visco_turb_au_format_K_Omega.valeur());
    }
  else
    fill_turbulent_viscosity_tab(n, tab_K_Omega, visco_turb);

  la_viscosite_turbulente_.changer_temps(temps);
  Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente la_viscosite_turbulente after", la_viscosite_turbulente_.valeurs());
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_K_Omega::fill_turbulent_viscosity_tab(const int n, const DoubleTab& tab_K_Omega, DoubleTab& turbulent_viscosity)
{
  for (int i = 0; i < n; ++i)
    {
      if (tab_K_Omega(i, 1) <= OMEGA_MIN_)
        turbulent_viscosity[i] = 0;
      else
        {
          if (model_variant_ == "SST")
            {
              // SST variant from the TurbModel group formulation
              const double tmpmax = std::max(CST_A1 * tab_K_Omega(i, 1), get_enstrophy()[i] * get_fieldF2()[i]);
              turbulent_viscosity[i] = tab_K_Omega(i, 0) * CST_A1 / tmpmax;
            }
          else
            turbulent_viscosity[i] = tab_K_Omega(i, 0) / tab_K_Omega(i, 1); // k/omega
        }
    }
}

/*! @brief Initialise three tabs when turbulence model variant is SST
 *
 *  The tables F1, F2 and enstrophy are used to computer a blending function.
 *  They are initialised with the total (real+virtual) DoF.
 *
 */
void Modele_turbulence_hyd_K_Omega::init_F1_F2_enstrophy()
{
  int const n = K_Omega().valeurs().dimension_tot(0);

  blenderF1_.resize(n);
  fieldF2_.resize(n);
  enstrophy_.resize(n);
}

int Modele_turbulence_hyd_K_Omega::preparer_calcul()
{
  eqn_transp_K_Omega().preparer_calcul();
  Modele_turbulence_hyd_RANS_K_Omega_base::preparer_calcul();
  if (model_variant_ == "SST")
    init_F1_F2_enstrophy();
  calculate_limit_viscosity<MODELE_TYPE::K_OMEGA>(K_Omega(), -123. /* unused */ );
  Debog::verifier("Modele_turbulence_hyd_K_Omega::preparer_calcul la_viscosite_turbulente", la_viscosite_turbulente_.valeurs());
  return 1;
}

bool Modele_turbulence_hyd_K_Omega::initTimeStep(double dt)
{
  return eqn_transport_K_Omega_.initTimeStep(dt);
}

/*! @brief Performs a time update of the turbulence model.
 *
 * Update the transport equation of k and omega, computes the wall function and the turbulence
 * viscosity.
 *
 * @param[in] (double temps) new time
 */
void Modele_turbulence_hyd_K_Omega::mettre_a_jour(double temps)
{
  Schema_Temps_base& sch = eqn_transp_K_Omega().schema_temps();
  eqn_transp_K_Omega().domaine_Cl_dis().mettre_a_jour(temps);
  if (!eqn_transp_K_Omega().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K_Omega());
  eqn_transp_K_Omega().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  Debog::verifier("Modele_turbulence_hyd_K_Omega::mettre_a_jour la_viscosite_turbulente before", la_viscosite_turbulente_.valeurs());
  calculate_limit_viscosity<MODELE_TYPE::K_OMEGA>(K_Omega(), -123. /* unused */ );
  Debog::verifier("Modele_turbulence_hyd_K_Omega::mettre_a_jour la_viscosite_turbulente after", la_viscosite_turbulente_.valeurs());
  statistiques().end_count(nut_counter_);
}
