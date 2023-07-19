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

#include <Modifier_nut_pour_fluide_dilatable.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Modele_turbulence_scal_base.h>
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
#include <Debog.h>
#include <Param.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Omega,
                        "Modele_turbulence_hyd_K_Omega",
                        Mod_turb_hyd_RANS_komega);
// XD k_omega mod_turb_hyd_rans k_omega -1 Turbulence model (k-omega).

/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Modele_turbulence_hyd_K_Omega::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

/*! @brief Simple appel a Mod_turb_hyd_RANS_komega::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Modele_turbulence_hyd_K_Omega::readOn(Entree& s)
{
  return Mod_turb_hyd_RANS_komega::readOn(s);
}

void Modele_turbulence_hyd_K_Omega::set_param(Param& param)
{
  Mod_turb_hyd_RANS_komega::set_param(param);
  param.ajouter("model_variant", &model_variant, Param::OPTIONAL); // XD_ADD_P Nom Model variant for k-omega (default value SST)
  param.ajouter_non_std("Transport_K_Omega", (this), Param::REQUIRED); // XD_ADD_P transport_k_omega Keyword to define the (k-omega) transportation equation.
  param.ajouter("PRANDTL_K", &Prandtl_K);
  param.ajouter("PRANDTL_Omega", &Prandtl_Omega);
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
    return Mod_turb_hyd_RANS_komega::lire_motcle_non_standard(mot, is);
  return 1;
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
  DoubleTab& visco_turb = la_viscosite_turbulente.valeurs();

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
  Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente la_viscosite_turbulente before", la_viscosite_turbulente.valeurs());
  Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente tab_K_Omega", tab_K_Omega);
  if (visco_turb.size() == n)
    non_prepare = 0.;
  non_prepare = mp_max(non_prepare);

  if (non_prepare == 1)
    {
      Champ_Inc visco_turb_au_format_K_Omega;

      visco_turb_au_format_K_Omega.typer(type);
      Champ_Inc_base& ch_visco_turb_K_Omega = visco_turb_au_format_K_Omega.valeur();
      ch_visco_turb_K_Omega.associer_domaine_dis_base(eqn_transp_K_Omega().domaine_dis().valeur());
      ch_visco_turb_K_Omega.nommer("diffusivite_turbulente");
      ch_visco_turb_K_Omega.fixer_nb_comp(1);
      ch_visco_turb_K_Omega.fixer_nb_valeurs_nodales(n);
      ch_visco_turb_K_Omega.fixer_unite("inconnue");
      ch_visco_turb_K_Omega.changer_temps(0.);
      DoubleTab& visco_turb_K_Omega = ch_visco_turb_K_Omega.valeurs();

      if(visco_turb_K_Omega.size() != n)
        {
          Cerr << "visco_turb_K_Omega size is " << visco_turb_K_Omega.size()
               << " instead of " << n << finl;
          exit();
        }

      // A la fin de cette boucle, le tableau visco_turb_K_Omega
      // contient les valeurs de la viscosite turbulente
      // au centre des faces du maillage.
      // Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente visco_turb_K_Omega before",visco_turb_K_Omega);
      fill_turbulent_viscosity_tab(tab_K_Omega, visco_turb_K_Omega);
      // Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente visco_turb_K_Omega after",visco_turb_K_Omega);

      // On connait donc la viscosite turbulente au centre des faces de chaque element
      // On cherche maintenant a interpoler cette viscosite turbulente au centre des
      // elements.
      la_viscosite_turbulente->affecter(visco_turb_au_format_K_Omega.valeur());
      Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente visco_turb_au_format_K_Omega", visco_turb_au_format_K_Omega.valeur());
    }
  else
    {
      fill_turbulent_viscosity_tab(tab_K_Omega, visco_turb);
    }

  la_viscosite_turbulente.changer_temps(temps);
  Debog::verifier("Modele_turbulence_hyd_K_Omega::calculer_viscosite_turbulente la_viscosite_turbulente after", la_viscosite_turbulente.valeurs());
  return la_viscosite_turbulente;
}

void Modele_turbulence_hyd_K_Omega::fill_turbulent_viscosity_tab(const DoubleTab& tab_K_Omega, DoubleTab& turbulent_viscosity)
{
  const int nb_elem = tab_K_Omega.dimension(0);
  for (int i = 0; i < nb_elem; ++i)
    {
      if (tab_K_Omega(i, 1) <= OMEGA_MIN)
        turbulent_viscosity[i] = 0;
      else
        turbulent_viscosity[i] = tab_K_Omega(i, 0)/tab_K_Omega(i, 1); // k/omega
    }
}

void imprimer_evolution_komega(const Champ_Inc& le_champ_K_Omega, const Schema_Temps_base& sch, int avant)
{
  if (sch.nb_pas_dt() == 0 || sch.limpr())
    {
      const DoubleTab& K_Omega = le_champ_K_Omega.valeurs();
      double k_min = DMAXFLOAT, omega_min = DMAXFLOAT, nut_min = DMAXFLOAT;
      double k_max = 0, omega_max = 0, nut_max = 0;
      int loc_k_min = -1, loc_omega_min = -1, loc_nut_min = -1;
      int loc_k_max = -1, loc_omega_max = -1, loc_nut_max = -1;
      int size = K_Omega.dimension(0);

      if (size < 0)
        {
          if (!sub_type(Champ_Inc_P0_base, le_champ_K_Omega.valeur()))
            {
              Cerr << "Unsupported K_Omega field in Modele_turbulence_hyd_K_Omega::imprimer_evolution_komega()" << finl;
              Process::exit(-1);
            }
          size = le_champ_K_Omega.valeur().equation().domaine_dis().domaine().nb_elem();
        }

      ConstDoubleTab_parts parts(le_champ_K_Omega.valeurs());
      for (int n = 0; n < size; n++)
        {
          const double k = K_Omega(n, 0);
          const double omega = K_Omega(n, 1);
          double nut = 0;

          if (omega > 0) nut = k/omega;

          if (k < k_min)
            {
              k_min = k;
              loc_k_min = n;
            }
          else if (k > k_max)
            {
              k_max = k;
              loc_k_max = n;
            }

          if (omega < omega_min)
            {
              omega_min = omega;
              loc_omega_min = n;
            }
          else if (omega > omega_max)
            {
              omega_max = omega;
              loc_omega_max = n;
            }

          if (nut < nut_min)
            {
              nut_min = nut;
              loc_nut_min = n;
            }
          else if (nut > nut_max)
            {
              nut_max = nut;
              loc_nut_max = n;
            }
        }
      /*
      k_min = Process::mp_min(k_min);
      eps_min = Process::mp_min(eps_min);
      nut_min = Process::mp_min(nut_min);
      k_max = Process::mp_max(k_max);
      eps_max = Process::mp_max(eps_max);
      nut_max = Process::mp_max(nut_max);
       */
      ArrOfDouble values(3);

      // min values
      values[0] = k_min;
      values[1] = omega_min;
      values[2] = nut_min;
      mp_min_for_each_item(values);
      k_min = values[0];
      omega_min = values[1];
      nut_min = values[2];

      // max values
      values[0] = k_max;
      values[1] = omega_max;
      values[2] = nut_max;
      mp_max_for_each_item(values);
      k_max = values[0];
      omega_max = values[1];
      nut_max = values[2];

      // ecriture
      if (Process::je_suis_maitre())
        {
          Cout << finl << "K_Omega evolution (" << (avant?"before":"after") << " law of the wall applies) at time " << le_champ_K_Omega.temps() << ":" << finl;
          Cout << "std::min(k)=" << k_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_k_min;
          Cout << finl;
          Cout << "std::min(omega)=" << omega_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_omega_min;
          Cout << finl;
          Cout << "std::min(nut)=" << nut_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_nut_min;
          Cout << finl;
          Cout << "std::max(k)=" << k_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_k_max;
          Cout << finl;
          Cout << "std::max(omega)=" << omega_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_omega_max;
          Cout << finl;
          Cout << "std::max(nut)=" << nut_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_nut_max;
          Cout << finl;
        }
    }
}

int Modele_turbulence_hyd_K_Omega::preparer_calcul()
{
  eqn_transp_K_Omega().preparer_calcul();
  Mod_turb_hyd_RANS_komega::preparer_calcul();

  // GF quand on demarre un calcul il est bon d'utliser la ldp
  // encore plus quand on fait une reprise !!!!!!!!
  Champ_Inc& ch_K_Omega = K_Omega();

  const Milieu_base& mil = equation().probleme().milieu();
  if (equation().probleme().is_dilatable())
    diviser_par_rho_si_dilatable(ch_K_Omega.valeurs(), mil);
  imprimer_evolution_komega(ch_K_Omega, eqn_transp_K_Omega().schema_temps(), 1);
  loipar.calculer_hyd(ch_K_Omega);
  eqn_transp_K_Omega().controler_K_Omega();
  calculer_viscosite_turbulente(ch_K_Omega.temps());
  limiter_viscosite_turbulente();
  // on remultiplie K_Omega par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Omega.valeurs(), mil);
      // correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Omega::preparer_calcul la_viscosite_turbulente", la_viscosite_turbulente.valeurs());
  imprimer_evolution_komega(ch_K_Omega, eqn_transp_K_Omega().schema_temps(), 0);
  return 1;
}

bool Modele_turbulence_hyd_K_Omega::initTimeStep(double dt)
{
  return eqn_transport_K_Omega.initTimeStep(dt);
}

/*! @brief Effectue une mise a jour en temps du modele de turbulence.
 *
 * Met a jour l'equation de transport K-Omega,
 *     calcule la loi de paroi et la viscosite turbulente
 *     au nouveau temps.
 *
 * @param (double temps) le temps de mise a jour
 */
void Modele_turbulence_hyd_K_Omega::mettre_a_jour(double temps)
{
  Champ_Inc& ch_K_Omega = K_Omega();
  Schema_Temps_base& sch = eqn_transp_K_Omega().schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base

  eqn_transp_K_Omega().domaine_Cl_dis().mettre_a_jour(temps);
  if (!eqn_transp_K_Omega().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K_Omega());
  eqn_transp_K_Omega().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  // cAlan : Mutualisable avec preparer ?
  const Milieu_base& mil = equation().probleme().milieu();
  Debog::verifier("Modele_turbulence_hyd_K_Omega::mettre_a_jour la_viscosite_turbulente before", la_viscosite_turbulente.valeurs());
  // on divise K_Omega par rho en QC pour revenir a K et Omega
  if (equation().probleme().is_dilatable())
    diviser_par_rho_si_dilatable(ch_K_Omega.valeurs(), mil);
  imprimer_evolution_komega(ch_K_Omega, eqn_transp_K_Omega().schema_temps(), 1);
  Debog::verifier("Modele_turbulence_hyd_K_Omega::mettre_a_jour loiparoi before", 0);
  loipar.calculer_hyd(ch_K_Omega);
  Debog::verifier("Modele_turbulence_hyd_K_Omega::mettre_a_jour loiparoi after", 1);
  eqn_transp_K_Omega().controler_K_Omega();
  calculer_viscosite_turbulente(ch_K_Omega.temps());
  limiter_viscosite_turbulente();
  // on remultiplie K_Omega par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Omega.valeurs(),mil);
      // correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Omega::mettre_a_jour la_viscosite_turbulente after", la_viscosite_turbulente.valeurs());
  imprimer_evolution_komega(ch_K_Omega, eqn_transp_K_Omega().schema_temps(), 0);
  statistiques().end_count(nut_counter_);
}

// cAlan : fonction pour le bicephale
// const Equation_base& Modele_turbulence_hyd_K_Omega::equation_k_omega(int i) const
// {
//   assert ((i==0));
//   return eqn_transport_K_Omega;
// }

// cAlan : templatable avec K_Eps
const Champ_base& Modele_turbulence_hyd_K_Omega::get_champ(const Motcle& nom) const
{
  try
    {
      return Mod_turb_hyd_RANS_komega::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }
  throw Champs_compris_erreur();
}

void Modele_turbulence_hyd_K_Omega::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Mod_turb_hyd_RANS_komega::get_noms_champs_postraitables(nom, opt);
}

void Modele_turbulence_hyd_K_Omega::verifie_loi_paroi()
{
  Nom lp = loipar.valeur().que_suis_je();
  Cerr << "We should probably do something here... Should be negligeable, isn't it?" << finl;
}
