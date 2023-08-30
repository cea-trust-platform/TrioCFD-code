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
#include <Schema_Temps.h>
#include <stat_counters.h>
#include <Param.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_2_Couches,"Modele_turbulence_hyd_K_Epsilon_2_Couches",Mod_turb_hyd_RANS);


/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Modele_turbulence_hyd_K_Eps_2_Couches::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}


/*! @brief Simple appel a Mod_turb_hyd_RANS::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Modele_turbulence_hyd_K_Eps_2_Couches::readOn(Entree& s )
{
  return Mod_turb_hyd_RANS::readOn(s);
}

void Modele_turbulence_hyd_K_Eps_2_Couches::set_param(Param& param)
{
  Mod_turb_hyd_RANS::set_param(param);
  param.ajouter_non_std("Transport_K_KEpsilon",(this),Param::REQUIRED);
}


int Modele_turbulence_hyd_K_Eps_2_Couches::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="Transport_K_KEpsilon")
    {
      eqn_transp_K_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_K_Eps();
      return 1;
    }
  else
    return Mod_turb_hyd_RANS::lire_motcle_non_standard(mot,is);
}

/*! @brief Calcule la viscosite turbulente au temps demande.
 *
 * @param (double temps) le temps auquel il faut calculer la viscosite
 * @return (Champ_Fonc&) la viscosite turbulente au temps demande
 * @throws erreur de taille de visco_turb_K_eps
 */
Champ_Fonc& Modele_turbulence_hyd_K_Eps_2_Couches::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK_Eps=eqn_transport_K_Eps.inconnue().valeur();
  Nom type=chK_Eps.que_suis_je();
  const DoubleTab& tab_K_Eps = chK_Eps.valeurs();
  DoubleTab& visco_turb =  la_viscosite_turbulente.valeurs();

  // K_Eps(i,0) = K au noeud i
  // K_Eps(i,1) = Epsilon au noeud i

  int n = tab_K_Eps.dimension(0);
  if (visco_turb.size() != n)
    {
      Champ_Inc visco_turb_au_format_K_eps;

      visco_turb_au_format_K_eps.typer(type);
      Champ_Inc_base& ch_visco_turb_K_eps=visco_turb_au_format_K_eps.valeur();
      ch_visco_turb_K_eps.associer_domaine_dis_base(eqn_transport_K_Eps.domaine_dis().valeur());
      ch_visco_turb_K_eps.nommer("diffusivite_turbulente");
      ch_visco_turb_K_eps.fixer_nb_comp(1);
      ch_visco_turb_K_eps.fixer_nb_valeurs_nodales(n);
      ch_visco_turb_K_eps.fixer_unite("inconnue");
      ch_visco_turb_K_eps.changer_temps(0.);

      DoubleTab& visco_turb_K_eps =  ch_visco_turb_K_eps.valeurs();

      if(visco_turb_K_eps.size() != n)
        {
          Cerr << "visco_turb_K_eps size is " << visco_turb_K_eps.size()
               << " instead of " << n << finl;
          exit();
        }
      for (int i=0; i<n; i++)
        {
          if (tab_K_Eps(i,1) <= LeEPS_MIN)
            visco_turb_K_eps[i] = 0;
          else
            visco_turb_K_eps[i] = CMU*tab_K_Eps(i,0)*tab_K_Eps(i,0)/tab_K_Eps(i,1);
        }
      la_viscosite_turbulente->affecter(visco_turb_au_format_K_eps.valeur());
    }
  else
    {
      for (int i=0; i<n; i++)
        {
          if (tab_K_Eps(i,1) <= LeEPS_MIN)
            visco_turb[i] = 0;
          else
            visco_turb[i] = CMU*tab_K_Eps(i,0)*tab_K_Eps(i,0)/tab_K_Eps(i,1);
        }
    }
  la_viscosite_turbulente.changer_temps(temps);
  return la_viscosite_turbulente;
}

int Modele_turbulence_hyd_K_Eps_2_Couches::preparer_calcul()
{
  eqn_transp_K_Eps().preparer_calcul();
  calculer_viscosite_turbulente(equation().schema_temps().temps_courant());
  Modele_turbulence_hyd_base::preparer_calcul();
  calculer_viscosite_turbulente(K_Eps().temps());
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
  return 1;
}

bool Modele_turbulence_hyd_K_Eps_2_Couches::initTimeStep(double dt)
{
  return eqn_transport_K_Eps.initTimeStep(dt);
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
  DoubleTab& visco_turb =  la_viscosite_turbulente.valeurs();
  Champ_Inc& ch_K_Eps = K_Eps();
  Schema_Temps_base& sch = eqn_transport_K_Eps.schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base
  eqn_transport_K_Eps.domaine_Cl_dis().mettre_a_jour(temps);
  sch.faire_un_pas_de_temps_eqn_base(eqn_transport_K_Eps);
  eqn_transport_K_Eps.mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  eqn_transport_K_Eps.controler_K_Eps();
  const Milieu_base& mil=equation().probleme().milieu();
  // on divise K_eps par rho en QC pour revenir a K et Eps
  if (equation().probleme().is_dilatable()) diviser_par_rho_si_dilatable(ch_K_Eps.valeurs(),mil);
  calculer_viscosite_turbulente(ch_K_Eps.temps());
  loipar->calculer_hyd(visco_turb,ch_K_Eps.valeurs());
  limiter_viscosite_turbulente();
  // on remultiplie K_eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Eps.valeurs(),mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente.valeurs().echange_espace_virtuel();
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

const Equation_base& Modele_turbulence_hyd_K_Eps_2_Couches::equation_k_eps(int i) const
{
  assert ((i==0));
  return eqn_transport_K_Eps;

}
