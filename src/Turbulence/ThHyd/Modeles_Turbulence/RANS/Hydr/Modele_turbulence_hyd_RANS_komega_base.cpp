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
// File:        Modele_turbulence_hyd_RANS_komega_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_RANS_komega_base.h>
#include <Transport_K_Omega_base.h>
#include <TRUSTTrav.h>
#include <Param.h>

Implemente_base_sans_constructeur(Modele_turbulence_hyd_RANS_komega_base,
                                  "Modele_turbulence_hyd_RANS_komega_base",
                                  Modele_turbulence_hyd_2_eq_base);
// XD mod_turb_hyd_rans_komega modele_turbulence_hyd_deriv mod_turb_hyd_rans_komega -1 Class for RANS turbulence model for Navier-Stokes equations.

Modele_turbulence_hyd_RANS_komega_base::Modele_turbulence_hyd_RANS_komega_base()
{
  Prandtl_K = 2.; // cAlan: careful, it is the inverse of the classical definition for code purpose
  Prandtl_Omega = 2.;
  OMEGA_MIN = 1.e-20; // OMEGA_MIN cannot be too small?
  OMEGA_MAX = 1.e+10;
  K_MIN = 1.e-20;
  model_variant = "SST";
  lquiet = 0;
}
/*! @brief Simple appel a Modele_turbulence_hyd_base::printOn(Sortie&)
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Modele_turbulence_hyd_RANS_komega_base::printOn(Sortie& is) const
{
  return Modele_turbulence_hyd_2_eq_base::printOn(is);
}


/*! @brief Simple appel a Modele_turbulence_hyd_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Modele_turbulence_hyd_RANS_komega_base::readOn(Entree& is)
{
  Modele_turbulence_hyd_2_eq_base::readOn(is);
  return is;
}

void Modele_turbulence_hyd_RANS_komega_base::set_param(Param& param)
{
  Modele_turbulence_hyd_base::set_param(param);
  param.ajouter("omega_min", &OMEGA_MIN); // XD_ADD_P double Lower limitation of omega (default value 1.e-10).
  param.ajouter("omega_max", &OMEGA_MAX); // XD_ADD_P double Upper limitation of omega (default value 1.e+10).
  param.ajouter("k_min", &K_MIN); // XD_ADD_P double Lower limitation of k (default value 1.e-10).
  param.ajouter_flag("quiet", &lquiet); // XD_ADD_P flag To disable printing of information about k and omega.
}

/*! @brief
 *
 */
void Modele_turbulence_hyd_RANS_komega_base::completer()
{
  eqn_transp_K_Omega().completer();
  verifie_loi_paroi();
}

void Modele_turbulence_hyd_RANS_komega_base::verifie_loi_paroi()
{
  Nom lp = loipar.valeur().que_suis_je();
  if (lp == "negligeable_VEF" || lp == "negligeable_VDF")
    {
      Cerr << "The turbulence model of type " << que_suis_je() << finl;
      Cerr << "must not be used with a wall law of type negligeable." << finl;
      Cerr << "Another wall law must be selected with this kind of turbulence model." << finl;
      exit();
    }
}

const Champ_base& Modele_turbulence_hyd_RANS_komega_base::get_champ(const Motcle& nom) const
{
  try
    {
      return Modele_turbulence_hyd_base::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }
  int nb_eq = nombre_d_equations();
  for (int i = 0; i < nb_eq; i++)
    {
      try
        {
          return equation_k_omega(i).get_champ(nom);
        }
      catch (Champs_compris_erreur)
        {
        }
    }
  throw Champs_compris_erreur();
  REF(Champ_base) ref_champ;
  return ref_champ;
}

void Modele_turbulence_hyd_RANS_komega_base::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_base::get_noms_champs_postraitables(nom, opt);

  for (int i = 0; i < nombre_d_equations(); i++)
    equation_k_omega(i).get_noms_champs_postraitables(nom, opt);
}

/*! @brief Sauvegarde le modele de turbulence sur un flot de sortie.
 *
 * (en vue d'une reprise)
 *     Sauvegarde le type de l'objet et
 *     l'equation de transport K-Omega associee.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (int) code de retour propage de: Transport_K_Omega::sauvegarder(Sortie&)
 */
int Modele_turbulence_hyd_RANS_komega_base::sauvegarder(Sortie& os) const
{
  Modele_turbulence_hyd_base::sauvegarder(os);
  return eqn_transp_K_Omega().sauvegarder(os);
}

/*! @brief Reprise du modele a partir d'un flot d'entree.
 *
 * Si l'equation portee par l'objet est non nulle
 *     on effectue une reprise "bidon".
 *
 * @param (Entree& is) un flot d'entree
 * @return (int) code de retour propage de: Transport_K_Omega::sauvegarder(Sortie&) ou 1 si la reprise est bidon.
 */
int Modele_turbulence_hyd_RANS_komega_base::reprendre(Entree& is)
{
  Modele_turbulence_hyd_base::reprendre(is);
  if (mon_equation.non_nul())
    {
      return eqn_transp_K_Omega().reprendre(is);
    }
  else
    {
      double dbidon;
      Nom bidon;
      DoubleTrav tab_bidon;
      is >> bidon >> bidon;
      is >> dbidon;
      tab_bidon.jump(is);
      return 1;
    }
}
