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
// File:        Transport_2eq_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_2eq_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <communications.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Zone_VF.h>
#include <Param.h>
#include <Debog.h>

Implemente_base_sans_constructeur(Transport_2eq_base, "Transport_2eq_base", Equation_base);

Transport_2eq_base::Transport_2eq_base()
{
  champs_compris_.ajoute_nom_compris("2eq_residu");
}

/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_2eq_base::printOn(Sortie& is) const
{
  return is << que_suis_je() << "\n";
}

/*! @brief Simple appel a Equation_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_2eq_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  return is;
}

void Transport_2eq_base::set_param(Param& param)
{
  Equation_base::set_param(param);
  param.ajouter_non_std("diffusion", (this));
  param.ajouter_non_std("convection", (this));
  param.ajouter_condition("is_read_diffusion", "The diffusion operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_condition("is_read_convection", "The convection operator must be read, select negligeable type if you want to neglect it.");
}

/*! @brief Associe un milieu physique a l'equation.
 *
 * @param (Milieu_base& un_milieu) le milieu physique a associer a l'equation
 */
void Transport_2eq_base::associer_milieu_base(const Milieu_base& un_milieu)
{
  le_fluide =  un_milieu;
}

/*! @brief Renvoie le milieu (fluide) associe a l'equation.
 *
 * @return (Milieu_base&) le milieu (fluide) associe a l'equation
 */
Milieu_base& Transport_2eq_base::milieu()
{
  if(!le_fluide.non_nul())
    {
      Cerr << "No fluid has been associated to the two equations Transport"
           << que_suis_je() << " equation." << finl;
      exit();
    }
  return le_fluide.valeur();
}

/*! @brief Renvoie le milieu (fluide) associe a l'equation.
 *
 * (version const)
 *
 * @return (Milieu_base&) le milieu (fluide) associe a l'equation
 */
const Milieu_base& Transport_2eq_base::milieu() const
{
  if(!le_fluide.non_nul())
    {
      Cerr << "No fluid has been associated to the Transport K_Epsilon"
           << que_suis_je() << " equation." << finl;
      exit();
    }
  return le_fluide.valeur();
}

void Transport_2eq_base::associer(const Equation_base& eqn_hydr)
{
  Equation_base::associer_pb_base(eqn_hydr.probleme());
  Equation_base::associer_sch_tps_base(eqn_hydr.schema_temps());
  Equation_base::associer_zone_dis(eqn_hydr.zone_dis());
}

double Transport_2eq_base::calculer_pas_de_temps() const
{
  // on prend le pas de temps de l'eq de NS.
  return probleme().equation(0).calculer_pas_de_temps();
}
