/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Convection_diffusion_turbulence_multiphase.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Equations
// Version:     /main/52
//
//////////////////////////////////////////////////////////////////////////////

#include <Convection_diffusion_turbulence_multiphase.h>
#include <Pb_Multiphase.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Avanc.h>
#include <Debog.h>
#include <Frontiere_dis_base.h>
#include <EcritureLectureSpecial.h>
#include <Champ_Uniforme.h>
#include <Matrice_Morse.h>
#include <Navier_Stokes_std.h>
#include <TRUSTTrav.h>
#include <Neumann_sortie_libre.h>
#include <Op_Conv_negligeable.h>
#include <Param.h>
#include <Schema_Implicite_base.h>
#include <SETS.h>
#include <EChaine.h>
#include <Scalaire_impose_paroi.h>
#include <Echange_global_impose.h>

#define old_forme

Implemente_base(Convection_diffusion_turbulence_multiphase,"Convection_diffusion_turbulence_multiphase",Convection_Diffusion_std);

/*! @brief Simple appel a: Convection_Diffusion_std::printOn(Sortie&)
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Convection_diffusion_turbulence_multiphase::printOn(Sortie& is) const
{
  return Convection_Diffusion_std::printOn(is);
}

/*! @brief Verifie si l'equation a une inconnue et un fluide associe et appelle Convection_Diffusion_std::readOn(Entree&).
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree& is) le flot d'entree modifie
 */
Entree& Convection_diffusion_turbulence_multiphase::readOn(Entree& is)
{
  assert(l_inco_ch.non_nul());
  assert(le_fluide.non_nul());
  Convection_Diffusion_std::readOn(is);
  return is;
}

void Convection_diffusion_turbulence_multiphase::completer()
{
  Convection_Diffusion_std::completer(); // en fait c'est Equation_base::completer() mais on sais pas un jour ...

  const Domaine_dis_base& zdis = domaine_dis();
  if (zdis.que_suis_je().debute_par("Domaine_VDF"))
    {
      // initialiser l'operateur grad SI VDF
      Op_Grad_.associer_eqn(*this);
      Op_Grad_.typer();
      Op_Grad_.l_op_base().associer_eqn(*this);
      const Domaine_Cl_dis_base& zcl = domaine_Cl_dis();
      const Champ_Inc& inco = inconnue();
      Op_Grad_->associer(zdis, zcl, inco);
    }
}

/*! @brief Associe un milieu physique a l'equation, le milieu est en fait caste en Fluide_base ou en Fluide_Ostwald.
 *
 * @param (Milieu_base& un_milieu)
 * @throws les proprietes physiques du fluide ne sont pas toutes specifiees
 */
void Convection_diffusion_turbulence_multiphase::associer_milieu_base(const Milieu_base& un_milieu)
{
  le_fluide = ref_cast(Fluide_base,un_milieu);
}

/*! @brief Renvoie le milieu physique de l'equation.
 *
 * (un Fluide_base upcaste en Milieu_base)
 *     (version const)
 *
 * @return (Milieu_base&) le Fluide_base upcaste en Milieu_base
 */
const Milieu_base& Convection_diffusion_turbulence_multiphase::milieu() const
{
  return le_fluide.valeur();
}


/*! @brief Renvoie le milieu physique de l'equation.
 *
 * (un Fluide_base upcaste en Milieu_base)
 *
 * @return (Milieu_base&) le Fluide_base upcaste en Milieu_base
 */
Milieu_base& Convection_diffusion_turbulence_multiphase::milieu()
{
  return le_fluide.valeur();
}

/*! @brief Impression des flux sur les bords sur un flot de sortie.
 *
 * Appelle Equation_base::impr(Sortie&)
 *
 * @param (Sortie& os) un flot de sortie
 * @return (int) code de retour propage
 */
int Convection_diffusion_turbulence_multiphase::impr(Sortie& os) const
{
  return Equation_base::impr(os);
}

/*! @brief Renvoie le nom du domaine d'application de l'equation.
 *
 * Ici "Turbulence".
 *
 * @return (Motcle&) le nom du domaine d'application de l'equation
 */
const Motcle& Convection_diffusion_turbulence_multiphase::domaine_application() const
{
  static Motcle mot("Turbulence");
  return mot;
}
