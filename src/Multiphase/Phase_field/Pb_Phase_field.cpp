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
// File:        Pb_Phase_field.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Phase_field/src
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Pb_Phase_field.h>
#include <Fluide_Incompressible.h>
#include <Constituant.h>
#include <Verif_Cl.h>

Implemente_instanciable(Pb_Phase_field,"Pb_Phase_field",Pb_Fluide_base);
// XD pb_phase_field Pb_base pb_phase_field -1 Problem to solve local instantaneous incompressible-two-phase-flows. Complete description of the Phase Field model for incompressible and immiscible fluids can be found into this PDF: TRUST_ROOT/doc/TRUST/phase_field_non_miscible_manuel.pdf
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible 0 The fluid medium associated with the problem.
// XD attr constituant constituant constituant 1 Constituents.
// XD attr navier_stokes_phase_field navier_stokes_phase_field navier_stokes_phase_field 1 Navier Stokes equation for the Phase Field problem.
// XD attr convection_diffusion_phase_field convection_diffusion_phase_field convection_diffusion_phase_field 1 Cahn-Hilliard equation of the Phase Field problem. The unknown of this equation is the concentration C.

/*! @brief Simple appel a: Pb_Fluide_base::printOn(Sortie&) Ecrit le probleme sur un flot de sortie.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Pb_Phase_field::printOn(Sortie& os) const
{
  return Pb_Fluide_base::printOn(os);
}


/*! @brief Simple appel a: Pb_Fluide_base::readOn(Entree&) Lit le probleme a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Pb_Phase_field::readOn(Entree& is)
{
  return Pb_Fluide_base::readOn(is);
}

/*! @brief Renvoie le nombre d'equation, Renvoie 2 car il y a 2 equations a un probleme
 *
 *     hydraulique avec transport:
 *       - l'equation de Navier Stokes
 *       - une equation de convection-diffusion (eventuellement vectorielle)
 *
 * @return (int) le nombre d'equation
 */
int Pb_Phase_field::nombre_d_equations() const
{
  return 2;
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_std si i=0 Renvoie l'equation de convection-diffusion de type
 *
 *     Convection_Diffusion_Concentration si i=1
 *     (l'equation de convection diffusion peut-etre vectorielle)
 *     (version const)
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation correspondante a l'index
 */
const Equation_base& Pb_Phase_field::equation(int i) const
{
  assert ((i==0) || (i==1));
  if (i == 0)
    {
      return eq_hydraulique;
    }
  else
    {
      return eq_concentration;
    }
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes si i=0 Renvoie l'equation de convection-diffusion de type
 *
 *     Convection_Diffusion_Concentration si i=1
 *     (l'equation de convection diffusion peut-etre vectorielle)
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation correspondante a l'index
 */
Equation_base& Pb_Phase_field::equation(int i)
{
  assert ((i==0) || (i==1));
  if (i == 0)
    {
      return eq_hydraulique;
    }
  else
    {
      return eq_concentration;
    }
}


/*! @brief Associe un milieu au probleme, Si le milieu est de type
 *
 *       - Fluide_Incompressible, il sera associe a l'equation de l'hydraulique
 *       - Constituant, il sera associe a l'equation de convection-diffusion
 *     Un autre type de milieu provoque une erreur
 *
 * @param (Milieu_base& mil) le milieu physique a associer au probleme
 * @throws mauvais type de milieu physique
 */
void Pb_Phase_field::associer_milieu_base(const Milieu_base& mil)
{
  if ( sub_type(Fluide_Incompressible,mil) )
    eq_hydraulique.associer_milieu_base(mil);
  else if ( sub_type(Constituant,mil) )
    eq_concentration.associer_milieu_base(mil);
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je() << " ne peut etre associe a " << finl;
      Cerr << "un probleme de type Pb_Phase_field " << finl;
      exit();
    }
}

void Pb_Phase_field::typer_lire_milieu(Entree& is)
{
  const int nb_milieu = 2;
  le_milieu_.resize(nb_milieu);
  for (int i = 0; i < nb_milieu; i++)
    {
      is >> le_milieu_[i]; // On commence par la lecture du milieu
      associer_milieu_base(le_milieu_[i].valeur()); // On l'associe a chaque equations (methode virtuelle pour chaque pb ...)
    }

  Probleme_base::discretiser_equations();

  // remontee de l'inconnue vers le milieu
  for (int i = 0; i < nombre_d_equations(); i++)
    {
      equation(i).associer_milieu_equation();
      equation(i).milieu().discretiser((*this), la_discretisation_.valeur());
    }
}

/*! @brief Teste la compatibilite des equations de convection-diffusion et de l'hydraulique.
 *
 * Le test se fait sur les conditions
 *     aux limites discretisees de chaque equation.
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_concentration(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return (int) code de retour propage
 */
int Pb_Phase_field::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_co = eq_concentration.domaine_Cl_dis();
  return tester_compatibilite_hydr_concentration(domaine_Cl_hydr,domaine_Cl_co);
}
