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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Pb_Hydraulique_Turbulent_ALE.cpp
// Directory : $ALE_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Hydraulique_Turbulent_ALE.h>
#include <Fluide_Incompressible.h>
#include <Verif_Cl.h>
#include <Verif_Cl_Turb.h>
#include <Les_mod_turb.h>
#include <Modele_turbulence_hyd_RANS_keps_base.h>
#include <Schema_Implicite_base.h>
#include <Implicite_ALE.h>
#include <Schema_Euler_explicite.h>
#include <Schema_Euler_explicite_ALE.h>
#include <Schema_Euler_Implicite.h>

Implemente_instanciable(Pb_Hydraulique_Turbulent_ALE,"Pb_Hydraulique_Turbulent_ALE",Pb_Fluide_base);
// XD Pb_Hydraulique_Turbulent_ALE Pb_base Pb_Hydraulique_Turbulent_ALE -1 Resolution of hydraulic turbulent problems for ALE
// XD  attr fluide_incompressible fluide_incompressible fluide_incompressible 0 The fluid medium associated with the problem.
// XD  attr Navier_Stokes_Turbulent_ALE Navier_Stokes_Turbulent_ALE Navier_Stokes_Turbulent_ALE 0 Navier-Stokes_ALE equations as well as the associated turbulence model equations on mobile domain (ALE)



/*! @brief Simple appel a: Pb_Fluide_base::printOn(Sortie&) Ecrit le probleme sur un flot de sortie.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Pb_Hydraulique_Turbulent_ALE::printOn(Sortie& os) const
{
  return Pb_Fluide_base::printOn(os);
}


/*! @brief Simple appel a: Pb_Fluide_base::readOn(Entree&) Lit le probleme a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Pb_Hydraulique_Turbulent_ALE::readOn(Entree& is)
{
  return Pb_Fluide_base::readOn(is);
}

/*! @brief Renvoie le nombre d'equation, Renvoie 1 car il y a seulement 1 equation a un probleme
 *
 *     hydraulique turbulent: l'equation de Navier Stokes avec Turbulence
 *
 * @return (int) le nombre d'equation
 */
int Pb_Hydraulique_Turbulent_ALE::nombre_d_equations() const
{
  return 1;
}

/*! @brief Renvoie l' equation d'hydraulique de type Navier_Stokes_Turbulent si i=0 sort (exit) sinon.
 *
 *     (version const)
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation d'hydraulique de type Navier_Stokes_Turbulent
 */
const Equation_base& Pb_Hydraulique_Turbulent_ALE::equation(int i) const
{
  if ( !( i==0 ) )
    {
      Cerr << "\nError in Pb_Hydraulique_Turbulent_ALE::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}

/*! @brief Renvoie l' equation d'hydraulique de type Navier_Stokes_Turbulent si i=0 sort (exit) sinon.
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation d'hydraulique de type Navier_Stokes_Turbulent
 */
Equation_base& Pb_Hydraulique_Turbulent_ALE::equation(int i)
{
  if ( !( i==0 ) )
    {
      Cerr << "\nError in Pb_Hydraulique_Turbulent_ALE::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}



/*! @brief Associe le milieu au probleme.
 *
 * Le milieu doit etre de type fluide incompressible.
 *
 * @param (Milieu_base& mil) le milieu physique a associer au probleme
 * @throws le milieu n'est pas du type Fluide_Incompressible
 */
void Pb_Hydraulique_Turbulent_ALE::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_Incompressible,mil))
    eq_hydraulique.associer_milieu_base(mil);
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je() << " ne peut etre associe a " << finl;
      Cerr << "un probleme de type Pb_Hydraulique_Turbulent_ALE " << finl;
      exit();
    }
}

/*! @brief Teste la compatibilite des equations de convection-diffusion et de l'hydraulique.
 *
 * Le test se fait sur les conditions
 *     aux limites discretisees de chaque equation.
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_turb(const Domaine_Cl_dis&)
 *
 * @return (int) code de retour propage
 */
int Pb_Hydraulique_Turbulent_ALE::verifier()
{
  const Domaine_Cl_dis& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();

  if ( sub_type(Modele_turbulence_hyd_RANS_keps_base, eq_hydraulique.get_modele(TURBULENCE).valeur() ))
    {
      const Modele_turbulence_hyd_RANS_keps_base& le_mod_RANS = ref_cast(Modele_turbulence_hyd_RANS_keps_base, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }

  return 1;
}

void Pb_Hydraulique_Turbulent_ALE::associer_sch_tps_base(const Schema_Temps_base& sch)
{
  Probleme_base::associer_sch_tps_base(sch);

  // Verify that user choosed adapted time scheme/solver
  if (sub_type(Schema_Euler_explicite,le_schema_en_temps_.valeur()))
    {
      if(!sub_type(Schema_Euler_explicite_ALE,le_schema_en_temps_.valeur()))
        {
          Cerr <<"Error: for Mobile domain (Domaine_ALE):  replace  " <<sch.que_suis_je()<<" with "<< sch.que_suis_je() <<"_ALE and restart!"<< finl;
          Process::exit();
        }
    }
  else if (sub_type(Schema_Euler_Implicite,le_schema_en_temps_.valeur()))
    {
      Schema_Euler_Implicite& sch_imp = ref_cast(Schema_Euler_Implicite,le_schema_en_temps_.valeur());
      if (!sub_type(Implicite_ALE,sch_imp.solveur().valeur()))
        {
          Cerr <<"Error: for Mobile domain (Domaine_ALE), you can only use Scheme_euler_implicit time scheme with solveur implicite_ALE "<<finl;
          Cerr <<"Replace  solveur implicite with solveur implicite_ALE and restart!"<< finl;
          Process::exit();
        }

    }
}
