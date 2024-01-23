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
// File:        Op_Diff_K_Omega_base.cpp
// Directory:   $TURBULENCE_ROOT/src/Kernel/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Omega_base.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Discretisation_base.h>
#include <stat_counters.h>
#include <Champ_Uniforme.h>
#include <Champ_Don.h>

Implemente_base(Op_Diff_K_Omega_base, "Op_Diff_K_Omega_base", Operateur_base);

Implemente_instanciable(Op_Diff_K_Omega_negligeable,
                        "Op_Diff_K_Omega_negligeable",
                        Op_Diff_K_Omega_base);

Implemente_instanciable(Op_Diff_K_Omega,
                        "Op_Diff_K_Omega",
                        DERIV(Op_Diff_K_Omega_base));


/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie& s) le flot de sortie modifie
 */
Sortie& Op_Diff_K_Omega_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie& s) le flot de sortie modifie
 */
Sortie& Op_Diff_K_Omega_negligeable::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie& s) le flot de sortie modifie
 */
Sortie& Op_Diff_K_Omega::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


/*! @brief NE FAIT RIEN A surcharger dans les classes derivees
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree
 */
Entree& Op_Diff_K_Omega_base::readOn(Entree& s )
{
  return s ;
}

/*! @brief NE FAIT RIEN A surcharger dans les classes derivees
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree
 */
Entree& Op_Diff_K_Omega_negligeable::readOn(Entree& s )
{
  return s ;
}

/*! @brief Simple appel a Operateur::lire(Entree&) Lit l'operateur a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Op_Diff_K_Omega::readOn(Entree& s )
{
  Operateur::lire(s);
  return s;
}


/*! @brief Type l'operateur s'associe a son equation
 *
 *     Associe la diffusivite turbulente a l'operateur base.
 *
 */
void Op_Diff_K_Omega::typer()
{
  if (typ=="negligeable")
    {
      DERIV(Op_Diff_K_Omega_base)::typer("Op_Diff_K_Omega_negligeable");
      valeur().associer_diffusivite_turbulente();
    }
  else
    {
      Nom nom_type="Op_Diff_K_Omega_";
      Nom discr = equation().discretisation().que_suis_je();

      // les operateurs de diffusion sont communs aux discretisations VEF et VEFP1B
      if(discr=="VEFPreP1B") discr="VEF";
      //Cerr <<" >>>>>>>>>>>>>>>>>>>>>>>>>>>" <<  equation().que_suis_je() << finl;


      Transport_K_Omega_base& leq = ref_cast(Transport_K_Omega_base,equation());
      Nom qc = leq.modele_turbulence().equation().que_suis_je();
      Cerr << ">>>>>>> Nom eq = " << qc << finl;
      if (qc=="Navier_Stokes_QC" || qc == "Navier_Stokes_Turbulent_QC")
        {
          nom_type+="QC_";
        }

      if (!sub_type(Champ_Uniforme,diffusivite()))
        {
          if (discr!="VEF")
            {
              nom_type+="var_";
            }
        }

      nom_type +=discr;
      Cerr << " >>>>>> Operateur = " << nom_type << finl;
      if (discr!="VDF_Hyper")
        {
          nom_type += "_";
          Nom type_inco=equation().inconnue()->que_suis_je();
          nom_type+=(type_inco.suffix("Champ_"));
          if (axi)
            nom_type += "_Axi";
        }
      DERIV(Op_Diff_K_Omega_base)::typer(nom_type);
      valeur().associer_eqn(equation());
      valeur().associer_diffusivite_turbulente();
      valeur().associer_diffusivite(diffusivite());
      Cerr << valeur().que_suis_je() << finl;
    }

}

/*! @brief Appel a l'objet sous-jacent.
 *
 * Ajoute la contribution de l'operateur au tableau
 *     passe en parametre
 *
 * @param (DoubleTab& inconnue) tableau contenant les donnees sur lesquelles on applique l'operateur.
 * @param (DoubleTab& resu) tableau auquel on ajoute la contribution de l'operateur
 * @return (DoubleTab&) le tableau contenant le resultat
 */
DoubleTab& Op_Diff_K_Omega::ajouter(const DoubleTab& inconnue, DoubleTab& resu) const
{
  statistiques().begin_count(diffusion_counter_);
  DoubleTab& tmp = valeur().ajouter(inconnue, resu);
  statistiques().end_count(diffusion_counter_);
  return tmp;
}


/*! @brief Appel a l'objet sous-jacent.
 *
 * Initialise le tableau passe en parametre avec la contribution
 *     de l'operateur.
 *
 * @param (DoubleTab& inconnue) tableau contenant les donnees sur lesquelles on applique l'operateur.
 * @param (DoubleTab& resu) tableau dans lequel stocke la contribution de l'operateur
 * @return (DoubleTab&) le tableau contenant le resultat
 */
DoubleTab& Op_Diff_K_Omega::calculer(const DoubleTab& inconnue, DoubleTab& resu) const
{
  statistiques().begin_count(diffusion_counter_);
  DoubleTab& tmp = valeur().calculer(inconnue, resu);
  statistiques().end_count(diffusion_counter_);
  return tmp;
}
