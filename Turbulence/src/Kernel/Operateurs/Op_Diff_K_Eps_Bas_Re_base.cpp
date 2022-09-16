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
// File:        Op_Diff_K_Eps_Bas_Re_base.cpp
// Directory:   $TURBULENCE_ROOT/src/Kernel/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Eps_Bas_Re_base.h>
#include <Champ_Don.h>
#include <Champ_Uniforme.h>
#include <Discretisation_base.h>

Implemente_base(Op_Diff_K_Eps_Bas_Re_base,"Op_Diff_K_Eps_Bas_Re_base",Operateur_base);
Implemente_instanciable(Op_Diff_K_Eps_Bas_Re_negligeable,"Op_Diff_K_Eps_Bas_Re_negligeable",Op_Diff_K_Eps_Bas_Re_base);
Implemente_deriv(Op_Diff_K_Eps_Bas_Re_base);
Implemente_instanciable(Op_Diff_K_Eps_Bas_Re,"Op_Diff_K_Eps_Bas_Re",DERIV(Op_Diff_K_Eps_Bas_Re_base));


/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie& s) le flot de sortie modifie
 */
Sortie& Op_Diff_K_Eps_Bas_Re_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie& s) le flot de sortie modifie
 */
Sortie& Op_Diff_K_Eps_Bas_Re_negligeable::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}



/*! @brief Associe la diffusivite a l'operateur.
 *
 * @param (Champ_Don& ch) le champ representant la diffusivite
 */
void Op_Diff_K_Eps_Bas_Re_negligeable::associer_diffusivite(const Champ_base& ch)
{
  la_diffusivite=ch;
}


/*! @brief Renvoie le champ representant la diffusivite.
 *
 * @return (Champ_Don_base&) le champ representant la diffusivite
 */
const Champ_base& Op_Diff_K_Eps_Bas_Re_negligeable::diffusivite() const
{
  return la_diffusivite;
}





/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie& s) le flot de sortie modifie
 */
Sortie& Op_Diff_K_Eps_Bas_Re::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


/*! @brief NE FAIT RIEN A surcharger dans les classes derivees
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree
 */
Entree& Op_Diff_K_Eps_Bas_Re_base::readOn(Entree& s )
{
  return s ;
}

/*! @brief NE FAIT RIEN A surcharger dans les classes derivees
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree
 */
Entree& Op_Diff_K_Eps_Bas_Re_negligeable::readOn(Entree& s )
{
  return s ;
}

/*! @brief Simple appel a Operateur::lire(Entree&) Lit l'operateur a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Op_Diff_K_Eps_Bas_Re::readOn(Entree& s )
{
  Operateur::lire(s);
  return s;
}



/*! @brief Type l'operateur s'associe a son equation
 *
 *     Associe la diffusivite turbulente a l'operateur base.
 *
 */
void Op_Diff_K_Eps_Bas_Re::typer()
{
  Cerr << "dans Op_Diff_K_Eps_Bas_Re::typer()  typ = " << typ << finl;
  if (typ=="negligeable")
    {
      DERIV(Op_Diff_K_Eps_Bas_Re_base)::typer("Op_Diff_K_Eps_Bas_Re_negligeable");
      valeur().associer_diffusivite_turbulente();
    }
  else
    {
      assert(la_diffusivite.non_nul());
      Nom nom_type;
      Cerr << "equation().valeur() " << equation() << finl;
      //if ( (!sub_type(Transport_K_Eps,equation().que_suis_je())) && (!sub_type(Transport_K_Eps_Bas_Re,equation().que_suis_je())) )
      nom_type="Op_Diff_K_Eps_Bas_Re_";
      Nom discr = equation().discretisation().que_suis_je();
      if (discr=="VEFPreP1B") discr = "VEF";
      Nom type_diff;
      nom_type +=discr;
      if (discr!="VDF_Hyper")
        {
          if(sub_type(Champ_Uniforme,diffusivite()))
            type_diff="_const_";
          else
            type_diff="_var_";
          nom_type+= type_diff;
          Nom type_inco=equation().inconnue()->que_suis_je();
          nom_type+=(type_inco.suffix("Champ_"));
          if (axi)
            nom_type += "_Axi";
        }
      Cerr << " type = " << nom_type << finl;
      DERIV(Op_Diff_K_Eps_Bas_Re_base)::typer(nom_type);
      valeur().associer_eqn(equation());
      valeur().associer_diffusivite_turbulente();
      valeur().associer_diffusivite(diffusivite());
    }
  Cerr << " fin du typage on a " << valeur().que_suis_je() << finl;
}

