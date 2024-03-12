/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Entree_fluide_temperature_imposee_H.cpp
// Directory:   $TRIO_U_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Entree_fluide_temperature_imposee_H.h>
#include <Fluide_Quasi_Compressible.h>
#include <Equation_base.h>
#include <Motcle.h>

Implemente_instanciable(Entree_fluide_temperature_imposee_H,"Entree_temperature_imposee_H",Entree_fluide_temperature_imposee);

/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Entree_fluide_temperature_imposee_H::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

/*! @brief Simple appel a: Cond_lim_base::readOn(Entree& )
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree& s) le flot d'entree modifie
 */
Entree& Entree_fluide_temperature_imposee_H::readOn(Entree& s)
{
  return Entree_fluide_temperature_imposee::readOn(s);
}

/*! @brief Complete les conditions aux limites.
 *
 */
void Entree_fluide_temperature_imposee_H::completer()
{
  le_fluide = ref_cast(Fluide_Quasi_Compressible,mon_dom_cl_dis->equation().milieu());
}

/*! @brief Renvoie un booleen indiquant la compatibilite des conditions aux limites avec l'equation specifiee en parametre.
 *
 *     Des CL de type Entree_fluide_temperature_imposee sont compatibles
 *     avec une equation dont le domaine est la Thermique
 *     ou bien indetermine.
 *
 * @param (Equation_base& eqn) l'equation avec laquelle il faut verifier la compatibilite
 * @return (int) valeur booleenne, 1 si les CL sont compatibles avec l'equation 0 sinon
 */
int Entree_fluide_temperature_imposee_H::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Thermique="Thermique";
  if ( (dom_app==Thermique))
    return 1;
  else
    {
      err_pas_compatible(eqn);
      return 0;
    }
}
void Entree_fluide_temperature_imposee_H::mettre_a_jour(double temps)
{
  Entree_fluide_temperature_imposee::mettre_a_jour(temps);
  double Pth=le_fluide.valeur().pression_th();
  coef_P_sur_R_=le_fluide.valeur().loi_etat().valeur().calculer_masse_volumique(Pth,1.);

}
/*! @brief Renvoie la valeur imposee sur la i-eme composante du champ a la frontiere.
 *
 * @param (int i) indice suivant la premiere dimension du champ
 * @return (double) la valeur imposee sur la composante du champ specifiee
 * @throws deuxieme dimension du champ de frontiere superieur a 1
 */
double Entree_fluide_temperature_imposee_H::val_imp(int i) const
{
//Cerr<<" VAL IMP "<<modifier_val_imp<< finl;
  double T;
  if (le_champ_front.valeurs().size()==1)
    {
      T= (le_champ_front(0,0));

    }
  else if (le_champ_front.valeurs().dimension(1)==1)
    {
      T= le_champ_front(i,0);

    }
  else
    {
      Cerr << "Temperature_imposee_paroi_H::val_imp erreur" << finl;
      T=-1.;
      abort();
    }
  if (modifier_val_imp==1)
    {
      double rho=coef_P_sur_R_/T;
      //Cerr<<"iiiiiii "<<rho<<finl;
      return rho;
    }
  else
    return T;
}


/*! @brief Renvoie la valeur imposee sur la (i,j)-eme composante du champ a la frontiere.
 *
 * @param (int i) indice suivant la premiere dimension du champ
 * @param (int j) indice suivant la deuxieme dimension du champ
 * @return (double) la valeur imposee sur la composante du champ specifiee
 */
double Entree_fluide_temperature_imposee_H::val_imp(int i, int j) const
{
  if (le_champ_front.valeurs().dimension(0)==1)
    {
      if (modifier_val_imp==1)
        return le_fluide->calculer_H(le_champ_front(0,j));
      else
        return le_champ_front(0,j);

    }
  else
    {
      if (modifier_val_imp==1)
        return le_fluide->calculer_H(le_champ_front(i,j));
      else
        return le_champ_front(i,j);
    }
}
