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
// File:        Temperature_imposee_paroi_rayo_semi_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VEF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#include <Temperature_imposee_paroi_rayo_semi_transp.h>
#include <Champ_front_contact_rayo_semi_transp_VEF.h>
#include <Equation_base.h>

Implemente_instanciable(Temperature_imposee_paroi_rayo_semi_transp,"Paroi_temperature_imposee_rayo_semi_transp",Temperature_imposee_paroi);


/*! @brief
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Temperature_imposee_paroi_rayo_semi_transp::printOn(Sortie& os) const
{
  return os;
}



/*! @brief Simple appel a: Temperature_imposee_paroi::readOn(Entree&) Lit la CL a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Temperature_imposee_paroi_rayo_semi_transp::readOn(Entree& is)
{
  return Temperature_imposee_paroi::readOn(is);
}


void Temperature_imposee_paroi_rayo_semi_transp::completer()
{
  Temperature_imposee_paroi::completer();
}



/*! @brief Renvoie un booleen indiquant la compatibilite des conditions aux limites avec l'equation specifiee en parametre.
 *
 *     Des CL de type Temperature_imposee_paroi sont compatibles
 *     avec une equation dont le domaine est la Thermique
 *     ou bien indetermine.
 *
 * @param (Equation_base& eqn) l'equation avec laquelle il faut verifier la compatibilite
 * @return (int) valeur booleenne, 1 si les CL sont compatibles avec l'equation 0 sinon
 */
int Temperature_imposee_paroi_rayo_semi_transp::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Thermique="Thermique";
  Motcle indetermine="indetermine";
  if ( (dom_app==Thermique) || (dom_app==indetermine) )
    return 1;
  else
    {
      err_pas_compatible(eqn);
      return 0;
    }
}


const Cond_lim_base&  Temperature_imposee_paroi_rayo_semi_transp::la_cl() const
{
  return (*this);
}


Champ_front_base&  Temperature_imposee_paroi_rayo_semi_transp::temperature_bord()
{
  return le_champ_front;
}

void  Temperature_imposee_paroi_rayo_semi_transp::calculer_temperature_bord(double temps)
{
  if (sub_type(Champ_front_contact_rayo_semi_transp_VEF,le_champ_front.valeur()))
    {
      Champ_front_contact_rayo_semi_transp_VEF& Ch_contact
        = ref_cast(Champ_front_contact_rayo_semi_transp_VEF,le_champ_front.valeur());
      Ch_contact.calculer_temperature_bord(temps);
    }
  else
    {
      // La temperature de paroi etant directement donnee par le champ_front
      // associe a la condition a la limite, il n'y a rien a calculer ici
      ;
    }
}

void Temperature_imposee_paroi_rayo_semi_transp::associer_modele(const Modele_rayo_semi_transp& un_modele)
{
  Cond_Lim_rayo_semi_transp::associer_modele(un_modele);
  if (sub_type(Champ_front_contact_rayo_semi_transp_VEF,le_champ_front.valeur()))
    {
      Champ_front_contact_rayo_semi_transp_VEF& Ch_contact
        = ref_cast(Champ_front_contact_rayo_semi_transp_VEF,le_champ_front.valeur());
      Ch_contact.associer_modele_rayo(un_modele);
    }
}
