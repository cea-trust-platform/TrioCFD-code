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
// File:        Temperature_imposee_paroi_rayo_transp.cpp
// Directory:   $TRUST_ROOT/src/Rayonnement/VEF
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Temperature_imposee_paroi_rayo_transp.h>
#include <Equation_base.h>
#include <Zone_VEF.h>
#include <Champ_front_contact_rayo_transp_VEF.h>

Implemente_instanciable(Temperature_imposee_paroi_rayo_transp,"Paroi_temperature_imposee_rayo_transp",Temperature_imposee_paroi);


// printOn et readOn

Sortie& Temperature_imposee_paroi_rayo_transp::printOn(Sortie& is ) const
{
  return is;
}

Entree& Temperature_imposee_paroi_rayo_transp::readOn(Entree& s )
{
  return Temperature_imposee_paroi::readOn(s);
}


void Temperature_imposee_paroi_rayo_transp::completer()
{
  //  Cerr<<"Temperature_imposee_paroi_rayo_transp::completer debut"<<finl;
  Temperature_imposee_paroi::completer();
  preparer_surface(frontiere_dis(),zone_Cl_dis());
}

void Temperature_imposee_paroi_rayo_transp::mettre_a_jour(double temps)
{
  //  Cerr<<"Temperature_imposee_paroi_rayo_transp::mettre_a_jour"<<finl;
  calculer_Teta_i(temps);
}


void Temperature_imposee_paroi_rayo_transp::calculer_Teta_i(double temps)
{
  //  Cerr<<"Temperature_imposee_paroi_rayo_transp::calculer_Teta_i()"<<finl;
  if (sub_type(Champ_front_contact_rayo_transp_VEF,le_champ_front.valeur()))
    {
      Champ_front_contact_rayo_transp_VEF& Ch_contact
        = ref_cast(Champ_front_contact_rayo_transp_VEF,le_champ_front.valeur());
      Ch_contact.calculer_temperature_bord(temps);
    }
  else
    {
      // La temperature de paroi etant directement donnee par le champ_front
      // associe a la condition a la limite, il n'y a rien a calculer ici
      ;
    }

  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int nb_faces_bord = front_vf.nb_faces();
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    {
      Teta_i[numfa]=val_imp(numfa);
    }
}


void Temperature_imposee_paroi_rayo_transp::associer_modele_rayo(Modele_Rayonnement_base& mod)
{
  le_modele_rayo = ref_cast(Modele_Rayonnement_Milieu_Transparent, mod);

  if (sub_type(Champ_front_contact_rayo_transp_VEF,le_champ_front.valeur()))
    {
      Champ_front_contact_rayo_transp_VEF& Ch_contact
        = ref_cast(Champ_front_contact_rayo_transp_VEF,le_champ_front.valeur());
      Ch_contact.associer_modele_rayo(modele_rayo());
    }
}
