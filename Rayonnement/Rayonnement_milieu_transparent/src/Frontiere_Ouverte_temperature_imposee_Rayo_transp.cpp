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
// File:        Frontiere_Ouverte_temperature_imposee_Rayo_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#include <Frontiere_Ouverte_temperature_imposee_Rayo_transp.h>
#include <Equation_base.h>
#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Zone_VF.h>

Implemente_instanciable(Frontiere_Ouverte_temperature_imposee_Rayo_transp,"Frontiere_Ouverte_temperature_imposee_Rayo_transp",Entree_fluide_temperature_imposee);


// printOn et readOn

Sortie& Frontiere_Ouverte_temperature_imposee_Rayo_transp::printOn(Sortie& is ) const
{
  return is;
}

Entree& Frontiere_Ouverte_temperature_imposee_Rayo_transp::readOn(Entree& s )
{
  /*
    le_champ_front.typer("Champ_front_uniforme");
    Motcle motlu;
    Motcles les_motcles(1);
    {
    les_motcles[0] = "T_ext";
    }
    s >> motlu;
    int rang = les_motcles.search(motlu);
    switch(rang)
    {
    case 0:
    {
    s >> le_champ_ext;
    break;
    }
    default:
    {
    Cerr << "Erreur a la lecture de la condition aux limites de type: " << finl;
    Cerr << que_suis_je() << finl;
    Cerr << "On attendait " << les_motcles << " a la place de " <<  motlu << finl;
    exit();
    }
    }
    return s;
  */
  return Entree_fluide_temperature_imposee::readOn(s);
}


void Frontiere_Ouverte_temperature_imposee_Rayo_transp::completer()
{
  //  Paroi_Rayo_transp::completer();
  Entree_fluide_temperature_imposee::completer();
  preparer_surface(frontiere_dis(),zone_Cl_dis());
}

void Frontiere_Ouverte_temperature_imposee_Rayo_transp::mettre_a_jour(double temps)
{
  Entree_fluide_temperature_imposee::mettre_a_jour(temps);
  calculer_Teta_i();
}


void Frontiere_Ouverte_temperature_imposee_Rayo_transp::calculer_Teta_i()
{
  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int nb_faces_bord = front_vf.nb_faces();
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    Teta_i[numfa]=val_imp(numfa);
}

