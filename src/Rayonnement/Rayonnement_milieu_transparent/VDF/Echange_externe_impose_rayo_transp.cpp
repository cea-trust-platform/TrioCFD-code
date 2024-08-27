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
// File:        Echange_externe_impose_rayo_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src/VDF
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_externe_impose_rayo_transp.h>
//providoire
////#include <Paroi_2couches_scal_VDF.h>

Implemente_instanciable(Echange_externe_impose_rayo_transp,"Paroi_Echange_externe_impose_rayo_transp",Echange_externe_impose);


// printOn et readOn

Sortie& Echange_externe_impose_rayo_transp::printOn(Sortie& is ) const
{
  return is;
}



Entree& Echange_externe_impose_rayo_transp::readOn(Entree& is )
{

  Motcle motlu;
  Motcles les_motcles(2);
  {
    les_motcles[0] = "h_imp";
    les_motcles[1] = "T_ext";
  }

  int ind = 0;
  while (ind < 2)
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);

      switch(rang)
        {
        case 0:
          {
            is >> h_imp_;
            break;
          }
        case 1:
          {
            is >> T_ext();
            break;
          }
        default:
          {
            Cerr << "Erreur a la lecture de la condition aux limites de type "<<finl;
            Cerr << "Echange_externe_impose_rayo_transp " << finl;
            Cerr << "On attendait " << les_motcles << "a la place de " <<  motlu << finl;
            exit();
          }
        }
      ind++;
    }

  if(local_min_vect(h_imp_->valeurs()) < 1.e9)
    {
      Cerr<<"Erreur sur l'utilisation de la condition a la limite"<<finl;
      Cerr<<"Echange_externe_impose_rayo_transp. Celle ci ne peut"<<finl;
      Cerr<<"etre utilisee pour un probleme de rayonnement que pour "<<finl;
      Cerr<<"imposer une temperature sur une paroi"<<finl;
      exit();
    }

  return is;
}



void Echange_externe_impose_rayo_transp::completer()
{
  Echange_externe_impose::completer();
  preparer_surface(frontiere_dis(),domaine_Cl_dis());
}


void  Echange_externe_impose_rayo_transp::calculer_Teta_i()
{
  const Front_VF& front_vf = ref_cast(Front_VF,frontiere_dis());
  int nb_faces_bord = front_vf.nb_faces();
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    Teta_i[numfa]=T_ext(numfa);
}
void Echange_externe_impose_rayo_transp::mettre_a_jour(double temps)
{
  Echange_externe_impose::mettre_a_jour(temps);
  calculer_Teta_i();
}

