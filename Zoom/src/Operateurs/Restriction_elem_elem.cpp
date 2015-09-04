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
// File:        Restriction_elem_elem.cpp
// Directory:   $TRUST_ROOT/src/Zoom/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Restriction_elem_elem.h>

Implemente_instanciable(Restriction_elem_elem,"Restriction_elem_elem",Restriction_base);

//// printOn
//

Sortie& Restriction_elem_elem::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Restriction_elem_elem::readOn(Entree& s )
{
  return s ;
}




// Description:
//    calcul du nombre d'elements fins contenus dans chaque element grossier
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Restriction_elem_elem::calculer(const Zone_VF& zone_VFG,
                                     const Zone_VF& zone_VFF,
                                     const IntVect& connect)
{
  //Cerr<<"debut de  Restriction_elem_elem::calculer_nb_elemF_dans_chaque_elemG"<<finl;
  int nb_elemG = zone_VFG.nb_elem();
  int nbelemsF;
  int num_elemG;
  int nombre_total_elemF = connect.size_array();
  nb_elemF_.resize(nb_elemG);
  nb_elemF_ = 0;
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      //on recupere le numero de l'element G correspondant
      num_elemG = connect(nbelemsF);
      //on ajoute 1 dans le tableau nb_elemF dans la case correspondant a num_elemG
      nb_elemF_(num_elemG) = nb_elemF_(num_elemG) + 1;
    }

  //Cerr<<"fin de  Restriction_elem_elem::calculer_nb_elemF_dans_chaque_elemG"<<finl;
}




// Description:
//    Restriction de l'inconnue fine sur le domaine grossier
//    pour inconnue TEMPERATURE ET PRESSION
//    en VDF-VDF pour T du domG, en VDF-VDF, VDF-VEF, VEF-VEF
//    et VEF-VDF pour p du domG
//    remarque : pas de couplage!!
//    car valeur
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Restriction_elem_elem::restreindre(const Zone_VF& zone_VFG,
                                        const Zone_VF& zone_VFF,
                                        const IntVect& connect,
                                        DoubleTab& valG,
                                        const DoubleTab& valF,int nb_comp)
{
  // Cerr<<"debut de  Restriction_elem_elem::restreindre"<<finl;
  int nbG;
  int nbF;
  int num_elemG;
  int comp;
  int nb_elem_fin;
  int nb_elemG = zone_VFG.nb_elem();
  int nb_elemFb = zone_VFF.nb_elem();

  //Mettre a 0 les valeurs grossieres a corriger
  for (nbG=0; nbG<nb_elemG; nbG++)
    {
      if (nb_elemF_(nbG) > 1E-9)
        {
          if (nb_comp == 1)
            valG(nbG) = 0.;
          else
            for(comp=0; comp<nb_comp; comp++)
              {
                valG(nbG, comp) = 0.;
              }
        }
    }

  //on ajoute la valeur calculee au centre de gravite
  //de chaque element fin contenu dans l'element grossier a corriger
  for (nbF=0; nbF<nb_elemFb; nbF++)
    {
      num_elemG = connect(nbF);
      if (nb_comp == 1)
        {
          valG(num_elemG) = valG(num_elemG) + valF(nbF);
        }
      else
        for(comp=0; comp<nb_comp; comp++)
          {
            valG(num_elemG, comp) = valG(num_elemG, comp) + valF(nbF, comp);
          }
    }

  //Pour corriger l'inconnue grossiere,
  //on fait la moyenne
  for (nbG=0; nbG<nb_elemG; nbG++)
    {

      nb_elem_fin = nb_elemF_(nbG);

      if (nb_elemF_(nbG) > 1E-9)
        {
          if (nb_comp == 1)
            {
              valG(nbG) /= nb_elem_fin;
            }
          else
            for(comp=0; comp<nb_comp; comp++)
              {
                valG(nbG, comp) = valG(nbG, comp) / nb_elem_fin;
              }
        }
    }
  //Cerr<<"fin de Restriction_elem_elem::restreindre"<<finl;
}

