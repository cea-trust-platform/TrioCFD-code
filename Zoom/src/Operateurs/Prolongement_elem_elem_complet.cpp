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
// File:        Prolongement_elem_elem_complet.cpp
// Directory:   $TRUST_ROOT/src/Zoom/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_elem_elem_complet.h>
#include <Domaine.h>
#include <Zone_VF.h>

Implemente_instanciable(Prolongement_elem_elem_complet,"Prolongement_elem_elem_complet",Prolongement_base);

//// printOn
//

Sortie& Prolongement_elem_elem_complet::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_elem_elem_complet::readOn(Entree& s )
{
  return s ;
}




// Description:
//    Prolongement du centre de gravite des faces grossieres
//    au centre de gravite de TOUTES les faces fines
//    PROLONGEMENT POUR PRESSION ET TEMPERATURE
//    EN 2D ET 3D
//    EN VDF ET VEF pour la pression
//    EN VDF pour la temperature
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

void Prolongement_elem_elem_complet::prolonger(Zone_VF& zone_VFG,
                                               Zone_VF& zone_VFF,
                                               const Frontiere& frontF,
                                               IntVect& connect,
                                               const DoubleTab& valG, DoubleTab& tab,
                                               int nb_comp)
{

  //Cerr<<"debut de Prolongement_elem_elem_complet::prolonger"<<finl;
  int nb_elemF = connect.size_array();
  int nbelemsF;
  int num_elemG;
  int nbfaceG;
  Zone& zoneg = zone_VFG.zone();
  int nb_faceG = zoneg.nb_faces_elem();

  IntVect numerosG1;
  DoubleVect numerosG1_dist;
  IntVect numerosG2;
  DoubleVect numerosG2_dist;
  numerosG1.resize(nb_faceG+1);
  numerosG2.resize(nb_faceG+1);
  numerosG1_dist.resize(nb_faceG+1);
  numerosG2_dist.resize(nb_faceG+1);
  int numG, num_faceG;
  IntTab& num_face_elemG = zone_VFG.elem_faces();

  //coord des centres de gravite des elems
  tab = 0.;

  DoubleVect coord_somG;
  coord_somG.resize(dimension);
  ArrOfInt elem;
  int compo;

  int nb_elemG_pour_prolongement;

  //elemG voisins de chaque face grossiere
  IntTab& voisins = zone_VFG.face_voisins();

  //Pour chaque elemF
  for (nbelemsF=0; nbelemsF<nb_elemF; nbelemsF++)
    {
      nb_elemG_pour_prolongement = 0;
      // on recupere le numero de l'elemG correspondant
      num_elemG = connect(nbelemsF);
      //on ajoute la valeur valG(num_elemG)
      if (nb_comp == 1)
        {
          tab(nbelemsF) = tab(nbelemsF) + valG(num_elemG);
          nb_elemG_pour_prolongement = nb_elemG_pour_prolongement+1;
        }
      else
        {
          for(compo = 0; compo<nb_comp; compo++)
            tab(nbelemsF, compo) = tab(nbelemsF, compo) + valG(num_elemG, compo);
          nb_elemG_pour_prolongement = nb_elemG_pour_prolongement+1;
        }


      //ELEM VOISIN G CORRESPONDANT
      //Pour chaque faceG
      for (nbfaceG=0; nbfaceG<nb_faceG; nbfaceG++)
        {
          num_faceG = num_face_elemG(num_elemG, nbfaceG);

          //on trouve l'elemG voisin de num_elemG
          numG = voisins(num_faceG, 0);

          if (numG == num_elemG )
            {
              numG = voisins(num_faceG, 1);   // element exterieur
              if (numG != -1)
                {
                  if (nb_comp == 1)
                    {
                      tab(nbelemsF) = tab(nbelemsF) + valG(num_elemG);
                      nb_elemG_pour_prolongement = nb_elemG_pour_prolongement+1;
                    }
                  else
                    {
                      for(compo = 0; compo<nb_comp; compo++)
                        tab(nbelemsF, compo) = tab(nbelemsF, compo) + valG(num_elemG, compo);
                      nb_elemG_pour_prolongement = nb_elemG_pour_prolongement+1;
                    }
                }
            }
        }


      if (nb_comp == 1)
        {
          tab(nbelemsF) = tab(nbelemsF)/nb_elemG_pour_prolongement;
        }
      else
        {
          for(compo = 0; compo<nb_comp; compo++)
            tab(nbelemsF, compo) = tab(nbelemsF, compo)/nb_elemG_pour_prolongement;
        }
    }
  //Cerr<<"fin de Prolongement_elem_elem_complet::prolonger"<<finl;
}




//NE FAIT RIEN
void Prolongement_elem_elem_complet::calculer(Zone_VF& zonef,
                                              Zone_VF& zoneg,
                                              IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}
