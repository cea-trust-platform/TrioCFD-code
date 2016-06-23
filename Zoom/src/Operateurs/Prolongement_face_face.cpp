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
// File:        Prolongement_face_face.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_face_face.h>
#include <Zone_VDF.h>
#include <Zone_VEF.h>
Implemente_instanciable(Prolongement_face_face,"Prolongement_face_face",Prolongement_base);

//// printOn
//

Sortie& Prolongement_face_face::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_face_face::readOn(Entree& s )
{
  return s ;
}




// Description:
//    Prolongement de l'inconnue grossiere sur la frontiere fine
//    pour inconnue VITESSE
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
void Prolongement_face_face::prolonger(Zone_VF& zone_VFG, Zone_VF& zone_VFF,
                                       const Frontiere& frontF,
                                       IntVect& connect,
                                       const DoubleTab& valG, DoubleTab& tab,
                                       int nb_compo)
{

  int nbFacesF;
  int num_faceG;
  tab = 0.;
  int comp;

  const int prem_face_bord_fin  =  frontF.num_premiere_face();

  int nb_faces_front_fine = frontF.nb_faces();

  for (nbFacesF=0; nbFacesF<nb_faces_front_fine; nbFacesF++)
    {
      const int num_face=prem_face_bord_fin+nbFacesF;
      //on recupere le numero de la face G correspondante
      num_faceG = connect(num_face);
      if (nb_compo == 1)
        {
          tab(nbFacesF, 0) = valG(num_faceG);
        }
      else if(sub_type(Zone_VEF, zone_VFG))
        {
          for (comp=0; comp<nb_compo; comp++)
            {
              tab(nbFacesF, comp) = valG(num_faceG, comp);
            }
        }
      else if(sub_type(Zone_VDF, zone_VFG) && sub_type(Zone_VDF, zone_VFF))
        {
          Zone_VDF& zone_vdfG = ref_cast(Zone_VDF, zone_VFG);

          const IntVect& oriG = zone_vdfG.orientation();
          const IntTab& face_voisG = zone_vdfG.face_voisins();
          const Zone& zoneg = zone_vdfG.zone();
          const int nb_faces_elemG = zoneg.type_elem().nb_faces();
          const IntTab& elem_facesG =  zone_vdfG.elem_faces();
          int num_face_gros;
          int num_elemG;
          if(oriG(num_faceG) == 0)
            {
              //composante en X
              tab(nbFacesF, 0) = valG(num_faceG);
              //somme des composantes normales
              //MODIF POUR IMPOSER LES COMPOSANTES en Y
              num_elemG = face_voisG(num_faceG,0);
              int nb_faces_gros = 0;
              if(num_elemG != -1)
                {
                  //pour chaque faceG
                  for(int nbFacesG=0; nbFacesG<nb_faces_elemG; nbFacesG++)
                    {
                      num_face_gros = elem_facesG(num_elemG, nbFacesG);
                      if(oriG(num_face_gros) == 1)
                        {
                          tab(nbFacesF, 1) = tab(nbFacesF, 1)+ valG(num_face_gros);
                          nb_faces_gros = nb_faces_gros+1;
                        }
                    }
                }
              num_elemG = face_voisG(num_faceG,1);
              if(num_elemG != -1)
                {
                  //pour chaque faceG
                  for(int nbFacesG=0; nbFacesG<nb_faces_elemG; nbFacesG++)
                    {
                      num_face_gros = elem_facesG(num_elemG, nbFacesG);
                      if(oriG(num_face_gros) == 1)
                        {
                          tab(nbFacesF, 1) = tab(nbFacesF, 1) + valG(num_face_gros);
                          nb_faces_gros = nb_faces_gros+1;
                        }
                    }
                }
              tab(nbFacesF, 1) = tab(nbFacesF, 1) / nb_faces_gros;

              ////////////////////////////////// fin des modif pour les compo tgtielles
            }
          else
            {
              tab(nbFacesF, 1) = valG(num_faceG);
              //MODIF POUR IMPOSER LES COMPOSANTES en X
              int num_elemGb = face_voisG(num_faceG,0);
              int nb_faces_gros = 0;
              if(num_elemGb != -1)
                {
                  //pour chaque faceG
                  for(int nbFacesG=0; nbFacesG<nb_faces_elemG; nbFacesG++)
                    {
                      num_face_gros = elem_facesG(num_elemGb, nbFacesG);
                      if(oriG(num_face_gros) == 0)
                        {
                          tab(nbFacesF, 0) = tab(nbFacesF, 0) + valG(num_face_gros);
                          nb_faces_gros = nb_faces_gros+1;
                        }
                    }
                }
              num_elemGb = face_voisG(num_faceG,1);
              if(num_elemGb != -1)
                {
                  //pour chaque faceG
                  for(int nbFacesG=0; nbFacesG<nb_faces_elemG; nbFacesG++)
                    {
                      num_face_gros = elem_facesG(num_elemGb, nbFacesG);
                      if(oriG(num_face_gros) == 0)
                        {
                          tab(nbFacesF, 0) = tab(nbFacesF, 0) + valG(num_face_gros);
                          nb_faces_gros = nb_faces_gros+1;
                        }
                    }
                }
              tab(nbFacesF, 0) = tab(nbFacesF, 0) / nb_faces_gros;
              //////////////// fin des modif pour les compo normales
            }
        }
    }

  //Cerr<<"fin de Prolongement_face_face::prolonger"<<finl;
}




//NE FAIT RIEN
void Prolongement_face_face::calculer(Zone_VF& zonef,
                                      Zone_VF& zoneg,
                                      IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}
