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
// File:        Prolongement_elem_face.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_elem_face.h>
#include <Domaine.h>
#include <Domaine_VF.h>

Implemente_instanciable(Prolongement_elem_face,"Prolongement_elem_face",Prolongement_base);

//// printOn
//

Sortie& Prolongement_elem_face::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_elem_face::readOn(Entree& s )
{
  return s ;
}




/*! @brief Prolongement de l'inconnue grossiere sur la frontiere fine pour inconnue TEMPERATURE ET PRESSION
 *
 *     en VDF pour T du domG, en VDF et VEF pour p du domG
 *     en 2D et en 3D
 *     (avec une discretisation quelconque pour le domaine fin
 *     qui ne contient pas de domaine plus petit que lui-meme)
 *     remarque : pas de couplage!!
 *     car valeur
 *
 */
void Prolongement_elem_face::prolonger(Domaine_VF& domaine_VFG, Domaine_VF& domaine_VFF,
                                       const Frontiere& frontF,
                                       IntVect& connect,
                                       const DoubleTab& valG, DoubleTab& tab,
                                       int nb_comp)
{
  //Cerr<<"debut de Prolongement_elem_face::prolonger"<<finl;

  int nbFacesF;
  int num_faceG;
  int num_elemG;
  int nbElemG;
  double somme;
  double somme_dist;
  double d;

  const int prem_face_bord_fin  =  frontF.num_premiere_face();
  int nb_faces_front_fine = frontF.nb_faces();

  //elemG voisins de chaque face grossiere
  IntTab& voisins = domaine_VFG.face_voisins();

  //Pour chaque face fine
  for (nbFacesF=0; nbFacesF<nb_faces_front_fine; nbFacesF++)
    {
      int num_face=prem_face_bord_fin+nbFacesF;
      //on recupere le numero de la face G correspondante
      num_faceG = connect(num_face);
      //Pour chaque elemG voisin de cette faceG
      somme = 0;
      somme_dist = 0;
      for (nbElemG=0; nbElemG<2; nbElemG++)
        {
          num_elemG = voisins(num_faceG, nbElemG);
          if(num_elemG != -1)
            {
              d = distances_(num_face, nbElemG);
              somme = somme + valG(num_elemG) * d;
              somme_dist = somme_dist + d;
            }
        }
      tab(nbFacesF, 0) = somme / somme_dist;
    }
  //Cerr<<"fin de Prolongement_elem_face::prolonger"<<finl;
}




//calcul des distances du centre de gravite de chaque face fine
//au centre de gravite des elements grossiers correspondants
void Prolongement_elem_face::calculer(Domaine_VF& domainef,
                                      Domaine_VF& domaineg,
                                      IntVect& connect_ff)
{
  //Cerr<<"debut de Prolongement_elem_face::calculer_distances"<<finl;
  int nbFacesF;
  int nbElemG;
  int num_faceG;
  int num_elemG;
  int dim;
  double somme_carres;
  const int prem_face_bord_fin  =  domainef.premiere_face_bord();
  const int der_face_bord_fin  =  prem_face_bord_fin+domainef.nb_faces_bord();
  //nombre de faces fines
  int nb_lignes = connect_ff.size_array();
  distances_.resize(nb_lignes, 2);
  distances_ = -1;
  //elemG voisins de chaque face grossiere
  IntTab& voisins = domaineg.face_voisins();
  //centre de gravite des elemG
  const DoubleTab& cg_elem_grossier = domaineg.xp();
  //centre de gravite des faces fines
  const DoubleTab& cg_face_fine = domainef.xv();

  //Pour chaque face fine
  for (nbFacesF=prem_face_bord_fin; nbFacesF<der_face_bord_fin; nbFacesF++)
    {
      //on recupere le numero de la faceG correspondant a le face fine
      num_faceG = connect_ff(nbFacesF);
      if(num_faceG != -1)
        {
          //Pour chaque elemG voisin de cette faceG (il y en a au maximun 2!!)
          for (nbElemG=0; nbElemG<2; nbElemG++)
            {
              num_elemG = voisins(num_faceG, nbElemG);
              //calcul de la distance
              //si il n'y a pas d'elemG correspondant
              //(si la face fine est sur le bord grossier)
              if(num_elemG == -1)
                somme_carres = 0.;
              else
                {
                  somme_carres = 0.;
                  for(dim=0; dim<dimension; dim++)
                    somme_carres = somme_carres + (cg_face_fine(nbFacesF, dim)-cg_elem_grossier(num_elemG, dim))*(cg_face_fine(nbFacesF, dim)-cg_elem_grossier(num_elemG, dim));
                }
              distances_(nbFacesF, nbElemG) = sqrt(somme_carres);
            }
        }
    }
  //Cerr<<"fin de Prolongement_elem_face::calculer_distances"<<finl;
}
