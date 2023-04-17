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
// File:        Restriction_elem_elem_2.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Restriction_elem_elem_2.h>

Implemente_instanciable(Restriction_elem_elem_2,"Restriction_elem_elem_2",Restriction_base);

//// printOn
//

Sortie& Restriction_elem_elem_2::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Restriction_elem_elem_2::readOn(Entree& s )
{
  return s ;
}




/*! @brief calcul du nombre d'elements fins contenus dans chaque element grossier
 *
 */
void Restriction_elem_elem_2::calculer(const Domaine_VF& domaine_VFG,
                                       const Domaine_VF& domaine_VFF,
                                       const IntVect& connect)
{
  //Cerr<<"debut de  Restriction_elem_elem_2::calculer_nb_elemF_pour_chaque_elemG"<<finl;
  int nb_elemG = domaine_VFG.nb_elem();
  int nbelemsF;
  int num_elemG;
  int est_dedans;
  int nombre_total_elemF = connect.size_array();
  nb_elemF_.resize(nb_elemG);
  nb_elemF_ = 0;
  int nbfacesF;
  const Domaine& domainef = domaine_VFF.domaine();
  const Domaine& domaineg = domaine_VFG.domaine();
  int nombre_face_elemF = domainef.nb_faces_elem();
  const Elem_geom& elem_geomG = domaineg.type_elem();
  int num_elem_voisin;
  const IntTab& voisinsF = domaine_VFF.face_voisins();
  int num_face;
  const IntTab& num_face_elemF = domaine_VFF.elem_faces();
  const DoubleTab& cg_elem_fin = domaine_VFF.xp();
  DoubleVect coordF;
  coordF.resize(dimension);
  int dim;




  //NOMBRE D'ELEMENTS FINS INTERIEURS A CHAQUE ELEMENT GROSSIER
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      //on recupere le numero de l'element G correspondant
      num_elemG = connect(nbelemsF);
      //on ajoute 1 dans le tableau nb_elemF dans la case correspondant a num_elemG
      nb_elemF_(num_elemG) = nb_elemF_(num_elemG) + 1;
    }



  //ON COMPTE LES VOISINS EXTERIEURS A L'ELEM GROSSIER
  //POUR CHAQUE ELEM FIN DANS CET ELEM GROSSIER
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      num_elemG = connect(nbelemsF);
      //Pour chaque face de l'elemF
      for (nbfacesF=0; nbfacesF<nombre_face_elemF; nbfacesF++)
        {
          num_face = num_face_elemF(nbelemsF, nbfacesF);
          //voisin de la face
          num_elem_voisin = voisinsF(num_face, 0);
          if (num_elem_voisin == nbelemsF)
            num_elem_voisin = voisinsF(num_face, 1);
          //Ce voisin - s'il existe - est-il exterieur a l'elemG num_elemG?
          if(num_elem_voisin != -1)
            {
              for(dim=0; dim<dimension; dim++)
                coordF(dim) = cg_elem_fin(num_elem_voisin, dim);
              //regardons si son centre de gravite est dans l'elemG ou non
              est_dedans = elem_geomG.contient(coordF, num_elemG);
              if (est_dedans == 0)
                {
                  //s'il n'est pas dans l'elemG
                  //on le compte dans nb_elemF_
                  nb_elemF_(num_elemG) = nb_elemF_(num_elemG) + 1;
                }
            }
        }
    }
  //Cerr<<"fin de  Restriction_elem_elem_2::calculer_nb_elemF_pour_chaque_elemG"<<finl;
}




/*! @brief Restriction de l'inconnue fine sur le domaine grossier en utilisant les elemF interieurs a l'elemG ET
 *
 *     les elemF autour de l'elemG
 *
 */
void Restriction_elem_elem_2::restreindre(const Domaine_VF& domaine_VFG,
                                          const Domaine_VF& domaine_VFF,
                                          const IntVect& connect,
                                          DoubleTab& valG,
                                          const DoubleTab& valF,int nb_comp)
{
  //Cerr<<"debut de  Restriction_elem_elem_2::restreindre"<<finl;
  int nbG;
  int nbF;
  int num_elemG;
  int comp;
  int nb_elem_fin;
  int nb_elemG = domaine_VFG.nb_elem();
  int nb_elemFb = domaine_VFF.nb_elem();
  int nbelemsF;
  int nombre_total_elemF = connect.size_array();
  int nbfacesF;

  const Domaine& domainef = domaine_VFF.domaine();
  const Domaine& domaineg = domaine_VFG.domaine();
  int nombre_face_elemF = domainef.nb_faces_elem();
  const Elem_geom& elem_geomG = domaineg.type_elem();

  int num_face;
  const IntTab& num_face_elemF = domaine_VFF.elem_faces();
  int num_elem_voisin;
  const IntTab& voisinsF = domaine_VFF.face_voisins();
  int dim;
  DoubleVect coordF;
  coordF.resize(dimension);
  const DoubleTab& cg_elem_fin = domaine_VFF.xp();
  int est_dedans;




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



  //on ajoute la valeur calculee au centre de gravite
  //de chaque element fin situe autour de l'element grossier a corriger
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      num_elemG = connect(nbelemsF);
      //Pour chaque face de l'elemF
      for (nbfacesF=0; nbfacesF<nombre_face_elemF; nbfacesF++)
        {
          num_face = num_face_elemF(nbelemsF, nbfacesF);
          //voisin de la face
          num_elem_voisin = voisinsF(num_face, 0);
          if (num_elem_voisin == nbelemsF)
            num_elem_voisin = voisinsF(num_face, 1);
          //Ce voisin - s'il existe - est-il exterieur a l'elemG num_elemG?
          if(num_elem_voisin != -1)
            {
              for(dim=0; dim<dimension; dim++)
                coordF(dim) = cg_elem_fin(num_elem_voisin, dim);
              //regardons si son centre de gravite est dans l'elemG ou non
              est_dedans = elem_geomG.contient(coordF, num_elemG);
              if (est_dedans == 0)
                {
                  //s'il n'est pas dans l'elemG
                  //on le compte dans nb_elemF_
                  if (nb_comp == 1)
                    {
                      valG(num_elemG) = valG(num_elemG) + valF(num_elem_voisin);
                    }
                  else
                    for(comp=0; comp<nb_comp; comp++)
                      {
                        valG(num_elemG, comp) = valG(num_elemG, comp) + valF(num_elem_voisin, comp);
                      }
                }
            }
        }
    }


  //Pour corriger l'inconnue grossiere,
  //on fait la moyenne
  for (nbG=0; nbG<nb_elemG; nbG++)
    {
      nb_elem_fin = nb_elemF_(nbG);
      if (nb_elem_fin > 1E-9)
        {
          if (nb_comp == 1)
            {
              valG(nbG) = valG(nbG) / nb_elem_fin;
            }
          else
            for(comp=0; comp<nb_comp; comp++)
              {
                valG(nbG, comp) = valG(nbG, comp) / nb_elem_fin;
              }
        }
    }
  //Cerr<<"fin de Restriction_elem_elem_2::restreindre"<<finl;
}

