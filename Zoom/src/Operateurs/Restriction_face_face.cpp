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
// File:        Restriction_face_face.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Restriction_face_face.h>

Implemente_instanciable(Restriction_face_face,"Restriction_face_face",Restriction_base);

//// printOn
//

Sortie& Restriction_face_face::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Restriction_face_face::readOn(Entree& s )
{
  return s ;
}




// Description:
//    calcul du nombre de faces fines dont l'elem contient le centre
//    de gravite pour chaque face grossiere
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
void Restriction_face_face::calculer(const Zone_VF& zone_VFG,
                                     const Zone_VF& zone_VFF,
                                     const IntVect& connect)
{
  //Cerr<<"debut de  Restriction_face_face::calculer"<<finl;
  int nb_elemG = zone_VFG.nb_elem();
  int nb_faceG = zone_VFG.nb_faces();
  const DoubleTab& cg_facesG = zone_VFG.xv();
  int max1, max2, max3;
  int nbelemsF;
  IntVect nb_elemF_;
  IntTab tab_facesF;
  int num_elemG;
  int nbfacesF;
  int nbfacesG;
  const Zone& zoneF = zone_VFF.zone();
  int nb_facesF_par_elemF = zoneF.nb_faces_elem();
  int nb_facesG_par_elemG = zone_VFG.zone().nb_faces_elem();
  int nombre_total_elemF = connect.size_array();
  nb_elemF_.resize(nb_elemG);
  nb_faceF_.resize(nb_faceG+1);
  nb_faceF_ = 0;
  const IntTab& num_facesF = zone_VFF.elem_faces();
  const IntTab& num_facesG = zone_VFG.elem_faces();
  const Elem_geom& mon_elemF = zoneF.type_elem();
  int dim;
  DoubleVect coord_cgG;
  coord_cgG.resize(dimension);
  int num_faceG;

  //On va stocker dans nb_faceF_
  //le nb de faces fines dont l'element fin contient
  //le centre de gravite de la face de chaque elemG
  //et dans la derniere case : le nombre maximal de faces
  //fines qui puissent contenir un centre de gravite de face grossiere

  //ON COMPLETE LE VECTEUR nb_elemF_
  nb_elemF_ = 0;
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      //on recupere le numero de l'element G correspondant
      num_elemG = connect(nbelemsF);
      //on ajoute 1 dans le tableau nb_elemF dans la case correspondant a num_elemG
      nb_elemF_(num_elemG) = nb_elemF_(num_elemG) + 1;
    }




  //ON CALCULE UN MAX "LARGE" DU NB DE FACES FINES DANS UN ELEM GROSSIER
  max1 = 0;
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      //on recupere le numero de l'element G correspondant
      num_elemG = connect(nbelemsF);
      int nb = nb_elemF_(num_elemG) * nb_facesF_par_elemF;
      if (max1 < nb)
        max1 = nb;
    }
  //le centre de gravite grossier appartient au plus a 2 elements
  //ainsi les faces fines dans ces 2 elements peuvent intervenir
  max1 = max1 * 2;

  //ON COMPLETE L'ATTRIBUT nb_faceF_
  tab_facesF.resize(nb_faceG, max1);
  tab_facesF = -1;
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      //on recupere le numero de l'element G correspondant
      num_elemG = connect(nbelemsF);
      //Pour chaque face grossiere
      for(nbfacesG=0; nbfacesG<nb_facesG_par_elemG; nbfacesG++)
        {
          //numero global de la face grossiere
          num_faceG = num_facesG(num_elemG, nbfacesG);
          //On regarde si son centre de gravite est dans l'elemF
          for(dim=0; dim<dimension; dim++)
            coord_cgG(dim) = cg_facesG(num_faceG, dim);
          int est_dedans = mon_elemF.contient(coord_cgG, nbelemsF);
          //s'il est dedans
          if ( est_dedans == 1 )
            {
              //Pour chaque face fine
              for(nbfacesF=0; nbfacesF<nb_facesF_par_elemF; nbfacesF++)
                {
                  //numero de la face fine
                  int numF = num_facesF(nbelemsF, nbfacesF);
                  //A-t-elle deja ete comptee?
                  int i = 0;
                  int trouve = 0;
                  while( (tab_facesF(num_faceG, i) != -1) && (i<max1) && (trouve == 0) )
                    {
                      if (tab_facesF(num_faceG, i) == numF)
                        trouve = 1;
                      else
                        i = i+1;
                    }
                  if (trouve == 0)
                    {
                      tab_facesF(num_faceG, i) = numF;
                      nb_faceF_(num_faceG) = nb_faceF_(num_faceG)+1;
                    }
                }
            }
        }
    }


  //ON TROUVE LE MAX "EXACT" ET ON LE STOCKE DANS LA DERNIERE CASE DU TABLEAU nb_faceF_
  max2 = 0;
  //Pour chaque face grossiere
  for(nbfacesG=0; nbfacesG<nb_faceG; nbfacesG++)
    {
      max3 = nb_faceF_(nbfacesG);
      if ( max2 < max3 )
        max2 = max3;
    }
  nb_faceF_(nbfacesG) = max2;

  //Cerr<<"fin de  Restriction_face_face::calculer"<<finl;
}




// Description:
//    Restriction de l'inconnue fine sur le domaine grossier
//    pour inconnue TEMPERATURE
//    en VDF-VEF pour T du domG,
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
void Restriction_face_face::restreindre(const Zone_VF& zone_VFG,
                                        const Zone_VF& zone_VFF,
                                        const IntVect& connect,
                                        DoubleTab& valG,
                                        const DoubleTab& valF,int nb_comp)
{
  //Cerr<<"debut de  Restriction_face_face::restreindre"<<finl;
  const DoubleTab& cg_facesG = zone_VFG.xv();
  int nb_faceG = zone_VFG.nb_faces();
  const Zone& zoneF = zone_VFF.zone();
  int nb_facesF_par_elemF = zone_VFF.zone().nb_faces_elem();
  int nb_facesG_par_elemG = zone_VFG.zone().nb_faces_elem();
  int nb_faces_fines;
  IntTab tab_facesF;
  int nbG;
  int comp;
  int num_elemG;
  int nbfacesF;
  const IntTab& num_facesF = zone_VFF.elem_faces();
  const IntTab& num_facesG = zone_VFG.elem_faces();
  int nbfacesG;
  int num_faceG;
  int dim;
  DoubleVect coord_cgG;
  coord_cgG.resize(dimension);
  const Elem_geom& mon_elemF = zoneF.type_elem();
  int max1 = nb_faceF_(nb_faceG);
  tab_facesF.resize(nb_faceG, max1+1);
  tab_facesF = -1;



  int nbelemsF;
  int nombre_total_elemF =  zone_VFF.nb_elem();

  //Mettre a 0 les valeurs grossieres a corriger
  for (nbG=0; nbG<nb_faceG; nbG++)
    {
      if (nb_faceF_(nbG) > 1E-9)
        {
          if (nb_comp == 1)
            {
              valG(nbG) = 0.;
            }
          else
            for(comp=0; comp<nb_comp; comp++)
              {
                valG(nbG, comp) = 0.;
              }
        }
    }





  if(nb_comp == 1)
    {
      //on ajoute la valeur calculee au centre de gravite
      //de chaque face fine des elemF contenant le centre de gravite
      //pour chaque face grossiere a corriger

      ///////////////////////////////
      ///////////////////////////////
      //ESSAI !!!!!!!!!!!!!!!!
      /////////////////////////////////
      /////////////////////////////////
      //Pour chaque element fin
      for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
        {
          //on recupere le numero de l'element G correspondant
          num_elemG = connect(nbelemsF);
          //Pour chaque face grossiere
          for(nbfacesG=0; nbfacesG<nb_facesG_par_elemG; nbfacesG++)
            {
              //numero global de la face grossiere
              num_faceG = num_facesG(num_elemG, nbfacesG);
              //On regarde si son centre de gravite est dans l'elemF
              for(dim=0; dim<dimension; dim++)
                coord_cgG(dim) = cg_facesG(num_faceG, dim);
              int est_dedans = mon_elemF.contient(coord_cgG, nbelemsF);
              //s'il est dedans
              if ( est_dedans == 1 )
                {
                  //Pour chaque face fine
                  for(nbfacesF=0; nbfacesF<nb_facesF_par_elemF; nbfacesF++)
                    {
                      //numero de la face fine
                      int numF = num_facesF(nbelemsF, nbfacesF);
                      //A-t-elle deja ete comptee?
                      int i = 0;
                      int trouve = 0;

                      while( (tab_facesF(num_faceG, i) != -1) && (i<max1) && (trouve == 0) )
                        {
                          if (tab_facesF(num_faceG, i) == numF)
                            trouve = 1;
                          else
                            i = i+1;
                        }
                      if (trouve == 0)
                        {
                          tab_facesF(num_faceG, i) = numF;
                          if(nb_comp == 1)
                            valG(num_faceG) = valG(num_faceG) + valF(numF);
                          else
                            for(comp=0; comp<nb_comp; comp++)
                              {
                                valG(num_faceG, comp) = valG(num_faceG, comp) + valF(numF, comp);
                              }

                        }
                    }
                }
            }
        }



      //Pour corriger l'inconnue grossiere,
      //on fait la moyenne
      for (nbG=0; nbG<nb_faceG; nbG++)
        {
          nb_faces_fines = nb_faceF_(nbG);
          if (nb_faces_fines > 1E-9)
            {
              if (nb_comp == 1)
                {
                  valG(nbG) = valG(nbG) / nb_faces_fines;
                }
              else
                for(comp=0; comp<nb_comp; comp++)
                  {
                    valG(nbG, comp) = valG(nbG, comp) / nb_faces_fines;
                  }
            }
        }
    }



  else
    {
      //SI INCO VITESSE VEF
      Cerr<<"SI INCO VITESSE VEF A faire !!!!!!!!!!"<<finl;
      exit();
    }
  //Cerr<<"fin de Restriction_face_face::restreindre"<<finl;
}
