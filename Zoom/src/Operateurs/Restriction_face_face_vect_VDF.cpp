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
// File:        Restriction_face_face_vect_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Restriction_face_face_vect_VDF.h>

Implemente_instanciable(Restriction_face_face_vect_VDF,"Restriction_face_face_vect_VDF",Restriction_base);

//// printOn
//

Sortie& Restriction_face_face_vect_VDF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Restriction_face_face_vect_VDF::readOn(Entree& s )
{
  return s ;
}




// Description:
//    calcul du nombre de faces fines contenant le centre
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
void Restriction_face_face_vect_VDF::calculer(const Zone_VF& zone_VFG,
                                              const Zone_VF& zone_VFF,
                                              const IntVect& connect)
{
  //Cerr<<"debut de  Restriction_face_face_vect_VDF::calculer"<<finl;



  //coord du centre de gravite des faces
  const DoubleTab& cg_face_fine = zone_VFF.xv();
  const DoubleTab& cg_face_gros = zone_VFG.xv();
  int j;
  int orien;
  IntTab trouve_facef_incluse;


  int nb_elemG = zone_VFG.nb_elem();
  int nb_faceG = zone_VFG.nb_faces();

  trouve_facef_incluse.resize(nb_faceG);
  trouve_facef_incluse = 0;

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
  DoubleVect coord_cgG;
  coord_cgG.resize(dimension);
  int num_faceG;

  //ATTENTION, JE DOIS TYPER LES ZONES EN VDF
  const Zone_VDF& zone_VDFG = ref_cast( Zone_VDF,zone_VFG);
  const Zone_VDF& zone_VDFF = ref_cast( Zone_VDF,zone_VFF);
  const IntVect& oriG = zone_VDFG.orientation();
  const IntVect& oriF = zone_VDFF.orientation();


  //On va stocker dans nb_faceF_
  //pour chaque faceG,
  //le nb de faces fines
  //qui sont de meme orientation que la faceG
  //et dont l'elemF est dans l'elemG
  //MODIF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //S'IL EXISTE UNE FACE FINE INCLUSE DANS LA FACE G
  //ON NE COMPTERA QUE LES FACES FINES QUI SONT DANS CETTE FACE G

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
  //la face grossiere appartient au plus a 2 elements grossiers
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


          num_faceG = num_facesG(num_elemG, nbfacesG);
          //Pour chaque face fine
          for(nbfacesF=0; nbfacesF<nb_facesF_par_elemF; nbfacesF++)
            {
              //numero de la face fine
              int numF = num_facesF(nbelemsF, nbfacesF);
              if(trouve_facef_incluse(num_faceG) == 0)
                {
                  ;
                  //on regarde si la face fine est de meme orientation
                  if(oriF(numF) == oriG(num_faceG))
                    {
                      orien = oriF(numF);
                      //et si elle est incluse dans la face grossiere
                      if((fabs(cg_face_fine(numF, orien)-cg_face_gros(num_faceG,orien))) < 1.e-9)
                        {
                          trouve_facef_incluse(num_faceG) = 1;
                          for(j=0; j<max1; j++)
                            tab_facesF(num_faceG, j) = -1;
                          tab_facesF(num_faceG, 0) = numF;
                          nb_faceF_(num_faceG) = 1;
                        }
                      //si elle N'est  PAS incluse dans la face grossiere
                      else
                        {
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
              else
                {
                  //Cerr<<"si on a deja trouve une face incluse"<<finl;
                  //Si la faceF est de meme orientation
                  //et si elle est incluse, on la compte
                  if(oriF(numF) == oriG(num_faceG))
                    {
                      orien = oriF(numF);
                      if((fabs(cg_face_fine(numF, orien)-cg_face_gros(num_faceG,orien ))) < 1.e-9)
                        {
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

  //Cerr<<"fin de  Restriction_face_face_vect_VDF::calculer"<<finl;
}




// Description:
//    Restriction de l'inconnue fine sur le domaine grossier
//    pour inconnue VITESSE
//    en VDF-VDF
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
void Restriction_face_face_vect_VDF::restreindre(const Zone_VF& zone_VFG,
                                                 const Zone_VF& zone_VFF,
                                                 const IntVect& connect,
                                                 DoubleTab& valG,
                                                 const DoubleTab& valF,int nb_comp)
{
  //Cerr<<"debut de  Restriction_face_face_vect_VDF::restreindre"<<finl;


  //coord du centre de gravite des faces fines
  const DoubleTab& cg_face_fine = zone_VFF.xv();
  int j;
  int orien;
  IntTab trouve_facef_incluse;

  const DoubleTab& cg_facesG = zone_VFG.xv();
  int nb_faceG = zone_VFG.nb_faces();

  trouve_facef_incluse.resize(nb_faceG);
  trouve_facef_incluse = 0;

  int nb_facesF_par_elemF = zone_VFF.zone().nb_faces_elem();
  int nb_facesG_par_elemG = zone_VFG.zone().nb_faces_elem();
  int nb_faces_fines;
  IntTab tab_facesF;
  int nbG;
  int num_elemG;
  int nbfacesF;
  const IntTab& num_facesF = zone_VFF.elem_faces();
  const IntTab& num_facesG = zone_VFG.elem_faces();
  int nbfacesG;
  int num_faceG;

  int max1 = nb_faceF_(nb_faceG);
  tab_facesF.resize(nb_faceG, max1+1);
  tab_facesF = -1;


  //ATTENTION, JE DOIS TYPER LES ZONES EN VDF
  const Zone_VDF& zone_VDFG = ref_cast( Zone_VDF,zone_VFG);
  const Zone_VDF& zone_VDFF = ref_cast( Zone_VDF,zone_VFF);
  const IntVect& oriG = zone_VDFG.orientation();
  const IntVect& oriF = zone_VDFF.orientation();




  int nbelemsF;
  int nombre_total_elemF =  zone_VFF.nb_elem();

  //Mettre a 0 les valeurs grossieres a corriger
  for (nbG=0; nbG<nb_faceG; nbG++)
    {
      if (nb_faceF_(nbG) > 1E-9)
        {
          valG(nbG) = 0;
        }
    }



  //on ajoute la valeur de la vitesse au centre de gravite des faces
  //fines dont l'orientation est la meme que la face grossiere
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      //on recupere le numero de l'element G correspondant
      num_elemG = connect(nbelemsF);
      //Pour chaque face grossiere
      for(nbfacesG=0; nbfacesG<nb_facesG_par_elemG; nbfacesG++)
        {
          num_faceG = num_facesG(num_elemG, nbfacesG);
          //Pour chaque face fine
          for(nbfacesF=0; nbfacesF<nb_facesF_par_elemF; nbfacesF++)
            {
              //numero de la face fine
              int numF = num_facesF(nbelemsF, nbfacesF);
              if(trouve_facef_incluse(num_faceG) == 0)
                {
                  //on regarde si la face fine est de meme orientation
                  if(oriF(numF) == oriG(num_faceG))
                    {
                      orien = oriF(numF);
                      //et si elle est incluse dans la face grossiere
                      if(fabs(cg_face_fine(numF, orien)-cg_facesG(num_faceG,orien )) < 1.e-9)
                        {
                          trouve_facef_incluse(num_faceG) = 1;
                          for(j=0; j<max1; j++)
                            tab_facesF(num_faceG, j) = -1;
                          tab_facesF(num_faceG, 0) = numF;
                          valG(num_faceG) = 0.;
                          valG(num_faceG) = valG(num_faceG) + valF(numF);
                        }
                      //si elle N'est  PAS incluse dans la face grossiere
                      else
                        {
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
                              valG(num_faceG) = valG(num_faceG) + valF(numF);
                            }
                        }
                    }
                }
              else
                {
                  //Si la faceF est de meme orientation
                  //et si elle est incluse, on la compte
                  if(oriF(numF) == oriG(num_faceG))
                    {
                      orien = oriF(numF);
                      if(fabs(cg_face_fine(numF, orien)-cg_facesG(num_faceG,orien )) < 1.e-9)
                        {
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
                              valG(num_faceG) = valG(num_faceG) + valF(numF);
                            }
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
          valG(nbG) = valG(nbG) / nb_faces_fines;
        }
    }


  //Cerr<<"fin de Restriction_face_face_vect_VDF::restreindre"<<finl;
}
