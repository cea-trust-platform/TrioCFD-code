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
// File:        Restriction_face_elem.cpp
// Directory:   $TRUST_ROOT/src/Zoom/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Restriction_face_elem.h>

Implemente_instanciable(Restriction_face_elem,"Restriction_face_elem",Restriction_base);

//// printOn
//

Sortie& Restriction_face_elem::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Restriction_face_elem::readOn(Entree& s )
{
  return s ;
}




// Description:
//    calcul du nombre de faces fines contenues dans chaque element grossier
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
void Restriction_face_elem::calculer(const Zone_VF& zone_VFG,
                                     const Zone_VF& zone_VFF,
                                     const IntVect& connect)
{
  //Cerr<<"debut de  Restriction_face_elem::calculer_nb_faceF_dans_chaque_elemG"<<finl;
  int nb_elemG = zone_VFG.nb_elem();
  int max1, max2, max3;
  int nbelemsF;
  int nbelemsG;
  IntVect nb_elemF_;
  IntTab tab_facesF;
  int num_elemG;
  int nbfacesF;
  int nb_facesF_par_elemF = zone_VFF.zone().nb_faces_elem();
  int nombre_total_elemF = connect.size_array();
  nb_elemF_.resize(nb_elemG);
  nb_faceF_.resize(nb_elemG+1);
  nb_faceF_ = 0;
  const IntTab& num_facesF = zone_VFF.elem_faces();
  //On va stocker dans nb_faceF_
  //le nb de faces fines par element grossier
  //et dans la derniere case : le nombre maximal de faces
  //fines que l'on puisse avoir dans un element grossier



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


  //ON COMPLETE L'ATTRIBUT nb_faceF_
  tab_facesF.resize(nb_elemG, max1);
  tab_facesF = -1;
  //Pour chaque element fin
  for (nbelemsF=0; nbelemsF<nombre_total_elemF; nbelemsF++)
    {
      //on recupere le numero de l'element G correspondant
      num_elemG = connect(nbelemsF);
      //Pour chaque face fine
      for(nbfacesF=0; nbfacesF<nb_facesF_par_elemF; nbfacesF++)
        {
          //numero de la face fine
          int numF = num_facesF(nbelemsF, nbfacesF);
          //A-t-elle deja ete comptee?
          int i = 0;
          int trouve = 0;
          while( (tab_facesF(num_elemG, i) != -1) && (i<max1) && (trouve == 0) )
            {
              if (tab_facesF(num_elemG, i) == numF)
                trouve = 1;
              else
                i = i+1;
            }
          if (trouve == 0)
            {
              tab_facesF(num_elemG, i) = numF;
              nb_faceF_(num_elemG) = nb_faceF_(num_elemG)+1;
            }
        }
    }

  //ON TROUVE LE MAX "EXACT" ET ON LE STOCKE DANS LA DERNIERE CASE DU TABLEAU nb_faceF_
  max2 = 0;
  //Pour chaque element grossier
  for (nbelemsG=0; nbelemsG<nb_elemG; nbelemsG++)
    {
      max3 = nb_faceF_(nbelemsG);
      if ( max2 < max3 )
        max2 = max3;
    }
  nb_faceF_(nbelemsG) = max2;

  //Cerr<<"fin de  Restriction_face_elem::calculer_nb_faceF_dans_chaque_elemG"<<finl;
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
void Restriction_face_elem::restreindre(const Zone_VF& zone_VFG,
                                        const Zone_VF& zone_VFF,
                                        const IntVect& connect,
                                        DoubleTab& valG,
                                        const DoubleTab& valF,int nb_comp)
{
  //Cerr<<"debut de  Restriction_face_elem::restreindre"<<finl;
  int nb_facesF_par_elemF = zone_VFF.zone().nb_faces_elem();
  int nb_elemG = zone_VFG.nb_elem();
  int nb_elemF = zone_VFF.nb_elem();
  int nb_faces_fines;
  int max = nb_faceF_(nb_elemG);
  IntTab tab_facesF;
  tab_facesF.resize(nb_elemG, max+1);
  tab_facesF = -1;
  int nbG;
  int nbF;
  int comp;
  int num_elemG;
  int nbfacesF;
  const IntTab& num_facesF = zone_VFF.elem_faces();

  //Mettre a 0 les valeurs grossieres a corriger
  //  Cerr<<"Valeurs grossieres mise a zero"<<finl;
  for (nbG=0; nbG<nb_elemG; nbG++)
    {
      if (nb_faceF_(nbG) > 1E-9)
        {
          //Cerr<<"Numero de l'elemG pour mise a zero : "<<nbG<<finl;
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


  //on ajoute la valeur calculee au centre de gravite
  //de chaque face fine contenue dans l'element grossier a corriger
  for (nbF=0; nbF<nb_elemF; nbF++)
    {
      num_elemG = connect(nbF);
      //Pour chaque face fine
      for(nbfacesF=0; nbfacesF<nb_facesF_par_elemF; nbfacesF++)
        {
          //numero de la face fine
          int numF = num_facesF(nbF, nbfacesF);
          //A-t-elle deja ete comptee?
          int i = 0;
          int trouve = 0;
          while( (tab_facesF(num_elemG, i) != -1) && (i<max+1) && (trouve == 0) )
            {
              if (tab_facesF(num_elemG, i) == numF)
                trouve = 1;
              else
                i = i+1;
            }
          if (trouve == 0)
            {
              tab_facesF(num_elemG, i) = numF;
              if (nb_comp == 1)
                {
                  valG(num_elemG) = valG(num_elemG) + valF(numF);
                }
              else
                for(comp=0; comp<nb_comp; comp++)
                  {
                    valG(num_elemG, comp) = valG(num_elemG, comp) + valF(numF, comp);
                  }
            }
        }
    }


  //Pour corriger l'inconnue grossiere,
  //on fait la moyenne
  for (nbG=0; nbG<nb_elemG; nbG++)
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
  //Cerr<<"fin de Restriction_face_elem::restreindre"<<finl;
}
