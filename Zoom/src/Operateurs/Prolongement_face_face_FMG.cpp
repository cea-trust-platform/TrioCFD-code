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
// File:        Prolongement_face_face_FMG.cpp
// Directory:   $TRUST_ROOT/src/Zoom/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_face_face_FMG.h>
#include <Domaine.h>
#include <Zone_VF.h>

Implemente_instanciable(Prolongement_face_face_FMG,"Prolongement_face_face_FMG",Prolongement_base);

//// printOn
//

Sortie& Prolongement_face_face_FMG::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_face_face_FMG::readOn(Entree& s )
{
  return s ;
}




// Description:
//    Prolongement du centre de gravite des faces grossieres
//    au centre de gravite de TOUTES les faces fines
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
void Prolongement_face_face_FMG::prolonger(Zone_VF& zone_VFG,
                                           Zone_VF& zone_VFF,
                                           const Frontiere& frontF,
                                           IntVect& connect,
                                           const DoubleTab& valG, DoubleTab& tab,
                                           int nb_compo)
{
  //Cerr<<"debut de Prolongement_face_face_FMG::prolonger"<<finl;
  int nb_elemF = connect.size_array();
  Zone& zonef = zone_VFF.zone();
  Zone& zoneg = zone_VFG.zone();
  Domaine& domg = zoneg.domaine();
  int nb_faces_elem_fin = zonef.nb_faces_elem();
  int nb_faces_elem_gros = zoneg.nb_faces_elem();
  const DoubleTab& cg_face_fine = zone_VFF.xv();
  const DoubleTab& som_gros = domg.coord_sommets();
  const IntTab& faces_som = zone_VFG.face_sommets();
  int trouve=0, i, nbj;
  IntTab& num_face_elemF = zone_VFF.elem_faces();
  IntTab& num_face_elemG = zone_VFG.elem_faces();
  DoubleVect coord_vect_faceF;
  DoubleVect coord_vect_faceG;
  coord_vect_faceF.resize(2);
  coord_vect_faceG.resize(2);
  tab = 0.;
  int nbElemF, num_elemG, nbFacesF, j=0, compo;
  double det;
  double epsilon = 1.e-9;
  double sign;

  double p0, p1, p2;
  int som0, som1, som2;
  double prod;



  //Pour chaque elemF
  for(nbElemF=0; nbElemF<nb_elemF; nbElemF++)
    {
      num_elemG = connect(nbElemF);
      //Pour chaque faceF
      for(nbFacesF=0; nbFacesF<nb_faces_elem_fin; nbFacesF++)
        {
          i = num_face_elemF(nbElemF, nbFacesF);
          //on regarde si elle appartient a une faceG de num_elemG

          if(dimension == 2)
            {
              trouve = 0;
              for (nbj=0; nbj<nb_faces_elem_gros; nbj++)
                {
                  j = num_face_elemG(num_elemG, nbj);
                  det = (cg_face_fine(i,1)-som_gros(faces_som(j,0),1)) * (som_gros(faces_som(j,1),0)-som_gros(faces_som(j,0),0)) - (cg_face_fine(i,0)-som_gros(faces_som(j,0),0)) * (som_gros(faces_som(j,1),1)-som_gros(faces_som(j,0),1));


                  if(fabs(som_gros(faces_som(j,1),0)-som_gros(faces_som(j,0),0))<epsilon)
                    {
                      sign = (cg_face_fine(i,1)-som_gros(faces_som(j,0),1))/(som_gros(faces_som(j,1),1)-som_gros(faces_som(j,0),1));
                    }
                  else
                    {
                      sign = (cg_face_fine(i,0)-som_gros(faces_som(j,0),0))/(som_gros(faces_som(j,1),0)-som_gros(faces_som(j,0),0));
                    }
                  ///////////
                  if (fabs(det) < epsilon  &&  sign >0  && sign <1)
                    {
                      //on met directement la valG dans la valF
                      if (nb_compo == 1)
                        {
                          tab(i) = valG(j);
                        }
                      else
                        for(compo = 0; compo<nb_compo; compo++)
                          tab(i, compo) = valG(j, compo);
                      trouve = 1;
                      break;
                    }
                }
            }



          else if(dimension == 3)
            {
              som0 = faces_som(j,0);
              som1 = faces_som(j,1);
              som2 = faces_som(j,2);


              //Cerr<<"EST-IL DANS LE PREMIER TRIANGLE ??"<<finl;
              //Attention
              // Il faut donc determiner le sens (trigo ou anti trigo)
              //pour la numerotation :
              // Calcul de prod = 01 vectoriel 02 selon z
              // prod > 0 : sens trigo
              // prod < 0 : sens anti trigo
              prod = (som_gros(som1,0)-som_gros(som0,0))*(som_gros(som2,1)-som_gros(som0,1)) - (som_gros(som1,1)-som_gros(som0,1))*(som_gros(som2,0)-som_gros(som0,0));
              if (prod >= 0)
                prod = 1;
              else
                prod = -1;
              // Calcul de p0 = 0M vectoriel 1M selon z
              p0 = (cg_face_fine(i,0)-som_gros(som0,0))*(cg_face_fine(i,1)-som_gros(som1,1)) -
                   (cg_face_fine(i,1)-som_gros(som0,1))*(cg_face_fine(i,0)-som_gros(som1,0));
              p0 = p0 * prod;
              // Calcul de p1 = 1M vectoriel 2M selon z
              p1 = (cg_face_fine(i,0)-som_gros(som1,0))*(cg_face_fine(i,1)-som_gros(som2,1)) -
                   (cg_face_fine(i,1)-som_gros(som1,1))*(cg_face_fine(i,0)-som_gros(som2,0));
              p1 = p1 * prod;
              // Calcul de p2 = 2M vectoriel 0M selon z
              p2 = (cg_face_fine(i,0)-som_gros(som2,0))*(cg_face_fine(i,1)-som_gros(som0,1)) -
                   (cg_face_fine(i,1)-som_gros(som2,1))*(cg_face_fine(i,0)-som_gros(som0,0));
              p2 = p2 * prod;
              if (p0 >= 0 && p1 >= 0 && p2 >= 0)
                {
                  //Cerr<<"EST SUR UNE FACE DE BORD"<<finl;
                  if (nb_compo == 1)
                    {
                      tab(i) = valG(j);
                    }
                  else
                    for(compo = 0; compo<nb_compo; compo++)
                      tab(i, compo) = valG(j, compo);
                  trouve = 1;
                  break;
                }


              if (trouve == 0)
                {
                  //Cerr<<"EST-IL DANS LE DEUXIEME TRIANGLE ??"<<finl;
                  som0 = faces_som(j,3);
                  som1 = faces_som(j,2);
                  som2 = faces_som(j,1);


                  //Attention
                  // Il faut donc determiner le sens (trigo ou anti trigo)
                  //pour la numerotation :
                  // Calcul de prod = 01 vectoriel 02 selon z
                  // prod > 0 : sens trigo
                  // prod < 0 : sens anti trigo
                  prod = (som_gros(som1,0)-som_gros(som0,0))*(som_gros(som2,1)-som_gros(som0,1)) - (som_gros(som1,1)-som_gros(som0,1))*(som_gros(som2,0)-som_gros(som0,0));
                  if (prod >= 0)
                    prod = 1;
                  else
                    prod = -1;
                  // Calcul de p0 = 0M vectoriel 1M selon z
                  p0 = (cg_face_fine(i,0)-som_gros(som0,0))*(cg_face_fine(i,1)-som_gros(som1,1)) -
                       (cg_face_fine(i,1)-som_gros(som0,1))*(cg_face_fine(i,0)-som_gros(som1,0));
                  p0 = p0 * prod;
                  // Calcul de p1 = 1M vectoriel 2M selon z
                  p1 = (cg_face_fine(i,0)-som_gros(som1,0))*(cg_face_fine(i,1)-som_gros(som2,1)) -
                       (cg_face_fine(i,1)-som_gros(som1,1))*(cg_face_fine(i,0)-som_gros(som2,0));
                  p1 = p1 * prod;
                  // Calcul de p2 = 2M vectoriel 0M selon z
                  p2 = (cg_face_fine(i,0)-som_gros(som2,0))*(cg_face_fine(i,1)-som_gros(som0,1)) -
                       (cg_face_fine(i,1)-som_gros(som2,1))*(cg_face_fine(i,0)-som_gros(som0,0));
                  p2 = p2 * prod;
                  if (p0 >= 0 && p1 >= 0 && p2 >= 0)
                    {
                      //Cerr<<"EST SUR UNE FACE DE BORD"<<finl;
                      if (nb_compo == 1)
                        {
                          tab(i) = valG(j);
                        }
                      else
                        for(compo = 0; compo<nb_compo; compo++)
                          tab(i, compo) = valG(j, compo);
                      trouve = 1;
                      break;
                    }
                }
            }


          if (trouve == 0  &&  tab(i) == 0 )
            {
              for (nbj=0; nbj<nb_faces_elem_gros; nbj++)
                {
                  j = num_face_elemG(num_elemG, nbj);

                  if (nb_compo == 1)
                    tab(i) += valG(j);
                  else
                    for(compo = 0; compo<nb_compo; compo++)
                      tab(i, compo) += valG(j, compo);
                }
              tab(i) /= nb_faces_elem_gros;
            }
        }
    }
  //Cerr<<"fin de Prolongement_face_face_FMG::prolonger"<<finl;
}



//NE FAIT RIEN
void Prolongement_face_face_FMG::calculer(Zone_VF& zonef,
                                          Zone_VF& zoneg,
                                          IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}


