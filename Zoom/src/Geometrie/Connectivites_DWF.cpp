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
// File:        Connectivites_DWF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Geometrie
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Connectivites_DWF.h>

#define PRECISION 1.e-9

Implemente_instanciable(Connectivites_DWF,"Connectivites_DWF",Connectivites_base);

//// printOn
//

Sortie& Connectivites_DWF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Connectivites_DWF::readOn(Entree& s )
{
  return s ;
}


void Connectivites_DWF::calculer_connectivites_face_face(Zone_VF& zonef,
                                                         Zone_VF& zoneg,
                                                         Domaine& domg)
{
  //Cerr<<"debut de Connectivites_base::calculer_connectivites_face_face"<<finl;
  int nb_zone = domg.nb_zones();
  const Elem_geom& type_elem = domg.zone(nb_zone-1).type_elem();
  const int nb_faces_bordF = zonef.nb_faces_bord();
  const int nb_faces_bordG = zoneg.nb_faces_bord();

  //coord du centre de gravite des faces fines
  const DoubleTab& cg_face_fine = zonef.xv();

  const int prem_face_bord_gros  =  zoneg.premiere_face_bord();
  const int der_face_bord_gros  =  prem_face_bord_gros+zoneg.nb_faces_bord();
  const int prem_face_int_gros  =  zoneg.premiere_face_int();
  const int der_face_int_gros  = prem_face_int_gros + zoneg.nb_faces_internes();

  //coordonnees des sommets grossiers
  const DoubleTab& som_gros = domg.coord_sommets();
  const IntTab& faces_som = zoneg.face_sommets();

  const double epsilon  = PRECISION;



  connect_faceF_faceG.resize(nb_faces_bordF);
  connect_faceF_faceG = -1;
  connect_faceG_faceF.resize(nb_faces_bordG);
  connect_faceG_faceF = -1;

  double sign;
  double det;

  int trouve;
  double prod, p0, p1, p2;
  int som0;
  int som1;
  int som2;
  int j;

  if (dimension == 2)
    {
      for (int i=0; i<nb_faces_bordF; i++) // boucle faces fines de bords
        {
          //On regarde si le centre de gravite fin est un sommet grossier
          trouve = 0;
          j = prem_face_bord_gros;
          while((j<der_face_bord_gros) && (trouve == 0))
            {
              //On regarde si le centre de gravite fin est un sommet grossier
              //d'une face grossiere de bord
              if(((fabs(cg_face_fine(i,0)-som_gros(faces_som(j,0),0))<PRECISION) && (fabs(cg_face_fine(i,1)-som_gros(faces_som(j,0),1))<PRECISION)) || ((fabs(cg_face_fine(i,0)-som_gros(faces_som(j,1),0))<PRECISION) && (fabs(cg_face_fine(i,1)-som_gros(faces_som(j,1),1))<PRECISION)))
                {
                  trouve = 1;
                  connect_faceF_faceG(i)=j;
                  if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                  else connect_faceG_faceF(j) = i;
                }
              else
                j = j+1;
            }


          if(trouve == 0)
            {
              //si ce n'est pas le sommet d'une faceG de bord
              j = prem_face_int_gros;
              while((j<der_face_int_gros) && (trouve == 0))
                {
                  //On regarde si le centre de gravite fin est un sommet grossier
                  //d'une face grossiere interne
                  if(((fabs(cg_face_fine(i,0)-som_gros(faces_som(j,0),0))<PRECISION) && (fabs(cg_face_fine(i,1)-som_gros(faces_som(j,0),1))<PRECISION)) || ((fabs(cg_face_fine(i,0)-som_gros(faces_som(j,1),0))<PRECISION) && (fabs(cg_face_fine(i,1)-som_gros(faces_som(j,1),1))<PRECISION)))
                    {
                      trouve = 1;
                      connect_faceF_faceG(i)=j;
                      if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                      else connect_faceG_faceF(j) = i;
                    }
                  else
                    j = j+1;
                }
            }

          if(trouve == 0)
            {
              //si le centre de graviteF n'est pas un sommet grossier
              for (j=prem_face_bord_gros; j<der_face_bord_gros; j++) // boucle faces grosses de bords
                {
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
                      connect_faceF_faceG(i)=j;
                      if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                      else connect_faceG_faceF(j) = i;

                      trouve = 1;
                      break;
                    }
                }

              if (trouve == 0)
                {
                  for (j=prem_face_int_gros; j<der_face_int_gros; j++) // boucle faces grosses internes
                    {
                      det = (cg_face_fine(i,1)-som_gros(faces_som(j,0),1)) * (som_gros(faces_som(j,1),0)-som_gros(faces_som(j,0),0)) - (cg_face_fine(i,0)-som_gros(faces_som(j,0),0)) * (som_gros(faces_som(j,1),1)-som_gros(faces_som(j,0),1));

                      if(fabs(som_gros(faces_som(j,1),0)-som_gros(faces_som(j,0),0))<epsilon)
                        {
                          sign = (cg_face_fine(i,1)-som_gros(faces_som(j,0),1))/(som_gros(faces_som(j,1),1)-som_gros(faces_som(j,0),1));
                        }
                      else
                        {
                          sign = (cg_face_fine(i,0)-som_gros(faces_som(j,0),0))/(som_gros(faces_som(j,1),0)-som_gros(faces_som(j,0),0));
                        }
                      ////////
                      if (fabs(det) < epsilon  &&  sign >0  && sign <1)
                        {
                          connect_faceF_faceG(i)=j;
                          if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                          else connect_faceG_faceF(j) = i;

                          trouve = 1;
                          break;
                        }
                    }
                }
            }
        }
    }
  else if (dimension == 3)
    {
      DoubleVect pointM;
      pointM.resize(3);


      for (int i=0; i<nb_faces_bordF; i++) // boucle faces fines de bords
        {
          trouve = 0;
          for (j=prem_face_bord_gros; j<der_face_bord_gros; j++) // boucle faces grosses de bords
            {
              //Cerr<<"POUR CHAQUE FACE GROSSIERE DE BORD"<<finl;
              som0 = faces_som(j,0);
              som1 = faces_som(j,1);
              som2 = faces_som(j,2);

              // On regarde tout d'abord si le point cherche n'est pas un des
              // sommets du triangle

              if( ((fabs(som_gros(som0,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som0,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som0,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som1,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som1,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som1,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som2,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som2,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som2,2)-cg_face_fine(i,2))<PRECISION)) )
                {
                  connect_faceF_faceG(i)=j;
                  if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                  else connect_faceG_faceF(j) = i;

                  trouve = 1;
                  break;
                }
              if (trouve == 0)
                {

                  //Attention
                  // Il faut donc determiner le sens (trigo ou anti trigo)
                  //pour la numerotation :


                  //si les z sont les memes :
                  if((fabs(som_gros(som0,2)-som_gros(som1,2))<PRECISION) && (fabs(som_gros(som1,2)-som_gros(som2,2))<PRECISION))
                    {
                      //si le centre de gravite est dans le plan forme des 3 sommets
                      if(fabs(cg_face_fine(i,2)- som_gros(som0,2))<PRECISION)
                        {
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
                        }
                      else
                        {
                          p0 = -1;
                          p1 = -1;
                          p2 = -1;
                        }
                    }
                  else
                    //si les y sont les memes :
                    if((fabs(som_gros(som0,1) - som_gros(som1,1))<PRECISION) && (fabs(som_gros(som1,1) - som_gros(som2,1))<PRECISION))
                      {
                        //Cerr<<"les y sont les memes"<<finl;
                        //si le centre de gravite est dans le plan forme des 3 sommets
                        if(fabs(cg_face_fine(i,1) - som_gros(som0,1))<PRECISION)
                          {
                            // Calcul de prod = 01 vectoriel 02 selon y
                            // prod > 0 : sens trigo
                            // prod < 0 : sens anti trigo
                            prod = (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,0)-som_gros(som0,0)) - (som_gros(som1,0)-som_gros(som0,0))*(som_gros(som2,2)-som_gros(som0,2));
                            if (prod >= 0)
                              prod = 1;
                            else
                              prod = -1;
                            // Calcul de p0 = 0M vectoriel 1M selon y
                            p0 = (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,0)-som_gros(som1,0)) -
                                 (cg_face_fine(i,0)-som_gros(som0,0))*(cg_face_fine(i,2)-som_gros(som1,2));
                            p0 = p0 * prod;
                            // Calcul de p1 = 1M vectoriel 2M selon y
                            p1 = (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,0)-som_gros(som2,0)) -
                                 (cg_face_fine(i,0)-som_gros(som1,0))*(cg_face_fine(i,2)-som_gros(som2,2));
                            p1 = p1 * prod;
                            // Calcul de p2 = 2M vectoriel 0M selon y
                            p2 = (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,0)-som_gros(som0,0)) -
                                 (cg_face_fine(i,0)-som_gros(som2,0))*(cg_face_fine(i,2)-som_gros(som0,2));
                            p2 = p2 * prod;
                          }
                        else
                          {
                            p0 = -1;
                            p1 = -1;
                            p2 = -1;
                          }
                      }
                    else
                      //si les x sont les memes :
                      if((fabs(som_gros(som0,0) - som_gros(som1,0))<PRECISION) && (fabs(som_gros(som1,0) - som_gros(som2,0))<PRECISION))
                        {
                          //si le centre de gravite est dans le plan forme des 3 sommets
                          /*                           if(fabs(cg_face_fine(i,0)-0.6) < 0.01)
                                                       {
                                                       Cerr<<"centre de gravite "<<cg_face_fine(i,0) <<  finl;
                                                       Cerr<<"som gros "<<som_gros(som0,0) <<  finl;
                                                       }
                          */                          if(fabs(cg_face_fine(i,0) - som_gros(som0,0))<PRECISION)
                            {
                              //          Cerr<<"le centre de gravite est dans le plan"<<finl;
                              // Calcul de prod = 01 vectoriel 02 selon x
                              // prod > 0 : sens trigo
                              // prod < 0 : sens anti trigo
                              prod = (som_gros(som1,1)-som_gros(som0,1))*(som_gros(som2,2)-som_gros(som0,2)) -
                                     (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,1)-som_gros(som0,1));
                              if (prod >= 0)
                                prod = 1;
                              else
                                prod = -1;
                              // Calcul de p0 = 0M vectoriel 1M selon x
                              p0 = (cg_face_fine(i,1)-som_gros(som0,1))*(cg_face_fine(i,2)-som_gros(som1,2)) -
                                   (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,1)-som_gros(som1,1));
                              p0 = p0 * prod;
                              // Calcul de p1 = 1M vectoriel 2M selon x
                              p1 = (cg_face_fine(i,1)-som_gros(som1,1))*(cg_face_fine(i,2)-som_gros(som2,2)) -
                                   (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,1)-som_gros(som2,1));
                              p1 = p1 * prod;
                              // Calcul de p2 = 2M vectoriel 0M selon x
                              p2 = (cg_face_fine(i,1)-som_gros(som2,1))*(cg_face_fine(i,2)-som_gros(som0,2)) -
                                   (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,1)-som_gros(som0,1));
                              p2 = p2 * prod;
                            }
                          else
                            {
                              p0 = -1;
                              p1 = -1;
                              p2 = -1;
                            }
                        }
                      else
                        {
                          p0 = -1;
                          p1 = -1;
                          p2 = -1;
                        }





                  if (p0 >= 0 && p1 >= 0 && p2 >= 0)
                    {
                      //Cerr<<"EST SUR UNE FACE DE BORD"<<finl;
                      connect_faceF_faceG(i)=j;
                      if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                      else connect_faceG_faceF(j) = i;

                      trouve = 1;

                      break;
                    }



                  else  if ( type_elem->que_suis_je() != "Tetraedre" )
                    {
                      //Cerr<<"EST-IL DANS LE DEUXIEME TRIANGLE ??"<<finl;
                      som0 = faces_som(j,3);
                      som1 = faces_som(j,2);
                      som2 = faces_som(j,1);

                      // On regarde tout d'abord si le point cherche n'est pas un des
                      // sommets du triangle
                      if( ((fabs(som_gros(som0,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som0,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som0,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som1,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som1,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som1,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som2,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som2,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som2,2)-cg_face_fine(i,2))<PRECISION)) )
                        {
                          connect_faceF_faceG(i)=j;
                          if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                          else connect_faceG_faceF(j) = i;

                          trouve = 1;
                          break;
                        }
                      if (trouve == 0)
                        {

                          //Attention
                          // Il faut donc determiner le sens (trigo ou anti trigo)
                          //pour la numerotation :


                          //si les z sont les memes et si le centre de gravite est dans le plan forme des 3 sommets:
                          if((fabs(som_gros(som0,2) - som_gros(som1,2))<PRECISION) && (fabs(som_gros(som1,2) - som_gros(som2,2))<PRECISION))
                            {
                              //si le centre de gravite fin est dans le plan forme des 3 sommets
                              if(fabs(cg_face_fine(i,2) - som_gros(som0,2))<PRECISION)
                                {
                                  //Cerr<<" z sont les memes et centre de gravite fin est dans le plan forme des 3 sommets"<<finl;
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
                                }
                              else
                                {
                                  p0 = -1;
                                  p1 = -1;
                                  p2 = -1;
                                }
                            }
                          else
                            //si les y sont les memes :
                            if((fabs(som_gros(som0,1) - som_gros(som1,1))<PRECISION) && (fabs(som_gros(som1,1) - som_gros(som2,1))<PRECISION))
                              {
                                //si le centre de gravite fin est dans le plan forme des 3 sommets
                                if(fabs(cg_face_fine(i,1) - som_gros(som0,1))<PRECISION)
                                  {
                                    //Cerr<<" y sont les memes et centre de gravite fin est dans le plan forme des 3 sommets"<<finl;
                                    // Calcul de prod = 01 vectoriel 02 selon y
                                    // prod > 0 : sens trigo
                                    // prod < 0 : sens anti trigo
                                    prod = (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,0)-som_gros(som0,0)) - (som_gros(som1,0)-som_gros(som0,0))*(som_gros(som2,2)-som_gros(som0,2));
                                    if (prod >= 0)
                                      prod = 1;
                                    else
                                      prod = -1;
                                    // Calcul de p0 = 0M vectoriel 1M selon y
                                    p0 = (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,0)-som_gros(som1,0)) -
                                         (cg_face_fine(i,0)-som_gros(som0,0))*(cg_face_fine(i,2)-som_gros(som1,2));
                                    p0 = p0 * prod;
                                    // Calcul de p1 = 1M vectoriel 2M selon y
                                    p1 = (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,0)-som_gros(som2,0)) -
                                         (cg_face_fine(i,0)-som_gros(som1,0))*(cg_face_fine(i,2)-som_gros(som2,2));
                                    p1 = p1 * prod;
                                    // Calcul de p2 = 2M vectoriel 0M selon y
                                    p2 = (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,0)-som_gros(som0,0)) -
                                         (cg_face_fine(i,0)-som_gros(som2,0))*(cg_face_fine(i,2)-som_gros(som0,2));
                                    p2 = p2 * prod;
                                  }
                                else
                                  {
                                    p0 = -1;
                                    p1 = -1;
                                    p2 = -1;
                                  }
                              }
                            else
                              //si les x sont les memes :
                              if((fabs(som_gros(som0,0) - som_gros(som1,0))<PRECISION) && (fabs(som_gros(som1,0) - som_gros(som2,0))<PRECISION))
                                {
                                  //si le centre de gravite fin est dans le plan forme des 3 sommets
                                  if(fabs(cg_face_fine(i,0) - som_gros(som0,0))<PRECISION)
                                    {
                                      //Cerr<<" x sont les memes et centre de gravite fin est dans le plan forme des 3 sommets"<<finl;
                                      // Calcul de prod = 01 vectoriel 02 selon x
                                      // prod > 0 : sens trigo
                                      // prod < 0 : sens anti trigo
                                      prod = (som_gros(som1,1)-som_gros(som0,1))*(som_gros(som2,2)-som_gros(som0,2)) -
                                             (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,1)-som_gros(som0,1));
                                      if (prod >= 0)
                                        prod = 1;
                                      else
                                        prod = -1;
                                      // Calcul de p0 = 0M vectoriel 1M selon x
                                      p0 = (cg_face_fine(i,1)-som_gros(som0,1))*(cg_face_fine(i,2)-som_gros(som1,2)) -
                                           (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,1)-som_gros(som1,1));
                                      p0 = p0 * prod;
                                      // Calcul de p1 = 1M vectoriel 2M selon x
                                      p1 = (cg_face_fine(i,1)-som_gros(som1,1))*(cg_face_fine(i,2)-som_gros(som2,2)) -
                                           (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,1)-som_gros(som2,1));
                                      p1 = p1 * prod;
                                      // Calcul de p2 = 2M vectoriel 0M selon x
                                      p2 = (cg_face_fine(i,1)-som_gros(som2,1))*(cg_face_fine(i,2)-som_gros(som0,2)) -
                                           (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,1)-som_gros(som0,1));
                                      p2 = p2 * prod;
                                    }
                                  else
                                    {
                                      p0 = -1;
                                      p1 = -1;
                                      p2 = -1;
                                    }
                                }


                              else
                                {
                                  p0 = -1;
                                  p1 = -1;
                                  p2 = -1;
                                }




                          if (p0 >= 0 && p1 >= 0 && p2 >= 0)
                            {
                              //Cerr<<"EST SUR UNE FACE DE BORD"<<finl;
                              connect_faceF_faceG(i)=j;
                              if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                              else connect_faceG_faceF(j) = i;

                              trouve = 1;
                              break;
                            }
                        }
                    }
                }
            }

          if (trouve == 0)
            {
              for (j=prem_face_int_gros; j<der_face_int_gros; j++) // boucle faces grosses internes
                {
                  //Cerr<<"POUR CHAQUE FACE GROSSIERE INTERNE"<<finl;
                  som0 = faces_som(j,0);
                  som1 = faces_som(j,1);
                  som2 = faces_som(j,2);

                  // On regarde tout d'abord si le point cherche n'est pas un des
                  // sommets du triangle
                  if( ((fabs(som_gros(som0,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som0,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som0,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som1,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som1,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som1,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som2,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som2,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som2,2)-cg_face_fine(i,2))<PRECISION)) )
                    {
                      connect_faceF_faceG(i)=j;
                      if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                      else connect_faceG_faceF(j) = i;

                      trouve = 1;
                      break;
                    }
                  if (trouve == 0)
                    {


                      //Attention
                      // Il faut donc determiner le sens (trigo ou anti trigo)
                      //pour la numerotation :


                      //si les z sont les memes :
                      if((fabs(som_gros(som0,2) - som_gros(som1,2))<PRECISION) && (fabs(som_gros(som1,2) - som_gros(som2,2))<PRECISION))
                        {
                          //Cerr<<" z sont les memes"<<finl;
                          //si le centre de gravite est dans le plan forme par les 3 sommets
                          if(fabs(cg_face_fine(i,2) - som_gros(som0,2))<PRECISION)
                            {
                              //Cerr<<" z sont les memes et centre de gravite fin est dans le plan forme des 3 sommets"<<finl;
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
                            }
                          else
                            {
                              p0 = -1;
                              p1 = -1;
                              p2 = -1;
                            }
                        }
                      else
                        //si les y sont les memes :
                        if((fabs(som_gros(som0,1) - som_gros(som1,1))<PRECISION) && (fabs(som_gros(som1,1) - som_gros(som2,1))<PRECISION))
                          {
                            //Cerr<<" y sont les memes"<<finl;
                            //si le centre de gravite est dans le plan forme par les 3 sommets
                            if(fabs(cg_face_fine(i,1) - som_gros(som0,1))<PRECISION)
                              {
                                //Cerr<<" y sont les memes et centre de gravite fin est dans le plan forme des 3 sommets"<<finl;
                                // Calcul de prod = 01 vectoriel 02 selon y
                                // prod > 0 : sens trigo
                                // prod < 0 : sens anti trigo
                                prod = (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,0)-som_gros(som0,0)) - (som_gros(som1,0)-som_gros(som0,0))*(som_gros(som2,2)-som_gros(som0,2));
                                if (prod >= 0)
                                  prod = 1;
                                else
                                  prod = -1;
                                // Calcul de p0 = 0M vectoriel 1M selon y
                                p0 = (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,0)-som_gros(som1,0)) -
                                     (cg_face_fine(i,0)-som_gros(som0,0))*(cg_face_fine(i,2)-som_gros(som1,2));
                                p0 = p0 * prod;
                                // Calcul de p1 = 1M vectoriel 2M selon y
                                p1 = (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,0)-som_gros(som2,0)) -
                                     (cg_face_fine(i,0)-som_gros(som1,0))*(cg_face_fine(i,2)-som_gros(som2,2));
                                p1 = p1 * prod;
                                // Calcul de p2 = 2M vectoriel 0M selon y
                                p2 = (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,0)-som_gros(som0,0)) -
                                     (cg_face_fine(i,0)-som_gros(som2,0))*(cg_face_fine(i,2)-som_gros(som0,2));
                                p2 = p2 * prod;
                              }
                            else
                              {
                                p0 = -1;
                                p1 = -1;
                                p2 = -1;
                              }
                          }
                        else
                          //si les x sont les memes :
                          if((fabs(som_gros(som0,0) - som_gros(som1,0))<PRECISION) && (fabs(som_gros(som1,0) - som_gros(som2,0))<PRECISION))
                            {
                              //Cerr<<" x sont les memes"<<finl;
                              //si le centre de gravite est dans le plan forme par les 3 sommets
                              if(fabs(cg_face_fine(i,0) - som_gros(som0,0))<PRECISION)
                                {
                                  //Cerr<<" x sont les memes et centre de gravite fin est dans le plan forme des 3 sommets"<<finl;
                                  // Calcul de prod = 01 vectoriel 02 selon x
                                  // prod > 0 : sens trigo
                                  // prod < 0 : sens anti trigo
                                  prod = (som_gros(som1,1)-som_gros(som0,1))*(som_gros(som2,2)-som_gros(som0,2)) -
                                         (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,1)-som_gros(som0,1));
                                  if (prod >= 0)
                                    prod = 1;
                                  else
                                    prod = -1;
                                  // Calcul de p0 = 0M vectoriel 1M selon x
                                  p0 = (cg_face_fine(i,1)-som_gros(som0,1))*(cg_face_fine(i,2)-som_gros(som1,2)) -
                                       (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,1)-som_gros(som1,1));
                                  p0 = p0 * prod;
                                  // Calcul de p1 = 1M vectoriel 2M selon x
                                  p1 = (cg_face_fine(i,1)-som_gros(som1,1))*(cg_face_fine(i,2)-som_gros(som2,2)) -
                                       (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,1)-som_gros(som2,1));
                                  p1 = p1 * prod;
                                  // Calcul de p2 = 2M vectoriel 0M selon x
                                  p2 = (cg_face_fine(i,1)-som_gros(som2,1))*(cg_face_fine(i,2)-som_gros(som0,2)) -
                                       (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,1)-som_gros(som0,1));
                                  p2 = p2 * prod;
                                }
                              else
                                {
                                  p0 = -1;
                                  p1 = -1;
                                  p2 = -1;
                                }
                            }

                          else
                            {
                              p0 = -1;
                              p1 = -1;
                              p2 = -1;
                            }



                      if (p0 >= 0 && p1 >= 0 && p2 >= 0)
                        {
                          connect_faceF_faceG(i)=j;
                          if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                          else connect_faceG_faceF(j) = i;

                          trouve = 1;
                          break;
                        }


                      else if ( type_elem->que_suis_je() != "Tetraedre" )
                        {
                          som0 = faces_som(j,3);
                          som1 = faces_som(j,2);
                          som2 = faces_som(j,1);
                          //Cerr<<"est-il dans le deuxieme triangle?"<<finl;
                          // On regarde tout d'abord si le point cherche n'est pas un des
                          // sommets du triangle
                          if( ((fabs(som_gros(som0,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som0,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som0,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som1,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som1,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som1,2)-cg_face_fine(i,2))<PRECISION)) || ((fabs(som_gros(som2,0)-cg_face_fine(i,0))<PRECISION) && (fabs(som_gros(som2,1)-cg_face_fine(i,1))<PRECISION) && (fabs(som_gros(som2,2)-cg_face_fine(i,2))<PRECISION)) )
                            {
                              connect_faceF_faceG(i)=j;
                              if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                              else connect_faceG_faceF(j) = i;
                              trouve = 1;
                              break;
                            }
                          if (trouve == 0)
                            {

                              //Attention
                              // Il faut donc determiner le sens (trigo ou anti trigo)
                              //pour la numerotation :


                              //si les z sont les memes :
                              if((fabs(som_gros(som0,2) - som_gros(som1,2))<PRECISION) && (fabs(som_gros(som1,2) - som_gros(som2,2))<PRECISION))
                                {
                                  //si le centre de gravite est dans le plan defini par les 3 sommets
                                  if(fabs(cg_face_fine(i,2) - som_gros(som0,2))<PRECISION)
                                    {
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
                                    }
                                  else
                                    {
                                      p0 = -1;
                                      p1 = -1;
                                      p2 = -1;
                                    }
                                }
                              else
                                //si les y sont les memes :
                                if((fabs(som_gros(som0,1) - som_gros(som1,1))<PRECISION) && (fabs(som_gros(som1,1) - som_gros(som2,1))<PRECISION))
                                  {
                                    //si le centre de gravite est dans le plan defini par les 3 sommets
                                    if(fabs(cg_face_fine(i,1) - som_gros(som0,1))<PRECISION)
                                      {
                                        // Calcul de prod = 01 vectoriel 02 selon y
                                        // prod > 0 : sens trigo
                                        // prod < 0 : sens anti trigo
                                        prod = (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,0)-som_gros(som0,0)) - (som_gros(som1,0)-som_gros(som0,0))*(som_gros(som2,2)-som_gros(som0,2));
                                        if (prod >= 0)
                                          prod = 1;
                                        else
                                          prod = -1;
                                        // Calcul de p0 = 0M vectoriel 1M selon y
                                        p0 = (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,0)-som_gros(som1,0)) -
                                             (cg_face_fine(i,0)-som_gros(som0,0))*(cg_face_fine(i,2)-som_gros(som1,2));
                                        p0 = p0 * prod;
                                        // Calcul de p1 = 1M vectoriel 2M selon y
                                        p1 = (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,0)-som_gros(som2,0)) -
                                             (cg_face_fine(i,0)-som_gros(som1,0))*(cg_face_fine(i,2)-som_gros(som2,2));
                                        p1 = p1 * prod;
                                        // Calcul de p2 = 2M ve1ctoriel 0M selon y
                                        p2 = (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,0)-som_gros(som0,0)) -
                                             (cg_face_fine(i,0)-som_gros(som2,0))*(cg_face_fine(i,2)-som_gros(som0,2));
                                        p2 = p2 * prod;
                                      }
                                    else
                                      {
                                        p0 = -1;
                                        p1 = -1;
                                        p2 = -1;
                                      }
                                  }
                                else
                                  //si les x sont les memes :
                                  if((fabs(som_gros(som0,0) - som_gros(som1,0))<PRECISION) && (fabs(som_gros(som1,0) - som_gros(som2,0))<PRECISION))
                                    {
                                      //si le centre de gravite est dans le plan defini par les 3 sommets
                                      if(fabs(cg_face_fine(i,0) - som_gros(som0,0))<PRECISION)
                                        {
                                          // Calcul de prod = 01 vectoriel 02 selon x
                                          // prod > 0 : sens trigo
                                          // prod < 0 : sens anti trigo
                                          prod = (som_gros(som1,1)-som_gros(som0,1))*(som_gros(som2,2)-som_gros(som0,2)) -
                                                 (som_gros(som1,2)-som_gros(som0,2))*(som_gros(som2,1)-som_gros(som0,1));
                                          if (prod >= 0)
                                            prod = 1;
                                          else
                                            prod = -1;
                                          // Calcul de p0 = 0M vectoriel 1M selon x
                                          p0 = (cg_face_fine(i,1)-som_gros(som0,1))*(cg_face_fine(i,2)-som_gros(som1,2)) -
                                               (cg_face_fine(i,2)-som_gros(som0,2))*(cg_face_fine(i,1)-som_gros(som1,1));
                                          p0 = p0 * prod;
                                          // Calcul de p1 = 1M vectoriel 2M selon x
                                          p1 = (cg_face_fine(i,1)-som_gros(som1,1))*(cg_face_fine(i,2)-som_gros(som2,2)) -
                                               (cg_face_fine(i,2)-som_gros(som1,2))*(cg_face_fine(i,1)-som_gros(som2,1));
                                          p1 = p1 * prod;
                                          // Calcul de p2 = 2M vectoriel 0M selon x
                                          p2 = (cg_face_fine(i,1)-som_gros(som2,1))*(cg_face_fine(i,2)-som_gros(som0,2)) -
                                               (cg_face_fine(i,2)-som_gros(som2,2))*(cg_face_fine(i,1)-som_gros(som0,1));
                                          p2 = p2 * prod;
                                        }
                                      else
                                        {
                                          p0 = -1;
                                          p1 = -1;
                                          p2 = -1;
                                        }
                                    }

                                  else
                                    {
                                      p0 = -1;
                                      p1 = -1;
                                      p2 = -1;
                                    }



                              if (p0 >= 0 && p1 >= 0 && p2 >= 0)
                                {
                                  connect_faceF_faceG(i)=j;
                                  if (connect_faceG_faceF(j)>0) connect_faceG_faceF(j) = -1;
                                  else connect_faceG_faceF(j) = i;

                                  trouve = 1;
                                  break;
                                }
                            }
                        }
                    }
                }
            }

        }
    }
}




// Description:
//    Calcul des connectivites faces/faces
// Precondition:
// Parametre:Zone& zoneG, Zone& zoneF
//    Signification: zone discretisee grossiere et fine
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour:
//    Signification:
//    Contraintes:
// Exception: da
// Effets de bord:
// Postcondition: connect_faceF_faceG et connect_faceG_faceF sont remplis
void Connectivites_DWF::calculer_connectivites(Zone_VF& zonef,
                                               Zone_VF& zoneg,
                                               Domaine& domg)
{
  calculer_connectivites_face_face(zonef, zoneg, domg);
  calculer_connectivites_elem_elem(zonef, zoneg);
}

