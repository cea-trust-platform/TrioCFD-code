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
// File:        Prolongement_face_face_DWF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_face_face_DWF.h>

Implemente_instanciable(Prolongement_face_face_DWF,"Prolongement_face_face_DWF",Prolongement_base);

//// printOn
//

Sortie& Prolongement_face_face_DWF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_face_face_DWF::readOn(Entree& s )
{
  return s ;
}




// Description:
//    Prolongement de l'inconnue grossiere sur la frontiere fine
//    pour inconnue VITESSE
// Dans ce prolongement, la frontiere fine ne repose pas sur une frontiere grossiere
// !  chaque face fine de la frontiere fine est situee a l'interieur  d'un element grossier.
// Les connectivites a passer sont du type elemF->elemG
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
void Prolongement_face_face_DWF::prolonger(Zone_VF& zone_VFG, Zone_VF& zone_VFF,
                                           const Frontiere& frontF,
                                           IntVect& connect,
                                           const DoubleTab& valG, DoubleTab& tab,
                                           int nb_compo)
{
  tab = 0.;
  int comp;

  const int prem_face_bord_fin  =  frontF.num_premiere_face();
  const int nb_faces_front_fine = frontF.nb_faces();
  const IntTab& face_voisinsF = zone_VFF.face_voisins();


  const IntTab& elem_facesG =  zone_VFG.elem_faces();

  if(sub_type(Zone_VDF, zone_VFG) && sub_type(Zone_VDF, zone_VFF))
    {
      const Zone_VDF& zone_vdfF = ref_cast(Zone_VDF, zone_VFF);
      const IntVect& oriF = zone_vdfF.orientation();
      //const Zone_VDF& zone_vdfG = ref_cast(Zone_VDF, zone_VFG);
      // const IntVect& oriG = zone_vdfG.orientation();
      for (int ifaceF=0; ifaceF<nb_faces_front_fine; ifaceF++)
        {
          const int num_face=prem_face_bord_fin+ifaceF; // face fine courante
          int elemF = face_voisinsF(num_face,0);
          if (elemF==-1) elemF = face_voisinsF(num_face,1); // element fin correspondant
          const int elemG = connect(elemF);                 // element grossier correspondant

          if (nb_compo > 1)// cas vectoriel : on moyenne
            {
              const int ori_fine = oriF(ifaceF);

              if(ori_fine == 0)
                {
                  int fg1 = elem_facesG(elemG,1);
                  int fg2 = elem_facesG(elemG,1+dimension);

                  //composante en Y
                  tab(ifaceF, 1) = 0.5*(valG(fg1)+valG(fg2));

                  //composante en X

                  fg1 = elem_facesG(elemG,0);
                  fg2 = elem_facesG(elemG,dimension);
                  tab(ifaceF, 0) =  0.5*(valG(fg1)+valG(fg2));

                  if (dimension == 3)
                    {
                      fg1 = elem_facesG(elemG,2);
                      fg2 = elem_facesG(elemG,2+dimension);

                      //composante en Z
                      tab(ifaceF, 2) = 0.5*(valG(fg1)+valG(fg2));
                    }
                }
              else if (ori_fine == 1)
                {
                  int fg1 = elem_facesG(elemG,0);
                  int fg2 = elem_facesG(elemG,dimension);

                  //composante en X
                  tab(ifaceF, 0) = 0.5*(valG(fg1)+valG(fg2));

                  //composante en Y

                  fg1 = elem_facesG(elemG,1);
                  fg2 = elem_facesG(elemG,1+dimension);
                  tab(ifaceF, 1) =  0.5*(valG(fg1)+valG(fg2));

                  if (dimension == 3)
                    {
                      fg1 = elem_facesG(elemG,2);
                      fg2 = elem_facesG(elemG,2+dimension);

                      //composante en Z
                      tab(ifaceF, 2) = 0.5*(valG(fg1)+valG(fg2));
                    }
                }
              else if (ori_fine == 2)
                {
                  int fg1 = elem_facesG(elemG,0);
                  int fg2 = elem_facesG(elemG,dimension);

                  //composante en X
                  tab(ifaceF, 0) = 0.5*(valG(fg1)+valG(fg2));

                  fg1 = elem_facesG(elemG,1);
                  fg2 = elem_facesG(elemG,1+dimension);

                  //composante en Y
                  tab(ifaceF, 1) = 0.5*(valG(fg1)+valG(fg2));

                  //composante en Z

                  fg1 = elem_facesG(elemG,2);
                  fg2 = elem_facesG(elemG,2+dimension);
                  tab(ifaceF, 2) =  0.5*(valG(fg1)+valG(fg2));
                }

            }
        }
    }
  //if(sub_type(Zone_VEF, zone_VFG))
  {
    for (comp=0; comp<nb_compo; comp++)
      {
        //tab(ifaceF, comp) = valG(num_faceG, comp);
      }
  }

  //Cerr<<"fin de Prolongement_face_face_DWF::prolonger"<<finl;
}




//NE FAIT RIEN
void Prolongement_face_face_DWF::calculer(Zone_VF& zonef,
                                          Zone_VF& zoneg,
                                          IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}
