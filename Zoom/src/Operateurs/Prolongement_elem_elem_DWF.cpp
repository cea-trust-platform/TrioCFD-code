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
// File:        Prolongement_elem_elem_DWF.cpp
// Directory:   $TRUST_ROOT/src/Zoom/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_elem_elem_DWF.h>

Implemente_instanciable(Prolongement_elem_elem_DWF,"Prolongement_elem_elem_DWF",Prolongement_base);

//// printOn
//

Sortie& Prolongement_elem_elem_DWF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_elem_elem_DWF::readOn(Entree& s )
{
  return s ;
}




// Description:
//    Prolongement de l'inconnue grossiere sur la frontiere fine
//    pour inconnue TEMPERATURE ou PRESSION
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
void Prolongement_elem_elem_DWF::prolonger(Zone_VF& zone_VFG, Zone_VF& zone_VFF,
                                           const Frontiere& frontF,
                                           IntVect& connect,
                                           const DoubleTab& valG, DoubleTab& tab,
                                           int nb_compo)
{
  tab = 0.;

  const int prem_face_bord_fin  =  frontF.num_premiere_face();
  const int nb_faces_front_fine = frontF.nb_faces();
  const IntTab& face_voisinsF = zone_VFF.face_voisins();

  if(sub_type(Zone_VDF, zone_VFG) && sub_type(Zone_VDF, zone_VFF))
    {
      for (int ifaceF=0; ifaceF<nb_faces_front_fine; ifaceF++)
        {
          const int num_face=prem_face_bord_fin+ifaceF;        // face fine courante
          int elemF = face_voisinsF(num_face,0);
          if (elemF==-1) elemF = face_voisinsF(num_face,1);        // element fin correspondant
          const int elemG = connect(elemF);                        // element grossier correspondant
          tab(ifaceF, 0) = valG(elemG);
        }
    }
}



//NE FAIT RIEN
void Prolongement_elem_elem_DWF::calculer(Zone_VF& zonef,
                                          Zone_VF& zoneg,
                                          IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}
