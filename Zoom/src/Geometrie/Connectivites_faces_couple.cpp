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
// File:        Connectivites_faces_couple.cpp
// Directory:   $TRUST_ROOT/src/Zoom/Geometrie
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Connectivites_faces_couple.h>

Implemente_instanciable(Connectivites_faces_couple,"Connectivites_faces_couple",Connectivites_base);

//// printOn
//

Sortie& Connectivites_faces_couple::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Connectivites_faces_couple::readOn(Entree& s )
{
  return s ;
}


void Connectivites_faces_couple::calculer_nb_facesF(Zone_VF& zone_vfG)
{
  Cerr<<"debut de Connectivites_faces_couple::calculer_nb_facesF"<<finl;
  //Cerr<<"CALCUL DU NOMBRE DE FACES FINES PAR FACE GROSSIERE"<<finl;
  int nb_facesG = zone_vfG.nb_faces();
  nb_facesF_.resize(nb_facesG);
  int nb_faces_fines = connect_faceF_faceG.size_array();
  int face, num_faceG;

  //CALCUL DU NOMBRE DE FACES FINES PAR FACE GROSSIERE
  nb_facesF_ = 0;
  for(face=0; face<nb_faces_fines; face++)
    {
      //on recupere la face grossiere correspondante
      num_faceG = connect_faceF_faceG(face);
      if(num_faceG != -1)
        {
          nb_facesF_(num_faceG) = nb_facesF_(num_faceG)+1;
        }
    }
}



// Description:
//    Calcul des connectivites entre facesF et facesG
//    Calcul du nombre de faces fines par face grossiere
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
// Postcondition: connect_elemF_elemG et connect_face_face sont remplis
void Connectivites_faces_couple::calculer_connectivites(Zone_VF& zonef,
                                                        Zone_VF& zoneg,
                                                        Domaine& domg)
{
  Cerr<<"debut de Connectivites_faces_couple::calculer_connectivites"<<finl;
  Connectivites_base::calculer_connectivites_face_face(zonef, zoneg, domg);
  calculer_nb_facesF(zoneg);
}
