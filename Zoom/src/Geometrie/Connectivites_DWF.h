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
// File:        Connectivites_DWF.h
// Directory:   $TRUST_ROOT/src/Zoom/Geometrie
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Connectivites_DWF_included
#define Connectivites_DWF_included


#include <Zone_VF.h>
#include <Connectivites_base.h>

#include <Domaine.h>
//
// .DESCRIPTION class Connectivites_DWF
//
// .SECTION voir aussi


//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Connectivites_DWF
// Cette classe etend la classe de Connectivites_base. Elle ne calcule que les connectivites
// faces/faces et dans les deux sens Grossier->Fin et Fin->Grossier en supposant que
// 1 face grossiere = 1 face fine sinon stop.
/////////////////////////////////////////////////////////////////////////////////

class Connectivites_DWF : public Connectivites_base
{
  Declare_instanciable(Connectivites_DWF);

public:
  virtual void calculer_connectivites(Zone_VF& zonef, Zone_VF& zoneg,
                                      Domaine& domg);
  IntVect& connectivites_faceG_faceF() ;


private:
  IntVect connect_faceG_faceF;
  void calculer_connectivites_face_face(Zone_VF& zonef, Zone_VF& zoneg,
                                        Domaine& domg);

};

// Fonctions inline



inline IntVect& Connectivites_DWF::connectivites_faceG_faceF()
{
  return connect_faceG_faceF;
}


#endif
