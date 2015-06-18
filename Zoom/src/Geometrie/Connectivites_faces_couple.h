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
// File:        Connectivites_faces_couple.h
// Directory:   $TRUST_ROOT/src/Zoom/Geometrie
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Connectivites_faces_couple_included
#define Connectivites_faces_couple_included


#include <Zone_VF.h>
#include <Connectivites_base.h>

#include <Domaine.h>
//
// .DESCRIPTION class Connectivites_faces_couple
//
// .SECTION voir aussi


//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Connectivites_faces_couple
// Cette classe etend la classe de Connectivites_base. Elle fournit en plus un tableau
// indiquant le nombre de faces fines contenues dans chaque face grossiere
/////////////////////////////////////////////////////////////////////////////////

class Connectivites_faces_couple : public Connectivites_base
{
  Declare_instanciable(Connectivites_faces_couple);

public:
  virtual void calculer_connectivites(Zone_VF& zonef, Zone_VF& zoneg,
                                      Domaine& domg);
  inline IntVect& nb_facesF();
  inline const IntVect& nb_facesF() const;


private:
  IntVect nb_facesF_; // tableau d'entier de taille NFACE__FRONT_GROS

  void calculer_nb_facesF(Zone_VF& zoneg);

};

// Fonctions inline



inline IntVect& Connectivites_faces_couple::nb_facesF()
{
  return nb_facesF_;
}


inline const IntVect& Connectivites_faces_couple::nb_facesF() const
{
  return nb_facesF_;
}


#endif
