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
// File:        Connectivites_base.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Connectivites_base_included
#define Connectivites_base_included


#include <TRUSTVect.h>
class Domaine_VF;
class Domaine;
/*! @brief class Connectivites_base
 *
 */

//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Connectivites_base
//
/////////////////////////////////////////////////////////////////////////////////

class Connectivites_base : public Objet_U
{
  Declare_instanciable(Connectivites_base);

public:
  virtual void calculer_connectivites(Domaine_VF& domainef, Domaine_VF& domaineg,
                                      Domaine& domg);
  inline IntVect& connectivites_faceF_faceG();
  inline IntVect& connectivites_elemF_elemG();
  inline const IntVect& connectivites_faceF_faceG()const ;
  inline const IntVect& connectivites_elemF_elemG()const ;


protected:
  IntVect connect_faceF_faceG;
  IntVect connect_elemF_elemG;
  void calculer_connectivites_face_face(Domaine_VF& domainef, Domaine_VF& domaineg,
                                        Domaine& domg);
  void calculer_connectivites_elem_elem(Domaine_VF& domainef, Domaine_VF& domaineg);



};

// Fonctions inline

inline IntVect& Connectivites_base::connectivites_faceF_faceG()
{
  return connect_faceF_faceG;
}



inline IntVect& Connectivites_base::connectivites_elemF_elemG()
{
  return connect_elemF_elemG;
}


inline const IntVect& Connectivites_base::connectivites_faceF_faceG() const
{
  return connect_faceF_faceG;
}



inline const IntVect& Connectivites_base::connectivites_elemF_elemG() const
{
  return connect_elemF_elemG;
}


#endif
