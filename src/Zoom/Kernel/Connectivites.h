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
// File:        Connectivites.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Connectivites_included
#define Connectivites_included


#include <Connectivites_base.h>
#include <TRUST_Deriv.h>
#include <TRUSTTabs_forward.h>
#include <Domaine_forward.h>

class Domaine_VF;

/*! @brief class Connectivites
 *
 */

//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Connectivites
//
/////////////////////////////////////////////////////////////////////////////////

class Connectivites : public DERIV(Connectivites_base)
{
  Declare_instanciable(Connectivites);

public:
  inline void calculer_connectivites(Domaine_VF& domainef, Domaine_VF& domaineg, Domaine& domg);
  inline IntVect& connectivites_faceF_faceG();
  inline IntVect& connectivites_elemF_elemG();
  inline const IntVect& connectivites_faceF_faceG() const;
  inline const IntVect& connectivites_elemF_elemG() const;
};

// Fonctions inline

inline void Connectivites::calculer_connectivites(Domaine_VF& domainef, Domaine_VF& domaineg,
                                                  Domaine& domg)
{
  valeur().calculer_connectivites(domainef,domaineg, domg);
}

inline IntVect& Connectivites::connectivites_faceF_faceG()
{
  return valeur().connectivites_faceF_faceG();
}

inline const IntVect& Connectivites::connectivites_faceF_faceG() const
{
  return valeur().connectivites_faceF_faceG();
}

inline IntVect& Connectivites::connectivites_elemF_elemG()
{
  return valeur().connectivites_elemF_elemG();
}

inline const IntVect& Connectivites::connectivites_elemF_elemG() const
{
  return valeur().connectivites_elemF_elemG();
}

#endif
