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
// File:        Restriction_base.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Restriction_base_included
#define Restriction_base_included



#include <Objet_U.h>
#include <TRUSTTabs_forward.h>
class Zone_VF;
#include <TRUSTTabs_forward.h>

/*! @brief class Restriction_base
 *
 */

//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Restriction_base
//
/////////////////////////////////////////////////////////////////////////////////

class Restriction_base : public Objet_U
{
  Declare_base(Restriction_base);

public:

  virtual void restreindre(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF, const IntVect& connect,
                           DoubleTab& incoG,
                           const DoubleTab& incoF,int nbcomp) =0;
  virtual void calculer(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF, const IntVect& connect) =0;
  virtual void restreindre(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF,
                           const  IntVect& connect,
                           DoubleTab& incoG,
                           const DoubleTab& incoF, int nb_comp,
                           int num_prem_face_frontG) =0;

};

#endif
