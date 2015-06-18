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
// File:        Prolongement_base.h
// Directory:   $TRUST_ROOT/src/Zoom/Kernel
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Prolongement_base_included
#define Prolongement_base_included


#include <Objet_U.h>
class Zone_VF;
class Frontiere;
class IntVect;
class DoubleTab;

//
// .DESCRIPTION class Prolongement_base
//
// .SECTION voir aussi


//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Prolongement_base
//
/////////////////////////////////////////////////////////////////////////////////

class Prolongement_base : public Objet_U
{
  Declare_base(Prolongement_base);

public:

  virtual void calculer(Zone_VF& , Zone_VF& , IntVect&) =0;
  virtual void prolonger(Zone_VF& zone_VFG, Zone_VF& zone_VFF,
                         const Frontiere& frontF,
                         IntVect& connect,
                         const DoubleTab& incoG,
                         DoubleTab& tab, int nb_comp) =0;


};

#endif
