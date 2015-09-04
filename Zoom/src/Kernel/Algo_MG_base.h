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
// File:        Algo_MG_base.h
// Directory:   $TRUST_ROOT/src/Zoom/Kernel
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


// Faire 2 classes : algo_composite et algo_MG


#ifndef Algo_MG_base_included
#define Algo_MG_base_included

#include <Ref_Pb_MG.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     classe Algo_MG_base
//     Cette classe abstraite definit un algorithme de resolution
//     associe a un couplage (de type Pb_MG).
//     Il est suppose que la resolution se fait par avance en temps.
//
// .SECTION voir aussi
//////////////////////////////////////////////////////////////////////////////

class Algo_MG_base : public Objet_U
{
  Declare_base(Algo_MG_base);

public:

  virtual int initialiser()=0;
  virtual double computeTimeStep(bool& stop) const=0;
  virtual bool solveTimeStep()=0;


  // Methodes d'acces aux membres :

  virtual int associer_pb(const Pb_MG&);
  Pb_MG& pb()
  {
    return mon_probleme.valeur();
  }
  const Pb_MG& pb() const
  {
    return mon_probleme.valeur();
  }

protected :

  REF(Pb_MG) mon_probleme;

};

#endif
