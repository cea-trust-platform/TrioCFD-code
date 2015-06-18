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
// File:        Algo_Couple_1.h
// Directory:   $TRUST_ROOT/src/Zoom/Algos
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Algo_Couple_1_included
#define Algo_Couple_1_included



#include <Algo_MG_base.h>
class Pb_MG;
//
// .DESCRIPTION class Algo_Couple_1
//
// .SECTION voir aussi


//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Algo_Couple_1
//
/////////////////////////////////////////////////////////////////////////////////

class Algo_Couple_1 : public Algo_MG_base
{
  Declare_instanciable(Algo_Couple_1);

public:

  virtual int initialiser();

  // Les methodes suivantes sont specifiques a une resolution par
  // avance en temps.

  virtual double computeTimeStep(bool& stop) const;
  virtual bool solveTimeStep();

private:
  int dt_unif;
  double dt_;
};

#endif
