/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Navier_Stokes_std_ALE.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Navier_Stokes_std_ALE_included
#define Navier_Stokes_std_ALE_included


#include <Navier_Stokes_std.h>

class Navier_Stokes_std_ALE: public Navier_Stokes_std
{

  Declare_instanciable(Navier_Stokes_std_ALE);

public :
  void renewing_jacobians( DoubleTab& derivee ) override;
  void div_ale_derivative( DoubleTrav& deriveeALE, double timestep, DoubleTab& derivee, DoubleTrav& secmemP ) override;
  void update_pressure_matrix( void ) override;
  void discretiser() override;
  void mettre_a_jour(double) override;
  virtual int sauvegarder(Sortie&) const override;
  virtual int reprendre(Entree&) override;

protected :
  Champ_Inc ALEMeshVelocity_;
  Champ_Inc ALEMeshTotalDisplacement_;
  Champ_Inc ALEMeshStructuralPressure_;
  Champ_Inc ALEMeshStructuralVonMises_;
  Champ_Inc ALEMeshStructuralForces_;
};

#endif
