/****************************************************************************
* Copyright (c) 2019, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Navier_Stokes_std_sensibility.cpp
// Directory : $BALTIK_COUPLAGE_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////

#include <Navier_Stokes_std_sensibility.h>
#include <Probleme_base.h>
#include <Param.h>
#include <Debog.h>
#include <LecFicDiffuse.h>
#include <communications.h>

Implemente_instanciable( Navier_Stokes_std_sensibility, "Navier_Stokes_standard_sensibility", Navier_Stokes_std) ;

Sortie& Navier_Stokes_std_sensibility::printOn( Sortie& os ) const
{
  return Navier_Stokes_std::printOn( os );
}

Entree& Navier_Stokes_std_sensibility::readOn( Entree& is )
{
  Navier_Stokes_std::readOn(is);

  return is;
}

void Navier_Stokes_std_sensibility::set_param(Param& param)
{

  Navier_Stokes_std::set_param(param);
  param.ajouter_non_std("state",(this),Param::REQUIRED);

}


int Navier_Stokes_std_sensibility::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="state")
    {
      Cerr << "Reading and typing of the state : " << finl;
      return 1;
    }
  else
    return Navier_Stokes_std::lire_motcle_non_standard(mot, is);
}
