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
// File:        UtilitairesPrepro.h
// Directory:   $TRUST_ROOT/src/UtilitairesAssemblages
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////


#ifndef UtilitairesPrepro_included
#define UtilitairesPrepro_included

#include <Double.h>

//
// .DESCRIPTION class UtilitairesPrepro
//

//
//

class UtilitairesPrepro
{
public:

  static inline int hexa(const double& x,const double& y,const double& ep);
  static inline int coord_assemblage_M(const double& x,const double& y, int M0, const double& pas_y);
  static inline int coord_assemblage_N(const double& x,const double& y, int M0, const double& pas_y);

};


int UtilitairesPrepro::hexa(const double& x,const double& y,const double& ep)

{
  // convention des faces = CADET

  if ( !sup_ou_egal(y,-ep/2.) )         return 6;
  if ( !sup_ou_egal(y,-ep-x*sqrt(3.)) ) return 5;
  if ( !sup_ou_egal(y,-ep+x*sqrt(3.)) ) return 1;
  if ( !inf_ou_egal(y, ep/2.) )                return 3;
  if ( !inf_ou_egal(y, ep-x*sqrt(3.)) ) return 2;
  if ( !inf_ou_egal(y, ep+x*sqrt(3.)) ) return 4;

  return 0;  // IN
}

#ifndef _COMPILE_AVEC_GCC_3
#ifndef _COMPILE_AVEC_INTEL
// PL: round n'est pas toujours portable, utiliser plutot floor
inline double round(double value)
{
  if (value < 0)
    return -(floor(-value + 0.5));
  else
    return   floor( value + 0.5);
}
#endif
#endif

int UtilitairesPrepro::coord_assemblage_M(const double& x,const double& y, int M0, const double& pas_y)
{
  return M0 + (int) round((sqrt(3.)/3.*x-y)/pas_y) ;
}



int UtilitairesPrepro::coord_assemblage_N(const double& x,const double& y, int M0, const double& pas_y)
{
  return M0 + (int) round((sqrt(3.)/3.*x+y)/pas_y) ;
}


#endif
