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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <DoubleTrav.h>
#include <Zone_VDF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem,"Source_Transport_K_Eps_Bas_Reynolds_W_VDF_P0_VDF",Source_Transport_Bas_Reynolds_VDF_Elem_base);

Sortie& Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::readOn(Entree& is ) { return Source_Transport_Bas_Reynolds_VDF_Elem_base::readOn(is); }

void Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem::fill_resu_bas_reyn(const DoubleTrav& P, const DoubleTrav& D, const DoubleTrav& E, const DoubleTrav& F1, const DoubleTrav& F2,DoubleTab& resu) const
{
  const DoubleTab& K_eps_Bas_Re = eqn_keps_bas_re->inconnue().valeurs();
  const DoubleVect& volumes = la_zone_VDF->volumes(), &porosite_vol = la_zone_VDF->porosite_elem();

  for (int elem = 0; elem < la_zone_VDF->nb_elem(); elem++)
    {
      resu(elem,0) += (P(elem)-K_eps_Bas_Re(elem,1)-D(elem))*volumes(elem)*porosite_vol(elem);
      if (K_eps_Bas_Re(elem,0) >= 10.e-10)
        resu(elem,1) += (C1*F1(elem)*P(elem)- C2*F2(elem)*K_eps_Bas_Re(elem,1))*volumes(elem)*porosite_vol(elem)
                        *K_eps_Bas_Re(elem,1)/K_eps_Bas_Re(elem,0)+E(elem)*volumes(elem)*porosite_vol(elem);
    }
}
