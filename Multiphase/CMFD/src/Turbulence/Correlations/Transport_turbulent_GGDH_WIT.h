/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Transport_turbulent_GGDH_WIT.h
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_turbulent_GGDH_WIT_included
#define Transport_turbulent_GGDH_WIT_included
#include <TRUSTTab.h>
#include <Transport_turbulent_base.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Transport_turbulent_GGDH_WIT
//    Transport turbulent de type GGDH pour l'equation sur WIT:
//    < u'_i theta'> = - C_s * temps caract * <u'_i u'_j> * d_j theta
//////////////////////////////////////////////////////////////////////////////

class Transport_turbulent_GGDH_WIT : public Transport_turbulent_base
{
  Declare_instanciable(Transport_turbulent_GGDH_WIT);
public:
  virtual int dimension_min_nu() const override
  {
    return dimension * dimension; //anisotrope complet!
  }
  virtual int gradu_required() const override
  {
    return 1; /* on a besoin de grad u */
  };
  virtual void modifier_nu(const Convection_Diffusion_std& eq, const Viscosite_turbulente_base& visc_turb, DoubleTab& nu) const override;
private:
  double delta_ = 1.3 ; // rapport volume perturbe par la bulle / volume de la bulle
  double gamma_ = 1. ; //rapport d'aspect des bulles
  double limiteur_alpha_ = 1e-6;
  double C_s = 1.;
};

#endif
