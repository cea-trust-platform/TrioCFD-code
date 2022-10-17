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
// File:        Flux_parietal_adaptatif.h
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Flux_parietal_adaptatif_included
#define Flux_parietal_adaptatif_included
#include <TRUSTTab.h>
#include <Flux_parietal_base.h>
#include <Ref_Correlation.h>

/*! @brief classe Flux_parietal_adaptatif classe qui implemente une correlation de flux parietal monophasique
 *
 *       pour un ecoulement turbulent avec une loi de paroi adaptative
 *       (i.e. qui peut gerer la zone visqueuse comme la zone log en proche paroi)
 *
 *
 *
 */
class Flux_parietal_adaptatif : public Flux_parietal_base
{
  Declare_instanciable(Flux_parietal_adaptatif);
public:
  virtual void qp(const input_t& input, output_t& output) const override;

  virtual void completer() override;

protected :
  double calc_theta_plus(double y, double u_tau, double mu, double lambda, double rho, double Cp, double Diam_hyd_) const;

  REF(Correlation) correlation_loi_paroi_;
};

#endif
