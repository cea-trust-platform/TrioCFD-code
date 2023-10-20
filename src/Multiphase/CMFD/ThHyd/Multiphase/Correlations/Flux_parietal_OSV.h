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
// File:        Flux_parietal_OSV.h
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Flux_parietal_OSV_included
#define Flux_parietal_OSV_included
#include <TRUSTTab.h>
#include <Flux_parietal_base.h>
#include <Correlation.h>
#include <Param.h>

/*! @brief Flux_parietal_OSV class that implements an onset of significant void-based heat flux partition
 *
 * It requires the use of a single-phase correlation, usually the adaptive thermal wall law is used
 * It also requires a multi-phase total heat flux correlation, likes those of Thom, Jens & Lottes or Frost & Dzakowic.
 * Jens & Lottes is the default for now but others can be adjusted in the future.
 * See https://hal.science/hal-04524455 for details
 *
 */
class Flux_parietal_OSV : public Flux_parietal_base
{
  Declare_instanciable(Flux_parietal_OSV);
public:
  virtual void qp(const input_t& input, output_t& output) const override;

  virtual void completer() override;

  int calculates_bubble_nucleation_diameter() const override {return 1;} ;
  int needs_saturation() const override {return 1;} ;

  virtual int T_at_wall() const override {return 0;};

protected :
  Correlation correlation_monophasique_;
  REF(Correlation) correlation_loi_paroi_;
  double Dh_ = -1. ;
  int Water_  = 0;
  int R12_   = 0;
  int R22_   = 0;
  int R134a_ = 0;
  double fac_pression_ = 1.;
  double fac_puissance_ = 1.;
  double beta_ = -7.;
};

#endif
