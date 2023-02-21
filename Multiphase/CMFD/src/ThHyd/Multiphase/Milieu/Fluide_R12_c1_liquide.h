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

#ifndef Fluide_R12_c1_liquide_included
#define Fluide_R12_c1_liquide_included

#include <Fluide_reel_base.h>

class Fluide_R12_c1_liquide: public Fluide_reel_base
{
  Declare_instanciable(Fluide_R12_c1_liquide);

  /* bornes (p, T_l) de BOUNDRSR12.H */
  std::map<std::string, std::array<double, 2>> unknown_range() const override
  {
    return { { "pression" , { 1000., 4000000.}}, { "temperature", { -273., 150.}} };
  }

protected :
  void rho_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void dP_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void dT_rho_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void h_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void dP_h_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void dT_h_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void cp_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void beta_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void mu_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;
  void lambda_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int id = 0) const override;

private :
  int like_eos_ = 0;
};

#endif /* Fluide_eau_c3_liquide_included */
