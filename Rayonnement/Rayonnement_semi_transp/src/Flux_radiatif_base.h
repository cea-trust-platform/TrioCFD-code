/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Flux_radiatif_base.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Flux_radiatif_base_included
#define Flux_radiatif_base_included

#include <Neumann_paroi.h>
class Equation_base;

/*! @brief Decrire ici la classe Flux_radiatif_base
 *
 *  .FINHTML
 *  .FINEPS
 *
 *
 */
class Flux_radiatif_base : public Neumann_paroi
{
  Declare_base(Flux_radiatif_base);

public :
  int compatible_avec_eqn(const Equation_base&) const override { return 1; }
  inline Champ_front& emissivite();
  inline const Champ_front& emissivite() const;
  inline double& A();
  inline const double& A() const;

  virtual void calculer_flux_radiatif(const Equation_base& eq_temp)=0;
  inline Champ_front& flux_radiatif();
  inline const Champ_front& flux_radiatif() const;
  void completer() override;

  double flux_impose(int i) const override;
  double flux_impose(int i,int j) const override;

protected :
  double A_;
  Champ_front emissivite_;
  Champ_front flux_radiatif_;
};



inline Champ_front& Flux_radiatif_base::flux_radiatif()
{
  return flux_radiatif_;
}

inline const Champ_front& Flux_radiatif_base::flux_radiatif() const
{
  return flux_radiatif_;
}

inline double& Flux_radiatif_base::A()
{
  return A_;
}

inline const double& Flux_radiatif_base::A() const
{
  return A_;
}


inline Champ_front& Flux_radiatif_base::emissivite()
{
  return emissivite_;
}

inline const Champ_front& Flux_radiatif_base::emissivite() const
{
  return emissivite_;
}

#endif
