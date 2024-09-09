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
// File:        Dispersion_bulles_turbulente_Burns.h
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Dispersion_bulles_turbulente_Burns_included
#define Dispersion_bulles_turbulente_Burns_included

#include <Dispersion_bulles_base.h>
#include <Correlation_base.h>
#include <TRUST_Ref.h>

/*! @brief classe Dispersion_bulles_turbulente_Burns
 *     coefficients de dispersion selon le modele Burns et al 2004
 *
 */

class Dispersion_bulles_turbulente_Burns : public Dispersion_bulles_base
{
  Declare_instanciable(Dispersion_bulles_turbulente_Burns);
public:
  void coefficient(const input_t& input, output_t& output) const override;
  void completer() override;

protected:
  REF(Correlation) correlation_drag_;
  int n_l = -1; //phase liquide
  double Prt_ = .9 ; // Turbulent Prandtl number
  double minimum_ = -1.;
  double a_res_ = -1.;
  double g_ = 9.81;
  double C_lambda_ = 2.7;
  double gamma_ = 1.; // rapport d'aspect moyen des bulles ; mis sous forme de daouble pour le moment ; possibilite d'ajouter une fermeture a l'avenir
  double coefBIA_ = 0.;

};

#endif
