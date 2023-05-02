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
// File:        Source_Dissipation_energie_cin_turb.h
// Directory:   $TRUST_ROOT/src/Turbulence/Sources
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Dissipation_energie_cin_turb_included
#define Source_Dissipation_energie_cin_turb_included

#include <Sources_Multiphase_base.h>

class Convection_Diffusion_std;
/*! @brief class Source_Dissipation_energie_cin_turb
 *
 *  Terme de dissipation -beta * (alpha * rho * k) * omega dans l'equation d'energie cinetique turbulente
 *
 *  la phase dont la turbulence est decrite avec le modele k-tau doit etre ecrite en premier dans le bloc phases { } du jeu de donnees
 *  Actuellement k est necessairement un scalaire.
 *  Si cela est amene a evolue pour permettre de la turbulence dans plusieurs phases, il faudra alors revoir cette classe en iterant sur les id_composites des phases turbulentes.
 *  en l'etat, si plusieurs phases sont turbulentes et sont decrites par le modele k-tau, alors elles doivent se suivre dans le bloc phases { } du jeu de donnees
 *
 */
class Source_Dissipation_energie_cin_turb : public Sources_Multiphase_base
{

  Declare_base(Source_Dissipation_energie_cin_turb);

public:
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;
  void completer() override;

protected:
  double beta_k = 0.09; // Wilcox
};

#endif
