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
// File:        Convection_Diffusion_Phase_field.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Convection_Diffusion_Phase_field_included
#define Convection_Diffusion_Phase_field_included

#include <Convection_Diffusion_Concentration.h>
#include <Operateur_Grad.h>
#include <Debog.h>
#include <Champ_Fonc.h>
#include <Champ_Don.h>
#include <Ref_Champ_Don.h>
#include <Constituant.h>

/*! @brief classe Convection_Diffusion_Phase_field Cas particulier de Convection_Diffusion_Concentration
 *
 *      pour un ou plusieurs constituants.
 *      Dans le cas de plusieurs constituants les champs
 *      concentration et diffusivite sont vectoriels.
 *
 * @sa Convection_Diffusion_Concentration
 */
class Convection_Diffusion_Phase_field : public Convection_Diffusion_Concentration
{
  Declare_instanciable_sans_constructeur(Convection_Diffusion_Phase_field);
public:

  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  Convection_Diffusion_Phase_field();
  void discretiser() override;
  Champ_Fonc& set_mutilde_() { return ch_mutilde; }
  DoubleTab& set_mutilde() { return mutilde; }
  DoubleTab& set_mutilde_demi() { return mutilde_demi; }
  DoubleTab& set_c_demi() { return c_demi; }
  DoubleTab& set_div_alpha_rho_gradC() { return div_alpha_rho_gradC; }
  DoubleTab& set_div_alpha_gradC() { return div_alpha_gradC; }
  DoubleTab& set_alpha_gradC_carre() { return alpha_gradC_carre; }
  DoubleTab& set_pression_thermo() { return pression_thermo; }
  const Champ_Fonc& get_mutilde_() const { return ch_mutilde; }
  const DoubleTab& get_mutilde() const { return mutilde; }
  const DoubleTab& get_mutilde_demi() const { return mutilde_demi; }
  const DoubleTab& get_c_demi() const { return c_demi; }
  const DoubleTab& get_div_alpha_rho_gradC() const { return div_alpha_rho_gradC; }
  const DoubleTab& get_div_alpha_gradC() const { return div_alpha_gradC; }
  const DoubleTab& get_alpha_gradC_carre() const { return alpha_gradC_carre; }
  const DoubleTab& get_pression_thermo() const { return pression_thermo; }
  int preparer_calcul() override;
  Operateur_Grad& operateur_gradient();
  const Operateur_Grad& operateur_gradient() const;
  void completer() override;

  inline int& get_mutype_() { return mutype_; }
  inline int get_mutype_() const { return mutype_; }
protected:
  Champ_Fonc ch_mutilde;
  Operateur_Grad gradient;
  DoubleTab div_alpha_gradC, div_alpha_rho_gradC, alpha_gradC_carre;
  DoubleTab pression_thermo, mutilde, c_demi, mutilde_demi;

  int mutype_;
};

#endif
