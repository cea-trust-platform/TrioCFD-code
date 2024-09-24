/****************************************************************************
* Copyright (c) 2023, CEA
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
/*! @brief class Op_Diff_K_Omega_VEF_Face Cette classe represente l'operateur de diffusion turbulente
 *
 *    La discretisation est VEF
 *   Les methodes pour l'implicite sont codees.
 *
 */
#ifndef Op_Diff_K_Omega_VEF_base_included
#define Op_Diff_K_Omega_VEF_base_included

#define PRDT_K_DEFAUT 0.5
#define PRDT_OMEGA_DEFAUT 0.5

#include <Op_Diff_K_Omega_base.h>
#include <Op_VEF_Face.h>
#include <TRUST_Ref.h>





class Domaine_Cl_VEF;
class Domaine_VEF;
class Champ_P1NC;
class Champ_Don_base;
class Modele_turbulence_hyd_K_Omega;

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Diff_K_Omega_VEF_base
//
//////////////////////////////////////////////////////////////////////////////

class Op_Diff_K_Omega_VEF_base : public Op_Diff_K_Omega_base
{

  Declare_base_sans_constructeur(Op_Diff_K_Omega_VEF_base);

public:

  inline Op_Diff_K_Omega_VEF_base(double Prandt_K = PRDT_K_DEFAUT ,
                                  double Prandt_Omega = PRDT_OMEGA_DEFAUT );
  void completer() override;
  void associer_diffusivite(const Champ_base& ch_diff) override;
  void associer_diffusivite_turbulente() override = 0;
  inline void associer_diffusivite_turbulente(const Champ_Fonc_base&);
  inline void associer_Pr_K_Omega(double, double);
  inline const Champ_Fonc_base& diffusivite_turbulente() const;
  inline const Champ_base& diffusivite() const override
  {
    return diffusivite_.valeur();
  };

protected:
  double Prdt_K;
  double Prdt_Omega;

  // For SST model
  static constexpr double SIGMA_OMEGA1 = 0.5;
  static constexpr double SIGMA_OMEGA2 = 0.856;

  REF(Champ_Fonc_base) diffusivite_turbulente_;
  REF(Champ_Don_base) tmp;
  REF(Champ_base) diffusivite_;
  REF(Modele_turbulence_hyd_K_Omega) turbulence_model;

};

//
// Fonctions inline de la classe Op_Diff_K_Omega_VEF_base
//

inline Op_Diff_K_Omega_VEF_base::Op_Diff_K_Omega_VEF_base(double Prandt_K,
                                                          double Prandt_Omega)
  : Prdt_K(Prandt_K), Prdt_Omega(Prandt_Omega) {}

inline const Champ_Fonc_base& Op_Diff_K_Omega_VEF_base::diffusivite_turbulente() const
{
  return diffusivite_turbulente_.valeur();
}

inline void Op_Diff_K_Omega_VEF_base::associer_diffusivite_turbulente(const Champ_Fonc_base& ch)
{
  diffusivite_turbulente_=ch;
}

inline void Op_Diff_K_Omega_VEF_base::associer_Pr_K_Omega(double Pr_K, double Pr_Omega)
{
  Prdt_K = Pr_K;
  Prdt_Omega = Pr_Omega;
}

#endif
