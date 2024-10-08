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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Transport_2eq_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Transport_2eq_base_included
#define Transport_2eq_base_included

#include <Modele_turbulence_hyd_2_eq_base.h>
#include <Equation_base.h>
#include <TRUST_Ref.h>

class Modele_turbulence_hyd_2_eq_base;
class Champ_Inc_base;
class Milieu_base;

/*! @brief Classe Transport_2eq_base Classe de base pour les equations
 *
 *     de transport des modeles a deux equations type k-DISSIP
 *
 */
class Transport_2eq_base: public Equation_base
{
  Declare_base(Transport_2eq_base);

public:

  void set_param(Param&) override;
  void associer(const Equation_base&);
  inline void associer_vitesse(const Champ_base&);
  void associer_milieu_base(const Milieu_base&) override;
  Milieu_base& milieu() override;
  const Milieu_base& milieu() const override;
  double calculer_pas_de_temps() const override;

protected:

  OBS_PTR(Milieu_base) le_fluide;
  OBS_PTR(Champ_Inc_base) la_vitesse_transportante;
};

inline void Transport_2eq_base::associer_vitesse(const Champ_base& vit)
{
  la_vitesse_transportante = ref_cast(Champ_Inc_base, vit);
}

#endif
