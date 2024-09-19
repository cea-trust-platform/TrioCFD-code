/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Transport_Fluctuation_Temperature_W.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/UNUSED
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_Fluctuation_Temperature_W_included
#define Transport_Fluctuation_Temperature_W_included

#include <Convection_Diffusion_std.h>
#include <Op_Diff_Fluctu_Temp_Base.h>
#include <TRUST_Ref.h>

class Fluide_base;
class Modele_turbulence_scal_Fluctuation_Temperature_W;
class Motcle;

class Transport_Fluctuation_Temperature_W : public Convection_Diffusion_std
{

  Declare_instanciable_sans_constructeur(Transport_Fluctuation_Temperature_W);

public :

  Transport_Fluctuation_Temperature_W();
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void associer_milieu_base(const Milieu_base& ) override;
  void associer_modele_turbulence(const Modele_turbulence_scal_Fluctuation_Temperature_W& );
  void discretiser() override;
  void completer() override;
  inline const Milieu_base& milieu() const override;
  inline Milieu_base& milieu() override;
  inline const Modele_turbulence_scal_Fluctuation_Temperature_W& modele_turbulence() const;
  inline Modele_turbulence_scal_Fluctuation_Temperature_W& modele_turbulence();
  int nombre_d_operateurs() const override;
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  const Champ_Inc_base& inconnue() const override;
  Champ_Inc_base& inconnue() override;



  const Motcle& domaine_application() const override;
  int controler_grandeur();

protected :

  OWN_PTR(Champ_Inc_base) le_champ_Fluctu_Temperature;
  Op_Diff_Fluctu_Temp terme_diffusif;

  REF(Fluide_base)le_fluide;
  REF(Modele_turbulence_scal_Fluctuation_Temperature_W) mon_modele_Fluctu_Temp;



};

/*! @brief renvoie le champ inconnue.
 *
 */
inline Champ_Inc_base& Transport_Fluctuation_Temperature_W::inconnue()
{
  return le_champ_Fluctu_Temperature;
}

/*! @brief renvoie le champ inconnue.
 *
 */
inline const Champ_Inc_base& Transport_Fluctuation_Temperature_W::inconnue() const
{
  return le_champ_Fluctu_Temperature;
}

inline const Modele_turbulence_scal_Fluctuation_Temperature_W& Transport_Fluctuation_Temperature_W::modele_turbulence() const
{
  assert(mon_modele_Fluctu_Temp.non_nul());
  return mon_modele_Fluctu_Temp.valeur();
}

inline Modele_turbulence_scal_Fluctuation_Temperature_W& Transport_Fluctuation_Temperature_W::modele_turbulence()
{
  assert(mon_modele_Fluctu_Temp.non_nul());
  return mon_modele_Fluctu_Temp.valeur();
}

#endif


