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
// File:        Mod_turb_hyd_RANS_komega.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Mod_turb_hyd_RANS_komega_included
#define Mod_turb_hyd_RANS_komega_included

#include <Motcle.h>
#include <Mod_turb_hyd_RANS_2eq.h>
class Equation_base;
class Transport_K_Omega_base;


/*! @brief Classe Mod_turb_hyd_RANS_komega Classe de base des modeles de type RANS_komega
 *
 * @sa Mod_turb_hyd_base
 */
class Mod_turb_hyd_RANS_komega: public Mod_turb_hyd_RANS_2eq
{

  Declare_base_sans_constructeur(Mod_turb_hyd_RANS_komega);

public:

  Mod_turb_hyd_RANS_komega();
  void set_param(Param& param) override;
  virtual int nombre_d_equations() const = 0;
  virtual Transport_K_Omega_base& eqn_transp_K_Omega() = 0;
  virtual const Transport_K_Omega_base& eqn_transp_K_Omega() const = 0;
  virtual const Equation_base& equation_k_omega(int) const=0 ;
  void completer() override;

  virtual void verifie_loi_paroi();
  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;

  // virtual const Equation_base& equation_k_omega(int) const=0;

  inline double get_Prandtl_K() const;
  inline double get_Prandtl_Omega() const;
  inline double get_OMEGA_MIN() const;
  inline double get_OMEGA_MAX() const;
  inline double get_K_MIN() const;
  inline int get_lquiet() const;

  //Methodes de l interface des champs postraitables
  /////////////////////////////////////////////////////
  //Methode creer_champ non codee a surcharger si necessaire
  //virtual void creer_champ(const Motcle& motlu);
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  /////////////////////////////////////////////////////

protected:
  double Prandtl_K, Prandtl_Omega; // cAlan beware! rename and put in 2eq ?
  double OMEGA_MIN, OMEGA_MAX, K_MIN;
  int lquiet;
  Motcle model_variant; // default model will be k-omega SST.

};

inline double Mod_turb_hyd_RANS_komega::get_Prandtl_K() const
{
  return Prandtl_K;
}

inline double Mod_turb_hyd_RANS_komega::get_Prandtl_Omega() const
{
  return Prandtl_Omega;
}

inline double Mod_turb_hyd_RANS_komega::get_OMEGA_MIN() const
{
  return OMEGA_MIN;
}

inline double Mod_turb_hyd_RANS_komega::get_OMEGA_MAX() const
{
  return OMEGA_MAX;
}

inline double Mod_turb_hyd_RANS_komega::get_K_MIN() const
{
  return K_MIN;
}

inline int Mod_turb_hyd_RANS_komega::get_lquiet() const
{
  return lquiet;
}

#endif
