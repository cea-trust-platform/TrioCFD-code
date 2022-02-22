/****************************************************************************
* Copyright (c) 2022, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_VDF_Elem_base.h
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_VDF_Elem_base_included
#define Source_Transport_VDF_Elem_base_included

#include <Calcul_Production_K_VDF.h>
#include <Ref_Equation_base.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Ref_Zone_VDF.h>
#include <Source_base.h>
#include <Modele_turbulence_hyd_K_Eps.h>
class Modele_Fonc_Bas_Reynolds;
class Probleme_base;
class Zone_Cl_dis;
class Zone_dis;

class Source_Transport_VDF_Elem_base : public Source_base, public Calcul_Production_K_VDF
{
  Declare_base_sans_constructeur( Source_Transport_VDF_Elem_base );
public :
  Source_Transport_VDF_Elem_base() { }
  Source_Transport_VDF_Elem_base(double cs1, double cs2) : C1(cs2), C2(cs2) { }
  DoubleTab& calculer(DoubleTab& ) const;
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  virtual void associer_pb(const Probleme_base& );
  virtual DoubleTab& ajouter(DoubleTab& ) const = 0;
  inline virtual void mettre_a_jour(double t) { Calcul_Production_K_VDF::mettre_a_jour(t); }

protected :
  static constexpr double C1__ = 1.44, C2__ = 1.92, C3__ = 1.0; // Chabard et N3S
  double C1 = C1__, C2 = C2__;
  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;
  REF(Equation_base) eq_hydraulique;

  void verifier_pb_keps();
  DoubleTab& ajouter_keps(DoubleTab& ) const;

private:
  // methodes a surcharger sinon throw !!
  virtual const DoubleTab& get_visc_turb() const { return not_implemented<DoubleTab&>(__func__); }
  virtual const Modele_Fonc_Bas_Reynolds& get_modele_fonc_bas_reyn() const { return not_implemented<Modele_Fonc_Bas_Reynolds&>(__func__); }
  virtual void calculer_terme_production(const Champ_Face&, const DoubleTab& , const DoubleTab& , DoubleVect&) const { return not_implemented<void>(__func__); }
  virtual void calcul_D_E(const DoubleTab& , const DoubleTab& , const Champ_Don& , DoubleTab& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void calcul_F1_F2(const Champ_base& , DoubleTab& , DoubleTab& , DoubleTab& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void fill_resu_bas_rey(const DoubleVect& , const DoubleTab& , const DoubleTab& , const DoubleTab& , const DoubleTab& , DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void fill_resu(const DoubleVect& , DoubleTab& ) const { return not_implemented<void>(__func__); }
};

#endif /* Source_Transport_VDF_Elem_base_included */
