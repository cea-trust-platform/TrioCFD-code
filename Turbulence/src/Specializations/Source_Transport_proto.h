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
// File      : Source_Transport_proto.h
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_proto_included
#define Source_Transport_proto_included

#include <Ref_Convection_Diffusion_Concentration.h>
#include <Ref_Convection_Diffusion_Temperature.h>
#include <Ref_Champ_Don_base.h>
#include <Ref_Equation_base.h>
#include <Ref_Champ_Don.h>
#include <Source_base.h>

class Modele_Fonc_Realisable_base;
class Modele_Fonc_Bas_Reynolds;
class Fluide_base;

struct Source_Transport_proto
{
protected:
  Source_Transport_proto() { }
  Source_Transport_proto(double cs1, double cs2) : C1(cs1), C2(cs2) { }

  Entree& readOn_proto(Entree&, const Nom& );
  Entree& readOn_nothing(Entree& , const Nom& );
  Entree& readOn_anisotherme(Entree& , const Nom& );
  Entree& readOn_concen(Entree& , const Nom& );
  Entree& readOn_anisotherme_concen(Entree& , const Nom& );
  Entree& readOn_real(Entree& , const Nom& );
  Entree& readOn_anisotherme_real(Entree& , const Nom& );
  Entree& readOn_concen_real(Entree& , const Nom& );
  Entree& readOn_anisotherme_concen_real(Entree& , const Nom& );

  void verifier_pb_keps(const Probleme_base&, const Nom& );
  void verifier_pb_keps_anisotherme(const Probleme_base&, const Nom& );
  void verifier_pb_keps_concen(const Probleme_base&, const Nom& );
  void verifier_pb_keps_anisotherme_concen(const Probleme_base&, const Nom& );

  void verifier_milieu_anisotherme(const Probleme_base&, const Nom& );
  void verifier_milieu_concen(const Probleme_base&, const Nom& );
  void verifier_milieu_anisotherme_concen(const Probleme_base&, const Nom& );
  void verifier_beta_concen(const Fluide_base&);

  void associer_pb_proto(const Probleme_base& );
  void associer_pb_anisotherme(const Probleme_base& );
  void associer_pb_concen(const Probleme_base& );
  void associer_pb_anisotherme_concen(const Probleme_base& );

  static constexpr double C1__ = 1.44, C2__ = 1.92, C3__ = 1.0; // Chabard et N3S
  static constexpr double C21_R__ = 1.9, C3_R__ = 1.0 /* = C3__ */; // Pour realisable, Chabard et N3S
  static constexpr double C11__ = 1.55, C21__ = 2.; // pour Bas Re !
  static constexpr int interpolation_viscosite_turbulente__ = 0; // moyenne arithmetique par defaut
  double C1 = C1__, C2 = C2__, C3 = C3__;
  int _interpolation_viscosite_turbulente = interpolation_viscosite_turbulente__;

  REF(Champ_Don) beta_t, beta_c;
  REF(Champ_Don_base) gravite;
  REF(Equation_base) eq_hydraulique;
  REF(Convection_Diffusion_Temperature) eq_thermique;
  REF(Convection_Diffusion_Concentration) eq_concentration;
};

inline void error_keps(const Nom& source, const Nom& nom)
{
  Cerr << "Error ! You can't use the " << source << " source term with a " << nom << " problem/medium !!" << finl;
  Cerr << "Check the reference manual. It is may be another source term which should be used." << finl;
  Process::exit();
}

template <typename RETURN_TYPE>
RETURN_TYPE not_implemented(const char * nom_funct)
{
  Cerr << "The method " << nom_funct << " should be implemented in a derived class !" << finl;
  throw;
}

#endif /* Source_Transport_proto_included */
