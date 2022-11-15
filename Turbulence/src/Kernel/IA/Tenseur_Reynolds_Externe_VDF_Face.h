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
// File:        Tenseur_Reynolds_Externe_VDF_Face.h
// Directory:   $TURBULENCE_ROOT/src/Kernel/IA
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Tenseur_Reynolds_Externe_VDF_Face_included
#define Tenseur_Reynolds_Externe_VDF_Face_included

#include <Source_base.h>
#include <Terme_Source_Qdm.h>
#include <Ref_Zone_VDF.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Equation.h>
#include <Ref_Modele_turbulence_hyd_K_Eps.h>
#include <Ref_Navier_Stokes_Turbulent.h>
#include <Ref_Transport_K_Eps.h>
#include <TBNN.h>
class Probleme_base;

/*! @brief class Tenseur_Reynolds_Externe_VDF_Face
 *
 *
 *
 * @sa Source_base
 */
class Tenseur_Reynolds_Externe_VDF_Face : public Source_base, public Terme_Source_Qdm
{

  Declare_instanciable_sans_constructeur_ni_destructeur(Tenseur_Reynolds_Externe_VDF_Face);

public:
  Tenseur_Reynolds_Externe_VDF_Face();
  ~Tenseur_Reynolds_Externe_VDF_Face() override;

  void associer_pb(const Probleme_base& ) override;
  inline int has_interface_blocs() const override { return 1; }
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const override;
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const override {}
  DoubleTab& calculer(DoubleTab& ) const override;
  void mettre_a_jour(double ) override;
  void completer() override;

protected:
  void readNN();

  REF(Navier_Stokes_Turbulent)           eqn_NS_;
  REF(Modele_turbulence_hyd_K_Eps)       modele_K_Eps_;
  REF(Probleme_base)                     probleme_;
  REF(Transport_K_Eps)                   eqn_transport_K_Eps_;

  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& ) override;

  void Calcul_RSLambda();
  DoubleTab& Calcul_bij_TBNN(DoubleTab& resu) const;
  DoubleTab& Calcul_Tenseur_Reynolds( DoubleTab& ) const;

  DoubleTab lambda_1_etoile_;
  DoubleTab lambda_2_etoile_;
  DoubleTab lambda_3_etoile_;
  DoubleTab lambda_4_etoile_;
  DoubleTab lambda_5_etoile_;

  DoubleTab T1_etoile_;
  DoubleTab T2_etoile_;
  DoubleTab T3_etoile_;
  DoubleTab T4_etoile_;
  DoubleTab T5_etoile_;
  DoubleTab T6_etoile_;
  DoubleTab T7_etoile_;
  DoubleTab T8_etoile_;
  DoubleTab T9_etoile_;
  DoubleTab T10_etoile_;

  int nelem_;

  Nom nn_casename;                  // nom du reseau de neurones a charger
  TBNN *tbnn;                       // objet reseau de neurones
};

#endif
