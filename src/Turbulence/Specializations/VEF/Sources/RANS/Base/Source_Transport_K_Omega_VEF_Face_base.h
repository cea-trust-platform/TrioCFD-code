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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_K_Omega_VEF_Face_base.h
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS/Base
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Transport_K_Omega_VEF_Face_base_included
#define Source_Transport_K_Omega_VEF_Face_base_included

#include <Calcul_Production_K_VEF.h>
#include <Tensors_Computation_VEF.h>
#include <Source_Transport_proto.h>
#include <TRUST_Ref.h>

class Transport_K_Omega;
class Domaine_Cl_VEF;
class Domaine_VEF;
class Modele_turbulence_hyd_K_Omega;

class Source_Transport_K_Omega_VEF_Face_base: public Source_base,
  public Source_Transport_proto,
  public Calcul_Production_K_VEF,
  public Tensors_Computation_VEF
{
  Declare_base_sans_constructeur(Source_Transport_K_Omega_VEF_Face_base);
public :
  Source_Transport_K_Omega_VEF_Face_base() { }
  Source_Transport_K_Omega_VEF_Face_base(double cs1, double cs2): Source_Transport_proto(cs1, cs2) { }

  void associer_pb(const Probleme_base& pb) override;
  void associer_domaines(const Domaine_dis_base&, const Domaine_Cl_dis_base&) override;
  DoubleTab& calculer(DoubleTab&) const override;
  DoubleTab& ajouter(DoubleTab&) const override = 0; // XXX XXX XXX Elie Saikali : like that !!;

  void elem_to_face(const Domaine_VF&, const DoubleTab& grad_elems, DoubleTab& grad_faces) const;

  inline void mettre_a_jour(double temps) override { Calcul_Production_K_VEF::mettre_a_jour(temps); }
  inline void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override { /* Do Nothing */ }

  virtual void compute_cross_diffusion(DoubleTab& gradKgradOmega) const { return not_implemented<void>(__func__); };
  virtual void compute_blending_F1(DoubleTab& gradKgradOmega) const { return not_implemented<void>(__func__); };

protected :
  DoubleTab& ajouter_komega(DoubleTab&) const;

  OBS_PTR(Domaine_VEF) le_dom_VEF;
  OBS_PTR(Domaine_Cl_VEF) le_dom_Cl_VEF;
  OBS_PTR(Modele_turbulence_hyd_K_Omega) turbulence_model;

  // Constants for the classic k-omega model Wilcox 1988
  static constexpr double BETA_K = 0.09; // Cmu or BETA_STAR, but clearer with _K
  static constexpr double BETA_OMEGA = 3./40.; // BETA
  static constexpr double SIGMA_K = 0.5; // SIGMA_STAR
  static constexpr double SIGMA_OMEGA = 0.5; // SIGMA
  static constexpr double ALPHA_OMEGA = 5./9.; // ALPHA

  // Constants for the k-omega SST model
  static constexpr double SIGMA_K1 = 0.85;
  static constexpr double SIGMA_K2 = 1.0;
  static constexpr double SIGMA_OMEGA1 = 0.5;
  static constexpr double SIGMA_OMEGA2 = 0.856;
  static constexpr double BETA1 = 0.075;
  static constexpr double BETA2 = 0.0828;
  static constexpr double KAPPA = 0.41;
  double const GAMMA1 = BETA1/BETA_K - SIGMA_OMEGA1*KAPPA*KAPPA/sqrt(BETA_K);
  double const GAMMA2 = BETA2/BETA_K - SIGMA_OMEGA2*KAPPA*KAPPA/sqrt(BETA_K);

  std::map<std::string, DoubleTab> contrib_komega;

private:
  // methodes a surcharger sinon throw !!
  virtual const DoubleTab& get_visc_turb() const { return not_implemented<DoubleTab&>(__func__); }
  virtual const DoubleTab& get_cisaillement_paroi() const { return not_implemented<DoubleTab&>(__func__); }
  virtual const DoubleTab& get_K_pour_production() const { return not_implemented<DoubleTab&>(__func__); }
  virtual const Nom get_type_paroi() const { return not_implemented<Nom>(__func__); }
  virtual void calcul_tenseur_reyn(const DoubleTab&, const DoubleTab&, DoubleTab& ) const { return not_implemented<void>(__func__); }
  virtual void fill_resu(const DoubleVect&, const DoubleTrav&, const DoubleTab&, DoubleTab& ) const { return not_implemented<void>(__func__); }
};

#endif /* Source_Transport_K_Omega_VEF_Face_base_included */
