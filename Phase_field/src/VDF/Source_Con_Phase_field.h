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
// File:        Source_Con_Phase_field.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Phase_field/src/VDF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Con_Phase_field_included
#define Source_Con_Phase_field_included


#include <Source_Con_Phase_field_base.h>
#include <Matrice_Morse.h>
#include <Ref_Zone_VDF.h>
#include <Ref_Zone_Cl_VDF.h>
//class Probleme_base;

//
// .DESCRIPTION class Source_Con_Phase_field
//

//
// .SECTION voir aussi Source_base, SouConPhase_field
//
//

class Source_Con_Phase_field : public Source_Con_Phase_field_base
{

  Declare_instanciable(Source_Con_Phase_field);

public:
  void associer_pb(const Probleme_base& );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double );
  void premier_demi_dt();

  inline const DoubleVect& get_u_carre();
  inline const double& get_drhodc();

protected:
  int tpsaff;
  double rho0;
  DoubleTab accr;
  DoubleVect u_carre_;
  double drhodc_;
  DoubleTab prov_face_,prov_elem_;
  double alpha, beta, kappa;
  double rho1, rho2;
  double mu1, mu2;
  int implicitation_;
  int gmres_;
  double epsilon_;
  double eps_;
  int boussi_;
  int diff_boussi_;
  DoubleTab g_;
  int mutype_;
  int couplage_;
  int kappa_ind;
  int type_kappa_;
  int kappa_moy_;
  double mult_kappa;
  int nkr, nit;
  double rec_min, rec_max, epsGMRES;

  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  DoubleTab& laplacien(const DoubleTab&, DoubleTab&) const;
  DoubleTab& div_kappa_grad(const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  void calculer_alpha_gradC_carre(DoubleTab&) const;
  void calculer_div_alpha_rho_gradC(DoubleTab&) const;
  void calculer_div_alpha_gradC(DoubleTab&) const;
  void calculer_pression_thermo(DoubleTab&) const;
  virtual void assembler_matrice_point_fixe(Matrice_Morse&);
  virtual void calculer_point_fixe(const DoubleTab&, const DoubleTab&, const Matrice_Morse&, DoubleTab&, DoubleTab&);
  virtual void calculer_u2_elem(DoubleVect&);
  virtual void calculer_mutilde(DoubleTab&) const;
  virtual double dpsidc(const double&) const;
  virtual void drhodc();

  virtual void construire_systeme(const DoubleTab&, const Matrice_Morse&, DoubleTab&, const DoubleTab&);
  virtual void matvect(const DoubleTab&, const Matrice_Morse&, const DoubleTab&, const DoubleTab&, DoubleTab&);
  // Modifie par DJ
  //---------------
  virtual int non_lin_gmres(const DoubleTab&, const DoubleTab&, const Matrice_Morse&, DoubleTab&, DoubleTab&);
  /*   virtual int non_lin_gmres(const DoubleTab&, const DoubleTab&, const Matrice_Morse&, DoubleTab&); */
  //---------------

  REF(Probleme_base) le_probleme2;
  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;
};

inline const DoubleVect& Source_Con_Phase_field::get_u_carre()
{
  return u_carre_;
}

inline const double& Source_Con_Phase_field::get_drhodc()
{
  return drhodc_;
}

#endif
