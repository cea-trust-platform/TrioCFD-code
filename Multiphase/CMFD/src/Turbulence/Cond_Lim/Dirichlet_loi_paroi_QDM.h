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
// File:        Dirichlet_loi_paroi_QDM.h
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Dirichlet_loi_paroi_QDM_included
#define Dirichlet_loi_paroi_QDM_included

#include <TRUSTTab.h>
#include <Dirichlet_loi_paroi.h>
#include <TRUST_Ref.h>

class Correlation;

/*! @brief Classe Dirichlet_loi_paroi_QDM
 *
 */
class Dirichlet_loi_paroi_QDM : public Dirichlet_loi_paroi
{

  Declare_instanciable(Dirichlet_loi_paroi_QDM);

public :
  int compatible_avec_eqn(const Equation_base&) const override;
  virtual int initialiser(double temps) override;
  virtual int avancer(double temps) override {return 1;}; // Avancer ne fait rien car le champ est modifie dans mettre_a_jour
  void mettre_a_jour(double tps) override;
  virtual void completer() override;

  virtual double val_imp(int i) const override {return d_(i,0);};
  virtual double val_imp(int i, int j) const override {return d_(i,j);};
  virtual double val_imp_au_temps(double temps, int i) const override {Process::exit(que_suis_je() + " : You shouldn't go through val_imp_au_temps but through val_imp ! ") ; return 1.;};
  virtual double val_imp_au_temps(double temps, int i, int j) const override {Process::exit(que_suis_je() + " : You shouldn't go through val_imp_au_temps but through val_imp ! ") ; return 1.;};

protected :
  void me_calculer();

  DoubleTab d_;

  double von_karman_ = 0.41 ;
  double beta_omega = 0.075;
  double beta_k = 0.09;

  // Initialized at 0. and adapted in the readOn depending on turbulence selected
  double y_p_prod_k_ = -1.e8 ;
  double fac_prod_k_ = -1.e8 ;
  double y_p_prod_k_grand_ = -1.e8 ;
  double fac_prod_k_grand_ = -1.e8 ;

};

#endif
