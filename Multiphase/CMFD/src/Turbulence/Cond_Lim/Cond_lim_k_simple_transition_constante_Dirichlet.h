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
// File:        Cond_lim_k_simple_transition_constante_Dirichlet.h
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Cond_lim_k_simple_transition_constante_Dirichlet_included
#define Cond_lim_k_simple_transition_constante_Dirichlet_included

#include <TRUSTTab.h>
#include <Dirichlet_loi_paroi.h>
#include <Ref_Correlation.h>

/*! @brief Classe Cond_lim_tau_omega_simple_demi
 *
 */
class Cond_lim_k_simple_transition_constante_Dirichlet : public Dirichlet_loi_paroi
{

  Declare_instanciable(Cond_lim_k_simple_transition_constante_Dirichlet);

public :
  int compatible_avec_eqn(const Equation_base&) const override;
  virtual int initialiser(double temps) override;
  virtual int avancer(double temps) override {return 1;}; // Avancer ne fait rien car le champ est modifie dans mettre_a_jour
  void mettre_a_jour(double tps) override;
  double calc_k(double y, double u_tau, double visc);
  virtual void liste_faces_loi_paroi(IntTab&) override;
  virtual void completer() override;

  virtual double val_imp(int i) const override {return d_(i,0);};
  virtual double val_imp(int i, int j) const override {return d_(i,j);};
  virtual double val_imp_au_temps(double temps, int i) const override {Process::exit(que_suis_je() + " : You shouldn't go through val_imp_au_temps but through val_imp ! ") ; return 1.;};
  virtual double val_imp_au_temps(double temps, int i, int j) const override {Process::exit(que_suis_je() + " : You shouldn't go through val_imp_au_temps but through val_imp ! ") ; return 1.;};

protected :
  void me_calculer();

  DoubleTab d_;

  double von_karman_ = 0.41 ;
  double beta_k_ = 0.09;
};

#endif
