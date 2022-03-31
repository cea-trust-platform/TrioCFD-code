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
// File:        Neumann_loi_paroi_faible_k.h
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Neumann_loi_paroi_faible_k_included
#define Neumann_loi_paroi_faible_k_included

#include <DoubleTab.h>
#include <Neumann_loi_paroi.h>
#include <CL_loi_paroi.h>
#include <Cond_lim_base.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Loi_paroi_faible_k
//    Cette condition limite correspond a un flux impose pour une condition aux limites adaptative faible de l'equation de
//    transport de k pour la turbulence.
//    Le flux impose est calcule a partir de la correlation de loi de paroi adaptative.
// .SECTION voir aussi
//    Neumann
//////////////////////////////////////////////////////////////////////////////
class Neumann_loi_paroi_faible_k : public Neumann_loi_paroi
{

  Declare_instanciable(Neumann_loi_paroi_faible_k);

public :
  int compatible_avec_eqn(const Equation_base&) const override;
  virtual int initialiser(double temps) override;
  virtual int avancer(double temps) override {return 1;}; // Avancer ne fait rien car le champ est modifie dans mettre_a_jour
  void mettre_a_jour(double tps) override;
  double calc_dyplus_kplus(double yp);
  double deriv_u_plus_de_y_plus(double y_p) ;
  double deriv_u_plus_de_y_plus_2(double y_p) ;
  virtual double flux_impose(int i) const override;
  virtual double flux_impose(int i,int j) const override;
  virtual void liste_faces_loi_paroi(IntTab&) override;
  virtual void completer() override;

protected :
  void me_calculer();

  DoubleTab valeurs_flux_;
  double von_karman_ = 0.41 ;
  double beta_omega = 0.075;
  double beta_k = 0.09;
  double limiteur_y_p = 0.01; // To prevent numerical issues ; no consequence on the calculation, as it falls in the region where the blending function is zero
};

#endif
