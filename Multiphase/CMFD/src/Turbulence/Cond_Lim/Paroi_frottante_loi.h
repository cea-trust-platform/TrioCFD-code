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
// File:        Paroi_frottante_loi.h
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Paroi_frottante_loi_included
#define Paroi_frottante_loi_included

#include <DoubleTab.h>
#include <Frottement_global_impose.h>
#include <CL_loi_paroi.h>
#include <Cond_lim_base.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Paroi_frottante_loi
//    Cette condition limite correspond a un flux impose pour une condition aux limites adaptative faible de l'equation de
//    transport de QDM.
//    Le coefficient de frottement est calcule a partir de la correlation de loi de paroi adaptative.
// .SECTION voir aussi
//    Neumann
//////////////////////////////////////////////////////////////////////////////
class Paroi_frottante_loi : public Frottement_global_impose
{

  Declare_instanciable(Paroi_frottante_loi);

public :
  int compatible_avec_eqn(const Equation_base&) const;
  virtual int initialiser(double temps) ;
  virtual int avancer(double temps) {return 1;}; // Avancer ne fait rien car le champ est modifie dans mettre_a_jour
  void mettre_a_jour(double tps);
  void me_calculer();
  double calc_flux(double y, double u_tau, double visc);
  virtual double coefficient_frottement(int i) const;
  virtual double coefficient_frottement(int i,int j) const;
  virtual void liste_faces_loi_paroi(IntTab&) ;

protected :

  REF(Correlation) correlation_loi_paroi_;
  double deriv_u_plus_de_y_plus(double y_p);


// A reflechir
  DoubleTab valeurs_coeff_;
  double von_karman_ = 0.41 ;
  double beta_omega = 0.075;
  double beta_k = 0.09;
};

#endif
