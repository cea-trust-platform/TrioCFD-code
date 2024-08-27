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

#ifndef Paroi_frottante_simple_included
#define Paroi_frottante_simple_included

#include <TRUSTTab.h>
#include <Frottement_global_impose.h>
#include <Cond_lim_base.h>
#include <Correlation_base.h>
#include <TRUST_Ref.h>


/*! @brief Classe Paroi_frottante_simple Cette condition limite correspond a un flux impose pour une condition aux limites adaptative faible de l'equation de
 *
 *     transport de QDM.
 *     Le coefficient de frottement est calcule a partir de la correlation de loi de paroi adaptative.
 *
 * @sa Neumann
 */
class Paroi_frottante_simple : public Frottement_global_impose
{

  Declare_instanciable(Paroi_frottante_simple);

public :
  void completer() override ;
  virtual int initialiser(double temps) override;
  virtual int avancer(double temps) override {return 1;}; // Avancer ne fait rien car le champ est modifie dans mettre_a_jour
  void mettre_a_jour(double tps) override;
  void me_calculer();
  virtual double coefficient_frottement(int i) const override {return valeurs_coeff_(i,0);};
  virtual double coefficient_frottement(int i,int j) const override {return valeurs_coeff_(i,j);};
  virtual double coefficient_frottement_grad(int i) const override {return valeurs_coeff_grad_(i,0);};
  virtual double coefficient_frottement_grad(int i,int j) const override {return valeurs_coeff_grad_(i,j);};
  virtual void liste_faces_loi_paroi(IntTab&) override;

protected :

  virtual double fac_coeff_grad(double y_p) const { return 1.;};
  REF(Correlation) correlation_loi_paroi_;

  DoubleTab valeurs_coeff_;
  DoubleTab valeurs_coeff_grad_;

};

#endif
