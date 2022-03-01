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
// File:        Echange_global_impose.h
// Directory:   $TRUST_ROOT/src/ThSol
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Echange_global_impose_adaptatif_included
#define Echange_global_impose_adaptatif_included



#include <Cond_lim_base.h>
#include <Echange_impose_base.h>
#include <Echange_global_impose.h>
#include <Ref_Milieu_base.h>
#include <DoubleTab.h>
#include <Ref_Correlation.h>
#include <Correlation.h>
#include <Frontiere_dis_base.h>
#include <Ref_Frontiere_dis_base.h>
#include <Equation_base.h>


//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Echange_global_impose_adaptatif
//    Cette classe represente le cas particulier de la classe Echange_global_impose
//    pour une loi de paroi adaptative.
//    h_imp est calcule grace a la correlation de loi de paroi adaptative
//    Text est quelconque et donne par un champ front donne par l'utilisateur dans le jdd
//
//
// .SECTION voir aussi
//    Echange_impose_base Echange_externe_impose
//////////////////////////////////////////////////////////////////////////////
class Echange_global_impose_adaptatif : public Echange_global_impose
{

  Declare_instanciable(Echange_global_impose_adaptatif);
public:
  virtual void mettre_a_jour(double );
  virtual int initialiser(double temps);
  virtual int compatible_avec_eqn(const Equation_base&) const;
  int compatible_avec_discr(const Discretisation_base& discr) const {return 1;};
  virtual void liste_faces_loi_paroi(IntTab& tab);

  virtual double h_imp(int num) const       {return valeurs_h_(num, 0);};
  virtual double h_imp(int num,int k) const {return valeurs_h_(num, k);};
  virtual Champ_front& h_imp()              {Process::exit("Echange_global_impose_adaptatif : champ_front h_imp is empty !"); return h_imp_;};
  virtual const Champ_front& h_imp() const  {Process::exit("Echange_global_impose_adaptatif : champ_front h_imp is empty !"); return h_imp_;};
  virtual inline Frontiere_dis_base& frontiere_dis() {return la_frontiere_dis_;};
  virtual inline const Frontiere_dis_base& frontiere_dis() const {return la_frontiere_dis_;};

  // on va chercher les methodes Cond_lim_base qui operent sur le_champ_front mais pas sur valeurs_h_
  virtual void set_temps_defaut(double temps)                     {Cond_lim_base::set_temps_defaut(temps);};
  virtual void fixer_nb_valeurs_temporelles(int nb_cases)         {Cond_lim_base::fixer_nb_valeurs_temporelles(nb_cases);};
  virtual void changer_temps_futur(double temps,int i)            {Cond_lim_base::changer_temps_futur(temps,i);};
  virtual int avancer(double temps)                               {return Cond_lim_base::avancer(temps);};
  virtual int reculer(double temps)                               {return Cond_lim_base::reculer(temps);};
  virtual void associer_fr_dis_base(const Frontiere_dis_base& fr) {la_frontiere_dis_ = fr; le_champ_front.associer_fr_dis_base(fr);};

protected :
  DoubleTab valeurs_h_ ;
  REF(Correlation) correlation_loi_paroi_;
  REF(Frontiere_dis_base) la_frontiere_dis_;
  double mon_temps = -1;

  void me_calculer();
  double calc_theta_plus(double y, double u_tau, double visc, double cond);

  double von_karman_ = 0.41;
  double Pr_turb = 0.85;
  double Diam_hyd_ = -1;
};


#endif
