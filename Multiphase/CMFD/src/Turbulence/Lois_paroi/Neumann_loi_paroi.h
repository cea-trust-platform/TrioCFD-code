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
// File:        Neumann_loi_paroi.h
// Directory:   $TRUST_ROOT/src/CoviMAC/Cond_Lim
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Neumann_loi_paroi_included
#define Neumann_loi_paroi_included

#include <Neumann_paroi.h>
#include <Param.h>
#include <IntTab.h>
#include <Ref_Correlation.h>
#include <Correlation.h>
#include <Frontiere_dis_base.h>
#include <Ref_Frontiere_dis_base.h>


//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Neumann_loi_paroi
//    Classe de base pour les flux impose pour une condition aux limites adaptative faible des equations de turbulence
//    Le flux impose est calcule a partir de la correlation de loi de paroi adaptative.
// .SECTION voir aussi
//    Neumann
//////////////////////////////////////////////////////////////////////////////
class Neumann_loi_paroi : public Neumann_paroi
{

  Declare_base(Neumann_loi_paroi);

public:

  virtual void liste_faces_loi_paroi(IntTab&) =0;
  virtual void associer_correlation(const Correlation& corr)
  {

    correlation_loi_paroi_ = corr;

  };
  virtual void associer_fr_dis_base(const Frontiere_dis_base& fr) {la_frontiere_dis=fr;};
  virtual int initialiser(double temps) =0 ;
  virtual void associer_zone_cl_dis_base(const Zone_Cl_dis_base& zcl) { ma_zone_cl_dis=zcl;};

  // fonctions de cond_lim_base qui necessitent le champ_front qu'on met a zero car on fait abstraction du champ_front
  virtual void completer() {};
  virtual void fixer_nb_valeurs_temporelles(int nb_cases) {};
  virtual inline Frontiere_dis_base& frontiere_dis() {return la_frontiere_dis;};
  virtual inline const Frontiere_dis_base& frontiere_dis() const {return la_frontiere_dis;};
  virtual void changer_temps_futur(double temps,int i) {};
  virtual void set_temps_defaut(double temps) {};
  virtual void calculer_coeffs_echange(double temps) {};
  virtual void verifie_ch_init_nb_comp() {};


protected:
  REF(Correlation) correlation_loi_paroi_;
  REF(Frontiere_dis_base) la_frontiere_dis;
  double mon_temps = -1;
};

#endif
