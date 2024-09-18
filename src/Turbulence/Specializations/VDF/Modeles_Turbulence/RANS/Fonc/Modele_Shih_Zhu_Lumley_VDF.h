/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Modele_Shih_Zhu_Lumley_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Shih_Zhu_Lumley_VDF_included
#define Modele_Shih_Zhu_Lumley_VDF_included

#include <Modele_Fonc_Realisable_base.h>
#include <TRUSTTabs_forward.h>
#include <Domaine_Cl_dis.h>
#include <Param.h>
#include <Equation_base.h>
#include <TRUST_Ref.h>

#define BR_EPS 1.e-20

class Domaine_Cl_VDF;
class Domaine_VDF;



class Modele_Shih_Zhu_Lumley_VDF : public Modele_Fonc_Realisable_base
{

  Declare_instanciable(Modele_Shih_Zhu_Lumley_VDF);

public :

  virtual void set_param(Param& param);
  void mettre_a_jour(double) override;
  void Calcul_S(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse) override ;
  void Calcul_C1(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) override ;
  void Calcul_Cmu_et_S(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN) override  ;
  void Contributions_Sources(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) override ;
  void Contributions_Sources_Paroi(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                   const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) override ;

  void Calcul_C1_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override ;
  void Calcul_Cmu_et_S_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override  ;
  void Contributions_Sources_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override ;
  void Contributions_Sources_Paroi_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN,
                                       const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) override ;

  void associer(const Domaine_dis_base& , const Domaine_Cl_dis& ) override;
  void Initialisation(const Domaine_dis_base& domaine_dis) ;

protected:

  int nelem_;

  REF(Domaine_VDF) le_dom_VDF;
  REF(Domaine_Cl_VDF) le_dom_Cl_VDF;

  double A0_;
};

#endif



