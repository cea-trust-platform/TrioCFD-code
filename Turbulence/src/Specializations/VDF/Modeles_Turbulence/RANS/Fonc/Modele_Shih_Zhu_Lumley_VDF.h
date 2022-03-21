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
#include <Ref_Zone_VDF.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Zone_Cl_dis.h>
#include <Param.h>
#include <Equation_base.h>


#define BR_EPS 1.e-20

class Zone_dis;
class Zone_Cl_dis;
class DoubleVect;
class DoubleTab;
class Zone_Cl_VDF;



class Modele_Shih_Zhu_Lumley_VDF : public Modele_Fonc_Realisable_base
{

  Declare_instanciable(Modele_Shih_Zhu_Lumley_VDF);

public :

  virtual void set_param(Param& param);
  void mettre_a_jour(double) override;
  void Calcul_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse) override ;
  void Calcul_C1(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) override ;
  void Calcul_Cmu_et_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN) override  ;
  void Contributions_Sources(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) override ;
  void Contributions_Sources_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                   const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) override ;

  void Calcul_C1_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override ;
  void Calcul_Cmu_et_S_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override  ;
  void Contributions_Sources_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override ;
  void Contributions_Sources_Paroi_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN,
                                       const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) override ;

  void associer(const Zone_dis& , const Zone_Cl_dis& ) override;
  void Initialisation(const Zone_dis& zone_dis) ;

protected:

  int nelem_;

  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;

  double A0_;
};

#endif



