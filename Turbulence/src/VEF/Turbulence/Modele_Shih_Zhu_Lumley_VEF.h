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
// File:        Modele_Shih_Zhu_Lumley_VEF.h
// Directory:   $TRUST_ROOT/src/VEF/Turbulence
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Shih_Zhu_Lumley_VEF_included
#define Modele_Shih_Zhu_Lumley_VEF_included

#include <Modele_Fonc_Realisable_base.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>
#include <Zone_Cl_dis.h>
#include <Param.h>

#define BR_EPS 1.e-10

class Zone_dis;
class Zone_Cl_dis;
class DoubleVect;
class DoubleTab;
class Zone_Cl_VEF;



class Modele_Shih_Zhu_Lumley_VEF : public Modele_Fonc_Realisable_base
{

  Declare_instanciable(Modele_Shih_Zhu_Lumley_VEF);

public :

  virtual void set_param(Param& param);
  virtual void mettre_a_jour(double);
  virtual void Calcul_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse) ;
  virtual void Calcul_C1 (const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) ;
  virtual void Calcul_Cmu_et_S (const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN)  ;
  // Transformation tenseur elem vers faces
  virtual DoubleTab& calcul_tenseur_face(DoubleTab&, const DoubleTab&,const Zone_VEF&, const Zone_Cl_VEF&) const;
  virtual void Contributions_Sources(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) ;
  virtual void Contributions_Sources_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                           const DoubleTab& visco, const DoubleTab& visco_turb,const DoubleTab& loi_paroi,const int idt) ;

  void associer(const Zone_dis& , const Zone_Cl_dis& );
  void init_tenseur_elem(DoubleTab&, const Zone_VEF&, const int) const;
  void init_tenseur_face(DoubleTab&, const Zone_VEF&, const int) const;
  void init_tenseur_elem(DoubleTab&, const Zone_VEF&, const int) ;
  void init_tenseur_face(DoubleTab&, const Zone_VEF&, const int) ;
  void Initialisation(const Zone_dis& zone_dis) ;
  void Calcul_Tenseurs_S_et_R_elem(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse) ;
  void Calcul_Tenseurs_S_et_R_elem_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,
                                         const DoubleTab& visco_tab, const DoubleTab& visco_turb,
                                         const DoubleTab& tab_paroi,const int idt) ;

protected:

  int nfaces_;

  DoubleTab S_elem_;
  DoubleTab R_elem_;

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zone_Cl_VEF;

  double A0_;
};

#endif



