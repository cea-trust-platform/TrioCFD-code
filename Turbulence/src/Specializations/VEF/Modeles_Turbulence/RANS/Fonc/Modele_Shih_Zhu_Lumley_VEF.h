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
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Shih_Zhu_Lumley_VEF_included
#define Modele_Shih_Zhu_Lumley_VEF_included

#include <Modele_Fonc_Realisable_base.h>
#include <TRUSTTabs_forward.h>
#include <Domaine_Cl_dis.h>
#include <Param.h>
#include <TRUST_Ref.h>


#define BR_EPS 1.e-20

class Domaine_dis;
class Domaine_Cl_dis;
class Domaine_Cl_VEF;
class Domaine_VEF;



class Modele_Shih_Zhu_Lumley_VEF : public Modele_Fonc_Realisable_base
{

  Declare_instanciable(Modele_Shih_Zhu_Lumley_VEF);

public :

  virtual void set_param(Param& param);
  void mettre_a_jour(double) override;
  void Calcul_S(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse) override ;
  void Calcul_C1 (const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) override ;
  void Calcul_Cmu_et_S (const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN) override  ;
  // Transformation tenseur elem vers faces
  virtual DoubleTab& calcul_tenseur_face(DoubleTab&, const DoubleTab&,const Domaine_VEF&, const Domaine_Cl_VEF&) const;
  void Contributions_Sources(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) override ;
  void Contributions_Sources_Paroi(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                   const DoubleTab& visco, const DoubleTab& visco_turb,const DoubleTab& loi_paroi,const int idt) override ;

  void Calcul_C1_BiK(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override ;
  void Calcul_Cmu_et_S_BiK(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override  ;
  void Contributions_Sources_BiK(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) override ;
  void Contributions_Sources_Paroi_BiK(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN,
                                       const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) override ;

  void associer(const Domaine_dis& , const Domaine_Cl_dis& ) override;
  void init_tenseur_elem(DoubleTab&, const Domaine_VEF&, const int) const;
  void init_tenseur_face(DoubleTab&, const Domaine_VEF&, const int) const;
  void init_tenseur_elem(DoubleTab&, const Domaine_VEF&, const int) ;
  void init_tenseur_face(DoubleTab&, const Domaine_VEF&, const int) ;
  void Initialisation(const Domaine_dis& domaine_dis) ;
  void Calcul_Tenseurs_S_et_R_elem(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse) ;
  void Calcul_Tenseurs_S_et_R_elem_Paroi(const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& vitesse,
                                         const DoubleTab& visco_tab, const DoubleTab& visco_turb,
                                         const DoubleTab& tab_paroi,const int idt) ;

protected:

  int nfaces_;

  DoubleTab S_elem_;
  DoubleTab R_elem_;

  REF(Domaine_VEF) le_dom_VEF;
  REF(Domaine_Cl_VEF) le_dom_Cl_VEF;

  double A0_;
};

#endif



