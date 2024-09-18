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
// File:        Modele_Lam_Bremhorst_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Lam_Bremhorst_VEF_included
#define Modele_Lam_Bremhorst_VEF_included

#include <Modele_Fonc_Bas_Reynolds_Base.h>
#include <TRUSTTabs_forward.h>

#include <Param.h>
#include <TRUST_Ref.h>


#define A1 3.9
#define CNL1 0.8
#define CNL2 11.
#define CNL3 4.5
#define CNL4 1000.
#define CNL5 1.
/* NASA Contractor Report 4145
Two-Equation Low-Reynolds-Number
Turbulence Modeling of Transitional
Boundary Layer Flows Characteristic
of Gas Turbine Blades
Rodney C. Schmidt and Suhas V. Patankar
*/
#define BR_EPS 1.e-20

class Domaine_Cl_VEF;
class Domaine_VEF;

class Modele_Lam_Bremhorst_VEF : public Modele_Fonc_Bas_Reynolds_Base
{

  Declare_instanciable(Modele_Lam_Bremhorst_VEF);

public :


  virtual void set_param(Param& param);
  DoubleTab& Calcul_D(DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const override;
  DoubleTab& Calcul_E(DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const override ;
  DoubleTab& Calcul_Cmu(DoubleTab&,
                        const Domaine_dis_base&,  const Domaine_Cl_dis_base&,
                        const DoubleTab&, const DoubleTab&, const double) const override;
  DoubleTab& Calcul_Cmu_Paroi(DoubleTab&,
                              const Domaine_dis_base&, const Domaine_Cl_dis_base&,
                              const DoubleTab& , const DoubleTab& ,
                              const DoubleTab& ,const int,
                              const DoubleTab&, const DoubleTab&,
                              const double) const override;
  //  virtual DoubleTab& Calcul_F1(DoubleTab&, const Domaine_dis_base& ) const ;
  DoubleTab& Calcul_F1( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& P,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const override;
  int Calcul_is_Reynolds_stress_isotrope() const override;
  int Calcul_is_Cmu_constant() const override;
  void init_tenseur_elem(DoubleTab&, const Domaine_VEF&, const int) const;
  void init_tenseur_face(DoubleTab&, const Domaine_VEF&, const int) const;
  virtual DoubleTab calcul_norme_elem(const Domaine_VEF&, const DoubleTab) const;
  // Transformation tenseur elem vers faces
  virtual DoubleTab& calcul_tenseur_face(DoubleTab&, const DoubleTab&,
                                         const Domaine_VEF&, const Domaine_Cl_VEF&) const;

  bool calcul_tenseur_Re(const DoubleTab&, const DoubleTab&, DoubleTab&) const override;
  virtual DoubleTab calcul_tenseur_Re_elem(const Discretisation_base& dis,
                                           const Domaine_dis_base&,const DoubleTab&, const DoubleTab&, const DoubleTab&,
                                           const DoubleTab&, const Champ_base& K_Eps) const;
  DoubleTab& Calcul_F2(DoubleTab&, DoubleTab&,const Domaine_dis_base&,const DoubleTab&,const Champ_base&) const override ;
  DoubleTab& Calcul_Fmu ( DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const Champ_Don& )const override ;
  Entree& lire(const Motcle&, Entree&);
  // void associer_pb(const Probleme_base& );
  void associer(const Domaine_dis_base& , const Domaine_Cl_dis_base& ) override;
  void mettre_a_jour(double) override;
  void lire_distance_paroi( ) override;
  DoubleTab calcul_tenseur_Re_elem_shih(const Discretisation_base&,
                                        const Domaine_dis_base&,
                                        const DoubleTab&, const DoubleTab&,
                                        const DoubleTab&, const DoubleTab& ,const DoubleTab& ,
                                        const Champ_base& K_Eps) const;
  DoubleTab calcul_tenseur_Re_shih(const Discretisation_base& dis,
                                   const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,
                                   const DoubleTab& G,
                                   const Champ_base& K_Eps) const;



  DoubleTab& Calcul_Cmu_BiK(DoubleTab&,const Domaine_dis_base&, const Domaine_Cl_dis_base&, const DoubleTab&, const DoubleTab&, const DoubleTab&, const double) const override;
  DoubleTab& Calcul_Cmu_Paroi_BiK(DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&,
                                  const DoubleTab& , const DoubleTab& ,
                                  const DoubleTab& ,const int,
                                  const DoubleTab&, const DoubleTab&, const DoubleTab&,
                                  const double) const override;
  bool calcul_tenseur_Re_BiK(const DoubleTab&, const DoubleTab&, DoubleTab&) const override;
  virtual DoubleTab calcul_tenseur_Re_elem_BiK(const Discretisation_base& dis,
                                               const Domaine_dis_base&,const DoubleTab&, const DoubleTab&, const DoubleTab&,
                                               const DoubleTab&, const Champ_base& K, const Champ_base& Eps) const;

  DoubleTab& Calcul_D_BiK(DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const override;
  DoubleTab& Calcul_E_BiK(DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const override ;
  DoubleTab& Calcul_F1_BiK( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& P,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const override;
  DoubleTab& Calcul_F2_BiK(DoubleTab&, DoubleTab&,const Domaine_dis_base&,const DoubleTab&,const DoubleTab&,const Champ_base&) const override ;
  DoubleTab& Calcul_Fmu_BiK ( DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const Champ_Don& )const override ;

protected:

  REF(Domaine_VEF) le_dom_VEF;
  REF(Domaine_Cl_VEF) le_dom_Cl_VEF;
  /*   REF(Fluide_Incompressible) le_fluide; */
  /*   REF(Champ_Inc) la_vitesse_transportante; */
  /*   REF(Transport_K_Eps_Bas_Reynolds) eq_transport_K_Eps_Bas_Re; */
  /*   REF(Champ_Don) visco; */

};

#endif



