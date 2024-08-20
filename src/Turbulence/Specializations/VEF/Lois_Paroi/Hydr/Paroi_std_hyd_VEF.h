/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Paroi_std_hyd_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

// .NOM Paroi_std_hyd_VEF :
// .ENTETE TRUST VEF Turbulence
// .LIBRAIRIE libhydVEFturb
// .FILE Paroi_std_hyd_VEF.h
// .FILE Paroi_std_hyd_VEF.cpp

#ifndef Paroi_std_hyd_VEF_included
#define Paroi_std_hyd_VEF_included

#include <Paroi_hyd_base_VEF.h>
#include <distances_VEF.h>
#include <Paroi_log_QDM.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Modele_turbulence_hyd_K_Omega.h>

class Champ_Fonc_base;
class Domaine_dis;
class Domaine_Cl_dis;

/*! @brief CLASS: Paroi_std_hyd_VEF
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
void remplir_face_keps_imposee(int& flag_face_keps_imposee_,int methode_calcul_face_keps_impose_, IntVect& face_keps_imposee_, const Domaine_VEF& domaine_VEF,const REF(Domaine_Cl_VEF) le_dom_Cl_VEF,int is_champ_P1NC);


class Paroi_std_hyd_VEF : public Paroi_hyd_base_VEF, public Paroi_log_QDM
{

  Declare_instanciable_sans_constructeur(Paroi_std_hyd_VEF);

public:

  Paroi_std_hyd_VEF();
  void set_param(Param& param) override;
  int init_lois_paroi() override;
  int calculer_hyd(DoubleTab& ) override;
  int calculer_hyd_BiK(DoubleTab& , DoubleTab& ) override;
  int calculer_hyd(DoubleTab& , DoubleTab& ) override;

  void imprimer_ustar(Sortie& ) const override;
  double calculer_u_plus(const int ,const double ,const  double erugu );

  inline void check_turbulence_model();
  void compute_turbulent_quantities(double&, double&, double d_plus, double u_star, double d_visco, double dist);
  virtual int calculer_k_eps(double& , double& , double , double , double , double);

  //
  void compute_k(double& k, const double yp, const double u_star);
  void compute_epsilon(double& epsilon, const double yp, const double u_star, const double d_visco);
  void compute_omega(double& omega, const double yp, const double u_star, const double d_visco, const double dist);
  void compute_k_epsilon(double& k, double& epsilon, const double yplus, const double u_star, const double d_visco, const double dist);
  void compute_k_omega(double& k, double& omega, const double yplus, const double u_star, const double d_visco, const double dist);




protected:

  int methode_calcul_face_keps_impose_;  // 0 std avant: 1 toutes les faces accroches 2:: comme avant mais toutes les faces des elts accroches. 3: comme avant sans test si bord...  4: que les faces des elts acrroches

  DoubleVect uplus_;

  DoubleVect seuil_LP_;
  IntVect iterations_LP_;

  double u_star_impose_;
  int is_u_star_impose_;
  virtual int init_lois_paroi_hydraulique();

  int turbulence_model_type {0}; // To redirect the computation of the wall quantities

  static constexpr double BETA_OMEGA {0.075};
  static constexpr double BETA_K {0.09};  // equals to Cmu
};

KOKKOS_FUNCTION
double calculer_u_plus(const int ind_face, const double u_plus_d_plus, const double erugu, const double Kappa, DoubleArrView seuil_LP, IntArrView iterations_LP);

/*! @brief Returns an integer value depending on the turbulence model.
 *
 */
inline void Paroi_std_hyd_VEF::check_turbulence_model()
{
  turbulence_model_type = 1;
  // if (sub_type(Modele_turbulence_hyd_K_Eps, mon_modele_turb_hyd.valeur()))
  // turbulence_model_type = 1;
  if (sub_type(Modele_turbulence_hyd_K_Omega, mon_modele_turb_hyd.valeur()))
    turbulence_model_type = 2;
  // else
  // Process::exit("The turbulence model should either be K_Eps or K_Omega");
}

/*! @brief cette classe permet de specifier des options a la loi de paroi standard.
 *
 * Elle est reservee aux experts.
 *
 */
class Loi_expert_hydr_VEF : public Paroi_std_hyd_VEF
{
  Declare_instanciable(Loi_expert_hydr_VEF);
};

#endif
