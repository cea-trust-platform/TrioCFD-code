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
// File:        Paroi_std_hyd_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////



#ifndef Paroi_std_hyd_VDF_included
#define Paroi_std_hyd_VDF_included

#include <Paroi_hyd_base_VDF.h>
#include <Paroi_log_QDM.h>


class Champ_Fonc_base;
class Param;

/*! @brief CLASS: Paroi_std_hyd_VDF
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class Paroi_std_hyd_VDF : public Paroi_hyd_base_VDF, public Paroi_log_QDM
{
  Declare_instanciable(Paroi_std_hyd_VDF);

public:

  void set_param(Param& param) override;
  int init_lois_paroi() override;
  int calculer_hyd(DoubleTab& ) override; // for K_Epsilon model
  int calculer_hyd_BiK(DoubleTab&, DoubleTab& ) override; // for K_Epsilon model, splitted version
  int calculer_hyd(DoubleTab&, DoubleTab& ) override; // for K_Epsilon model
  int calculer_hyd(DoubleTab&, int isKeps, DoubleTab& ); // for K_Epsilon model
  int compute_law_komega(DoubleTab&); // for K_Omega model
  inline const DoubleVect& tab_u_plus() const { return uplus_; }
  inline double tau_tang(int face, int k) const;
  void imprimer_ustar(Sortie& ) const override;
  void calculer_uplus_dplus(DoubleVect&, DoubleVect&, DoubleVect&, int, double, double, double ) ;

protected:
  DoubleVect uplus_;

  virtual int init_lois_paroi_hydraulique();
  virtual int preparer_calcul_hyd(DoubleTab&);

  // = K_Omega model functions
  virtual int initialize_wall_law_komega(DoubleTab&);
  int compute_layer_selection(double, double, DoubleTab&, double, double, int, int);
  int compute_viscous_layer(DoubleTab& field_komega, double dist_y, double viscosity, int elem);
  int compute_buffer_layer(DoubleTab&, double, double, int, int);
  int compute_log_layer(DoubleTab&, double, int, int);
  void set_yplus_komega();

  // = Viscous sub-layer
  int calculer_u_star_sous_couche_visq(double, double, double, int);
  double calculer_u_star_sous_couche_visq(double, double, double);
  int calculer_sous_couche_visq(DoubleTab&, double, int, double, int);
  int calculer_sous_couche_visq(DoubleTab&, DoubleTab&, int);

  // = Buffer layer
  int calculer_u_star_sous_couche_tampon(double&, double, double, double, int);
  double calculer_u_star_sous_couche_tampon(double, double, double, double&);
  int calculer_sous_couche_tampon(DoubleTab&, double, double, int, int);
  int calculer_sous_couche_tampon(DoubleTab&, DoubleTab&, double, double, int, int);

  // = Inertial sub-layer
  int calculer_u_star_sous_couche_log(double, double, double, int);
  virtual double calculer_u_star(double&, double&, double&);
  int calculer_sous_couche_log(DoubleTab&, double, int, int);
  int calculer_sous_couche_log(DoubleTab&, DoubleTab&, double, int, int);

  // Identify in which part of the boundary layer we are
  int calculer_local(double, double, DoubleTab&, double, double, int, int); // k-eps version
  int calculer_local(double, double, DoubleTab&, DoubleTab&, double, double, int, int); // nu_t, tab_k
  double calculer_local(double, double, double, double, int&, double&);

  void calculer_moyennes_parois(double&, double&, double&, double&, double&, double&);
  void modifs_valeurs_turb(int, int, double, double, double, double, DoubleTab&, DoubleTab&);

  // Constant
  static constexpr double BETA1 {0.075};
  const double Cmu025 {std::pow(0.09, 0.25)}; // Cmu to the power 0.25
  const double sCmu {std::sqrt(0.09)}; // Square root of Cmu
  double ypluslam {0};
  int blended_ {0};
};

inline double Paroi_std_hyd_VDF::tau_tang(int face, int k) const
{
  if(face >= Cisaillement_paroi_.dimension(0))
    face -= le_dom_VDF->nb_faces_internes();

  if(face >= Cisaillement_paroi_.dimension_tot(0))
    {
      Cerr << "Erreur dans tau_tang " << finl;
      Cerr << "dimension : " << Cisaillement_paroi_.dimension(0) << finl;
      Cerr << "dimension_tot : " << Cisaillement_paroi_.dimension_tot(0) << finl;
      Cerr << "face : " << face << finl;
      exit();
    }
  return Cisaillement_paroi_(face, k);
}

/*! @brief cette classe permet de specifier des options a la loi de paroi standard.
 *
 * Elle est reservee aux experts.
 *
 */
class Loi_expert_hydr_VDF : public Paroi_std_hyd_VDF
{
  Declare_instanciable(Loi_expert_hydr_VDF);
};
#endif
