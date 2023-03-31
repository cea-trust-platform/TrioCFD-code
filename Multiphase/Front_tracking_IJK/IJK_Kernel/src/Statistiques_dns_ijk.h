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

#ifndef Statistiques_dns_ijk_H
#define Statistiques_dns_ijk_H
#include <FixedVector.h>
#include <IJK_Field.h>

#include <Objet_U.h>
#include <TRUSTArrays.h>
#include <Noms.h>
#include <Param.h>

//
class IJK_Grid_Geometry;
class Statistiques_dns_ijk : public Objet_U
{
  Declare_instanciable(Statistiques_dns_ijk);
public:
  virtual Sortie& completer_print(Sortie& os) const
  {
    // porte d'entree pour completer dans la classe fille
    return os;
  };
  virtual void completer_read(Param& param)
  {
    Cerr << "Mot lu nom compris par " << que_suis_je() <<finl;
    Process::exit();
  };
  void postraiter(Sortie&, int flag_valeur_instantanee = 0) const;
  void postraiter_k(Sortie&, int flag_valeur_instantanee = 0) const; //modif AT 20/06/2013
  void update_stat(const FixedVector<IJK_Field_double, 3>& vitesse,
                   const IJK_Field_double& pression,
                   const IJK_Field_double& temperature,
                   const IJK_Field_double& masse_vol,
                   const IJK_Field_double& champ_mu,
                   const IJK_Field_double& champ_lambda,
                   const ArrOfDouble_with_ghost& delta_z_local_pour_delta,
                   const bool flag_nu_anisotropic,
                   const int flag_turbulent_viscosity,
                   const IJK_Field_double& champ_turbulent_mu_xx,
                   const IJK_Field_double& champ_turbulent_mu_xy,
                   const IJK_Field_double& champ_turbulent_mu_xz,
                   const IJK_Field_double& champ_turbulent_mu_yy,
                   const IJK_Field_double& champ_turbulent_mu_yz,
                   const IJK_Field_double& champ_turbulent_mu_zz,
                   const bool flag_kappa_anisotropic,
                   const int flag_turbulent_diffusivity,
                   const IJK_Field_double& champ_turbulent_kappa_x,
                   const IJK_Field_double& champ_turbulent_kappa_y,
                   const IJK_Field_double& champ_turbulent_kappa_z,
                   const int flag_structural_uu,
                   const FixedVector<IJK_Field_double, 6>& structural_uu_tensor,
                   const int flag_structural_uscalar,
                   const FixedVector<IJK_Field_double, 3>& structural_uscalar_vector,
                   const int flag_formulation_favre,
                   const int flag_formulation_velocity,
                   const double cp_gaz,
                   const double pression_thermodynamique,
                   double dt);
  //modif AT 20/06/2013
  void update_stat_k(const FixedVector<IJK_Field_double, 3>& vitesse,
                     const IJK_Field_double& pression,
                     //const IJK_Field_double &temperature,
                     const IJK_Field_double& masse_vol,
                     const IJK_Field_double& champ_mu,
                     //const IJK_Field_double &champ_lambda,
                     const double pression_thermodynamique,
                     const double terme_source_acceleration,
                     double dt);

  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  virtual void initialize(const IJK_Grid_Geometry&)
  {
    Cerr << "On est pas suppose pouvoir etre ici : Statistiques_dns_ijk::initialize(const IJK_Grid_Geometry &)" << finl;
    Process::exit();
  }
  virtual void initialize(const IJK_Grid_Geometry& ,const double T_KMAX, const double T_KMIN, const double constante_specifique_gaz);
  const double& t_integration() const
  {
    return t_integration_;
  }
  double t_integration_k() const
  {
    return t_integration_k_;//modif AT 20/06/2013
  }
  int check_converge() const
  {
    return check_converge_;//modif AT 20/06/2013
  }
  // accesseur pour recuperer dans un autre
  inline  VECT(ArrOfDouble) vitesse_moyenne() const
  {
    return vit_moy_;
  }
  inline ArrOfDouble masse_volumique_moyenne() const // DD 16/10/2015
  {
    return rho_moy_;
  }
  inline ArrOfDouble viscosite_cinematique_moyenne() const // DD 16/10/2015
  {
    return nu_moy_;
  }
  inline int is_converge() const
  {
    return check_converge_;
  }
  void compute_and_store_gradU_cell(const IJK_Field_double& vitesse_i,
                                    const IJK_Field_double& vitesse_j,
                                    const IJK_Field_double& vitesse_k,
                                    /* Et les champs en sortie */
                                    IJK_Field_double& dudx, IJK_Field_double& dvdy, IJK_Field_double& dwdx,
                                    IJK_Field_double& dudz, IJK_Field_double& dvdz, IJK_Field_double& dwdz);

  void cell_to_cell_gradient(const int i, const int j, const int k,
                             const IJK_Field_double& dudx, const IJK_Field_double& dvdy, const IJK_Field_double& dwdx,
                             const IJK_Field_double& dudz, const IJK_Field_double& dvdz, const IJK_Field_double& dwdz,
                             /* Et les outputs en ref aussi!! */
                             double& ddudxy, double& ddudxz, double& ddudyz,
                             double& ddvdxy, double& ddvdxz, double& ddvdyz,
                             double& ddwdxy, double& ddwdxz, double& ddwdyz) const ;

  double face_to_cell_gradient(const IJK_Field_double& vitesse_i,
                               const IJK_Field_double& vitesse_j,
                               const IJK_Field_double& vitesse_k,
                               const int i, const int j, const int k,
                               const double dz,
                               double& duidx,
                               double& dujdx,
                               double& dukdx,
                               double& duidy,
                               double& dujdy,
                               double& dukdy,
                               double& duidz,
                               double& dujdz,
                               double& dukdz,
                               const bool on_the_first_cell,
                               const bool on_the_last_cell,
                               const int bc_type) const;
  // Calcul le gradient et la derivee seconde de la vitesse :
  double calculer_gradients_vitesse(const IJK_Field_double& vitesse_i,
                                    const IJK_Field_double& vitesse_j,
                                    const IJK_Field_double& vitesse_k,
                                    const int i, const int j, const int k,
                                    const double dz,
                                    double& duidx, double&   dujdx, double&   dukdx,
                                    double& duidy, double&   dujdy, double&   dukdy,
                                    double& duidz, double&   dujdz, double&   dukdz,
                                    double& dduidxx, double& ddujdxx, double& ddukdxx,
                                    double& dduidyy, double& ddujdyy, double& ddukdyy,
                                    double& dduidzz, double& ddujdzz, double& ddukdzz,
                                    const bool on_the_first_cell,
                                    const bool on_the_last_cell) const;
  double calculer_vraie_dissipation(const double& pseudo_dissip,
                                    const double& duidx, const double& duidy, const double& duidz,
                                    const double& dujdx, const double& dujdy, const double& dujdz,
                                    const double& dukdx, const double& dukdy, const double& dukdz) const;

  double calculer_produit_scalaire_faces_to_center(const IJK_Field_double& ui,
                                                   const IJK_Field_double& uj,
                                                   const IJK_Field_double& uk,
                                                   const IJK_Field_double& vi,
                                                   const IJK_Field_double& vj,
                                                   const IJK_Field_double& vk,
                                                   const int i,
                                                   const int j,
                                                   const int k
                                                  );
  IJK_Field_double compute_and_store_scalar_product_face_to_face(
    const IJK_Field_double& v1_i,
    const IJK_Field_double& v1_j,
    const IJK_Field_double& v1_k,
    const IJK_Field_double& v2_i,
    const IJK_Field_double& v2_j,
    const IJK_Field_double& v2_k);
  void compute_vecA_minus_vecB_in_vecA(FixedVector<IJK_Field_double, 3>& vecA, const FixedVector<IJK_Field_double, 3>& vecB);
protected:
  // Z coordinates of statistics points
  ArrOfDouble elem_coord_;
  //taille des mailles;
  double dx_;//modif AT 20/06/2013
  double dy_;//modif AT 20/06/2013
  //  double dz_;//modif AT 20/06/2013
  ArrOfDouble tab_dz_;//modif AT 20/06/2013

  // Last instantaneous value of the space average (only on processor 0)
  VECT(ArrOfDouble) moyenne_spatiale_instantanee_;
  VECT(ArrOfDouble) moyenne_spatiale_ec_; // FA 17/03/2014 pour spatiales de Ec
  // Temporal integral of statistics variables
  VECT(ArrOfDouble) integrale_temporelle_;
  // Integration time
  double t_integration_;
  VECT(Nom) noms_moyennes_;
  VECT(Nom) noms_k_;
  int nval_;
  int kval_;
  double t_integration_k_;//modif AT 20/06/2013
  int check_converge_;//modif AT 20/06/2013
  VECT(ArrOfDouble) integrale_k_;//modif AT 20/06/2013
  VECT(ArrOfDouble) vit_moy_;//modif AT 20/06/2013
  ArrOfDouble rho_moy_;// DD 16/10/2015
  ArrOfDouble nu_moy_;// DD 16/10/2015
  // F.A CL pour derivee
  double TCL_kmax_;
  double TCL_kmin_;
  double constante_specifique_gaz_;
};
#endif
