/****************************************************************************
 * Copyright (c) 2019, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Corrige_flux_FT.cpp
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <Corrige_flux_FT.h>
#include <DebogIJK.h>
#include <IJK_FT.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Thermique.h>
#include <Intersection_Interface_ijk.h>
#include <Param.h>
#include <stat_counters.h>

void Corrige_flux_FT::initialize(const IJK_Splitting& splitting,
                                 const IJK_Field_double& field,
                                 const IJK_Interfaces& interfaces,
                                 const IJK_FT_double& ijk_ft)
{
  // TODO: est-ce que les pointeurs restent bien toujours les meme ?
  // Si oui je peux ne les initialiser qu'une fois, sinon il faudra
  // les mettre à jour.
  interfaces_ = &interfaces;
  field_ = &field;
  splitting_ = &splitting;
  ref_ijk_ft_ = ijk_ft;
}

void Corrige_flux_FT_temperature_conv::initialize(
  const IJK_Splitting& splitting, const IJK_Field_double& field,
  const IJK_Interfaces& interfaces, const IJK_FT_double& ijk_ft,
  const double rhocpl, const double rhocpv, const double ldal,
  const double ldav)
{
  Corrige_flux_FT::initialize(splitting, field, interfaces, ijk_ft);
  intersection_ijk_face_.initialize(splitting, interfaces);
  intersection_ijk_cell_.initialize(splitting, interfaces);
  temp_interface_face_.set_smart_resize(1);
  temp_interface_cell_.set_smart_resize(1);
  temperature_barys_.set_smart_resize(1);
  temperature_ghost_.set_smart_resize(1);
  temp_interface_face_.resize(2*intersection_ijk_face_.n());
  temp_interface_cell_.resize(intersection_ijk_cell_.n());
  temperature_barys_.resize(intersection_ijk_face_.n(), 2);
  temperature_ghost_.resize(intersection_ijk_cell_.n(), 2);
  rhocp_l_ = rhocpl;
  rhocp_v_ = rhocpv;
  lda_l_ = ldal;
  lda_v_ = ldav;
  rho_cp_.initialize(rhocpl, rhocpv);
  lda_.initialize(ldal, ldav);
}

void Corrige_flux_FT_temperature_conv::update()
{

  ArrOfDouble temp_vap, temp_liqu;
  temp_vap.set_smart_resize(1);
  temp_liqu.set_smart_resize(1);

  // On commence par calculer les temperatures aux faces mouillées
  intersection_ijk_face_.maj_interpolation_coo_on_interfaces();

  temp_vap.resize(intersection_ijk_face_.n());
  temp_liqu.resize(intersection_ijk_face_.n());

  // On lance l'interpolation sur l'interface,
  calcul_temp_flux_interf_pour_bary_face(temp_vap, temp_liqu);

  // puis l'interpolation retour au bary de la face mouillée.
  interp_back_to_bary_faces(temp_vap, temp_liqu);

  // Puis des températures ghost pour les flux à proximité de l'interface
  intersection_ijk_cell_.maj_interpolation_coo_on_interfaces(
    ref_ijk_ft_->itfce().I());

  temp_vap.resize(intersection_ijk_cell_.n());
  temp_liqu.resize(intersection_ijk_cell_.n());

  // On lance l'interpolation sur l'interface,
  calcul_temp_flux_interf_pour_bary_cell(temp_vap, temp_liqu);

  // puis l'interpolation retour au bary de la face mouillée.
  update_temperature_ghost(temp_vap, temp_liqu);
}

double Corrige_flux_FT_temperature_conv::quick(
  const double Tim1,
  const double Ti,
  const double Tip1,
  const double Tip2,
  const double velocity) const
{

  const double dminfloat = std::pow(10.,-15.);
  const double dx_squared_over_8 = 0.1; //TODO: a changer evidemment
  const double curv0 = Tip1 - 2* Ti + Tim1;
  const double curv1 = Tip2 - 2* Tip1 + Ti;
  double delta0 = std::max(Tip1, Tim1) - std::min(Tip1, Tim1);
  double delta1 = std::max(Tip2, Ti) - std::min(Tip2, Ti);
  double fram0, fram1;
  if (std::abs(delta0) < dminfloat)
    {
      fram0 = 0.;
    }
  else
    {
      fram0 = std::pow((Ti - std::min(Tip1, Tim1)) / delta0 * 2. - 1., 3.);
      fram0 = std::min(fram0, 1.);
    }
  if (std::abs(delta1) < dminfloat)
    {
      fram1 = 0.;
    }
  else
    {
      fram1 = std::pow((Tip1 - std::min(Tip2, Ti)) / delta1 * 2. - 1., 3.);
      fram1 = std::min(fram1, 1.);
    }

  const double fram = std::max(fram0, fram1);
  const double curv = select_double(velocity, 0., curv1, curv0);
  const double T_amont = select_double(velocity, 0., Ti /* if velocity < 0 */, Tip1 /* if velocity > 0 */);
  double T_interp = (Ti + Tip1) * 0.5 - dx_squared_over_8 * curv;
  T_interp = (1. - fram) * T_interp + fram * T_amont;
  return T_interp;
}


void Corrige_flux_FT_temperature_conv::remplace_flux_par_quick_ghost_amont_1(
  const double frac_liquide,
  const double s_face,
  IJK_Field_local_double * const flux
) const
{
  double T_interp;

  const auto elem_i = parcours_.elem(0);
  const double velocity = ref_ijk_ft_->get_velocity(
                          )[parcours_.face()](elem_i[0], elem_i[1], elem_i[2]);
  double decal = -0.5;
  if (velocity > 0.)
    {
      decal = 0.5;
    }
  const bool le_flux_est_juste_aval_interface = is_flux_upwind_from_interface(decal);
  if (le_flux_est_juste_aval_interface)
    {
      T_interp = extrapolation_amont_1_depuis_l_interface(frac_liquide, decal);
    }
  else
    {
      T_interp = interpolation_quick_avec_1_ghost(frac_liquide, decal);
    }
  (*flux)(elem_i[0], elem_i[1], 0) = T_interp * rho_cp_(frac_liquide) * s_face * velocity;
}

bool Corrige_flux_FT_temperature_conv::is_flux_upwind_from_interface(const double decal) const
{
  //TODO
  return false;
}

double Corrige_flux_FT_temperature_conv::extrapolation_amont_1_depuis_l_interface(
  const double frac_liquide,
  const double decal) const
{
  const FixedVector<int, 3> elem = parcours_.elem(1 + (int)decal);
  const int& i = elem[0];
  const int& j = elem[1];
  const int& k = elem[2];
  const auto i_diph = intersection_ijk_cell_(i,j,k);
  // TODO: attention c'est completement faux, la température à l'interface dans le
  // tableau correspond aux températures d'interface aux projection des barycentre
  // des faces mouillées.
  const double Ti = temp_interface_cell_(i_diph);
  const double qi = q_interface_cell_(i_diph);
  // const double d = intersection_ijk_cell_.dist_interf()(i_diph);
  const double d = 0.5 * splitting_->get_grid_geometry().get_constant_delta(parcours_.face());
  const Vecteur3 norm_interf {intersection_ijk_cell_.norm_interf()(i_diph, 0), intersection_ijk_cell_.norm_interf()(i_diph, 1), intersection_ijk_cell_.norm_interf()(i_diph, 2)} ;
  const Vecteur3 norm_face {(double)parcours_.get_normale_vec()[0], (double)parcours_.get_normale_vec()[1], (double)parcours_.get_normale_vec()[2]};
  const double lda = frac_liquide * lda_l_ + (1.-frac_liquide) * lda_v_;
  return Ti + qi/lda * d * Vecteur3::produit_scalaire(norm_face, norm_interf);
}


double Corrige_flux_FT_temperature_conv::interpolation_quick_avec_1_ghost(
  const double frac_liquide,
  const double decal) const
{

  const bool is_temp_liquide = std::abs(frac_liquide - 1.) < EPS_;

  const auto elem_im1 = parcours_.elem(-1);
  const double temperature_im1 = get_ghost_temp_if_cell_is_diph(elem_im1, is_temp_liquide);

  const auto elem_i = parcours_.elem(0);
  const double temperature_i = get_ghost_temp_if_cell_is_diph(elem_i, is_temp_liquide);

  const auto elem_ip1 = parcours_.elem(1);
  const double temperature_ip1 = get_ghost_temp_if_cell_is_diph(elem_ip1, is_temp_liquide);

  const auto elem_ip2 = parcours_.elem(2);
  const double temperature_ip2 = get_ghost_temp_if_cell_is_diph(elem_ip2, is_temp_liquide);

  const double T_interp = quick(
                            temperature_im1,
                            temperature_i,
                            temperature_ip1,
                            temperature_ip2,
                            decal);
  return T_interp;
}

bool Corrige_flux_FT::test_if_stencil_inclut_bout_interface_liquide() const
{
  const int& dir= parcours_.face();

  const double velocity = ref_ijk_ft_->get_velocity(
                          )[dir](parcours_.i(), parcours_.j(), parcours_.k());
  double stencil_inclut_interface = 1.;
  int decal = 0;
  if ( velocity > 0.)
    decal = -1;

  for (int c = 0; c < 3; c++)
    {
      const auto c_elem = parcours_.elem(c + decal);
      stencil_inclut_interface *= ref_ijk_ft_->itfce().I(
                                    c_elem[0],
                                    c_elem[1],
                                    c_elem[2]);
    }

  // On considère que ce n'est pas grave s'il n'y a qu'un tout petit bout d'interface
  return std::abs(1. - stencil_inclut_interface) > 0.05;
}

bool Corrige_flux_FT::test_if_stencil_inclut_bout_interface_vapeur() const
{
  const int& dir= parcours_.face();

  const double velocity = ref_ijk_ft_->get_velocity(
                          )[dir](parcours_.i(), parcours_.j(), parcours_.k());
  double stencil_inclut_interface = 1.;
  int decal = 0;
  if ( velocity > 0.)
    decal = -1;

  for (int c = 0; c < 3; c++)
    {
      const auto c_elem = parcours_.elem(c + decal);
      stencil_inclut_interface *=  (
                                     1. - ref_ijk_ft_->itfce().I(
                                       c_elem[0],
                                       c_elem[1],
                                       c_elem[2]));
    }

  // On considère que ce n'est pas grave s'il n'y a qu'un tout petit bout d'interface
  return std::abs(1. - stencil_inclut_interface) > 0.05;
}

void Corrige_flux_FT_temperature_conv::multiplie_par_rho_cp_de_la_face_monophasique(
  const double frac_liquide,
  IJK_Field_local_double * const flux
) const
{
  const auto elem  = parcours_.elem(0);
  const double rho_cp_face = rho_cp_(frac_liquide);
  // Cerr << "Rho cp de la face : " << rho_cp_face << finl;
  (*flux)(elem[0], elem[1], 0) *= rho_cp_face;
  // flux_ij *= rho_cp_face;
}


void ParcoursIJKDir::set_indices_to_keep()
{
  if (dir_ == 0)
    {
      indices_to_keep_[0] = 1;
      indices_to_keep_[1] = 2;
    }
  else if (dir_ == 1)
    {
      indices_to_keep_[0] = 2;
      indices_to_keep_[1] = 0;
    }
  else
    {
      indices_to_keep_[0] = 0;
      indices_to_keep_[1] = 1;
    }
}

void ParcoursIJKDir::set_next_elem()
{
  if (dir_ == 0)
    {
      next_elem_[0] = 1;
      next_elem_[1] = 0;
      next_elem_[2] = 0;
    }
  else if (dir_ == 1)
    {
      next_elem_[0] = 0;
      next_elem_[1] = 1;
      next_elem_[2] = 0;
    }
  else
    {
      next_elem_[0] = 0;
      next_elem_[1] = 0;
      next_elem_[2] = 1;
    }
}

const double& Corrige_flux_FT_temperature_conv::get_ghost_temp_if_cell_is_diph(
  const FixedVector<int, 3>& elem,
  const bool from_liqu_phase) const
{
  const int& i = elem[0];
  const int& j = elem[1];
  const int& k = elem[2];
  const IJK_Field_double& indic = ref_ijk_ft_->itfce().I();
  const bool cell_is_diph = (indic(i,j,k) * (1. - indic(i,j,k)) > DMINFLOAT);
  if (cell_is_diph)
    {
      const int i_diph = intersection_ijk_cell_(i,j,k);
      return temperature_ghost_(i_diph,from_liqu_phase);
    }
  else
    {
      return (*field_)(i,j,k);
    }
}

void Corrige_flux_FT_temperature_conv::corrige_flux_faceIJ(
  IJK_Field_local_double * const flux,
  const int k_layer,
  const int dir)
{
  // on applique pour de vrai la correction au flux
  assert((dir >= 0) && (dir < 3));
  parcours_.set_dir(dir);

  const auto& surfaces = interfaces_->get_surface_vapeur_par_face();
  const double s_face = parcours_.calculer_surface_face(splitting_->get_grid_geometry());

  const int ni = field_->ni();
  const int nj = field_->nj();
  for (int i = 0; i < ni; i++)
    for (int j = 0; j < nj; j++)
      {
        parcours_.set_elem(i, j, k_layer);

        const double s_vap = surfaces[dir](i, j, k_layer);
        const double frac_vapeur = s_vap / s_face;
        const double frac_liquide = 1.- frac_vapeur;

        const bool face_monophasique = frac_liquide * frac_vapeur < EPS_;
        const bool stencil_liquide = test_if_stencil_inclut_bout_interface_liquide();
        const bool stencil_vapeur = test_if_stencil_inclut_bout_interface_vapeur();
        const bool stencil_inclut_interface = stencil_liquide & stencil_vapeur;

        multiplie_par_rho_cp_de_la_face_monophasique(frac_liquide, flux);
        // if ((face_monophasique) && (not stencil_inclut_interface))
        //   {
        //     multiplie_par_rho_cp_de_la_face_monophasique(frac_liquide, flux);
        //   }
        // else if (face_monophasique) {
        //   remplace_flux_par_quick_ghost_amont_1(frac_liquide, s_face, flux);
        // } else {
        //   remplace_flux_par_somme_rhocpf_Tf_v_Sf(frac_liquide, s_face, flux);
        // }
        // // remplace_flux_par_somme_rhocpf_Tf_v_Sf(frac_liquide, s_face, flux);
      }
}

void Corrige_flux_FT_temperature_conv::remplace_flux_par_somme_rhocpf_Tf_v_Sf(
  const double frac_liquide,
  const double s_face,
  IJK_Field_local_double  * const flux) const
{
  // La je suis bien en train de regarder une face mouillee.
  // On remet a 0 le flux, puis on ajoute les valeurs sur chaque partie
  // trouvée.
  // Cerr << "Cet elem est diphasique" << endl;
  const int& i = parcours_.i();
  const int& j = parcours_.j();
  const int& k_layer = parcours_.k();
  const int& dir = parcours_.face();
  const double velocity = ref_ijk_ft_->get_velocity()[dir](i, j, k_layer);
  const int i_diph = intersection_ijk_face_(i, j, k_layer, dir);

  // Maintenant j'ajoute les valeurs pour chaque phase (liquide et vapeur)
  (*flux)(i, j, 0) = (rhocp_l_ * temperature_barys_(i_diph, 0) * frac_liquide +
                      rhocp_v_ * temperature_barys_(i_diph, 1) * (1. - frac_liquide)) * velocity * s_face ;
}

void Corrige_flux_FT_temperature_conv::calcul_temp_flux_interf_pour_bary_face(ArrOfDouble& temp_vap, ArrOfDouble& temp_liqu)
{
  const double ldal = lda_l_;
  const double ldav = lda_v_;
  const auto& geom = splitting_->get_grid_geometry();
  const double dist = 1.52 * std::pow((
                                        std::pow(geom.get_constant_delta(0), 2.) +
                                        std::pow(geom.get_constant_delta(1), 2.) +
                                        std::pow(geom.get_constant_delta(2), 2.)), 0.5);

  DoubleTab coo_liqu1;
  DoubleTab coo_vap1;

  calcul_temperature_flux_interface(
    *field_, ldal, ldav, dist, intersection_ijk_face_.pos_interf(),
    intersection_ijk_face_.norm_interf(), temp_interface_face_, q_interface_face_,
    temp_liqu, temp_vap, coo_liqu1, coo_vap1);
}

void Corrige_flux_FT_temperature_conv::calcul_temperature_flux_interface(
  const IJK_Field_double& temperature, const double ldal, const double ldav,
  const double dist, const DoubleTab& positions, const DoubleTab& normal_on_interf,
  ArrOfDouble& temperature_interp,
  ArrOfDouble& flux_normal_interp,
  ArrOfDouble& temp_liqu,
  ArrOfDouble& temp_vap,
  DoubleTab& coo_liqu,
  DoubleTab& coo_vap)
{
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, dist, coo_liqu);
  Intersection_Interface_ijk_face::get_position_interpolation_normal_interf(
    positions, normal_on_interf, -dist, coo_vap);
  // TODO mettre une valeur acceptable dans le cas ou la temperature n est pas
  // connue
  temp_liqu.resize_array(coo_liqu.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_liqu, temp_liqu, 1.e10);
  temp_vap.resize(coo_vap.dimension(0));
  ijk_interpolate_skip_unknown_points(temperature, coo_vap, temp_vap, 1.e10);
  const int n_point_interp = positions.dimension(0);
  temperature_interp.resize_array(n_point_interp);
  flux_normal_interp.resize_array(n_point_interp);
  for (int i_point_interp = 0; i_point_interp < n_point_interp; i_point_interp++)
    {
      const double Ti =
        (temp_liqu(i_point_interp) * ldal + temp_vap(i_point_interp) * ldav) / (ldal + ldav);
      // Cerr << "Position liqu : " << coo_liqu(i_point_interp, 0) << ", " <<
      // coo_liqu(i_point_interp, 1) << ", " << coo_liqu(i_point_interp, 2) << finl; Cerr <<
      // "Position vap : " << coo_vap(i_point_interp, 0) << ", " << coo_vap(i_point_interp, 1) <<
      // ", " << coo_vap(i_point_interp, 2) << finl; Cerr << "Position interface : " <<
      // positions(i_point_interp, 0) << ", " << positions(i_point_interp, 1) << ", " <<
      // positions(i_point_interp, 2) << finl; Cerr << "Tl : " << temp_liqu(i_point_interp) << ",
      // Tv : " << temp_vap(i_point_interp) << ", Ti : " << Ti << finl;
      if (temp_liqu(i_point_interp) > 1.e9)
        {
          Cerr << "Problem temperature interface" << finl;
          temp_liqu(i_point_interp) = 1.e9;
        }
      temperature_interp(i_point_interp) = Ti;
      flux_normal_interp(i_point_interp) = ldav * (Ti - temp_vap(i_point_interp)) / dist;
    }
}

void Corrige_flux_FT_temperature_conv::interp_back_to_bary_faces(const ArrOfDouble& temp_vap, const ArrOfDouble& temp_liqu)
{
  const int n_diph = intersection_ijk_face_.n();
  const int n_point_interp = temp_vap.size_array();
  Cerr << "N diph " << n_diph << "n point inter = 2*n_diph " << n_point_interp << finl;
  assert(2*n_diph == n_point_interp);
  const auto& geom = splitting_->get_grid_geometry();
  const double d1 = 1.52 * std::pow((
                                      std::pow(geom.get_constant_delta(0), 2.) +
                                      std::pow(geom.get_constant_delta(1), 2.) +
                                      std::pow(geom.get_constant_delta(2), 2.)), 0.5);
  temperature_barys_.resize(n_diph, 2);
  // On réalise une interpolation proportionelle a la distance entre la
  // tempertaure d'interface et la température de la mm phase qui est au
  // point un peu éloignée utilisé précédemment pour calculer la temperature
  // d'interface.
  for (int i_diph = 0; i_diph < n_diph; i_diph++)
    {
      const double di_vap = std::abs(intersection_ijk_face_.dist_interf()(2*i_diph));
      const double di_liqu = std::abs(intersection_ijk_face_.dist_interf()(2*i_diph+1));
      assert(d1 - di_vap > 0.);
      assert(d1 - di_liqu > 0.);
      // La distance entre le point a l'interface et le point d'interpolation de
      // la temperature monophasique vaut bien d1 - di. Aucun cas n'est censé
      // donner une valeur négative.
      const double d1_vap_inv = 1. / (d1 - di_vap + EPS_);
      const double di_vap_inv = 1. / (di_vap + EPS_);
      const double d1_liqu_inv = 1. / (d1 - di_liqu + EPS_);
      const double di_liqu_inv = 1. / (di_liqu + EPS_);
      temperature_barys_(i_diph, 1) =
        (temp_interface_face_(2*i_diph) * di_vap_inv + temp_vap(2*i_diph) * d1_vap_inv) /
        (di_vap_inv + d1_vap_inv);
      // liquide
      temperature_barys_(i_diph, 0) =
        (temp_interface_face_(2*i_diph+1) * di_vap_inv + temp_liqu(2*i_diph+1) * d1_vap_inv) /
        (di_vap_inv + d1_vap_inv);
    }
}

void Corrige_flux_FT_temperature_conv::calcul_temp_flux_interf_pour_bary_cell(ArrOfDouble& temp_vap, ArrOfDouble& temp_liqu)
{
  const double ldal = lda_l_;
  const double ldav = lda_v_;
  const auto& geom = splitting_->get_grid_geometry();
  const double dist = 1.52 * std::pow((
                                        std::pow(geom.get_constant_delta(0), 2.) +
                                        std::pow(geom.get_constant_delta(1), 2.) +
                                        std::pow(geom.get_constant_delta(2), 2.)), 0.5);
  DoubleTab coord_vap, coord_liqu;
  calcul_temperature_flux_interface(
    *field_, ldal, ldav, dist, intersection_ijk_cell_.pos_interf(),
    intersection_ijk_cell_.norm_interf(), temp_interface_cell_, q_interface_cell_, temp_liqu,
    temp_vap, coord_liqu, coord_vap);
}

void Corrige_flux_FT_temperature_conv::update_temperature_ghost(const ArrOfDouble& temp_vap, const ArrOfDouble& temp_liqu)
{
  // intersection_ijk_cell_.maj_interpolation_coo_on_interfaces(
  //   ref_ijk_ft_->itfce().I());
  const double ldal = lda_l_;
  const double ldav = lda_v_;
  const int n_diph = temp_interface_cell_.size_array();
  assert(n_diph == intersection_ijk_cell_.n());
  temperature_ghost_.resize(n_diph, 2);
  // On réalise une interpolation proportionelle a la distance entre la
  // temperature d'interface au centre de la cellule diph avec tantot le
  // gradient liquide ou vapeur.
  for (int i_diph = 0; i_diph < n_diph; i_diph++)
    {
      const double di = intersection_ijk_cell_.dist_interf()(i_diph);
      temperature_ghost_(i_diph, 0) =
        temp_interface_cell_(i_diph) + di * q_interface_cell_(i_diph) / ldav;
      temperature_ghost_(i_diph, 1) =
        temp_interface_cell_(i_diph) + di * q_interface_cell_(i_diph) / ldal;
    }
}
