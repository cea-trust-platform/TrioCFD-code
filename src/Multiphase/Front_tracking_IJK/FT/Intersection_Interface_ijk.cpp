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
// File      : Intersection_Interface_ijk.cpp
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <Intersection_Interface_ijk.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Interfaces.h>

/*
 * Intersection_Interface_ijk
 */

Implemente_base_sans_constructeur(Intersection_Interface_ijk, "Intersection_Interface_ijk", Objet_U);

Sortie& Intersection_Interface_ijk::printOn(Sortie& os) const
{
  return os;
}

Entree& Intersection_Interface_ijk::readOn(Entree& is)
{
  return is;
}

void Intersection_Interface_ijk::projete_interface(
  const Vecteur3& normale,
  const Vecteur3& bary,
  const Vecteur3& bary_vap,
  Vecteur3& posit_proj)
{
  // On a : xp = x - n (n . (x-x0))
  posit_proj = normale;
  posit_proj *= -Vecteur3::produit_scalaire(normale, bary_vap - bary);
  posit_proj += bary_vap;
}

void Intersection_Interface_ijk::distance_point_point(
  const Vecteur3& point_1,
  const Vecteur3& point_2,
  double& distance)
{
  // On a : xp = x - n (n . (x-x0))
  distance = Vecteur3::produit_scalaire(point_2 - point_1, point_2 - point_1);
  distance = std::sqrt(distance);
}

void Intersection_Interface_ijk::distance_point_point_signe(
  const Vecteur3& point_1,
  const Vecteur3& point_2,
  const Vecteur3& normal_interf,
  double& distance)
{
  // On a : xp = x - n (n . (x-x0))
  distance_point_point(point_1, point_2, distance);
  Vecteur3 interface_to_point = point_2 - point_1;
  const int sign = signbit(Vecteur3::produit_scalaire(normal_interf, interface_to_point));
  if (sign)
    distance *= -1;
}

void Intersection_Interface_ijk::get_position_interpolation_normal_interf(
  const DoubleTab& position_on_interf, const DoubleTab& normal_on_interf,
  const double dist, DoubleTab& positions)
{
  const int n_size = position_on_interf.dimension(0);
  positions.resize(n_size, 3);
  // positions = 0.;
  for (int i_diph = 0; i_diph < n_size; i_diph++)
    {
      double norm_normale = 0.;
      for (int c = 0; c < 3; c++)
        norm_normale += std::pow(normal_on_interf(i_diph, c), 2.);
      norm_normale = std::pow(norm_normale, 0.5);
      for (int c = 0; c < 3; c++)
        {
          // Cerr << "Position on interf : " << position_on_interf(i_diph, c) << ",
          // dist : " << dist << ", normal_on_interf : " << normal_on_interf(i_diph,
          // c) << ", norm_normale : " << norm_normale << finl;
          positions(i_diph, c) = position_on_interf(i_diph, c) +
                                 dist * normal_on_interf(i_diph, c) / norm_normale;
        }
    }
}

void Intersection_Interface_ijk::get_mean_interface_cell(
  const int elem, Vecteur3& normale, Vecteur3& bary) const
{
  Int3 ijk = splitting_->convert_packed_to_ijk_cell(elem);
  bool no_pb = false;
  normale = interfaces_->nI(ijk[0], ijk[1], ijk[2]);
  bary = interfaces_->xI(ijk[0], ijk[1], ijk[2]);
  for (int c=0; c < 3; c++)
    {
      if (std::abs(normale[c]) > EPS_)
        no_pb = true;
    }
  if (!no_pb)
    Cerr << "Attention, il y a un decalage entre l'indicatrice et les valeurs moyennes calculee dans la cellule" << finl;
  assert(no_pb);
}

/*
 * Intersection_Interface_ijk_face
 */

Implemente_instanciable_sans_constructeur(Intersection_Interface_ijk_face, "Intersection_Interface_ijk_face", Intersection_Interface_ijk);

Sortie& Intersection_Interface_ijk_face::printOn(Sortie& os) const
{
  return os;
}

Entree& Intersection_Interface_ijk_face::readOn(Entree& is)
{
  return is;
}

int Intersection_Interface_ijk_face::initialize(const IJK_Splitting& splitting, const IJK_Interfaces& interfaces)
{
  // TODO: est-ce que les pointeurs restent bien toujours les meme ?
  // Si oui je peux ne les initialiser qu'une fois, sinon il faudra
  // les mettre à jour.
  interfaces_ = &interfaces;
  splitting_ = &splitting;
  allocate_velocity(idiph_ijk_, splitting, 2);
//	idiph_ijk_[0].allocate(splitting, IJK_Splitting::FACES_I, 2);
//	idiph_ijk_[1].allocate(splitting, IJK_Splitting::FACES_J, 2);
//	idiph_ijk_[2].allocate(splitting, IJK_Splitting::FACES_K, 2);
  return 1;
}

void Intersection_Interface_ijk_face::compute_mean_interface_face(
  const int i, const int j, const int k, const int dir,
  Vecteur3& normale, Vecteur3& bary) const
{
  Vecteur3 bary_facettes_dans_elem0;
  Vecteur3 normale0;
  const int elem0 = splitting_->convert_ijk_cell_to_packed(i, j, k);
  get_mean_interface_cell(elem0, normale0, bary_facettes_dans_elem0);

  // elem1 est l'element avant elem0 dans la direction dir
  int elem1;
  assert((dir >= 0) && (dir < 3));
  if (dir == 0)
    elem1 = splitting_->convert_ijk_cell_to_packed(i - 1, j, k);
  else if (dir == 1)
    elem1 = splitting_->convert_ijk_cell_to_packed(i, j - 1, k);
  else
    elem1 = splitting_->convert_ijk_cell_to_packed(i, j, k - 1);

  Vecteur3 bary_facettes_dans_elem1;
  Vecteur3 normale1;
  get_mean_interface_cell(elem1, normale1, bary_facettes_dans_elem1);

  // Ici on fait une moyenne centree, mais on pourait pondérer par la distance
  // du barycentre a la face par exemple.
  normale = normale0;
  normale += normale1;
  normale *= .5;
  bary = bary_facettes_dans_elem0;
  bary += bary_facettes_dans_elem1;
  bary *= .5;
}

void Intersection_Interface_ijk_face::calcul_projection_bary_face_mouillee_interface_moy(
  DoubleTab& positions,
  IntTab& indices,
  DoubleTab& normales_de_la_proj,
  DoubleTab& distance_barys_interface
)
{
  const auto& surfaces = interfaces_->get_surface_vapeur_par_face();
  const auto& barys = interfaces_->get_barycentre_vapeur_par_face();
  const int ni = surfaces[0].ni();
  const int nj = surfaces[0].nj();
  const int nk = surfaces[0].nk();
  const double eps = 1.e-6;
  double surf_vap;
  double surf_liqu;
  int i_diph = 0;
  n_diph_ = interfaces_->get_nb_face_mouillees();
  positions.resize(2 * n_diph_, 3);
  indices.resize(2 * n_diph_, 5);
  normales_de_la_proj.resize(2 * n_diph_, 3);
  distance_barys_interface.resize(2*n_diph_,1);

  for (int dir = 0; dir < 3; dir++)
    {
      double surf_cell = 1.;
      for (int c = 0; c < 3; c++)
        {
          if (c != dir)
            {
              surf_cell *= splitting_->get_grid_geometry().get_constant_delta(c);
            }
        }
      for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++)
          for (int k = 0; k < nk; k++)
            {
              // S'il y a une face traversée par l'interface
              // Attention, je ne trouve visiblement aucune interface...
              surf_vap = surfaces[dir](i, j, k);
              if ((surf_vap/surf_cell > 0.) and (surf_vap/surf_cell < 1.))
                Cerr << "frac surf vap " << surf_vap / surf_cell << finl;
              if ((surf_vap > eps * surf_cell) and (surf_vap < (1. - eps) * surf_cell))
                {
                  // On renseigne le numéro de cette face dans le tableau des faces
                  // diphasiques
                  idiph_ijk_[dir](i, j, k) = i_diph;

                  surf_liqu = surf_cell - surf_vap;
                  // FIXME: l'assert ne passe pas, il faut trouver pourquoi ma surface
                  // est trop grande. On a dans un test surf_vap = dx et surf_cell =
                  // dx^2.
                  assert(surf_liqu > 0.);
                  // Je calcule l'équation de l'interface moyenne entre les
                  // deux cellules (de chaque coté de la face).
                  Vecteur3 bary {0., 0., .0};
                  Vecteur3 normale {0., 0., .0};
                  compute_mean_interface_face(i, j, k, dir, normale, bary);
                  // on copie normale dans le tableau final, c'etait peut etre pas
                  // absolument necessaire d'avoir recours a une copie mais tant pis.
                  // Le tableau aurait pu être 2 fois plus petit, mais dans le cas ou
                  // on utilise pas la mm normale pour les deux phases (pas
                  // implementé) cette taille est la bonne.
                  for (int c = 0; c < 3; c++)
                    {
                      normales_de_la_proj(i_diph, c) = normale[c];
                      normales_de_la_proj(i_diph + 1, c) = normale[c];
                    }

                  // On remplit le tableau des indices pour le reparcourir plus tard.
                  // Ca aurait ete pratique de créer un objet interface field
                  // avec les indices des cellules diphasiques, des faces coupées par
                  // l'interfaces etc, mais ce sera pour une autre fois.
                  for (int phase = 1; phase > 0; phase--)
                    {
                      indices(i_diph, 0) = i;
                      indices(i_diph, 1) = j;
                      indices(i_diph, 2) = k;
                      indices(i_diph, 3) = dir;
                      indices(i_diph, 4) = phase;
                    }

                  // Je calcule la projection du bary vapeur sur l'interface.
                  const Vecteur3 bary_vap
                  {
                    barys[0][dir](i, j, k), barys[1][dir](i, j, k),
                          barys[2][dir](i, j, k)
                  };
                  Vecteur3 bary_face;
                  if (dir == 0)
                    bary_face = splitting_->get_coords_of_dof(i,j,k, IJK_Splitting::FACES_I);
                  else if (dir == 1)
                    bary_face = splitting_->get_coords_of_dof(i,j,k, IJK_Splitting::FACES_J);
                  else if (dir == 2)
                    bary_face = splitting_->get_coords_of_dof(i,j,k, IJK_Splitting::FACES_K);
                  else
                    Process::exit();
                  Vecteur3 bary_liqu = bary_face - bary_vap;
                  bary_liqu *= (surf_vap / surf_liqu);
                  bary_liqu += bary_face;
                  Vecteur3 position;
                  projete_interface(normale, bary, bary_vap, position);
                  for (int c = 0; c < 3; c++)
                    positions(i_diph, c) = position[c];
                  distance_point_point(position, bary_vap, distance_barys_interface(i_diph, 0));
                  projete_interface(normale, bary, bary_liqu, position);
                  for (int c = 0; c < 3; c++)
                    positions(i_diph + 1, c) = position[c];
                  distance_point_point(position, bary_liqu, distance_barys_interface(i_diph+1, 0));
                  i_diph += 2;
                }
              else
                idiph_ijk_[dir](i, j, k) = -1;
            }
    }
  if (i_diph <2)
    Cerr << "Aucune face coupee rencontree" << finl;
}

void Intersection_Interface_ijk_face::maj_interpolation_coo_on_interfaces()
{
  if (!projected_on_interface_flag_)
    {
      calcul_projection_bary_face_mouillee_interface_moy(positions_on_interf_,
                                                         ijkf_interfaces_,
                                                         normal_on_interf_,
                                                         dist_to_interf_);
      projected_on_interface_flag_ = true;
    }
  //TODO: maj de dist_interf_
}
/*
 * Intersection_Interface_ijk_cell
 */

Implemente_instanciable_sans_constructeur(Intersection_Interface_ijk_cell, "Intersection_Interface_ijk_cell", Intersection_Interface_ijk);

Sortie& Intersection_Interface_ijk_cell::printOn(Sortie& os) const
{
  return os;
}

Entree& Intersection_Interface_ijk_cell::readOn(Entree& is)
{
  return is;
}

int Intersection_Interface_ijk_cell::initialize(
  const IJK_Splitting& splitting, const IJK_Interfaces& interfaces)
{
  // TODO: est-ce que les pointeurs restent bien toujours les meme ?
  // Si oui je peux ne les initialiser qu'une fois, sinon il faudra
  // les mettre à jour.
  interfaces_ = &interfaces;
  splitting_ = &splitting;
  projected_on_interface_flag_ = false;
  // TODO: attention, tant qu on ne passe pas dans la boucle
  // calcul_projection_centre_...
  // n_diph_ n est pas initialise
  n_diph_ = 1;
  idiph_ijk_.allocate(splitting, IJK_Splitting::ELEM, 2);
  return 1;
}

void Intersection_Interface_ijk_cell::calcul_projection_centre_sur_interface_moy(
  const IJK_Field_double& indicatrice,
  DoubleTab& positions,
  IntTab& indices,
  DoubleTab& normales_de_la_proj,
  DoubleTab& distance_centre_interface)
{
  // TODO: peut être verifier que la methode d'acces est judicieuse
  // const auto& barys = interfaces_->get_barycentre_vapeur_par_face();
  const int ni = indicatrice.ni();
  const int nj = indicatrice.nj();
  const int nk = indicatrice.nk();
  // ArrOfInt liste_composantes_connexes_dans_element;
  // liste_composantes_connexes_dans_element.set_smart_resize(1);
  // const auto &mesh = interfaces_->maillage_ft_ijk();
  Vecteur3 centre {0., 0., .0};
  Vecteur3 bary_interf {0., 0., .0};
  Vecteur3 normale_interf {0., 0., .0};
  Vecteur3 position {0., 0., .0};
  Vecteur3 distance {0., 0., .0};
  int i_diph = 0;
  // TODO: attention il faut recuperer le nombre de cellule diphasique ici.
  n_diph_ = 0;
  for (int i = 0; i < ni; i++)
    for (int j = 0; j < nj; j++)
      for (int k = 0; k < nk; k++)
        {
          //if (indicatrice(i, j, k) * (1. - indicatrice(i, j, k)) > LOCAL_EPS)
          if (fabs(indicatrice(i, j, k)) > VAPOUR_INDICATOR_TEST && fabs(indicatrice(i, j, k)) < LIQUID_INDICATOR_TEST)
            n_diph_++;
        }
  positions.resize(n_diph_, 3);
  indices.resize(n_diph_, 3);
  normales_de_la_proj.resize(n_diph_, 3);
  distance_centre_interface.resize(n_diph_,1);

  for (int i = 0; i < ni; i++)
    for (int j = 0; j < nj; j++)
      for (int k = 0; k < nk; k++)
        {
          // S'il y a une cellule traversée par l'interface
          // if (indicatrice(i, j, k) * (1. - indicatrice(i, j, k)) > LOCAL_EPS)
          if (fabs(indicatrice(i, j, k)) > VAPOUR_INDICATOR_TEST && fabs(indicatrice(i, j, k)) < LIQUID_INDICATOR_TEST)
            {
              // On renseigne le numéro de cette cellule dans le tableau des
              // cellules diphasiques
              idiph_ijk_(i, j, k) = i_diph;

              const int elem = splitting_->convert_ijk_cell_to_packed(i, j, k);
              get_mean_interface_cell(elem, normale_interf, bary_interf);
              // on copie normale dans le tableau final, c'etait peut etre pas
              // absolument necessaire d'avoir recours a une copie mais tant pis.
              for (int c = 0; c < 3; c++)
                normales_de_la_proj(i_diph, c) = normale_interf[c];

              // On remplit le tableau des indices pour le reparcourir plus tard.
              // Ca aurait ete pratique de créer un objet interface field
              // avec les indices des cellules diphasiques, des faces coupées par
              // l'interfaces etc, mais ce sera pour une autre fois.
              indices(i_diph, 0) = i;
              indices(i_diph, 1) = j;
              indices(i_diph, 2) = k;

              // Je calcule la projection du centre sur l'interface.
              centre = splitting_->get_coords_of_dof(i, j, k, IJK_Splitting::ELEM);
              projete_interface(normale_interf, bary_interf, centre, position);
              // distance_point_point(position, centre, distance_centre_interface(i_diph, 0));
              distance_point_point_signe(position, centre, normale_interf, distance_centre_interface(i_diph, 0));
              for (int c = 0; c < 3; c++)
                positions(i_diph, c) = position[c];
              i_diph++;
            }
          else
            idiph_ijk_(i, j, k) = -1;
        }
  assert(i_diph == n_diph_);
}

void Intersection_Interface_ijk_cell::calcul_projection_centre_faces_sur_interface_moy(const IJK_Field_double& indicatrice,
                                                                                       const IntTab& indices,
                                                                                       const DoubleTab& normales_interf,
                                                                                       DoubleTab& positions,
                                                                                       IntTab& indices_voisins,
                                                                                       FixedVector<DoubleVect, 3>& indices_faces_corrections,
                                                                                       DoubleTab& distance_centre_faces_interface) const
{
  Vecteur3 bary_interf {0., 0., .0};
  Vecteur3 normale_interf {0., 0., .0};
  Vecteur3 bary_face {0., 0., .0};
  Vecteur3 position {0., 0., .0};
  const int nb_diph = n_diph_; // ijk_interfaces_.size();
  positions.resize(nb_diph, 3, 6);
  indices_voisins.resize(nb_diph, 6);
  distance_centre_faces_interface.resize(nb_diph, 6);
  int neighbours_i[6] = NEIGHBOURS_I;
  int neighbours_j[6] = NEIGHBOURS_J;
  int neighbours_k[6] = NEIGHBOURS_K;
  int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
  int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
  int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;
  int nb_faces_to_correct = 0.;
  for (int i_diph=0; i_diph<nb_diph ; i_diph++)
    {
      const int i = indices(i_diph, 0);
      const int j = indices(i_diph, 1);
      const int k = indices(i_diph, 2);
      /*
       * Get back the mean interface parameters
       * FIXME: get_mean_interface_cell(elem, normale_interf, bary_interf);
       * Avoid this calculation again
       */
      const int elem = splitting_->convert_ijk_cell_to_packed(i, j, k);
      get_mean_interface_cell(elem, normale_interf, bary_interf);
      for (int l=0; l<6; l++)
        {
          const int ii = neighbours_i[l];
          const int jj = neighbours_j[l];
          const int kk = neighbours_k[l];
          const int ii_f = neighbours_faces_i[l];
          const int jj_f = neighbours_faces_j[l];
          const int kk_f = neighbours_faces_k[l];
          indices_voisins(i_diph, l) = 0;
          /*
           * Check if the neighbours is a pure liquid cell !
           */
          // if (fabs(1.-indicatrice(i+ii, j+jj, k+kk)) < LOCAL_EPS)
          if (fabs(indicatrice(i+ii, j+jj, k+kk)) > LIQUID_INDICATOR_TEST)
            {
              indices_voisins(i_diph, l) = l+1;
              nb_faces_to_correct++;
              if (ii)
                bary_face = splitting_->get_coords_of_dof(i+ii_f, j+jj_f, k+kk_f, IJK_Splitting::FACES_I);
              if (jj)
                bary_face = splitting_->get_coords_of_dof(i+ii_f, j+jj_f, k+kk_f, IJK_Splitting::FACES_J);
              if (kk)
                bary_face = splitting_->get_coords_of_dof(i+ii_f, j+jj_f, k+kk_f, IJK_Splitting::FACES_K);
              const double nx = normales_interf(i_diph, 0);
              const double ny = normales_interf(i_diph, 1);
              const double nz = normales_interf(i_diph, 2);
              normale_interf = {nx, ny, nz};
              projete_interface(normale_interf, bary_interf, bary_face, position);
              for (int c = 0; c < 3; c++)
                positions(i_diph, c, l) = position[c];
              distance_point_point_signe(position, bary_face, normale_interf, distance_centre_faces_interface(i_diph, l));
            }
        }
    }
  for (int dir=0; dir<3; dir++)
    indices_faces_corrections[dir].resize(nb_faces_to_correct);
}

void Intersection_Interface_ijk_cell::compute_face_to_correct()
{
  /*
   * FIXME: Can we use an append_array of something ?
   */
  const int nb_diph = n_diph_;
  int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
  int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
  int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;
  int nb_faces_to_correct = 0.;
  for (int i_diph=0; i_diph<nb_diph ; i_diph++)
    {
      const int i = ijk_interfaces_(i_diph, 0);
      const int j = ijk_interfaces_(i_diph, 1);
      const int k = ijk_interfaces_(i_diph, 2);
      for (int l=0; l<6; l++)
        {
          const int ii_f = neighbours_faces_i[l];
          const int jj_f = neighbours_faces_j[l];
          const int kk_f = neighbours_faces_k[l];
          if (ijk_pure_face_neighbours_(i_diph, l))
            {
              DoubleVect& i_pure_face_to_correct = ijk_pure_face_to_correct_[0];
              DoubleVect& j_pure_face_to_correct = ijk_pure_face_to_correct_[1];
              DoubleVect& k_pure_face_to_correct = ijk_pure_face_to_correct_[2];
              i_pure_face_to_correct[nb_faces_to_correct] = (i + ii_f);
              j_pure_face_to_correct[nb_faces_to_correct] = (j + jj_f);
              k_pure_face_to_correct[nb_faces_to_correct] = (k + kk_f);
              nb_faces_to_correct++;
              // i_pure_face_to_correct.append_array(i + ii_f);
              // j_pure_face_to_correct.append_array(j + jj_f);
              // k_pure_face_to_correct.append_array(k + kk_f);
            }
        }
    }
}

void Intersection_Interface_ijk_cell::update_interpolations_cell_centres_on_interface(const IJK_Field_double& interfaces)
{
  if (!projected_on_interface_flag_)
    {
      calcul_projection_centre_sur_interface_moy(interfaces,
                                                 positions_on_interf_,
                                                 ijk_interfaces_,
                                                 normal_on_interf_,
                                                 dist_to_interf_);
      projected_on_interface_flag_ = true;
    }
}

void Intersection_Interface_ijk_cell::update_interpolations_cell_faces_on_interface(const IJK_Field_double& interfaces)
{
  if (!face_centres_projected_on_interface_flag_)
    {
      calcul_projection_centre_faces_sur_interface_moy(interfaces,
                                                       ijk_interfaces_,
                                                       normal_on_interf_,
                                                       positions_pure_faces_on_interf_,
                                                       ijk_pure_face_neighbours_,
                                                       ijk_pure_face_to_correct_,
                                                       dist_pure_faces_to_interf_);
      compute_face_to_correct();
      face_centres_projected_on_interface_flag_ = true;
    }
}

void Intersection_Interface_ijk_cell::update_interpolations_cell_centres_on_interface()
{
  update_interpolations_cell_centres_on_interface(interfaces_->I());
}

void Intersection_Interface_ijk_cell::update_interpolations_cell_faces_on_interface()
{
  update_interpolations_cell_faces_on_interface(interfaces_->I());
}

void Intersection_Interface_ijk_cell::update_interpolations_cell_centres_on_interface_new()
{
  update_interpolations_cell_centres_on_interface(interfaces_->In());
}

void Intersection_Interface_ijk_cell::update_interpolations_cell_faces_on_interface_new()
{
  update_interpolations_cell_faces_on_interface(interfaces_->In());
}
