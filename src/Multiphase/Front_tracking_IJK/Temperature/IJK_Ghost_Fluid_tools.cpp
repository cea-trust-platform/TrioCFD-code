/****************************************************************************
* Copyright (c) 2023, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Ghost_Fluid_tools.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Ghost_Fluid_tools.h>
#include <Probleme_base.h>
#include <DebogIJK.h>
#include <stat_counters.h>
#include <Operateur_IJK_elem_diff.h>

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >>6;
  return num_bulle;
}

static void extrapolate(const Domaine_VF& domaine_vf,
                        const IJK_Field_double& interfacial_area,
                        const IJK_Field_double& distance,
                        IJK_Field_double& field,
                        const int stencil_width)
{
  const double invalid_test = -1.e30;
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nb_faces_elem = elem_faces.dimension(1);
  const int nb_elem = elem_faces.dimension(0);
  IJK_Field_double field_old;
  const IJK_Splitting& splitting_distance = distance.get_splitting();
  /*
   * n_iterations = stencil_width is the minimum to get a propagation of information from the interface to the border
   * of the extrapolation. But doing more will lead to smoother values... And it probably costs close to nothing
   */
  const double n_iterations = 5 * stencil_width;
  for (int iteration = 0; iteration < n_iterations; iteration++)
    {
      // Copy the old field value as we do not want to use the current iteration values.
      field_old = field;
      // La valeur sur un element est la moyenne des valeurs sur les elements voisins
      for (int i_elem = 0; i_elem < nb_elem; i_elem++)
        {
          // Do not touch field in interfacial cells.
          // Iterate on other values.
          const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(i_elem);
          const double d = distance(num_elem_ijk[DIRECTION_I],
                                    num_elem_ijk[DIRECTION_J],
                                    num_elem_ijk[DIRECTION_K]);
          const double interfacial_area_elem = interfacial_area(num_elem_ijk[DIRECTION_I],
                                                                num_elem_ijk[DIRECTION_J],
                                                                num_elem_ijk[DIRECTION_K]);
          if ((d > invalid_test) && (interfacial_area_elem < invalid_test))
            {
              double sum_field = 0.;
              double coeff = 0.;
              for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                {
                  const int face = elem_faces(i_elem, i_face);
                  const int neighbour = faces_elem(face, 0) + faces_elem(face, 1) - i_elem;
                  if (neighbour >= 0)
                    {
                      // Not a boundary...
                      const Int3 num_elem_neighbour_ijk = splitting_distance.convert_packed_to_ijk_cell(neighbour);
                      double field_neighbour = field_old(num_elem_neighbour_ijk[DIRECTION_I],
                                                         num_elem_neighbour_ijk[DIRECTION_J],
                                                         num_elem_neighbour_ijk[DIRECTION_K]);
                      const double distance_neighbour = distance(num_elem_neighbour_ijk[DIRECTION_I],
                                                                 num_elem_neighbour_ijk[DIRECTION_J],
                                                                 num_elem_neighbour_ijk[DIRECTION_K]);
                      if (distance_neighbour > invalid_test)
                        {
                          // Give more weight in the smoothing to values closer to the interface:
                          if (fabs(distance_neighbour) < -1.e30)
                            {
                              Cerr << "Distance is very much at zero whereas interfacial_area is zero too... Pathological case to be looked into closely. " << finl;
                              Cerr << "Is it from a Break-up or coalescence? " << finl;
                              Cerr << "see Convection_Diffusion_Temperature_FT_Disc and static void extrapolate" << finl;
                              Cerr << "Contact TRUST support." << finl;
                              Process::exit();
                            }
                          /*
                           * TODO: Check the difference between extrapolate_champ_elem
                           * and extrapolate from Convection_Diffusion_Temperature_FT_Disc.cpp
                           */
//                          const double inv_distance_squared = 1./ (distance_neighbour * distance_neighbour);
//                          sum_field += field_neighbour * inv_distance_squared;
//                          coeff += inv_distance_squared;
                          sum_field += field_neighbour;
                          coeff++;
                        }
                    }
                }
              if (coeff > 0.)
                field(num_elem_ijk[DIRECTION_I],
                      num_elem_ijk[DIRECTION_J],
                      num_elem_ijk[DIRECTION_K]) = sum_field / coeff;
            }
        }
      field.echange_espace_virtuel(field.ghost());
    }
}

void compute_eulerian_normal_distance_field(const IJK_Interfaces& interfaces, //  ref_problem_ft_disc,
                                            IJK_Field_double& distance_field,
                                            FixedVector<IJK_Field_double, 3>& normal_vect,
                                            const int& n_iter)
{
  /*
   * Compute the normal distance to the interface
   */
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "compute_eulerian_normal_distance_field");
  statistiques().begin_count(stat_counter);

  static const double invalid_distance_value = -1.e30;
  const int dim = 3; // in IJK

  // Vertex coordinates of the eulerian domain
  const Domaine_dis_base& mon_dom_dis = interfaces.get_domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const DoubleTab& centre_element = domaine_vf.xp();

  // Approximation of the normal distance
  distance_field.data() = invalid_distance_value * 1.1;
  for (int dir = 0; dir < dim; dir++)
    normal_vect[dir].data() = 0.;

  const int nb_elem = mon_dom_dis.domaine().nb_elem();

  const Maillage_FT_IJK& maillage = interfaces.maillage_ft_ijk();

  // Same splitting for the normal vector field
  const IJK_Splitting& splitting_distance = distance_field.get_splitting();

  /*
   * M.G: Copy of B.M from Transport_Interfaces_FT_Disc::calculer_distance_interface
   * Distance calculation for the thickness 0 (vertices of the elements crossed by
   * the interface). For each element, we calculate the plane intersecting the
   * barycentre of the facets portions. The normal vector to the plane corresponds to
   * the average normal vector weighting by the surface portions. The distance
   * interface/element is the distance between this plane and the centre of the element.
   */
  {
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_elem = intersections.index_elem();
    const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
    const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
    const IntTab& facettes = maillage.facettes();
    const DoubleTab& sommets = maillage.sommets();
    // Loop on the elements
    for (int elem = 0; elem < nb_elem; elem++)
      {
        int index = index_elem[elem];
        // Moyenne ponderee des normales aux facettes qui traversent l'element
        double normale[3] = {0., 0., 0.};
        // Barycentre of the facets/element intersections
        double centre[3] = {0., 0., 0.};
        // Sum of the weights
        double surface_tot = 0.;
        // Loop on the facets which cross the element
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data =
              intersections.data_intersection(index);

            const int num_facette = data.numero_facette_;
#ifdef AVEC_BUG_SURFACES
            const double surface = data.surface_intersection_;
#else
            const double surface = data.fraction_surface_intersection_ * surface_facettes[num_facette];
#endif
            surface_tot += surface;
            for (int i = 0; i < dim; i++)
              {
                normale[i] += surface * normale_facettes(num_facette, i);
                // Barycentre calculation of the facets/element intersections
                double g_i = 0.; // i-component of the barycentre coordinate
                for (int j = 0; j < dim; j++)
                  {
                    const int som   = facettes(num_facette, j);
                    const double coord = sommets(som, i);
                    const double coeff = data.barycentre_[j];
                    g_i += coord * coeff;
                  }
                centre[i] += surface * g_i;
              }
            index = data.index_facette_suivante_;
          }
        if (surface_tot > 0.)
          {
            /*
             * The stored vector is not normed : norm ~ surface portion
             * centre = sum(centre[facet] * surface) / sum(surface)
             */
            const double inverse_surface_tot = 1. / surface_tot;
            double norme = 0.;
            int j;
            for (j = 0; j < dim; j++)
              {
                norme += normale[j] * normale[j];
                const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
                normal_vect[j](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = normale[j];
                centre[j] *= inverse_surface_tot;
              }
            if (norme > 0)
              {
                double i_norme = 1./sqrt(norme);
                double distance = 0.;
                for (j = 0; j < dim; j++)
                  {
                    double n_j = normale[j] * i_norme; // normal vector normed
                    distance += (centre_element(elem, j) - centre[j]) * n_j;
                  }
                const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
                distance_field(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = distance;
              }
          }
      }
    distance_field.echange_espace_virtuel(distance_field.ghost());
    normal_vect.echange_espace_virtuel();
    for (int dir = 0; dir < dim; dir++)
      DebogIJK::verifier("IJK_Ghost_Fluid_tools::compute_eulerian_normal_distance_field", normal_vect[dir]);
    DebogIJK::verifier("IJK_Ghost_Fluid_tools::compute_eulerian_normal_distance_field", distance_field);
  }

  FixedVector<IJK_Field_double, 3> terme_src(normal_vect);
  FixedVector<IJK_Field_double, 3> tmp(normal_vect);

  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces   = domaine_vf.elem_faces();
  const int nb_elem_voisins = elem_faces.line_size();

  // Normal vector calculation at the element location:
  int iteration;
  for (iteration = 0; iteration < n_iter; iteration++)
    {
      /*
       * Smoothing iterator, in theory:
       * normal = normal + (source_term - laplacian(normal)) * factor
       * converge towards laplacian(normal) = source_term
       * in practice:
       * normal = average(normal on neighbours) + source_term
       */
      const double un_sur_ncontrib = 1. / (1. + nb_elem_voisins);
      int elem, i, k;
      for (elem = 0; elem < nb_elem; elem++)
        {
          // Averaging the normal vector on the neighbours
          double n[3] = {0., 0., 0.};
          const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
          for (i = 0; i < dim; i++)
            {
              n[i] = normal_vect[i](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
            }
          for (k = 0; k < nb_elem_voisins; k++)
            {
              // We look for the neighbour by face k
              const int face = elem_faces(elem, k);
              const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
              const Int3 num_elem_voisin_ijk = splitting_distance.convert_packed_to_ijk_cell(e_voisin);
              if (e_voisin >= 0) // Not on a boundary
                for (i = 0; i < dim; i++)
                  n[i] += normal_vect[i](num_elem_voisin_ijk[DIRECTION_I],num_elem_voisin_ijk[DIRECTION_J],num_elem_voisin_ijk[DIRECTION_K]);
            }
          for (i = 0; i < dim; i++)
            {
              tmp[i](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) =
                terme_src[i](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) + n[i] * un_sur_ncontrib;
            }
        }
      normal_vect = tmp;
      normal_vect.echange_espace_virtuel();
    }
  /*
   * We normalise the normal vector and we create a list of elements
   * for whom the normal is known
   */
  ArrOfIntFT liste_elements;
  {
    int elem;
    for (elem = 0; elem < nb_elem; elem++)
      {
        const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
        double nx = normal_vect[0](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
        double ny = normal_vect[1](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
        double nz = normal_vect[2](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
        double norme2 = nx*nx + ny*ny + nz*nz;
        if (norme2 > 0.)
          {
            double i_norme = 1. / sqrt(norme2);
            normal_vect[0](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = nx * i_norme;
            normal_vect[1](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = ny * i_norme;
            normal_vect[2](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = nz * i_norme;
            liste_elements.append_array(elem);
          }
      }
    normal_vect.echange_espace_virtuel(); // This swap is essential
  }
  // Distance calculation at the interface
  IJK_Field_double terme_src_dist(distance_field);
  IJK_Field_double tmp_dist(distance_field);

  for (iteration = 0; iteration < n_iter; iteration++)
    {
      int i_elem, elem;
      const int liste_elem_size = liste_elements.size_array();
      for (i_elem = 0; i_elem < liste_elem_size; i_elem++)
        {
          elem = liste_elements[i_elem];
          const Int3 num_elem_ijk = splitting_distance.convert_packed_to_ijk_cell(elem);
          if (terme_src_dist(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) > invalid_distance_value)
            {
              // For all the element already crossed by the interface, the value is not computed again
              tmp_dist(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = distance_field(num_elem_ijk[DIRECTION_I],
                                                                                                                         num_elem_ijk[DIRECTION_J],
                                                                                                                         num_elem_ijk[DIRECTION_K]);
            }
          else
            {
              // For the others, we compute a distance value per neighbour
              double ncontrib = 0.;
              double somme_distances = 0.;
              int k;
              for (k = 0; k < nb_elem_voisins; k++)
                {
                  // Look for a neighbour by the face k
                  const int face = elem_faces(elem, k);
                  const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                  const Int3 num_elem_voisin_ijk = splitting_distance.convert_packed_to_ijk_cell(e_voisin);
                  if (e_voisin >= 0) // Not on a boundary
                    {
                      const double distance_voisin = distance_field(num_elem_voisin_ijk[DIRECTION_I],
                                                                    num_elem_voisin_ijk[DIRECTION_J],
                                                                    num_elem_voisin_ijk[DIRECTION_K]);
                      if (distance_voisin > invalid_distance_value)
                        {
                          // Average normal distance between an element and its neighbours
                          double nx = normal_vect[0](num_elem_ijk[DIRECTION_I],
                                                     num_elem_ijk[DIRECTION_J],
                                                     num_elem_ijk[DIRECTION_K]) +
                                      normal_vect[0](num_elem_voisin_ijk[DIRECTION_I],
                                                     num_elem_voisin_ijk[DIRECTION_J],
                                                     num_elem_voisin_ijk[DIRECTION_K]);
                          double ny = normal_vect[1](num_elem_ijk[DIRECTION_I],
                                                     num_elem_ijk[DIRECTION_J],
                                                     num_elem_ijk[DIRECTION_K]) +
                                      normal_vect[1](num_elem_voisin_ijk[DIRECTION_I],
                                                     num_elem_voisin_ijk[DIRECTION_J],
                                                     num_elem_voisin_ijk[DIRECTION_K]);
                          double nz = normal_vect[2](num_elem_ijk[DIRECTION_I],
                                                     num_elem_ijk[DIRECTION_J],
                                                     num_elem_ijk[DIRECTION_K]) +
                                      normal_vect[2](num_elem_voisin_ijk[DIRECTION_I],
                                                     num_elem_voisin_ijk[DIRECTION_J],
                                                     num_elem_voisin_ijk[DIRECTION_K]);
                          double norm2 = nx*nx + ny*ny + nz*nz;
                          if (norm2 > 0.)
                            {
                              double i_norm = 1./sqrt(norm2);
                              nx *= i_norm;
                              ny *= i_norm;
                              nz *= i_norm;
                            }
                          // Element to neighbour vector calculation
                          double dx = centre_element(elem, 0) - centre_element(e_voisin, 0);
                          double dy = centre_element(elem, 1) - centre_element(e_voisin, 1);
                          double dz = centre_element(elem, 2) - centre_element(e_voisin, 2);
                          double d = nx * dx + ny * dy + nz * dz + distance_voisin;
                          somme_distances += d;
                          ncontrib++;
                        }
                    }
                }
              // Averaging the distances obtained from neighbours
              if (ncontrib > 0.)
                {
                  double d = somme_distances / ncontrib;
                  tmp_dist(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = d;
                }
            }
        }
      distance_field = tmp_dist;
      distance_field.echange_espace_virtuel(distance_field.ghost());
    }
  statistiques().end_count(stat_counter);
}

//void compute_eulerian_normal_vector_field(const IJK_Field_double& distance, const FixedVector<IJK_Field_double, 3>& normal_vect)
//{
//  /*
//   * Compute the gradient of the normal distance field
//   */
//
//}

void compute_eulerian_curvature_field_from_distance_field(const IJK_Field_double& distance,
                                                          IJK_Field_double& curvature,
                                                          const IJK_Field_local_double& boundary_flux_kmin,
                                                          const IJK_Field_local_double& boundary_flux_kmax)
{
  /*
   * Compute the divergence of the normal vector field or the laplacian of the eulerian distance field
   */
  // Laplacian operator
  Operateur_IJK_elem_diff laplacian_distance;
  laplacian_distance.typer_diffusion_op("uniform");
  // Initialise with unit lambda
  laplacian_distance.initialize(distance.get_splitting());
  const double lambda = 1.;
  laplacian_distance.set_uniform_lambda(lambda);
  // Calculate Laplacian(dist)
  laplacian_distance.calculer(distance, curvature, boundary_flux_kmin, boundary_flux_kmax);
  const IJK_Grid_Geometry& geom = curvature.get_splitting().get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const double dz = geom.get_constant_delta(DIRECTION_K);
  const double vol = dx * dy * dz;
  const int nx = curvature.ni();
  const int ny = curvature.nj();
  const int nz = curvature.nk();
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        curvature(i,j,k) /= vol;
  curvature.echange_espace_virtuel(curvature.ghost());
}

void compute_eulerian_curvature_field_from_normal_vector_field(const FixedVector<IJK_Field_double, 3>& normal_vect,
                                                               IJK_Field_double& curvature)
{

}

void compute_eulerian_curvature_field_from_interface(const FixedVector<IJK_Field_double, 3>& normal_vect,
                                                     const IJK_Interfaces& interfaces,
                                                     IJK_Field_double& interfacial_area,
                                                     IJK_Field_double& curvature,
                                                     const int& n_iter,
                                                     const int igroup)
{
  /*
   * From IJK_Interfaces::calculer_normales_et_aires_interfaciales
   * Called in update_stat_ft IJK_FT through an instance of IJK_FT_Post !!!!!
   */

  // Vertex coordinates of the eulerian domain
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "compute_eulerian_curvature_field_from_interface");
  statistiques().begin_count(stat_counter);

  static const double invalid_curvature_value = -1.e30;

  interfacial_area.data() = invalid_curvature_value * 1.1;
  curvature.data() = invalid_curvature_value * 1.1;

  // Vertex coordinates of the eulerian domain
  const Domaine_dis_base& mon_dom_dis = interfaces.get_domaine_dis().valeur();
  const int nb_elem = mon_dom_dis.domaine().nb_elem();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  const IJK_Splitting& splitting_curvature = normal_vect.get_splitting();

  const Maillage_FT_IJK& maillage = interfaces.maillage_ft_ijk();
  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const IntTab& facettes = maillage.facettes();
  const ArrOfDouble& courbure = maillage.get_update_courbure_sommets();

  const int n_fa7 = maillage.nb_facettes();
  // Calculate the curvature in the cells crossed by the interface
  const ArrOfInt& compo_facette = maillage.compo_connexe_facettes();
  ArrOfInt compo_to_group = interfaces.get_compo_to_group();
  for (int fa7 = 0; fa7 < n_fa7; fa7++)
    {
      int icompo = compo_facette[fa7];
      if (icompo<0)
        {
          // Portion d'interface ghost. On recherche le vrai numero
          icompo = decoder_numero_bulle(-icompo);
        }
      if ((compo_to_group[icompo] != igroup) && (igroup != -1))
        continue;
      const double sf = surface_facettes[fa7];
      int index = intersections.index_facette()[fa7];
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_elem = data.numero_element_;
          const Int3 ijk = splitting_curvature.convert_packed_to_ijk_cell(num_elem);
          const double surf = data.fraction_surface_intersection_ * sf;
          for (int isom = 0; isom < 3; isom++)
            {
              const int num_som = facettes(fa7, isom);
              const double kappa = courbure[num_som];
              // No volume consideration
              if (curvature(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) < invalid_curvature_value)
                curvature(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) = kappa * surf / 3.;
              else
                curvature(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) += kappa * surf / 3.;
            }
          if (interfacial_area(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) < invalid_curvature_value)
            interfacial_area(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) = surf;
          else
            interfacial_area(ijk[DIRECTION_I], ijk[DIRECTION_J], ijk[DIRECTION_K]) += surf;
          index = data.index_element_suivant_;
        }
    }
  interfacial_area.echange_espace_virtuel(interfacial_area.ghost());
  curvature.echange_espace_virtuel(curvature.ghost());
  {
    const int nx = curvature.ni();
    const int ny = curvature.nj();
    const int nz = curvature.nk();
    for (int k=0; k < nz ; k++)
      for (int j=0; j< ny; j++)
        for (int i=0; i < nx; i++)
          if (curvature(i,j,k) > invalid_curvature_value)
            curvature(i,j,k) /= interfacial_area(i,j,k);
  }
  interfacial_area.echange_espace_virtuel(interfacial_area.ghost());
  curvature.echange_espace_virtuel(curvature.ghost());

  // Get back the cells filled with non-zero normal vectors
  ArrOfIntFT liste_elements;
  {
    for (int elem = 0; elem < nb_elem; elem++)
      {
        const Int3 num_elem_ijk = splitting_curvature.convert_packed_to_ijk_cell(elem);
        double nx = normal_vect[0](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
        double ny = normal_vect[1](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
        double nz = normal_vect[2](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
        double norme2 = nx*nx + ny*ny + nz*nz;
        if (norme2 > 0.)
          liste_elements.append_array(elem);
      }
  }
  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces   = domaine_vf.elem_faces();
  const int nb_elem_voisins = elem_faces.line_size();

  // Curvature calculation at the interface
  IJK_Field_double terme_src_curv(curvature);
  IJK_Field_double tmp_curv(curvature);

  for (int iteration = 0; iteration < n_iter; iteration++)
    {
      int i_elem, elem;
      const int liste_elem_size = liste_elements.size_array();
      for (i_elem = 0; i_elem < liste_elem_size; i_elem++)
        {
          elem = liste_elements[i_elem];
          const Int3 num_elem_ijk = splitting_curvature.convert_packed_to_ijk_cell(elem);
          if (terme_src_curv(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) > invalid_curvature_value)
            {
              // For all the element already crossed by the interface, the value is not computed again
              tmp_curv(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = curvature(num_elem_ijk[DIRECTION_I],
                                                                                                                    num_elem_ijk[DIRECTION_J],
                                                                                                                    num_elem_ijk[DIRECTION_K]);
            }
          else
            {
              // For the others, we compute a distance value per neighbour
              double ncontrib = 0.;
              double sum_kappa = 0.;
              int k;
              for (k = 0; k < nb_elem_voisins; k++)
                {
                  // Look for a neighbour by the face k
                  const int face = elem_faces(elem, k);
                  const int e_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                  const Int3 num_elem_voisin_ijk = splitting_curvature.convert_packed_to_ijk_cell(e_voisin);
                  if (e_voisin >= 0) // Not on a boundary
                    {
                      const double curvature_voisin = curvature(num_elem_voisin_ijk[DIRECTION_I],
                                                                num_elem_voisin_ijk[DIRECTION_J],
                                                                num_elem_voisin_ijk[DIRECTION_K]);
                      if (curvature_voisin > invalid_curvature_value)
                        {
                          // Average normal distance between an element and its neighbours

                          sum_kappa += curvature_voisin;
                          ncontrib++;
                        }
                    }
                }
              // Averaging the distances obtained from neighbours
              if (ncontrib > 0.)
                {
                  double kappa = sum_kappa / ncontrib;
                  tmp_curv(num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]) = kappa;
                }
            }
        }
      curvature = tmp_curv;
      curvature.echange_espace_virtuel(curvature.ghost());
    }
}

void compute_eulerian_normal_temperature_gradient_interface(const IJK_Field_double& distance,
                                                            const IJK_Field_double& indicator,
                                                            const	IJK_Field_double& temperature,
                                                            IJK_Field_double& grad_T_interface)
{
  /*
   * Compute the normal temperature gradient at the bubble interface
   */
  const int ni = temperature.ni();
  const int nj = temperature.nj();
  const int nk = temperature.nk();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        grad_T_interface(i,j,k) = temperature(i,j,k) / distance(i,j,k);
  /*
   * Check if indicatrice of the neighbours is zero + interfacial_area to locate the mixed cells
   */
}

void propagate_eulerian_normal_temperature_gradient_interface(const IJK_Interfaces& interfaces,
                                                              const IJK_Field_double& interfacial_area,
                                                              const IJK_Field_double& distance,
                                                              IJK_Field_double& temperature,
                                                              const int stencil_width)
{
  /*
   * Propagate value of grad_T_int stored in pure liquid phase towards the vapour phase and mixed cells
   */
  const Domaine_dis_base& mon_dom_dis = interfaces.get_domaine_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, mon_dom_dis);
  extrapolate(domaine_vf, interfacial_area, distance, temperature, stencil_width);
//  ArrOfIntFT liste_elements;
//  {
//    for (int elem = 0; elem < nb_elem; elem++)
//      {
//        const Int3 num_elem_ijk = splitting_curvature.convert_packed_to_ijk_cell(elem);
//        double nx = normal_vect[0](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
//        double ny = normal_vect[1](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
//        double nz = normal_vect[2](num_elem_ijk[DIRECTION_I], num_elem_ijk[DIRECTION_J], num_elem_ijk[DIRECTION_K]);
//        double norme2 = nx*nx + ny*ny + nz*nz;
//        if (norme2 > 0.)
//          liste_elements.append_array(elem);
//      }
//  }
}

void compute_eulerian_extended_temperature(const IJK_Field_double& grad_T_int,
                                           const IJK_Field_double& temperature,
                                           const IJK_Field_double& curvature,
                                           const IJK_Field_double& distance,
                                           int taylor_expansion_order)
{
  /*
   * Compute the extended temperature field using prppagated values of the temperature gradient
   */

}
