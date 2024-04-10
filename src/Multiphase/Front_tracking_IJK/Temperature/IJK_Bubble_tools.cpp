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
// File      : IJK_Bubble_tools.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Bubble_tools.h>


static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >> 6;
  return num_bulle;
}

std::vector<int> arg_sort_array(const ArrOfDouble& array_to_sort)
{
  const int n = array_to_sort.size_array();
  // IntVect indices(n);
  std::vector<int> indices(n);
  for (int j=0; j<n; j++)
    indices[j]=j;
  std::sort(indices.begin(), indices.end(), [&array_to_sort](int i, int j) {return array_to_sort[i] < array_to_sort[j];});
  return indices;
}

std::vector<int> arg_sort_array_phi(const ArrOfDouble& angle_incr, const ArrOfDouble& first_angle, const ArrOfDouble& array_to_sort)
{
  const int n = array_to_sort.size_array();
  const double constant_angle_incr = abs(angle_incr[1] - angle_incr[0]);
  std::vector<int> indices;
  std::vector<int> indices_subarray;
  ArrOfDouble sub_array;
  for (int i=0; i<angle_incr.size_array(); i++)
    {
      const double angle_min = angle_incr(i) - constant_angle_incr / 2;
      const double angle_max = angle_incr(i) + constant_angle_incr / 2;
      for (int j=0; j<first_angle.size_array(); j++)
        if (first_angle(j) >= angle_min && first_angle(j) < angle_max)
          {
            indices_subarray.push_back(j);
            sub_array.append_array(array_to_sort(j));
          }
      std::vector<int> indices_subarray_sorted = arg_sort_array(sub_array);
      const int size_subarray_sorted = (int) indices_subarray_sorted.size();
      for (int k=0; k<size_subarray_sorted; k++)
        {
          const int index = indices_subarray_sorted[k];
          indices.push_back((int) indices_subarray[index]);
        }
      indices_subarray.clear();
      sub_array.reset();
    }

  assert((int) indices.size() == n);
  if (!((int) indices.size() == n))
    Process::exit();

  return indices;
}

/* FROM void IJK_Interfaces::calculer_volume_bulles
 * L'index de la bulle ghost est (entre -1 et -nbulles_ghost):
 * const int idx_ghost = get_ghost_number_from_compo(compo);
 * // On la place en fin de tableau :
 * compo = nbulles_reelles - 1 - idx_ghost;
 */

void compute_bounding_box_fill_compo(const IJK_Interfaces& interfaces,
                                     DoubleTab& bounding_box,
                                     DoubleTab& min_max_larger_box,
                                     IJK_Field_double& eulerian_compo_connex,
                                     IJK_Field_double& eulerian_compo_connex_ghost,
                                     DoubleTab& bubbles_barycentre)
{
  /*
  * bounding_box(b, dir, m) :
  * b -> Numero de la composante connexe de la bulle.
  * dir -> Direction i,j ou k.
  * m   -> min (0) ou max (1)
  */
  interfaces.calculer_bounding_box_bulles(bounding_box);
  int nb_bubbles = interfaces.get_nb_bulles_reelles();
  int nb_ghost_bubbles = interfaces.get_nb_bulles_ghost();
  eulerian_compo_connex.data() = -1;
  eulerian_compo_connex.echange_espace_virtuel(eulerian_compo_connex.ghost());
  eulerian_compo_connex_ghost.data() = -1;
  eulerian_compo_connex_ghost.echange_espace_virtuel(eulerian_compo_connex_ghost.ghost());
  IntTab ghost_to_real_bubble(nb_ghost_bubbles);
  for (int l = 0; l < nb_ghost_bubbles; l++)
    {
      const int ighost = interfaces.ghost_compo_converter(l);
      const int ibulle_reelle = decoder_numero_bulle(-ighost);
      ghost_to_real_bubble(l) = ibulle_reelle;
    }

  ArrOfDouble bubbles_volume;
  interfaces.calculer_volume_bulles(bubbles_volume, bubbles_barycentre);

  /*
   * Considered a constant grid spacing
   */
  const IJK_Splitting& splitting =eulerian_compo_connex.get_splitting();
  const IJK_Grid_Geometry& geometry = splitting.get_grid_geometry();
  double dx = geometry.get_constant_delta(DIRECTION_I);
  double dy = geometry.get_constant_delta(DIRECTION_J);
  double dz = geometry.get_constant_delta(DIRECTION_K);
  double delta_xyz[3] = {dx, dy, dz};

  /*
   * Look for a larger bounding box (3D)
   */
  double geom_origin_x = geometry.get_origin(DIRECTION_I);
  double geom_origin_y = geometry.get_origin(DIRECTION_J);
  double geom_origin_z = geometry.get_origin(DIRECTION_K);
  double origin_x = geom_origin_x + splitting.get_offset_local(DIRECTION_I) * dx;
  double origin_y = geom_origin_y + splitting.get_offset_local(DIRECTION_J) * dy;
  double origin_z = geom_origin_z + splitting.get_offset_local(DIRECTION_K) * dz;
  double geom_origin[3] = {geom_origin_x, geom_origin_y, geom_origin_z};
  double origin[3] = {origin_x, origin_y, origin_z};
  //
  // DoubleTab min_max_larger_box(nb_bubbles, 3, 2);
  min_max_larger_box.resize(nb_bubbles, 3, 2);
  DoubleTab min_max_larger_box_absolute(nb_bubbles, 3, 2);
  for (int ibubble = 0; ibubble < nb_bubbles; ibubble++)
    {
      for (int dir = 0; dir < 3; dir++)
        {
          min_max_larger_box(ibubble, dir, 0) = origin[dir] + trunc((bounding_box(ibubble, dir, 0) - geom_origin[dir]) / delta_xyz[dir])
                                                * delta_xyz[dir];
          min_max_larger_box_absolute(ibubble, dir, 0) = min_max_larger_box(ibubble, dir, 0) - bubbles_barycentre(ibubble, dir);
        }
      for (int dir = 0; dir < 3; dir++)
        {
          min_max_larger_box(ibubble, dir, 1) = origin[dir] + trunc((bounding_box(ibubble, dir, 1) - geom_origin[dir] + delta_xyz[dir]) / delta_xyz[dir]) * delta_xyz[dir];
          min_max_larger_box_absolute(ibubble, dir, 1) = min_max_larger_box(ibubble, dir, 1) - bubbles_barycentre(ibubble, dir);
        }
    }
  /*
   * FT fields
   */
  const int nk = eulerian_compo_connex.nk();
  const int nj = eulerian_compo_connex.nj();
  const int ni = eulerian_compo_connex.ni();
  // const IJK_Field_double& indic = interfaces.I_ft();
  const IJK_Field_double& indic = interfaces.In_ft();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        for (int ibubble = 0; ibubble < (nb_bubbles + nb_ghost_bubbles); ibubble++)
          {
            const double cell_pos_x = origin_x + (i + 0.5) * delta_xyz[0];
            const double cell_pos_y = origin_y + (j + 0.5) * delta_xyz[1];
            const double cell_pos_z = origin_z + (k + 0.5) * delta_xyz[2];
            double cell_pos[3] = {cell_pos_x, cell_pos_y, cell_pos_z};
            const double chi_l = indic(i,j,k);
            int cell_pos_bool = true;
            int bubble_index;
            for (int dir = 0; dir < 3; dir++)
              {
                double min_box;
                double max_box;
                if (ibubble < nb_bubbles)
                  {
                    bubble_index = ibubble;
                    min_box = min_max_larger_box(bubble_index, dir, 0);
                    max_box = min_max_larger_box(bubble_index, dir, 1);
                  }
                else
                  {
                    bubble_index = ghost_to_real_bubble(ibubble - nb_bubbles);
                    min_box = min_max_larger_box_absolute(bubble_index, dir, 0) + bubbles_barycentre(ibubble, dir);
                    max_box = min_max_larger_box_absolute(bubble_index, dir, 1) + bubbles_barycentre(ibubble, dir);
                  }
                cell_pos_bool = (cell_pos_bool && cell_pos[dir] > min_box && cell_pos[dir] < max_box);
              }
            if (fabs(chi_l) < LIQUID_INDICATOR_TEST)
              {
                if (cell_pos_bool)
                  {
                    eulerian_compo_connex(i,j,k) = bubble_index;
                    eulerian_compo_connex_ghost(i,j,k) = ibubble;
                  }
              }
          }
}

void compute_interfacial_compo_fill_compo(const IJK_Interfaces& interfaces, IJK_Field_double& eulerian_compo_connex)
{

}

void compute_rising_velocity(const FixedVector<IJK_Field_double, 3>& velocity, const IJK_Interfaces& interfaces,
                             const IJK_Field_int& eulerian_compo_connex_ns, const int& gravity_dir,
                             ArrOfDouble& rising_velocities, DoubleTab& rising_vectors,
                             Vecteur3& liquid_velocity)
{
  /*
   * Constant cell volume
   */
  const int nk = eulerian_compo_connex_ns.nk();
  const int nj = eulerian_compo_connex_ns.nj();
  const int ni = eulerian_compo_connex_ns.ni();
  const IJK_Field_double& indic = interfaces.I();
  int nb_bubbles = interfaces.get_nb_bulles_reelles();
  DoubleTab sum_indicator(nb_bubbles);
  DoubleTab sum_velocity_x_indicator(nb_bubbles);
  DoubleTab sum_velocity_y_indicator(nb_bubbles);
  DoubleTab sum_velocity_z_indicator(nb_bubbles);
  liquid_velocity = 0.;
  double sum_indicator_liquid = 0.;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          const double chi_v = (1. - indic(i,j,k));
          const double vel_x = velocity[0](i,j,k);
          const double vel_y = velocity[1](i,j,k);
          const double vel_z = velocity[2](i,j,k);
          int compo_connex = eulerian_compo_connex_ns(i,j,k);
          if (compo_connex >= 0)
            {
              sum_indicator(compo_connex) += chi_v;
              sum_velocity_x_indicator(compo_connex) += chi_v * vel_x;
              sum_velocity_y_indicator(compo_connex) += chi_v * vel_y;
              sum_velocity_z_indicator(compo_connex) += chi_v * vel_z;
            }
          if (chi_v > VAPOUR_INDICATOR_TEST)
            {
              Vecteur3 liquid_velocity_local = {vel_x, vel_y, vel_z};
              liquid_velocity_local *= (1-chi_v);
              sum_indicator_liquid += (1-chi_v);
              liquid_velocity += liquid_velocity_local;
            }
        }
  for (int ibubble = 0; ibubble < nb_bubbles; ibubble++)
    {
      sum_indicator(ibubble) = Process::mp_sum(sum_indicator(ibubble));
      sum_velocity_x_indicator(ibubble) = Process::mp_sum(sum_velocity_x_indicator(ibubble));
      sum_velocity_y_indicator(ibubble) = Process::mp_sum(sum_velocity_y_indicator(ibubble));
      sum_velocity_z_indicator(ibubble) = Process::mp_sum(sum_velocity_z_indicator(ibubble));
      sum_velocity_x_indicator(ibubble) /= sum_indicator(ibubble);
      sum_velocity_y_indicator(ibubble) /= sum_indicator(ibubble);
      sum_velocity_z_indicator(ibubble) /= sum_indicator(ibubble);
      rising_velocities(ibubble) = sqrt( sum_velocity_x_indicator(ibubble) * sum_velocity_x_indicator(ibubble)
                                         + sum_velocity_y_indicator(ibubble) * sum_velocity_y_indicator(ibubble)
                                         + sum_velocity_z_indicator(ibubble) * sum_velocity_z_indicator(ibubble));
      if (rising_velocities(ibubble) > DMINFLOAT)
        {
          rising_vectors(ibubble, 0) = sum_velocity_x_indicator(ibubble) / rising_velocities(ibubble);
          rising_vectors(ibubble, 1) = sum_velocity_y_indicator(ibubble) / rising_velocities(ibubble);
          rising_vectors(ibubble, 2) = sum_velocity_z_indicator(ibubble) / rising_velocities(ibubble);
        }
      else
        {
          assert(gravity_dir >=0);
          rising_vectors(ibubble, gravity_dir) = 1.;
        }
    }
  sum_indicator_liquid = Process::mp_sum(sum_indicator_liquid);
  double liquid_velocity_x = liquid_velocity[0];
  double liquid_velocity_y = liquid_velocity[1];
  double liquid_velocity_z = liquid_velocity[2];
  liquid_velocity_x = Process::mp_sum(liquid_velocity_x);
  liquid_velocity_y = Process::mp_sum(liquid_velocity_y);
  liquid_velocity_z = Process::mp_sum(liquid_velocity_z);
  liquid_velocity = {liquid_velocity_x, liquid_velocity_y, liquid_velocity_z};
  liquid_velocity *= (1 / (1e-30 + sum_indicator_liquid));
}

void fill_rising_velocity_double(const IJK_Field_double * eulerian_compo_connex_ns,
                                 const ArrOfDouble& rising_velocities,
                                 IJK_Field_double& eulerian_rising_velocity)
{
  const int nk = (*eulerian_compo_connex_ns).nk();
  const int nj = (*eulerian_compo_connex_ns).nj();
  const int ni = (*eulerian_compo_connex_ns).ni();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          double compo_connex = (*eulerian_compo_connex_ns)(i,j,k);
          int int_compo_connex = (int) compo_connex;
          if (int_compo_connex >= 0)
            {
              eulerian_rising_velocity(i,j,k) = rising_velocities(int_compo_connex);
            }
        }
}

void fill_rising_velocity_int(const IJK_Field_int& eulerian_compo_connex_ns,
                              const ArrOfDouble& rising_velocities,
                              IJK_Field_double& eulerian_rising_velocity)
{
  const int nk = eulerian_compo_connex_ns.nk();
  const int nj = eulerian_compo_connex_ns.nj();
  const int ni = eulerian_compo_connex_ns.ni();
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          int int_compo_connex = eulerian_compo_connex_ns(i,j,k);
          if (int_compo_connex >= 0)
            {
              eulerian_rising_velocity(i,j,k) = rising_velocities(int_compo_connex);
            }
        }
}
