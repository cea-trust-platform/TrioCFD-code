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
// File      : IJK_Composantes_Connex.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Composantes_Connex.h>
#include <IJK_FT.h>
#include <IJK_Interfaces.h>
#include <IJK_Bubble_tools.h>

Implemente_instanciable( IJK_Composantes_Connex, "IJK_Composantes_Connex", Objet_U ) ;

static int decoder_numero_bulle(const int code)
{
  const int num_bulle = code >> 6;
  return num_bulle;
}

Sortie& IJK_Composantes_Connex::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_Composantes_Connex::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

int IJK_Composantes_Connex::initialize(const IJK_Splitting& splitting,
                                       const IJK_Interfaces& interfaces,
                                       const bool is_switch)
{
  int nalloc = 0;
  if (!is_switch)
    {
      interfaces_ = &interfaces;

      if (Process::nproc() == 1)
        {
          eulerian_compo_connex_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
          nalloc += 1;
          eulerian_compo_connex_ft_.data() = -1.;
          eulerian_compo_connex_ft_.echange_espace_virtuel(eulerian_compo_connex_ft_.ghost());

          eulerian_compo_connex_ghost_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
          nalloc += 1;
          eulerian_compo_connex_ghost_ft_.data() = -1.;
          eulerian_compo_connex_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_ghost_ft_.ghost());

          eulerian_compo_connex_ns_.allocate(splitting, IJK_Splitting::ELEM, 0);
          nalloc += 1;
          eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());

          eulerian_compo_connex_ghost_ns_.allocate(splitting, IJK_Splitting::ELEM, 0);
          nalloc += 1;
          eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
        }
      eulerian_compo_connex_from_interface_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 0);
      nalloc += 1;
      eulerian_compo_connex_from_interface_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ft_.ghost());

      eulerian_compo_connex_from_interface_ns_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 1;
      eulerian_compo_connex_from_interface_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ns_.ghost());

      eulerian_compo_connex_from_interface_ghost_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 0);
      nalloc += 1;
      eulerian_compo_connex_from_interface_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ft_.ghost());

      eulerian_compo_connex_from_interface_ghost_ns_.allocate(splitting, IJK_Splitting::ELEM, 0);
      nalloc += 1;
      eulerian_compo_connex_from_interface_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ns_.ghost());

      eulerian_compo_connex_from_interface_int_ns_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 1;
      eulerian_compo_connex_from_interface_int_ns_.data() = -1;
      eulerian_compo_connex_from_interface_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());

      eulerian_compo_connex_from_interface_ghost_int_ns_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 1;
      eulerian_compo_connex_from_interface_ghost_int_ns_.data() = -1;
      eulerian_compo_connex_from_interface_ghost_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_int_ns_.ghost());

      eulerian_compo_connex_valid_compo_field_.allocate(splitting, IJK_Splitting::ELEM, 1);
      nalloc += 1;
      eulerian_compo_connex_valid_compo_field_.data() = 0;
      eulerian_compo_connex_valid_compo_field_.echange_espace_virtuel(eulerian_compo_connex_valid_compo_field_.ghost());


    }
  return nalloc;
}

void IJK_Composantes_Connex::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
}

void IJK_Composantes_Connex::initialise_bubbles_params()
{
  interfaces_->calculer_volume_bulles(bubbles_volume_, bubbles_barycentre_);
}

int IJK_Composantes_Connex::associate_rising_velocities_parameters(const IJK_Splitting& splitting,
                                                                   const int& compute_rising_velocities,
                                                                   const int& fill_rising_velocities)
{
  int nalloc = 0;
  compute_rising_velocities_ = compute_rising_velocities;
  fill_rising_velocities_ = fill_rising_velocities;
  if (fill_rising_velocities_)
    {
      eulerian_rising_velocities_.allocate(splitting, IJK_Splitting::ELEM, 0);
      eulerian_rising_velocities_.data() = 0;
      nalloc += 1;
      eulerian_rising_velocities_.echange_espace_virtuel(eulerian_rising_velocities_.ghost());
    }
  return nalloc;
}

void IJK_Composantes_Connex::compute_bounding_box_fill_compo_connex()
{
  if (Process::nproc() != 1)
    interfaces_->calculer_bounding_box_bulles(bounding_box_);
  else
    {
      compute_bounding_box_fill_compo(ref_ijk_ft_->itfce(),
                                      bounding_box_,
                                      min_max_larger_box_,
                                      eulerian_compo_connex_ft_,
                                      eulerian_compo_connex_ghost_ft_,
                                      bubbles_barycentre_);
      eulerian_compo_connex_ft_.echange_espace_virtuel(eulerian_compo_connex_ft_.ghost());
      eulerian_compo_connex_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_ghost_ft_.ghost());

      eulerian_compo_connex_ns_.data() = -1;
      eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());
      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_compo_connex_ft_, eulerian_compo_connex_ns_);
      eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());

      eulerian_compo_connex_ghost_ns_.data() = -1;
      eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
      ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_compo_connex_ghost_ft_, eulerian_compo_connex_ghost_ns_);
      eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
    }

}

void IJK_Composantes_Connex::compute_compo_connex_from_interface()
{
  interfaces_->calculer_volume_bulles(bubbles_volume_, bubbles_barycentre_);

  fill_mixed_cell_compo();

  const IJK_Splitting& splitting = eulerian_compo_connex_from_interface_int_ns_.get_splitting();
  int neighours_i[6] = NEIGHBOURS_I;
  int neighours_j[6] = NEIGHBOURS_J;
  int neighours_k[6] = NEIGHBOURS_K;

  const int nx = eulerian_compo_connex_from_interface_int_ns_.ni();
  const int ny = eulerian_compo_connex_from_interface_int_ns_.nj();
  const int nz = eulerian_compo_connex_from_interface_int_ns_.nk();
  //ArrOfInt elems_vap;
  ArrOfInt elems_valid;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        if(eulerian_compo_connex_valid_compo_field_(i,j,k))
          elems_valid.append_array(splitting.convert_ijk_cell_to_packed(i,j,k));


  int elems_valid_size = elems_valid.size_array();
  while (elems_valid_size > 0)
    {
      ArrOfInt elems_valid_copy = elems_valid;
      elems_valid.reset();
      for (int elem=0; elem<elems_valid_copy.size_array(); elem++)
        {
          const Int3 num_elem_ijk = splitting.convert_packed_to_ijk_cell(elems_valid_copy[elem]);
          const int i = num_elem_ijk[DIRECTION_I];
          const int j = num_elem_ijk[DIRECTION_J];
          const int k = num_elem_ijk[DIRECTION_K];
          const int num_compo = eulerian_compo_connex_from_interface_int_ns_(i,j,k);
          const int num_compo_ghost = eulerian_compo_connex_from_interface_ghost_int_ns_(i,j,k);
          for (int l = 0; l < 6; l++)
            {
              const int ii = neighours_i[l];
              const int jj = neighours_j[l];
              const int kk = neighours_k[l];
              if((i + ii < 0 || j + jj < 0 || k + kk < 0) || (i + ii >= nx || j + jj >= ny || k + kk >= nz))
                break;
              const int num = eulerian_compo_connex_from_interface_int_ns_(i + ii,j + jj,k + kk);
              const double indic_neighbour =  interfaces_->In()(i + ii,j + jj,k + kk);
              if (num == -1 && indic_neighbour < VAPOUR_INDICATOR_TEST)
                {
                  const int num_elem = splitting.convert_ijk_cell_to_packed(i + ii,j + jj,k + kk);
                  elems_valid.append_array(num_elem);
                  eulerian_compo_connex_from_interface_int_ns_(i + ii,j + jj,k + kk) = num_compo;
                  eulerian_compo_connex_from_interface_ghost_int_ns_(i + ii,j + jj,k + kk) = num_compo_ghost;
                }
            }
        }
      elems_valid_size = elems_valid.size_array();
    }
  eulerian_compo_connex_from_interface_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());
  eulerian_compo_connex_from_interface_ghost_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());
}

void IJK_Composantes_Connex::fill_mixed_cell_compo()
{

  const Domaine_dis_base& mon_dom_dis = interfaces_->get_domaine_dis().valeur();
  const int nb_elem = mon_dom_dis.domaine().nb_elem();
  const Maillage_FT_IJK& maillage = interfaces_->maillage_ft_ijk();

  // Same splitting for the normal vector field
  const IJK_Splitting& splitting_ft = ref_ijk_ft_->get_splitting_ft();

  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  eulerian_compo_connex_from_interface_ft_.data() = -1;
  eulerian_compo_connex_from_interface_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ft_.ghost());
  eulerian_compo_connex_from_interface_ghost_ft_.data() = -1;
  eulerian_compo_connex_from_interface_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ft_.ghost());
  // Loop on the elements
  const ArrOfInt& compo_facettes = maillage.compo_connexe_facettes();
  ArrOfInt compo_per_cell;
  ArrOfInt compo_ghost_per_cell;
  ArrOfInt count_compo_per_cell;
  ArrOfInt count_compo_ghost_per_cell;
  int counter, counter_ghost;
  FixedVector<IJK_Field_double *,2> compo_connex_non_ghost_ghost;
  compo_connex_non_ghost_ghost[0] = &eulerian_compo_connex_from_interface_ft_;
  compo_connex_non_ghost_ghost[1] = &eulerian_compo_connex_from_interface_ghost_ft_;
  const int nbulles_reelles = interfaces_->get_nb_bulles_reelles();
  int l;
  for (int elem = 0; elem < nb_elem; elem++)
    {
      int index = index_elem[elem];
      compo_per_cell.reset();
      compo_ghost_per_cell.reset();
      count_compo_per_cell.reset();
      count_compo_ghost_per_cell.reset();
      // Loop on the facets which cross the element
      while (index >= 0)
        {
          const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
          const int num_facette = data.numero_facette_;
          // const int n = mesh.nb_facettes();
          const int compo_facet = compo_facettes[num_facette];
          const int compo_size_array = compo_per_cell.size_array();
          const int compo_ghost_size_array = compo_ghost_per_cell.size_array();
          int compo_real;
          int compo_ghost;
          if (compo_facet < 0)
            {
              compo_real = decoder_numero_bulle(-compo_facet);
              compo_ghost = compo_real + nbulles_reelles;
            }
          else
            {
              compo_real = compo_facet;
              compo_ghost = compo_facet;
            }
          counter = 0;
          counter_ghost = 0;
          for (l=0; l<compo_size_array; l++)
            {
              if (compo_real == compo_per_cell(l))
                {
                  count_compo_per_cell(l) += 1;
                  break;
                }
              counter++;
            }
          for (l=0; l<compo_ghost_size_array; l++)
            {
              if (compo_ghost == compo_ghost_per_cell(l))
                {
                  count_compo_ghost_per_cell(l) += 1;
                  break;
                }
              counter_ghost++;
            }
          if (counter == compo_size_array)
            {

              compo_per_cell.append_array(compo_real);
              count_compo_per_cell.append_array(1);
            }
          if (counter_ghost == compo_ghost_size_array)
            {

              compo_ghost_per_cell.append_array(compo_ghost);
              count_compo_ghost_per_cell.append_array(1);
            }
          index = data.index_facette_suivante_;
        }
      const int n = compo_per_cell.size_array();
      const int n_ghost = compo_ghost_per_cell.size_array();
      if (n > 0)
        {
          std::vector<int> indices(n);
          for (int j=0; j<n; j++)
            indices[j]=j;
          std::sort(indices.begin(), indices.end(), [&count_compo_per_cell](int i, int j) {return count_compo_per_cell[i] < count_compo_per_cell[j];});
          const int max_compo_per_cell = compo_per_cell(indices[n-1]);
          const Int3 num_elem_ijk = splitting_ft.convert_packed_to_ijk_cell(elem);
          eulerian_compo_connex_from_interface_ft_(num_elem_ijk[DIRECTION_I],num_elem_ijk[DIRECTION_J],num_elem_ijk[DIRECTION_K]) = (double) max_compo_per_cell;
        }
      if (n_ghost > 0)
        {
          std::vector<int> indices(n);
          for (int j=0; j<n; j++)
            indices[j]=j;
          std::sort(indices.begin(), indices.end(), [&count_compo_ghost_per_cell](int i, int j) {return count_compo_ghost_per_cell[i] < count_compo_ghost_per_cell[j];});
          const int max_compo_ghost_per_cell = compo_ghost_per_cell(indices[n-1]);
          const Int3 num_elem_ijk = splitting_ft.convert_packed_to_ijk_cell(elem);
          eulerian_compo_connex_from_interface_ghost_ft_(num_elem_ijk[DIRECTION_I],num_elem_ijk[DIRECTION_J],num_elem_ijk[DIRECTION_K]) = (double) max_compo_ghost_per_cell;
        }
    }
  eulerian_compo_connex_from_interface_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ft_.ghost());
  eulerian_compo_connex_from_interface_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ft_.ghost());
  eulerian_compo_connex_from_interface_ns_.data() = -1;
  eulerian_compo_connex_from_interface_ghost_ns_.data() = -1;
  ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_compo_connex_from_interface_ft_, eulerian_compo_connex_from_interface_ns_);
  ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_compo_connex_from_interface_ghost_ft_, eulerian_compo_connex_from_interface_ghost_ns_);
  eulerian_compo_connex_from_interface_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ns_.ghost());
  eulerian_compo_connex_from_interface_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_ns_.ghost());
  const int nx = eulerian_compo_connex_from_interface_int_ns_.ni();
  const int ny = eulerian_compo_connex_from_interface_int_ns_.nj();
  const int nz = eulerian_compo_connex_from_interface_int_ns_.nk();
  eulerian_compo_connex_valid_compo_field_.data() = 0;
  eulerian_compo_connex_from_interface_int_ns_.data() = -1;
  eulerian_compo_connex_from_interface_ghost_int_ns_.data() = -1;
  for (int k=0; k < nz ; k++)
    for (int j=0; j< ny; j++)
      for (int i=0; i < nx; i++)
        {
          eulerian_compo_connex_from_interface_int_ns_(i,j,k) = (int) eulerian_compo_connex_from_interface_ns_(i,j,k);
          eulerian_compo_connex_from_interface_ghost_int_ns_(i,j,k) = (int) eulerian_compo_connex_from_interface_ghost_ns_(i,j,k);
          if (eulerian_compo_connex_from_interface_int_ns_(i,j,k) >= 0)
            eulerian_compo_connex_valid_compo_field_(i,j,k) = 1;
        }
  eulerian_compo_connex_from_interface_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_int_ns_.ghost());
  eulerian_compo_connex_from_interface_ghost_int_ns_.echange_espace_virtuel(eulerian_compo_connex_from_interface_ghost_int_ns_.ghost());
  eulerian_compo_connex_valid_compo_field_.echange_espace_virtuel(eulerian_compo_connex_valid_compo_field_.ghost());
}

void IJK_Composantes_Connex::compute_rising_velocities()
{
  if (compute_rising_velocities_)
    {
      int nb_bubbles = ref_ijk_ft_->itfce().get_nb_bulles_reelles();
      rising_velocities_ = ArrOfDouble(nb_bubbles);
      rising_vectors_ = DoubleTab(nb_bubbles, 3);
      compute_rising_velocity(ref_ijk_ft_->get_velocity(), ref_ijk_ft_->itfce(),
                              eulerian_compo_connex_from_interface_int_ns_, ref_ijk_ft_->get_direction_gravite(),
                              rising_velocities_, rising_vectors_,
                              liquid_velocity_);
      if (fill_rising_velocities_)
        {
          eulerian_rising_velocities_.data() = 0.;
          eulerian_rising_velocities_.echange_espace_virtuel(eulerian_rising_velocities_.ghost());
          fill_rising_velocity_int(eulerian_compo_connex_from_interface_int_ns_, rising_velocities_, eulerian_rising_velocities_);
        }
    }
  else
    Cerr << "Don't compute the ghost temperature field" << finl;
}
