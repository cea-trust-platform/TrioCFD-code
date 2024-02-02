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
// File      : IJK_Thermals.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermals.h>
#include <IJK_FT.h>
#include <IJK_switch_FT.h>
#include <IJK_switch.h>

Implemente_instanciable( IJK_Thermals, "IJK_Thermals", LIST(IJK_Thermal) ) ;

IJK_Thermals::IJK_Thermals(const IJK_FT_double& ijk_ft) : IJK_Thermals()
{
  ref_ijk_ft_ = ijk_ft;
}

Sortie& IJK_Thermals::printOn( Sortie& os ) const
{
  return os;
}

Entree& IJK_Thermals::readOn( Entree& is )
{
  if ((*this).size())
    vide();
  LIST(IJK_Thermal)::readOn(is);
  return is;
}

void IJK_Thermals::set_fichier_reprise(const char *lataname)
{
  for (auto& itr : *this)
    itr.set_fichier_reprise(lataname);
}

const Nom& IJK_Thermals::get_fichier_reprise()
{
  return (*this)[0].get_fichier_reprise();
}

void IJK_Thermals::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  associer_post(ijk_ft.get_post());
  associer_interface_intersections(ijk_ft.itfce().get_intersection_ijk_cell(), ijk_ft.itfce().get_intersection_ijk_face());
  for (auto& itr : *this)
    {
      itr.associer(ijk_ft);
      itr.associer_ghost_fluid_fields(ghost_fluid_fields_);
    }
  ghost_fluid_fields_.associer(ijk_ft);
  retrieve_ghost_fluid_params();
}

void IJK_Thermals::retrieve_ghost_fluid_params()
{
  int compute_distance = 1;
  int compute_curvature = 1;
  int n_iter_distance = 6;
  IJK_Field_local_double boundary_flux_kmin;
  IJK_Field_local_double boundary_flux_kmax;
  (*this)[0].get_boundary_fluxes(boundary_flux_kmin, boundary_flux_kmax);
  for (auto& itr : *this)
    itr.retrieve_ghost_fluid_params(compute_distance,
                                    compute_curvature,
                                    n_iter_distance);
  ghost_fluid_fields_.retrieve_ghost_fluid_params(compute_distance,
                                                  compute_curvature,
                                                  n_iter_distance,
                                                  boundary_flux_kmin,
                                                  boundary_flux_kmax);
}

void IJK_Thermals::associer_post(const IJK_FT_Post& ijk_ft_post)
{
  ref_ijk_ft_post_ = ijk_ft_post;
  for (auto& itr : *this)
    itr.associer_post(ijk_ft_post);
}

void IJK_Thermals::associer_switch(const Switch_FT_double& ijk_ft_switch)
{
  ref_ijk_ft_switch_ = ijk_ft_switch;
  for (auto& itr : *this)
    itr.associer_switch(ref_ijk_ft_switch_);
}

void IJK_Thermals::associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                    const Intersection_Interface_ijk_face& intersection_ijk_face)
{
  ref_intersection_ijk_cell_ = intersection_ijk_cell;
  ref_intersection_ijk_face_ = intersection_ijk_face;
  for (auto& itr : *this)
    itr.associer_interface_intersections(intersection_ijk_cell, intersection_ijk_face);
}

double IJK_Thermals::get_modified_time()
{
  double modified_time = ref_ijk_ft_->get_current_time();
  for (auto& itr : *this)
    modified_time = std::max(modified_time, itr.get_modified_time());
  return modified_time;
}

void IJK_Thermals::get_rising_velocities_parameters(int& compute_rising_velocities, int& fill_rising_velocities)
{
  for (auto& itr : *this)
    itr.get_rising_velocities_parameters(compute_rising_velocities, fill_rising_velocities);
}

void IJK_Thermals::sauvegarder_temperature(Nom& lata_name,
                                           const int& stop)
{
  int idth = 0;
  for (auto& itr : *this)
    itr.sauvegarder_temperature(lata_name, idth, stop);
  idth++;
}

void IJK_Thermals::sauvegarder_thermals(SFichier& fichier)
{
  int flag_list_not_empty_th = 0;
  if ((*this).size() > 0)
    {
      fichier << " thermals {\n" ;
      flag_list_not_empty_th = 1;
    }
  for(auto itr = (*this).begin(); itr != (*this).end(); )
    {
      fichier << *itr ;
      ++itr;
      if (itr != (*this).end())
        fichier << ", \n" ;
      else
        fichier << "\n" ;
    }
  if (flag_list_not_empty_th)
    fichier << " } \n" ;
}

void IJK_Thermals::compute_timestep(double& dt_thermals, const double dxmin)
{
  for (auto& itr : (*this))
    {
      const double dt_th = itr.compute_timestep(dt_thermals, dxmin);
      // We take the most restrictive of all thermal problems and use it for all:
      dt_thermals = std::min(dt_thermals, dt_th);
    }
}

void IJK_Thermals::initialize(const IJK_Splitting& splitting, int& nalloc)
{
  ghost_fluid_fields_.initialize(nalloc, splitting);
  int idth =0;
  Nom thermal_outputs_rank_base = Nom("thermal_outputs_rank_");
  const int max_digit = 3;
  for (auto& itr : (*this))
    {
      nalloc += itr.initialize(splitting, idth);
      if (!ref_ijk_ft_->disable_diphasique())
        itr.update_thermal_properties();
      const int max_rank_digit = idth < 1 ? 1 : (int) (log10(idth) + 1);
      thermal_rank_folder_.add(thermal_outputs_rank_base
                               + Nom(std::string(max_digit - max_rank_digit, '0')) + Nom(idth));
      idth++;
    }
  if (!ref_ijk_ft_->get_disable_diphasique())
    {
      overall_bubbles_quantities_folder_ = Nom("overall_bubbles_quantities");
      interfacial_quantities_thermal_probes_folder_ = Nom("interfacial_quantities_thermal_probes");
      local_quantities_thermal_probes_folder_ = Nom("local_quantities_thermal_probes");
      local_quantities_thermal_probes_time_index_folder_ = Nom("local_quantities_thermal_probes_time_index_");
    }
  for (auto& itr : (*this))
    {
      lata_step_reprise_.push_back(itr.valeur().get_latastep_reprise());
      lata_step_reprise_ini_.push_back(itr.valeur().get_latastep_reprise_ini());
    }
}

void IJK_Thermals::recompute_temperature_init()
{
  for (auto& itr : (*this))
    itr.recompute_temperature_init();
}

int IJK_Thermals::size_thermal_problem(Nom thermal_problem)
{
  int size=0;
  for (auto& itr : (*this))
    {
      if (thermal_problem == itr.get_thermal_problem_type())
        size++;
    }
  return size;
}

void IJK_Thermals::update_thermal_properties()
{
  for (auto& itr : (*this))
    itr.update_thermal_properties();
}

void IJK_Thermals::euler_time_step(const double timestep)
{
  for (auto& itr : (*this))
    itr.euler_time_step(timestep);
  ghost_fluid_fields_.enforce_distance_curvature_values_for_post_processings();
}

void IJK_Thermals::euler_rustine_step(const double timestep)
{
  for (auto& itr : (*this))
    if (itr.get_thermal_problem_type() == Nom("onefluid"))
      {
        itr.update_thermal_properties();
        if (itr.get_conserv_energy_global())
          {
            const double dE = itr.get_E0() - itr.compute_global_energy();
            itr.euler_rustine_step(timestep, dE);
          }
      }
}

void IJK_Thermals::rk3_sub_step(const int rk_step, const double total_timestep, const double time)
{
  for (auto& itr : (*this))
    {
      int thermal_rank = itr.get_thermal_rank();
      switch (thermal_rank)
        {
        case 0:
          Cerr << "RK3 Time scheme is not implemented yet with" << itr.get_thermal_words()[thermal_rank] << finl;
          break;
        case 1:
          Cerr << "RK3 Time scheme is not implemented yet with" << itr.get_thermal_words()[thermal_rank] << finl;
          break;
        case 2:
          itr.rk3_sub_step(rk_step, total_timestep, time);
          Cerr << "RK3 Time scheme is implemented with" << itr.get_thermal_words()[thermal_rank] << finl;
          break;
        case 3:
          Cerr << "RK3 Time scheme is not implemented  with" << itr.get_thermal_words()[thermal_rank] << finl;
          break;
        default:
          Process::exit();
        }
    }
}

void IJK_Thermals::rk3_rustine_sub_step(const int rk_step, const double total_timestep,
                                        const double fractionnal_timestep, const double time)
{
  for (auto& itr : (*this))
    if (itr.get_thermal_problem_type() == Nom("onefluid") )
      {
        itr.update_thermal_properties();
        if (itr.get_conserv_energy_global())
          {
            const double dE = itr.get_E0() - itr.compute_global_energy();
            itr.rk3_rustine_sub_step(rk_step, total_timestep, fractionnal_timestep, time, dE);
          }
      }
}

void IJK_Thermals::posttraiter_tous_champs_thermal(Motcles& liste_post_instantanes_)
{
  int idx_th = 0;
  for (auto &itr : (*this))
    {
      itr.posttraiter_tous_champs_thermal(liste_post_instantanes_, idx_th);
      ++idx_th;
    }
}

void IJK_Thermals::posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes,
                                                          const char *lata_name,
                                                          const int latastep,
                                                          const double current_time,
                                                          int& n)
{
  Cerr << "Post-process Eulerian fields related to the temperature resolution" << finl;
  int idx_th = 0;
  for (auto &itr : (*this))
    {
      int nb = itr.posttraiter_champs_instantanes_thermal(liste_post_instantanes,
                                                          lata_name,
                                                          latastep,
                                                          current_time,
                                                          idx_th);
      // Interfacial thermal fields :
      if (!(ref_ijk_ft_->disable_diphasique()))
        nb += itr.posttraiter_champs_instantanes_thermal_interface(liste_post_instantanes,
                                                                   lata_name,
                                                                   latastep,
                                                                   current_time,
                                                                   idx_th);
      if (idx_th == 0)
        n -= nb; // On compte comme "un" tous les CHAMPS_N (ou N est la longueur de la liste)
      ++idx_th;
    }
}

int IJK_Thermals::init_switch_thermals(const IJK_Splitting& splitting)
{
  int nb_allocated_arrays=0;
  int idx =0;
  for (auto& itr : (*this))
    {
      Cout << "Reading the old temperature field from " << Nom(itr.get_fichier_sauvegarde())
           << " to fill the (*this) field."<< finl;
      nb_allocated_arrays += itr.initialize_switch(splitting, idx);
      idx++;
    }
  return nb_allocated_arrays;
}

void IJK_Thermals::prepare_thermals(const char *lata_name)
{
  for (auto& itr : (*this))
    itr.set_fichier_sauvegarde(lata_name);
}

void IJK_Thermals::ecrire_fichier_reprise(SFichier& fichier, const char *lata_name)
{
  Cerr << "  potentially saving temperature fields... " << finl;
  int flag_list_not_empty = 0;
  if ((int) (*this).size() > 0)
    {
      fichier << " thermals {\n" ;
      flag_list_not_empty = 1;
    }
  int idx =0;
  for (auto itr = (*this).begin(); itr != (*this).end(); )
    {
      (*itr).set_fichier_sauvegarde(lata_name);
      fichier << *itr ;
      ++itr;
      if (itr != (*this).end() )
        fichier << ", \n" ;
      else
        fichier << "\n" ;
      Cerr << "  end of temperature field #" << idx << "... " << finl;
      ++idx;
    }
  if (flag_list_not_empty)
    fichier << " } \n" ;
}

int IJK_Thermals::ghost_fluid_flag()
{
  int ghost_fluid = 0;
  for (auto& itr : (*this))
    {
      ghost_fluid = itr.get_ghost_fluid_flag();
      if (ghost_fluid)
        return ghost_fluid;
    }
  return ghost_fluid;
}

void IJK_Thermals::compute_ghost_cell_numbers_for_subproblems(const IJK_Splitting& splitting, int ghost_init)
{
  for (auto& itr : (*this))
    itr.compute_ghost_cell_numbers_for_subproblems(splitting, ghost_init);
}

int IJK_Thermals::get_probes_ghost_cells(int ghost_init)
{
  int ghost_cells = ghost_init;
  for (auto& itr : (*this))
    {
      const int itr_ghost_cells = itr.get_probes_ghost_cells();
      if (itr_ghost_cells > ghost_cells)
        ghost_cells = itr_ghost_cells;
    }
  return ghost_cells;
}

void IJK_Thermals::update_intersections()
{
  for (auto& itr : (*this))
    itr.update_intersections();
}

void IJK_Thermals::clean_ijk_intersections()
{
  for (auto& itr : (*this))
    itr.clean_ijk_intersections();
}

void IJK_Thermals::compute_eulerian_distance()
{
  ghost_fluid_fields_.compute_eulerian_distance();
}

void IJK_Thermals::compute_eulerian_curvature()
{
  ghost_fluid_fields_.compute_eulerian_curvature();
}

void IJK_Thermals::compute_eulerian_curvature_from_interface()
{
  ghost_fluid_fields_.compute_eulerian_curvature_from_interface();
}

void IJK_Thermals::compute_eulerian_distance_curvature()
{
  compute_eulerian_distance();
  compute_eulerian_curvature_from_interface();
}

int IJK_Thermals::get_disable_post_processing_probes_out_files() const
{
  int disable_post_processing_probes_out_files = 1;
  for (auto& itr : (*this))
    disable_post_processing_probes_out_files = (disable_post_processing_probes_out_files && itr.get_disable_post_processing_probes_out_files());
  return disable_post_processing_probes_out_files;
}

void IJK_Thermals::set_latastep_reprise(const bool stop)
{
  if (stop)
    for (auto& itr : (*this))
      itr.valeur().set_latastep_reprise(ref_ijk_ft_->get_tstep() + 1);
}

void IJK_Thermals::thermal_subresolution_outputs(const int& dt_post_thermals_probes)
{
  const int disable_post_processing_probes_out_files = get_disable_post_processing_probes_out_files();
  if (!disable_post_processing_probes_out_files && post_pro_first_call_)
    {
      if(!ini_folder_out_files_)
        {
          create_folders_for_probes();
          ini_folder_out_files_ = 1;
        }
      int rank = 0;
      for (auto& itr : (*this))
        {
          const int last_time = ref_ijk_ft_->get_tstep() + lata_step_reprise_ini_[rank];
          const int max_digit_time = 8;
          const int nb_digit_tstep = last_time < 1 ? 1 : (int) (log10(last_time) + 1);
          Nom local_quantities_thermal_probes_time_index_folder = thermal_rank_folder_[rank] + "/"
                                                                  + local_quantities_thermal_probes_folder_ + "/"
                                                                  + local_quantities_thermal_probes_time_index_folder_
                                                                  + Nom(std::string(max_digit_time - nb_digit_tstep, '0')) + Nom(last_time);
          Nom overall_bubbles_quantities = thermal_rank_folder_[rank] + "/" + overall_bubbles_quantities_folder_;
          Nom interfacial_quantities_thermal_probes = thermal_rank_folder_[rank] + "/" + interfacial_quantities_thermal_probes_folder_;

          create_folders(local_quantities_thermal_probes_time_index_folder);


          itr.thermal_subresolution_outputs(interfacial_quantities_thermal_probes,
                                            overall_bubbles_quantities,
                                            local_quantities_thermal_probes_time_index_folder);
          // .sauv written before the post-processing on probes
          int latastep_reprise = lata_step_reprise_ini_[rank] + ref_ijk_ft_->get_tstep() + 2;
          const int nb_dt_max = ref_ijk_ft_->get_nb_timesteps();
          if ((ref_ijk_ft_->get_tstep() + dt_post_thermals_probes) >= nb_dt_max)
            latastep_reprise = nb_dt_max + 1;
          itr.valeur().set_latastep_reprise(latastep_reprise);
          rank++;
        }
    }
  post_pro_first_call_++;
}

void IJK_Thermals::create_folders_for_probes()
{
  for (int idth = 0; idth < (*this).size(); idth++)
    {
      create_folders(thermal_rank_folder_[idth]);
      create_folders(thermal_rank_folder_[idth] + "/" + overall_bubbles_quantities_folder_);
      create_folders(thermal_rank_folder_[idth] + "/" + interfacial_quantities_thermal_probes_folder_);
      create_folders(thermal_rank_folder_[idth] + "/" + local_quantities_thermal_probes_folder_);
    }
}

void IJK_Thermals::create_folders(Nom folder_name_base)
{
  if (Process::je_suis_maitre())
    {
      Nom spacing = " ";
      Nom folder_name = "mkdir -p";
      folder_name = folder_name + spacing + folder_name_base.getString(); // donothing";
      Cerr << "Shell command executed: " << folder_name << finl;
      int error = system(folder_name);
      assert(!error);
      if (error)
        Process::exit();
    }
}

void IJK_Thermals::set_first_step_thermals_post(int& first_step_thermals_post)
{
  first_step_thermals_post = 0;
  for (int idth = 0; idth < (*this).size(); idth++)
    first_step_thermals_post = (first_step_thermals_post || (*this)[idth].get_first_step_thermals_post());
}

void IJK_Thermals::set_temperature_ini()
{
  for (auto& itr : (*this))
    itr.compute_temperature_init();
}

void IJK_Thermals::recompute_interface_smoothing()
{
  set_temperature_ini();
  set_post_pro_first_call();
}

void IJK_Thermals::compute_new_thermal_field(Switch_FT_double& switch_double_ft,
                                             const IJK_Splitting& new_mesh,
                                             const Nom& lata_name,
                                             DoubleTab& coeff_i,
                                             IntTab Indice_i,
                                             DoubleTab& coeff_j,
                                             IntTab Indice_j,
                                             DoubleTab& coeff_k,
                                             IntTab Indice_k)
{
  IJK_Field_double new_thermal_field;
  if ((*this).size() > 0)
    {
      switch_double_ft.calculer_coords_elem();
      switch_double_ft.calculer_coeff(coeff_i,Indice_i,coeff_j,Indice_j,coeff_k,Indice_k);
      new_thermal_field.allocate(new_mesh /* it is in fact a splitting */, IJK_Splitting::ELEM, 0);
    }
  int idth = 0;
  for (auto& itr : (*this))
    {
      switch_double_ft.switch_scalar_field(itr.get_temperature(),
                                           new_thermal_field,
                                           coeff_i, Indice_i,
                                           coeff_j ,Indice_j,
                                           coeff_k ,Indice_k);

      Cout << "Writing " << Nom("TEMPERATURE_") + Nom(idth) << " into " << lata_name << finl;
      dumplata_scalar(lata_name, Nom("TEMPERATURE_") + Nom(idth), new_thermal_field, 0 /*we store a 0 */);
      ++idth;
    }
}
