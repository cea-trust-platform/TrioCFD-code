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
// File      : IJK_Thermal.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Thermal_included
#define IJK_Thermal_included

#include <IJK_Thermal_base.h>
#include <TRUST_Deriv.h>
#include <IJK_Ghost_Fluid_Fields.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Thermal
//
// <Description of class IJK_Thermal>
//
/////////////////////////////////////////////////////////////////////////////

class IJK_FT_double;
class Switch_FT_double;

class IJK_Thermal: public DERIV(IJK_Thermal_base)
{
  Declare_instanciable( IJK_Thermal );

public:

  inline Nom& get_thermal_problem_type() { return thermal_problem_type_; }
  inline int& get_thermal_rank() { return thermal_rank_; }
  inline Motcles& get_thermal_words() { return thermal_words_; }
  inline Motcles& get_thermal_suffix() { return lata_suffix_; }

  inline void associer(const IJK_FT_double &ijk_ft);
  inline void associer_post(const IJK_FT_Post &ijk_ft_post);
  inline void associer_switch(const Switch_FT_double &ijk_ft_switch);
  inline void associer_interface_intersections(const Intersection_Interface_ijk_cell &intersection_ijk_cell, const Intersection_Interface_ijk_face &intersection_ijk_face);

  void posttraiter_tous_champs_thermal(Motcles &liste, const int idx) const;
  int posttraiter_champs_instantanes_thermal(const Motcles &liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx);
  int posttraiter_champs_instantanes_thermal_interface(const Motcles &liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx);
  int posttraiter_champs_instantanes_thermal_interface_ref(const Motcles &liste_post_instantanes, const char *lata_name, const int latastep, const double current_time, const int idx);
  Entree& typer_thermal(Entree &is);
  void thermal_subresolution_outputs(const Nom &interfacial_quantities_thermal_probes, const Nom &overall_bubbles_quantities, const Nom &local_quantities_thermal_probes_time_index_folder);

protected:
  int thermal_rank_;
  Nom thermal_problem_type_;
  Nom prefix_;
  Motcles thermal_words_;
  Motcles lata_suffix_;

  REF(IJK_FT_double) ref_ijk_ft_;
  REF(IJK_FT_Post) ref_ijk_ft_post_;
  REF(Switch_FT_double) ref_ijk_ft_switch_;
  REF(Intersection_Interface_ijk_cell) ref_intersection_ijk_cell_;
  REF(Intersection_Interface_ijk_face) ref_intersection_ijk_face_;

};

inline void IJK_Thermal::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
  valeur().associer(ijk_ft);
}

inline void IJK_Thermal::associer_post(const IJK_FT_Post& ijk_ft_post)
{
  ref_ijk_ft_post_ = ijk_ft_post;
  valeur().associer_post(ijk_ft_post);
}

inline void IJK_Thermal::associer_switch(const Switch_FT_double& ijk_ft_switch)
{
  ref_ijk_ft_switch_ = ijk_ft_switch;
  valeur().associer_switch(ref_ijk_ft_switch_);
}

inline void IJK_Thermal::associer_interface_intersections(const Intersection_Interface_ijk_cell& intersection_ijk_cell, const Intersection_Interface_ijk_face& intersection_ijk_face)
{
  ref_intersection_ijk_cell_ = intersection_ijk_cell;
  ref_intersection_ijk_face_ = intersection_ijk_face;
  valeur().associer_interface_intersections(intersection_ijk_cell, intersection_ijk_face);
}

#endif /* IJK_Thermal_included */
