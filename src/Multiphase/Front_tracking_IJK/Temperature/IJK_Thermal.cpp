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
// File      : IJK_Thermal.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Thermal.h>
#include <IJK_FT.h>

Implemente_instanciable_sans_constructeur( IJK_Thermal, "IJK_Thermal", DERIV(IJK_Thermal_base) ) ;

IJK_Thermal::IJK_Thermal()
{
  thermal_rank_ = 0;
  thermal_problem_type_ = "subresolution";
  prefix_="IJK_Thermal_";
  thermal_words_ = Motcles(4);
  {
    thermal_words_[0] = "subresolution";
    thermal_words_[1] = "multiplesubresolutions";
    thermal_words_[2] = "onefluid";
    thermal_words_[3] = "onefluidenergy";
  }
  lata_suffix_ = Motcles(4);
  {
    lata_suffix_[0] = "SUBRES_";
    lata_suffix_[1] = "MSUBRES_";
    lata_suffix_[2] = "OF_";
    lata_suffix_[3] = "OFE_";
  }
}

Sortie& IJK_Thermal::printOn( Sortie& os ) const
{
  DERIV(IJK_Thermal_base)::printOn( os );
  return os;
}

Entree& IJK_Thermal::readOn( Entree& is )
{
  return typer_thermal(is);
}

Entree& IJK_Thermal::typer_thermal( Entree& is )
{
  Cerr << "Read and Cast IJK_Thermal :" << finl;
  Motcle word;
  is >> word;
  Nom type = "";
  thermal_rank_ = thermal_words_.search(word);
  type += prefix_;
  switch(thermal_rank_)
    {
    case 0 :
      {
        type += "Subresolution";
        break;
      }
    case 1 :
      {
        type += "Multiple_Subresolutions";
        break;
      }
    case 2 :
      {
        type += "Onefluid";
        break;
      }
    case 3 :
      {
        type += "OnefluidEnergy";
        break;
      }
    default :
      {
        Cerr << "ERROR : Thermal problems that are already implemented are:" << finl;
        Cerr << thermal_words_ << finl;
        abort();
      }
    }
  thermal_problem_type_=thermal_words_[thermal_rank_];
  typer(type);
  is >> valeur(); // Call the readOn
  return is;
}

void IJK_Thermal::posttraiter_tous_champs_thermal(Motcles& liste, const int idx) const
{
  /*
   * thermal_words_[0] = "subresolution";
   * thermal_words_[1] = "multiplesubresolutions";
   * thermal_words_[2] = "onefluid";
   * thermal_words_[3] = "onefluidenergy";
   */
  liste.add("TEMPERATURE");
  liste.add("TEMPERATURE_ANA");
  liste.add("ECART_T_ANA");
  liste.add("DIV_LAMBDA_GRAD_T_VOLUME");
  liste.add("GRAD_T");
  //
  liste.add("DISTANCE");
  liste.add("CURVATURE");
  if (thermal_rank_==0 || thermal_rank_==1)
    {
      /*
       * TODO: ADD SOME PARTICULAR FIELDS OR DO SWITCH CASE(thermal_rank)
       */
    }
  else
    {
      /*
       * TODO: CHECK IF GRAD_T0 MUST STILL BE POST-PROCESSED
       */
      liste.add("CP");
      liste.add("LAMBDA");
      //
      liste.add("SOURCE_TEMPERATURE");
      liste.add("TEMPERATURE_ADIM_BULLES");
      liste.add("TEMPERATURE_PHYSIQUE_T");
      liste.add("TEMPERATURE_ADIMENSIONNELLE_THETA");
      liste.add("SOURCE_TEMPERATURE_ANA");
      liste.add("ECART_SOURCE_TEMPERATURE_ANA");
      //
      liste.add("GRAD_T0");
      liste.add("GRAD_T1");
      liste.add("GRAD_T2");
      //
      liste.add("T_RUST");
      liste.add("DIV_RHO_CP_T_V");
    }
}

int IJK_Thermal::posttraiter_champs_instantanes_thermal(const Motcles& liste_post_instantanes,
                                                        const char * lata_name,
                                                        const int latastep,
                                                        const double current_time,
                                                        const int idx)
{
  Cerr << liste_post_instantanes << finl;
  int n = 0; // number of post-processed field
  std::ostringstream oss;
  Motcle lata_suffix = lata_suffix_[thermal_rank_];
  /*
   * TEMPERATURE
   */
  oss << "TEMPERATURE_" << lata_suffix << idx;
  Nom nom_temp(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
    {
      n++, dumplata_scalar(lata_name, nom_temp, get_temperature(), latastep);
    }
  oss.str("");

  oss << "TEMPERATURE_BEFORE_EXTRAP_" << lata_suffix << idx;
  Nom nom_temp_before_extrap(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_BEFORE_EXTRAP")) || (liste_post_instantanes.contient_(nom_temp_before_extrap)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_before_extrap, get_temperature_before_extrapolation(), latastep);
    }
  oss.str("");

  oss << "TEMPERATURE_CELL_NEIGHBOURS_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_CELL_NEIGHBOURS")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours, get_temperature_cell_neighbours(), latastep);
    }
  oss.str("");

  oss << "TEMPERATURE_CELL_NEIGHBOURS_DEBUG_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours_debug(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_CELL_NEIGHBOURS_DEBUG")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours_debug)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours_debug, get_temperature_cell_neighbours_debug(), latastep);
    }
  oss.str("");

  oss << "NEIGHBOURS_CORRECTED_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours_corrected(oss.str().c_str());
  if ((liste_post_instantanes.contient_("NEIGHBOURS_CORRECTED")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours_corrected)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours_corrected, get_cell_neighbours_corrected(), latastep);
    }
  oss.str("");

  oss << "NEIGHBOURS_COLINEARITY_" << lata_suffix << idx;
  Nom nom_temp_cell_neighbours_colinearity(oss.str().c_str());
  if ((liste_post_instantanes.contient_("NEIGHBOURS_COLINEARITY")) || (liste_post_instantanes.contient_(nom_temp_cell_neighbours_colinearity)))
    {
      n++, dumplata_scalar(lata_name, nom_temp_cell_neighbours_colinearity, get_neighbours_temperature_colinearity_weighting(), latastep);
    }
  oss.str("");

  /*
   * TEMPERATURE_ANA
   */
  oss << "TEMPERATURE_ANA_" << lata_suffix << idx;
  Nom nom_ana(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_ANA")) || (liste_post_instantanes.contient_(nom_ana)))
    {
      //set_field_data(itr.temperature_ana_, itr.expression_T_ana_, current_time);
      set_field_T_ana();
      n++, dumplata_scalar(lata_name, nom_ana, get_temperature_ana(), latastep);
    }
  oss.str("");

  /*
   * ECART_T_ANA
   */
  oss << "ECART_T_ANA_" << lata_suffix << idx;
  Nom nom_ecart_ana(oss.str().c_str());
  if ((liste_post_instantanes.contient_("ECART_T_ANA") || liste_post_instantanes.contient_(nom_ecart_ana)))
    {
      calculer_ecart_T_ana();
      n++, dumplata_scalar(lata_name, nom_ecart_ana, get_ecart_t_ana(), latastep);
    }
  oss.str("");

  /*
   * DIV_LAMBDA_GRAD_T_VOLUME
   */
  oss << "DIV_LAMBDA_GRAD_T_VOLUME_" << lata_suffix << idx;
  Nom nom_div_lambda_grad_T(oss.str().c_str());
  if ((liste_post_instantanes.contient_("DIV_LAMBDA_GRAD_T_VOLUME") || liste_post_instantanes.contient_(nom_div_lambda_grad_T)))
    {
      n++, dumplata_scalar(lata_name, nom_div_lambda_grad_T, get_div_lambda_grad_T(), latastep);
    }
  oss.str("");

  /*
   * GRAD_T
   */
  oss << "GRAD_T_" << lata_suffix << idx;
  Nom nom_grad(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T") || (liste_post_instantanes.contient_(nom_grad)))
    {
      const FixedVector<IJK_Field_double, 3>& grad_T = get_grad_T();
      n++, dumplata_vector(lata_name, nom_grad, grad_T[0], grad_T[1], grad_T[2], latastep);
    }
  oss.str("");

  /*
   * DISTANCE
   */
  oss << "DISTANCE_" << lata_suffix << idx;
  Nom nom_distance(oss.str().c_str());
  if ((liste_post_instantanes.contient_("DISTANCE") || liste_post_instantanes.contient_(nom_distance)))
    {
      n++, dumplata_scalar(lata_name, nom_distance, get_eulerian_distance_ns(), latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_X
   */
  oss << "NORMAL_VECTOR_X_" << lata_suffix << idx;
  Nom nom_normal_vector_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_X") || liste_post_instantanes.contient_(nom_normal_vector_x))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_x, get_normal_vector_ns()[0], latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_Y
   */
  oss << "NORMAL_VECTOR_Y_" << lata_suffix << idx;
  Nom nom_normal_vector_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_Y") || liste_post_instantanes.contient_(nom_normal_vector_y))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_y, get_normal_vector_ns()[1], latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_Z
   */
  oss << "NORMAL_VECTOR_Z_" << lata_suffix << idx;
  Nom nom_normal_vector_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_Z") || liste_post_instantanes.contient_(nom_normal_vector_z))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_z, get_normal_vector_ns()[2], latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_X_NORMED
   */
  oss << "NORMAL_VECTOR_X_NORMED_" << lata_suffix << idx;
  Nom nom_normal_vector_x_normed(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_X_NORMED") || liste_post_instantanes.contient_(nom_normal_vector_x_normed))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_x_normed, get_normal_vector_ns_normed()[0], latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_Y_NORMED
   */
  oss << "NORMAL_VECTOR_Y_NORMED_" << lata_suffix << idx;
  Nom nom_normal_vector_y_normed(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_Y_NORMED") || liste_post_instantanes.contient_(nom_normal_vector_y_normed))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_y_normed, get_normal_vector_ns_normed()[1], latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_Z_NORMED
   */
  oss << "NORMAL_VECTOR_Z_NORMED_" << lata_suffix << idx;
  Nom nom_normal_vector_z_normed(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_Z_NORMED") || liste_post_instantanes.contient_(nom_normal_vector_z_normed))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_z_normed, get_normal_vector_ns_normed()[2], latastep);
    }
  oss.str("");

  /*
   * CURVATURE
   */
  oss << "CURVATURE_" << lata_suffix << idx;
  Nom nom_curvature(oss.str().c_str());
  if ((liste_post_instantanes.contient_("CURVATURE") || liste_post_instantanes.contient_(nom_curvature)))
    {
      n++, dumplata_scalar(lata_name, nom_curvature, get_eulerian_curvature_ns(), latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL AREA
   */
  oss << "INTERFACIAL_AREA_" << lata_suffix << idx;
  Nom nom_interfacial_area(oss.str().c_str());
  if ((liste_post_instantanes.contient_("INTERFACIAL_AREA") || liste_post_instantanes.contient_(nom_interfacial_area)))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_area, get_interfacial_area_ns(), latastep);
    }
  oss.str("");

  /*
   * GRAD_T_INTERFACE
   */
  oss << "GRAD_T_INTERFACE_" << lata_suffix << idx;
  Nom nom_grad_T_interface(oss.str().c_str());
  if ((liste_post_instantanes.contient_("GRAD_T_INTERFACE") || liste_post_instantanes.contient_(nom_grad_T_interface)))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_interface, get_grad_T_interface_ns(), latastep);
    }
  oss.str("");

  /*
   * DISTANCE_FT
   */
  oss << "DISTANCE_FT_" << lata_suffix << idx;
  Nom nom_distance_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("DISTANCE_FT") || liste_post_instantanes.contient_(nom_distance_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_distance_ft, get_eulerian_distance_ft(), latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_X_FT
   */
  oss << "NORMAL_VECTOR_X_FT_" << lata_suffix << idx;
  Nom nom_normal_vector_x_ft(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_X_FT") || liste_post_instantanes.contient_("NORMAL_VECTOR_X_FT") || liste_post_instantanes.contient_(nom_normal_vector_x_ft))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_x_ft, get_normal_vector_ft()[0], latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_Y_FT
   */
  oss << "NORMAL_VECTOR_Y_FT_" << lata_suffix << idx;
  Nom nom_normal_vector_y_ft(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_Y_FT") || liste_post_instantanes.contient_("NORMAL_VECTOR_Y_FT") || liste_post_instantanes.contient_(nom_normal_vector_y_ft))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_y_ft, get_normal_vector_ft()[1], latastep);
    }
  oss.str("");

  /*
   * NORMAL_VECTOR_Z_FT
   */
  oss << "NORMAL_VECTOR_Z_FT_" << lata_suffix << idx;
  Nom nom_normal_vector_z_ft(oss.str().c_str());
  if (liste_post_instantanes.contient_("NORMAL_VECTOR_Z_FT") || liste_post_instantanes.contient_("NORMAL_VECTOR_Z_FT") || liste_post_instantanes.contient_(nom_normal_vector_z_ft))
    {
      n++, dumplata_scalar(lata_name, nom_normal_vector_z_ft, get_normal_vector_ft()[2], latastep);
    }
  oss.str("");

  /*
   * CURVATURE_FT
   */
  oss << "CURVATURE_FT_" << lata_suffix << idx;
  Nom nom_curvature_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("CURVATURE_FT") || liste_post_instantanes.contient_(nom_curvature_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_curvature_ft, get_eulerian_curvature_ft(), latastep);
    }
  oss.str("");

  /*
   * INTERFACIAL_AREA_FT
   */
  oss << "INTERFACIAL_AREA_FT_" << lata_suffix << idx;
  Nom nom_interfacial_area_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("INTERFACIAL_AREA_FT") || liste_post_instantanes.contient_(nom_interfacial_area_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_interfacial_area_ft, get_interfacial_area_ft(), latastep);
    }
  oss.str("");

  /*
   * GRAD_T_INTERFACE_FT
   */
  oss << "GRAD_T_INTERFACE_FT_" << lata_suffix << idx;
  Nom nom_grad_T_interface_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("GRAD_T_INTERFACE_FT") || liste_post_instantanes.contient_(nom_grad_T_interface_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_interface_ft, get_grad_T_interface_ft(), latastep);
    }
  oss.str("");

  /*
   * TEMPERATURE_FT
   */
  oss << "TEMPERATURE_FT_" << lata_suffix << idx;
  Nom nom_temperature_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("TEMPERATURE_FT") || liste_post_instantanes.contient_(nom_temperature_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_temperature_ft, get_temperature_ft(), latastep);
    }
  oss.str("");

  /*
   * INDICATOR_FT
   */
  oss << "INDICATOR_FT_" << lata_suffix << idx;
  Nom nom_indicator_ft(oss.str().c_str());
  if ((liste_post_instantanes.contient_("INDICATOR_FT") || liste_post_instantanes.contient_(nom_indicator_ft)))
    {
      n++, dumplata_scalar(lata_name, nom_indicator_ft, ref_ijk_ft_->itfce().I_ft(), latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_X_ELEM
   */
  oss << "GRAD_T_DIR_X_ELEM_" << lata_suffix << idx;
  Nom nom_grad_T_dir_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_X_ELEM") || liste_post_instantanes.contient_("GRAD_T_ELEM")
      || liste_post_instantanes.contient_(nom_grad_T_dir_x))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_x, get_gradient_temperature_elem()[0], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Y_ELEM
   */
  oss << "GRAD_T_DIR_Y_ELEM_" << lata_suffix << idx;
  Nom nom_grad_T_dir_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Y_ELEM_") || liste_post_instantanes.contient_("GRAD_T_ELEM")
      || liste_post_instantanes.contient_(nom_grad_T_dir_y))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_y, get_gradient_temperature_elem()[1], latastep);
    }
  oss.str("");

  /*
   * GRAD_T_DIR_Z_ELEM
   */
  oss << "GRAD_T_DIR_Z_ELEM_" << lata_suffix << idx;
  Nom nom_grad_T_dir_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("GRAD_T_DIR_Z_ELEM_") || liste_post_instantanes.contient_("GRAD_T_ELEM")
      || liste_post_instantanes.contient_(nom_grad_T_dir_z))
    {
      n++, dumplata_scalar(lata_name, nom_grad_T_dir_z, get_gradient_temperature_elem()[2], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_XX_ELEM
   */
  oss << "HESS_T_DIR_XX_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_xx(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_XX_ELEM") || liste_post_instantanes.contient_("HESS_DIAG_T_ELEM")
      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_xx))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xx, get_hessian_diag_temperature_elem()[0], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_YY_ELEM
   */
  oss << "HESS_T_DIR_YY_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_yy(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_YY_ELEM") || liste_post_instantanes.contient_("HESS_DIAG_T_ELEM")
      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_yy))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_yy, get_hessian_diag_temperature_elem()[1], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_ZZ_ELEM
   */
  oss << "HESS_T_DIR_ZZ_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_zz(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_ZZ_ELEM") || liste_post_instantanes.contient_("HESS_DIAG_T_ELEM")
      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_zz))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_zz, get_hessian_diag_temperature_elem()[2], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_XY_YX_ELEM
   */
  oss << "HESS_T_DIR_XY_YX_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_xy_yx(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_XY_YX_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_xy_yx))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xy_yx, get_hessian_cross_temperature_elem()[2], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_XZ_ZX_ELEM
   */
  oss << "HESS_T_DIR_XZ_ZX_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_xz_zx(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_XZ_ZX_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_xz_zx))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xz_zx, get_hessian_cross_temperature_elem()[1], latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_YZ_ELEM
   */
  oss << "HESS_T_DIR_YZ_ZY_ELEM_" << lata_suffix << idx;
  Nom nom_hess_T_dir_yz_zy(oss.str().c_str());
  if (liste_post_instantanes.contient_("HESS_T_DIR_YZ_ZY_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_yz_zy))
    {
      n++, dumplata_scalar(lata_name, nom_hess_T_dir_yz_zy, get_hessian_cross_temperature_elem()[0], latastep);
    }
  oss.str("");

//  /*
//   * HESS_T_DIR_XY_ELEM
//   */
//  oss << "HESS_T_DIR_XY_ELEM_" << lata_suffix << idx;
//  Nom nom_hess_T_dir_xy(oss.str().c_str());
//  if (liste_post_instantanes.contient_("HESS_T_DIR_XY_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
//      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_xy))
//    {
//      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xy, get_hessian_cross_temperature_elem()[2], latastep);
//    }
//  oss.str("");
//
//  /*
//   * HESS_T_DIR_YX_ELEM
//   */
//  oss << "HESS_T_DIR_YX_ELEM_" << lata_suffix << idx;
//  Nom nom_hess_T_dir_yx(oss.str().c_str());
//  if (liste_post_instantanes.contient_("HESS_T_DIR_YX_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
//      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_yx))
//    {
//      n++, dumplata_scalar(lata_name, nom_hess_T_dir_yx, get_hessian_cross_temperature_elem()[2], latastep);
//    }
//  oss.str("");

  /*
   * HESS_T_DIR_XZ_ELEM
   */
//  oss << "HESS_T_DIR_XZ_ELEM_" << lata_suffix << idx;
//  Nom nom_hess_T_dir_xz(oss.str().c_str());
//  if (liste_post_instantanes.contient_("HESS_T_DIR_XZ_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
//      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_xz))
//    {
//      n++, dumplata_scalar(lata_name, nom_hess_T_dir_xz, get_hessian_cross_temperature_elem()[1], latastep);
//    }
//  oss.str("");

  /*
   * HESS_T_DIR_ZX_ELEM
   */
//  oss << "HESS_T_DIR_ZX_ELEM_" << lata_suffix << idx;
//  Nom nom_hess_T_dir_zx(oss.str().c_str());
//  if (liste_post_instantanes.contient_("HESS_T_DIR_ZX_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
//      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_zx))
//    {
//      n++, dumplata_scalar(lata_name, nom_hess_T_dir_zx, get_hessian_cross_temperature_elem()[1], latastep);
//    }
//  oss.str("");

  /*
   * HESS_T_DIR_YZ_ELEM
   */
//  oss << "HESS_T_DIR_YZ_ELEM_" << lata_suffix << idx;
//  Nom nom_hess_T_dir_yz(oss.str().c_str());
//  if (liste_post_instantanes.contient_("HESS_T_DIR_YZ_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
//      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_yz))
//    {
//      n++, dumplata_scalar(lata_name, nom_hess_T_dir_yz, get_hessian_cross_temperature_elem()[0], latastep);
//    }
//  oss.str("");

  /*
   * HESS_T_DIR_YX_ELEM
   */
//  oss << "HESS_T_DIR_ZY_ELEM_" << lata_suffix << idx;
//  Nom nom_hess_T_dir_zy(oss.str().c_str());
//  if (liste_post_instantanes.contient_("HESS_T_DIR_ZY_ELEM") || liste_post_instantanes.contient_("HESS_CROSS_T_ELEM")
//      || liste_post_instantanes.contient_("HESS_T_ELEM") || liste_post_instantanes.contient_(nom_hess_T_dir_zy))
//    {
//      n++, dumplata_scalar(lata_name, nom_hess_T_dir_zy, get_hessian_cross_temperature_elem()[0], latastep);
//    }
//  oss.str("");

  /*
   * BARY_X
   */
  oss << "BARY_X_" << lata_suffix << idx;
  Nom nom_bary_x(oss.str().c_str());
  if (liste_post_instantanes.contient_("BARY") || liste_post_instantanes.contient_("BARY_X") || liste_post_instantanes.contient_(nom_bary_x))
    {
      n++, dumplata_scalar(lata_name, nom_bary_x, get_bary()[0], latastep);
    }
  oss.str("");

  /*
   * BARY_Y
   */
  oss << "BARY_Y_" << lata_suffix << idx;
  Nom nom_bary_y(oss.str().c_str());
  if (liste_post_instantanes.contient_("BARY") || liste_post_instantanes.contient_("BARY_Y") || liste_post_instantanes.contient_(nom_bary_y))
    {
      n++, dumplata_scalar(lata_name, nom_bary_y, get_bary()[1], latastep);
    }
  oss.str("");

  /*
   * BARY_Z
   */
  oss << "BARY_Z_" << lata_suffix << idx;
  Nom nom_bary_z(oss.str().c_str());
  if (liste_post_instantanes.contient_("BARY") || liste_post_instantanes.contient_("BARY_Z") || liste_post_instantanes.contient_(nom_bary_z))

    {
      n++, dumplata_scalar(lata_name, nom_bary_z, get_bary()[2], latastep);
    }
  oss.str("");

  /*
   * EULERIAN_COMPO
   */
  oss << "EULERIAN_COMPO_FT_" << lata_suffix << idx;
  Nom nom_eulerian_compo_ft(oss.str().c_str());
  if (liste_post_instantanes.contient_("EULERIAN_COMPO_FT") || liste_post_instantanes.contient_(nom_eulerian_compo_ft))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_compo_ft, get_eulerian_compo_connex_ft(), latastep);
    }
  oss.str("");

  /*
   * EULERIAN_COMPO
   */
  oss << "EULERIAN_COMPO_" << lata_suffix << idx;
  Nom nom_eulerian_compo(oss.str().c_str());
  if (liste_post_instantanes.contient_("EULERIAN_COMPO") || liste_post_instantanes.contient_(nom_eulerian_compo))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_compo, get_eulerian_compo_connex_ns(), latastep);
    }
  oss.str("");

  /*
   * EULERIAN_COMPO
   */
  oss << "EULERIAN_COMPO_GHOST_FT_" << lata_suffix << idx;
  Nom nom_eulerian_compo_ghost_ft(oss.str().c_str());
  if (liste_post_instantanes.contient_("EULERIAN_COMPO_GHOST_FT") || liste_post_instantanes.contient_(nom_eulerian_compo_ghost_ft))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_compo_ghost_ft, get_eulerian_compo_connex_ghost_ft(), latastep);
    }
  oss.str("");

  /*
   * EULERIAN_COMPO
   */
  oss << "EULERIAN_COMPO_GHOST_" << lata_suffix << idx;
  Nom nom_eulerian_compo_ghost(oss.str().c_str());
  if (liste_post_instantanes.contient_("EULERIAN_COMPO_GHOST") || liste_post_instantanes.contient_(nom_eulerian_compo_ghost))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_compo_ghost, get_eulerian_compo_connex_ghost_ns(), latastep);
    }
  oss.str("");

  /*
   * EULERIAN_COMPO
   */
  oss << "EULERIAN_COMPO_INTERFACE_FT_" << lata_suffix << idx;
  Nom nom_eulerian_compo_from_interface_ft(oss.str().c_str());
  if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_FT") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_ft))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_ft, get_eulerian_compo_connex_from_interface_ft(), latastep);
    }
  oss.str("");

  /*
   * EULERIAN_COMPO
   */
  oss << "EULERIAN_COMPO_INTERFACE_NS_" << lata_suffix << idx;
  Nom nom_eulerian_compo_from_interface_ns(oss.str().c_str());
  if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_NS") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_ns))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_ns, get_eulerian_compo_connex_from_interface_ns(), latastep);
    }
  oss.str("");

  /*
   * EULERIAN_COMPO
   */
  oss << "EULERIAN_COMPO_INTERFACE_INT_" << lata_suffix << idx;
  Nom nom_eulerian_compo_from_interface_int(oss.str().c_str());
  if (liste_post_instantanes.contient_("EULERIAN_COMPO_INTERFACE_INT") || liste_post_instantanes.contient_(nom_eulerian_compo_from_interface_int))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_compo_from_interface_int, get_eulerian_compo_connex_int_from_interface_ns(), latastep);
    }
  oss.str("");

  /*
   * RISING_VELOCITIES
   */
  oss << "RISING_VELOCITIES_" << lata_suffix << idx;
  Nom nom_eulerian_rising_velocities(oss.str().c_str());
  if (liste_post_instantanes.contient_("RISING_VELOCITIES") || liste_post_instantanes.contient_(nom_eulerian_rising_velocities))

    {
      n++, dumplata_scalar(lata_name, nom_eulerian_rising_velocities, get_eulerian_rising_velocities(), latastep);
    }
  oss.str("");

  /*
   * HESS_T_DIR_YZ_ELEM
   */
  if (get_debug())
    {
      oss << "DEBUG_LRS_CELLS_" << lata_suffix << idx;
      Nom nom_debug_lrs_cells(oss.str().c_str());
      if (liste_post_instantanes.contient_("DEBUG_LRS_CELLS") || liste_post_instantanes.contient_(nom_debug_lrs_cells))
        {
          n++, dumplata_scalar(lata_name, nom_debug_lrs_cells, get_debug_lrs_cells(), latastep);
        }
      oss.str("");
    }
  return n;
}

int IJK_Thermal::posttraiter_champs_instantanes_thermal_interface(const Motcles& liste_post_instantanes,
                                                                  const char *lata_name,
                                                                  const int latastep,
                                                                  const double current_time,
                                                                  const int idx)
{
  int n = 0; // nombre de champs postraites
  if (thermal_rank_==0 || thermal_rank_==1)
    {
      /*
       * TODO: COMPUTE INTERFACIAL GRADIENT
       */
      //  Cerr << liste_post_instantanes << finl;
      //  int n = 0; // nombre de champs postraites
      //  std::ostringstream oss;
      //
      //  /*
      //   * NORMAL_T_GRADIENT
      //   */
      //  oss << "NORMAL_T_GRADIENT" << idx;
      //  Nom nom_normal_temp(oss.str().c_str());
      //  if (liste_post_instantanes_.contient_("NORMAL_T_GRADIENT") || (liste_post_instantanes.contient_(nom_normal_temp)))
      //    {
      //
      //    }
      //  oss.str("");
      //
      //  /*
      //   *
      //   */
      //
      //  /*
      //   *
      //   */
    }
  else
    {
      n = posttraiter_champs_instantanes_thermal_interface_ref(liste_post_instantanes, lata_name, latastep, current_time, idx);
    }
  return n;
}

int IJK_Thermal::posttraiter_champs_instantanes_thermal_interface_ref(const Motcles& liste_post_instantanes,
                                                                      const char *lata_name,
                                                                      const int latastep,
                                                                      const double current_time,
                                                                      const int idx)
{
  Cerr << liste_post_instantanes << finl;
  int n = 0; // nombre de champs postraites
  Motcle lata_suffix = lata_suffix_[thermal_rank_];

  std::ostringstream oss;
  oss << "INTERFACE_TEMPERATURE_" << lata_suffix << idx;
  Nom nom_temp(oss.str().c_str());

  std::ostringstream oss2;
  oss2 << "INTERFACE_PHIN_" << lata_suffix << idx;
  Nom nom_phin(oss2.str().c_str());

  if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE"))
      || (liste_post_instantanes.contient_("INTERFACE_PHIN"))
      || (liste_post_instantanes.contient_(nom_temp))
      || (liste_post_instantanes.contient_(nom_phin)))
    {
      //  Computing interfacial temperature at fa7 centre :
      const Maillage_FT_IJK& mesh = ref_ijk_ft_->get_maillage_ft_ijk(); // ref_ijk_ft_post_->interfaces_.maillage_ft_ijk();
      const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
      const int nb_facettes = mesh.nb_facettes();
      ArrOfDouble interfacial_temperature;
      ArrOfDouble interfacial_phin;
      // To transfer the field to FT splitting (because interfaces are there...) !!! NEEDED for compute_interfacial_temperature
      IJK_Field_double& temperature_ft = get_temperature_ft();
      ref_ijk_ft_->redistribute_to_splitting_ft_elem(get_temperature(), temperature_ft);
      temperature_ft.echange_espace_virtuel(temperature_ft.ghost());
      // results are prop to the area :
      //itr.compute_interfacial_temperature(interfacial_temperature, interfacial_phin, itr.get_storage());
      compute_interfacial_temperature2(interfacial_temperature, interfacial_phin);
      for (int fa7 = 0; fa7 < nb_facettes; fa7++)
        {
          const double sf = surface_facettes[fa7];
          interfacial_temperature[fa7] /= sf;
          interfacial_phin[fa7] /= sf;
        }
      if ((liste_post_instantanes.contient_("INTERFACE_TEMPERATURE")) || (liste_post_instantanes.contient_(nom_temp)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_temp, "ELEM", interfacial_temperature, latastep);
      if ((liste_post_instantanes.contient_("INTERFACE_PHIN")) || (liste_post_instantanes.contient_(nom_phin)))
        n++, dumplata_ft_field(lata_name, "INTERFACES", nom_phin, "ELEM", interfacial_phin, latastep);
    }
  oss.str("");
  return n;
}

void IJK_Thermal::thermal_subresolution_outputs(SFichier& fic)
{
  if (thermal_rank_==0 || thermal_rank_==1)
    valeur().set_thermal_subresolution_outputs(fic);
}

