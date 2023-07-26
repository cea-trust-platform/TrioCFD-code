/****************************************************************************
* Copyright (c) 2019, CEA
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
// File      : Triple_Line_Model_FT_Disc.h
// Directory : $TCL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Triple_Line_Model_FT_Disc_included
#define Triple_Line_Model_FT_Disc_included

#include <Objet_U.h>
#include <Domaine_VDF.h>
#include <FTd_tools.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <TRUST_Ref.h>

class Navier_Stokes_FT_Disc;
class Convection_Diffusion_Temperature_FT_Disc;

//#include <Maillage_FT_Disc.h>
class Param;
class Parcours_interface;
class Equation_base;
#define TCL_MODEL 1

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Triple_Line_Model_FT_Disc
//
// <Description of class Triple_Line_Model_FT_Disc>
//
/////////////////////////////////////////////////////////////////////////////

class Triple_Line_Model_FT_Disc : public Objet_U
{

  Declare_instanciable_sans_constructeur( Triple_Line_Model_FT_Disc ) ;

public :
  Triple_Line_Model_FT_Disc();

  void initialize();
  void completer();
  int associer_(Objet_U& ob) override;
  void set_param(Param& p);
  double get_Qtcl()
  {
//   return Qtcl_;
    //  Qtcl_ = coeffb_ + coeffa_*(theta_app_*22/(7*180));
//   Cerr << "coeffb " << coeffb_ << " coeffa " << coeffa_ << "theta " << (theta_app_) << finl;
    //  return (coeffb_ + coeffa_*(theta_app_))*(6.2/(theta_app_));
    return Qtcl_;
  };
  const double& get_lv() const
  {
    return lv_;
  };
  const double& get_ym() const
  {
    return ym_;
  };
  double get_sm()
  {
    return sm_;
  };
  const double& get_initial_CL_xcoord() const
  {
    return initial_CL_xcoord_;
  };
  const double& get_rhocpl() const
  {
    return rhocpl_;
  };
  const bool& is_activated() const
  {
    return activated_;
  };
  const int& is_capillary_activated() const
  {
    return capillary_effect_on_theta_activated_;
  };

  const int& tag_tcl() const
  {
    return tag_tcl_;
  };
  double compute_capillary_number() const;
  void get_in_out_coords(const Domaine_VDF& zvdf, const int elem,
                         const double dt,
                         DoubleTab& in_out, FTd_vecteur3& norm_elem, double& surface_tot);

  double compute_local_cos_theta(const Parcours_interface& parcours, const int num_face, const FTd_vecteur3& norm_elem) const;

  double compute_Qint(const DoubleTab& in_out, const double theta_app_locs,
                      const int num_face, double& Qmeso) const;

  void compute_TCL_fluxes_in_all_boundary_cells(ArrOfInt& elems_with_CL_contrib,
                                                ArrOfDouble& mpoint_from_CL, ArrOfDouble& Q_from_CL);

  void corriger_mpoint(DoubleTab& mpoint) const;
  void corriger_secmem(const double coef, DoubleTab& secmem2) const;
  void correct_wall_adjacent_temperature_according_to_TCL_fluxes(DoubleTab& temperature) const;
  void set_wall_adjacent_temperature_according_to_TCL_model(DoubleTab& temperature) const;
  void correct_TCL_energy_evolution(DoubleTab& temperature) const;
  void set_theta_app(double val)
  {
    theta_app_= val;
    Cerr << "theta_app inside set_theta_app " << theta_app_ << finl;
  };

  void associer_eq_temperature(const Convection_Diffusion_Temperature_FT_Disc& eq_temp);
  void associer_eq_ns(const Equation_base& eq);
  void associer_eq_interf(const Equation_base& eq);


  inline ArrOfInt& elems();
  inline const ArrOfInt& elems() const;
  inline ArrOfInt& boundary_faces();
  inline const ArrOfInt& boundary_faces() const;
  inline ArrOfDouble& mp();
  inline const ArrOfDouble& mp() const;
  inline ArrOfDouble& Q();
  inline const ArrOfDouble& Q() const;
//  double initial_CL_xcoord_;
protected :

  bool activated_;
  int capillary_effect_on_theta_activated_;
  int TCL_energy_correction_;
  int tag_tcl_;
  double coeffa_; //
  double coeffb_; //
  double Qtcl_; //
  double lv_; // The length...
  double theta_app_;
  double x_cl_; // position of the CL.
  // The end of micro-region:
  double ym_;
  double sm_;
// double old_xcl_ ; // Former position of the contact line. Not nice, but needed to copute the CL velocity
  double initial_CL_xcoord_; // Former position of the contact line. Not nice, but needed to copute the CL velocity
  double kl_cond_; // We store the liquid conductivity for easy access.
  double rhocpl_;

  // Information on the TCL region :
  // Note that the same elem may appear twice in the list, once for the micro contribution, once for the meso.
  ArrOfInt elems_; // The list of elements containing either the TCL itself, or the meso domaine.
  ArrOfInt boundary_faces_; // corresponding list of faces on the wall boundary.
  ArrOfDouble mp_; // corresponding value to be assigned to each cell (as a mass flux m [unit?])
  ArrOfDouble Q_;  // corresponding value to be assigned to each cell (as a thermal flux Q [unit?])

  // Variables for post-processing :
  double integration_time_;
  double instant_mmicro_evap_;
  double instant_mmeso_evap_;
  double instant_vmicro_evap_;
  double instant_vmeso_evap_;
  double integrated_m_evap_;
  double integrated_vmicro_evap_;
  double integrated_vmeso_evap_;
  double vevap_int_;

  REF(Convection_Diffusion_Temperature_FT_Disc)  ref_eq_temp_;
  REF(Transport_Interfaces_FT_Disc)  ref_eq_interf_;
  REF(Navier_Stokes_FT_Disc) ref_ns_;
};


inline ArrOfInt& Triple_Line_Model_FT_Disc::elems()
{
#ifdef DEBUG
  {
    const Transport_Interfaces_FT_Disc& eq_transport = ref_eq_interface_.valeur();
    const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
    assert(tcl.tag_tcl() == maillage.get_mesh_tag());
    toto code non lu...
  }
#endif
  const Transport_Interfaces_FT_Disc& eq_transport = ref_eq_interf_.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  if (tag_tcl() == maillage.get_mesh_tag())
    {
      Cerr << "Why do you need to access elems_ in a non-const maner when tag_cl_= " << tag_tcl()
           << " and maillage.get_mesh_tag()= " << maillage.get_mesh_tag()
           << " are equal (which means by-the-way that the model should be up-to-date already!)" << finl;
      Process::exit();
    }
  return elems_;
}

inline const ArrOfInt& Triple_Line_Model_FT_Disc::elems() const
{
  assert(tag_tcl() == ref_eq_interf_.valeur().maillage_interface().get_mesh_tag());
  return elems_;
}

inline ArrOfInt& Triple_Line_Model_FT_Disc::boundary_faces()
{
  return boundary_faces_;
}

inline const ArrOfInt& Triple_Line_Model_FT_Disc::boundary_faces() const
{
  assert(tag_tcl() == ref_eq_interf_.valeur().maillage_interface().get_mesh_tag());
  return boundary_faces_;
}

inline ArrOfDouble& Triple_Line_Model_FT_Disc::mp()
{
  return mp_;
}

inline const ArrOfDouble& Triple_Line_Model_FT_Disc::mp() const
{
  return mp_;
}

inline ArrOfDouble& Triple_Line_Model_FT_Disc::Q()
{
  return Q_;
}

inline const ArrOfDouble& Triple_Line_Model_FT_Disc::Q() const
{
  return Q_;
}

#endif /* Triple_Line_Model_FT_Disc_included */
