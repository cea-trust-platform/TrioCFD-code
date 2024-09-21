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

#include <Transport_Interfaces_FT_Disc.h>
#include <Domaine_VDF.h>
#include <FTd_tools.h>
#include <TRUST_Ref.h>

class Convection_Diffusion_Temperature_FT_Disc;
class Navier_Stokes_FT_Disc;
class Parcours_interface;
class Equation_base;
class Param;

#define TCL_MODEL 1

class Triple_Line_Model_FT_Disc : public Objet_U
{
  Declare_instanciable_sans_constructeur( Triple_Line_Model_FT_Disc ) ;
public :
  Triple_Line_Model_FT_Disc();

  void associer_pb(const Probleme_base&);
  void initialize();
  void completer();
  void set_param(Param& p);
  double get_Qtcl(const int num_face);
  const double& get_lv() const { return lv_; }
  const double& get_ym() const { return ym_; }
  double get_sm() { return sm_; }
  const double& get_initial_CL_xcoord() const { return initial_CL_xcoord_; }
  const double& get_rhocpl() const { return rhocpl_; }
  const bool& is_activated() const { return activated_; }
  const int& is_capillary_activated() const { return capillary_effect_on_theta_activated_; }
  const int& is_read_via_file() const { return read_via_file_; }
  const int& Rc_tcl_GridN() const { return Rc_tcl_GridN_; }
  const double& Rc_inject() const { return Rc_inject_; }
  const double& thetaC_tcl() const { return thetaC_tcl_; }
  const int& reinjection_tcl() const { return reinjection_tcl_; }
  const int& distri_first_facette() const { return distri_first_facette_; }
  const bool& ready_inject_tcl() const { return ready_inject_tcl_; }
  bool& ready_inject_tcl() { return ready_inject_tcl_; }
  const double& tempC_tcl() const { return tempC_tcl_; }
  const int& tag_tcl() const { return tag_tcl_; }

  double compute_capillary_number() const;
  int get_any_tcl_face() const;

  // Using the exact in/out intersections of the interface within the cell (works only in 2D, and may assume the wall in y-)
  void get_in_out_coords(const Domaine_VDF& zvdf, const int elem, const double dt, DoubleTab& in_out, FTd_vecteur3& norm_elem, double& surface_tot);

  // computes an approximate interface position (equivalent mean-plane)
  // returns in_out table and the mean normal in elem, and the interface surface
  void compute_approximate_interface_inout(const Domaine_VDF& zvdf, const int korient, const int elem, const int num_face, DoubleTab& in_out, FTd_vecteur3& norm_elem, double& surface_tot) const;

  double compute_local_cos_theta(const Parcours_interface& parcours, const int num_face, const FTd_vecteur3& norm_elem) const;

  // Computes the integral flux in a cell in the meso zone.
  // in and out are the interface position on both sides of the cell.
  double compute_Qint(const DoubleTab& in_out, const double theta_app_locs, const int num_face, double& Qmeso) const;

  void compute_TCL_fluxes_in_all_boundary_cells(ArrOfInt& elems_with_CL_contrib, ArrOfInt& faces_with_CL_contrib, ArrOfDouble& mpoint_from_CL, ArrOfDouble& Q_from_CL);

  void corriger_mpoint(DoubleTab& mpoint) const;
  void corriger_secmem(const double coef, DoubleTab& secmem2) const;
  void correct_wall_adjacent_temperature_according_to_TCL_fluxes(DoubleTab& temperature) const;
  void set_wall_adjacent_temperature_according_to_TCL_model(DoubleTab& temperature) const;
  void correct_TCL_energy_evolution(DoubleTab& temperature) const;
  double get_theta_app(const int num_face);

  enum InoutMethod { EXACT, APPROX, BOTH };

  //InoutMethod inout_method() const { return inout_method_; }
  inline ArrOfInt& elems();

  inline const ArrOfInt& elems() const
  {
    assert(tag_tcl() == ref_eq_interf_->maillage_interface().get_mesh_tag());
    return elems_;
  }

  inline ArrOfInt& boundary_faces() { return boundary_faces_; }

  inline const ArrOfInt& boundary_faces() const
  {
    assert(tag_tcl() == ref_eq_interf_->maillage_interface().get_mesh_tag());
    return boundary_faces_;
  }

  inline ArrOfDouble& mp() { return mp_; }
  inline const ArrOfDouble& mp() const { return mp_; }
  inline ArrOfDouble& Q() { return Q_; }
  inline const ArrOfDouble& Q() const { return Q_; }
//  double initial_CL_xcoord_;
protected:

  bool activated_;
  int deactivate_;
  int capillary_effect_on_theta_activated_;
  int TCL_energy_correction_;
  int n_ext_meso_; // number of layers in the extension of the meso-zone.
  int tag_tcl_;
  //double coeffa_; //
  //double coeffb_; //
  double Qtcl_; //
  double lv_; // The length...
  double theta_app_;
  double x_cl_; // position of the CL.
  // The end of micro-region:
  double ym_;
  double sm_;
  // End of meso region
  double ymeso_;
// double old_xcl_ ; // Former position of the contact line. Not nice, but needed to copute the CL velocity
  double initial_CL_xcoord_; // Former position of the contact line. Not nice, but needed to copute the CL velocity
  double kl_cond_; // We store the liquid conductivity for easy access.
  double rhocpl_;
  int inout_method_ = 0;
  // Tabulation TCL model
  int read_via_file_;
  Nom Nom_ficher_tcl_;
  DoubleTab tab_Mtcl_;
  // total nb of colons
  int nb_colon_tab_;
  // num -colon for temperature
  int num_colon_tem_;
  // num -colon for velocity
  // int num_colon_vel_;
  // num -colon for theta app
  int num_colon_app_;
  // num -colon for Qtcl
  int num_colon_qtcl_;

  // parameters to force numerical breakup
  // Radius of nucleate site; [in number of grids]
  // Distance below which the theta_app will be replaced by a big contact angle thetaC_tcl_ [in degree]
  // to force bubble pinching / necking
  int Rc_tcl_GridN_;
  double Rc_inject_;
  double thetaC_tcl_;
  // parameters injection seed nucleate
  // Size: Using the same Radius defined above (Rc_tcl_gridN_), [in number of grids]
  // Use the average temperature in the same zone to tell if nucleate site is activted and inject bubbles
  // Temperature beyond which the site is activeted [in K] T-Tsat
  // in this case,  we need the initial value for contact angle, thetaC_tcl_ will be used
  // reinjection: if a re-injection of interface if nessaire
  int reinjection_tcl_;
  double tempC_tcl_;
  bool ready_inject_tcl_;
  // summaries:
  //        - when reinjection_tcl_ is 0;  i.e., by default, the pinching breakup with be used,
  //           Rc_tcl_GridN_ thetaC_tcl_ control the breakup, thetaC_tcl_ should be large
  //        - when reinjection_tcl_ is 1;  i.e., there is no the pinching breakup. The bubble goes away directly, we reinject seed
  //           Rc_tcl_GridN_ thetaC_tcl_ control the initial shape, thetaC_tcl_ should be small...

  // if to distribute the Qtcl into the first facette
  int distri_first_facette_;

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

  Nom nom_eq_hydr_, nom_eq_therm_, nom_eq_interf_;

  REF(Probleme_base) pb_base_;
  REF(Convection_Diffusion_Temperature_FT_Disc) ref_eq_temp_;
  REF(Transport_Interfaces_FT_Disc) ref_eq_interf_;
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

#endif /* Triple_Line_Model_FT_Disc_included */
