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
* OR BUSINESS INTelERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Triple_Line_Model_FT_Disc.cpp
// Directory : $TCL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Triple_Line_Model_FT_Disc.h>
#include <Param.h>
#include <Convection_Diffusion_Temperature_FT_Disc.h>
// passed into the .h for assert in inline method... #include <Transport_Interfaces_FT_Disc.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Domaine_VF.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <LecFicDiffuse.h>
#include <Domaine_Cl_VDF.h>
#include <Echange_contact_VDF_FT_Disc.h>


static const double TSAT_CONSTANTE = 0.;

Implemente_instanciable_sans_constructeur( Triple_Line_Model_FT_Disc, "Triple_Line_Model_FT_Disc", Objet_U ) ;
// XD Triple_Line_Model_FT_Disc objet_u Triple_Line_Model_FT_Disc -1 Triple Line Model (TCL)
// XD ref navier_stokes_ft_disc navier_stokes_ft_disc
// XD ref transport_interfaces_ft_disc transport_interfaces_ft_disc
// XD ref convection_diffusion_temperature_ft_disc convection_diffusion_temperature_ft_disc


Triple_Line_Model_FT_Disc::Triple_Line_Model_FT_Disc() :
  // Values initialized:
  activated_(false),
  deactivate_(0),
  capillary_effect_on_theta_activated_(0),
  TCL_energy_correction_(0),
  n_ext_meso_(2),
  tag_tcl_(-1),
  //coeffa_(0),
  //coeffb_(0),
  Qtcl_(0),
  lv_(0.),
  theta_app_(0.),
  x_cl_(0.),
  ym_(0.),
  sm_(0.),
  ymeso_(0.),
//  old_xcl_(0.),
  initial_CL_xcoord_(0.),
  kl_cond_(-1.),  // Invalid
  rhocpl_(-1.),  // Invalid
  read_via_file_(0),
  nb_colon_tab_(5),
  num_colon_tem_(0),
  num_colon_app_(1),
  num_colon_qtcl_(2),
  Rc_tcl_GridN_(4),
  Rc_inject_(8e-05),
  thetaC_tcl_(150.),
  reinjection_tcl_(0),
  tempC_tcl_(10.),
  ready_inject_tcl_(0),
  distri_first_facette_(0),
  integration_time_(0.),
  instant_mmicro_evap_(0.),
  instant_mmeso_evap_(0.),
  instant_vmicro_evap_(0.),
  instant_vmeso_evap_(0.),
  integrated_m_evap_(0.),
  integrated_vmicro_evap_(0.),
  integrated_vmeso_evap_(0.),
  vevap_int_(0.)
{
}

Sortie& Triple_Line_Model_FT_Disc::printOn( Sortie& os ) const
{
  //Objet_U::printOn(os);
  /*  os << "{\n"
  		     << "   lv " << lv_ << "\n"
  		     << "   Qtcl " << Qtcl_ << "\n"
       << "}\n" ;
  */
  Param p(que_suis_je());
  // Je peux pas car c'est const... set_param(p);
  p.print(os);
  return os;
}

Entree& Triple_Line_Model_FT_Disc::readOn( Entree& is )
{
  Param p(que_suis_je());
  set_param(p);
  p.lire_avec_accolades_depuis(is);
  if (!deactivate_)
    activated_ = true;
  if (read_via_file_)
    {
      LecFicDiffuse fic;
      if(!fic.ouvrir(Nom_ficher_tcl_))
        {
          Cerr << " File " <<Nom_ficher_tcl_<< " doesn't exist. To generate it, please, check JDD " << finl;
          exit();
        }




      Cerr << "Reading of the ascii file : " << Nom_ficher_tcl_ << finl;

      const int nb_colons = nb_colon_tab_; // nb colons in the readed file

      for (int k = 0; k < nb_colons; k++)
        {
          Nom name_colon_;
          fic >> name_colon_;
          Cerr << "--Name of the " << k + 1 << "-th colon : " << name_colon_
               << finl;
        }


      double tmp;
      ArrOfDouble list_val;

      while (fic.get(&tmp, 1))
        {
          list_val.append_array(tmp);
        }

      int nb_lines = list_val.size_array () / nb_colons;
      tab_Mtcl_.resize (nb_colons, nb_lines);

      for (int i = 0; i < nb_lines; i++)
        for (int k = 0; k < nb_colons; k++)
          {
            tab_Mtcl_ (k, i) = list_val[i * nb_colons + k];
          }
    }
  return is;
}

// Description: Methode appelee par readOn. Declaration des parametres a lire dans le .data
void Triple_Line_Model_FT_Disc::set_param(Param& p)
{
  p.ajouter("Qtcl", &Qtcl_); // XD_ADD_P floattant Heat flux contribution to micro-region [W/m]
  p.ajouter("lv", &lv_); // XD_ADD_P floattant Slip length (unused)
  //p.ajouter("coeffa", &coeffa_); // XD_ADD_P floattant not_set
  //p.ajouter("coeffb", &coeffb_); // XD_ADD_P floattant not_set
  p.ajouter("theta_app", &theta_app_); // XD_ADD_P floattant Apparent contact angle (Cox-Voinov)
  //p.ajouter("ylim", &ym_); // XD_ADD_P floattant not_set
  p.ajouter("ym", &ym_); // XD_ADD_P floattant Wall distance of the point M delimiting micro/meso transition [m]
  p.ajouter("sm", &sm_,Param::REQUIRED); // XD_ADD_P floattant Curvilinear abscissa of the point M delimiting micro/meso transition [m]
  p.ajouter("ymeso", &ymeso_); // XD_ADD_P floattant Meso region extension in wall-normal direction [m]
  p.ajouter("n_extend_meso", &n_ext_meso_); // XD_ADD_P entier Meso region extension in number of cells [-]
  p.ajouter("initial_CL_xcoord", &initial_CL_xcoord_); // XD_ADD_P floattant Initial interface position (unused)
  p.ajouter("read_via_file", &read_via_file_);
  p.ajouter("Rc_tcl_GridN", &Rc_tcl_GridN_); // XD_ADD_P floattant Radius of nucleate site; [in number of grids]
  // p.ajouter("Rc_inject", &Rc_inject_);
  p.ajouter("thetaC_tcl", &thetaC_tcl_); // XD_ADD_P floattant imposed contact angle [in degree] to force bubble pinching / necking once TCL entre nucleate site
  p.ajouter_flag("reinjection_tcl", &reinjection_tcl_); // XD_ADD_P flag This flag activates the automatic injection of a new nucleate seed with a specified shape when the temperature in the nucleation site becomes higher than a certain threshold (tempC_tcl). The shape of the seed is determined by the radius Rc_tcl_GridN and the contact angle thetaC_tcl. The nucleation site is considered free when there are no bubbles present. The site size is defined by Rc_tcl_GridN. This temperature threshold, termed tempC_tcl, is the activation temperature. Setting this temperature implies a wall temperature, therefore, activating reinjection_tcl is ONLY possible for a simulation coupled with solid conduction. NL2 When reinjection_tcl is activated, the values of tempC_tcl (default 10K), Rc_tcl_GridN (default 4 grid sizes), and thetaC_tcl (default 150 degrees) should be provided. Unless (STRONGLY not recommended), the default values (indicated in parentheses) will be used. NL2 If reinjection_tcl is not activated (by default), the mechanism of Numerically forcing bubble pinching/necking will be used for multi-cycle simulation. Once the Triple Contact Line (TCL) enters the nucleation site, a big contact angle thetaC_tcl is imposed to initiate bubble pinching/necking. After the bubble pinching ends, the large bubble above will depart, leaving the remaining part to serve as the nucleate seed. This process is equivalent to immediately inserting a new seed with a prescribed shape (determined by the nucleation site size and contact angle) once a bubble departs. Site size is defined by Rc_tcl_GridN (default 4 grid sizes). Contact angle thetaC_tcl (default 150 degrees). Useful for a standalone (not coupling with solid conduction) simulation.
  p.ajouter_flag("distri_first_facette", &distri_first_facette_); // XD_ADD_P flag This flag determines whether to distribute the Qtcl into all grids occupied by the first facette according to their area proportions. When set, the flux is redistributed into all grids occupied by the first facette based on their area proportions. Default value is 0, the flux is distributed differently: similar to the Meso zone, it is only distributed to grids within the Micro-zone (where the height of the front y is smaller than the size of Micro ym). The distribution of this flux is logarithmically proportional to y between 5.6nm (here interpreted as the value 0 in logarithm) and ym. In practice, in most cases, it will distribute all the flux locally in the first grid.
  p.ajouter("tempC_tcl", &tempC_tcl_); // see in _.h file
  p.ajouter("file_name", &Nom_ficher_tcl_); // XD_ADD_P floattant Input file to set TCL model
  p.ajouter("nb_colon", &nb_colon_tab_);
  p.ajouter("num_colon_tem", &num_colon_tem_);
  // p.ajouter("num_colon_vel", &num_colon_vel_);
  p.ajouter("num_colon_app", &num_colon_app_);
  p.ajouter("num_colon_qtcl", &num_colon_qtcl_);
  p.ajouter_flag("enable_energy_correction", &TCL_energy_correction_);
  p.ajouter_flag("capillary_effect_on_theta", &capillary_effect_on_theta_activated_);
  p.ajouter_flag("deactivate", &deactivate_); // XD_ADD_P flag Simple way to disable completely the TCL model contribution
  p.ajouter("inout_method", &inout_method_); // XD_ADD_P chaine(into=["exact","approx","both"]) Type of method for in out calc. By defautl, exact method is used
  p.dictionnaire("exact", EXACT);
  p.dictionnaire("approx", APPROX);
  p.dictionnaire("both", BOTH); // Makes both and compare them
}

void Triple_Line_Model_FT_Disc::associer_eq_temperature(const Convection_Diffusion_Temperature_FT_Disc& eq_temp)
{
  if (! sub_type(Convection_Diffusion_Temperature_FT_Disc, eq_temp))
    {
      Cerr << "Error for associer_eq_temperature is not given what we expect\n"
           << "A Convection_Diffusion_Temperature_FT_Disc medium was expected." << finl;
      Process::exit();
    }
  ref_eq_temp_ = ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq_temp);
}

void Triple_Line_Model_FT_Disc::associer_eq_interf(const Equation_base& eq)
{
  if (!sub_type(Transport_Interfaces_FT_Disc,eq))
    {
      Cerr<<"No interface transport equation has been associated to TCL model"<<finl;
      Process::exit();
    }
  ref_eq_interf_ = ref_cast(Transport_Interfaces_FT_Disc, eq);
}

void Triple_Line_Model_FT_Disc::associer_eq_ns(const Equation_base& eq)
{
  if (!sub_type(Navier_Stokes_FT_Disc,eq))
    {
      Cerr<<"No NS equation has been associated to TCL model"<<finl;
      Process::exit();
    }
  ref_ns_ = ref_cast(Navier_Stokes_FT_Disc, eq);
}

int Triple_Line_Model_FT_Disc::associer_(Objet_U& ob)
{
  if (sub_type(Transport_Interfaces_FT_Disc, ob))
    {
      associer_eq_interf(ref_cast(Transport_Interfaces_FT_Disc, ob));
      return 1;
    }
  else if (sub_type(Navier_Stokes_FT_Disc,ob))
    {
      associer_eq_ns(ref_cast(Navier_Stokes_FT_Disc,ob));
      return 1;
    }
  else if (sub_type(Convection_Diffusion_Temperature_FT_Disc,ob))
    {
      associer_eq_temperature(ref_cast(Convection_Diffusion_Temperature_FT_Disc,ob));
      return 1;
    }
  else
    {
      Cerr << "Error in TCL association to obj " << ob.que_suis_je() << finl;
      Process::exit();
    }
  return 0;
}

void Triple_Line_Model_FT_Disc::initialize()
{
  if (!ref_eq_interf_.non_nul())
    {
      Cerr << "No interfacial transport has been associated to TCL model" << finl;
      Cerr << "Please use Associate [TCL] [Transport_Interface_FT_Disc]" << finl;
      Process::exit();
    }
  if (!ref_eq_temp_.non_nul())
    {
      Cerr << "No temperature equation has been associated to TCL model" << finl;
      Cerr << "Please use Associate [TCL] [Convection_Diffusion_FT_Disc]" << finl;
      Process::exit();
    }
  if (!ref_ns_.non_nul())
    {
      Cerr << "No NS equation has been associated to TCL model" << finl;
      Cerr << "Please use Associate [TCL] [Navier_Stokes_FT_Disc]" << finl;
      Process::exit();
    }
  const Transport_Interfaces_FT_Disc& transport = ref_eq_interf_.valeur();
  const Maillage_FT_Disc& maillage = transport.maillage_interface();
  tag_tcl_ = maillage.get_mesh_tag();

  // Issue: we don't know which phase is for the temperature at the initialize step.
  // It is required to store kl_cond_
  // So we do it in the completer.

  // how to access fluid diphasique? Through (eq_ns or pb)? We have ns.
  // const Milieu_base& milieu = ref_eq_temp_.valeur().milieu();
  // const Fluide_Diphasique& fluid_dipha = ref_cast(Fluide_Diphasique, milieu); -> no it's not a Fluide_diphasique
  // Or maybe we can give access to fluide_dipha_ from Convection_Diffusion_Temperature_FT_Disc.
  //
  // int phase = ref_eq_temp_.valeur().get_phase();
  // L_vap_ = fluide_dipha.valeur().chaleur_latente();

  // We may compute the true position here for old_xcl_
  // All of (integration_time_, instant_m_evap_, instant_vmicro_evap_, instant_vmeso_evap_,
  // integrated_m_evap_, integrated_vmicro_evap_, integrated_vmeso_evap_)
  // are set to zero by construction. For restart, something more tricky should be done.

  // Information on the TCL region :
  // Resize the tables :
  elems_.resize_array(0);
  boundary_faces_.resize_array(0);
  mp_.resize_array(0);
  Q_.resize_array(0);
}

int Triple_Line_Model_FT_Disc::get_any_tcl_face() const
{
  int num_bc_face=-1;
  const Transport_Interfaces_FT_Disc& transport = ref_eq_interf_.valeur();
  const Maillage_FT_Disc& maillage = transport.maillage_interface();
  for (int s = 0; s < maillage.nb_sommets(); s++)
    {
      if (maillage.sommet_ligne_contact(s) and (not(maillage.sommet_virtuel(s))))
        {
          num_bc_face = maillage.sommet_face_bord(s);
          break;
        }
    }
  return num_bc_face;
}

void Triple_Line_Model_FT_Disc::completer()
{
  // Via the temperature transport equation, we directly get access to a Fluide_Incompressible:
  const Milieu_base& milieu = ref_eq_temp_.valeur().milieu();
  kl_cond_ = milieu.conductivite()(0,0);
  rhocpl_ =  milieu.masse_volumique()(0,0) * milieu.capacite_calorifique()(0,0);

  if ((n_ext_meso_ != 1)and(ymeso_>DMINFLOAT))
    {
      Cerr << "Check your datafile. n_extend_meso and ymeso should not but used simultaneously." << finl;
      Cerr << "n_extend_meso = " << n_ext_meso_<< finl;
      Cerr << "ymeso = " << ymeso_ <<finl;
      Process::exit();
    }
  if (n_ext_meso_ <2 )
    {
      Cerr << "n_extend_meso = " << n_ext_meso_<< finl;
      Cerr << "Check your datafile. n_extend_meso should be >= 2" << finl;
      Process::exit();
    }

  if (ymeso_==0.)
    {
      // ymeso has not been initialised. n_ext_meso is used instead:
      const Navier_Stokes_FT_Disc& ns = ref_ns_.valeur();
      const Domaine_Cl_dis_base& zcldis = ns.domaine_Cl_dis().valeur();
      // get wall BC frontiere nb
      int num_bord = -1;
      for (int i=0; i<zcldis.nb_cond_lim(); i++)
        {
          const Cond_lim& la_cl = zcldis.les_conditions_limites(i);
          // For each BC, we check its type to see if it's a wall:
          // BC for hydraulic equation
          bool is_wall = sub_type(Dirichlet_paroi_fixe,la_cl.valeur())
                         || sub_type(Dirichlet_paroi_defilante,la_cl.valeur());
          if (is_wall)
            {
              num_bord = i;
              break;
            }
        }
      if (num_bord<0)
        Process::exit( "[TCL: ymeso]: !!! NO WALL-TYPE BOUNDARY WAS FOUND IN THE DOMAINE, PLESEASE CHECK JDD" );

      const Frontiere& fr=zcldis.les_conditions_limites(num_bord).frontiere_dis().frontiere();
      const int nb_face = fr.nb_faces();
      const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, ns.domaine_dis().valeur());
      // get TCL cells
      if (nb_face>0)
        {
          const int num_face= fr.num_premiere_face();;
          int elem_voisin = zvdf.face_voisins(num_face, 0) + zvdf.face_voisins(num_face, 1) +1;
          const double cell_height = 2.*std::fabs(zvdf.dist_face_elem0(num_face,elem_voisin));
          ymeso_ = cell_height * n_ext_meso_;
        }
      ymeso_ = Process::mp_max(ymeso_);
      Cerr << "[TCL] ymeso is set to " << ymeso_ << " according to provided n_ext_meso_ " << n_ext_meso_ << finl;
    }

  // ymeso_ -= Objet_U::precision_geom;
  ymeso_ = std::fmax(ymeso_,DMINFLOAT);

  if (reinjection_tcl())
    {
      const Navier_Stokes_FT_Disc& ns = ref_ns_.valeur();
      const Domaine_Cl_dis_base& zcldis = ns.domaine_Cl_dis().valeur();
      // get wall BC frontiere nb
      int num_bord = -1;
      for (int i=0; i<zcldis.nb_cond_lim(); i++)
        {
          const Cond_lim& la_cl = zcldis.les_conditions_limites(i);
          // For each BC, we check its type to see if it's a wall:
          // BC for hydraulic equation
          bool is_wall = sub_type(Dirichlet_paroi_fixe,la_cl.valeur())
                         || sub_type(Dirichlet_paroi_defilante,la_cl.valeur());
          if (is_wall)
            {
              num_bord = i;
              break;
            }
        }
      if (num_bord<0)
        Process::exit( "[TCL: Re-injection]: !!! NO WALL-TYPE BOUNDARY WAS FOUND IN THE DOMAINE, PLESEASE CHECK JDD" );

      const Frontiere& fr=zcldis.les_conditions_limites(num_bord).frontiere_dis().frontiere();
      const int nb_face = fr.nb_faces();
      const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, ns.domaine_dis().valeur());
      // supposing deltax = deltay in the first thermal layer
      if (nb_face>0)
        {
          const int num_face=  fr.num_premiere_face();;
          int elem_voisin = zvdf.face_voisins(num_face, 0) + zvdf.face_voisins(num_face, 1) +1;
          const double cell_height = 2.*std::fabs(zvdf.dist_face_elem0(num_face,elem_voisin));
          Rc_inject_ = cell_height * Rc_tcl_GridN_;
        }
      Rc_inject_ = Process::mp_max(Rc_inject_);
      Cerr << "[TCL] Rc_inject_ is set to " << Rc_inject_ << " according to provided Rc_tcl_GridN_ " << Rc_tcl_GridN_ << finl;
    }
  // Rc_tcl_GridN_ -= Objet_U::precision_geom;
  Rc_inject_ = std::fmax(Rc_inject_,DMINFLOAT);

  // how to access fluid diphasique? Through (eq_ns or pb)? We have neither so far.
  // const Fluide_Diphasique& fluid_dipha = ref_cast(Fluide_Diphasique, milieu); -> no it's not a Fluide_diphasique
  // Or maybe we can give access to fluide_dipha_ from Convection_Diffusion_Temperature_FT_Disc.
  //
  // int phase = ref_eq_temp_.valeur().get_phase();
  // L_vap_ = fluide_dipha.valeur().chaleur_latente();

  // We may compute the true position here for old_xcl_
  // All of (integration_time_, instant_m_evap_, instant_vmicro_evap_, instant_vmeso_evap_,
  // integrated_m_evap_, integrated_vmicro_evap_, integrated_vmeso_evap_)
  // are set to zero by construction. For restart, something more tricky should be done.
}

double Triple_Line_Model_FT_Disc::compute_capillary_number() const
{


  return 0.;
}

// Description : For a given elem in the domaine, computes the crossing points of the interface and the averaged normal.
//               The method is not const because it changes the value of old_xcl_.
// 				 I believe there is a convention that y_in <= y_out
// Parameters :
// elem -> the id of the cell of interest in the domaine
// norm_elem : the averaged normal into the element.
// surface_tot : The local interfacial area [m2].
// in_out : Real coordinates of the interface intersection with the cell
//
// in_out(0,0) -> x_in
// in_out(0,1) -> y_in
// in_out(1,0) -> x_out
// in_out(1,1) -> y_out
void Triple_Line_Model_FT_Disc::get_in_out_coords(const Domaine_VDF& zvdf, const int elem,
                                                  const double dt,
                                                  DoubleTab& in_out, FTd_vecteur3& norm_elem, double& surface_tot)
{
  const int dim = Objet_U::dimension;
  const Transport_Interfaces_FT_Disc& transport = ref_eq_interf_.valeur();
  const Maillage_FT_Disc& maillage = transport.maillage_interface();
  const IntTab& facettes = maillage.facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const DoubleTab& sommets = maillage.sommets();

  in_out.resize(2, dim); // The first param '2' is for in and out
  if (dim != 2)
    {
      Cerr << que_suis_je() << "::get_in_out_coords() only coded in 2 dimensions" << finl;
      Process::exit();
    }

  // We get the barycenter of the boundary face:
  //     const double xG = zvdf.xv(num_face, 0);
  //     const double yG = zvdf.xv(num_face, 1);
  //     double zG = 0.;
  //     if (dim3)
  //       zG = zvdf.xv(num_face, 2);

  //...VP-modifications for calculating the intersections.....
  // We get the barycenter of the center of the elem:
  //        const double fraction = data.fraction_surface_intersection_;
  const double xP = zvdf.xp(elem, 0);
  const double yP = zvdf.xp(elem, 1);
  //Cerr << "Barycentre coordinates are " << " xP = " << xP << "and yP =" << yP << finl;

  // We get the width of the elem:
  const double delta_x = zvdf.dim_elem(elem, 0);
  const double delta_y = zvdf.dim_elem(elem, 1);

  //x coordinates of the barycenter of the right and left faces:
  double rface_dis = xP + delta_x/2;
  double lface_dis = xP - delta_x/2;
  //y coordinates of the barycenter of the top and bottom faces:
  double tface_dis = yP + delta_y/2;
  double bface_dis = yP - delta_y/2;
  DoubleTab coo_faces(2,dim);
  coo_faces(0,0) = lface_dis;
  coo_faces(1,0) = rface_dis;
  coo_faces(0,1) = bface_dis;
  coo_faces(1,1) = tface_dis;

  FTd_vecteur3 bary_facettes_dans_elem = {0., 0., 0.} ;
  // Boucle sur les facettes qui traversent l'element elem :
  double yl = 0.;
  double yr = 0.;
  double xl = 0.;
  double xr = 0.;
  // We anticipate for more than 2 intersections because with tolerances, when the FT sommet is
  // close to the element face, 2 intersections will be found instead of one.
  double xint[4] = {}; // X coords of the intersections.
  double yint[4] = {}; // Y coords of the intersections.
  int nint = 0;
  const double eps = Objet_U::precision_geom;
  int index=intersections.index_elem()[elem];
  while (index >= 0)
    {
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
      const int fa7 = data.numero_facette_;
      const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
      const double surface_facette = surface_facettes[fa7];
      const double surf = data.fraction_surface_intersection_ * surface_facette;
      surface_tot +=surf;
      FTd_vecteur3 coord_barycentre_fraction = {0., 0., 0.};
      FTd_vecteur3 coord_sommet = {0., 0., 0.};
      for (int dir = 0; dir< Objet_U::dimension; dir++)
        {
          const double nx = normale_facettes(fa7,dir);
          norm_elem[dir] += nx * surf;
        }
      int inner_nodes = 0;

      for (int isom = 0; isom< Objet_U::dimension; isom++)
        {
          const int num_som = facettes(fa7, isom); // numero du sommet dans le tableau sommets
          // Cerr << " node-number #" << num_som << "face-number #" << fa7 << finl;

          const double bary_som = data.barycentre_[isom];
          for (int dir = 0; dir < Objet_U::dimension; dir++)
            {
              //modifications
              coord_sommet[dir] = sommets(num_som,dir);
              coord_barycentre_fraction[dir] += bary_som * sommets(num_som,dir) * surf;
              //  coord_barycentre_fraction[dir] += bary_som * sommets(num_som,dir);
            }

          //  Cerr << "coordinate x = " << coord_sommet[0] <<  " coordinate y = " << coord_sommet[1] << finl;
          //            if(coord_sommet[0] < rface_dis && coord_sommet[0] > lface_dis &&
          //                coord_sommet[1] < tface_dis && coord_sommet[1] >= bface_dis)
          if(coord_sommet[0] - rface_dis< eps && coord_sommet[0] - lface_dis > -eps &&
              coord_sommet[1]-tface_dis < eps && coord_sommet[1] - bface_dis > -eps)
            {
              //  Cerr << " This node lies inside the element #" << elem << finl;
              inner_nodes += 1;
            }
        }
      // Cerr << " nb of inner nodes = " << inner_nodes << " in elem " << elem
      //      << " when runnning fa7= " << fa7 << finl;
      // Cerr << " inside " << inside[0] << " "<< inside[1] << " "<< inside[2] << " "<< finl;

      if (inner_nodes >= 2)
        {
          // Cerr << " This facette lies totally inside the element #" << elem << finl;
          // Still we check if one intersection is close to the border :
          for (int isom = 0; isom< Objet_U::dimension; isom++)
            {
              const int node_number = facettes(fa7, isom); // numero du sommet dans le tableau sommets
              // Cerr << " node-number #" << node_number << "face-number #" << fa7 << finl;
              for (int dir = 0; dir < Objet_U::dimension; dir++)
                {
                  coord_sommet[dir] = sommets(node_number,dir);
                }
              if((std::abs(coord_sommet[0]-lface_dis) <= eps)
                  ||(std::abs(coord_sommet[0]-rface_dis) <= eps)
                  ||(std::abs(coord_sommet[1]-bface_dis) <= eps)
                  ||(std::abs(coord_sommet[1]-tface_dis) <= eps))
                {
                  xint[nint] = coord_sommet[0];
                  yint[nint] = coord_sommet[1];
                  nint += 1;
                  Cerr<< "Intersection found at node " << node_number << " belonging to " <<fa7 << finl;
                }
            }
          // Cerr<< "nint = "<< nint << finl;
        }
      else if((inner_nodes == 1)||(inner_nodes == 0))
        {
          // TODO : If this method works also when there are 2 intersections, we
          // could remove the for loop before and the if (....) ... else if ... and just keep below:
          //  Cerr << " find the intersection points of the facettes and elem-faces "<< finl;
          // One point in, one point out, search  1 intersect
          // OR Two points out, search 2 intersect
          const int s0 = facettes(fa7, 0);
          const int s1 = facettes(fa7, 1);
          const double x0 = sommets(s0,0);
          const double x1 = sommets(s1,0);
          const double y0 = sommets(s0,1);
          const double y1 = sommets(s1,1);
          for (int jj=0; jj< 2; jj++)
            {
              const double xF = coo_faces(jj,0);
              if ((x0-xF)*(x1-xF)<eps*(std::fabs(x0-x1)+eps))
                {
                  // points s0 and s1 lies on both sides of face jj.
                  double tmp_yint = 0.;
                  if (std::fabs(x0-x1)< eps)
                    {
                      tmp_yint = y0;
                      Cerr << "Almost // to the face, which intersect?" <<finl;
                      Cerr << "Should we count two extrems as intersect" << finl;
                      xint[nint] = xF;
                      yint[nint] = std::max(bface_dis,std::min(tface_dis, y0));
                      nint += 1;
                      xint[nint] = xF;
                      yint[nint] = std::max(bface_dis,std::min(tface_dis, y1));
                      nint += 1;
                      Process::exit();
                    }
                  else
                    tmp_yint = y0-(x0-xF)*(y0-y1)/(x0-x1);
                  // Is the intersection inside the cell:
                  if ((bface_dis-tmp_yint)*(tface_dis-tmp_yint)<eps*delta_y)
                    {
                      // lface or rface is cut:
                      xint[nint] = xF;
                      yint[nint] = tmp_yint;
                      nint += 1;
                    }
                }
              const double yF = coo_faces(jj,1);
              if ((y0-yF)*(y1-yF)<eps*(std::fabs(y0-y1)+eps))
                {
                  // points s0 and s1 lies on both sides of face jj.
                  double tmp_xint = 0.;
                  if (std::fabs(y0-y1)< eps)
                    {
                      tmp_xint = x0;
                      Cerr << "Almost // to the face, which intersect?" <<finl;
                      Cerr << "Should we count two extrems as intersect" << finl;
                      yint[nint] = yF;
                      xint[nint] = std::max(lface_dis,std::min(rface_dis, x0));
                      nint += 1;
                      yint[nint] = yF;
                      xint[nint] = std::max(lface_dis,std::min(rface_dis, x1));
                      nint += 1;
                      Cerr << "Code to be validated" << finl;
                      Process::exit();
                    }
                  else
                    tmp_xint = x0-(y0-yF)*(x0-x1)/(y0-y1);
                  // Is the intersection inside the cell:
                  if ((lface_dis-tmp_xint)*(rface_dis-tmp_xint)<eps*delta_x)
                    {
                      // bface or tface is cut:
                      yint[nint] = yF;
                      xint[nint] = tmp_xint;
                      nint += 1;
                    }
                }
            }
        }
      else
        {
          Cerr << "WTF inner_nodes : " << inner_nodes << finl;
          Process::exit();
        }

      for (int dir = 0; dir < Objet_U::dimension; dir++)
        bary_facettes_dans_elem[dir] += coord_barycentre_fraction[dir];
      index = data.index_facette_suivante_;
    }

  // removing two closing point
  if (nint > 2)
    {
      for (int i=0; i < nint; i++)
        {
          for (int j=0; j < nint; j++)
            {
              if (i < j)
                {
                  const double dist = sqrt(pow((xint[i]-xint[j]), 2) + pow((yint[i]-yint[j]), 2));
                  if (dist < Objet_U::precision_geom)
                    {
                      if (j!=nint-1)
                        {

                          xint[j] = xint[nint-1];
                          yint[j] = yint[nint-1];
                          nint--;
                          j--;
                        }
                      else
                        {
                          nint--;


                        }
                    }
                }

            }
        }
    }
  // We should have found 2 intersections in the end:
  if (nint != 2)
    {
      if (nint > 2)
        {
          Cerr << "Too many intersections " << nint << " -> trying to suppress duplicates" << finl;
          Cerr << "Too many points, maybe a FT sommet was close to a border and is counted twice" << finl;
          Cerr << "Intersections are: " << finl;
          for (int i=0; i<nint; i++)
            Cerr << xint[i] << " " << yint[i] << finl;
          Cerr << "In cell: " << elem << " at : " << finl;
          Cerr << "x in : " << lface_dis << " " << lface_dis << finl;
          Cerr << "y in : " << bface_dis << " " << tface_dis << finl;
          for (int i=0; i<nint-1; i++)
            {
              for (int i2=0; i2<nint; i2++)
                {
                  double dist = std::sqrt(pow(xint[i]-xint[i2],2) +pow(yint[i]-yint[i2],2));
                  if (dist < Objet_U::precision_geom)
                    {
                      // Ce sont les memes points :
                      if (i2!= nint)
                        {
                          // Ce n'est pas le dernier point. On met le dernier point a la place de i2 puis on refait i2:
                          xint[i2] = xint[nint-1];
                          yint[i2] = yint[nint-1];
                        }
                      // On supprime donc le dernier:
                      nint--;
                    }
                }
            }
        }
      if (nint > 2)
        {
          Cerr << "Wrong number of intersections " << nint << finl;
          Cerr << "STILL Too many points, maybe a FT sommet was close to a border and is counted twice" << finl;
          Cerr << "Intersections are: " << finl;
          for (int i=0; i<nint; i++)
            Cerr << xint[i] << " " << yint[i] << finl;
          Cerr << "In cell: " << elem << "at : " << finl;
          Cerr << "x in : " << lface_dis << " " << lface_dis << finl;
          Cerr << "y in : " << bface_dis << " " << tface_dis << finl;
          // TODO : le tri n'a pas fonctionne? trier, supprimer le/les doublons et decrementer le nint.
          Process::exit();
        }
    }

  if (surface_tot > 0.)
    {
      for (int dir = 0; dir< Objet_U::dimension; dir++)
        norm_elem[dir] *= 1./surface_tot;
    }
// Cerr<< "nint = "<< nint << finl;
//  Cerr << " xint[0] = " << xint[0] << " xint[1] = " << xint[1] << finl;
  if(xint[1] > xint[0])
    {
      xr = xint[1];
      xl = xint[0];
      yr = yint[1];
      yl = yint[0];
    }
  else
    {
      xr = xint[0];
      xl = xint[1];
      yr = yint[0];
      yl = yint[1];
    }
//  Cerr << "Elem #" << elem << " yr = " << yr << " xr = " << xr
//       << " yl = " << yl << " xl = " << xl <<finl;
// Fill in the table :
  in_out(0,0) = xl;
  in_out(0,1) = yl;
  in_out(1,0) = xr;
  in_out(1,1) = yr;
}

// const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ref_eq_temp_.valeur().domaine_dis().valeur());
// const IntTab& faces_elem = domaine_vf.face_voisins();
// // One of the neighbours doesnot exist so it has "-1". We get the other elem by:
// const int elem = faces_elem(num_face, 0) + faces_elem(num_face, 1) +1;
//  const double d = compute_distance(ref_cast(Domaine_VF, ref_eq_temp_.valeur().domaine_dis().valeur()),
//                   num_face, elem);
//
// Computes the distance between a face centre and an element centre.
double compute_distance(const Domaine_VF& domaine_vf, const int num_face, const int elem)
{
  double d=0;
  double P[3] = {0.,0.,0.}, xyz_face[3] = {0.,0.,0.};
  xyz_face[0] =  domaine_vf.xv(num_face,0);
  xyz_face[1] =  domaine_vf.xv(num_face,1);
  P[0] = domaine_vf.xp(elem, 0);
  P[1] = domaine_vf.xp(elem, 1);
  if (Objet_U::dimension == 3)
    {
      xyz_face[2] =  domaine_vf.xv(num_face,2);
      P[2] = domaine_vf.xp(elem, 2);
    }

  for (int i=0; i<3; i++)
    d += (xyz_face[i] - P[i])*(xyz_face[i] - P[i]);
  d= sqrt(d);
  return d;
}

// Q_int = Q_meso*circum_dis is W
double Triple_Line_Model_FT_Disc::compute_Qint(const DoubleTab& in_out, const double theta_app_loc,
                                               const int num_wall_face, double& Q_meso) const
{
  double circum_dis = 1.;
  if (Objet_U::bidim_axi)
    {
      const double xl = in_out(0,0);
      const double xr = in_out(1,0);
      double radial_dis = (xr + xl)/2;
      const double angle_bidim_axi = Maillage_FT_Disc::angle_bidim_axi();
      circum_dis = angle_bidim_axi*radial_dis;
      Cerr << "[TCL: MESO:] FILLING LIST Qmeso [bidim_axi] circum_dis = " << circum_dis << finl;
    }

  const double yl = in_out(0,1);
  const double yr = in_out(1,1);
  double Twall = 0.;
  double flux = 0.;
  ref_eq_temp_.valeur().get_flux_and_Twall(num_wall_face, flux, Twall);

  const double ytop = std::fmax(yr,yl);
  const double ybot = std::fmin(yr,yl);
  // Meso region is between ym_ and ymeso_:

  // for the last meso cell, if Front do not cut the right face of cell, only part of Meso contri is considered
  // the corresponding flux will be underestimated
  // we should ONLY consider cells before this last one
  // To do this, here we put 0 for the this flux, it will then be eliminated in the subsequent part of the compute_tcl_flux.
  double ln_y = 0.;
  if (!est_egal(ytop, ymeso_))
    ln_y = log(std::fmin(ymeso_,std::fmax(ytop,ym_))/std::fmin(ymeso_,std::fmax(ybot,ym_)));

  assert(kl_cond_>0.);
//  Cerr << "ln_y = " << ln_y << " time_total = " << temps << " Theta_app_local = " << theta_app_loc
//       <<" delT = "<< Twall << " kl = "<< kl_cond_ << finl;
  Q_meso = kl_cond_*(Twall/theta_app_loc)*ln_y; // Twall here is (Wall temperature - saturation temperature)..unit of Q_meso is W/m

  double Q_int = Q_meso*circum_dis; // unit of Q_int is W
//  Cerr << "Q_meso = " << Q_meso << finl;
  return Q_int;
}

// in_out are real absolute coordinates
void Triple_Line_Model_FT_Disc::compute_approximate_interface_inout(const Domaine_VDF& zvdf, const int korient,
                                                                    const int elem, const int num_face,
                                                                    DoubleTab& in_out, FTd_vecteur3& normale, double& surface_tot) const
{
  const int dim = Objet_U::dimension;
  if (Objet_U::dimension == 3)
    {
      Cerr << "not coded yet " <<finl;
      Process::exit();
    }
  in_out.resize(2,dim);
  const Maillage_FT_Disc& maillage = ref_eq_interf_.valeur().maillage_interface();
  const IntTab& facettes = maillage.facettes();
  const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const DoubleTab& sommets = maillage.sommets();

  // Boucle sur les facettes qui traversent l'element elem :
  // Moyenne ponderee des normales aux facettes qui traversent l'element
  //normale = {0., 0., 0.};
  normale[0] = 0.;
  normale[1] = 0.;
  normale[2] = 0.;
  // Centre de gravite de l'intersection facettes/element
  double centre[3] = {0., 0., 0.};
  // Somme des poids
  surface_tot =0.;
  int index=intersections.index_elem()[elem];
  while (index >= 0)
    {
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
      const int fa7 = data.numero_facette_;
      const DoubleTab& normale_facettes = maillage.get_update_normale_facettes();
      const double surface_facette = surface_facettes[fa7];
      const double surf = data.fraction_surface_intersection_ * surface_facette;
      surface_tot +=surf;
      for (int i = 0; i< dim; i++)
        {
          normale[i] += surf * normale_facettes(fa7, i);
          // Calcul du centre de gravite de l'intersection facette/element
          double g_i = 0.; // Composante i de la coordonnee du centre de gravite
          for (int j = 0; j < dim; j++)
            {
              const int som   = facettes(fa7, j);
              const double coord = sommets(som, i);
              const double coeff = data.barycentre_[j];
              g_i += coord * coeff;
            }
          centre[i] += surf * g_i;
        }
      index = data.index_facette_suivante_;
    }
  if (surface_tot > 0.)
    {
      const double inverse_surface_tot = 1. / surface_tot;
      double norme = 0.;
      for (int j = 0; j < dim; j++)
        {
          norme += normale[j] * normale[j];
          centre[j] *= inverse_surface_tot;
        }
      if (norme > 0)
        {
          double i_norme = 1./sqrt(norme);
          for (int j = 0; j < dim; j++)
            {
              double n_j = normale[j] * i_norme; // normale normee
              normale[j] = n_j;
            }
        }
    }

  // L'interface moyenne passe donc par "centre" (P) et de normale unitaire "normale" n=(a,b,c).
  // The interfacial plane equation is :
  // a * (x-xP) + b*(y-yP) + c * (z-zP) = 0

  // The cell width is : zvdf.dim_elem(elem, k);
  // In the tangential direction:
  const double delta_tangential = zvdf.dim_elem(elem, 1-korient); // 1-korient is the direction // to wall in 2D.
  const double xc = zvdf.xp(elem,1-korient);              // Cell-center in the direction // to the wall in 2D.
  // So left and right faces are at :
  double xL = xc - 0.5*delta_tangential;
  double xR = xc + 0.5*delta_tangential;

  // In the wall normal direction:
  const double delta_normal = zvdf.dim_elem(elem, korient); // korient is the direction perp to wall in 2D.
  const double yc = zvdf.xp(elem,korient);
  const double ymin = yc - 0.5*delta_normal;
  const double ymax = yc + 0.5*delta_normal;

  // And the normal to the faces is such that nface[1-korient] = \pm 1.

  // So left and right faces are at :
  const double a = -normale[1-korient];
  const double b = normale[korient];

  const double yG = centre[korient];
  const double xG = centre[1-korient];
  double yL = 0.;
  double yR = 0.;
  if (std::fabs(b)< DMINFLOAT)
    {
      // vertical interface:
      yL = ymin;
      yR = ymax;
      xL = xG;
      xR = xG;
    }
  else
    {
      const double coeff_directeur = a/b;
      yL = yG + coeff_directeur * (xL - xG);
      yR = yG + coeff_directeur * (xR - xG);
      yL = std::min(std::max(ymin,yL),ymax);
      yR = std::min(std::max(ymin,yR),ymax);

      if (std::fabs(a)< DMINFLOAT)
        {
          // horizontal interface:
          // in and out are at left and right faces so we already have :
          // xL = xL; xR = xR;
        }
      else
        {
          xL = xG + (yL-yG)/coeff_directeur;
          xR = xG + (yR-yG)/coeff_directeur;
        }
    }

  // Fill in the table :
  in_out(0,0) = xL;
  in_out(0,1) = yL;
  in_out(1,0) = xR;
  in_out(1,1) = yR;
}

double Triple_Line_Model_FT_Disc::compute_local_cos_theta(const Parcours_interface& parcours, const int num_face, const FTd_vecteur3& norm_elem) const
{
  FTd_vecteur3 nface = {0., 0., 0.} ;
  // GF is the vector from the cell to the interface barycenters:
  //        FTd_vecteur3 GF = {bary_facettes_dans_elem[0]-xG, bary_facettes_dans_elem[1]-yG, bary_facettes_dans_elem[2]-zG} ;
  // The 3 zeros are useless parameters...
  parcours.calculer_normale_face_bord(num_face, 0.,0.,0., nface[0], nface[1], nface[2]) ;

  double nl = 0.;
  double nl_fac = 0.;
  double n_dot = 0.;
  for (int dir = 0; dir< Objet_U::dimension; dir++)
    {
      nl += pow(norm_elem[dir],2);
      nl_fac += pow(nface[dir],2);
      n_dot += norm_elem[dir]*nface[dir];
    }
  double cos_theta = n_dot/(sqrt(nl)*sqrt(nl_fac));
//  Cerr << " cos_theta = " << cos_theta << finl;
  double theta_app_loc = acos(-cos_theta);
  return theta_app_loc;
}


/*
double Triple_Line_Model_FT_Disc::get_update_AAA()
{
	if (tag_tcl_ == maillage.get_mesh_tag())
	{
		return AAA;
	}

	// update the tcl model :
	tag_tcl_ = maillage.get_mesh_tag();
	compute_AAA();
}
*/

// Search for the numero of the wall_face in a given direction
int wall_face_towards(const int iface, const int starting_elem, const int num_bord, const Domaine_VDF& zvdf)
{
  const IntTab& face_voisins = zvdf.face_voisins();
  const IntTab& elem_faces = zvdf.elem_faces();
  //const Front_VF& the_wall = zvdf.front_VF(num_bord);
  int num_face_wall = -1;
  int elem = starting_elem;
  // Conventions TRUST VDF :
  //static const int NUM_FACE_GAUCHE = 0;
  //static const int NUM_FACE_BAS = 1;
  //static const int NUM_FACE_HAUT = 3;
  //static const int NUM_FACE_DROITE = 2;
  int count = 0;
  int is_wall=0;
  while (not(is_wall))
    {
      num_face_wall = elem_faces(elem,iface);
      if (zvdf.face_numero_bord(num_face_wall) == num_bord)
        {
          // The face num_face_wall is actually on the same boundary (wall) as num_face
          is_wall = 1;
          break;
        }
      else
        {
          // We move to the next elem:
          elem = face_voisins(num_face_wall, 0) + face_voisins(num_face_wall, 1) - elem;
          count++;
        }
      if (count>30) //n_ext_meso_)
        {
          Cerr << "Looking too far and wall not found" << finl;
          Process::exit();
        }
      if (elem==-1)
        {
          Cerr << "We reached a boundary but not the expected wall?"
               << " border reached :" <<  zvdf.face_numero_bord(num_face_wall)
               << " is not within num_bord :" << num_bord << finl;
          Process::exit();
        }
    }
  return num_face_wall;
}

bool is_in_list(const ArrOfInt& list, const int elem)
{
  const int nb_mixed = list.size_array();
  int j=0;
  for(j=0; j<nb_mixed; j++)
    {
      const int elemj = list[j];
      if (elem == elemj)
        {
          // yes, then we hit the "continue" in the next if and the loop continues with the next element
          break;
        }
    }
  if (j!=nb_mixed)
    return true;

  return false;
}

// Arguments :
//    elems_with_CL_contrib : The list of elements containing either the TCL itself, or the meso domaine.
//                            In all of them, the flux_evap has to be modified.
//    num_faces : Corresponding wall faces to which the flux should be assigned.
//    mpoint_from_CL : corresponding value to be assigned to each cell (as a surface mass flux [kg/(m2s)] of total interface area within the cell)
//    Q_from_CL : corresponding value to be assigned to each cell (as a thermal flux Q [unit?])
//
// Note that the same elem may appear twice (or more when contribution has to be displaced...) in the list,
// once (or more) for the micro contribution, once for the meso.
// It is not a problem, it will be considered by "+=" later
//
// num_faces is a locally temporary storages that in the end, is copied back into the class attribute boundary_faces_.
void Triple_Line_Model_FT_Disc::compute_TCL_fluxes_in_all_boundary_cells(ArrOfInt& elems_with_CL_contrib,
                                                                         ArrOfInt& num_faces,
                                                                         ArrOfDouble& mpoint_from_CL, ArrOfDouble& Q_from_CL)
{
  const Navier_Stokes_FT_Disc& ns = ref_ns_.valeur();
  const double dt = ns.schema_temps().pas_de_temps();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ns.domaine_dis().valeur());
  const DoubleVect& volumes = domaine_vf.volumes();
  const Fluide_Diphasique& fluide = ns.fluide_diphasique();
  const Fluide_Incompressible& phase_0 = fluide.fluide_phase(0);
  const Fluide_Incompressible& phase_1 = fluide.fluide_phase(1);
  const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique().valeurs();
  const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique().valeurs();
  const double rho_phase_0 = tab_rho_phase_0(0,0);
  const double rho_phase_1 = tab_rho_phase_1(0,0);
  const double jump_inv_rho = 1./rho_phase_1 - 1./rho_phase_0;
  const double Lvap = fluide.chaleur_latente();
  //const double coef = jump_inv_rho/Lvap;
  // const int dim3 = Objet_U::dimension == 3;
  const Transport_Interfaces_FT_Disc& eq_transport = ref_eq_interf_.valeur();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const IntTab& facettes = maillage.facettes();
  const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
  const ArrOfInt& index_facette_elem = intersections.index_facette();
  const DoubleTab& sommets = maillage.sommets();
  //      REF(Transport_Interfaces_FT_Disc) & refeq_transport =
  //      variables_internes().ref_eq_interf_proprietes_fluide;
  //      const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const bool first_call_for_this_mesh = (tag_tcl_ != maillage.get_mesh_tag());
  double instant_Qmicro = 0., instant_Qmeso = 0.;
  if (first_call_for_this_mesh)
    {
      tag_tcl_ = maillage.get_mesh_tag(); // update the TCL tag to match the mesh that was used.
      // reset instantaneous values :
      instant_mmicro_evap_ = 0.;
      instant_Qmicro = 0.;
      instant_vmicro_evap_ = 0.;
      instant_mmeso_evap_ = 0.;
      instant_vmeso_evap_ = 0.;
      instant_Qmeso = 0.;
      // Update the integration time :
      integration_time_ += dt;
    }
  else
    {
      Cerr << "Now that elems_, mp_and  Q_ are stored as members to the TCL class, we shouldn't need to call this function "
           << "several times for the same timestep. Please check or contact support." << finl;
      Process::exit();
    }

  const Parcours_interface& parcours = eq_transport.parcours_interface();
  // const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();

  // We store the elements (eulerian) crossed by a contact line.
  // As the specific source term is introduced here, a contribution of mdot should not be included afterwards.

  elems_with_CL_contrib.resize_array(0);

  num_faces.resize_array(0);

  mpoint_from_CL.resize_array(0);

  Q_from_CL.resize_array(0);

  // List to be completed with the meso region:
  ArrOfInt list_meso_elems;
  ArrOfInt list_meso_faces;

  // List to be completed with the micro region:
  ArrOfInt list_micro_elems;
  ArrOfInt list_micro_faces;
  ArrOfInt list_micro_indexs;

  // list of Qtot W/m from TCL model for printing
  ArrOfDouble list_micro_fraction;

  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, ns.domaine_dis().valeur());
  const IntTab& face_voisins = zvdf.face_voisins();
  const IntVect& orientation = zvdf.orientation();
  // const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis);
  const IntTab& elem_faces = zvdf.elem_faces();




  // 1. First loop on contact line cells:
  // Micro cells: 1, at wall BC + 2, 0<indicatrice <1 3, height < ym
  {
    // const double temps = ns.schema_temps().temps_courant();
    Cerr <<  " ****************************** TCL ****************************** " << finl;
    Cerr << que_suis_je() << "::compute_TCL_fluxes_in_all_boundary_cells() searching TCL cells" <<finl;
    const Domaine_Cl_dis_base& zcldis = ns.domaine_Cl_dis().valeur();
    // get wall BC frontiere nb
    int num_bord = -1;
    int nbwall_found = 0;
    for (int i=0; i<zcldis.nb_cond_lim(); i++)
      {
        const Cond_lim& la_cl = zcldis.les_conditions_limites(i);
        const Nom& bc_name = la_cl.frontiere_dis().le_nom();
        // For each BC, we check its type to see if it's a wall:
        // BC for hydraulic equation
        bool is_wall = sub_type(Dirichlet_paroi_fixe,la_cl.valeur())
                       || sub_type(Dirichlet_paroi_defilante,la_cl.valeur());
        if (is_wall)
          {
            Cerr << "[TCL]: Wall-type BC found at " <<
                 bc_name <<finl;
            num_bord = i;
            nbwall_found = nbwall_found +1;
//             break;
          }
      }
    if (nbwall_found != 1)
      Process::exit(Nom(nbwall_found) + " wall-type BC found!!!!!  Check JDD, pls" );
    if (num_bord<0)
      Process::exit( "[TCL]: !!! NO WALL-TYPE BOUNDARY WAS FOUND IN THE DOMAINE, PLESEASE CHECK JDD" );

    const Frontiere& fr=zcldis.les_conditions_limites(num_bord).frontiere_dis().frontiere();
    const int nb_first_face = fr.num_premiere_face();
    const int nb_face = fr.nb_faces();
    const DoubleTab& indica = eq_transport.inconnue ().valeurs ();

    // get TCL cells
    for (int ii = 0; ii < nb_face; ii++)
      {
        const int num_face = ii + nb_first_face;
        const int elemi = face_voisins (num_face, 0)
                          + face_voisins (num_face, 1) + 1;
        // if (is_in_list(list_micro_elems,elemi))
        //  continue;

        bool is_interf = (indica(elemi,0)>0.) && (indica(elemi,0)<1.) ;
        if (is_interf)
          {
            int index=intersections.index_elem()[elemi]; // index > 0 if Front presents
            if (index < 0)
              Process::exit( "[TCL]: !!! INDICATTICE and FT are NOT CONSISTENT " );
            else
              Cerr << "[TCL] FOUND interface cell at #" << elemi <<" with face number #" << num_face << finl;

            const int korient = orientation(num_face); // 0:x 1:y or 2:z (gives the direction it is normal to the wall)
            const double xwall = zvdf.xv(num_face, korient);
            double nwall[3] = { 0., 0., 0. }; // Away from the wall by convention here
            (zvdf.xp(elemi, korient) > xwall) ? nwall[korient] = 1. : nwall[korient] = -1.;
            // "-nwall " because we want to go to the wall

            int a = (-nwall[korient] < 0) ? 0 : Objet_U::dimension;
            int b = korient; // if the normal is along y, korient=1 and we want 1 for x
            int iface = a + b;


            if (distri_first_facette())
              {
                while (index >= 0)       // loop over index inside an elem, loop effective when one facette has go through serveral mesh
                  {
                    const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                    const int fa7 = data.numero_facette_;
                    const int sommet0 = facettes(fa7, 0);
                    const int sommet1 = facettes(fa7, 1);

                    if ((maillage.sommet_ligne_contact(sommet0)) || (maillage.sommet_ligne_contact(sommet1)))
                      {
                        int indextmp=index_facette_elem[fa7];
                        while (indextmp >= 0)
                          {
                            const Intersections_Elem_Facettes_Data& datatmp = intersections.data_intersection(indextmp);
                            const int elem = datatmp.numero_element_;
                            if (!is_in_list(list_micro_elems,elem))
                              {
                                // Cerr << "[TCL: MICRO] FOUND MICRO cell at #" << elem <<" also part of first facette #" << finl;
                                list_micro_elems.append_array(elem);

                                // Cerr << "[TCL: MICRO]: SERCHING corresponding face number at wall, DIRECTION = " << iface << finl;
                                // Cerr << "[TCL: MICRO]: (2D case: 0 LEFT; 1 DOWN; 2 Right; 3 UP.)" << finl;
                                const int num_face_wall = wall_face_towards(iface, elem, num_bord, zvdf);
                                list_micro_faces.append_array(num_face_wall);
                                list_micro_indexs.append_array(indextmp);
                                // Cerr << "[TCL: MICRO]: Corresponding face number at wall " << num_face_wall << finl;
                                Cerr << "[TCL: MICRO]: elem " << elem << " | face number at wall " << num_face_wall << finl;
                              }
                            indextmp = datatmp.index_element_suivant_;
                          }

                        if (Objet_U::bidim_axi)
                          {
                            int num_som = maillage.sommet_ligne_contact(sommet0) ? sommet0: sommet1;
                            x_cl_ = sommets(num_som,0);
                            Cerr << "[TCL: MICRO] position of the CL is set to be " << x_cl_ << finl;
                          }
                      }
                    index = data.index_facette_suivante_;
                  }
              }
            else
              {

                while (index >= 0)       // loop over index inside an elem, loop effective when one facette has go through serveral mesh
                  {
                    const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                    const int fa7 = data.numero_facette_;
                    const int sommet0 = facettes(fa7, 0);
                    const int sommet1 = facettes(fa7, 1);

                    if ((maillage.sommet_ligne_contact(sommet0)) || (maillage.sommet_ligne_contact(sommet1)))
                      {
                        int indextmp=index_facette_elem[fa7];
                        if (indextmp >= 0)
                          {
                            const Intersections_Elem_Facettes_Data& datatmp = intersections.data_intersection(indextmp);
                            const int elem = datatmp.numero_element_;
                            if (!is_in_list(list_micro_elems,elem))
                              {
                                // Cerr << "[TCL: MICRO] FOUND MICRO cell at #" << elem <<" also part of first facette #" << finl;
                                list_micro_elems.append_array(elem);

                                // Cerr << "[TCL: MICRO]: SERCHING corresponding face number at wall, DIRECTION = " << iface << finl;
                                // Cerr << "[TCL: MICRO]: (2D case: 0 LEFT; 1 DOWN; 2 Right; 3 UP.)" << finl;
                                const int num_face_wall = wall_face_towards(iface, elem, num_bord, zvdf);
                                list_micro_faces.append_array(num_face_wall);
                                list_micro_indexs.append_array(indextmp);
                                // Cerr << "[TCL: MICRO]: Corresponding face number at wall " << num_face_wall << finl;
                                Cerr << "[TCL: MICRO]: elem " << elem << " | face number at wall " << num_face_wall << finl;
                              }
                          }

                        if (Objet_U::bidim_axi)
                          {
                            int num_som = maillage.sommet_ligne_contact(sommet0) ? sommet0: sommet1;
                            x_cl_ = sommets(num_som,0);
                            Cerr << "[TCL: MICRO] position of the CL is set to be " << x_cl_ << finl;
                          }
                      }
                    index = data.index_facette_suivante_;

                  }



              }

            Cerr << "[TCL] Check if this CELL A MICRO and/or MESO cell by h" << finl;

            const double dist1 = std::fabs(zvdf.dist_face_elem0(num_face,elemi)); // h
            const double half_cell_height1 = std::fabs(zvdf.dist_face_elem0(elem_faces(elemi,iface),elemi)); // half distanct to furface parallel to wall
            const double hcell_top = dist1+half_cell_height1;
            const double hcell_bot = dist1-half_cell_height1;

            if (hcell_top <= ym_ + Objet_U::precision_geom && !is_in_list(list_micro_elems,elemi))
              {
                if (!is_in_list(list_micro_elems,elemi) && !distri_first_facette())
                  {

                    list_micro_elems.append_array(elemi);
                    list_micro_faces.append_array(num_face);
                    Cerr << "[TCL: MICRO]: elem " << elemi << " | face number at wall " << num_face << finl;
                  }
              }
            else if (hcell_bot <= ym_ + Objet_U::precision_geom)
              {
                // Cerr << "[TCL: MICRO+MESO] FOUND cell can contribue both MICRO and MESO #" << elemi << finl;
                if (!distri_first_facette() && !is_in_list(list_micro_elems,elemi))
                  {
                    list_micro_elems.append_array(elemi);
                    list_micro_faces.append_array(num_face);
                  }

                if (!is_in_list(list_meso_elems,elemi))
                  {
                    list_meso_faces.append_array(num_face);
                    list_meso_elems.append_array(elemi);
                  }

                Cerr << "[TCL: MICRO+MESO]: elem " << elemi << " | face number at wall " << num_face << finl;

              }
            else if (hcell_bot <= ymeso_ + Objet_U::precision_geom && !is_in_list(list_meso_elems,elemi) )
              {
                list_meso_elems.append_array(elemi);
                list_meso_faces.append_array(num_face);
                Cerr << "[TCL: MESO]: elem " << elemi << " | face number at wall " << num_face << finl;
              }



            // Cerr << "[TCL] Searching MICRO CELLS (INTERFACE + 1.e-10m<Y<YM) among Neighboorhood cells of iElem #" << elemi << finl;
            // Cerr << "[TCL] Searching MESO CELLS (INTERFACE + YM<Y<YMeso) among Neighboorhood cells of iElem #" << elemi << finl;
            // loop to find all cells with interface
            Cerr << "[TCL] Searching neighboor cells from Elem #" << elemi << " and face #" << num_face << finl;

            ArrOfInt future_new_elems;

            future_new_elems.append_array(elemi); // Initialize the list with the current micro cells
            while (future_new_elems.size_array()) // There are some elements to deal with
              {
                ArrOfInt new_elems(future_new_elems);
                future_new_elems.resize_array(0); // Empty the future list to be constructed
                for (int i=0; i< new_elems.size_array(); i++)
                  {
                    const int ielem = new_elems[i];
                    // Cerr << "[TCL ] Searching MICRO (INTERFACE + 1.e-10m<Y<YM)  or MESO (INTERFACE + YM<Y<YMeso) cells among Neighboorhood of iElem #" << ielem << finl;
                    const int nb_voisins = elem_faces.dimension(1);
                    for (int iv=0; iv < nb_voisins; iv++)
                      {
                        const int num_face2 = elem_faces(ielem,iv);
                        int elem_voisin = face_voisins(num_face2, 0) + face_voisins(num_face2, 1) - ielem;
                        if (elem_voisin<0)
                          continue; // We hit a border, there is no neighbour here
                        const int nbelemt= intersections.index_elem().size_array();
                        // set to -1 if the next elem is out of range, i.e., elem_voisin is a virtual element (paralisation)
                        int index2= (elem_voisin < nbelemt) ? intersections.index_elem()[elem_voisin]:-1;
                        if (index2 < 0)
                          continue; // This element is not crossed by the interface. Go to the next face.
                        // We won't want to do ielem twice (unless it was a micro elem)
                        // Have we seen this element already in the TCL list?
                        if (is_in_list(list_micro_elems,elem_voisin) && is_in_list(list_meso_elems,elem_voisin))
                          continue;

                        // Cerr << "[TCL ] FOUND interface cell #" << elem_voisin <<" among Neighboorhood cells of iElem #" << ielem << finl;
                        // Here, we are with a new element elem_voisin to append to both lists


                        const double dist = std::fabs(zvdf.dist_face_elem0(num_face,elem_voisin)); // h
                        const double half_cell_height = std::fabs(zvdf.dist_face_elem0(elem_faces(elem_voisin,iface),elem_voisin)); // half distanct to furface parallel to wall
                        const double hcell_dn = dist-half_cell_height;
                        const double hcell_up = dist+half_cell_height;
                        if (hcell_up <= ym_+Objet_U::precision_geom)
                          {
                            if (!distri_first_facette() && !is_in_list(list_micro_elems,elem_voisin))
                              {
                                // Cerr << "[TCL: MICRO] FOUND MICRO cell at #" << elem_voisin <<" among Neighboorhood cells of iElem #" << ielem << finl;
                                list_micro_elems.append_array(elem_voisin);
                                // Cerr << "[TCL: MICRO]: SERCHING corresponding face number at wall, DIRECTION = " << iface << finl;
                                // Cerr << "[TCL: MICRO]: (2D case: 0 LEFT; 1 DOWN; 2 Right; 3 UP.)" << finl;
                                const int num_face_wall = wall_face_towards(iface, elem_voisin, num_bord, zvdf);
                                list_micro_faces.append_array(num_face_wall);
                                // Cerr << "[TCL: MICRO]: Corresponding face number at wall " << num_face_wall << finl;
                                Cerr << "[TCL: MICRO]: elem " << elem_voisin << " | face number at wall " << num_face_wall << finl;
                                future_new_elems.append_array(elem_voisin);
                              }
                          }
                        else if (hcell_dn <= ym_+Objet_U::precision_geom)
                          {
                            // Cerr << "[TCL: MICRO+MESO] FOUND cell can contribue both MICRO and MESO #" << elem_voisin <<" among Neighboorhood cells of iElem #" << ielem << finl;
                            // Cerr << "[TCL: MICRO+MESO]: SERCHING corresponding face number at wall, DIRECTION = " << iface << finl;
                            // Cerr << "[TCL: MICRO+MESO]: (2D case: 0 LEFT; 1 DOWN; 2 Right; 3 UP.)" << finl;
                            const int num_face_wall = wall_face_towards(iface, elem_voisin, num_bord, zvdf);


                            if (!distri_first_facette() && !is_in_list(list_micro_elems,elem_voisin))
                              {
                                list_micro_elems.append_array(elem_voisin);
                                list_micro_faces.append_array(num_face_wall);
                                Cerr << "[TCL: MICRO]: elem " << elem_voisin << " | face number at wall " << num_face_wall << finl;
                                future_new_elems.append_array(elem_voisin);
                              }

                            if (!is_in_list(list_meso_elems,elem_voisin))
                              {
                                list_meso_elems.append_array(elem_voisin);
                                list_meso_faces.append_array(num_face_wall);
                                Cerr << "[TCL: MESO]: elem " << elem_voisin << " | face number at wall " << num_face_wall << finl;
                                if(!is_in_list(future_new_elems,elem_voisin))
                                  future_new_elems.append_array(elem_voisin);
                              }
                            // Cerr << "[TCL: MICRO+MESO]: Corresponding face number at wall " << num_face_wall << finl;

                          }
                        else if (hcell_dn<=ymeso_+Objet_U::precision_geom)
                          {
                            if (!is_in_list(list_meso_elems,elem_voisin))
                              {
                                // Cerr << "[TCL: MESO] FOUND  MESO cell  at #" << elem_voisin <<" among Neighboorhood cells of iElem #" << ielem << finl;
                                list_meso_elems.append_array(elem_voisin);
                                // Cerr << "[TCL: MESO]: SERCHING corresponding face number at wall, DIRECTION = " << iface << finl;
                                // Cerr << "[TCL: MESO]: (2D case: 0 LEFT; 1 DOWN; 2 Right; 3 UP.)" << finl;
                                const int num_face_wall = wall_face_towards(iface, elem_voisin, num_bord, zvdf);
                                list_meso_faces.append_array(num_face_wall);
                                Cerr << "[TCL: MESO]: elem " << elem_voisin << " | face number at wall " << num_face_wall << finl;
                                future_new_elems.append_array(elem_voisin);
                              }
                          }
                        else
                          {
                            // Cerr << "[TCL] HOWEVER ## elem " << elem_voisin << " NOT part of MICRO or Macro Region" << finl;
                            break;
                          }

                      }
                  }
              }
          }
      } // ==========================fin==============loop 1===================

    // Cerr << "[TCL: MICRO]: list_micro_elems number of elems: " << list_micro_elems << finl;
    // Cerr << "[TCL: MICRO]: list_micro_faces number of faces: " << list_micro_faces << finl;
    // Cerr <<  " *************************************** " << finl;
    // Cerr << "[TCL: MESO]: list_meso_elems number of elems: " << list_meso_elems << finl;
    // Cerr << "[TCL: MESO]: list_meso_faces number of faces: " << list_meso_faces << finl;


    // 2. send loop to fill-in micro contribution
    // We have our list of cells, we can complete their contribution and add them to the list
    {
      const double nn = list_micro_elems.size_array();
      for (int idx=0; idx<nn; idx++)
        {
          const int elem = list_micro_elems[idx];
          const int num_face = list_micro_faces[idx];

          double surface_tot = 0.;
          DoubleTab in_out(2,Objet_U::dimension); // The coords of left and right points where the interface cross the cell boundary
          FTd_vecteur3 norm_elem = {0., 0., 0.};
          if (inout_method_ == EXACT)
            {
              get_in_out_coords(zvdf, elem, dt, in_out, norm_elem, surface_tot);
            }
          else if (inout_method_ == APPROX)
            {
              const int korient= zvdf.orientation(num_face);
              compute_approximate_interface_inout(zvdf, korient,
                                                  elem, num_face,
                                                  in_out, norm_elem, surface_tot);
            }
          else
            {
              // both methods are used and compared :
              get_in_out_coords(zvdf, elem, dt, in_out, norm_elem, surface_tot);

              DoubleTab in_out_approx(2,Objet_U::dimension);
              FTd_vecteur3 norm_elem_approx = {0., 0., 0.};
              double surface_tot_approx = 0.;
              const int korient= zvdf.orientation(num_face);
              compute_approximate_interface_inout(zvdf, korient,
                                                  elem, num_face,
                                                  in_out_approx, norm_elem_approx, surface_tot_approx);
              Cerr << "comparison of approx to exact method in elem " << elem << finl;
              // in_out(0,0) -> x_in
              // in_out(0,1) -> y_in
              // in_out(1,0) -> x_out
              // in_out(1,1) -> y_out
              // Sorting points in/out with point 'in' closer to origin than 'ou' for the comparison:
              const double d0 = std::sqrt(in_out(0,0)*in_out(0,0)+in_out(0,1)*in_out(0,1));
              const double d1 = std::sqrt(in_out(1,0)*in_out(1,0)+in_out(1,1)*in_out(1,1));
              const double x_in = (d0<d1) ? in_out(0,0) : in_out(1,0);
              const double y_in = (d0<d1) ? in_out(0,1) : in_out(1,1);
              const double x_ou = (d0<d1) ? in_out(1,0) : in_out(0,0);
              const double y_ou = (d0<d1) ? in_out(1,1) : in_out(0,1);
              // approx points
              const double da0 = std::sqrt(in_out_approx(0,0)*in_out_approx(0,0)+in_out_approx(0,1)*in_out_approx(0,1));
              const double da1 = std::sqrt(in_out_approx(1,0)*in_out_approx(1,0)+in_out_approx(1,1)*in_out_approx(1,1));
              const double xa_in = (da0<da1) ? in_out_approx(0,0) : in_out_approx(1,0);
              const double ya_in = (da0<da1) ? in_out_approx(0,1) : in_out_approx(1,1);
              const double xa_ou = (da0<da1) ? in_out_approx(1,0) : in_out_approx(0,0);
              const double ya_ou = (da0<da1) ? in_out_approx(1,1) : in_out_approx(0,1);
              if ((std::fabs(x_in-xa_in)>DMINFLOAT) or (std::fabs(x_ou-xa_ou)>DMINFLOAT))
                Cerr << "inout x -- exact "
                     << x_in << " "
                     << x_ou << " approx "
                     << xa_in << " "
                     << xa_ou << " delta "
                     << x_in-xa_in << " "
                     << x_ou-xa_ou  << finl;
              if ((std::fabs(y_in-ya_in)>DMINFLOAT) or (std::fabs(y_ou-ya_ou)>DMINFLOAT))
                Cerr << "inout y -- exact "
                     << y_in << " "
                     << y_ou << " approx "
                     << ya_in << " "
                     << ya_ou << " delta "
                     << y_in-ya_in << " "
                     << y_ou-ya_ou  << finl;
              if ((std::fabs(norm_elem[0]-norm_elem_approx[0])>DMINFLOAT) or
                  (std::fabs(norm_elem[1]-norm_elem_approx[1])>DMINFLOAT))
                Cerr << "normal_elem -- exact " << norm_elem[0] << " " << norm_elem[1]
                     << " approx " << norm_elem_approx[0] << " " << norm_elem_approx[1]
                     << " delta " << norm_elem[0]-norm_elem_approx[0] << " " << norm_elem[1]-norm_elem_approx[1]
                     << finl;
              if (std::fabs(surface_tot-surface_tot_approx)>DMINFLOAT)
                Cerr << "surface -- exact "
                     << surface_tot << " approx "
                     << surface_tot_approx  << " delta "
                     << surface_tot-surface_tot_approx << finl;
              Cerr << "End of comparison" << finl;
            }


          // The offset distance should be accounted for because the methods in_out counts the origin at the cell
          {
            // in_out are real absolute coordinates.
            // We should make the "y" a normal distance to TCL
            const int korient = zvdf.orientation(num_face);
            const double xwall =zvdf.xv(num_face, korient);
            in_out(0,korient) -= xwall;
            in_out(1,korient) -= xwall;
            //in_out(0,1) += dist;
            //in_out(1,1) += dist;
            if (zvdf.dist_face_elem0(num_face,elem)>0)
              {
                // The face is on top of the element, so distance will be negative.
                // We want to count it positive by convention so :
                in_out(0,korient) *= -1;
                in_out(1,korient) *= -1;
              }
            assert(in_out(0,korient) >=-Objet_U::precision_geom);
            assert(in_out(1,korient) >=-Objet_U::precision_geom);
          }
          // const double xl = in_out(0,0);
          const double yl = in_out(0,1);
          // const double xr = in_out(1,0);
          const double yr = in_out(1,1);
          // Cerr << "xl,yl, xr,yr = " << xl << " " <<  yl << " " <<  xr << " " <<  yr << finl;



          double circum_dis = 1.;
          if (Objet_U::bidim_axi)
            {
              const double angle_bidim_axi = Maillage_FT_Disc::angle_bidim_axi();
              circum_dis = angle_bidim_axi*x_cl_;
              Cerr << "[TCL: MICRO] FILLING LIST MICRO [bidim_axi] circum_dis = " << circum_dis << finl;
            }


          const double ytop = std::fmax(yr,yl);
          const double ybot = std::fmin(yr,yl);

          // Micro region is between 1.e-10 and ym_:
          const double ymin = exp(-19);

          double Qtot = get_Qtcl(num_face);
          double Qmicro = 0.;
          double fraction = 0.;

          if(!distri_first_facette())
            {
              Qmicro = log(std::fmin(ym_,std::fmax(ytop,ymin))/ std::fmin(ym_,std::fmax(ybot,ymin))) / (log(ym_) +19.)*Qtot;
              fraction = log(std::fmin(ym_,std::fmax(ytop,ymin))/ std::fmin(ym_,std::fmax(ybot,ymin))) / (log(ym_) +19.);
            }
          else
            {
              // const int index=intersections.index_elem()[elem]; // index > 0 if Front presents
              const int index= list_micro_indexs[idx];
              const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
              fraction = data.fraction_surface_intersection_;
              Qmicro = fraction*Qtot;
            }

          double Q_int = Qmicro*circum_dis; // unit of Q_int is W

          if (std::fabs(Q_int)<DMINFLOAT)
            continue;

          if (first_call_for_this_mesh)
            {
              // We update values and integrals only at the first call
              instant_mmicro_evap_ += (Q_int/Lvap)*dt;
              instant_vmicro_evap_ = instant_mmicro_evap_*jump_inv_rho; // '=' not '+=' because the '+=' is already in instant_mmicro_evap_
            }


          elems_with_CL_contrib.append_array(elem);
          num_faces.append_array(num_face);
          mpoint_from_CL.append_array(Q_int/(surface_tot*Lvap));
          Q_from_CL.append_array(Q_int);
          list_micro_fraction.append_array(fraction);


          instant_Qmicro += Qmicro;

          Cerr.precision(16);
          Cerr << "[TCL: MICRO]:" << " Instantaneous tcl-evaporation = "<< instant_vmicro_evap_ << finl;
          Cerr << "[TCL: MICRO] time  = " << integration_time_ << " Qmicro= " << Qmicro << " Qint " << Q_int << finl;
          Cerr.precision(7);


        }
    }

    instant_Qmicro = mp_sum(instant_Qmicro);

    if (elems_with_CL_contrib.size_array()>0)
      {
        Cerr << " nb_contact_line_contribution " << Q_from_CL.size_array() << finl;
        Cerr << " ********* VALUES OF MICRO CONTRIBUTION *************** " << finl;
        Cerr << "elem\t| face\t| Q\t\t| mp\t\t| fraction" << finl;
        for (int idx = 0; idx < Q_from_CL.size_array(); idx++)
          {
            Cerr << elems_with_CL_contrib[idx] << "\t| " <<  num_faces[idx] << "\t| "
                 << Q_from_CL[idx] << "\t| " <<   mpoint_from_CL[idx] << "\t| "
                 <<list_micro_fraction[idx] << finl;
          }
        Cerr << " ************************************************************* " << finl;
      }

// 3. Third loop to fill-in meso contribution
// We have our list of cells, we can complete their contribution and add them to the list
    {
      const double nn = list_meso_elems.size_array();
      for (int idx=0; idx<nn; idx++)
        {
          const int elem = list_meso_elems[idx];
          const int num_face = list_meso_faces[idx];


          double surface_tot = 0.;
          DoubleTab in_out(2,Objet_U::dimension); // The coords of left and right points where the interface cross the cell boundary
          FTd_vecteur3 norm_elem = {0., 0., 0.};
          if (inout_method_ == EXACT)
            {
              get_in_out_coords(zvdf, elem, dt, in_out, norm_elem, surface_tot);
            }
          else if (inout_method_ == APPROX)
            {
              const int korient= zvdf.orientation(num_face);
              compute_approximate_interface_inout(zvdf, korient,
                                                  elem, num_face,
                                                  in_out, norm_elem, surface_tot);
            }
          else
            {
              // both methods are used and compared :
              get_in_out_coords(zvdf, elem, dt, in_out, norm_elem, surface_tot);

              DoubleTab in_out_approx(2,Objet_U::dimension);
              FTd_vecteur3 norm_elem_approx = {0., 0., 0.};
              double surface_tot_approx = 0.;
              const int korient= zvdf.orientation(num_face);
              compute_approximate_interface_inout(zvdf, korient,
                                                  elem, num_face,
                                                  in_out_approx, norm_elem_approx, surface_tot_approx);
              Cerr << "comparison of approx to exact method in elem " << elem << finl;
              // in_out(0,0) -> x_in
              // in_out(0,1) -> y_in
              // in_out(1,0) -> x_out
              // in_out(1,1) -> y_out
              // Sorting points in/out with point 'in' closer to origin than 'ou' for the comparison:
              const double d0 = std::sqrt(in_out(0,0)*in_out(0,0)+in_out(0,1)*in_out(0,1));
              const double d1 = std::sqrt(in_out(1,0)*in_out(1,0)+in_out(1,1)*in_out(1,1));
              const double x_in = (d0<d1) ? in_out(0,0) : in_out(1,0);
              const double y_in = (d0<d1) ? in_out(0,1) : in_out(1,1);
              const double x_ou = (d0<d1) ? in_out(1,0) : in_out(0,0);
              const double y_ou = (d0<d1) ? in_out(1,1) : in_out(0,1);
              // approx points
              const double da0 = std::sqrt(in_out_approx(0,0)*in_out_approx(0,0)+in_out_approx(0,1)*in_out_approx(0,1));
              const double da1 = std::sqrt(in_out_approx(1,0)*in_out_approx(1,0)+in_out_approx(1,1)*in_out_approx(1,1));
              const double xa_in = (da0<da1) ? in_out_approx(0,0) : in_out_approx(1,0);
              const double ya_in = (da0<da1) ? in_out_approx(0,1) : in_out_approx(1,1);
              const double xa_ou = (da0<da1) ? in_out_approx(1,0) : in_out_approx(0,0);
              const double ya_ou = (da0<da1) ? in_out_approx(1,1) : in_out_approx(0,1);
              if ((std::fabs(x_in-xa_in)>DMINFLOAT) or (std::fabs(x_ou-xa_ou)>DMINFLOAT))
                Cerr << "inout x -- exact "
                     << x_in << " "
                     << x_ou << " approx "
                     << xa_in << " "
                     << xa_ou << " delta "
                     << x_in-xa_in << " "
                     << x_ou-xa_ou  << finl;
              if ((std::fabs(y_in-ya_in)>DMINFLOAT) or (std::fabs(y_ou-ya_ou)>DMINFLOAT))
                Cerr << "inout y -- exact "
                     << y_in << " "
                     << y_ou << " approx "
                     << ya_in << " "
                     << ya_ou << " delta "
                     << y_in-ya_in << " "
                     << y_ou-ya_ou  << finl;
              if ((std::fabs(norm_elem[0]-norm_elem_approx[0])>DMINFLOAT) or
                  (std::fabs(norm_elem[1]-norm_elem_approx[1])>DMINFLOAT))
                Cerr << "normal_elem -- exact " << norm_elem[0] << " " << norm_elem[1]
                     << " approx " << norm_elem_approx[0] << " " << norm_elem_approx[1]
                     << " delta " << norm_elem[0]-norm_elem_approx[0] << " " << norm_elem[1]-norm_elem_approx[1]
                     << finl;
              if (std::fabs(surface_tot-surface_tot_approx)>DMINFLOAT)
                Cerr << "surface -- exact "
                     << surface_tot << " approx "
                     << surface_tot_approx  << " delta "
                     << surface_tot-surface_tot_approx << finl;
              Cerr << "End of comparison" << finl;
            }
          double factor = surface_tot/pow(volumes[elem], 2./3.);
          if (factor < 1.e-6)
            {
              Cerr << "We are skipping cell " << elem << "because the interface in it has negligible surface!!" << finl;
              continue; // skip cells crossed by almost no interface.
            }

          // The offset distance should be accounted for because the methods in_out counts the origin at the cell
          {
            // in_out are real absolute coordinates.
            // We should make the "y" a normal distance to TCL
            const int korient = zvdf.orientation(num_face);
            const double xwall =zvdf.xv(num_face, korient);
            in_out(0,korient) -= xwall;
            in_out(1,korient) -= xwall;
            //in_out(0,1) += dist;
            //in_out(1,1) += dist;
            if (zvdf.dist_face_elem0(num_face,elem)>0)
              {
                // The face is on top of the element, so distance will be negative.
                // We want to count it positive by convention so :
                in_out(0,korient) *= -1;
                in_out(1,korient) *= -1;
              }
            assert(in_out(0,korient) >=-Objet_U::precision_geom);
            assert(in_out(1,korient) >=-Objet_U::precision_geom);
          }
          const double xl = in_out(0,0);
          const double yl = in_out(0,1);
          const double xr = in_out(1,0);
          const double yr = in_out(1,1);

          // Cerr << "xl,yl, xr,yr = " << xl << " " <<  yl << " " <<  xr << " " <<  yr << finl;
          double theta_app_loc = compute_local_cos_theta(parcours, num_face, norm_elem);
          double Q_meso = 0.;
          double Q_int = compute_Qint(in_out, theta_app_loc, num_face, Q_meso);
          //  const double value = jump_inv_rho*Q_int/Lvap;
          //  Cerr << "[TCL-meso] Local source in elem=["  << elem << "] with value= " << value << "[m2.s-1]" << finl;
          //  Cerr << "Qint= " << Q_int << finl;

          if (std::fabs(Q_int)<DMINFLOAT)
            continue;

          if (first_call_for_this_mesh)
            {
              // We update values and integrals only at the first call
              instant_mmeso_evap_ += (Q_int/Lvap)*dt;                   // total mass of liquid evaporated
              instant_vmeso_evap_ = (instant_mmeso_evap_*jump_inv_rho); // total volume of liquid evaporated
              integrated_vmeso_evap_ = instant_vmeso_evap_;
            }
          Cerr.precision(16);
          Cerr << "[TCL: MESO] time  = " << integration_time_ << " Qmeso= " << Q_meso << " Qint " << Q_int << finl;
          Cerr.precision(7);
          // Cerr << "time = " << integration_time_ << " instantaneous mass-meso-evaporation = " << instant_mmeso_evap_ << " instantaneous meso-evaporation = " << instant_vmeso_evap_ << finl;
          elems_with_CL_contrib.append_array(elem);
          //mpoint_from_CL.append_array((Q_meso*surface_tot/v)/Lvap); // VP: replaced Q_int by Q_meso....unit--kg/m^2.s
          double x1=0., y1=0.;
          double x2=0., y2=0.;
          // point 0 will be the left if (left lower than right)&&(M lower than left)
          // otherwise, it should be point M if M is in the elem
          double x0=0., y0=0.;
          if (yl<yr)
            {
              x1=xl;
              y1=yl;
              x2=xr;
              y2=yr;
            }
          else
            {
              x2=xl;
              y2=yl;
              x1=xr;
              y1=yr;
            }
          if (ym_<=y1)
            {
              // micro not in that cell:
              x0=x1;
              y0=y1;
            }
          else if(ym_<=y2)
            {
              // point M is in that cell:
              // ie PART of the interface is in micro.
              y0=ym_;
              const double f = (y0-y1)/(y2-y1);
              x0=x1+f*(x2-x1);
              Cerr << "[x0 vs x_cl] t= "<< integration_time_ << " elem= " << elem << " x0/x_cl " << x0<< " " <<x_cl_<< finl;
            }
          else
            {
              // That cell is fully micro ie ym>y2>y1
              // in that case, smeso should be 0. :
              x0=x2;
              y0=y2;
            }
          const double dx = x2-x0;
          const double dy = y2-y0;
          const double s_meso = sqrt(dx*dx+dy*dy);
          if (s_meso<DMINFLOAT)
            {
              // Qmeso is probably zero too
              mpoint_from_CL.append_array(0.);
            }
          else
            {
              mpoint_from_CL.append_array(Q_int/(surface_tot*Lvap)); // LW 2023.27.09 coherance with Qint used in Q_from_CL to compute a mean weighted by surface_tot
            }
          instant_Qmeso += Q_meso;
          {
            // comparison s_meso as dist_LR to the reconstruction based on surface_tot
            // Si la CL est dans la maille...
            if ((Objet_U::bidim_axi) && ((xr-x_cl_)*(xl-x_cl_)<=0.))
              {
                double S_micro = Maillage_FT_Disc::angle_bidim_axi()*x_cl_*sm_;
                const double rm = 0.5*(xr+xl);
                const double s_meso_assessed = (surface_tot-S_micro)/(Maillage_FT_Disc::angle_bidim_axi()*rm);
                Cerr.precision(16);
                Cerr << "[TCL: MESO: Assessment-S_meso] t= "<< integration_time_ << " elem= " << elem << " LR= " ;
                Cerr << s_meso <<  " Ameso/2pir= " << s_meso_assessed << " err(%)= " << 50.*std::abs(s_meso_assessed-s_meso)/(s_meso_assessed+s_meso) << finl;
                Cerr.precision(7);
              }
          }
          Q_from_CL.append_array(Q_int);
          num_faces.append_array(num_face);
          // instant_mmeso_evap_ += Q_int/Lvap;                      // total mass of liquid evaporated
          // instant_vmeso_evap_ += instant_mmeso_evap_/rho_phase_1; // total volume of liquid evaporated
        }
    }
    instant_Qmeso = mp_sum(instant_Qmeso);

    const double verbosity=1;
    if (verbosity)
      {
        const int nb_contact_line_contribution = elems_with_CL_contrib.size_array();
        if (nb_contact_line_contribution>0 )
          {
            const double t = ns.schema_temps().temps_courant();
            Cerr << " time " << t << finl;
            Cerr << " nb_contact_line_contribution " << nb_contact_line_contribution << finl;
            Cerr << " ********* VALUES OF CONTACT LINE CONTRIBUTION *************** " << finl;
            Cerr << "elem\t| face\t| Q\t\t| mp" << finl;
            for (int idx = 0; idx < nb_contact_line_contribution; idx++)
              {
                Cerr << elems_with_CL_contrib[idx] << "\t| " <<  num_faces[idx] << "\t| "
                     << Q_from_CL[idx] << "\t| " << mpoint_from_CL[idx]  << finl;
              }
            Cerr << " ************************************************************* " << finl;
            Cerr << "[instant-Q] time micro meso " << t << " " <<
                 instant_Qmicro << " " << instant_Qmeso << finl;
            Cerr << " ************************************************************* " << finl;
          }
      }

    // double micro_mevap = instant_mmicro_evap_ + instant_mmeso_evap_;
    // You want the total value on all processors:
    // micro_mevap = Process::mp_sum(micro_mevap);

    // Temporary code : Later, we will either (i) fill elems_ directly without using a temporary elems_with_CL_contrib
    // or (ii) make a list without repetition (summing up values for micro and meso together)
    // Store :

    elems_.resize_array(0); // to empty the list
    boundary_faces_.resize_array(0);
    mp_.resize_array(0);
    Q_.resize_array(0);
    const int nn = elems_with_CL_contrib.size_array();
    assert(num_faces.size_array() == nn);
    assert(mpoint_from_CL.size_array() == nn);
    assert(Q_from_CL.size_array() == nn);
    for (int idx = 0; idx < nn; idx++)
      {
        const int elem = elems_with_CL_contrib[idx];
        elems_.append_array(elem);
        boundary_faces_.append_array(num_faces[idx]);
        mp_.append_array(mpoint_from_CL[idx]);
        Q_.append_array(Q_from_CL[idx]);
      }
    // GB 2023.04.07: a time integral. Formulae modified but not validated.
    const double collect_vevap = Process::mp_sum(instant_vmicro_evap_ + instant_vmeso_evap_)*dt;
    if (Process::je_suis_maitre())
      {
        vevap_int_ +=collect_vevap;
        Cerr << "[TCL: MICRO+MESO] vevap_int_(micro+meso) = " << vevap_int_ << " time = " << integration_time_ << finl;
      }
  }
}

// correct mpoint with values from TCL model.
void Triple_Line_Model_FT_Disc::corriger_mpoint(DoubleTab& mpoint) const
{
  const int nb_contact_line_contribution = elems_.size_array();
  Cerr << " nb_contact_line_contribution #" << nb_contact_line_contribution << finl;

  // I believe it is required to put zero in all cells first to remove any possible macroscopic contribution.
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      mpoint(elem) = 0.;
    }

  // Then, we can sum our micro and meso contributions :
  // Unit: kg/(m2s)
  // The area on which mp_ should be applied is by definition the total interface area within that cell.
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      mpoint(elem) += mp_[idx];
      Cerr << "[TCL] contrib added to mpoint for cell " << elem << " with value= " << mp_[idx] << finl;
    }
}

// correct secmem2 with values from TCL model.
void Triple_Line_Model_FT_Disc::corriger_secmem(const double coef, DoubleTab& secmem2) const
{
  const int nb_contact_line_contribution = elems_.size_array();
  // const double max_val_before = max_array(secmem2);
  // const double min_val_before = min_array(secmem2);
  // I believe it is required to put zero in all cells with TCL model first to remove any possible macroscopic contribution.
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      secmem2(elem) =0.;
    }

  // Then, we can sum our micro and meso contributions :
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      const double Q_val = Q_[idx];
      const double value = coef*Q_val;
      secmem2(elem) += value; // '+=' because several contribs are in the same cell element.
    }
  // Cerr << "[secmem2] max before/after TCL model: " << max_val_before << " / " << max_array(secmem2) << finl;
  // Cerr << "[secmem2] min before/after TCL model: " << min_val_before << " / " << min_array(secmem2) << finl;
}

// In cells where the TCL model is activated, the mean liquid temperature in each cell should be
// computed in agreement with the meso-scale model. Based on the quasi-steady solution of diffusion in
// a wedge, the temperature is assumed to depend only on the local angle (in polar coordinates) as :
//    T = Twall + DeltaT * theta/alpha,
//    where DeltaT = Tsat-Twall, alpha is the wedge angle (current interface slope at local coord 's'
//    (ie within the current cell), and theta is the angular coordinate.
// Based on this temperature profile, the mean value in cell 'j' can be readily calculated as :
//    <T>(j) = Twall + 0.5 * DeltaT
//
void Triple_Line_Model_FT_Disc::set_wall_adjacent_temperature_according_to_TCL_model(DoubleTab& temperature) const
{
  for (int idx = 0; idx < elems_.size_array(); idx++)
    {
      const int elem = elems_[idx];
      const int num_face = boundary_faces_[idx];
      const double Twall =ref_eq_temp_.valeur().get_Twall_at_face(num_face);
      const double DeltaT = TSAT_CONSTANTE-Twall;
      const double Tbefore =  temperature(elem);
      temperature(elem) =Twall + 0.5*DeltaT;
      Cerr << "[Temperature correction] elem= " << elem << " Tbefore/Tafter: " << Tbefore << " / " << temperature(elem) << finl;
    }
}

// In the TCL region, each cell recieve a wall-flux phi_w assumed in direct transfer with the interface, that is:
//    Qwall = A_w phi_w = A_i phi_i = Qinterf = Qmicro + Qmeso.
// Besides, internal (diffusive) fluxes are set to zero.
// And the cell mean-temperature is set according to the model to :
//    <T>(j) = Twall + 0.5 * DeltaT
// Consequently, if either (Twall is time-dependent) or (a cell gets in/out of the TCL region), or (there is a
// convective flux at the matching face), then there is an imbalance between the time evolution of the internal
// energy in the TCL region and the fluxes.
// Several options :
//   (i)   We store it to correct the interfacial flux (for instance the micro) at the next time-step?
//   (ii)  The disequilibrium is used to modify the wall flux (complicated, as fluxes have already been computed)?
//   (iii) We correct the closest fully liquid neighbor with this energy difference.
//
// Option (iii) seems the simplest, so we go for it.
double compute_global_TCL_energy_disequilibrium()
{
  return 0.;
}

void Triple_Line_Model_FT_Disc::correct_TCL_energy_evolution(DoubleTab& temperature) const
{
  if (TCL_energy_correction_ == 0)
    return;

  // In Joules :
  const double disequilibrium = compute_global_TCL_energy_disequilibrium();
  const int elem = 0; //closest_liquid_neighbour has to be determined;
  Cerr << "Code unfinished in Triple_Line_Model_FT_Disc::correct_TCL_energy_evolution" << finl;
  Process::exit();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ref_eq_temp_.valeur().domaine_dis().valeur());
  const double vol = domaine_vf.volumes(elem);
  // Minus sign because we want to compensate the disequilibrium
  const double deltaT = -disequilibrium/(rhocpl_*vol);
  temperature(elem) += deltaT;
}

void Triple_Line_Model_FT_Disc::correct_wall_adjacent_temperature_according_to_TCL_fluxes(DoubleTab& temperature) const
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ref_eq_temp_.valeur().domaine_dis().valeur());
  const int nb_contact_line_contribution = elems_.size_array();
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      const int num_face = boundary_faces_[idx];
      double flux=0., Twall=0.;
      ref_eq_temp_.valeur().get_flux_and_Twall(num_face,
                                               flux, Twall);
      temperature(elem) =Twall;
    }
  // MEMO : remove the previous loop if num_faces is computed here.
  // Then, we can sum our micro and meso contributions :
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      const int num_face = boundary_faces_[idx];
      // Get the distance between the center of the elem and the given face:
      const double d = compute_distance(domaine_vf, num_face, elem);
      double flux=0., Twall=0.;
      ref_eq_temp_.valeur().get_flux_and_Twall(num_face, flux, Twall);
      // It is the total flux, so It should simply be :
      temperature(elem) = Twall - flux*d/kl_cond_; // Can be done several times, no problem.
    }
}

double Triple_Line_Model_FT_Disc:: get_Qtcl(const int num_face)
{
//   return Qtcl_;
  //  Qtcl_ = coeffb_ + coeffa_*(theta_app_*22/(7*180));
//   Cerr << "coeffb " << coeffb_ << " coeffa " << coeffa_ << "theta " << (theta_app_) << finl;
  //  return (coeffb_ + coeffa_*(theta_app_))*(6.2/(theta_app_));
  if (read_via_file_)
    {
      const double Twall = ref_eq_temp_.valeur ().get_Twall_at_face (
                             num_face);

      const int num_col = num_colon_qtcl_; // colon number corresponding to Qtcl
      const int num_tem = num_colon_tem_;

      int ind = 0;
      while (ind < tab_Mtcl_.dimension (1) && tab_Mtcl_ (num_tem, ind) < Twall)
        {
          ++ind;
        }
      if (ind == 0)
        {
          if (!est_egal (tab_Mtcl_ (num_tem, ind), 0))
            Qtcl_ =  tab_Mtcl_ (num_col, ind) *Twall / tab_Mtcl_ (num_tem, ind);
          else
            Qtcl_ = 0. ;
        }
      else if (ind == tab_Mtcl_.dimension (1))
        {
          Qtcl_ = tab_Mtcl_ (num_col, ind - 1);
        }
      else
        {
          double x1 = tab_Mtcl_ (num_tem, ind - 1);
          double x2 = tab_Mtcl_ (num_tem, ind);
          double y1 = tab_Mtcl_ (num_col, ind - 1);
          double y2 = tab_Mtcl_ (num_col, ind);
          assert(!est_egal (x1, x2));
          Qtcl_ = y1 + (y2 - y1) * (Twall - x1) / (x2 - x1);
        }
    }
  return Qtcl_;
}

double Triple_Line_Model_FT_Disc:: get_theta_app(const int num_face)
{
  if (read_via_file_)
    {

      double Twall = 0.;
      const Domaine_Cl_dis_base& zcldis = ref_cast(
                                            Domaine_Cl_dis_base, ref_eq_temp_->domaine_Cl_dis ().valeur ());
      const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, zcldis);
      const Cond_lim& la_cl = zclvdf.la_cl_de_la_face (num_face);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      if (sub_type(Echange_contact_VDF_FT_Disc, la_cl.valeur ()))
        {
          const Echange_contact_VDF_FT_Disc& la_cl_typee = ref_cast(
                                                             Echange_contact_VDF_FT_Disc, la_cl.valeur ());
          const double T_imp = la_cl_typee.Ti_wall (num_face - ndeb);
          Twall = T_imp;
        }
      else if (sub_type(Echange_impose_base, la_cl.valeur ()))
        {
          const Echange_impose_base& la_cl_typee = ref_cast(Echange_impose_base,
                                                            la_cl.valeur ());
          const double T_imp = la_cl_typee.T_ext (num_face - ndeb);
          Twall = T_imp;
        }
      else
        {
          Cerr
              << "How can we set Twall when temperature(elem) is not valid? Or is it?"
              << finl;
          Process::exit ();
        }

      // const double Twall = ref_eq_temp_.valeur ().get_Twall_at_face (
      //                       num_face);

      const int num_col = num_colon_app_; // colon number corresponding to theta_app
      const int num_tem = num_colon_tem_;

      int ind = 0;
      while (ind < tab_Mtcl_.dimension (1) && tab_Mtcl_ (num_tem, ind) < Twall)
        {
          ++ind;
        }
      if (ind == 0)
        {
          theta_app_ = tab_Mtcl_ (num_col, 0);
        }
      else if (ind == tab_Mtcl_.dimension (1))
        {
          theta_app_ = tab_Mtcl_ (num_col, ind - 1);
        }
      else
        {
          double x1 = tab_Mtcl_ (num_tem, ind - 1);
          double x2 = tab_Mtcl_ (num_tem, ind);
          double y1 = tab_Mtcl_ (num_col, ind - 1);
          double y2 = tab_Mtcl_ (num_col, ind);
          assert(!est_egal (x1, x2));
          theta_app_ = y1 + (y2 - y1) * (Twall - x1) / (x2 - x1);
        }
      // if num <= 1, i.e., =0 or =1 we do not change the theta_app_:No numerical forcing breakup
      // else replaced by a big value for forcing numerical breakup
      if ( Rc_tcl_GridN() >= 1 && !reinjection_tcl())
        {
          const Navier_Stokes_FT_Disc& ns = ref_ns_.valeur ();
          const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, ns.domaine_dis ().valeur ());

          const IntVect& orientation = zvdf.orientation();
          const IntTab& elem_faces = zvdf.elem_faces();

          int elem_voisin = zvdf.face_voisins (num_face, 0) + zvdf.face_voisins (num_face, 1) + 1;
          int numfacex = -1;

          const int nb_voisins = elem_faces.dimension(1);
          for (int iv=0; iv < nb_voisins; iv++)
            {
              const int num_face2 = elem_faces(elem_voisin,iv);
              if (orientation(num_face2) == 0 )
                {
                  numfacex = num_face2;
                  break;
                }
            }
          assert(numfacex>=0);
          const double cell_size = 2. * std::fabs (zvdf.dist_face_elem0 (numfacex, elem_voisin));
          const double xwall_ = zvdf.xv(num_face, 0);
          if (xwall_ <= cell_size*Rc_tcl_GridN() )
            theta_app_ = thetaC_tcl();
        }
    }
  return theta_app_;
}


