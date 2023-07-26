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

static const double TSAT_CONSTANTE = 0.;

Implemente_instanciable_sans_constructeur( Triple_Line_Model_FT_Disc, "Triple_Line_Model_FT_Disc", Objet_U ) ;



Triple_Line_Model_FT_Disc::Triple_Line_Model_FT_Disc() :
  // Values initialized:
  activated_(false),
  capillary_effect_on_theta_activated_(0),
  TCL_energy_correction_(0),
  tag_tcl_(-1),
  coeffa_(0),
  coeffb_(0),
  Qtcl_(0),
  lv_(0.),
  theta_app_(0.),
  x_cl_(0.),
  ym_(0.),
  sm_(0.),
//  old_xcl_(0.),
  initial_CL_xcoord_(0.),
  kl_cond_(-1.),  // Invalid
  rhocpl_(-1.),  // Invalid
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
  activated_ = true;
  return is;
}

// Description: Methode appelee par readOn. Declaration des parametres a lire dans le .data
void Triple_Line_Model_FT_Disc::set_param(Param& p)
{
  p.ajouter("Qtcl", &Qtcl_);
  p.ajouter("lv", &lv_);
  p.ajouter("coeffa", &coeffa_);
  p.ajouter("coeffb", &coeffb_);
  p.ajouter("theta_app", &theta_app_);
  p.ajouter("ylim", &ym_);
  p.ajouter("ym", &ym_);
  p.ajouter("sm", &sm_,Param::REQUIRED);
  p.ajouter("initial_CL_xcoord", &initial_CL_xcoord_);
  p.ajouter_flag("enable_energy_correction", &TCL_energy_correction_);
  p.ajouter_flag("capillary_effect_on_theta", &capillary_effect_on_theta_activated_);
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
  elems_.set_smart_resize(1);
  elems_.resize_array(0);
  boundary_faces_.set_smart_resize(1);
  boundary_faces_.resize_array(0);
  mp_.set_smart_resize(1);
  mp_.resize_array(0);
  Q_.set_smart_resize(1);
  Q_.resize_array(0);
}


void Triple_Line_Model_FT_Disc::completer()
{
  // Via the temperature transport equation, we directly get access to a Fluide_Incompressible:
  const Milieu_base& milieu = ref_eq_temp_.valeur().milieu();
  kl_cond_ = milieu.conductivite()(0,0);
  rhocpl_ =  milieu.masse_volumique()(0,0) * milieu.capacite_calorifique()(0,0);

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
//
// in_out(0,0) -> x_in
// in_out(0,1) -> y_in
// in_out(1,0) -> x_out
// in_out(1,1) -> y_out
// norm_elem : the averaged normal into the element.
// surface_tot : The local interfacial area.
void Triple_Line_Model_FT_Disc::get_in_out_coords(const Domaine_VDF& zvdf, const int elem,
                                                  const double dt,
                                                  DoubleTab& in_out, FTd_vecteur3& norm_elem, double& surface_tot)
{
//	  const double& dt = schema_temps().pas_de_temps();
  const int dim = Objet_U::dimension;
  const int dim3 = Objet_U::dimension == 3;
  const Transport_Interfaces_FT_Disc& transport = ref_eq_interf_.valeur();
  const Maillage_FT_Disc& maillage = transport.maillage_interface();
//  const double temps = maillage.temps();
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
  double zP = 0.;
  if (dim3)
    zP = zvdf.xp(elem, 2);

  // We get the width of the elem:
  const double delta_x = zvdf.dim_elem(elem, 0);
  const double delta_y = zvdf.dim_elem(elem, 1);
  double delta_z = 0.;
  if (dim3)
    delta_z = zvdf.dim_elem(elem, 2);

  //x coordinates of the barycenter of the right and left faces:
  double rface_dis = xP + delta_x/2;
  double lface_dis = xP - delta_x/2;
  //y coordinates of the barycenter of the top and bottom faces:
  double tface_dis = yP + delta_y/2;
  double bface_dis = yP - delta_y/2;

  //    const double tol_dis = 1.0e-8; // tolerance-distance of node from the wall (not sure if we need it)
// Cerr << " rface_dis = "<< rface_dis << " lface_dis = " << lface_dis << finl;
  //double surface_tot = 0.;
  FTd_vecteur3 bary_facettes_dans_elem = {0., 0., 0.} ;
  FTd_vecteur3 norm = {0., 0., 0.}; // added vp
  //FTd_vecteur3 norm_elem = {0., 0., 0.};
  // Boucle sur les facettes qui traversent l'element elem :
  //       int sng_node = 0.;
  double x_sommet0 = 0.;
  double y_sommet0 = 0.;
  //     double x_sommetr = 0.;
  //     double y_sommetr = 0.;
  //     double x_sommetl = 0.;
  //     double y_sommetl = 0.;
  //     double x_node = 0.;
  //     double y_node = 0.;
  //    double nx_ref = 0.;
  //    double ny_ref = 0.;
  //double nx0 = 0.;
// double ny0 = 0.;
  //    double nxr = 0.;
  //    double nyr = 0.;
  //    double nxl = 0.;
  //    double nyl = 0.;
  double x_sommet1 = 0.;
// double y_sommet1 = 0.;
// double nx1 = 0.;
// double ny1 = 0.;
  double slope = 0.;
  double slopei = 0.;
  double yl = 0.;
  double yr = 0.;
  double xl = 0.;
  double xr = 0.;
  double ytest = 0.;
  double xtest = 0.;
  double xint[2] = {};
  double yint[2] = {};
  int nint = 0;
  const double eps = 1.e-10;
  ytest = delta_z*zP;  //dummy
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
          norm[dir] = nx;
        }

      slope = -norm[0]/norm[1]; //slope of the facette
      slopei = -norm[1]/norm[0];
      //    Cerr << "i am here" << finl;
      int inner_nodes = 0;
      int int_face = 0;

      for (int isom = 0; isom< Objet_U::dimension; isom++)
        {
          const int num_som = facettes(fa7, isom); // numero du sommet dans le tableau sommets
          // Cerr << " node-number #" << num_som << "face-number #" << fa7 << finl;

          //          const double sommet_d = sommets(num_som, isom);
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
          if(coord_sommet[0] < rface_dis && coord_sommet[0] > lface_dis &&
              coord_sommet[1] < tface_dis && std::fabs(coord_sommet[1]) >= bface_dis)
            {
              //  Cerr << " This node lies inside the element #" << elem << finl;
              inner_nodes += 1;
            }
        }
      Cerr << " no of inner nodes = " << inner_nodes << finl;

      if (inner_nodes >= 2)
        {
          // Cerr << " This facette lies totally inside the element #" << elem << finl;
          for (int isom = 0; isom< Objet_U::dimension; isom++)
            {
              const int node_number = facettes(fa7, isom); // numero du sommet dans le tableau sommets
              // Cerr << " node-number #" << node_number << "face-number #" << fa7 << finl;
              for (int dir = 0; dir < Objet_U::dimension; dir++)
                {
                  coord_sommet[dir] = sommets(node_number,dir);
                }
              //    if(coord_sommet[1] == bface_dis || coord_sommet[1] <= 1.e-10)
              if(std::abs(coord_sommet[1]-bface_dis) <= eps)
                {
                  xint[nint] = coord_sommet[0];
                  //   yint[nint] = coord_sommet[1];
                  yint[nint] = 0.;
                  //  Cerr << " time_base= " << temps << " xint[0]= " << xint[nint] << " yint[0]= " << yint[nint] << finl;
                  // nint += 1;
                  //   Cerr<< "nint = "<< nint << finl;
                  //     const double U_cl = (xint[nint]-initial_CL_xcoord_)/dt;
                  //      old_xcl_ = xint[nint];
                  //     Cerr << " time_cl= " << temps << " CL_velocity= " << U_cl << " xcl= " << old_xcl_ << " time-step= " << dt << finl;
                  int_face += 1;
                }

              x_sommet0 = coord_sommet[0];
              y_sommet0 = coord_sommet[1];
            }
          // Cerr<< "nint = "<< nint << finl;
        }
      if(inner_nodes == 1)
        {
          // Cerr<< "nint = "<< nint << finl;
          //  Cerr << " find the intersection points of the facettes and elem-faces "<< finl;
          for (int isom = 0; isom< Objet_U::dimension; isom++)
            {
              const int node_number = facettes(fa7, isom); // numero du sommet dans le tableau sommets
              //  Cerr << " node-number #" << node_number << "face-number #" << fa7 << finl;
              for (int dir = 0; dir < Objet_U::dimension; dir++)
                {
                  coord_sommet[dir] = sommets(node_number,dir);
                }
              //           if(coord_sommet[0] < rface_dis && coord_sommet[0] > lface_dis &&
              //             coord_sommet[1] < tface_dis && coord_sommet[1] >= bface_dis)
              if(coord_sommet[0] < rface_dis && coord_sommet[0] > lface_dis &&
                  coord_sommet[1] < tface_dis && std::fabs(coord_sommet[1]) >= bface_dis)
                {
                  // Cerr << " sng_node = " << sng_node << finl;
                  //   if(std::abs(coord_sommet[0]-rface_dis) > std::abs(coord_sommet[0]-lface_dis))
                  //  {
                  //    Cerr << " This node is closer to the left-face and therefore we check for "
                  //    		"the intersection with left or bottom face " << finl;
                  //     if(sng_node == 0.)

                  x_sommet0 = coord_sommet[0];
                  y_sommet0 = coord_sommet[1];
                  // nx0 = norm[0];
                  //  ny0 = norm[1];
                  //   if(coord_sommet[1] == bface_dis || coord_sommet[1] <= 1.e-10)
                  if(std::abs(coord_sommet[1]-bface_dis) <= eps)
                    {
                      xint[nint] = x_sommet0;
                      //     yint[nint] = y_sommet0;
                      yint[nint] = 0.;
                      //  Cerr << " time_base= " << temps << " xint[0]= " << xint[nint] << " yint[0]= " << yint[nint] << finl;
                      //   const double U_cl = (xint[nint]-old_xcl_)/dt;
                      //   Cerr << "xint_old = " << old_xcl_ << finl;
                      //   old_xcl_ = xint[nint];
                      //    Cerr << " time_cl= " << temps << " CL_velocity= " << U_cl << " xcl= " << old_xcl_ << finl;
                      //  Cerr << "time_cl= " << temps << " CL_velocity= " << U_cl << " time-step= " << dt << finl;
                      nint += 1;
                      int_face += 1;
                    }
                }
              else
                {
                  x_sommet1 = coord_sommet[0];
                  //   y_sommet1 = coord_sommet[1];
                  //  nx1 = norm[0];
                  //  ny1 = norm[1];
                }
            }
          //  Cerr << " x_sommet0 = " << x_sommet0 << " y_sommet0 = "
          //     << y_sommet0 << " normalx = " << nx0 << " normaly = " << ny0 << finl;

          //   Cerr << " x_sommet1 = " << x_sommet1 << " y_sommet1 = "
          //      << y_sommet1 << " normalx = " << nx1 << " normaly = " << ny1 << finl;
          // we look for the intersection of this facette with a face of the element
          if(x_sommet1 >= rface_dis) // intersection is either at right or at top
            {
              // Cerr << "x_sommet1 is greater than rface_dis" << finl;
              // test for top-face
              ytest = delta_y;
              xtest = slopei*(ytest-y_sommet0) + x_sommet0;
              //  Cerr << "xtest = " << xtest << " slopei = "<< slopei << finl;
              if(xtest <= rface_dis && xtest >= lface_dis)
                {
                  xint[nint] = xtest;
                  yint[nint] = delta_y;
                  int_face += 1;
                  ///   Cerr << "intersection is at the top-face" << finl;
                }
              else
                {
                  xtest = rface_dis;
                  ytest = slope*(xtest-x_sommet0) + y_sommet0;
                  //  Cerr<< "nint = "<< nint << finl;
                  xint[nint] = xtest;
                  yint[nint] = ytest;
                  int_face += 1;
                  //  Cerr << "yint = " << ytest << "xint = "<< xint[nint]<<" slope = "<< slope << finl;
                  //   Cerr << "intersection is at the right-face" << finl;
                }
              //      Cerr << "Elem #" << elem << " ymax = " << yr << " xmax = " << xr <<finl;
            }
          if(x_sommet1 <= lface_dis) // intersection is either at left or at top
            {
              //  Cerr << "x_sommet1 is less than lface_dis" << finl;
              //test for left-face
              xtest = lface_dis;
              ytest = slope*(xtest-x_sommet0) + y_sommet0;
              if(ytest <= tface_dis && ytest >= bface_dis)
                {
                  xint[nint] = xtest;
                  yint[nint] = ytest;
                  int_face += 1;
                  //   Cerr << "intersection is at the left-face" << finl;
                }
              else
                // test for top-face
                {
                  ytest = delta_y;
                  xtest = slopei*(ytest-y_sommet0) + x_sommet0;
                  xint[nint] = xtest;
                  yint[nint] = delta_y;
                  int_face += 1;
                  //   Cerr << "intersection is at the top-face" << finl;
                }
            }
          if(x_sommet1 < rface_dis && x_sommet1 > lface_dis)
            {
              //    Cerr << "x_sommet1 is less than rface_dis and greater than lface_dis" << finl;
              //   Cerr << "intersection is at the top-face" << finl;
              // test for top-face
              ytest = delta_y;
              xtest = slopei*(ytest-y_sommet0) + x_sommet0;
              if(x_sommet0 > x_sommet1)
                {
                  xtest = -slopei*(ytest-y_sommet0) + x_sommet0;
                }
              yint[nint] = ytest;
              xint[nint] = xtest;
              int_face += 1;
            }
        }

      if(inner_nodes == 0.)
        {
          for (int isom = 0; isom< Objet_U::dimension; isom++)
            {
              const int num_som = facettes(fa7, isom); // numero du sommet dans le tableau sommets
              const double bary_som = data.barycentre_[isom];
              for (int dir = 0; dir < Objet_U::dimension; dir++)
                {
                  coord_sommet[dir] = sommets(num_som,dir);
                  coord_barycentre_fraction[dir] += bary_som * sommets(num_som,dir) * surf;
                }
              if(isom == 0.)
                {
                  x_sommet0 = coord_sommet[0];
                  y_sommet0 = coord_sommet[1];
                  //   nx0 = norm[0];
                  //   ny0 = norm[1];
                }
              if(isom == 1.)
                {
                  x_sommet1 = coord_sommet[0];
                  //  y_sommet1 = coord_sommet[1];
                  //  nx1 = norm[0];
                  //  ny1 = norm[1];
                }
            }
          //     Cerr << "outer node --- x_sommet0 = " << x_sommet0 << "y_sommet0 = "
          //          << y_sommet0 << " normalx = " << nx0 << " normaly = " << ny0 << finl;

          //     Cerr << "outer node --- x_sommet1 = " << x_sommet1 << "y_sommet1 = "
          //         << y_sommet1 << " normalx = " << nx1 << " normaly = " << ny1 << finl;

          // test for top
          ytest = delta_y;
          xtest = slopei*(ytest-y_sommet0) + x_sommet0;
          //     Cerr << " xtest_top = " << xtest << finl;
          if(xtest <= rface_dis && xtest >= lface_dis)
            {
              xint[nint] = xtest;
              yint[nint] = ytest;
              int_face += 1;
              //       Cerr << "intersection is at the top-face" << finl;
              nint += 1;
            }
          else if(std::fabs(xtest-lface_dis) < eps || std::fabs(xtest-rface_dis) < eps)
            {
              xint[nint] = xtest;
              yint[nint] = ytest;
              int_face += 1;
              //       Cerr << "intersection is at the top-face" << finl;
              nint += 1;
            }
          // test for left
          xtest = lface_dis;
          ytest = slope*(xtest-x_sommet0) + y_sommet0;
          //     Cerr << " ytest_left = " << ytest << finl;
          if(ytest <= delta_y && ytest >= bface_dis)
            {
              xint[nint] = xtest;
              yint[nint] = ytest;
              int_face += 1;
              //        Cerr << "intersection is at the left-face" << finl;
              nint += 1;
            }
          else if(std::fabs(ytest-tface_dis) < eps || std::fabs(ytest-bface_dis) < eps)
            {
              xint[nint] = xtest;
              yint[nint] = ytest;
              int_face += 1;
              //       Cerr << "intersection is at the left-face" << finl;
              nint += 1;
            }
          // test for right
          xtest = rface_dis;
          ytest = slope*(xtest-x_sommet0) + y_sommet0;
          //       Cerr << " ytest_right = " << ytest << finl;
          if(ytest <= delta_y && ytest >= bface_dis)
            {
              xint[nint] = xtest;
              yint[nint] = ytest;
              int_face += 1;
              //          Cerr << "intersection is at the right-face" << finl;
              nint += 1;
            }
          else if(std::fabs(ytest-tface_dis) < eps || std::fabs(ytest-bface_dis) < eps)
            {
              xint[nint] = xtest;
              yint[nint] = ytest;
              int_face += 1;
              //        Cerr << "intersection is at the right-face" << finl;
              nint += 1;
            }
        }
      if (int_face > 0)
        nint += 1;

      for (int dir = 0; dir < Objet_U::dimension; dir++)
        bary_facettes_dans_elem[dir] += coord_barycentre_fraction[dir];
      index = data.index_facette_suivante_;
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

double Triple_Line_Model_FT_Disc::compute_Qint(const DoubleTab& in_out, const double theta_app_loc,
                                               const int num_face, double& Q_meso) const
{
  double radial_dis = 0.;
  const double xl = in_out(0,0);
  const double xr = in_out(1,0);
  radial_dis = (xr + xl)/2;
  double circum_dis = 1.;
  if (Objet_U::bidim_axi)
    {
      const double angle_bidim_axi = Maillage_FT_Disc::angle_bidim_axi();
      circum_dis = angle_bidim_axi*radial_dis;
    }
  Cerr << "circum_dis = " << circum_dis << finl;
  //              const Echange_impose_base& la_cl_typee = ref_cast(Echange_impose_base, la_cl_face.valeur());
  //              const double T_imp = la_cl_typee.T_ext(num_face-ndeb);
  //              Cerr << "T_imp= " << T_imp << finl;

// const Transport_Interfaces_FT_Disc& transport = ref_eq_interf_.valeur();
// const Maillage_FT_Disc& maillage = transport.maillage_interface();
// const double temps = maillage.temps();
  const double yl = in_out(0,1);
  const double yr = in_out(1,1);
  double Twall = 0.;
  //Twall = ref_eq_temp_.valeur().get_Twall(num_face);
  //Cerr << Twall << finl;
  double flux = 0.;
  double d=0;
  {
    const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ref_eq_temp_.valeur().domaine_dis().valeur());
    const IntTab& faces_elem = domaine_vf.face_voisins();
    // On of the neighbours doesnot exist so it has "-1". We get the other elem by:
    const int elem = faces_elem(num_face, 0) + faces_elem(num_face, 1) +1;

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
  }
  ref_eq_temp_.valeur().get_flux_and_Twall(num_face, d, flux, Twall);

// Cerr << "Flux/Twall : " << flux << " " << Twall << finl;
  //Process::exit();

  const double ytop = std::fmax(yr,yl);
  const double ybot = std::fmin(yr,yl);
  double ln_y = log(std::fmax(ytop,ym_)/std::fmax(ybot,ym_));
  assert(kl_cond_>0.);
//  Cerr << "ln_y = " << ln_y << " time_total = " << temps << " Theta_app_local = " << theta_app_loc
//       <<" delT = "<< Twall << " kl = "<< kl_cond_ << finl;
  Q_meso = kl_cond_*(Twall/theta_app_loc)*ln_y; // Twall here is (Wall temperature - saturation temperature)..unit of Q_meso is W/m

  double Q_int = Q_meso*circum_dis; // unit of Q_int is W
//  Cerr << "Q_meso = " << Q_meso << finl;
  return Q_int;
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

// Arguments :
//    elems_with_CL_contrib : The list of elements containing either the TCL itself, or the meso domaine.
//                            In all of them, the flux_evap has to be modified.
//    mpoint_from_CL : corresponding value to be assigned to each cell (as a mass flux m [unit?])
//    Q_from_CL : corresponding value to be assigned to each cell (as a thermal flux Q [unit?])
//
// Note that the same elem may appear twice (or more when contribution has to be displaced...) in the list,
// once (or more) for the micro contribution, once for the meso.
// It is not a problem, it will be considered by "+=" later
//
// num_faces is a locally temporary storages that in the end, is copied back into the class attribute boundary_faces_.
void Triple_Line_Model_FT_Disc::compute_TCL_fluxes_in_all_boundary_cells(ArrOfInt& elems_with_CL_contrib,
                                                                         ArrOfDouble& mpoint_from_CL, ArrOfDouble& Q_from_CL)
{
  ArrOfInt num_faces;
  num_faces.set_smart_resize(1);
  num_faces.resize_array(0);
  //  ArrOfInt elems_with_CL_contrib;
  // ArrOfDouble mpoint_from_CL;
  // ArrOfDouble Q_from_CL;
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
  const int dim3 = Objet_U::dimension == 3;
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
  //const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();

  // We store the elements (eulerian) crossed by a contact line.
  // As the specific source term is introduced here, a contribution of mdot should not be included afterwards.

  elems_with_CL_contrib.set_smart_resize(1);
  elems_with_CL_contrib.resize_array(0);

  mpoint_from_CL.set_smart_resize(1);
  mpoint_from_CL.resize_array(0);

  Q_from_CL.set_smart_resize(1);
  Q_from_CL.resize_array(0);

//  DoubleVect dI_dt;
//  ns.domaine_dis().valeur().domaine().creer_tableau_elements(dI_dt);
//  ns.calculer_dI_dt(dI_dt);

  // First loop on contact line cells:
  {
    int fa7;
    const int nb_facettes = maillage.nb_facettes();
    for (fa7 = 0; fa7 < nb_facettes; fa7++)
      {
        const int sommet0 = facettes(fa7, 0);
        const int sommet1 = facettes(fa7, 1);
        int sommet2 = -1;
        int tcl_s2 = 1; // 1 if s2 is on the contact line.
        if (dim3)
          {
            sommet2 = facettes(fa7, 2);
            tcl_s2 = maillage.sommet_ligne_contact(sommet2);
          }
        int tcl_s0 = maillage.sommet_ligne_contact(sommet0)*(1- maillage.sommet_virtuel(sommet0)); // We don't want to stop on virtual contact line
        int tcl_s1 = maillage.sommet_ligne_contact(sommet1)*(1- maillage.sommet_virtuel(sommet1));
        if (tcl_s0+tcl_s1+tcl_s2 >= 2)
          {
            // We have a contact line in this cell...
            // In most cases, contact line will come to a sum of 2,
            // but in corners, or a 1-cell channel, it can be 3.
            //
            int num_bc_face=-1;
            int num_som = 0;
            if (tcl_s0)
              {
                num_som = sommet0;
                num_bc_face = maillage.sommet_face_bord(sommet0);
              }
            if (tcl_s1)
              {
                num_som = sommet1;
                num_bc_face = maillage.sommet_face_bord(sommet1);
              }
            // We need to check on what boundary we are:
            Journal() << "iii " << num_bc_face << finl;
            const int num_bord = domaine_vf.face_numero_bord(num_bc_face);
            Journal() << "jjj " << num_bord << finl;
            assert(num_bord>=0); // Otherwise we're inside and there's a big problem
            //const Frontiere_dis_base& la_cl =  domaine_vf.frontiere_dis(num_bord);
            const Domaine_Cl_dis_base& zcldis = ns.domaine_Cl_dis().valeur();
            const Cond_lim& la_cl = zcldis.les_conditions_limites(num_bord);
            const Nom& bc_name = la_cl.frontiere_dis().le_nom();

            Cerr << que_suis_je() << "::derivee_en_temps_inco() computing contrib at TCL. CL "
                 << la_cl.valeur();// << finl;
            // For each BC, we check its type to see if it's a wall:
            if ( sub_type(Dirichlet_paroi_fixe,la_cl.valeur())
                 || sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) )
              {
                Cerr << "  -> for this BC [" <<
                     bc_name <<"], we compute a specific phase-change rate at TCL." << finl;
              }
            else
              {
                Cerr << "  -> nothing for this BC [" << bc_name <<"] at TCL." << finl;
                // We're on the symmetry axis, or something else but not at the wall...
                continue;
              }
// 2D axisymetric case :
            double contact_line_length = 1.0;
            if (Objet_U::bidim_axi)
              {
                const double angle_bidim_axi = Maillage_FT_Disc::angle_bidim_axi();
                x_cl_ = sommets(num_som,0);
                contact_line_length = angle_bidim_axi*x_cl_;
              }
            else
              {
                // Not-bidim-axi :
                if (Objet_U::dimension == 2)
                  {
                    contact_line_length = 1.; // Assume 1m in the 3rd dimension
                    x_cl_ = 0.;
                    Process::exit();
                  }
                else
                  {
                    Cerr << " 3D case, not axi-symetric, to be coded later. " << finl;
                    Process::exit();
                  }
              }
            /* double surface_tot = 0.;
            int index = index_facette_elem[fa7];
            while (index >= 0)
              {
                const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                const double surface_facette = surface_facettes[fa7];
                const double surf = data.fraction_surface_intersection_ * surface_facette;
                surface_tot +=surf;
                index = data.index_element_suivant_;
              } */
            int index = index_facette_elem[fa7];
            while (index >= 0)
              {
                const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                double fraction = data.fraction_surface_intersection_; // this gives a weight... It is the fraction of the facette in contact with the wall (ie the CL facette)
                int elem = data.numero_element_;
                int num_face=-1;
                {
                  // If the contact angle is high enough, it may happen that the first front element (that
                  // contains a node of the contact line) ends above the row of cells directly in contact with
                  // the wall. In that case, the contribution of the corresponding fraction of facette
                  // Should be transfered to the previous
                  // const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ref_eq_temp_.valeur().domaine_dis().valeur());
                  const IntTab& elem_faces = domaine_vf.elem_faces();
                  const IntTab& faces_elem = domaine_vf.face_voisins();
                  const int nb_faces_voisins = elem_faces.dimension(1);
                  // Struggle to get the boundary face
                  int i;
                  for (i=0; i<nb_faces_voisins; i++)
                    {
                      num_face = elem_faces(elem,i);
                      // If it's a boundary face, one of the neighbours doesnot exist so it has "-1".
                      // We detect a boundary that way:
                      const int elemb = faces_elem(num_face, 0) + faces_elem(num_face, 1) +1;
                      if (elem == elemb)
                        {
                          break;
                        }
                    }
                  if (i==nb_faces_voisins)
                    {
                      Cerr << "No boundary face found in this element "<< elem
                           <<  " because the first front-element is steep and long enough to go all the way up "
                           << "to the second layer of cells... " << finl;
                      Cerr << "We have to adapt and transfer/report the contribution "
                           << "to the previous near wall-cell. " << finl;

                      const int nb_elems_so_far = elems_with_CL_contrib.size_array();
                      if (nb_elems_so_far==0)
                        {
                          Cerr << "Let's cry... what should we do? how can I report it to the next element?" << finl;
                          Process::exit();
                        }
                      else
                        {
                          const int previous_elem = elems_with_CL_contrib[nb_elems_so_far-1];
                          Cerr << "The TCL-micro contribution is transferred from " << elem << " to " << previous_elem << finl;
                          // On reporte totalement la contribution sur cet autre elem:
                          elem = previous_elem;
                          num_face = num_faces[nb_elems_so_far-1];
                          fraction = 1.;
                        }
                    }
                }
                const double temps = ns.schema_temps().temps_courant();
                //   Cerr << "temps= "<< temps << " dI_dt[TCL]= " << dI_dt(elem) << finl;
                //const double surface_facette = surface_facettes[fa7];
                //const double surf = data.fraction_surface_intersection_ * surface_facette;
                //    surface_tot +=surf;
                {
                  // const double v = volumes[elem];
                  //    Cerr << "Elem " << elem << " xp= " << domaine_vf.xp(elem,0) << " vol= " << v << finl;
                  double Qtot = get_Qtcl();
//                     if (dI_dt(elem) > 0.)
//                       {
//                         Qtot = 0.;
//                         Cerr << "Qtot =" << Qtot <<finl;
//                         Cerr << "  -> advancing contact line : disable TCL contribution to phase change" << finl;
//                       }

//                      else
//                        {
//                          Cerr << "  -> advancing contact line : disable TCL contribution to phase change" << finl;
//                          //    continue;
//                          const double Qtot = 0.;
                  //   Cerr << "Qtot=" << Qtot << " time_Qtot= " << temps << finl;
//                        }

// NEW COMMENT
                  if (dim3)
                    {
                      // We look for the two nodes on the contact line:
                      int s0,s1;
                      if (!tcl_s0)
                        {
                          s0 = sommet1;
                          s1 = sommet2;
                        }
                      else if (!tcl_s1)
                        {
                          s0 = sommet0;
                          s1 = sommet2;
                        }
                      else
                        {
                          // 2 or 3 points are contact line. Let's take the first two of which we're sure.
                          s0 = sommet0;
                          s1= sommet1;
                        }
                      double s0s1[3], xyz_f0[3], xyz_f1[3], f0f1[3];
                      s0s1[0] = sommets(s1,0) - sommets(s0,0);
                      s0s1[1] = sommets(s1,1) - sommets(s0,1);
                      s0s1[2] = sommets(s1,2) - sommets(s0,2);

                      int num_bc_face0 = maillage.sommet_face_bord(s0);
                      int num_bc_face1 = maillage.sommet_face_bord(s1);

                      if (num_bc_face0 == num_bc_face1)
                        {
                          // This Front element (facette) is in contact with only one bc_face:
                          contact_line_length = sqrt(s0s1[0]*s0s1[0] + s0s1[1]*s0s1[1] + s0s1[2]*s0s1[2]);
                        }
                      else
                        {
// The contact line is on several elements. We assume it's two, but it can be more!
                          // How can we know?
                          // We need to compute The fraction of TCL within each cell
                          xyz_f0[0] =  domaine_vf.xv(num_bc_face0,0);
                          xyz_f0[1] =  domaine_vf.xv(num_bc_face0,1);
                          xyz_f0[2] =  domaine_vf.xv(num_bc_face0,2);
                          xyz_f1[0] =  domaine_vf.xv(num_bc_face1,0);
                          xyz_f1[1] =  domaine_vf.xv(num_bc_face1,1);
                          xyz_f1[2] =  domaine_vf.xv(num_bc_face1,2);
                          f0f1[0] = xyz_f1[0] -  xyz_f0[0];
                          Cerr << "f0f1 " << f0f1[0] << finl;
                          // ...
                          // Need for geometry formulae to be exact...
                          // Let's do a simplification : assume we put half of the contribution
                          // to each of the two elem containing the FT contact line vertex:
                          contact_line_length = 0.5*sqrt(s0s1[0]*s0s1[0] + s0s1[1]*s0s1[1] + s0s1[2]*s0s1[2]);

                        }
                      //contact_line_length = pow(v, 1./3.); // if a cube... we should do better later.
                    }

                  // we multiply by the contact_line_length to represent :  int_v delta^i delta^wall dv homogeneous to [m]
                  //  double value = 0.;
                  //   const double value = jump_inv_rho*contact_line_length*Qtot/Lvap;
                  //  Cerr << "[TCL-micro] Local source in elem=["  << elem << "] with value= " << value << "[m3.s-1]" << finl;
                  //   Cerr << "contact_line_length= " << contact_line_length << " "  << v << " "  << Qtot << " "  << Lvap << " " << jump_inv_rho << finl;

//                      Cerr << "Secmem= "<< secmem(elem)*schema_temps().pas_de_temps() ;
//                      Cerr << " Secmem2 before= "<< secmem2(elem) ;
                  // @Vinod : See meso-scale model and comments. Should it be "+=" to be added on top of mpoint or simply setting "="?
                  //secmem2(elem) += fraction*value;
                  //   Cerr << "surface for tcl = " << surface_tot << " Volume of tcl cell = " << v << finl;

                  // BEWARE : If the front-element containing the TCL is overlapping two cells,
                  // then, the Qtot contribution will be counted twice. So, it is essential to use fraction
                  // to distribute the micro-contribution on 2 cells.
                  //
                  if (first_call_for_this_mesh)
                    {
                      // We update values and integrals only at the first call
                      instant_mmicro_evap_ += fraction*Qtot*contact_line_length/Lvap*dt;
                      instant_vmicro_evap_ = instant_mmicro_evap_*jump_inv_rho; // VP - I think it should be '=' not '+='
                    }
                  // no_elem += 1;
                  //      }
                  // instant_vmicro_evap_ += value*dt;
                  elems_with_CL_contrib.append_array(elem);
                  // mpoint_from_CL.append_array(fraction*Qtot*surface_tot/v/Lvap);  // GB 2019.11.11 Correction as Lvap was missing.
                  mpoint_from_CL.append_array(fraction*Qtot/(sm_*Lvap));  // GB 2021.05.27 Correction surface_tot/v replaced by 1/sm.
                  instant_Qmicro += fraction*Qtot;
                  Q_from_CL.append_array(fraction*Qtot*contact_line_length); // I'm not sure, but I believe "*contact_line_length" should be in it, as it is now.
                  num_faces.append_array(num_face);
                  Cerr << "temps = " << temps << " Instantaneous tcl-evaporation = "<< instant_vmicro_evap_ << finl;
                }
                index = data.index_element_suivant_;
              }
          }
      }
  }
  instant_Qmicro = mp_sum(instant_Qmicro);

  // Second loop to find meso cells:
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, ns.domaine_dis().valeur());
  const Domaine_Cl_dis_base& zcldis = ns.domaine_Cl_dis().valeur();

  //const Domaine_dis& aa= zcldis.domaine_dis();
  const Domaine_dis_base& domaine_dis = ns.domaine_dis().valeur();
  const IntTab& face_voisins = domaine_dis.face_voisins();
  //const IntTab& face_voisins = ns.le_dom_dis.valeur().valeur().face_voisins();
  //       const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, zcldis);
  // Boucle sur les bords pour traiter les conditions aux limites
  //            double m_evap = 0.;
  //            double instant_vmeso_evap_ = 0.;
  //            double Q_int = 0.;
  //    const DoubleVect& volumes = domaine_vf.volumes();
  int ndeb, nfin, num_face;
  for (int n_bord=0; n_bord<zvdf.nb_front_Cl(); n_bord++)
    {
      // Cerr << "I entered the loop" << finl;
      // break;
      //  const double v = volumes[elem];
      const Cond_lim& la_cl = zcldis.les_conditions_limites(n_bord);
      //         const Cond_lim& la_cl_face = zclvdf.la_cl_de_la_face(num_face);
      const Nom& bc_name = la_cl.frontiere_dis().le_nom();
      Cerr << que_suis_je() << "::derivee_en_temps_inco() computing near TCL contrib "
           << "at CL " << la_cl.valeur();// << finl;
      // For each BC, we check its type to see if it's a wall:
      if ( sub_type(Dirichlet_paroi_fixe,la_cl.valeur())
           || sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) )
        {
          Cerr << "  -> for this BC [" <<
               bc_name <<"], we compute a specific phase-change rate near TCL." << finl;
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          // Vecteur3& normale = 0.; // added vp

          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              // Tricky way to find the element (one of the neighbour doesn't exists and contains -1...)
              const int elem = -face_voisins(num_face, 0) * face_voisins(num_face, 1);
              //   const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(elem);
              //   const double fraction = data.fraction_surface_intersection_; // this gives a weight...
              const double vol = volumes[elem];
              int index=intersections.index_elem()[elem];
              if (index < 0)
                continue; // This element is not crossed by the interface. Go to the next face.

              //  Cerr << "Elem #" << elem << " is adjacent to the wall " << bc_name
              //      << " and also contains an interface. So it belongs to meso domaine." << finl;

              double surface_tot = 0.;
              // const double v = volumes[elem];
              DoubleTab in_out; // The coords of left and right points where the interface cross the cell boundary
              FTd_vecteur3 norm_elem = {0., 0., 0.};
              get_in_out_coords(zvdf, elem, dt, in_out, norm_elem, surface_tot);
              double factor = surface_tot/pow(vol, 2./3.);
              //   Cerr << "factor= " << factor << finl;
              if (factor < 1.e-6)
                {
                  Cerr << "We are skipping cell " << elem << "because the interface in it has negligible surface!!" << finl;
                  continue; // skip cells crossed by almost no interface.
                }
              // /*
              const double xl = in_out(0,0);
              const double yl = in_out(0,1);
              const double xr = in_out(1,0);
              const double yr = in_out(1,1);
              Cerr << "xl,yl, xr,yr = " << xl << " " <<  yl << " " <<  xr << " " <<  yr << finl;
              // */
              double theta_app_loc = compute_local_cos_theta(parcours, num_face, norm_elem);
              double Q_meso = 0.;
              double Q_int = compute_Qint(in_out, theta_app_loc, num_face, Q_meso);
              // const double v = volumes[elem];
              //  const double value = jump_inv_rho*Q_int/Lvap;
              //     Cerr << "[TCL-meso] Local source in elem=["  << elem << "] with value= " << value << "[m2.s-1]" << finl;

              // @Vinod : to analyze closely... At first, I thought it should be "+=", but it seems we want to replace the mpoint
              // contribution totally. Is it enough to change from "+=" to simply setting secmem2 (=) in that cell?
              // That way, mpoint remains has it is from the extrapolation and is ok for the neighbors but is not used in that cell...
              // secmem2(elem) += value;
              //     Cerr << "surface_total for meso cell = " << surface_tot << " Volume of meso cell = " << v << finl;
              //     Cerr << "Qint= " << Q_int << finl;

              if (first_call_for_this_mesh)
                {
                  // We update values and integrals only at the first call
                  //  instant_mmeso_evap_ += (Q_int*surface_tot/Lvap)*dt;          // total mass of liquid evaporated
                  instant_mmeso_evap_ += (Q_int/Lvap)*dt;          // total mass of liquid evaporated
                  instant_vmeso_evap_ = (instant_mmeso_evap_*jump_inv_rho);   // total volume of liquid evaporated
                  // VP - '+=' or '='?
                  integrated_vmeso_evap_ = instant_vmeso_evap_; // ?*dt;
                  //  Cerr << "Qint= " << Q_int << "dt = " << dt << finl;
                  //    Cerr << "time = " << integration_time_ << " Qmeso= " << Q_meso << " Qint " << Q_int << finl;
                }
              Cerr.precision(16);
              Cerr << "time = " << integration_time_ << " Qmeso= " << Q_meso << " Qint " << Q_int << finl;
              Cerr.precision(7);
              //   Cerr << "time = " << integration_time_ << " instantaneous mass-meso-evaporation = " << instant_mmeso_evap_ << " instantaneous meso-evaporation = " << instant_vmeso_evap_ << finl;
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
                  // Cerr << Q_meso << finl;
                  // Process::exit();
                  mpoint_from_CL.append_array(0.);
                }
              else
                mpoint_from_CL.append_array(Q_meso/(s_meso*Lvap));  // GB 2021.05.27 Correction surface_tot/v replaced by 1/s_meso.
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
                    Cerr << "[Assessment-S_meso] t= "<< integration_time_ << " elem= " << elem << " LR= " ;
                    Cerr << s_meso <<  " Ameso/2pir= " << s_meso_assessed << " err(%)= " << 50.*std::abs(s_meso_assessed-s_meso)/(s_meso_assessed+s_meso) << finl;
                    Cerr.precision(7);
                  }
              }
              // mpoint_from_CL.append_array((Q_int*surface_tot/v)/Lvap);
              //    Q_from_CL.append_array(Q_meso*surface_tot); //VP: corrections (1. Q_int replaced by Q_meso and multiplied by surface_tot)...unit--W
              Q_from_CL.append_array(Q_int);
              num_faces.append_array(num_face);
            }

          //      instant_vmeso_evap_ += value*dt
          Cerr << "time = " << integration_time_ << " instantaneous mass-meso-evaporation = " << instant_mmeso_evap_ << " instantaneous meso-evaporation = " << instant_vmeso_evap_ << finl;
        }
      else
        {
          Cerr << "Nothing specific for this type of BC " << bc_name << " near TCL. " << finl;
        }
//                  instant_mmeso_evap_ += Q_int/Lvap;          // total mass of liquid evaporated
//                  instant_vmeso_evap_ += instant_mmeso_evap_/rho_phase_1;   // total volume of liquid evaporated
    }
  instant_Qmeso = mp_sum(instant_Qmeso);

  const double verbosity=1;
  if (verbosity)
    {
      const int nb_contact_line_contribution = elems_with_CL_contrib.size_array();
      const double t = ns.schema_temps().temps_courant();
      Cerr << " ********* VALUES OF CONTACT LINE CONTRIBUTION *************** " << finl;
      Cerr << " time " << t << finl;
      Cerr << " nb_contact_line_contribution " << nb_contact_line_contribution << finl;
      Cerr << "# elem \t face \t Q \t\t mp" << finl;
      for (int idx = 0; idx < nb_contact_line_contribution; idx++)
        {
          Cerr << elems_with_CL_contrib[idx] << " \t " <<  num_faces[idx] << " \t "
               << Q_from_CL[idx] << " \t " << mpoint_from_CL[idx] << " \t " << finl;
//          Cerr << elems_with_CL_contrib[idx] << " \t " <<  Q_from_CL[idx] << " \t " << mpoint_from_CL[idx] << " \t " << finl;
        }
      Cerr << " ************************************************************* " << finl;
      Cerr << "[instant-Q] time micro meso " << t << " " <<
           instant_Qmicro << " " << instant_Qmeso << finl;
      Cerr << " ************************************************************* " << finl;
    }

//  if (first_call_for_this_mesh)
//    {

//                      integrated_vmicro_evap_ += instant_vmicro_evap_; //?
//   }

  // double micro_mevap = instant_mmicro_evap_ + instant_mmeso_evap_;
  vevap_int_ += instant_vmicro_evap_ + instant_vmeso_evap_;
//    vevap_int += vevap;
  // You want the total value on all processors:
  // micro_mevap = Process::mp_sum(micro_mevap);

  // Temporary code : Later, we will either (i) fill elems_ directly without using a temporary elems_with_CL_contrib
  // or (ii) make a list without repetition (summing up values for micro and meso together)
  // Store :

  elems_.resize_array(0); // to empty the list
  boundary_faces_.resize_array(0);
  mp_.resize_array(0);
  Q_.resize_array(0);
  for (int idx = 0; idx < elems_with_CL_contrib.size_array(); idx++)
    {
      const int elem = elems_with_CL_contrib[idx];
      elems_.append_array(elem);
      boundary_faces_.append_array(num_faces[idx]);
      mp_.append_array(mpoint_from_CL[idx]);
      Q_.append_array(Q_from_CL[idx]);
    }
  // vevap = Process::mp_sum(vevap);
//  Cerr << "vevap(micro+meso) = " << vevap << " time = " << integration_time_ << finl;
//  vevap_int_ = Process::mp_sum(vevap_int_);
  Cerr << "vevap_int_(micro+meso) = " << vevap_int_ << " time = " << integration_time_ << finl;
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
      //secmem2(elem) =0.;
      mpoint(elem) = 0.;
    }

  // Then, we can sum our micro and meso contributions :
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      //const double Q = Q_from_CL[idx];
      //const double value = jump_inv_rho/Lvap*Q;
      //secmem2(elem) += value;
      mpoint(elem) += mp_[idx];
      Cerr << "contrib added to mpoint for cell " << elem << " with value= " << mp_[idx] << finl;
    }
}

// correct secmem2 with values from TCL model.
void Triple_Line_Model_FT_Disc::corriger_secmem(const double coef, DoubleTab& secmem2) const
{
//  const Navier_Stokes_FT_Disc& ns = ref_ns_.valeur();
//  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ns.domaine_dis().valeur());
//  const DoubleVect& volumes = domaine_vf.volumes();
  const int nb_contact_line_contribution = elems_.size_array();
// const double max_val_before = max_array(secmem2);
//  const double min_val_before = min_array(secmem2);
  // I believe it is required to put zero in all cells first to remove any possible macroscopic contribution.
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
//      const double v = volumes[elem];
      const double value = coef*Q_val;
      secmem2(elem) += value; // '+=' or '='? - VP
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
      //const double Twall =ref_eq_temp_.valeur().get_Twall_at_elem(elem);
      const int num_face = boundary_faces_[idx];
      const double Twall =ref_eq_temp_.valeur().get_Twall_at_face(num_face);
      const double DeltaT = TSAT_CONSTANTE-Twall;
      const double Tbefore =  temperature(elem);
      temperature(elem) =Twall + 0.5*DeltaT;
      Cerr << "[Temperature correction] elem= " << elem << " before/after: " << Tbefore << " / " << temperature(elem) << finl;
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
                                               -1. /* distance_wall_interface  is not needed to get Twall */,
                                               flux, Twall);
      temperature(elem) =Twall;
    }
  // MEMO : remove the previous loop if num_faces is computed here.
  // Then, we can sum our micro and meso contributions :
  for (int idx = 0; idx < nb_contact_line_contribution; idx++)
    {
      const int elem = elems_[idx];
      const int num_face = boundary_faces_[idx];
      double d=0;
      {
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
      }
      double flux=0., Twall=0.;
      ref_eq_temp_.valeur().get_flux_and_Twall(num_face, d, flux, Twall);
      // It is the total flux, so It should simply be :
      temperature(elem) = Twall - flux*d/kl_cond_; // Can be done several times, no problem.
    }
}
