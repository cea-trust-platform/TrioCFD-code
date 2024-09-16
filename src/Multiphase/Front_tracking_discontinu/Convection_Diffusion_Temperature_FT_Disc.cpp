/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Convection_Diffusion_Temperature_FT_Disc.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/25
//
//////////////////////////////////////////////////////////////////////////////
#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Domaine_VF.h>
#include <Discretisation_base.h>
#include <Probleme_FT_Disc_gen.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Domaine.h>
#include <Sous_Domaine.h>
#include <Param.h>
#include <Champ_Uniforme.h>
#include <Echange_impose_base.h>
#include <Echange_contact_VDF_FT_Disc.h>
#include <Domaine_Cl_VDF.h>
#include <Neumann_paroi.h>
#include <Neumann_paroi_adiabatique.h>
#include <SChaine.h>
#include <Entree.h>
#include <EChaine.h>
#include <Interprete_bloc.h>
#include <Domaine_VDF.h>
#include <stat_counters.h>
#include <TRUST_Ref.h>
#include <Parametre_implicite.h>

static const double TSAT_CONSTANTE = 0.;

Implemente_instanciable_sans_constructeur(Convection_Diffusion_Temperature_FT_Disc,"Convection_Diffusion_Temperature_FT_Disc",Convection_Diffusion_Temperature);

Convection_Diffusion_Temperature_FT_Disc::Convection_Diffusion_Temperature_FT_Disc()
{
  phase_ = -1;
  correction_courbure_ordre_=0; // Par defaut, pas de prise en compte de la courbure pour corriger le champ etendu delta_vitesse
  stencil_width_ = 8; //GB : Valeur par defaut de stencil_width_
  temp_moy_ini_ = 0.; //GB : Valeur par defaut de la temperature moyenne initiale
  nom_sous_domaine_ = "unknown_sous_domaine"; //GB : Valeur par defaut de la sous-domaine de moyenne
  maintien_temperature_ = false;
  is_prescribed_mpoint_ = false;
  prescribed_mpoint_ = -1.e30;
  mixed_elems_.resize_array(0);
  lost_fluxes_.resize_array(0);
  derivee_energy_.resize_array(0);
  mixed_elems_diffu_.resize_array(0);
  lost_fluxes_diffu_.resize_array(0);
  mixed_elems_conv_.resize_array(0);
  lost_fluxes_conv_.resize_array(0);
  divergence_free_velocity_extension_=0; // Default set to historical behavior : velocity extension is NOT divergence-free
}

Sortie& Convection_Diffusion_Temperature_FT_Disc::printOn(Sortie& os) const
{
  return Convection_Diffusion_Temperature::printOn(os);
}
/*! @brief cf Convection_Diffusion_std::readOn(Entree&).
 *
 */
Entree& Convection_Diffusion_Temperature_FT_Disc::readOn(Entree& is)
{
  // Ne pas faire assert(fluide non nul)
  Convection_Diffusion_std::readOn(is);
  solveur_masse->set_name_of_coefficient_temporel("rho_cp_comme_T");
  Nom num=inconnue()->le_nom(); // On prevoir le cas d'equation de scalaires passifs
  num.suffix("temperature_thermique");
  Nom nom="Convection_chaleur";
  nom+=num;
  terme_convectif.set_fichier(nom);
  terme_convectif.set_description((Nom)"Convective heat transfer rate=Integral(-rho*cp*T*u*ndS) [W] if SI units used");
  nom="Diffusion_chaleur";
  nom+=num;
  terme_diffusif.set_fichier(nom);
  terme_diffusif.set_description((Nom)"Conduction heat transfer rate=Integral(lambda*grad(T)*ndS) [W] if SI units used");

  // If keyword correction_mpoint_diff_conv_energy_ has not been read, the table should be properly set:
  if (correction_mpoint_diff_conv_energy_.size_array() != 3)
    {
      Cerr << "We account for no energy correction!!" << finl;
      correction_mpoint_diff_conv_energy_.resize_array(3);
      correction_mpoint_diff_conv_energy_[0] = 0;
      correction_mpoint_diff_conv_energy_[1] = 0;
      correction_mpoint_diff_conv_energy_[2] = 0;
    }
  return is;
}

void Convection_Diffusion_Temperature_FT_Disc::set_param(Param& param)
{
  Convection_Diffusion_Temperature::set_param(param);
  param.ajouter("phase",&phase_);
  param.ajouter_condition("(value_of_phase_eq_0)_or_(value_of_phase_eq_1)","phase must be set to 0 or 1");
  param.ajouter("stencil_width",&stencil_width_);
  param.ajouter("correction_courbure_ordre", &correction_courbure_ordre_);
  param.ajouter_non_std("equation_interface",(this),Param::REQUIRED);
  param.ajouter_non_std("maintien_temperature",(this));
  param.ajouter_non_std("equation_navier_stokes",(this),Param::REQUIRED);
  param.ajouter_non_std("prescribed_mpoint", (this));
  param.ajouter("correction_mpoint_diff_conv_energy", &correction_mpoint_diff_conv_energy_);
  param.ajouter_flag("divergence_free_velocity_extension", &divergence_free_velocity_extension_, Param::OPTIONAL);
  param.ajouter_non_std("solveur_pression_fictive",(this),Param::OPTIONAL);
  param.ajouter("bc_opening_pressure",&name_bc_opening_pressure_,Param::OPTIONAL);
}

int Convection_Diffusion_Temperature_FT_Disc::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="equation_interface")
    {
      const Probleme_FT_Disc_gen& pb = ref_cast(Probleme_FT_Disc_gen,probleme());
      Motcle nom_eq;
      is >> nom_eq;
      if (Process::je_suis_maitre())
        Cerr << " Interface equation for the temperature convection diffusion equation :"
             << nom_eq << finl;
      ref_eq_interface_ = pb.equation_interfaces(nom_eq);
      return 1;
    }
  else if (mot=="solveur_pression_fictive")
    {
      divergence_free_velocity_extension_=1;
      Cerr << "Reading and typing of fictitious pressure solver (for velocity extension) : " << finl;
      is >> solveur_pression_;
      Cerr<<"Fake Pressure solver type : "<<solveur_pression_->que_suis_je()<< finl;
      solveur_pression_.nommer("solveur_pression_fictive");
      return 1;
    }
  else if (mot=="diffusion")
    {
      if (je_suis_maitre())
        Cerr << " "<< que_suis_je() <<"::lire diffusivite" << finl;
      // Need phase to know which diffusivity to use:
      if (phase_ < 0)
        {
          barrier();
          Cerr << " Error: phase must be specified before diffusion" << finl;
          exit();
        }
      Convection_Diffusion_std::lire_motcle_non_standard(mot,is);
      return 1;
    }
  else if (mot=="convection")
    {
      if (je_suis_maitre())
        Cerr <<  " "<< que_suis_je() <<"::lire convection" << finl;
      if (!ref_eq_ns_.non_nul())
        {
          barrier();
          Cerr << " Error: equation_navier_stokes must be specified before convection" << finl;
          exit();
        }
      Convection_Diffusion_std::lire_motcle_non_standard(mot,is);
      return 1;
    }
  else if (mot=="maintien_temperature")
    {
      if (je_suis_maitre())
        Cerr << " Equation_Concentration_FT::lire maintien_temperature" << finl;
      maintien_temperature_ = true;
      is >> nom_sous_domaine_;
      is >> temp_moy_ini_;
      return 1;
    }
  else if (mot=="prescribed_mpoint")
    {
      if (je_suis_maitre())
        Cerr << que_suis_je() <<"::lire prescribed_mpoint" << finl;
      is_prescribed_mpoint_ = true; // set flag to true.
      is >> prescribed_mpoint_;
      return 1;
    }
  else if (mot=="equation_navier_stokes")
    {
      const Probleme_FT_Disc_gen& pb = ref_cast(Probleme_FT_Disc_gen,probleme());
      Motcle nom_eq;
      is >> nom_eq;
      if (Process::je_suis_maitre())
        Cerr << " Navier_stokes equation for the temperature convection diffusion equation :"
             << nom_eq << finl;
      const Navier_Stokes_std& ns = pb.equation_hydraulique(nom_eq);
      ref_eq_ns_ = ns;
      return 1;
    }
  else
    return Convection_Diffusion_Temperature::lire_motcle_non_standard(mot,is);
}

const Champ_base& Convection_Diffusion_Temperature_FT_Disc::vitesse_pour_transport() const
{
  // return ref_eq_ns_->vitesse();
  // GB 2019.03.07 (code lost between 2022.02.16 and 2023.06.20!!)
  // We will go through all the effort of creating a continous velocity from the unknown of NS!
  // We might as well use it! ;-)
  return vitesse_convection_;
}

void Convection_Diffusion_Temperature_FT_Disc::preparer_pas_de_temps(void)
{
}

int Convection_Diffusion_Temperature_FT_Disc::get_phase() const
{
  return phase_;
}

void Convection_Diffusion_Temperature_FT_Disc::corriger_pas_de_temps(double dt)
{
}

static void extrapolate(const Domaine_VF&    domaine_vf,
                        const DoubleTab& interfacial_area,
                        const int stencil_width,
                        const DoubleTab& distance,
                        DoubleTab&        field)
{
  const double invalid_test = -1.e30;
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nb_faces_elem = elem_faces.dimension(1);
  const int nb_elem       = elem_faces.dimension(0);
  DoubleTab field_old;
  // n_iterations = stencil_width is the minimum to get a propagation of information from the interface to the border
  // of the extrapolation. But doing more will lead to smoother values... And it probably costs close to nothing
  const double n_iterations = 5*stencil_width;
  for (int iteration = 0; iteration < n_iterations; iteration++)
    {
      // Copy the old field value as we do not want to use the current iteration values.
      field_old = field;
      // La valeur sur un element est la moyenne des valeurs sur les elements voisins
      for (int i_elem = 0; i_elem < nb_elem; i_elem++)
        {
          // Do not touch field in interfacial cells.
          // Iterate on other values.
          const double d = distance[i_elem];
          if (( d > invalid_test) && (interfacial_area[i_elem]<DMINFLOAT))
            {
              double somme = 0.;
              double coeff = 0.;
              for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                {
                  const int face = elem_faces(i_elem, i_face);
                  const int voisin = faces_elem(face, 0) + faces_elem(face, 1) - i_elem;
                  if (voisin >= 0)
                    {
                      // Not a boundary...
                      double mp = field_old[voisin];
                      const double dvois = distance[voisin];
                      if (dvois > invalid_test)
                        {
                          // Give more weight in the smoothing to values closer to the interface:
                          if (fabs(dvois)<DMINFLOAT)
                            {
                              Cerr << "distance is very much at zero whereas interfacial_area is zero too... Pathological case to be looked into closely. " << finl;
                              Cerr << "Is it from a Break-up or coalescence? " << finl;
                              Cerr << "see Convection_Diffusion_Temperature_FT_Disc and static void extrapolate" << finl;
                              Cerr << "Contact TRUST support." << finl;
                              Process::exit();
                            }
                          const double inv_d2 = 1./(dvois*dvois);
                          somme += mp*inv_d2;
                          coeff+=inv_d2;
                        }
                    }
                }
              if (coeff > 0.)
                field[i_elem] = somme / coeff;
            }
        }
      field.echange_espace_virtuel();
    }
}

static void extrapoler_dans_phase(DoubleTab&        gradient,
                                  const DoubleTab& indicatrice,
                                  const Domaine_VF&    domaine_vf,
                                  const double invalid_test,
                                  const double indic_phase, const int nb_iter)
{
  DoubleTab gradient_old;
  for (int iteration = 0; iteration < nb_iter; iteration++)
    {
      // Copie de la valeur du gradient: on ne veut pas utiliser les valeurs
      // calculees lors de l'iteration courante
      gradient_old = gradient;
      // La valeur sur un element est la moyenne des valeurs sur les elements voisins
      for (int i_elem = 0; i_elem < domaine_vf.elem_faces().dimension(0); i_elem++)
        {
          if (indicatrice[i_elem] != indic_phase)
            {
              // Ne pas toucher au gradient de la phase "phase".
              // Iterer sur les autres valeurs.
              double somme = 0.;
              int coeff = 0;
              for (int i_face = 0; i_face < domaine_vf.elem_faces().dimension(1); i_face++)
                {
                  const int face = domaine_vf.elem_faces(i_elem, i_face);
                  const int voisin = domaine_vf.face_voisins(face, 0) + domaine_vf.face_voisins(face, 1) - i_elem;
                  if (voisin >= 0)
                    {
                      // Not a boundary...
                      const double grad = gradient_old[voisin];
                      if (grad > invalid_test)
                        {
                          somme += grad;
                          coeff++;
                        }
                    }
                }
              if (coeff > 0)
                gradient[i_elem] = somme / coeff;
            }
        }
      gradient.echange_espace_virtuel();
    }
}
// A partir des valeurs du "champ" dans la phase "phase", calcule
// les valeurs du "champ" dans un voisinage d'epaisseur "stencil_width"
// en extrapolant lineairement et en supposant que la valeur a l'interface
// vaut "interfacial_value".
static void extrapoler_champ_elem(const Domaine_VF&    domaine_vf,
                                  const DoubleTab& indicatrice,
                                  const DoubleTab& distance_interface,
                                  const DoubleTab& normale_interface,
                                  const DoubleTab& champ_div_n,
                                  const int   phase,
                                  const int   stencil_width,
                                  const double   interfacial_value,
                                  DoubleTab&        champ,
                                  DoubleTab&        gradient,
                                  const double temps)
{
  const IntTab& elem_faces = domaine_vf.elem_faces();
  //const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nb_faces_elem = elem_faces.dimension(1);
  const int nb_elem       = elem_faces.dimension(0);

  const DoubleTab& centre_gravite_elem = domaine_vf.xp();
  const DoubleTab& centre_gravite_face = domaine_vf.xv();
  const double invalid_test = -1.e30;
  const double invalid_value = -2.e30;
  // Calcul de la composante normale du gradient du champ:
  //  gradient normal = (champ - valeur_interface) / distance.
  assert(gradient.dimension(0) == distance_interface.dimension(0));
  gradient = invalid_value;
  const double indic_phase = (phase == 0) ? 0. : 1.;
  int i_elem;
  int err_count = 0;
  for (i_elem = 0; i_elem < nb_elem; i_elem++)
    {
      double d = distance_interface[i_elem];
      if (indicatrice[i_elem] == indic_phase && d > invalid_test)
        {

          // Test pour remedier a une eventuelle erreur de calcul de la
          // fonction distance dans les mailles voisines des mailles traversees
          // par l'interface. On sait que la distance est superieure a delta_x/2
          // A faire: calculer distance_min = delta_x/2
          //       if (std::fabs(d) < distance_min) {
          //         d = (d>0) ? distance_min : -distance_min;
          //         err_count++;
          //       }

          // Test pour savoir si la distance a l'interface est bien inferieure
          // au rayon du cercle inscrit a l'element
          // Si c'est le cas, on remplace la distance par plus ou moins ce rayon
          // suivant la phase dans laquelle se situe l'element
          double dist_elem_face_min = 1e30;
          for (int face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              double dist_elem_face = 0;
              const int face_glob = elem_faces(i_elem,face_loc);
              for (int i_dim=0; i_dim<Objet_U::dimension; i_dim++)
                {
                  double centre_elem_i = centre_gravite_elem(i_elem,i_dim);
                  double centre_face_i = centre_gravite_face(face_glob,i_dim);
                  dist_elem_face += (centre_elem_i-centre_face_i)*(centre_elem_i-centre_face_i);
                }
              dist_elem_face = sqrt(dist_elem_face);
              dist_elem_face_min = (dist_elem_face < dist_elem_face_min) ? dist_elem_face : dist_elem_face_min;
            }

          // Codage initial valable uniquement pour une discretisation VDF
          //       double dx = domaine_vf.volumes(i_elem)/domaine_vf.face_surfaces(elem_faces(i_elem,0));
          //       double dy = domaine_vf.volumes(i_elem)/domaine_vf.face_surfaces(elem_faces(i_elem,1));
          //       double dz = domaine_vf.volumes(i_elem)/domaine_vf.face_surfaces(elem_faces(i_elem,2));
          //       double dist_elem_face_min = ( ((dx<dy) ? dx : dy)<dz ? ((dx<dy) ? dx : dy) : dz )/2;

          // Si la distance entre le centre de l'element et l'interface est plus petite
          // que la plus petite distance entre le centre de l'element et le centre de ses faces
          // on ecrase la valeur de la distance a l'interface qui est manifestement fausse.
          // On choisit par defaut la plus petite distance entre le centre de l'element et le
          // centre de ses faces.  Le signe de la distance est determine en fonction de la phase
          // de l'element (la distance etant fausse, son signe a toutes les chances de l'etre
          // aussi).
          if (std::fabs(d) < dist_elem_face_min)
            {
              Cerr << "Time = " << temps << "; extrapoler_champ_elem: distance lower than dx/2" << finl;
              Cerr << "        Element position:" << finl;
              for  (int i_dim=0; i_dim<Objet_U::dimension; i_dim++)
                {
                  Cerr << "          x(" << i_dim << ") = " << centre_gravite_elem(i_elem,i_dim) << finl;
                }

              d = (phase == 0) ? -dist_elem_face_min : dist_elem_face_min;

              err_count++;
            }
          const double v = champ[i_elem];
          // Ceci est le gradient evalue a une distance d/2 de l'interface
          // (difference finie centree)
          const double grad = (v - interfacial_value) / d;
          // Correction du gradient pour trouver la valeur a l'interface :
          // on suppose que le transfert est stationnaire, et flux normal
          // a l'interface.
          const double div_n = champ_div_n[i_elem];
          gradient[i_elem] = grad * (1. + div_n * (d*0.5));
          //gradient[i_elem] = grad * (1. + div_n * (d*0.5) + 0.5 * div_n * div_n * (d*0.5) * (d*0.5));
        }
    }
  if (err_count)
    Cerr << "extrapoler_champ_elem: errcount=" << err_count << finl;
  gradient.echange_espace_virtuel();
  // On etale ce gradient par continuite sur une epaisseur de "stencil_value"
  // Les iterations de cet algorithme convergent vers une sorte de laplacien=0
  // avec condition aux limites de Dirichlet sur les elements de la phase "phase".
  extrapoler_dans_phase(gradient, indicatrice, domaine_vf, invalid_test, indic_phase, stencil_width/* nb_iter*/ );

  // 2023. Pour la vitesse de convection extrapolee, on reextrapole dans l'autre phase :
  const double autre_phase = 1.-indic_phase;
  extrapoler_dans_phase(gradient, indicatrice, domaine_vf, invalid_test, autre_phase, stencil_width/* nb_iter*/ );

  // On calcule la valeur extrapolee:
  for (i_elem = 0; i_elem < nb_elem; i_elem++)
    {
      const double d    = distance_interface[i_elem];
      const double grad = gradient[i_elem];
      if (indicatrice[i_elem] != indic_phase
          && d > invalid_test
          && grad > invalid_test)
        {
          // Extrapolation parabolique tenant compte de la courbure de l'interface :
          const double div_n = champ_div_n[i_elem];
          champ[i_elem] = d * grad * (1. - 0.5 * div_n * d) + interfacial_value;
          // TODO: Assess Higher order Taylor series by the inclusion of a new parameter in the datafile. (Formulae below to be checked)
          //  champ[i_elem] = 0.5 * d * grad * (1. - 0.5 * div_n * d) + interfacial_value; // multiplied by 0.5 (may be it was missing...we have to check)
          //  champ[i_elem] = d * grad * (1. - 0.5 * div_n * d + div_n * div_n * d * d / 6.) + interfacial_value;
        }
    }
  champ.echange_espace_virtuel();
}

/*! @brief met a jour le champ grad_t en fonction du champ inconnue.
 *
 * Attention, l'inconnue est modifiee (on etend le champ de temperature dans la phase
 *   opposee.
 *
 */
void Convection_Diffusion_Temperature_FT_Disc::calculer_grad_t()
{
  Transport_Interfaces_FT_Disc& eq_interface_ = ref_eq_interface_.valeur();

  const DoubleTab& indicatrice = eq_interface_.get_update_indicatrice().valeurs();
  const DoubleTab& distance_interface = eq_interface_.get_update_distance_interface().valeurs();
  const DoubleTab& normale_interface = eq_interface_.get_update_normale_interface().valeurs();
  //GB : Augmenter la constante de l'epaisseur
  //GB : const int stencil_width = 8;
  const int stencil_width = stencil_width_;
  const double interfacial_value = TSAT_CONSTANTE;

  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());

  const double temps = schema_temps().temps_courant();

  // Extrapolation lineaire de la temperature des mailles diphasiques a partir
  // de la temperature des mailles de la "phase".
  DoubleTab& temperature = inconnue()->valeurs();

  Navier_Stokes_FT_Disc& eq_navier_stokes = ref_cast(Navier_Stokes_FT_Disc, ref_eq_ns_.valeur());
  const DoubleTab& div_n = eq_navier_stokes.calculer_div_normale_interface().valeurs();

  extrapoler_champ_elem(domaine_vf, indicatrice, distance_interface, normale_interface, div_n,
                        phase_, stencil_width, interfacial_value,
                        temperature,
                        grad_t_->valeurs(),
                        temps);
}

void Convection_Diffusion_Temperature_FT_Disc::calculer_mpoint()
{
  calculer_mpoint(mpoint_);
}

void Convection_Diffusion_Temperature_FT_Disc::calculer_mpoint(Champ_base& mpoint)
{
  const double invalid_test = -1.e25;
  calculer_grad_t();

  if (is_prescribed_mpoint_)
    {
      mpoint.valeurs() = prescribed_mpoint_;
      return;
    }
  mpoint.valeurs() = grad_t_->valeurs();

  DoubleTab& val = mpoint.valeurs();
  const int n = val.size();
  for (int i = 0; i < n; i++)
    if (val[i] < invalid_test)
      val[i] = 0;

  const double k = fluide_dipha_->fluide_phase(phase_).conductivite()->valeurs()(0,0);
  const double L = fluide_dipha_->chaleur_latente();
  // L est la chaleur latente de changement de phase pour passer de
  // la phase 0 a la phase 1.
  double f = k / L;
  // Si on est dans la phase 0 et que L > 0, on doit avoir mpoint positif pour
  //  gradient(T) scalaire n negatif.
  if (phase_ == 0)
    f = -f;

  mpoint.valeurs() *= f;

  // Pour deverminage : on impose une variation de mpoint lineaire en z
  //const Domaine_VF & domaine_vf = ref_cast(Domaine_VF, domaine_dis());
  //for (int i = 0; i < n ; i++)
  //  mpoint.valeurs()[i] = 500./8957.*(1.+0.1*(domaine_vf.xp()(i,3)-0.0005)/0.001);

#if TCL_MODEL
  // Here, it is too early to correct the table mpoint.valeurs() because it has not been used yet to build the extended velocities.
  // If we were correcting it here with the TCL contribution, we would generate tangential artificial delta_u velocities that are not desired.
  // Consequently, the TCL effect will be explicitely applied to secmem2 and mpoint will be explicitely corrected too (but later,
  // just before post-processing in fact).
#endif
  mpoint.valeurs().echange_espace_virtuel();
}

// The list mixed_elems_ contains elems several times (once per operator conv/diff, once per face to a pure phase_ neighbour)
//                                                    (that is four times in most of 2D cells when convection+diffusion are applied)
// We reduce the list to unique occurences of elements.
static void collect_into_unique_occurence(ArrOfInt& mixed_elems, ArrOfDouble& lost_fluxes)
{
  const int nb_elem_with_duplicates = mixed_elems.size_array();
  int nb_elem = 0;
  // Cerr << "Algo may be optimized? It is in NxN = " << nb_elem_with_duplicates << " x "
  //     << nb_elem_with_duplicates << " = " << nb_elem_with_duplicates*nb_elem_with_duplicates <<finl;
  for (int i=0; i<nb_elem_with_duplicates; i++)
    {
      const int elemi =  mixed_elems[i];
      // Have we seen this element already?
      int j=0;
      for(j=0; j<i; j++)
        {
          const int elemj = mixed_elems[j];
          if (elemi == elemj)
            {
              // yes, then we hit the "continue" in the next if and the loop continues with the next element
              break;
            }
        }
      if (j!=i)
        continue;

      // Here, we are with a new element that will remain in the list (may be moved to the position nb_elem)
      lost_fluxes[nb_elem] = lost_fluxes[i]; // We store the first lost_flux for this elem
      for (j=i+1; j<nb_elem_with_duplicates; j++)
        {
          const int elemj = mixed_elems[j];
          if (elemi == elemj)
            {
              lost_fluxes[nb_elem] += lost_fluxes[j]; // and pile-up others...
            }
        }
      mixed_elems[nb_elem] = elemi;
      nb_elem++;
    }
  lost_fluxes.resize_array(nb_elem);
  mixed_elems.resize_array(nb_elem);
}

void Convection_Diffusion_Temperature_FT_Disc::correct_mpoint()
{
  Cerr << "Work in progress and widly incorrect. " << finl;
  Cerr << "Wait and see for further improvements. Exit" << finl;
  Process::exit();
  if (inconnue()->nb_valeurs_temporelles() == 1)
    {
      Cerr << "You need at least 2 positions to the wheel... Contact TRUST support. " << finl;
      Process::exit();
    }
  Transport_Interfaces_FT_Disc& eq_interface_ = ref_eq_interface_.valeur();
  //const Champ_base& ch_indic = ref_cast(Champ_Inc,eq_interface_.get_update_indicatrice());
  //const DoubleTab& indicatrice = ch_indic.valeurs();
  const DoubleTab& indicatrice = eq_interface_.inconnue()->valeurs();
  const DoubleTab& indicatrice_passe = eq_interface_.inconnue()->passe();
  const double& dt = schema_temps().pas_de_temps();
  //const DoubleTab& indicatrice_passe = ch_indic.passe();
  DoubleTab& temperature = inconnue()->valeurs();
  const DoubleTab& temperature_passe = inconnue()->passe();
  const double rhocp = fluide_dipha_->fluide_phase(phase_).masse_volumique()->valeurs()(0,0)
                       * fluide_dipha_->fluide_phase(phase_).capacite_calorifique()->valeurs()(0,0);

  const int nb_elem = mixed_elems_.size_array();
  //assert(mixed_elems_diffu_.size_array()==nb_elem);
  //assert(mixed_elems_conv_.size_array()==nb_elem); // all lists should now have the same size. Or maybe not due to BC?
  //                                                    But still, mixed_elems_ should be the longest and the other should be included
  derivee_energy_.resize_array(nb_elem);
  for(int i=0; i<nb_elem; i++)
    {
      const int elemi =  mixed_elems_[i];
      derivee_energy_[i] = (temperature[elemi] * indicatrice[elemi] - temperature_passe[elemi] * indicatrice_passe[elemi])* rhocp;
    }

  Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, ref_eq_ns_.valeur());
  const double Lvap  = fluide_dipha_->chaleur_latente();
  const DoubleVect& volume = ref_cast(Domaine_VF, domaine_dis()).volumes();
  DoubleTab& mp = mpoint_->valeurs();
  const DoubleTab& ai = ns.get_interfacial_area();
  const double temps = schema_temps().temps_courant();
  {
    double total_flux_lost = 0.;
    double total_derivee_energy = 0.;
    double mpai_tot = 0.;
    double mp_sum_before = 0.;
    double int_ai_before = 0.;
    for (int nd=0 ; nd<nb_elem ; nd++)
      {
        const int elembe = mixed_elems_[nd];
        int_ai_before += ai[elembe];
        //  Cerr << " elembe= " << elembe << finl;
        //   total_flux_lost -= lost_fluxes_(nd)*rhocp; // multiplied by rhocp (vp)  // "-" because depends on convention (GB)
        total_flux_lost -= lost_fluxes_(nd); // in the present TRUST version it should not be multiplied by rhocp
        total_derivee_energy += derivee_energy_(nd)*volume[nd];
        mpai_tot += mp[elembe]*ai[elembe]*Lvap;// multiplied by Lvap (vp)
        mp_sum_before += mp[elembe]*ai[elembe];
      }
    if (int_ai_before>DMINFLOAT)
      {
        const double mean_mp_before = mp_sum_before/int_ai_before;
        Cerr << " mp_sum_before= " << mp_sum_before << " mean_mp_before= " << mean_mp_before << " time= " << temps << finl;
      }
    total_flux_lost = mp_sum(total_flux_lost);
    mpai_tot = mp_sum(mpai_tot);
    total_derivee_energy = mp_sum(total_derivee_energy)/dt;
    Cerr << "[Basic-Mixed-cells-Energy-Balance] Time= " << temps << " nb_elems= " << nb_elem
         << " phi(positive if towards mixed cells)= " << total_flux_lost
         << " mp*ai*Lvap= " << mpai_tot
         << " dE/dt= " << total_derivee_energy
         << " imbalance= " << total_derivee_energy-total_flux_lost+mpai_tot << finl;
  }

  // Energy balance correction. Loop on mixed elems only :
  const int option=-1; // To disable that
  {
    const int account_for_diff = correction_mpoint_diff_conv_energy_[0];
    const int account_for_conv = correction_mpoint_diff_conv_energy_[1];
    const int account_for_mixed_cell_energy = correction_mpoint_diff_conv_energy_[2];
    double int_dmp_ai = 0.; // The integral over the interface of delta_mp
    double int_ai = 0.; // The interface area
    double int_mp_ai = 0.;
    for(int i=0; i<nb_elem; i++)
      {
        const int elem = mixed_elems_[i];
        // The convention is that phi_lost is viewed from the liquid side point-of-view.
        // (negative for evap as it's leaving the pure liquid)
        // So when we consider mixed cells, the incoming flux is "-phi_lost"
        const double phi_in_mixed_cell = lost_fluxes_[i];
        double phi_conv_lost_by_mixed_cell = 0.;
        double phi_diffu_added_to_mixed_cell = 0.;
        if (mixed_elems_diffu_[i] != elem)
          {
            Cerr << "Search for a solution? What case is it? diffu, BC? " << finl;
            Process::exit();
          }
        else
          {
            phi_diffu_added_to_mixed_cell = lost_fluxes_diffu_[i];
            //   Cerr << " lost-flux-diff= " << phi_diffu_added_to_mixed_cell << " at-time= " << temps << " in-elem= " << elem << finl;
          }
        if (i>=mixed_elems_conv_.size_array() || mixed_elems_conv_[i] != elem)
          {
            phi_conv_lost_by_mixed_cell =0.;
            int j=0;
            for (j=0; j<mixed_elems_conv_.size_array(); j++)
              {
                if (mixed_elems_conv_[j] == elem)
                  {
                    phi_conv_lost_by_mixed_cell = lost_fluxes_conv_[j];
                    break;
                  }
              }
            if (j==mixed_elems_conv_.size_array())
              {
                Cerr << "conv. The end of the list is reached, "
                     << "mixed_elems_["<< i<<"]= " << elem << " was not found in mixed_elem_conv_"<< finl;
                Cerr << "mixed_elems_conv_= " << mixed_elems_conv_ << finl;
                Cerr << "mixed_elems_= " << mixed_elems_ << finl;
                Cerr << "WE ASSUME IT IS BECAUSE IT IS A BC? any solution? WE IGNORE IT"<< finl;
                // Process::exit();
              }

          }
        else
          {
            phi_conv_lost_by_mixed_cell = lost_fluxes_conv_[i];
            //     Cerr << " lost-flux-conv= " << phi_conv_lost_by_mixed_cell << " at-time= " << temps << " in-elem= " << elem << finl;
          }

        double VdrhocpT_dt = 0.;
        if (account_for_mixed_cell_energy)
          {
            VdrhocpT_dt=(rhocp*volume[elem])*(temperature[elem] * indicatrice[elem] - temperature_passe[elem] * indicatrice_passe[elem])/dt;
            //   Cerr << "temp_current= " << temperature[elem] << " temp_previous= " << temperature_passe[elem] << finl;
          }

        if (option ==1)
          {
            // OPtion 1 : Correction of T
            const double value_before = temperature[elem];
            temperature[elem] = 1/indicatrice[elem] * (temperature_passe[elem] * indicatrice_passe[elem] + dt/(rhocp*volume[elem])*(mp[elem]*ai[elem]*Lvap - phi_in_mixed_cell*rhocp));
            // multiplied phi_lost by rhocp (vp)
            Cerr << "[Delta-T-due-to-Energy-correction] elem "<< elem <<" Tnew-Tuncorrected " << temperature[elem]-value_before << finl;
            //
            const double tmp = (temperature[elem] * indicatrice[elem] - temperature_passe[elem] * indicatrice_passe[elem])* rhocp;
            if (derivee_energy_[i] != tmp)
              {
                Cerr << "New derivee_energy after correction (before/after) : ( "<< derivee_energy_[i] <<" / " << tmp << " )." << finl;
                derivee_energy_[i] = tmp;
              }
          }
        else if (option ==2)
          {
            // Option 2 : Correction of mp: local corrections in all mixed cells
            if (ai[elem]>DMINFLOAT)
              {
                const double old_mp = mp[elem];

                double delta_mp = 0. ;
                if (account_for_diff)
                  {
                    delta_mp += (1/(ai[elem]*Lvap))*((mp[elem]*ai[elem]*Lvap + phi_diffu_added_to_mixed_cell)); // for Stefan it is negative sign before phi_difu_...
                    // (we need to make it uniform so that it works for all with the same sign convention....may depend on the normal)
                    Cerr << " delta-mp-diff= " << delta_mp << finl;
                  }
                if (account_for_conv)
                  {
                    // To be checked :
                    //   - the minus sign?
                    //  - what surface should we use? ai or Sface?
                    //  - The interpolation of T to the face?
                    //  double Tface = temperature[elem];
                    //   double rhov = fluide_dipha_->fluide_phase(1-phase_).masse_volumique()(0,0);

                    delta_mp += (1/(ai[elem]*Lvap))*(-phi_conv_lost_by_mixed_cell);
                    //         delta_mp += (1/(ai[elem]*Lvap))*(-(mp[elem]*ai[elem]*rhocp/rhov*Tface - phi_conv_lost_by_mixed_cell));
                    Cerr << " delta-mp-conv= " << delta_mp << finl;
                  }
                if (account_for_mixed_cell_energy)
                  delta_mp += (1/(ai[elem]*Lvap))*(VdrhocpT_dt);
                // divided by Lvap (now consistent in units-vp) and multiplied phi_lost by rhocp
                mp[elem] +=delta_mp;
                Cerr << " delta_mp= " << delta_mp << " old_mp= " << old_mp << " new_mp= " << mp[elem] << finl;
                if (fabs(old_mp)>DMINFLOAT)
                  Cerr << "relative correction of mp at time "<< temps << " elem= " << elem << " delta= " <<delta_mp/old_mp*100. << "%." <<finl;
              }
          }
        else if (option==3)
          {
            // Option 3 : Global Correction of mp:
            if ((ai[elem]>DMINFLOAT) and (temps>DMINFLOAT))
              {
                int_ai += ai[elem];
                int_mp_ai += mp[elem]*ai[elem];
                //  Cerr << " area-elem= " << ai[elem] << " temps= " << temps << finl;
                // Before we had : const double tmp = (1/(Lvap))*(VdrhocpT_dt-(mp[elem]*ai[elem]*Lvap - phi_in_mixed_cell*rhocp));

                double delta_mp = 0. ;
                if (account_for_diff)
                  {
                    //  delta_mp += (1/(ai[elem]*Lvap))*(-(mp[elem]*ai[elem]*Lvap - phi_diffu_added_to_mixed_cell));
                    delta_mp += (-1/Lvap)*((mp[elem]*ai[elem]*Lvap - phi_diffu_added_to_mixed_cell));// for Stefan it is negative sign before phi_difu_...
                    // (we need to make it uniform so that it works for all with the same sign convention....may depend on the normal)
                    //  delta_mp += mp[elem]*ai[elem] + (1/Lvap)*(phi_diffu_added_to_mixed_cell);
                    //                Cerr << " mpaiL= " << mp[elem]*ai[elem]*Lvap << " phi_diff= " << phi_diffu_added_to_mixed_cell << finl;
                    //               Cerr << " delta-mp-diff= " << delta_mp << finl;
                  }

                if (account_for_conv)
                  {
                    // To be checked :
                    //   - the minus sign?
                    //  - what surface should we use? ai or Sface?
                    //  - The interpolation of T to the face?
                    //   double Tface = temperature[elem];
                    //   double rhov = fluide_dipha_->fluide_phase(1-phase_).masse_volumique()(0,0);

                    //  delta_mp += (1/(ai[elem]*Lvap))*(-(mp[elem]*ai[elem]*rhocp/rhov*Tface - phi_conv_lost_by_mixed_cell));
                    // delta_mp += (1/(ai[elem]*Lvap))*((mp[elem]*ai[elem]*rhocp/rhov*Tface - phi_conv_lost_by_mixed_cell));
                    // delta_mp += (1/Lvap)*(mp[elem]*ai[elem]*rhocp/rhov*Tface);
                    // Cerr << " conv-new= " << (1/Lvap)*(mp[elem]*ai[elem]*rhocp/rhov*Tface) << finl;
                    delta_mp += (-1/Lvap)*(phi_conv_lost_by_mixed_cell);
                    //  Cerr << " conv-new= " << phi_conv_lost_by_mixed_cell << " in elem= " << elem << finl;
                    //  Cerr << " delta-mp-conv= " << delta_mp << finl;
                  }
                if (account_for_mixed_cell_energy)
                  {
                    delta_mp += (1/Lvap)*(VdrhocpT_dt);
//                   Cerr << " delta-mp-energy= " << (1/Lvap)*(VdrhocpT_dt) << finl;
                  }
                int_dmp_ai += delta_mp;
//               Cerr << " delta-mp-tot= " << int_dmp_ai << " time= " << temps << finl;
              }
          }
        // End of options
      }
    if ((option==3) and (int_ai>DMINFLOAT))
      //  if ((option==3))
      {
        const double mean_dmp = int_dmp_ai/int_ai;
        const double mean_mp = int_mp_ai/int_ai;
        double rel= 0.;
        if (fabs(mean_mp)>DMINFLOAT)
          rel=mean_dmp/mean_mp*100.;
        Cerr << "correction of mp at time "<< temps << " mean_mp= " << mean_mp << " relative= " << rel <<  "%." <<finl;
        /*  Bad correction:
         *
         for(int i=0; i<nb_elem; i++)
            {
               const int elem = mixed_elems_[i];
               mp[elem] +=mean_dmp;
               Cerr << " newmp= " << mp[elem] << finl;
             }*/
        // New correction truely based on AI:
        // (to be in perfect agreement with the condition (on untouched variables) used later for the extrapolation:
        /* for(int i=0; i<nb_elem; i++)
          {
            const int elem = mixed_elems_[i];
            if ((ai[elem]>DMINFLOAT) and (temps>DMINFLOAT))
                mp[elem] +=mean_dmp;
          } */

        // GB 18/10/2020. Application of the correction EVERYWHERE:
        mp +=mean_dmp;
      }

    double total_flux_lost = 0.;
    double total_derivee_energy = 0.;
    double mpai_tot = 0.;
//   double mp_tot = 0.;
    for (int nd=0 ; nd<nb_elem ; nd++)
      {
        total_flux_lost -= lost_fluxes_(nd)*rhocp; // multiplied by rhocp (vp) // "-" because depends on convention (GB)
        total_derivee_energy += derivee_energy_(nd)*volume[nd];
        const int elem = mixed_elems_[nd];
        mpai_tot += mp[elem]*ai[elem]*Lvap; // multiplied by Lvap (vp)
        //  mp_tot += mp[elem];
      }
    if ((option==3) and (int_ai>DMINFLOAT))
      {
        const double mean_mp_corr = mpai_tot/(int_ai*Lvap);
        Cerr << " mean_mp_corr= " << mean_mp_corr << " time= " << temps << finl;
      }
    total_flux_lost = mp_sum(total_flux_lost);
    mpai_tot = mp_sum(mpai_tot);
    total_derivee_energy = mp_sum(total_derivee_energy)/dt;
    Cerr << "[Corrected-balance] Time= " << temps << " nb_elems= " << nb_elem
         << " phi(positive if towards mixed cells)= " << total_flux_lost
         << " mp*ai*Lvap= " << mpai_tot
         << " dE/dt= " << total_derivee_energy
         << " imbalance= " << total_derivee_energy-total_flux_lost+mpai_tot << finl;
  }

  // If mpoint has been modified in interfacial cells (mixed cells). We want to extend it back into both phases
  // using the same procedure
  if ((option==0) or (option==3) or (option==2))
    {
      double mmax = 0.;
      const int n = ai.size_array();
      for(int i=0; i<n; i++)
        {
          if ((ai[i]>DMINFLOAT) && (fabs(mp[i])>mmax))
            {
              mmax =fabs(mp[i]);
            }
        }
      mmax = Process::mp_max(mmax);
      //  Cerr << "[Maximum-mp] Time= " << temps << " max(abs(mp))= " << mmax << finl;
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
      const DoubleTab& distance_interface = eq_interface_.get_update_distance_interface().valeurs();
      extrapolate(domaine_vf, ai, stencil_width_, distance_interface, mp);
    }
}

// With this approach, calculer_delta_u_interface is not called,
// so the velocity should be extended.
// Fills the field : vitesse_convection_
void Convection_Diffusion_Temperature_FT_Disc::compute_divergence_free_velocity_extension()
{
  Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, ref_eq_ns_.valeur());
  vitesse_convection_->valeurs() = ns.inconnue()->valeurs();
  // Projection of the convective field :
  //SolveurSys solveur_pression(ns.get_solveur_pression());
  Solveur_Masse solveur_masse_fictitious(ns.solv_masse()); // Copy the operator to change the coeff
  solveur_masse_fictitious->set_name_of_coefficient_temporel("no_coeff");

  //On utilise un operateur de divergence temporaire et pas celui porte par l equation
  //pour ne pas modifier les flux_bords_ rempli au cours de ...::mettre_a_jour
  Operateur_Div div_tmp;
  div_tmp.associer_eqn(ns);// It's slightly tricky but associating it to (*this) is not OK, because we want to apply it to a velocity.
  // As delta_u is not the unknown of this equation (no more than vitesse_convection_), we resort to ns to
  // make it work. Basically, it also means that delta_u gets the BC from u_ns?
  // Problem with div_tmp.typer() that get the type of equation to make a type.
  // Hence, it will see a conv/diff equation for a scalar... and will type div_tmp badly (not for a vector velocity field!)
  // Instead we associer_eqn to ns
  div_tmp.typer();
  Cerr << "[Conv_diff_FT_Disc] div_tmp recieved the type " << div_tmp.type() << finl;
  Cerr << "[Conv_diff_FT_Disc] div_tmp que suis je " << div_tmp.que_suis_je() << finl;
  if (0)
    {
      Equation_base& eqn=ref_eq_ns_.valeur();
      Nom inut;
      Nom nom_type=eqn.discretisation().get_name_of_type_for(ns.que_suis_je(),inut,ns);
      //  const Probleme_FT_Disc_gen& pb = ref_cast(Probleme_FT_Disc_gen,probleme());
      // ref_cast
      //  OWN_PTR(Operateur_Div_base)::typer(nom_type);
      Cerr << "[Conv_diff_FT_Disc] Velocity correction. Construction of the divergence operator type : " << nom_type << finl;
      //   Cerr << valeur().que_suis_je() << finl ;
    }

  div_tmp.l_op_base().associer_eqn(ns);//*this);
  //div_tmp.typer();
  div_tmp->completer();
  // En VDF, la methode 'completer()' est surchargee et fait
  //     1. Operateur_base::completer();
  //     2. iter.completer_();
  // donc si on fait juste :
  //     div_tmp.l_op_base().associer(domaine_dis(), zcl_fictitious_, vitesse_convection_);
  // On rate l'etape 2 qui a ete faite avec le completer plus haut...

  // WARNING : Si on prend les cl de l'ope de pression pour cette divergence, il y a une vitesse sur les sorties de NS
  //           qui se trouve face a des conditions limites de symetrie dans zcl_fictitious...
  //           Du coup, ca fait un gros div! Ce n'est pas ce qu'on veut.
  //           PAR CONTRE, pour le calcul du gradient ci-dessous, c'est bien les cl de zcl_fictitious qu'on veut,
  //           sinon, sortie libre va mettre un grad(P).n non nul au bord!!
  //		   DONC EN CONCLUSION, ON VEUT UN MIX!
  //

  // RHS for this pressure solveur : div(delta_u)
  // We need to create a table sized as temperature to store the rhs :
  DoubleTab& secmem = divergence_delta_U->valeurs();
  //secmem.copy(inconnue()->valeurs(), RESIZE_OPTIONS::NOCOPY_NOINIT);
  div_tmp->calculer(vitesse_convection_->valeurs(),secmem);

  // On ne conserve que la divergence des elements proches de l'interface, et supprime quand la distance est invalide
  if (0) // Cela pourrait etre utile si on construisait le saut de vitesse delta u_0 mais
    // cela semble presenter que peu d'interet...
    {
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
      const IntTab& face_voisins = domaine_vf.face_voisins();
      const IntTab& elem_faces = domaine_vf.elem_faces();
      const int nb_elem = secmem.size_array();
      const int   nb_faces_elem = elem_faces.dimension(1);
      const Transport_Interfaces_FT_Disc& eq_transport = ref_eq_interface_.valeur();
      // Distance a l'interface discretisee aux elements:
      const DoubleTab& distance = eq_transport.get_update_distance_interface().valeurs();
      //     const int nb_elem = secmem2.dimension(0);
      for (int elem = 0; elem < nb_elem; elem++)
        {
          const double dist = distance(elem);
          int i_face = -1;
          if (dist < -1e20)
            {
              // Distance invalide: on est loin de l'interface
              // invalide ou voisin invalide, on supprime la div :
              secmem(elem) = 0.;
            }
          else
            {
              // Y a-t-il un voisin pour lequel la distance est invalide
              for (i_face = 0; i_face < nb_faces_elem; i_face++)
                {
                  const int face = elem_faces(elem, i_face);
                  const int voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                  if (voisin >= 0)
                    {
                      const double d = distance(voisin);
                      if (d < -1e20)
                        {
                          // invalide ou voisin invalide, on supprime la div :
                          secmem(elem) = 0.;
                          break; // Un voisin invalide
                        }
                    }
                }
            }
        }
    }
  // Correction du second membre d'apres les conditions aux limites :
  assembleur_pression_->modifier_secmem(secmem);
  // Ajout pour la sauvegarde au premier pas de temps si reprise
  la_pression->changer_temps(schema_temps().temps_courant());

  // Resolution du systeme en pression : calcul de la_pression
  solveur_pression_.resoudre_systeme(matrice_pression_.valeur(),
                                     secmem,
                                     la_pression->valeurs()
                                    );
  assembleur_pression_->modifier_solution(la_pression->valeurs());
  // Calcul d(u)/dt = vpoint + 1/rho*grad(P)

  DoubleTab& gradP = gradient_pression_->valeurs();
  // Can I re-use a gradient from NS?
  Operateur_Grad gradient_tmp;
  gradient_tmp.associer_eqn(ns);//*this);
  gradient_tmp.typer();
  gradient_tmp.l_op_base().associer_eqn(ns);//*this);
  gradient_tmp->completer();
  // It is important to associate the gradient to the good BC via zcl,
  // Otherwise a normal gradient appears at the outlet.
  gradient_tmp.l_op_base().associer(domaine_dis(), zcl_fictitious_, vitesse_convection_); // Quel champ inco pour le gradP? A quoi sert-elle?
  gradient_tmp.calculer(la_pression->valeurs(), gradP);
  // I don't wanna divide by rho_face :
  solveur_masse_fictitious->appliquer(gradP); // divide by cell_volume

  // Correction of vitesse : "+=" seems the good sign, though I don't understand why...
  {
    int i, j;
    DoubleTab& vc = vitesse_convection_->valeurs();
    const int n = vc.dimension(0);
    if (vc.nb_dim() == 1)
      {
        // VDF
        for (i = 0; i < n; i++)
          vc(i) += gradP(i);
      }
    else
      {
        //VEF
        const int m = vc.dimension(1);
        for (i = 0; i < n; i++)
          for (j = 0; j < m; j++)
            vc(i,j) += gradP(i,j);
      }
    vc.echange_espace_virtuel();
  }

}
DoubleTab& Convection_Diffusion_Temperature_FT_Disc::derivee_en_temps_inco(DoubleTab& derivee)
{
  // We start the timestep by computing :
  // STEP 1. Extend velocity and temperature (via calculer_grad_t which does a linear extrapolation
  // from the values within "phase_")
  static const Stat_Counter_Id count2 = statistiques().new_counter(1, "calculer_mpoint", 0);
  statistiques().begin_count(count2);
  calculer_mpoint();
  statistiques().end_count(count2);

  Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, ref_eq_ns_.valeur());
  // Emptying lists for new timestep.
  mixed_elems_.resize_array(0);
  lost_fluxes_.resize_array(0);

  // STEP 2. Compute the continuous extension of velocity from phase_ into the other.
  // This velocity will be usefull for the convection operator.
  // Now, the step calculer_delta_u_interface calls get_mpoint instead of
  // calculer_mpoint (which used to in turns call: calculer_grad_t)
  // Therefore, it is important to have called explicitely the temperature
  // extension before at the begining of the function.
  static const Stat_Counter_Id count1 = statistiques().new_counter(1, "extend_ui", 0);
  statistiques().begin_count(count1);
  const Champ_Inc& vitesse_ns = ns.inconnue();
  if (!divergence_free_velocity_extension_)
    {
      ns.calculer_delta_u_interface(vitesse_convection_, phase_, correction_courbure_ordre_);
      //vitesse_convection_.valeurs() += vitesse_ns.valeurs();
      vitesse_convection_->valeurs() += vitesse_ns->futur(); // avec les jeux d'avancer reculer et les equations, il faut prendre futur pour avoir u^n ici
      //Cerr << inconnue()->temps()  << " =! " << vitesse_ns.temps() << " " << vitesse_convection_.temps() << finl;
      //Process::exit();
    }
  else
    {
      // Compute vitesse_convection_ which is equal to ns.inconnue() in phase_ and is extended (divergence-free) in the other phase:
      compute_divergence_free_velocity_extension();
    }
  statistiques().end_count(count1);

  // Application of 2 operators: convection & diffusion
  // STEP 3: Convection operator
  // STEP 4: Diffusion operator
  // Convection and diffusion are dealt with in the standard way;
  // diffusion can be implicited.
  // rhoCp is used for convection.
  // vitesse_convection_ is used for convection.
  DoubleTab& temperature = inconnue()->valeurs();
  derivee = 0.;

  // STEP 3: Diffusion operator
  //  static const Stat_Counter_Id count2 = statistiques().new_counter(1, "terme_diffusif_T", 0);
  //  statistiques().begin_count(count2);
  bool calcul_explicite = false;
  if (parametre_equation_.non_nul() && sub_type(Parametre_implicite, parametre_equation_.valeur()))
    {
      Parametre_implicite& param2 = ref_cast(Parametre_implicite, parametre_equation_.valeur());
      calcul_explicite = param2.calcul_explicite();
    }

  if (schema_temps().diffusion_implicite() && !calcul_explicite)
    {
      // do it later?
    }
  else
    {
      terme_diffusif.ajouter(temperature, derivee);
    }
  //  statistiques().end_count(count2);
  const int nb_diffu=mixed_elems_.size_array();
  // const double temps = schema_temps().temps_courant();
  lost_fluxes_diffu_.resize_array(nb_diffu);
  mixed_elems_diffu_.resize_array(nb_diffu);
  {
    double total_flux_lost = 0.;
    for (int nd=0 ; nd<nb_diffu ; nd++)
      {
        total_flux_lost += lost_fluxes_(nd);
        lost_fluxes_diffu_(nd) = lost_fluxes_(nd);
        mixed_elems_diffu_(nd) = mixed_elems_(nd);
      }
    total_flux_lost = mp_sum(total_flux_lost);
    //  Cerr << "[Lost-fluxes-diffusif] Time= " << temps << " nb_faces= " << nb_diffu
    //     << " Lost= " << total_flux_lost << finl;
  }

  // STEP 4: Convection operator
  //  static const Stat_Counter_Id count3 = statistiques().new_counter(1, "convection_T", 0);
  //  statistiques().begin_count(count3);
  {
    const DoubleTab& rhoCp = get_champ("rho_cp_comme_T").valeurs();
    DoubleTab derivee_tmp(derivee);
    derivee_tmp = 0.;
    terme_convectif.ajouter(temperature, derivee_tmp);
    derivee_tmp *= rhoCp;
    derivee += derivee_tmp;
#if TCL_MODEL
    const int nb_conv=mixed_elems_.size_array()-nb_diffu;
    lost_fluxes_conv_.resize_array(nb_conv);
    mixed_elems_conv_.resize_array(nb_conv);
    double total_flux_conv_lost = 0.;
    for (int nd=0 ; nd<nb_conv ; nd++)
      {
        const int elem = mixed_elems_(nb_diffu+nd);
        const double flux = lost_fluxes_(nb_diffu+nd)*rhoCp[elem];
        total_flux_conv_lost += flux;
        lost_fluxes_conv_(nd) = flux;
        mixed_elems_conv_(nd) = elem;
      }

    total_flux_conv_lost = mp_sum(total_flux_conv_lost);
    // Cerr << "[Lost-fluxes-convection] Time= " << temps << " nb_faces= " << nb_conv
    //    << " Lost= " << total_flux_conv_lost << finl;
#endif
  }
  //  statistiques().end_count(count3);

  solveur_masse->appliquer(derivee);
  if (schema_temps().diffusion_implicite() && !calcul_explicite)
    {
      Convection_Diffusion_Temperature::derivee_en_temps_inco(derivee);
      return derivee;
    }
  else
    {
      // TODO : J'aimerais bien deplacer la diffusion explicite la aussi.
      //        Il faut juste revoir les boucles pour deplacer le lost_fluxes associes ici aussi
      // terme_diffusif.ajouter(temperature, derivee);
    }

  // To remove duplicates in the list of mixed_elems (some elements are there several times, twice (conv+diff) for each face-to-pure-liquid)
  collect_into_unique_occurence(mixed_elems_, lost_fluxes_);
  collect_into_unique_occurence(mixed_elems_diffu_, lost_fluxes_diffu_);
  collect_into_unique_occurence(mixed_elems_conv_, lost_fluxes_conv_);

  {
    const int nb=mixed_elems_.size_array();
    double total_flux_lost = 0.;
    for (int nd=0 ; nd<nb ; nd++)
      total_flux_lost += lost_fluxes_(nd);

    total_flux_lost = mp_sum(total_flux_lost);
    //  Cerr << "[Lost-fluxes-conv+diff] Time= " << temps << " nb_faces= " << nb
    //     << " Lost= " << total_flux_lost << finl;
  }

  // STEP 5: (OPTIONNAL?) Attempt to correct mpoint to reach an energy conservation.
  // Apply a correction to mpoint in order to locally satisfy the energy balance in mixed cells (consistency between incoming fluxes that are lost
  // and the phase change mp*ai, all with regards to the liquid energy evolution in the cell: d(rho*cp*chi*T)/dt.
  {
    mpoint_uncorrected_->valeurs() = mpoint_->valeurs();
    DoubleTab& mp = mpoint_->valeurs();
    const DoubleTab& ai = ns.get_interfacial_area();
    const double temps = schema_temps().temps_courant();
    const int nb_elemi = mixed_elems_.size_array();
    double mp_sum_before_corr = 0.;
    double int_ai_before = 0.;
    for (int nd=0 ; nd<nb_elemi ; nd++)
      {
        const int elembi = mixed_elems_[nd];
        int_ai_before += ai[elembi];
        //     total_flux_lost -= lost_fluxes_(nd)*rhocp; // multiplied by rhocp (vp)  // "-" because depends on convention (GB)
        //     total_derivee_energy += derivee_energy_(nd)*volume[nd];
        //     mpai_tot += mp[elemb]*ai[elemb]*Lvap;// multiplied by Lvap (vp)
        //    Cerr << " ai= " << ai[elembi] << " int_ai= " << int_ai_before << finl;
        mp_sum_before_corr += mp[elembi]*ai[elembi];
      }
    if (int_ai_before>DMINFLOAT)
      {
        const double mean_mp_before_corr = mp_sum_before_corr/int_ai_before;
        Cerr << " mp_sum_before_corr= " << mp_sum_before_corr << " mean_mp_before_corr= " << mean_mp_before_corr << " time= " << temps << finl;
      }
    if ((!is_prescribed_mpoint_) && (max_array(correction_mpoint_diff_conv_energy_)))
      {
        // WARNING!!  We don't correct mpoint if correction_mpoint_diff_conv_energy_ are all set to 0
        // this method is not related to TCL model but rather to attempting energy conserving Ghost Fluid Method.
        correct_mpoint();
      }
  }

  const Probleme_FT_Disc_gen& pb = ref_cast(Probleme_FT_Disc_gen,probleme());
  const Triple_Line_Model_FT_Disc& tcl = pb.tcl();
  // const double max_val_before = max_array(temperature);
  // const double min_val_before = min_array(temperature);
  if (tcl.is_activated())
    {
      // No need to recompute TCL here, it should be up-to-date.
#ifdef DEBUG
      {
        const Transport_Interfaces_FT_Disc& eq_transport = ref_eq_interface_.valeur();
        const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
        assert(tcl.tag_tcl() == maillage.get_mesh_tag());
        toto ce code n est pas lu
      }
#endif
      // Ne compile pas, pourquoi?
      // assert(tcl.tag_tcl() == ref_eq_interface_->maillage_interface().get_mesh_tag());

      // We are late enough within the timestep for the correction not to affect the velocities reconstructions (delta_u)
      // because we expect it has been done before. Nonetheless, we need to do it before the post-pro to see our corrected fields
      // in the lata
      //
      // Correct the field mpoint in wall-adjacent cells to account for TCL model:
      // ---> It cannot be done before because it would disturbs the solution via delta u!
      tcl.corriger_mpoint(mpoint_->valeurs());

      // Correct the phase-change in wall adjacent cells:
      // The mean cell-temperature is simply derived from the TCL solution.
      tcl.set_wall_adjacent_temperature_according_to_TCL_model(temperature);
      temperature.echange_espace_virtuel();
    }
  // Cerr << "[temperature] max before/after TCL model: " << max_val_before << " / " << max_array(temperature) << finl;
  // Cerr << "[temperature] min before/after TCL model: " << min_val_before << " / " << min_array(temperature) << finl;

  return derivee;
}

void Convection_Diffusion_Temperature_FT_Disc::mettre_a_jour (double temps)
{
  Convection_Diffusion_Temperature::mettre_a_jour(temps);
  vitesse_convection_->mettre_a_jour(temps);
  // GB : Debut du maintien artificiel de la temperature.
  if (maintien_temperature_)
    {
      const Nom nom_sous_domaine = nom_sous_domaine_;
      const double temp_moy_ini = temp_moy_ini_;
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
      const Domaine& dom = domaine_vf.domaine();
      const Sous_Domaine& ss_domaine = dom.ss_domaine(nom_sous_domaine);

      Transport_Interfaces_FT_Disc& eq_interface_ = ref_eq_interface_.valeur();
      const DoubleTab& indicatrice = eq_interface_.get_update_indicatrice().valeurs();
      DoubleTab& temperature = inconnue()->valeurs();
      const DoubleVect& volume = ref_cast(Domaine_VF, domaine_dis()).volumes();
      const int nb_elem = domaine_vf.domaine().nb_elem();

      // Calcul de la moyenne sous domaine, et du facteur : fac = temp_moy_ini/temp_moy_ss_domaine
      double temp_moy_ss_domaine = 0.;
      double vol_liq = 0.;
      for(int i=0; i<ss_domaine.nb_elem_tot() /*methode d'acces au nombre d'elements*/; i++)
        {
          int index = ss_domaine[i];
          // assert(index < domaine_vf.domaine().nb_elem());
          if (index < nb_elem)
            {
              temp_moy_ss_domaine += temperature[index] * indicatrice[index] * volume[index];
              vol_liq += indicatrice[index] * volume[index];
            }
        }
      vol_liq = mp_sum(vol_liq);
      temp_moy_ss_domaine = mp_sum(temp_moy_ss_domaine);
      if (vol_liq == 0. || temp_moy_ss_domaine == 0.)
        {
          // ne fien faire
        }
      else
        {
          temp_moy_ss_domaine /= vol_liq ;
          double fac = temp_moy_ini / temp_moy_ss_domaine;

          // Correction du champ de temperature :
          temperature *= fac; // Methode de multiplication d'un tableau
        }
    }
  // GB : Fin.
}

void Convection_Diffusion_Temperature_FT_Disc::discretiser()
{
  phase_ = 1;

  const Discretisation_base& dis = discretisation();
  const double temps = schema_temps().temps_courant();
  const Domaine_dis_base& un_domaine_dis = domaine_dis();
  LIST(REF(Champ_base)) & champs_compris = liste_champs_compris_;
  const int nb_valeurs_temps = schema_temps().nb_valeurs_temporelles();

  Nom nom;

  nom = Nom("temperature_") + le_nom();
  dis.discretiser_champ("temperature", un_domaine_dis, nom, "K", 1 /* composantes */, nb_valeurs_temps, temps, la_temperature);
  champs_compris.add(la_temperature.valeur());
  champs_compris_.ajoute_champ(la_temperature);

  nom = Nom("temperature_grad_") + le_nom();
  dis.discretiser_champ("temperature", un_domaine_dis, nom, "K/m", 1 /* composantes */, temps, grad_t_);
  champs_compris.add(grad_t_.valeur());
  champs_compris_.ajoute_champ(grad_t_);

  nom = Nom("mpoint_") + le_nom();
  dis.discretiser_champ("temperature", un_domaine_dis, nom, "kg/(m2s)", 1 /* composante */, temps, mpoint_);
  champs_compris.add(mpoint_.valeur());
  champs_compris_.ajoute_champ(mpoint_);

  nom = Nom("mpoint_uncorrected_") + le_nom();
  dis.discretiser_champ("temperature", un_domaine_dis, nom, "kg/(m2s)", 1 /* composante */, temps, mpoint_uncorrected_);
  champs_compris.add(mpoint_uncorrected_.valeur());
  champs_compris_.ajoute_champ(mpoint_uncorrected_);

  nom = Nom("vitesse_conv_") + le_nom();
  dis.discretiser_champ("vitesse", un_domaine_dis, nom, "m/s", -1 /* nb composantes par defaut */, 1 /* valeur temporelle */, temps, vitesse_convection_);
  champs_compris.add(vitesse_convection_.valeur());
  champs_compris_.ajoute_champ(vitesse_convection_);

  // TODO: I'd like the flag to work... But it is not updated yet. Please Help me, at Adrien Bruneton
  if (1 || divergence_free_velocity_extension_) // Problem for ABN : On n'a pas encore lu le param.ajouter qui passe mon flag a 1.
    {
      Cerr << "Fake Pressure gradient discretization" << finl;
      nom = Nom("gradient_pression_") + le_nom();
      dis.discretiser_champ("vitesse", un_domaine_dis,
                            nom, "",
                            -1 /* nb composantes par defaut */, 1 /* valeur temporelle */, temps,
                            gradient_pression_);
      champs_compris.add(gradient_pression_.valeur());
      champs_compris_.ajoute_champ(gradient_pression_);

      Cerr << "Fake Pressure discretization" << finl;
      nom = Nom("fake_pressure_") + le_nom();
      dis.discretiser_champ("pression",un_domaine_dis,nom,"Pa.m3/kg",1,1,temps,la_pression);
      champs_compris.add(la_pression.valeur());
      champs_compris_.ajoute_champ(la_pression);

      Cerr << "Velocity (delta_u) divergence discretization" << finl;
      nom = Nom("divergence_delta_U_") + le_nom();
      dis.discretiser_champ("divergence_vitesse" /*directive */,un_domaine_dis, nom, "m3/s", 1,1,temps,divergence_delta_U);
      champs_compris.add(divergence_delta_U.valeur());
      champs_compris_.ajoute_champ(divergence_delta_U);

      discretiser_assembleur_pression();
    }
  Equation_base::discretiser();
}

void Convection_Diffusion_Temperature_FT_Disc::discretiser_assembleur_pression()
{
  Nom type = "Assembleur_P_";
  type += discretisation().que_suis_je();
  //type += "_homogene";
  Cerr << "Navier_Stokes_std::discretiser_assembleur_pression : type="<< type << finl;
  assembleur_pression_.typer(type);
  assembleur_pression_->associer_domaine_dis_base(domaine_dis());
}

void Convection_Diffusion_Temperature_FT_Disc::completer()
{
  Convection_Diffusion_Temperature::completer();
  const Transport_Interfaces_FT_Disc& ft = ref_eq_interface_.valeur();
  const int ndist = ft.get_n_iterations_distance();
  if (stencil_width_ <= ndist)
    {
      Cerr << "The value of the stencil to compute mpoint should at least be strictly larger than "
           << "the width on which the distance is computed (given by n_iterations_distance in " << ft.que_suis_je() << ")" << finl;
      Cerr << "Current options : stencil= " << stencil_width_ << " and n_iterations_distance= " << ndist << " are invalid! Exiting." << finl;
      Process::exit();
    }
  if (divergence_free_velocity_extension_)
    {
      if (!solveur_pression_.non_nul())
        {
          Cerr << "You are trying to make the extension of velocity divergence free. " << finl;
          Cerr << "A poisson solver is then required and should be defined in your equation by the addition of the optional keyword solveur_pression_fictive " << finl;
          Cerr << "as in e.g.: solveur_pression_fictive  GCP { precond ssor { omega 1.5 } seuil 1e-15 impr }  " << finl;
          Process::exit();
        }
      const Navier_Stokes_std& ns = ref_eq_ns_.valeur();

      //
      // Build dummy Domaine_cl_dis object to pass to the assembleur_pression_ object:
      //
      // Step 0: Build string corresponding to the list of CLs we want to pass:
      SChaine instructions;
      instructions << "{" << finl;
      const Domaine& ladomaine=ns.domaine_dis().domaine();
      int nfront = ladomaine.nb_front_Cl();
      // The idea is to open the pressure on boundaries in contact with the other phase ie "(1-phase_)"
      // The goal is to let some fictitious pressure out (to accomodate for an (\int div(u_conv) dv !=0)
      for (int ifront=0; ifront<nfront; ifront++)
        {
          const Nom& nom_front = ladomaine.frontiere(ifront).le_nom();
          if (name_bc_opening_pressure_.contient_(nom_front))
            instructions << "    " <<nom_front << " sortie_libre_rho_variable champ_front_uniforme 1 0" << finl;
          else
            instructions << "    " <<nom_front << " symetrie" << finl;
        }
      instructions << "}" << finl;
      Cerr << "Interpretation de la chaine suivante:" << finl << instructions.get_str();
      EChaine is(instructions.get_str());

      // Step 1: discretise (see Equation_base::discretiser())
      Cerr << "Discretisation of fictitious CL ..." << finl;
      zcl_fictitious_.typer("Domaine_Cl_VDF");
      Domaine_Cl_VDF& zcl_fictitious_vdf = ref_cast(Domaine_Cl_VDF, zcl_fictitious_.valeur());
      Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis());
      zcl_fictitious_vdf.associer(domaine_vdf);
      zcl_fictitious_->associer_eqn(ns);
      // zcl_fictitious_->associer_inconnue(inconnue()); // Useless

      // Step 2: read input (see Equation_base::lire_cl())
      Cerr << "Interpreting input string ..." << finl;
      is >> zcl_fictitious_vdf ;

      // Step 3: completer (see Equation_base::completer())
      Cerr << "Completing fictitious CL ..." << finl;
      zcl_fictitious_->completer();

      // Associate to the newly created zcl :
      // required because It's not a good plan to use ns.domaine_Cl_dis() because of outlet_BC
      assembleur_pression_->associer_domaine_cl_dis_base(zcl_fictitious_.valeur());
      assembleur_pression_->completer(ns); // Should it be associated to (*this) or ns?
      //                                        I think it does not matter because Assembleur_base::completer is called and does nothing
      // On assemble la matrice de pression une seule et unique fois (puisqu'elle ne depend pas de rho...).
      assembleur_pression_->assembler(matrice_pression_);
      // Informe le solveur que la matrice a change :
      solveur_pression_->reinit();
    }
}

// Pour que milieu().mettre_a_jour(temps) ne plante pas...
// Et d'autres comme diffusivite_pour_transport...
Milieu_base& Convection_Diffusion_Temperature_FT_Disc::milieu()
{
  if (!fluide_dipha_.non_nul())
    {
      Cerr << "You forgot to associate the diphasic fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  // Cast non const cause acces const a fluide_phase!!
  return ref_cast_non_const(Milieu_base,fluide_dipha_->fluide_phase(phase_));
}
const Milieu_base& Convection_Diffusion_Temperature_FT_Disc::milieu() const
{
  if (!fluide_dipha_.non_nul())
    {
      Cerr << "You forgot to associate the diphasic fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return fluide_dipha_->fluide_phase(phase_);
}

void Convection_Diffusion_Temperature_FT_Disc::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (! sub_type(Fluide_Diphasique, un_milieu))
    {
      Cerr << "Erreur dans Convection_Diffusion_Temperature_FT_Disc::associer_milieu_base\n"
           << " On attendait un fluide diphasique" << finl;
      Cerr << "Error for Convection_Diffusion_Temperature_FT_Disc::associer_milieu_base\n"
           << "A Fluide_Diphasique medium was expected." << finl;
      exit();
    }
  fluide_dipha_ = ref_cast(Fluide_Diphasique, un_milieu);
}

/*! @brief Methode appelee par Transport_Interfaces_xxx::test_suppression_interfaces_sous_domaine() lorqu'une interfaces disparait.
 *
 * Il faut remettre la temperature de saturation dans
 *   l'inclusion supprimee.
 *
 */
void Convection_Diffusion_Temperature_FT_Disc::suppression_interfaces(const IntVect& num_compo,
                                                                      const ArrOfInt& flags_compo_a_supprimer,
                                                                      int nouvelle_phase)
{
  // Si la nouvelle phase n'est pas la phase resolue, ne rien faire
  if (nouvelle_phase != phase_)
    return;

  const int n = domaine_dis().domaine().nb_elem();
  assert(num_compo.size() == n);
  const int nb_valeurs_temporelles = inconnue()->nb_valeurs_temporelles();
  // Il faut traiter toutes les cases temporelles:
  //  selon l'ordre des equations dans le probleme, la roue a deja ete tournee ou pas...
  // Note B.M. est ce que c'est compatible avec la spec de ICOCO ? (modif
  //  du temps n autorisee ???)
  for (int t = 0; t < nb_valeurs_temporelles; t++)
    {
      DoubleTab& temp = inconnue()->futur(t);
      for (int i = 0; i < n; i++)
        {
          const int c = num_compo[i];
          if (c >= 0 && flags_compo_a_supprimer[c])
            temp[i] = TSAT_CONSTANTE;
        }
      temp.echange_espace_virtuel();
    }
}

int Convection_Diffusion_Temperature_FT_Disc::preparer_calcul()
{
  return Equation_base::preparer_calcul();
}

// A few methods for TCL only:
double Convection_Diffusion_Temperature_FT_Disc::get_flux_to_face(const int num_face) const
{
  double interfacial_flux = 0.;
  const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis();
  if (!sub_type(Domaine_Cl_VDF, zcldis))
    {
      Cerr << "Woops! Not VDF";
      Process::exit();
    }
  const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, zcldis);
  const Cond_lim& la_cl = zclvdf.la_cl_de_la_face(num_face);
  const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
  const Nom& bc_name = la_cl->frontiere_dis().le_nom();
  const int ndeb = le_bord.num_premiere_face();
//  Cerr <<  " BC: " << la_cl.valeur() << " name: " << bc_name << finl;
//  Cerr << "Dealing with face " << num_face <<  " belonging to BC " << la_cl.valeur();
  if ( sub_type(Neumann_paroi_adiabatique,la_cl.valeur()) )
    {
      Cerr << "paroi_adiabatique" << finl;
      return interfacial_flux;
    }
  else if ( sub_type(Neumann_paroi,la_cl.valeur()) )
    {
      Cerr << "paroi_flux_impose" << finl;
      const Neumann_paroi& la_cl_typee = ref_cast(Neumann_paroi, la_cl.valeur());
      double phi_imp = la_cl_typee.champ_front()->valeurs()(num_face-ndeb); // Should it be valeurs instead of ?
      return phi_imp;
    }
  else if (sub_type(Echange_impose_base, la_cl.valeur()))
    {
      Cerr << "paroi_temperature_imposee (among other possibilities)" << finl;
      /* Le terme de flux calcule a partir du couple(h_imp,T_ext) s'ecrit :
      //                           h_t(T_ext - T_entier)*Surf
      //                          avec h_t : coefficient d'echange global.
       * */
      const Echange_impose_base& la_cl_typee = ref_cast(Echange_impose_base, la_cl.valeur());
      //  Cerr <<  "   Face " << num_face << " is actually #" << num_face-ndeb << " on this boundary." << finl;
      const double h = la_cl_typee.h_imp(num_face-ndeb);
      const double T_imp = la_cl_typee.T_ext(num_face-ndeb);
      //  Cerr <<  "   Reading coefficients h= " << h << " and T_imp= " << T_imp << " for heat flux evaluation." << finl;
      // The flux is between the wall and Tsat :
      interfacial_flux = h*(T_imp - TSAT_CONSTANTE); // What about the area? surf should be the interfacial area or come from the face area?
      //                                                it is taken into account after this function.
      return interfacial_flux;
    }
  /*  else if ( sub_type(paroi_contact,la_cl.valeur()) )
      {
        Cerr << que_suis_je() << "::get_flux_to_face() thermal BC " << la_cl.valeur()
             << " at face " << num_face << " not supported yet." << finl;
        Cerr << "The BC type " << la_cl << " for boundary "<< la_cl->champ_front().le_nom()
             << " is not supported yet."<< finl;
        Process::exit();
      } */
  else
    {
      Cerr << que_suis_je() << "::get_flux_to_face(). " ;
      Cerr << "The BC type " << la_cl.valeur() << " for boundary "<< bc_name
           << " is not supported yet."<< finl;
      Process::exit();
    }
  return interfacial_flux;
}

double Convection_Diffusion_Temperature_FT_Disc::get_Twall_at_face(const int num_face) const
{
  double flux=0., Twall=0.;
  get_flux_and_Twall(num_face,
                     flux, Twall);
  return Twall;
}

double Convection_Diffusion_Temperature_FT_Disc::get_Twall_at_elem(const int elem) const
{
  //ArrOfInt num_faces;
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& faces_elem = domaine_vf.face_voisins();
  const int nb_faces_voisins = elem_faces.dimension(1);
  // Struggle to get the boundary face
  int num_face=-1;
  int i;
  for (i=0; i<nb_faces_voisins; i++)
    {
      num_face = elem_faces(elem,i);
      // If it's a boundary face, one of the neighbours doesnot exist so it has "-1".
      // We detect a boundary that way:
      const int elemb = faces_elem(num_face, 0) + faces_elem(num_face, 1) +1;
      if (elem == elemb)
        {
          //num_faces[idx] = num_face;
          break;
        }
    }
  if (i==nb_faces_voisins)
    {
      Cerr << "Error. No boundary face found in this element "<< elem << finl;
      Process::exit();
    }
  double flux=0., Twall=0.;
  get_flux_and_Twall(num_face,
                     flux, Twall);
  return Twall;
}

void Convection_Diffusion_Temperature_FT_Disc::get_flux_and_Twall(const int num_face,
                                                                  double& flux, double& Twall) const
{
  flux = 0.;
  const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis();
  if (!sub_type(Domaine_Cl_VDF, zcldis))
    {
      Cerr << "Woops! Not VDF";
      Process::exit();
    }
  const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, zcldis);
  const Cond_lim& la_cl = zclvdf.la_cl_de_la_face(num_face);
  const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
  const Nom& bc_name = la_cl->frontiere_dis().le_nom();
  const int ndeb = le_bord.num_premiere_face();
// Cerr <<  " BC: " << la_cl.valeur() << " name: " << bc_name << finl;
// Cerr << "Dealing with face " << num_face <<  " belonging to BC " << la_cl.valeur();
  if ( sub_type(Neumann_paroi_adiabatique,la_cl.valeur()) )
    {
      Cerr << "paroi_adiabatique" << finl;
      flux = 0.;
      //
      Cerr << "How can we set Twall when temperature(elem) is not valid? Or is it?" << finl;
      Process::exit();
    }
  else if ( sub_type(Neumann_paroi,la_cl.valeur()) )
    {
      Cerr << "paroi_flux_impose" << finl;
      const Neumann_paroi& la_cl_typee = ref_cast(Neumann_paroi, la_cl.valeur());
      flux = la_cl_typee.champ_front()->valeurs()(num_face-ndeb); // Should it be valeurs instead of ?
      //
      Cerr << "How can we set Twall when temperature(elem) is not valid? Or is it?" << finl;
      Process::exit();
    }
  else if (sub_type(Echange_impose_base, la_cl.valeur()))
    {

      if (sub_type(Echange_contact_VDF_FT_Disc, la_cl.valeur()))
        {
          // Cerr << "Echange_contact_VDF_FT_Disc" << finl;
          /* Le terme de flux calcule a partir du (phi_ext) s'ecrit :
          //                           phi_ext_*Surf
          //                           avec phi_ext : heat flux density at the surface.
           * */
          const Echange_contact_VDF_FT_Disc& la_cl_typee = ref_cast(Echange_contact_VDF_FT_Disc, la_cl.valeur());
          //  Cerr <<  "   Face " << num_face << " is actually #" << num_face-ndeb << " on this boundary." << finl;
          const double h = la_cl_typee.h_imp(num_face-ndeb);
          const double T_imp = la_cl_typee.Ti_wall(num_face-ndeb);
          // Cerr <<  "   Reading coefficients h= " << h << " and T_imp= " << T_imp << " for heat flux evaluation." << finl;
          // We keep the term + h*(T_imp - TSAT_CONSTANTE) as h = 0.;
          //
          flux = la_cl_typee.flux_exterieur_impose(num_face-ndeb) + h*(T_imp - TSAT_CONSTANTE);
          Twall = T_imp;
        }
      else
        {
          Cerr << "paroi_temperature_imposee (among other possibilities)" << finl;
          /* Le terme de flux calcule a partir du couple(h_imp,T_ext) s'ecrit :
          //                           h_t(T_ext - T_entier)*Surf
          //                          avec h_t : coefficient d'echange global.
           * */
          const Echange_impose_base& la_cl_typee = ref_cast(Echange_impose_base, la_cl.valeur());
          //  Cerr <<  "   Face " << num_face << " is actually #" << num_face-ndeb << " on this boundary." << finl;
          const double h = la_cl_typee.h_imp(num_face-ndeb);
          const double T_imp = la_cl_typee.T_ext(num_face-ndeb);

          // Cerr <<  "   Reading coefficients h= " << h << " and T_imp= " << T_imp << " for heat flux evaluation." << finl;
          // The flux is between the wall and Tsat :
          flux = h*(T_imp - TSAT_CONSTANTE); // What about the area? surf should be the interfacial area or come from the face area?
          //                                                it is taken into account after this function.
          Twall = T_imp;

        }
    }
  /*  else if ( sub_type(paroi_contact,la_cl.valeur()) )
      {
        Cerr << que_suis_je() << "::get_flux_to_face() thermal BC " << la_cl.valeur()
             << " at face " << num_face << " not supported yet." << finl;
        Cerr << "The BC type " << la_cl << " for boundary "<< la_cl->champ_front().le_nom()
             << " is not supported yet."<< finl;
        Process::exit();
      } */
  else
    {
      Cerr << que_suis_je() << "::get_flux_and_Twall(). " ;
      Cerr << "The BC type " << la_cl.valeur() << " for boundary "<< bc_name
           << " is not supported yet."<< finl;
      Process::exit();
    }
}


double Convection_Diffusion_Temperature_FT_Disc::get_Twall(const int num_face) const
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
  const IntTab& faces_elem = domaine_vf.face_voisins();
  // On of the neighbours doesnot exist so it has "-1". We get the other elem by:
  const int elem = faces_elem(num_face, 0) + faces_elem(num_face, 1) +1;
  const DoubleTab& temperature = inconnue()->valeurs();

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

  double d=0;
  for (int i=0; i<3; i++)
    d += (xyz_face[i] - P[i])*(xyz_face[i] - P[i]);
  d= sqrt(d);
  const double flux = get_flux_to_face(num_face);
  const double k = fluide_dipha_->fluide_phase(phase_).conductivite()->valeurs()(0,0);
//  Cerr << "flux/d/k" <<  flux << " " << d << " " << k <<  finl;
  // flux is incoming. So "-flux" is needed.
  const double Twall = temperature(elem) - d/k*flux;
// Cerr << "We have Twall = "<< Twall << " at face= " << num_face << " elem= " << elem << finl;
  return Twall;
}
