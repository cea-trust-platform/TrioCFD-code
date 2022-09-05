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
#include <Ref_Convection_Diffusion_Temperature_FT_Disc.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Zone_VF.h>
#include <Discretisation_base.h>
#include <Probleme_FT_Disc_gen.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Domaine.h>
#include <Param.h>
#include <Champ_Uniforme.h>

static const double TSAT_CONSTANTE = 0.;

Implemente_instanciable_sans_constructeur(Convection_Diffusion_Temperature_FT_Disc,"Convection_Diffusion_Temperature_FT_Disc",Convection_Diffusion_Temperature);

Implemente_ref(Convection_Diffusion_Temperature_FT_Disc);

Convection_Diffusion_Temperature_FT_Disc::Convection_Diffusion_Temperature_FT_Disc()
{
  phase_ = -1;
  correction_courbure_ordre_=0; // Par defaut, pas de prise en compte de la courbure pour corriger le champ etendu delta_vitesse
  stencil_width_ = 8; //GB : Valeur par defaut de stencil_width_
  temp_moy_ini_ = 0.; //GB : Valeur par defaut de la temperature moyenne initiale
  nom_sous_zone_ = "unknown_sous_zone"; //GB : Valeur par defaut de la sous-zone de moyenne
  maintien_temperature_ = false;
}

Sortie& Convection_Diffusion_Temperature_FT_Disc::printOn(Sortie& os) const
{
  return Convection_Diffusion_Temperature::printOn(os);
}
// Description:
//  cf Convection_Diffusion_std::readOn(Entree&).
Entree& Convection_Diffusion_Temperature_FT_Disc::readOn(Entree& is)
{
  // Ne pas faire assert(fluide non nul)
  Convection_Diffusion_std::readOn(is);
  solveur_masse->set_name_of_coefficient_temporel("rho_cp_comme_T");
  Nom num=inconnue().le_nom(); // On prevoir le cas d'equation de scalaires passifs
  num.suffix("temperature_thermique");
  Nom nom="Convection_chaleur";
  nom+=num;
  terme_convectif.set_fichier(nom);
  terme_convectif.set_description((Nom)"Convective heat transfer rate=Integral(-rho*cp*T*u*ndS) [W] if SI units used");
  nom="Diffusion_chaleur";
  nom+=num;
  terme_diffusif.set_fichier(nom);
  terme_diffusif.set_description((Nom)"Conduction heat transfer rate=Integral(lambda*grad(T)*ndS) [W] if SI units used");
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
  else if (mot=="diffusion")
    {
      if (je_suis_maitre())
        Cerr << " Equation_Concentration_FT::lire diffusivite" << finl;
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
        Cerr << " Equation_Concentration_FT::lire convection" << finl;
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
      is >> nom_sous_zone_;
      is >> temp_moy_ini_;
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
  return 1;
}

const Champ_base& Convection_Diffusion_Temperature_FT_Disc::vitesse_pour_transport() const
{
  return ref_eq_ns_.valeur().vitesse();
}

void Convection_Diffusion_Temperature_FT_Disc::preparer_pas_de_temps(void)
{
}

void Convection_Diffusion_Temperature_FT_Disc::corriger_pas_de_temps(double dt)
{
}

// A partir des valeurs du "champ" dans la phase "phase", calcule
// les valeurs du "champ" dans un voisinage d'epaisseur "stencil_width"
// en extrapolant lineairement et en supposant que la valeur a l'interface
// vaut "interfacial_value".
static void extrapoler_champ_elem(const Zone_VF&    zone_vf,
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
  const IntTab& elem_faces = zone_vf.elem_faces();
  const IntTab& faces_elem = zone_vf.face_voisins();
  const int nb_faces_elem = elem_faces.dimension(1);
  const int nb_elem       = elem_faces.dimension(0);

  const DoubleTab& centre_gravite_elem = zone_vf.xp();
  const DoubleTab& centre_gravite_face = zone_vf.xv();
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
          //       double dx = zone_vf.volumes(i_elem)/zone_vf.face_surfaces(elem_faces(i_elem,0));
          //       double dy = zone_vf.volumes(i_elem)/zone_vf.face_surfaces(elem_faces(i_elem,1));
          //       double dz = zone_vf.volumes(i_elem)/zone_vf.face_surfaces(elem_faces(i_elem,2));
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
  DoubleTab gradient_old;
  for (int iteration = 0; iteration < stencil_width; iteration++)
    {
      // Copie de la valeur du gradient: on ne veut pas utiliser les valeurs
      // calculees lors de l'iteration courante
      gradient_old = gradient;
      // La valeur sur un element est la moyenne des valeurs sur les elements voisins
      for (i_elem = 0; i_elem < nb_elem; i_elem++)
        {
          if (indicatrice[i_elem] != indic_phase)
            {
              // Ne pas toucher au gradient de la phase "phase".
              // Iterer sur les autres valeurs.
              double somme = 0.;
              double coeff = 0.;
              for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                {
                  const int face = elem_faces(i_elem, i_face);
                  const int voisin = faces_elem(face, 0) + faces_elem(face, 1) - i_elem;
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
              if (coeff > 0.)
                gradient[i_elem] = somme / coeff;
            }
        }
      gradient.echange_espace_virtuel();
    }
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
          //champ[i_elem] = d * grad * (1. - 0.5 * div_n * d + div_n * div_n * d * d / 6.) + interfacial_value;
        }
    }
  champ.echange_espace_virtuel();
}

// Description: met a jour le champ grad_t en fonction du champ inconnue.
//  Attention, l'inconnue est modifiee (on etend le champ de temperature dans la phase
//  opposee.
void Convection_Diffusion_Temperature_FT_Disc::calculer_grad_t()
{
  Transport_Interfaces_FT_Disc& eq_interface = ref_eq_interface_.valeur();

  const DoubleTab& indicatrice = eq_interface.get_update_indicatrice().valeurs();
  const DoubleTab& distance_interface = eq_interface.get_update_distance_interface().valeurs();
  const DoubleTab& normale_interface = eq_interface.get_update_normale_interface().valeurs();
  //GB : Augmenter la constante de l'epaisseur
  //GB : const int stencil_width = 8;
  const int stencil_width = stencil_width_;
  const double interfacial_value = TSAT_CONSTANTE;

  const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());

  const double temps = schema_temps().temps_courant();

  // Extrapolation lineaire de la temperature des mailles diphasiques a partir
  // de la temperature des mailles de la "phase".
  DoubleTab& temperature = inconnue().valeur().valeurs();

  Navier_Stokes_FT_Disc& eq_navier_stokes = ref_cast(Navier_Stokes_FT_Disc, ref_eq_ns_.valeur());
  const DoubleTab& div_n = eq_navier_stokes.calculer_div_normale_interface().valeurs();

  extrapoler_champ_elem(zone_vf, indicatrice, distance_interface, normale_interface, div_n,
                        phase_, stencil_width, interfacial_value,
                        temperature,
                        grad_t_.valeur().valeurs(),
                        temps);

}

void Convection_Diffusion_Temperature_FT_Disc::calculer_mpoint(Champ_base& mpoint)
{
  const double invalid_test = -1.e25;

  calculer_grad_t();

  mpoint.valeurs() = grad_t_.valeur().valeurs();

  DoubleTab& val = mpoint.valeurs();
  const int n = val.size();
  for (int i = 0; i < n; i++)
    if (val[i] < invalid_test)
      val[i] = 0;

  const double k = fluide_dipha_.valeur().fluide_phase(phase_).conductivite()(0,0);
  const double L = fluide_dipha_.valeur().chaleur_latente();
  // L est la chaleur latente de changement de phase pour passer de
  // la phase 0 a la phase 1.
  double f = k / L;
  // Si on est dans la phase 0 et que L > 0, on doit avoir mpoint positif pour
  //  gradient(T) scalaire n negatif.
  if (phase_ == 0)
    f = -f;

  mpoint.valeurs() *= f;

  // Pour deverminage
  // mpoint.valeurs() = 10.;

  // Pour deverminage : on impose une variation de mpoint lineaire en z
  //const Zone_VF & zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
  //for (int i = 0; i < n ; i++)
  //  mpoint.valeurs()[i] = 500./8957.*(1.+0.1*(zone_vf.xp()(i,3)-0.0005)/0.001);

  mpoint.valeurs().echange_espace_virtuel();
}

DoubleTab& Convection_Diffusion_Temperature_FT_Disc::derivee_en_temps_inco(DoubleTab& derivee)
{
  // Application des deux operateurs: convection et diffusion
  DoubleTab& temperature = inconnue().valeur().valeurs();

  // Extrapolation lineaire de la temperature des mailles diphasiques a partir
  // de la temperature des mailles de la "phase".
  //  static const Stat_Counter_Id count = statistiques().new_counter(1, "ghost_T", 0);
  //  statistiques().begin_count(count);
  calculer_grad_t();
  //  statistiques().end_count(count);

  derivee = 0.;

  //  static const Stat_Counter_Id count2 = statistiques().new_counter(1, "terme_diffusif_T", 0);
  //  statistiques().begin_count(count2);
  terme_diffusif.ajouter(temperature, derivee);
  //  statistiques().end_count(count2);

  // Calcul du champ de vitesse de convection dans la "phase_"
  Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, ref_eq_ns_.valeur());
  ns.calculer_delta_u_interface(vitesse_convection_, phase_, correction_courbure_ordre_ /* ordre de la correction en courbure */);
  const Champ_Inc& vitesse_ns = ns.inconnue();
  vitesse_convection_.valeurs() += vitesse_ns.valeurs();

  //  static const Stat_Counter_Id count3 = statistiques().new_counter(1, "convection_T", 0);
  //  statistiques().begin_count(count3);
  {
    const DoubleTab& rhoCp = get_champ("rho_cp_comme_T").valeurs();
    DoubleTab derivee_tmp(derivee);
    derivee_tmp = 0.;
    terme_convectif.ajouter(temperature, derivee_tmp);
    derivee_tmp *= rhoCp;
    derivee += derivee_tmp;
  }
  //  statistiques().end_count(count3);


  solveur_masse.appliquer(derivee);

  return derivee;
}

void Convection_Diffusion_Temperature_FT_Disc::mettre_a_jour (double temps)
{
  Convection_Diffusion_Temperature::mettre_a_jour(temps);

  // GB : Debut du maintien artificiel de la temperature.
  if (maintien_temperature_)
    {
      const Nom nom_sous_zone = nom_sous_zone_;
      const double temp_moy_ini = temp_moy_ini_;
      const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
      const Domaine& dom = zone_vf.zone().domaine();
      const Sous_Zone& ss_zone = dom.ss_zone(nom_sous_zone);

      Transport_Interfaces_FT_Disc& eq_interface = ref_eq_interface_.valeur();
      const DoubleTab& indicatrice = eq_interface.get_update_indicatrice().valeurs();
      DoubleTab& temperature = inconnue().valeur().valeurs();
      const DoubleVect& volume = ref_cast(Zone_VF, zone_dis().valeur()).volumes();
      const int nb_elem = zone_vf.zone().nb_elem();

      // Calcul de la moyenne sous zone, et du facteur : fac = temp_moy_ini/temp_moy_ss_zone
      double temp_moy_ss_zone = 0.;
      double vol_liq = 0.;
      for(int i=0; i<ss_zone.nb_elem_tot() /*methode d'acces au nombre d'elements*/; i++)
        {
          int index = ss_zone[i];
          // assert(index < zone_vf.zone().nb_elem());
          if (index < nb_elem)
            {
              temp_moy_ss_zone += temperature[index] * indicatrice[index] * volume[index];
              vol_liq += indicatrice[index] * volume[index];
            }
        }
      vol_liq = mp_sum(vol_liq);
      temp_moy_ss_zone = mp_sum(temp_moy_ss_zone);
      if (vol_liq == 0. || temp_moy_ss_zone == 0.)
        {
          // ne fien faire
        }
      else
        {
          temp_moy_ss_zone /= vol_liq ;
          double fac = temp_moy_ini / temp_moy_ss_zone;

          // DEBUG : Imprime la valeure du facteur dans le .err
          int num_proc = Process::me();
          if (num_proc == 0)
            {
              const double temps_sch = schema_temps().temps_courant();
              Cerr << "GB:fac t= " << temps_sch << " f= " << fac << finl;
              // ofstream file("fac-maintien-temp.txt");
              // file << "GB:fac t= " << temps << " f= " << fac << finl;
            }

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
  const Zone_dis_base& une_zone_dis = zone_dis().valeur();
  LIST(REF(Champ_base)) & champs_compris = liste_champs_compris_;
  const int nb_valeurs_temps = schema_temps().nb_valeurs_temporelles();

  Nom nom;

  nom = Nom("temperature_") + le_nom();
  dis.discretiser_champ("temperature", une_zone_dis, nom, "K", 1 /* composantes */, nb_valeurs_temps, temps, la_temperature);
  champs_compris.add(la_temperature.valeur());
  champs_compris_.ajoute_champ(la_temperature);

  nom = Nom("temperature_grad_") + le_nom();
  dis.discretiser_champ("temperature", une_zone_dis, nom, "K/m", 1 /* composantes */, temps, grad_t_);
  champs_compris.add(grad_t_.valeur());
  champs_compris_.ajoute_champ(grad_t_);

  nom = Nom("vitesse_conv_") + le_nom();
  dis.discretiser_champ("vitesse", une_zone_dis, nom, "m/s", -1 /* nb composantes par defaut */, 1 /* valeur temporelle */, temps, vitesse_convection_);
  champs_compris.add(vitesse_convection_.valeur());
  champs_compris_.ajoute_champ(vitesse_convection_);
  Equation_base::discretiser();
}


// Pour que milieu().mettre_a_jour(temps) ne plante pas...
Milieu_base& Convection_Diffusion_Temperature_FT_Disc::milieu()
{
  if (!fluide_dipha_.non_nul())
    {
      Cerr << "You forgot to associate the diphasic fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  // Cast non const cause acces const a fluide_phase!!
  return ref_cast_non_const(Milieu_base,fluide_dipha_.valeur().fluide_phase(phase_));
}
const Milieu_base& Convection_Diffusion_Temperature_FT_Disc::milieu() const
{
  if (!fluide_dipha_.non_nul())
    {
      Cerr << "You forgot to associate the diphasic fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return fluide_dipha_.valeur().fluide_phase(phase_);
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

// Description:
//  Methode appelee par Transport_Interfaces_xxx::test_suppression_interfaces_sous_zone()
//  lorqu'une interfaces disparait. Il faut remettre la temperature de saturation dans
//  l'inclusion supprimee.
void Convection_Diffusion_Temperature_FT_Disc::suppression_interfaces(const IntVect& num_compo,
                                                                      const ArrOfInt& flags_compo_a_supprimer,
                                                                      int nouvelle_phase)
{
  // Si la nouvelle phase n'est pas la phase resolue, ne rien faire
  if (nouvelle_phase != phase_)
    return;

  const int n = zone_dis().zone().nb_elem();
  assert(num_compo.size() == n);
  const int nb_valeurs_temporelles = inconnue().nb_valeurs_temporelles();
  // Il faut traiter toutes les cases temporelles:
  //  selon l'ordre des equations dans le probleme, la roue a deja ete tournee ou pas...
  // Note B.M. est ce que c'est compatible avec la spec de ICOCO ? (modif
  //  du temps n autorisee ???)
  for (int t = 0; t < nb_valeurs_temporelles; t++)
    {
      DoubleTab& temp = inconnue().futur(t);
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
