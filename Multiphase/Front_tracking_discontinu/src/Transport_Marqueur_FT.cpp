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
// File:        Transport_Marqueur_FT.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_Marqueur_FT.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Discretisation_base.h>
#include <TRUSTTab.h>
#include <Fluide_Incompressible.h>
#include <Champ_Uniforme.h>
#include <Fluide_Diphasique.h>
#include <Zone_VF.h>
#include <Source_Reaction_Particules.h>
#include <Zone.h>
#include <Source_Masse_Ajoutee.h>
#include <communications.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Transport_Marqueur_FT,"Transport_Marqueur_FT",Transport_Interfaces_FT_Disc);
Implemente_ref(Transport_Marqueur_FT);

Transport_Marqueur_FT::Transport_Marqueur_FT()
{
  t_debut_integr_ = -1.;
  t_debut_inject_ = 1.e32;
  t_derniere_inject_ = -1.;
  dt_inject_ = -1.;
  t_debut_transfo_ = 1.e32;
  t_derniere_transfo_ = -1;
  dt_transfo_ = -1.;
  nb_it_ = 1;
  implicite_ = 0;
  diametre_min_ = -1.;
  beta_transfo_ = -1.;
  methode_calcul_vp_ = INTERPOLEE;
  methode_couplage_ = SUIVI;
  phase_marquee_ = -1;
  contrib_one_way = 1;
  champs_compris_.ajoute_nom_compris("densite_particules");
  champs_compris_.ajoute_nom_compris("volume_particules");

}


/*! @brief
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_Marqueur_FT::printOn(Sortie& os) const
{
  return os;
}

/*! @brief Lecture d'une equation sur un flot d'entree.
 *
 * Voir Equation_base::readOn()
 *   Exception : le bloc de lecture des conditions limites doit etre vide
 *                 Conditions_limites { }
 *
 */
Entree& Transport_Marqueur_FT::readOn(Entree& is)
{
  Equation_base::readOn(is);

  if (t_debut_integr_==-1.)
    t_debut_integr_ = schema_temps().temps_init();

  if (!est_egal(dt_inject_,-1.))
    {
      if (est_egal(t_debut_inject_,1.e32))
        t_debut_inject_ = t_debut_integr_;
      t_derniere_inject_ = t_debut_inject_-2.*dt_inject_;
    }

  if (!est_egal(dt_transfo_,-1))
    {
      if (est_egal(t_debut_transfo_,1.e32))
        t_debut_transfo_ = t_debut_integr_;
      t_derniere_transfo_ = t_debut_transfo_-2.*dt_transfo_;
    }

  const Navier_Stokes_std& eq_ns = ref_cast(Navier_Stokes_std,probleme().equation(0));
  const Milieu_base& mil = eq_ns.milieu();
  if ((sub_type(Fluide_Diphasique,mil)) && (phase_marquee_==-1))
    phase_marquee_ = 0;

  int ok_impl = 0;
  if (implicite_==1)
    {
      for (auto& itr : les_sources)
        {
          if (sub_type(Source_Masse_Ajoutee,itr.valeur()))
            ok_impl = 1;
        }
      if (!ok_impl)
        {
          Cerr<<"Implicite option cannot be used if there is no added weight source term."<<finl;
          exit();
        }
    }

  return is;
}

//Lecture d options specifique a l equation
//-methode_transport     : type de la methode de transport                  : vitesse_interpolee ou vitesse_particules
//-methode_couplage      : type de la methode de couplage                   : suivi ou one_way_coupling ou two_way_coupling
//-phase_marquee         : numero de la phase marquee                             : numero de la phase marquee (cas diphasique)
//-injection                 : informations specifiques a une injection          : cf lire_injection(...)
//-transformation_bulles : informations specifiques a une transformation : cf lire_transformation(...)
//-contribution_one_way  : pour activer (1) ou desactiver (0) l effet du fluide sur les particules
//-nb_iterations         : nombre de sous pas de temps pour resoudre l equation de qdm des particules (1 par defaut)
//-implicite             : pour impliciter (1) ou non (0) le schema en temps
void Transport_Marqueur_FT::set_param(Param& param)
{
  Equation_base::set_param(param);
  param.ajouter_non_std("methode_transport",(this));
  param.ajouter_non_std("methode_couplage",(this));
  param.ajouter("phase_marquee",&phase_marquee_);
  param.ajouter_non_std("injection",(this));
  param.ajouter_non_std("transformation_bulles",(this));
  param.ajouter("contribution_one_way",&contrib_one_way);
  param.ajouter_condition("(value_of_contribution_one_way_eq_0)_or_(value_of_contribution_one_way_eq_1)","contribution_one_way must be set to 0 or 1");
  param.ajouter("nb_iterations",&nb_it_);
  param.ajouter("implicite",&implicite_);
  param.ajouter_condition("(value_of_implicite_eq_0)_or_(value_of_implicite_eq_1)","implicite must be set to 0 or 1");
}

int Transport_Marqueur_FT::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  Motcle motlu;
  if (mot=="methode_transport")
    {
      is>>motlu;
      if (motlu=="vitesse_interpolee")
        methode_calcul_vp_ = INTERPOLEE;
      else if (motlu=="vitesse_particules")
        methode_calcul_vp_ = BILAN_QDM;
      else
        {
          Cerr<<"The keyword "<<motlu<<" is not recognized as a calculation method "<<finl;
          Cerr<<"of the velocity particles."<<finl;
          exit();
        }
      return 1;
    }
  else if (mot=="methode_couplage")
    {
      is>>motlu;
      if (motlu=="suivi")
        methode_couplage_ = SUIVI;
      else if (motlu=="one_way_coupling")
        methode_couplage_ = ONE_WAY_COUPLING;
      else if (motlu=="two_way_coupling")
        {
          methode_couplage_ = TWO_WAY_COUPLING;
          Navier_Stokes_std& eq_ns = ref_cast(Navier_Stokes_std,probleme().equation(0));
          Nom disc = probleme().discretisation().que_suis_je();
          Source source_part;
          Sources& les_sources_ns = eq_ns.sources();
          Source& source_particule=les_sources_ns.add(source_part);
          Cerr << "Creation of the Force Particule term for the two_way_coupling case."<< finl;
          Nom type_so = "Two_Way_Coupling";
          type_so += "_";
          if (Motcle(disc)=="VEFPreP1b")
            disc ="VEF";
          type_so += disc;
          Cerr<<"type_so = "<<type_so<<finl;
          source_particule.typer_direct(type_so);
          source_particule->associer_eqn(eq_ns);
          Source_Reaction_Particules& so_reac = ref_cast(Source_Reaction_Particules,source_particule.valeur());
          so_reac.associer_eq_marqueur(*this);
        }
      else
        {
          Cerr<<"The keyword "<<motlu<<" is not recognized as a fluid/particles coupling method."<<finl;
          exit();
        }
      return 1;
    }
  else if (mot=="injection")
    {
      lire_injection(is);
      return 1;
    }
  else if (mot=="transformation_bulles")
    {
      lire_transformation(is);
      return 1;
    }
  else
    return Equation_base::lire_motcle_non_standard(mot,is);
}


/*! @brief Lecture des conditions initiales dans un flot d'entree.
 *
 * Le format de lecture est le suivant:
 *      {
 *       ensemble_points { }                lecture de la positions des particules                   cf Maillage_FT_Disc::readOn
 *         [ proprietes_particules { } ]   lecture des proprietes materielles des particules cf Prprietes_part_vol::readOn
 *         [ t_debut_integration value ]
 *      }
 *
 */
Entree& Transport_Marqueur_FT::lire_cond_init(Entree& is)
{
  Param param_init(que_suis_je());
  param_init.ajouter("ensemble_points",&maillage_interface());
  param_init.ajouter("proprietes_particules",&proprietes_particules());
  param_init.ajouter("t_debut_integration",&t_debut_integr_);
  param_init.lire_avec_accolades_depuis(is);
  return is;
}

//La methode est surchargee pour ne pas avoir a lire des conditions limites
//specifiques pour cette equation car le traitement des particules aux CL
//depend de celles imposees a l hydraulique
Entree& Transport_Marqueur_FT::lire_cl(Entree& is)
{
  Param param_cl(que_suis_je());
  param_cl.lire_avec_accolades_depuis(is);
  return is;
}

//Lecture des informations specifiques a l injection
//-ensemble_points       : lecture de la positions des particules a injecter cf Maillage_FT_Disc::readOn
//-proprietes_particules : lecture des proprietes des particules a injecter cf Prprietes_part_vol::readOn
//-t_debut_injection         : Instant de la premiere injection
//-dt_injection                 : Periodique d injection
Entree& Transport_Marqueur_FT::lire_injection(Entree& is)
{
  dt_inject_= probleme().schema_temps().pas_temps_min();
  Param param_inj(que_suis_je());
  param_inj.ajouter("ensemble_points",&maillage_inject());
  param_inj.ajouter("proprietes_particules",&proprietes_inject());
  param_inj.ajouter("t_debut_injection",&t_debut_inject_);
  param_inj.ajouter("dt_injection",&dt_inject_);
  param_inj.lire_avec_accolades_depuis(is);
  return is;
}

//Lecture des informations specifiques a la transformation
//-localisation       : nb_sz nom_sz1 ... (lecture des sous zones ou sont declenchees la transformation)
//-taille               : diametre_min value ou beta_transfo value (diamtre minimum ou coefficient multiplicatif du volume de controle)
//-t_debut_transfo    : Instant de la premiere injection
//-dt_transfo              : Periodique de transformation
//-interface               : Lecture du nom de l equation d interface separant les deux phases
Entree& Transport_Marqueur_FT::lire_transformation(Entree& is)
{
  const Navier_Stokes_std& eq_ns = ref_cast(Navier_Stokes_std,probleme().equation(0));
  const Milieu_base& mil = eq_ns.milieu();
  if (!sub_type(Fluide_Diphasique,mil))
    {
      Cerr<<"Transformation of bubbles into particles can be selected"<<finl;
      Cerr<<"only for two phases calculation."<<finl;
      exit();
    }

  dt_transfo_ = probleme().schema_temps().pas_temps_min();
  Param param_transfo(que_suis_je());
  param_transfo.ajouter("localisation",&nom_sz_transfo_);
  param_transfo.ajouter("diametre_min",&diametre_min_);
  param_transfo.ajouter("beta_transfo",&beta_transfo_);
  param_transfo.ajouter("t_debut_transfo",&t_debut_transfo_);
  param_transfo.ajouter("dt_transfo",&dt_transfo_);
  param_transfo.ajouter("interface",&nom_eq_interf_,Param::REQUIRED);
  param_transfo.lire_avec_accolades_depuis(is);
  return is;
}

//Discretisation du champ inconnu (champ_bidon_)
void Transport_Marqueur_FT::discretiser(void)
{
  const Probleme_base&  pb = get_probleme_base();
  const Discretisation_base& discr = pb.discretisation();
  //Discretisation d un champ inconnu bidon
  Nom fieldname;
  Nom suffix("_");
  suffix += le_nom();
  fieldname = "CHAMP_BIDON";
  fieldname += suffix;
  const Zone_dis_base& une_zone_dis = zone_dis().valeur();
  const double temps = schema_temps().temps_courant();
  const int nb_valeurs_temps = schema_temps().nb_valeurs_temporelles();

  discr.discretiser_champ("champ_elem", une_zone_dis,
                          fieldname, "-",
                          1 /* composantes */, nb_valeurs_temps,
                          temps,
                          champ_bidon_);
  champ_bidon_.associer_eqn(*this);

  Transport_Interfaces_FT_Disc::discretiser();

  const Zone& zone = une_zone_dis.zone();
  maillage_interface().associer_zone(zone);
}

/*! @brief Complete la construction (initialisation) des objets : maillage_interface et proprietes_particules_
 *
 *         maillage_inject_ et proprietes_inject_
 *
 */
void Transport_Marqueur_FT::completer()
{
  les_sources.completer();
  ////la_zone_Cl_dis->completer();
  const Zone& zone = la_zone_dis.valeur().zone();

  Maillage_FT_Disc& ens_points = maillage_interface();
  ens_points.associer_zone(zone);
  ens_points.generer_structure();
  ens_points.associer_equation_transport(*this);

  proprietes_particules().completer();
  const ArrOfInt& som_init =  ens_points.som_init_util();
  proprietes_particules().nettoyer(som_init);

  Maillage_FT_Disc& ens_points_inject = maillage_inject();
  ens_points_inject.associer_zone(zone);
  ens_points_inject.generer_structure();
  ens_points_inject.associer_equation_transport(*this);

  proprietes_inject().completer();
  const ArrOfInt&   som_init_inject =  ens_points_inject.som_init_util();
  proprietes_inject().nettoyer(som_init_inject);

  // Pas d'appel a Equation_base::completer donc:
  initialise_residu();
}


//Initialisation de delta_v et dimensionnement de source_stockage_
//dans le cas ou l on resoud le bila de qdm des particules
int Transport_Marqueur_FT::preparer_calcul()
{

  if (methode_calcul_vp_==BILAN_QDM)
    {
      int nb_part = proprietes_particules().nb_particules();
      int dim = Objet_U::dimension;
      update_delta_v(0,nb_part,maillage_interface());
      source_stockage_.resize(nb_part,dim);
    }

  return 1;
}

const Champ_Inc& Transport_Marqueur_FT::inconnue(void) const
{
  return champ_bidon_;
}

Champ_Inc& Transport_Marqueur_FT::inconnue(void)
{
  return champ_bidon_;
}


//Preparation du pas de temps pour l equataion
//-injection et transformation
//-dans le cas ou l on resoud le bilan de qdm des particules
// calcul des proprietes du fluide aux positions des particules
// calcul des sources intrevenant dans le bialn de qdm
// actualisation de delta_v
bool Transport_Marqueur_FT::initTimeStep(double dt)
{
  Equation_base::initTimeStep(dt);

  double temps = probleme().schema_temps().temps_courant();

  injecter(temps);
  transformer(temps);

  if ((sup_ou_egal(temps+dt,t_debut_integr_)) || ((methode_couplage_==TWO_WAY_COUPLING) && (!contrib_one_way)))
    {
      if (methode_calcul_vp_==BILAN_QDM)
        {
          const int dim = Objet_U::dimension;
          const int nb_particules = proprietes_particules().nb_particules();
          const Maillage_FT_Disc& ens_points = maillage_interface();
          DoubleTab& vitesse_p = proprietes_particules().vitesse_particules();
          dt_p_ = dt/nb_it_;
          calculer_proprietes_fluide_pos_particules(ens_points);
          imposer_cond_lim();
          source_stockage_.resize(nb_particules,dim);

          for (int i=0; i<nb_it_; i++)
            {
              source_stockage_=0.;
              les_sources.ajouter(source_stockage_);
              //Stockage de delta_v pour le prochain pas de temps
              update_delta_v(0,nb_particules,ens_points,0);
              resoudre_edo(vitesse_p,source_stockage_,dt_p_);
            }
        }
    }

  return true;
}

//-Transport des particules
//-Nettoyage : suppression des particules virtuelles et de celles sortant du domaine
void Transport_Marqueur_FT::mettre_a_jour(double temps)
{
  integrer_ensemble_lagrange(temps);
  maillage_interface().nettoyer_noeuds_virtuels_et_frontieres();
}

/*! @brief
 *
 */
int Transport_Marqueur_FT::sauvegarder(Sortie& os) const
{
  return 0;
}
int Transport_Marqueur_FT::reprendre(Entree& is)
{
  return 1;
}

//Integration des trajectoires de particules
// La methode est une version simplifiee de deplacer_maillage_ft_v_fluide(...)
// -Interpolation du champ de vitesse de l equation de qdm ou resolution du bilan de qdm des particules
// -Calcul du deplacement a appliquer aux particules
// -Transport des particules en fonction du deplacement calcule
// -Mise a jour du temps de l objet ensemble de points meme si t<t_debut_integr
// Rq : t_debut_integr est fixe par defaut a t_init mais peut etre specifie par l utilisateur

void Transport_Marqueur_FT::integrer_ensemble_lagrange(const double temps)
{
  const double t_debut_integr = temps_debut_integration();
  Maillage_FT_Disc& ens_points = maillage_interface();
  if (sup_ou_egal(temps,t_debut_integr))
    {
      DoubleTab&   deplacement = deplacement_som();
      calcul_vitesse_p(deplacement);

      const double delta_t = probleme_base_->schema_temps().pas_de_temps();
      // Calcul du deplacement :
      deplacement *= delta_t;

      //preparation des tableaux contenant les proprietes volumiques
      Proprietes_part_vol& proprietes = proprietes_particules();
      preparer_tableaux_avant_transport(ens_points,proprietes);
      ens_points.transporter_simple(deplacement);
      //update des tableaux contenant les proprietes volumiques
      update_tableaux_apres_transport(ens_points,proprietes);

    } //Fin if t_debut_integr

  //Rq : La modification du temps devrait etre faite dans la mise a jour de l ensemble de points ...
  ens_points.changer_temps(temps);
}

//Evaluation des proprietes du fluide aux positions des particules
//vit_fluide_som_, grad_P_som_, rho_fluide_som_ et visco_dyn_fluide_som_ sont interpoles
void Transport_Marqueur_FT::calculer_proprietes_fluide_pos_particules(const Maillage_FT_Disc& ens_points)
{
  const Navier_Stokes_std& eq_ns = ref_cast(Navier_Stokes_std,probleme().equation(0));
  const Champ_base& champ_vitesse =  eq_ns.inconnue();
  const Champ_base& champ_pression =  eq_ns.pression();

  const Solveur_Masse& le_solveur_masse_ns = eq_ns.solv_masse();
  const Operateur_Grad& gradient = eq_ns.operateur_gradient();

  Champ_Inc champ_grad_p(eq_ns.grad_P());
  gradient.calculer(champ_pression.valeurs(), champ_grad_p.valeurs());
  le_solveur_masse_ns.appliquer(champ_grad_p.valeurs());
  champ_grad_p.valeurs().echange_espace_virtuel();
  DoubleTabFT tab_som;
  const DoubleTab& pos = ens_points.sommets();
  const ArrOfInt& elem = ens_points.sommet_elem();
  const int nb_pos_tot = pos.dimension(0);
  const int dim = Objet_U::dimension;

  int nb_positions;
  nb_positions=0;
  for (int i = 0; i < nb_pos_tot; i++)
    {
      const int num_elem = elem[i];
      if (num_elem >= 0)
        {
          nb_positions++;
        }
    }
  vit_fluide_som_.resize(nb_positions, dim);
  grad_P_som_.resize(nb_positions,dim);
  rho_fluide_som_.resize(nb_positions);
  visco_dyn_fluide_som_.resize(nb_positions);

  calculer_vitesse_transport_interpolee(champ_vitesse,ens_points,tab_som,1);
  nb_positions =0;
  for (int i = 0; i < nb_pos_tot; i++)
    {
      if (elem[i] >= 0)
        {
          for (int j = 0; j < dim; j++)
            vit_fluide_som_(nb_positions, j) = tab_som(i, j);
          nb_positions++;
        }
    }

  calculer_vitesse_transport_interpolee(champ_grad_p,ens_points,tab_som,1);
  nb_positions =0;
  for (int i = 0; i < nb_pos_tot; i++)
    {
      if (elem[i] >= 0)
        {
          for (int j = 0; j < dim; j++)
            grad_P_som_(nb_positions, j) = tab_som(i, j);
          nb_positions++;
        }
    }


  if (sub_type(Navier_Stokes_FT_Disc,eq_ns))
    {
      const Milieu_base& mil = eq_ns.milieu();
      if (sub_type(Fluide_Diphasique,mil))
        {
          const Fluide_Diphasique& fluide_diph = ref_cast(Fluide_Diphasique,mil);
          const Fluide_Incompressible& fluide = fluide_diph.fluide_phase(phase_marquee_);
          const DoubleTab& rho = fluide.masse_volumique().valeurs();
          const DoubleTab& visco_dyn = fluide.viscosite_dynamique().valeurs();
          rho_fluide_som_ = rho(0,0);
          visco_dyn_fluide_som_ = visco_dyn(0,0);
        }
      else
        {
          const Fluide_base& fluide = ref_cast(Fluide_base,mil);
          const DoubleTab& rho = fluide.masse_volumique().valeurs();
          const DoubleTab& visco_dyn = fluide.viscosite_dynamique().valeurs();
          rho_fluide_som_ = rho(0,0);
          visco_dyn_fluide_som_ = visco_dyn(0,0);
        }
      const Navier_Stokes_FT_Disc& eq_ns_FT = ref_cast(Navier_Stokes_FT_Disc,eq_ns);
      if (!eq_ns_FT.is_terme_gravite_rhog()
          && eq_ns_FT.milieu().a_gravite())
        {
          const DoubleTab& gravite = eq_ns_FT.milieu().gravite().valeurs();
          if (gravite.nb_dim() != 2 || gravite.dimension(1) != dimension)
            {
              Cerr << "Error for method calculer_proprietes_fluide_pos_particules\n";
              Process::exit();
            }

          for (int i=0; i<grad_P_som_.dimension(0); i++)
            for (int j=0; j<grad_P_som_.dimension(1); j++)
              grad_P_som_(i,j) += rho_fluide_som_(i)*gravite(0,j);

        }
    }
  else
    {
      const Fluide_base& fluide = ref_cast(Fluide_base,eq_ns.milieu());
      const Champ_base& champ_masse_vol =  fluide.masse_volumique().valeur();
      const Champ_base& champ_visco_dyn =  fluide.viscosite_dynamique().valeur();

      DoubleTab& les_positions = tableaux_positions();
      IntVect& les_elements = vecteur_elements();
      les_positions.resize(nb_pos_tot, dim);
      les_elements.resize(nb_pos_tot);
      //int i, j;
      nb_positions = 0;
      for (int i = 0; i < nb_pos_tot; i++)
        {
          const int num_elem = elem[i];
          if (num_elem >= 0)
            {
              for (int j = 0; j < dim; j++)
                les_positions(nb_positions, j) = pos(i, j);
              les_elements(nb_positions) = num_elem;
              nb_positions++;
            }
        }


      if (sub_type(Champ_Uniforme,champ_masse_vol))
        {
          double rho = champ_masse_vol.valeurs()(0,0);
          rho_fluide_som_ = rho;
          const DoubleTab& gravite = milieu().gravite().valeurs();

          //Dans le cas N_S et Fluide_Inc gradient_P contient grad(P/rho+gz) avec rho uniforme
          //Par consequent on multiplie gradient_P par rho et on retranche rho*g
          for (int i=0; i<grad_P_som_.dimension(0); i++)
            for (int j=0; j<grad_P_som_.dimension(1); j++)
              {
                grad_P_som_(i,j) *= rho;
                grad_P_som_(i,j) -= rho*gravite(0,j);
              }
        }
      else
        calculer_scalaire_interpole(champ_masse_vol,ens_points,rho_fluide_som_,1);

      if (sub_type(Champ_Uniforme,champ_visco_dyn))
        visco_dyn_fluide_som_ = champ_visco_dyn.valeurs()(0,0);
      else
        champ_visco_dyn.valeur_aux_elems(les_positions,les_elements,visco_dyn_fluide_som_);
    }

}

//Appel a la methode d injection de particules si instant d injection
void Transport_Marqueur_FT::injecter(double temps)
{
  if (linject(temps))
    {
      const Maillage_FT_Disc& maill_inject = maillage_inject();
      const Proprietes_part_vol& propri_inject = proprietes_inject();
      injection(maill_inject,propri_inject);
      t_derniere_inject_ = temps;
    }
}

//Injection de particules
//-ajout d un esnemble de positions (et d un ensemble de proprietes si particules materielles)
//-actualisation de delta_v pour la partie injectee
void Transport_Marqueur_FT::injection(const Maillage_FT_Disc& maill_inject, const Proprietes_part_vol& propr_inject)
{
  int nb_part_deb = proprietes_particules().nb_particules();
  ajouter_points(maill_inject);
  proprietes_particules().ajouter_proprietes(propr_inject);
  int nb_part_fin = proprietes_particules().nb_particules();
  update_delta_v(nb_part_deb,nb_part_fin,maill_inject,1);

}

//Suppression de particules entrant dans la phase non marquee
//Appel a la methode de transformation si instant de transformation
void Transport_Marqueur_FT::transformer(double temps)
{
  transformer_particules();

  if (ltransfo(temps))
    {
      Maillage_FT_Disc maillage_transfo;
      Proprietes_part_vol proprietes_transfo;
      transformation(maillage_transfo,proprietes_transfo);
      injection(maillage_transfo,proprietes_transfo);
      t_derniere_transfo_ = temps;
    }


}

//-Identification des groupes connexes
//-calcul des proprietes geometriques equivalentes (positions, volumes et diametres) de ces goupes
//-Detection des groupes a supprimer (situe dans une sous zone specifie par utilisateur
//                                      ou de diametre inferieur a un diametre specifie par l utilsateur)
//-Construction de l objet maillage (ensemble de points) et proprietes (proprietes materielles) a injecter
//-Suppression des groupes connexes identifies a supprimer

void Transport_Marqueur_FT::transformation(Maillage_FT_Disc& marqueurs,Proprietes_part_vol& proprietes)
{
  //Appel de la methode permettant de detecter les groupes connexes
  //et de deteminer les proprietes des particules resultant de la transformation
  const Equation_base& eq = probleme_base_->get_equation_by_name(nom_eq_interf_);
  Transport_Interfaces_FT_Disc& eq_interf = ref_cast_non_const(Transport_Interfaces_FT_Disc,eq);
  const Equation_base& eqn_hydraulique = probleme_base_->equation(0);
  const Champ_base& champ_vitesse = eqn_hydraulique.inconnue().valeur();
  DoubleTab& vitesse_p = proprietes.vitesse_particules();

  IntVect num_compo;
  const int nb_compo = eq_interf.calculer_composantes_connexes_pour_suppression(num_compo);

  ArrOfDouble volumes;
  DoubleTab positions;
  const DoubleTab& indic = eq_interf.get_update_indicatrice().valeurs();
  calcul_proprietes_geometriques(num_compo, nb_compo, indic, volumes, positions);
  ArrOfInt flags_compos_a_supprimer;
  detection_groupes_a_supprimer(volumes, positions, flags_compos_a_supprimer);
  construction_ensemble_proprietes(num_compo, nb_compo, marqueurs, proprietes, flags_compos_a_supprimer, positions, volumes);
  eq_interf.calculer_vitesse_transport_interpolee(champ_vitesse,marqueurs,vitesse_p,1);
  eq_interf.suppression_interfaces(num_compo, flags_compos_a_supprimer);
}

//Suppression des particules situees dans la phase non marquee
void Transport_Marqueur_FT::transformer_particules()
{
  if ((phase_marquee_!=-1) && (nom_eq_interf_!="??"))
    {
      const Milieu_base& mil = probleme().equation(0).milieu();
      if (!sub_type(Fluide_Diphasique,mil))
        {
          Cerr<<"Particles transformation can be activated only for a Fluide_Diphasique medium."<<finl;
          exit();
        }
      Maillage_FT_Disc& ens_points = maillage_interface();
      ens_points.nettoyer_phase(nom_eq_interf_,phase_marquee_);
    }
}

//calcul des proprietes geometriques equivalentes aux groupes connexes (positions, volumes, diametres)
//par integration de la quantite (volume ou position) pondere par l indicatrice
//Construction du tableau elems_pos
//-elems_pos(el,0) contient l element dans lequel est localise le centre de gravite G
//-elems_pos(el,1) contient un element contenu dans le groupe connexe de centre G
//car le centre de gravite du groupe connexe peut ne pas etre contenu dans le groupe

void Transport_Marqueur_FT::calcul_proprietes_geometriques(const IntVect&        num_compo,
                                                           const int         nb_compo,
                                                           const DoubleTab&        indic,
                                                           ArrOfDouble&         volumes,
                                                           DoubleTab&                 positions)
{
  const Zone_dis_base& zdis_base = zone_dis().valeur();
  const Zone_VF& zvf = ref_cast(Zone_VF,zdis_base);
  const DoubleVect& vol_elem = zvf.volumes();
  const DoubleTab& xp_elem = zvf.xp();
  int nb_elem = zvf.nb_elem();
  int dim = xp_elem.dimension(1);

  volumes.resize_array(nb_compo, Array_base::NOCOPY_NOINIT);
  positions.resize(nb_compo, dim, Array_base::NOCOPY_NOINIT);
  volumes = 0.;
  positions = 0.;

  //Calcul du volume - du diametre - du centre de gravite pour chaque groupe connexe
  for (int elem=0; elem<nb_elem; elem++)
    {
      const int num_group = num_compo(elem);
      if (num_group>=0)
        {
          const double indic_phase = phase_marquee_+(1.-2*phase_marquee_)*indic(elem);
          const double v = indic_phase * vol_elem(elem);
          volumes[num_group] += v;
          for (int j=0; j<dim; j++)
            positions(num_group,j) += xp_elem(elem,j) * v;
        }
    }

  mp_sum_for_each_item(volumes);
  mp_sum_for_each_item(positions);

  for (int num_group=0; num_group<nb_compo; num_group++)
    {
      const double inv_v = 1. / volumes[num_group];
      for (int j=0; j<dim; j++)
        positions(num_group,j) *= inv_v;
    }
}

//Parmi l ensemble des groupes connexes identifies, on detecte ceux qui sont a supprimer
//Les criteres de suppression retenus sont :
// -le centre de gravite du groupe est dans une sous zone definie par l utilisateur
// (pour un groupe connexe donne, on parcours les sous zones pour tester si l une d entre elles
//  contient le centre de gravite du groupe)
// -le diametre caracteristique du groupe connexe considere est plus petit qu un diametre_min
// defini par l utilisateur ou plus petite que beta fois le volume elementaire contenant le
// centre de gravite du groupe (beta est specifie par l utilisateur)

void Transport_Marqueur_FT::detection_groupes_a_supprimer(const ArrOfDouble& volumes,
                                                          const DoubleTab&    positions,
                                                          ArrOfInt&               flags_compo_a_supprimer)
{
  const Zone_VF& zvf = ref_cast(Zone_VF, zone_dis().valeur());
  const Zone& zone_geom = zvf.zone();
  const DoubleVect& volume_elem = zvf.volumes();

  // Recherche de l'element contenant le centre de gravite de chaque composante connexe
  const int nb_compo = volumes.size_array();
  ArrOfInt elems;

  zone_geom.chercher_elements(positions, elems);

  flags_compo_a_supprimer.resize_array(nb_compo, Array_base::NOCOPY_NOINIT);
  flags_compo_a_supprimer = 0;

  // Construction d'un tableau de flags des elements reels contenus dans les sous-zones de destruction
  const int nb_elem = zone_geom.nb_elem();
  ArrOfBit flag_sous_zone(nb_elem);
  flag_sous_zone = 0;
  {
    const int nb_sous_zones = nom_sz_transfo_.size();
    for (int i_sous_zone = 0; i_sous_zone < nb_sous_zones; i_sous_zone++)
      {
        const Sous_Zone& sous_zone = zone_geom.ss_zone(nom_sz_transfo_[i_sous_zone]);
        const int nb_elem_sous_zone = sous_zone.nb_elem_tot();
        for (int i = 0; i < nb_elem_sous_zone; i++)
          {
            const int ielem = sous_zone[i];
            // La sous_zone contient des indices d'elements virtuels
            if (ielem < nb_elem)
              flag_sous_zone.setbit(sous_zone[i]);
          }
      }
  }

  // Volume max d'une inclusion a supprimer (conversion du diametre min en volume)
  const double volume_max = 4. / 3. * M_PI * pow(diametre_min_ * 0.5, 3); // 4/3*PI*R^3

  for (int i = 0; i < nb_compo; i++)
    {
      // Test sur le centre de gravite: suppression des interfaces dont le centre de gravite
      // est dans la sous-zone de suppression.
      // Attention, seul le processeur qui a le centre de gravite peut mettre le flag a 1.
      // On synchronise les flags ensuite.
      const int num_elem = elems[i];
      if (num_elem >= 0 && num_elem < nb_elem)
        {
          // L'element est chez moi. Est-il dans la sous_zone ?
          if (flag_sous_zone[num_elem])
            {
              // Test sur le volume:
              //  (facteur 1/6 pour compatibilite avec codage precedent.
              //   Il n'y a pas de raison particuliere, c'est comme ca)
              const double volume_local_max = beta_transfo_ / 6. * volume_elem[num_elem];
              const double v = volumes[i];
              // Je remplace inf_strict par <. De toutes facons on compare maintenant les volumes pas les diametres
              // cela change le epsilon. De plus, un test "inf_strict" n'a pas d'interet ici.
              if (v < volume_max || v < volume_local_max)
                flags_compo_a_supprimer[i] = 1;
            }
        }
    }
  // Synchronisation des flags entre processeurs
  mp_max_for_each_item(flags_compo_a_supprimer);
}

//Construction de la structure maillage (ensemble de points) et proprietes (proprietes materielles) a injecter
// en remplacement des groupes connexes
//Pour la structure maillage :   ens_points.associer_zone(...)
//                                   ens_points.generer_structure(...)
//                                   ens_points.associer_equation_transport(...)

//Pour la structure proprietes : dimensionnement remplissage des differentes proprietes
//                                 proprietes.fixer_nb_particules(...)
//                                 proprietes.completer(...)
//                                 proprietes.nettoyer(...)

void Transport_Marqueur_FT::construction_ensemble_proprietes(const IntVect&                num_compo,
                                                             const int                nb_compo,
                                                             Maillage_FT_Disc&                ens_points,
                                                             Proprietes_part_vol&        propri,
                                                             const ArrOfInt&                flags_compos_a_supprimer,
                                                             const DoubleTab&                positions,
                                                             const ArrOfDouble&                volumes)
{
  int phase_transfo = (phase_marquee_==0?1:0);
  const Fluide_Diphasique& fluide_diph = ref_cast(Fluide_Diphasique,probleme().equation(0).milieu());
  const DoubleTab& rho = fluide_diph.fluide_phase(phase_transfo).masse_volumique().valeurs();
  double rho_val = rho(0,0);

  const int dim = positions.dimension(1);

  int size_new = 0;
  for (int i = 0; i < nb_compo; i++)
    if (flags_compos_a_supprimer[i])
      size_new++;

  //remplir sommets_lu et proprietes
  const Zone& zone = la_zone_dis->zone();
  ens_points.associer_zone(zone);

  DoubleTab&   soms_tmp =  ens_points.sommets_lu();
  DoubleTab& vitesse_tmp = propri.vitesse_particules();
  DoubleTab& temp_tmp = propri.temperature_particules();
  DoubleTab& dv_tmp = propri.delta_v();
  DoubleTab& mvol_tmp =  propri.masse_vol_particules();
  DoubleTab& vol_tmp = propri.volume_particules();
  DoubleTab& diam_tmp = propri.diametre_particules();

  soms_tmp.resize(size_new,dim);
  vitesse_tmp.resize(size_new,dim);
  dv_tmp.resize(size_new,dim);
  temp_tmp.resize(size_new,1);
  mvol_tmp.resize(size_new,1);
  vol_tmp.resize(size_new,1);
  diam_tmp.resize(size_new,1);

  int compt=0;
  for (int i = 0; i<nb_compo; i++)
    {
      if (flags_compos_a_supprimer[i]>0)
        {
          for (int j=0; j<dim; j++)
            {
              soms_tmp(compt,j) =  positions(i,j);
              vitesse_tmp(compt,j) = 0.;
              dv_tmp(compt,j) = 0.;
            }
          temp_tmp(compt,0) = 273.;
          mvol_tmp(compt,0) = rho_val;
          vol_tmp(compt,0) = volumes[i];
          diam_tmp(compt,0) = pow(volumes[i]*3./(4.*M_PI), 1./3) * 2.;
          compt++;
        }
    }
  propri.fixer_nb_particules(size_new);

  ens_points.generer_structure();
  ens_points.associer_equation_transport(*this);

  propri.completer();
  const ArrOfInt& som_init =  ens_points.som_init_util();
  propri.nettoyer(som_init);

}

//Calcul de la vitesse des particules vitesse_p
//-methode_calcul_vp_==INTERPOLEE : interpolation du champ de vitesse du fluide
//-methode_calcul_vp_==BILAN_QDM  : resolution du bilan de dqm des particules
//Calcul du deplacement = vitesse_p*dt

void Transport_Marqueur_FT::calcul_vitesse_p(DoubleTab& deplacement) const
{
  if (methode_calcul_vp_==INTERPOLEE)
    {
      const Champ_base *u0_ptr = 0;

      //Pour une equation de transport d interface solide, la REF refequation_vitesse_transport est vide
      //On utilise directement l equation de N_S et non pas efequation_vitesse_transport.valeur()
      const Equation_base& eqn_hydraulique = probleme_base_->equation(0);
      const Champ_base& champ_vitesse = eqn_hydraulique.inconnue().valeur();
      calculer_vitesse_transport_interpolee(champ_vitesse,
                                            maillage_interface(),
                                            deplacement, 1 /* recalculer le champ de vitesse L2 */);

      if (sub_type(Navier_Stokes_FT_Disc, eqn_hydraulique))
        {
          // On recupere le saut de vitesse a l'interface (changement de phase)
          const Navier_Stokes_FT_Disc& ns = ref_cast(Navier_Stokes_FT_Disc, eqn_hydraulique);
          u0_ptr = ns.get_delta_vitesse_interface();
          if (u0_ptr)
            {
              const Champ_base& u0 = *u0_ptr;
              DoubleTabFT d2(deplacement);
              calculer_vitesse_transport_interpolee(u0,
                                                    maillage_interface(),
                                                    d2, 1 /* recalculer le champ de vitesse L2 */);
              const int n = d2.dimension(0);
              const int dim = d2.dimension(1);

              for (int i = 0; i < n; i++)
                {
                  for (int j = 0; j < dim; j++)
                    {
                      const double depl2 = d2(i,j);
                      deplacement(i,j) -= depl2;
                    }
                }
            }
        }
    }
  else if (methode_calcul_vp_==BILAN_QDM)
    {
      const DoubleTab& vitesse_p = proprietes_particules().vitesse_particules();
      deplacement = vitesse_p;
    }
}

//Actualisation de delta_v (vitesse_p-vit_fluide_som_) pour les particules allant de n_deb a n_fin-1
//Les proprietes du fluide sont re-estimees sur ens_points si calc=1

void Transport_Marqueur_FT::update_delta_v(int n_deb,int n_fin,const Maillage_FT_Disc& ens_points,int calc)
{
  assert(n_fin-n_deb==ens_points.sommets().dimension(0));
  DoubleTab& delta_vit = proprietes_particules().delta_v();
  const DoubleTab& vitesse_p = proprietes_particules().vitesse_particules();
  int dim = Objet_U::dimension;
  if (mp_sum(n_fin-n_deb))
    {
      if (calc)
        calculer_proprietes_fluide_pos_particules(ens_points);

      if (!implicite_)
        {
          for (int i=n_deb; i<n_fin; i++)
            for (int j=0; j<dim; j++)
              delta_vit(i,j) = vit_fluide_som_(i-n_deb,j)-vitesse_p(i,j);
        }
      else
        {
          for (int i=n_deb; i<n_fin; i++)
            for (int j=0; j<dim; j++)
              delta_vit(i,j) = vit_fluide_som_(i-n_deb,j);
        }
    }

}

//Methode de resolution du bilan de qdm des particules
//vitesse_p = vitesse_p + sources*dt
//distinction pour explicite et implicite
void Transport_Marqueur_FT::resoudre_edo(DoubleTab& vitesse_p, DoubleTab& une_source_stockage,const double delta_t)
{

  const int nb_particules = proprietes_particules().nb_particules();
  const int dim = Objet_U::dimension;
  const DoubleTab& masse_vol_p = proprietes_particules().masse_vol_particules();
  const DoubleTab& volume_p = proprietes_particules().volume_particules();

  if (!implicite_)
    {
      for (int i=0; i<nb_particules; i++)
        for (int j=0; j<dim; j++)
          vitesse_p(i,j) += (1./(masse_vol_p(i,0)*volume_p(i,0)))*contrib_one_way*source_stockage_(i,j)*delta_t;
    }
  else
    {
      for (int i=0; i<nb_particules; i++)
        for (int j=0; j<dim; j++)
          vitesse_p(i,j) += (1./((masse_vol_p(i,0)+0.5*rho_fluide_som_(i))*volume_p(i,0)))*contrib_one_way*source_stockage_(i,j)*delta_t;

      //Evaluation du terme source en tenant compte de la partice imlplicite de la masse ajoutee
      for (int i=0; i<nb_particules; i++)
        for (int j=0; j<dim; j++)
          source_stockage_(i,j) -= (0.5*rho_fluide_som_(i)/(masse_vol_p(i,0)+0.5*rho_fluide_som_(i)))*source_stockage_(i,j);
    }
}

//Modifie le champ de vitesse pour les particules situees sur un bord qui porte :
//une condition limite de paroi_fixe ou de symetrie pour l hydraulique
//La modification consiste en une reflexion du champ de vitesse (condition de rebond)
//Appel a appliquer_reflexion_vitesse(...)

void Transport_Marqueur_FT::imposer_cond_lim()
{
  const Zone_VF& zone_vf = ref_cast(Zone_VF,zone_dis().valeur());
  DoubleTab& vitesse_p =  proprietes_particules().vitesse_particules();
  const Maillage_FT_Disc& maillage = maillage_interface();
  const DoubleTab& pos = maillage.sommets();
  const ArrOfInt& elem = maillage.sommet_elem();
  const ArrOfInt& som_face = maillage.sommet_face_bord();
  const int nb_pos_tot = pos.dimension(0);
  const int dim = Objet_U::dimension;

  //  int nb_positions;
  int face_bord;
  //nb_positions=0;
  double x = 0.,y= 0.,z= 0.;

  for (int som = 0; som < nb_pos_tot; som++)
    {
      const int num_elem = elem[som];
      if (num_elem >= 0)
        {
          if (maillage.sommet_ligne_contact(som))
            {
              x= pos(som,0);
              y= pos(som,1);
              if (dim==3)
                z= pos(som,2);
              face_bord = som_face[som];
              appliquer_reflexion_vitesse(x,y,z,som,face_bord,zone_vf,vitesse_p);
            }
        }
    }

}

//-Construction  d une base de vecteur vecteur_base(0), vecteur_base(1), vecteur_base(2)
// ou vecteur_base(0) est normal a la face ou l on veut modifier le vecteur vitesse
//-Inversion du signe de la composante normale du vecteur vitesse dans cette base
//-Calcul des nouvelles composantes du vecteur vitesse dans la base (e0,e1,e2)

void Transport_Marqueur_FT::appliquer_reflexion_vitesse(const double x, const double y, const double z,
                                                        const int som,int& face_bord,
                                                        const Zone_VF& zone_vf,
                                                        DoubleTab& vitesse_p)
{

  int dim = Objet_U::dimension;
  const DoubleTab& cg_faces = zone_vf.xv();
  DoubleTab vecteur_base(dim,dim);
  DoubleVect vit_base1(dim);

  //calcul de la norme du vecteur normal a la face vecteur_base(0)
  double norme = 0.;
  for (int j=0; j<dim; j++)
    norme += zone_vf.face_normales(face_bord,j)*zone_vf.face_normales(face_bord,j);
  norme = sqrt(norme);
  if (!est_egal(norme,0.))
    {
      for (int j=0; j<dim; j++)
        vecteur_base(0,j) = zone_vf.face_normales(face_bord,j)/norme;
    }
  else
    {
      Cerr<<"Problem encountered wih the norme calculation for the method "<<finl;
      Cerr<<"Maillage_FT_Disc::modifier_proprietes_particules"<<finl;
    }
  //construction du premier vecteur tangentiel a la paroi vecteur_base(1)
  //on prend un vecteur orthogonal construit a partir du point x y z
  //situe sur la face et du centre de gravite de la face
  vecteur_base(1,0) = x-cg_faces(face_bord,0);
  vecteur_base(1,1) = y-cg_faces(face_bord,1);
  if (dim==3)
    vecteur_base(1,2) = z-cg_faces(face_bord,2);
  norme = 0.;
  for (int j=0; j<dim; j++)
    norme += vecteur_base(1,j)*vecteur_base(1,j);
  norme = sqrt(norme);
  if (!est_egal(norme,0.))
    {
      for (int j=0; j<dim; j++)
        vecteur_base(1,j) = vecteur_base(1,j)/norme;
    }
  else
    {
      Cerr<<"Problem encountered wih the tan1 norme calculation for the method"<<finl;
      Cerr<<"Maillage_FT_Disc::modifier_proprietes_particules"<<finl;
    }
  //cas dimension 3
  //construction du second vecteur tangentiel a la paroi vecteur_base(2)
  //on prend le produit vectoriel entre les deux vecteurs norm et tan1
  if (dim==3)
    {
      vecteur_base(2,0)=vecteur_base(0,1)*vecteur_base(1,2)-vecteur_base(0,2)*vecteur_base(1,1);
      vecteur_base(2,1)=vecteur_base(0,2)*vecteur_base(1,0)-vecteur_base(0,0)*vecteur_base(1,2);
      vecteur_base(2,2)=vecteur_base(0,0)*vecteur_base(1,1)-vecteur_base(0,1)*vecteur_base(1,0);
    }

  //Calcul des composantes du vecteur vitesse dans la base vecteur_base
  for (int j=0; j<dim; j++)
    {
      for (int k=0; k<dim; k++)
        vit_base1(j) += vitesse_p(som,k)*vecteur_base(j,k);

    }

  //Inversion du signe de la composante normale
  vit_base1(0)  *= -1.;

  //Calcul des nouvelles composantes du vecteur vitesse dans la base (e0,e1,e2)
  for (int j=0; j<dim; j++)
    {
      vitesse_p(som,j) = 0.;
      for (int k=0; k<dim; k++)
        vitesse_p(som,j) += vit_base1(k)*vecteur_base(k,j);
    }

}

//Ajout de points a l objet maillage
void Transport_Marqueur_FT::ajouter_points(const Maillage_FT_Disc& un_maillage_inject)
{
  Maillage_FT_Disc& maillage = maillage_interface();
  int skip_facettes = 1;
  maillage.ajouter_maillage(un_maillage_inject,skip_facettes);
}


//Preparation des tableaux contenant les proprietes (vitesse_p, ...)
//pour pouvoir assurer le passage d un sous domaine a un autre dans le cas d un calcul parallele
void Transport_Marqueur_FT::preparer_tableaux_avant_transport(Maillage_FT_Disc& maillage,
                                                              Proprietes_part_vol& proprietes)
{
  int nb_particules = proprietes.nb_particules();
  if (mp_sum(nb_particules))
    {
      const Desc_Structure_FT& descripteur = maillage.desc_sommets();
      DoubleTab& vitesse_som = proprietes.vitesse_particules();
      DoubleTab& delta_vit_som = proprietes.delta_v();
      DoubleTab& temp_som = proprietes.temperature_particules();
      DoubleTab& masse_vol_som = proprietes.masse_vol_particules();
      DoubleTab& diametre_som = proprietes.diametre_particules();
      DoubleTab& volume_som = proprietes.volume_particules();

      maillage.preparer_tableau_avant_transport(vitesse_som,descripteur);
      maillage.preparer_tableau_avant_transport(delta_vit_som,descripteur);
      maillage.preparer_tableau_avant_transport(temp_som,descripteur);
      maillage.preparer_tableau_avant_transport(masse_vol_som,descripteur);
      maillage.preparer_tableau_avant_transport(diametre_som,descripteur);
      maillage.preparer_tableau_avant_transport(volume_som,descripteur);
    }
}

//Actualisation des tableaux contenant les proprietes (vitesse_p, ...)
//dans le cas ou une particule a change de sous domaine dans le cas d un calcul parallele
//Cette operation doit modifier l attribut nb_particules_ de proprietes
//dans le cas ou le nombre de particules a varie pour le sous domaine
void Transport_Marqueur_FT::update_tableaux_apres_transport(Maillage_FT_Disc& maillage,
                                                            Proprietes_part_vol& proprietes)
{
  int nb_particules = proprietes.nb_particules();
  if (mp_sum(nb_particules))
    {
      const Desc_Structure_FT& descripteur = maillage.desc_sommets();
      int nb_som = maillage.nb_sommets();
      DoubleTab& vitesse_som = proprietes.vitesse_particules();
      DoubleTab& delta_vit_som = proprietes.delta_v();
      DoubleTab& temp_som = proprietes.temperature_particules();
      DoubleTab& masse_vol_som = proprietes.masse_vol_particules();
      DoubleTab& diametre_som = proprietes.diametre_particules();
      DoubleTab& volume_som = proprietes.volume_particules();

      maillage.update_tableau_apres_transport(vitesse_som,nb_som,descripteur);
      maillage.update_tableau_apres_transport(delta_vit_som,nb_som,descripteur);
      maillage.update_tableau_apres_transport(temp_som,nb_som,descripteur);
      maillage.update_tableau_apres_transport(masse_vol_som,nb_som,descripteur);
      maillage.update_tableau_apres_transport(diametre_som,nb_som,descripteur);
      maillage.update_tableau_apres_transport(volume_som,nb_som,descripteur);

      proprietes.fixer_nb_particules(nb_som);
    }
}


void Transport_Marqueur_FT::creer_champ(const Motcle& motlu)
{
  Equation_base::creer_champ(motlu);
  if (motlu == "densite_particules")
    {
      if (!densite_particules_.non_nul())
        {
          //const & discr = ref_cast(Discret_Thyd, discretisation());
          const Discretisation_base& discr = probleme().discretisation();
          const Zone_dis_base& une_zone_dis = zone_dis().valeur();
          const double temps = schema_temps().temps_courant();
          Nom nom="densite_particules";
          Nom unite="sans_dimension";
          discr.discretiser_champ("champ_elem",une_zone_dis,nom,unite,1,temps,densite_particules_);
          champs_compris_.ajoute_champ(densite_particules_.valeur());
        }
    }

  if (motlu == "volume_particules")
    {
      if (!volume_particules_.non_nul())
        {
          //const Discret_Thyd& discr = ref_cast(Discret_Thyd, discretisation());
          const Discretisation_base& discr = probleme().discretisation();
          const Zone_dis_base& une_zone_dis = zone_dis().valeur();
          const double temps = schema_temps().temps_courant();
          Nom nom="volume_particules";
          Nom unite="sans_dimension";
          discr.discretiser_champ("champ_elem",une_zone_dis,nom,unite,1,temps,volume_particules_);
          champs_compris_.ajoute_champ(volume_particules_.valeur());
        }
    }
}

const Champ_base& Transport_Marqueur_FT::get_champ(const Motcle& nom) const
{
  if (nom=="densite_particules")
    {
      double temps = schema_temps().temps_courant();
      Champ_Fonc_base& ch_densite=ref_cast_non_const(Champ_Fonc_base,densite_particules_.valeur());
      DoubleTab& val_densite =  ch_densite.valeurs();
      calculer_valeurs_densite(val_densite);
      ch_densite.changer_temps(temps);
      return champs_compris_.get_champ(nom);
    }

  if (nom=="volume_particules")
    {
      double temps = schema_temps().temps_courant();
      Champ_Fonc_base& ch_vol=ref_cast_non_const(Champ_Fonc_base,volume_particules_.valeur());
      DoubleTab& val_vol =  ch_vol.valeurs();
      calculer_valeurs_volumes(val_vol);
      ch_vol.changer_temps(temps);
      return champs_compris_.get_champ(nom);
    }

  try
    {
      return Transport_Interfaces_FT_Disc::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }
  throw Champs_compris_erreur();
  REF(Champ_base) ref_champ;
  return ref_champ;
}


//Calcul des valeurs du champ postraitable densite_particules_
//Estimation du nombre de particules par maille
const DoubleTab& Transport_Marqueur_FT::calculer_valeurs_densite(DoubleTab& val_densite) const
{
  const Maillage_FT_Disc& ens_points = maillage_interface();
  if (ens_points.type_statut())
    {
      const ArrOfInt& sommet_elem = ens_points.sommet_elem();
      val_densite = 0.;
      int nb_som = ens_points.nb_sommets();
      int elem;
      for (int som=0; som<nb_som; som++)
        {
          elem = sommet_elem[som];
          val_densite(elem) += 1.;
        }
    }
  return val_densite;
}

//Calcul des valeurs du champ postraitable volume_particules_
//Estimation du volume occupe par les particules dans une maille
const DoubleTab& Transport_Marqueur_FT::calculer_valeurs_volumes(DoubleTab& val_volume) const
{
  const Maillage_FT_Disc& ens_points = maillage_interface();

  if (ens_points.type_statut())
    {
      const ArrOfInt& sommet_elem = ens_points.sommet_elem();
      const Proprietes_part_vol&  proprietes = proprietes_particules();
      const DoubleTab& volume_part = proprietes.volume_particules();

      val_volume = 0.;
      int nb_som = ens_points.nb_sommets();
      int elem;

      for (int som=0; som<nb_som; som++)
        {
          elem = sommet_elem[som];
          val_volume(elem) += volume_part(som, 0);
        }
    }
  return val_volume;
}

// Methode de travail de remplissage d'une FloatTab par un DoubleTabFT
inline void remplissage(const DoubleTab& tab, DoubleTab *ftab)
{
  const int nb_noeuds = tab.dimension(0);
  const int nb_compo = tab.dimension(1);
  ftab->resize(nb_noeuds, nb_compo);
  for (int som=0 ; som<nb_noeuds ; som++)
    for (int k=0 ; k<nb_compo ; k++)
      (*ftab)(som,k) = tab(som,k);
}
/*! @brief Cherche le champ discret aux interfaces dont le nom est "champ", et verifie qu'il peut etre postraite a la localisation demandee (loc).
 *
 *   Si oui on renvoie 1 et, si ftab est non nul, on remplit le champ ftab
 *   avec le champ demande.
 *   Si non, on renvoie 0.
 *   (la fonction est appelee avec ftab=0 lors de la lecture du postraitement,
 *    car on n'a pas besoin de la valeur du champ, on veut seulement verifier
 *    qu'il existe).
 *
 */
int Transport_Marqueur_FT::get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, DoubleTab *ftab) const
{
  int res = 1;

  const Motcle som = "sommets";            //postraitement possible uniquement aux sommets
  const Motcle elem = "elements";          //postraitement possible uniquement aux elements
  const Motcle bi = "elements et sommets"; //postraitement possible aux sommets et aux elements
  int nb_champs = 2;
  int proprietes_postraitables=1;
  // Condition pour postraiter les proprietes des particules
  // en attendant que les marqueurs aient des proprietes simples:
  // diametre nul
  // vitesse
  // etc...
  if (methode_calcul_vp_ == INTERPOLEE)
    proprietes_postraitables=0;

  if (proprietes_postraitables)
    nb_champs+=5;
  Motcles les_champs(nb_champs);
  {
    les_champs[0] = Postraitement_base::demande_description;
    les_champs[1] = "vitesse";
    if (proprietes_postraitables)
      {
        les_champs[2] = "delta_v";
        les_champs[3] = "temperature";
        les_champs[4] = "masse_volumique";
        les_champs[5] = "diametre";
        les_champs[6] = "volume";
      }
  }
  Motcles localisations(nb_champs);
  {
    localisations[0] = bi;
    // Pour des particules, on n'accepte que le postraitement aux sommets (1 particule=1 sommet!):
    for (int i=1; i<nb_champs; i++)
      localisations[i] = som;
  }

  int rang=les_champs.search(champ);

  if (rang==0)
    {
      Cerr<<"The real fields to be post-processed are :"<<finl;
      for (int i=1 ; i<nb_champs ; i++)
        {
          Cerr << " Fields("<<i<<") : "<< les_champs[i] << " # Localisations : " << localisations[i] << finl;
        }
      res = 0;
    }
  else if (rang==-1)     //test champ existe ?
    {
      //champ inexistant
      res = 0;
    }
  else if (! (localisations[rang]==bi
              || (localisations[rang]==som && loc==Postraitement_base::SOMMETS)
              || (localisations[rang]==elem && loc==Postraitement_base::ELEMENTS)) )   //test localisation autorisee ?
    {
      //localisation non autorisee
      res = 0;
    }
  else
    {
      if (ftab) // Si pointeur nul : ne pas calculer la valeur du champ.
        switch(rang)
          {
          case 1:
            {
              if (methode_calcul_vp_ == INTERPOLEE)
                {
                  DoubleTabFT vitesse;
                  calcul_vitesse_p(vitesse);
                  remplissage(vitesse, ftab);
                }
              else
                remplissage(proprietes_particules().vitesse_particules(), ftab);
              break;
            }
          case 2:
            {
              remplissage(proprietes_particules().delta_v(), ftab);
              break;
            }
          case 3:
            {
              remplissage(proprietes_particules().temperature_particules(), ftab);
              break;
            }
          case 4:
            {
              remplissage(proprietes_particules().masse_vol_particules(), ftab);
              break;
            }
          case 5:
            {
              remplissage(proprietes_particules().diametre_particules(), ftab);
              break;
            }
          case 6:
            {
              remplissage(proprietes_particules().volume_particules(), ftab);
              break;
            }
          default:
            Cerr << "Error for the method Transport_Marqueur_FT::get_champ_post_FT" << finl;
            exit();
          }
      res = 1;
    }

  return res;
}
/*! @brief Voir l'autre get_champ_post_FT.
 *
 * Cette fonction est specifique aux champs d'entiers.
 *
 */
int Transport_Marqueur_FT::get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, IntTab *itab) const
{
  return 0;
}
