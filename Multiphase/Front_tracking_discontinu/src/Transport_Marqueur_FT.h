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
// File:        Transport_Marqueur_FT.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_Marqueur_FT_included
#define Transport_Marqueur_FT_included

#include <Transport_Interfaces_FT_Disc.h>
#include <Marqueur_FT.h>

enum Methode_calcul_vp { INTERPOLEE, BILAN_QDM };
enum Methode_couplage { SUIVI, ONE_WAY_COUPLING, TWO_WAY_COUPLING };

//Description
//Classe chargee de la gestion des particules (ponctuelles ou materielles)

//La presence de particules dans l ecoulement peut resulter :
//- d une distribution en condition initiale
//- d une injection periodique en temps
//- d une transformation d une inclusion d une phase dans une autre (cas diphasique)

//L integration des trajectoires des particules est faite a partir d un champ de vitesse resultant soit :
//- d une interpolation du champ de vitesse du fluide aux positions des particules (methode_calcul_vp_=INTERPOLEE
//  ou  methode_couplage_=SUIVI)
//- d un bilan de q.d.m. des particules (methode_calcul_vp_=BILAN_QDM) pour ce dernier cas on peut prendre en compte
//  uniquement l action du fluide sur les particules (methode_couplage_=ONE_WAY_COUPLING) ou tenir compte en plus
//  de la reaction des particules (methode_couplage_=TWO_WAY_COUPLING)

//La classe declenche ou realise les operations de lecture, d actualisation ou de postraitement
//concernant les attributs maillage_interface et proprietes_particules_ ou encore
//maillage_inject_ et proprietes_inject_ (portes par Transport_Interfaces_FT_Disc_interne)
//
//maillage_interface (ou maillage_inject_) designe ici un ensemble de points sans notion de facettes
//proprietes_particules_ (ou proprietes_inject_) est un objet portant les proprietes physiques
//                                                 des particules (voir Proprietes_part_vol)

//Rq : Dans la version actuelle, les particules ponctuelles sont distinguees des particules materielles
//     par l attribut nb_particules_ de Proprietes_part_vol (nb_particules_ = 0 dans le cas de particules ponctuelles)

//La classe porte aussi les proprietes du fluide (vit_fluide_som_ ...) estimees aux positions de particules.
//Ces quantites sont necessaires pour estimer les sources intervenant dans le bilan de q.d.m. des particules
//(methode_calcul_vp_=BILAN_QDM)

//Description succinte des operations mises en jeu pour realiser le transport des particules
//Actualisation :
//-initTimestep()  : injection, transformation
//                     evaluation des proprietes du fluide
//                     mise a jour des conditions limites
//                     calcul des sources (Trainee, ...)
//                     actualisation de delta_v
//-mettre_a_jour() : calcul de la vitesse des particules
//                     transport des particules



class Transport_Marqueur_FT : public Transport_Interfaces_FT_Disc
{
  Declare_instanciable_sans_constructeur(Transport_Marqueur_FT);
public:

  Transport_Marqueur_FT();
  void set_param(Param& titi) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  Entree& lire_cond_init(Entree& is) override;
  Entree& lire_cl(Entree&) override;
  Entree& lire_injection(Entree& is);
  Entree& lire_transformation(Entree& is);

  void    discretiser(void) override;
  int  preparer_calcul() override;
  void    completer() override;

  void  mettre_a_jour(double temps) override;
  bool initTimeStep(double dt) override;
  void imposer_cond_lim();
  static void appliquer_reflexion_vitesse(const double x, const double y, const double z,
                                          const int som,int& face_bord,
                                          const Zone_VF& zone_vf,
                                          DoubleTab& vitesse);

  //Effectue l integration d un ensemble de points (sans notion de facettes)
  void integrer_ensemble_lagrange(const double temps) override;

  //Calcul des proprietes du fluide aux positions des particule
  void calculer_proprietes_fluide_pos_particules(const Maillage_FT_Disc& ens_points);
  //Calcul de la vitesse des particules (interpolation ou bilan de qdm)
  void calcul_vitesse_p(DoubleTab& deplacement) const;

  //Resolution du bilan de qdm des particules
  void resoudre_edo(DoubleTab& vitesse_p, DoubleTab& source_stockage, const double delta_t);

  //Actualiser d une partie ou tout le tableau delta_v
  void update_delta_v(int n_deb,int n_fin,const Maillage_FT_Disc& ens_points,int calc=1);
  //Methodes pour le parallelisme des proprietes de particules
  void preparer_tableaux_avant_transport(Maillage_FT_Disc& maillage,
                                         Proprietes_part_vol& proprietes);
  void update_tableaux_apres_transport(Maillage_FT_Disc& maillage,
                                       Proprietes_part_vol& proprietes);

  inline int linject(double temps) const;
  void injecter(double temps);
  void injection(const Maillage_FT_Disc& maill_inject,const Proprietes_part_vol& propr_inject);
  void ajouter_points(const Maillage_FT_Disc& maillage_inject);

  inline int ltransfo(double temps) const;
  void transformer(double temps);
  void transformer_particules();

  void transformation(Maillage_FT_Disc& maillage,Proprietes_part_vol& proprietes);

  void calcul_proprietes_geometriques(const IntVect&        num_compo,
                                      const int        nb_compo_glob,
                                      const DoubleTab&        indic,
                                      ArrOfDouble&        volumes,
                                      DoubleTab&        positions);

  void detection_groupes_a_supprimer(const ArrOfDouble& volumes,
                                     const DoubleTab&    positions,
                                     ArrOfInt&               flags_compo_a_supprimer);

  void construction_ensemble_proprietes(const IntVect&                num_compo,
                                        const int                nb_compo,
                                        Maillage_FT_Disc&                ens_points,
                                        Proprietes_part_vol&        propri,
                                        const ArrOfInt&                flags_compos_a_supprimer,
                                        const DoubleTab&                positions,
                                        const ArrOfDouble&                volumes);

  //On surcharge les deux methodes suivantes pour quelle ne fasse rien
  int  sauvegarder(Sortie& ) const override;
  int  reprendre(Entree&) override;

  //Methodes de l interface des champs postraitables (champs euleriens)
  /////////////////////////////////////////////////////
  //Methode creer_champ pas codee a surcharger si necessaire
  void creer_champ(const Motcle& motlu) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  /////////////////////////////////////////////////////

  //methodes utilisees pour le post-traitement des quantites lagrangiennes
  int get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, FloatTab *ftab = 0) const override;
  int get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, IntTab *itab = 0) const override;

  const DoubleTab& calculer_valeurs_densite(DoubleTab& val_densite) const;
  const DoubleTab& calculer_valeurs_volumes(DoubleTab& val_volume) const;


  const Champ_Inc& inconnue(void) const override;     //renvoie un champ bidon
  Champ_Inc&        inconnue(void) override;

  //Methodes d acces aux tableaux contenant les proprietes du fluide
  inline const DoubleTab& vitesse_fluide() const;
  inline const DoubleTab& rho_fluide() const;
  inline const DoubleTab& visco_dyn_fluide() const;
  inline const DoubleTab& grad_pression() const;

  inline const double& temps_debut_integration() const;
  inline const double& dela_t() const;
  inline const int& resol_implicite() const;

  //Methode d acces au terme source stocke (somme de chacune des sources)
  inline const DoubleTab& source_stockage() const;

protected:

  Champ_Fonc densite_particules_;   //Exprime le nombre de particules par maille
  Champ_Fonc volume_particules_;    //Exprime le volume des particules par maille

  double t_debut_integr_;            //Instant de demarrage de l integration
  double t_debut_inject_;            //Instant de la premiere injection
  double t_debut_transfo_;          //Instant de la premiere transformation
  double t_derniere_inject_;            //Instant de la derniere injection realisee
  double t_derniere_transfo_;            //Instant de la derniere transformation realisee
  double dt_inject_;                    //periode d injection
  double dt_transfo_;                    //periode de transformation

  Noms nom_sz_transfo_;             //Noms des sous_zones ou l on impose eventuellement un transformation
  double diametre_min_;             //Diametre minimum pour declencher une transformation
  double beta_transfo_;             //Parametre multiplicatif du volume de controle

  double dt_p_;                            //Pas de temps pour la resolution du bilan de q.d.m.
  int nb_it_;                    //dt_p_ = dt_/nb_it_

  int implicite_;                    //Activation de la resolution implicite (pour masse ajoutee)

  // L'inconnue de cette equation n a pas de sens
  //Creation d un champ bidon
  Champ_Inc champ_bidon_;

  int phase_marquee_;                    //numero de la phase marquee par des particules (-1 pour cas monophasique)
  Nom nom_eq_interf_;                        //Pour recuperer l equation d interface entre les deux phases
  Methode_calcul_vp methode_calcul_vp_; //Interpolation (INTERPOLEE par defaut) ou resolution du bilan de qdm des particules (BILAN_QDM)
  Methode_couplage methode_couplage_;   //Pas de couplage (SUIVI par defaut)
  //Action des particules sur le fluide (ONE_WAY_COUPLING)
  //+action des particules sur le fluide (TWO_WAY_COUPLING)

  //Tableaux contenant les proprietes du fluide aux positions des particules
  DoubleTab vit_fluide_som_;
  DoubleTab rho_fluide_som_;
  DoubleTab visco_dyn_fluide_som_;
  DoubleTab grad_P_som_;

  int contrib_one_way;                //Pour supprimer l action du fluide sur les particules (contrib_one_way=0)
  //par defaut contrib_one_way=1

  DoubleTab source_stockage_;                //Tableau qui sert a stocker la contribution des termes sources



};


inline const double& Transport_Marqueur_FT::temps_debut_integration() const
{
  return t_debut_integr_;
}

inline const double& Transport_Marqueur_FT::dela_t() const
{
  return dt_p_;
}

inline const int& Transport_Marqueur_FT::resol_implicite() const
{
  return implicite_;
}


inline const DoubleTab& Transport_Marqueur_FT::vitesse_fluide() const
{
  return vit_fluide_som_;
}

inline const DoubleTab& Transport_Marqueur_FT::rho_fluide() const
{
  return rho_fluide_som_;
}

inline const DoubleTab& Transport_Marqueur_FT::visco_dyn_fluide() const
{
  return visco_dyn_fluide_som_;
}

inline const DoubleTab& Transport_Marqueur_FT::grad_pression() const
{
  return grad_P_som_;
}

inline const DoubleTab& Transport_Marqueur_FT::source_stockage() const
{
  return source_stockage_;
}

inline int Transport_Marqueur_FT::linject(double temps) const
{
  if ((sup_ou_egal(temps,t_debut_inject_)) && (sup_strict(temps-t_derniere_inject_,dt_inject_)))
    return 1;
  else
    return 0;
}

inline int Transport_Marqueur_FT::ltransfo(double temps) const
{
  if ((sup_ou_egal(temps,t_debut_transfo_)) && (sup_strict(temps-t_derniere_transfo_,dt_transfo_)))
    return 1;
  else
    return 0;
}

#endif
