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
// File:        Remaillage_FT.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/patch_168/1
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Remaillage_FT_included
#define Remaillage_FT_included

#include <Objet_U.h>
#include <Zone_dis.h>
#include <Ref_Zone_VF.h>
#include <FTd_tools.h>

class Maillage_FT_Disc;
class ArrOfIntFT;
class IntTabFT;
class ArrOfDoubleFT;
class DoubleTabFT;
class Zone_VF;
class Param;

// ====================================================================
// .DESCRIPTION        : class Remaillage_FT
//  Cette classe implemente les procedures de remaillage des interfaces pour le Front-Tracking :
//
// .SECTION voir aussi
//  Transport_Interfaces_FT_Disc Maillage_FT_Disc

class Remaillage_FT : public Objet_U
{
  Declare_instanciable_sans_constructeur(Remaillage_FT);
public:
  Remaillage_FT();
  void set_param(Param& p);
  void associer_zone(const Zone_dis& zone_dis);

  int supprimer_facettes_nulles(Maillage_FT_Disc& maillage) const;
  int traite_adherence(Maillage_FT_Disc& maillage) const;
  int traite_decollement(Maillage_FT_Disc& maillage, const DoubleTab& deplacement) const;

  int a_remailler(double temps, const Maillage_FT_Disc& maillage) const;
  int a_lisser(double temps) const;

  // Remaillage local avec changement des connectivites:
  void remaillage_local_interface(double temps, Maillage_FT_Disc& maillage);

  // Regularisations du maillage (deplacement des sommets sans changer la connectivite)
  void corriger_volume(Maillage_FT_Disc& maillage, ArrOfDouble& var_volume);
  void barycentrer_lisser_systematique(double temps, Maillage_FT_Disc& maillage);
  void barycentrer_lisser_apres_remaillage(Maillage_FT_Disc& maillage, ArrOfDouble& var_volume);

  int  sauvegarder(Sortie& ) const;
  int  reprendre(Entree&);

  double calculer_variation_volume(const Maillage_FT_Disc& maillage,
                                   const DoubleTab& position_initiale,
                                   ArrOfDouble& varVolume) const;

  void lisser_dvolume(const Maillage_FT_Disc& maillage,
                      ArrOfDouble& var_volume,
                      const int nb_iterations) const;

  virtual void regulariser_courbure(Maillage_FT_Disc& maillage,
                                    const double coeff,
                                    ArrOfDouble& dvolume) const;

protected:

  int tester_a_remailler(const Maillage_FT_Disc& maillage) const;

  double regulariser_maillage(Maillage_FT_Disc& maillage,
                              ArrOfDouble& var_volume,
                              const double facteur_barycentrage_tangent,
                              const double coeff_lissage,
                              const int nb_iter_barycentrage,
                              const int nb_iter_lissage,
                              const int max_nb_iter_correction_volume,
                              const double seuil_dvolume) const;

  double redistribuer_sommets(Maillage_FT_Disc&   maillage,
                              const double        relaxation_direction_tangente,
                              const double        relaxation_direction_normale,
                              ArrOfDouble& var_volume_impose,
                              ArrOfDouble& var_volume_obtenu) const;

  int calculer_barycentre_facettes_voisines(const Maillage_FT_Disc& maillage,
                                            DoubleTab& barycentres) const;

  int calculer_connectivites_sommetFacettes(const Maillage_FT_Disc& maillage,ArrOfInt& fa7VoisinesSom_index, IntTab& fa7VoisinesSom_data) const;
  int calculer_correction_deplacement(DoubleTab& deplacement,const ArrOfDouble& varVolume,const DoubleTab& deplacement_varVolume, const ArrOfDouble& norme2_deplacement_varVolume) const;

  int calculer_differentielle_volume(const Maillage_FT_Disc& maillage,
                                     DoubleTab& differentielle_volume) const;

  double calculer_variation_volume_facette_2D(int fa7, const Maillage_FT_Disc& maillage,const DoubleTab& position_initiale) const;
  double calculer_variation_volume_facette_3D(int fa7, const Maillage_FT_Disc& maillage,const DoubleTab& position_initiale) const;

  virtual double calculer_longueurIdeale2_arete(const Maillage_FT_Disc& maillage, int som0, double x, double y, double z) const;
  //  double calculer_longueurIdeale2_arete(const Maillage_FT_Disc& maillage, int som0, int som1) const;
  int supprimer_petites_aretes(Maillage_FT_Disc& maillage, ArrOfDouble& varVolume) const;
  int diviser_grandes_aretes(Maillage_FT_Disc& maillage) const;
  //int marquer_sommets_petites_aretes(Maillage_FT_Disc& maillage, ArrOfInt & tab_somSupp) const;
  int marquer_aretes(Maillage_FT_Disc& maillage, IntTab& tab_aretesMarquees, ArrOfInt& tab_somD, DoubleTab& tab_deplacement_somD, int drap) const;
  int inserer_tab_aretes(int& nb_tab_aretes,IntTab& tab_aretes,DoubleTab& tab_criteres,int pe0, int numOwner0, int pe1, int numOwner1, int face_bord1,int peRequete, int fa7_peR, int iarete_fa7_peR) const;
  int chercher_arete_tab(int tmp, const ArrOfInt& tab_index, const IntTab& tab_aretes,int pe0, int numOwner0, int pe1, int numOwner1) const;
  double calculer_volume_sommets_supprimes(const Maillage_FT_Disc& maillage, const ArrOfInt& tab_somSupp,ArrOfDouble& varVolume) const;

  int supprimer_facettes_bord(Maillage_FT_Disc& maillage) const;
  int supprimer_doublons_facettes(Maillage_FT_Disc& maillage) const;
  int permuter_aretes(Maillage_FT_Disc& maillage) const;
  double qualiteTriangle(const FTd_vecteur3& som0, const FTd_vecteur3& som1, const FTd_vecteur3& som2, double& aire) const;

  int nettoyer_maillage(Maillage_FT_Disc& maillage) const;

  REF(Zone_VF) refzone_VF_;

  double temps_;
  double temps_dernier_remaillage_;
  double temps_dernier_lissage_;

  double dt_remaillage_;
  double dt_lissage_;

  int nb_iter_remaillage_;
  int nb_iter_barycentrage_;
  int nb_iter_bary_volume_seul_;
  double seuil_dvolume_residuel_;
  double relax_barycentrage_;

  // Criteres de longueur ideale des aretes
  double critere_arete_;
  int impr_;
  // L'un de ces deux vaut -1, l'autre doit etre positif:
  double valeur_longueur_fixe_;
  double facteur_longueur_ideale_;
  int equilateral_;

  double variation_volume_;
  double surface_interface_;


  // Coefficient de "diffusion" (en quelque sorte la CFL du schema
  // de lissage).
  double lissage_courbure_coeff_;
  // Nombre d'iterations de lissage a faire lors des operations de lissage (dt_lissage_)
  int lissage_courbure_iterations_systematique_;
  // Nombre d'iterations de lissage a faire en cas de remaillage local ou global
  //  (declanche uniquement si le maillage est effectivement modifie)
  int lissage_courbure_iterations_si_remaillage_;
  // Ancien parametre de lissage (voir readOn)
  int lissage_courbure_iterations_old_;
};


#endif
