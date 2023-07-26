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
// File:        Marching_Cubes.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Marching_Cubes_included
#define Marching_Cubes_included

#include <TRUSTTabs.h>
#include <Maillage_FT_Disc.h>
#include <TRUSTArray.h>
#include <TRUSTTabFT_forward.h>
#include <TRUSTTabs_forward.h>
#include <TRUST_Ref.h>

class Domaine_VF;
class ArrOfBit;
class Desc_Structure_FT;
class Domaine;

class Marching_Cubes : public Objet_U
{
  Declare_instanciable_sans_constructeur(Marching_Cubes);
public:
  Marching_Cubes();

  void associer_domaine_vf(const Domaine_VF& domaine_vf);

  int construire_iso(const DoubleVect& valeurs_sommets,
                     double isovaleur,
                     Maillage_FT_Disc& maillage,
                     DoubleVect& indicatrice_approchee,
                     const Maillage_FT_Disc::AjoutPhase phase,
                     int ignorer_collision = 0) const;

  int construire_iso(const Nom& expression, double isovaleur,
                     Maillage_FT_Disc& maillage,
                     DoubleVect& indicatrice_approchee,
                     const Maillage_FT_Disc::AjoutPhase phase,
                     DoubleTab& eval_expression_sommets,
                     int ignorer_collision = 0) const;

protected:
  void remplir_data_marching_cubes(const Domaine& domaine);

  // Ces deux fonctions seraient mieux a leur place dans Joint ou Domaine...
  void remplir_renum_virt_loc(const Domaine& domaine);
  void renum_sommets_dist_loc(const int pe_voisin,
                              ArrOfInt& num_sommets) const;

  void calculer_signe(const DoubleVect& valeurs_sommets,
                      const double isovaleur,
                      ArrOfBit& signe) const;

  int construire_noeuds_et_facettes(const ArrOfBit& signe,
                                    IntTab& def_noeud,
                                    IntTab& facettes,
                                    DoubleVect& indicatrice_approchee,
                                    const Maillage_FT_Disc::AjoutPhase phase) const;

  void construire_noeuds_liste_faces(const ArrOfBit& signe,
                                     const IntTab& faces_sommets,
                                     const int nb_faces_a_traiter,
                                     const int numero_PE,
                                     IntTab& def_noeud) const;

  void construire_noeuds_joints(const ArrOfBit& signe,
                                IntTab& def_noeud) const;

  void trier_les_noeuds(IntTab& def_noeud) const;

  void construire_noeuds_uniques(IntTab& def_noeud,
                                 Maillage_FT_Disc& maillage) const;

  void calculer_coord_noeuds(const DoubleVect& valeurs_sommets,
                             const double isovaleur,
                             const IntTab& def_noeud,
                             Maillage_FT_Disc& maillage) const;

  void correspondance_espaces_distant_virtuel(const IntTab& def_noeud,
                                              Desc_Structure_FT& desc) const;

  // Champs remplis lors de l'association avec la domaine
  // *********** DEBUT
  REF(Domaine_VF) ref_domaine_vf_;

  int nb_sommets_element;
  // Definition des aretes de l'element de base eulerien
  // (exemple pour des triangles 2D {{0,1},{1,2},{2,0}}
  int nb_aretes_element;
  IntTabFT mcubes_def_aretes;
  // Definition des aretes des faces
  int nb_sommets_par_face;
  int nb_aretes_faces;
  IntTabFT mcubes_def_aretes_faces;
  // Definition des facettes a creer en fonction du cas marching_cubes
  // voir Marching_cubes_data.h
  int nb_sommets_facette;
  ArrOfIntFT mcubes_index_facettes;
  ArrOfIntFT mcubes_facettes;
  ArrOfIntFT mcubes_nb_facettes;

  // La structure suivante est utilisee dans renum_sommets_dist_loc :
  VECT(IntTab) renum_virt_loc_;
  // Pour i=0..nproc(), indice_joint[i] est l'indice du joint eulerien avec le proc i
  // (-1 si pas de joint avec le proc i)
  ArrOfInt indice_joint_;
  // *********** FIN

  // Taille du tableau a la derniere construction d'une iso...
  mutable int last_def_noeud_size;

private:
  Marching_Cubes(const Marching_Cubes&);  // Interdit !
  const Marching_Cubes& operator=(const Marching_Cubes&);   // Interdit !
};

#endif
