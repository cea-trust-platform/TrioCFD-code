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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Maillage_FT_IJK.h
// Directory : $IJK_ROOT/src/FT/front_tracking
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Maillage_FT_IJK_included
#define Maillage_FT_IJK_included

#include <Objet_U.h>
#include <Maillage_FT_Disc.h>
#include <Linear_algebra_tools.h>
#include <IJK_Splitting.h>
#include <TRUSTTab.h>
class Domaine_dis;
class Parcours_interface;
/*! @brief : class Maillage_FT_IJK
 *
 *  <Description of class Maillage_FT_IJK>
 *
 *
 *
 */
class Maillage_FT_IJK : public Maillage_FT_Disc
{

  Declare_instanciable(Maillage_FT_IJK) ;

public:
  Maillage_FT_IJK(const Maillage_FT_IJK&) = default;
  void initialize(const IJK_Splitting&, const Domaine_dis&, const Parcours_interface&);
  const ArrOfInt&   compo_connexe_facettes() const
  {
    return compo_connexe_facettes_;
  };
  ArrOfInt&   compo_connexe_facettes_non_const()
  {
    return compo_connexe_facettes_;
  };

// Surcharge de Maillage_FT_Disc:
  void nettoyer_maillage() override;
  void transporter(const DoubleTab& deplacement) override;
  void deplacer_sommets(const ArrOfInt& liste_sommets_initiale,
                        const DoubleTab& deplacement_initial,
                        ArrOfInt& liste_sommets_sortis,
                        ArrOfInt& numero_face_sortie, int skip_facettes=0) override;
// Surcharge de maillage_ft_disc pour conserver les composantes connexes.
  void recopie(const Maillage_FT_Disc& source_mesh, Statut_Maillage niveau_copie) override;

// Surcharge de Maillage_FT_Disc :
// surcharge: initialise les composantes connexes et appelle la methode ajouter_maillage_IJK.
  void ajouter_maillage(const Maillage_FT_Disc& maillage_tmp,int skip_facettes=0) override;
// ajout et gestion du tableau des compo_connexes a partir du tableau compo_connex du maillage source
  void ajouter_maillage_IJK(const Maillage_FT_IJK& added_mesh);

  void lire_maillage_ft_dans_lata(const char *filename_with_path, int tstep,
                                  const char *geometryname);

// Ces methodes de conversion index <-> (i,j,k) ont ete deplaces dans la classe splitting.
#if 0
  int convert_ijk_cell_to_packed(int i, int j, int k) const
  {
    assert((i < 0 && j < 0 && k < 0)
           || (i >= 0 && i < nbmailles_euler_i_
               && j >= 0 && j < nbmailles_euler_j_
               && k >= 0 && k < nbmailles_euler_k_));
    if (i < 0)
      return -1;
    else
      return (k * nbmailles_euler_j_ + j) * nbmailles_euler_i_ + i;
  }
  Int3 convert_packed_to_ijk_cell(int index) const
  {
    Int3 ijk;
    if (index < 0)
      {
        ijk[0] = ijk[1] = ijk[2] = -1;
      }
    else
      {
        ijk[0] = index % nbmailles_euler_i_;
        index /= nbmailles_euler_i_;
        ijk[1] = index % nbmailles_euler_j_;
        index /= nbmailles_euler_j_;
        ijk[2] = index;
        assert(index < nbmailles_euler_k_);
      }
    return ijk;
  }

  void set_ijk_cell_index(int num_sommet, Int3 ijk)
  {
    sommet_elem_[num_sommet] = convert_ijk_cell_to_packed(ijk[0],ijk[1],ijk[2]);
  }

#endif
  void set_ijk_cell_index(int num_sommet, Int3 ijk)
  {
    const IJK_Splitting& s = ref_splitting_.valeur();
    sommet_elem_[num_sommet] = s.convert_ijk_cell_to_packed(ijk[0],ijk[1],ijk[2]);
  }
  Int3 get_ijk_cell_index(int num_sommet) const
  {
    const IJK_Splitting& s = ref_splitting_.valeur();
    int index = sommet_elem_[num_sommet];
    return s.convert_packed_to_ijk_cell(index);
  }

  void initialize_processor_neighbourhood();
// Surcharges de Maillage_FT_Disc :
  int check_sommets(int error_is_fatal = 1) const override;
  int check_mesh(int error_is_fatal = 1, int skip_facette_owner = 0, int skip_facettes = 0) const override;

  void calculer_compo_connexe_sommets(ArrOfIntFT& compo_connexe_sommets) const;

  void recopie_force_compo(const Maillage_FT_IJK& source_mesh, const int icompo);
// Methode qui modifie l'attribut compo_connexe_facettes_
  void set_composante_connexe(const int i_facette, const int icompo)
  {
    compo_connexe_facettes_[i_facette] = icompo;
  };

  double minimum_longueur_arrete() const;
  const REF(IJK_Splitting) ref_splitting() const
  {
    return ref_splitting_;
  }
protected:
  // Surcharge de Maillage_FT_Disc :
  void   calculer_costheta_minmax(DoubleTab& costheta) const override;

  const Maillage_FT_IJK& operator=(const Maillage_FT_IJK&)
  {
    Cerr << "Maillage_FT_IJK& operator=" << finl;
    Process::exit();
    return *this;
  }    // Interdit !
  void creer_facettes_virtuelles(const ArrOfInt& liste_facettes,
                                 const ArrOfInt& liste_pe,
                                 const ArrOfInt& facettes_send_pe_list,
                                 const ArrOfInt& facettes_recv_pe_list) override;

  REF(IJK_Splitting) ref_splitting_;

// Taille du domaine IJK sur chaque proc:
  int nbmailles_euler_i_;
  int nbmailles_euler_j_;
  int nbmailles_euler_k_;
// Tableaux des processeurs voisins :
  ArrOfInt voisinage_processeur_;
  ArrOfInt liste_processeurs_voisins_faces_;
  ArrOfInt liste_processeurs_voisins_aretes_;
  ArrOfInt liste_processeurs_voisins_coins_;

  // Pour chaque facette d'interface, numero de la composante connexe (numero de la bulle)
  ArrOfIntFT compo_connexe_facettes_;

};
#endif /* Maillage_FT_IJK_included */
