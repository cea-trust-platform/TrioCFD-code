/****************************************************************************
 * Copyright (c) 2019, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Intersection_Interface_ijk.h
// Directory : $IJK_ROOT/src/FT/
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Interface_ijk_included
#define Interface_ijk_included

#include <Objet_U.h>
#include <IJK_Field.h>
#include <IJK_Splitting.h>
#include <IJK_Field_forward.h>
#define NEIGHBOURS_I {-1, 1, 0, 0, 0, 0}
#define NEIGHBOURS_J {0, 0, -1, 1, 0, 0}
#define NEIGHBOURS_K {0, 0, 0, 0, -1, 1}
#define NEIGHBOURS_FACES_I {0, 1, 0, 0, 0, 0}
#define NEIGHBOURS_FACES_J {0, 0, 0, 1, 0, 0}
#define NEIGHBOURS_FACES_K {0, 0, 0, 0, 0, 1}

#define LOCAL_EPS 1e-12
/*! @brief : class Interface_ijk_face
 *
 *  Cette classe regroupe les information concernant les mailles qui sont
 *  coupées par l'interface (leur numéro par exemple).
 */

class IJK_Interfaces;

class Intersection_Interface_ijk : public Objet_U
{
  Declare_base_sans_constructeur(Intersection_Interface_ijk);
public:
  Intersection_Interface_ijk()
  {
    projected_on_interface_flag_ = false;
    n_diph_ = 1;
    interfaces_=nullptr;
    splitting_=nullptr;
  };
  virtual int initialize(const IJK_Splitting& splitting,
                         const IJK_Interfaces& interfaces) = 0;

  // Cette méthode permet de récupérer à partir des positions projetées sur
  // l'interface des positions à une certaine distance, dans la direction
  // normale. Ensuite on interpolera le field à ces positions de manière
  // monophasique car on sera loin de l'interface (et on ne prend pas en compte
  // le cas ou les bulles sont proches).
  static void get_position_interpolation_normal_interf(
    const DoubleTab& position_on_interf, const DoubleTab& normal_on_interf,
    const double dist, DoubleTab& positions);

  int n() const { return n_diph_; };

  const DoubleTab& pos_interf() const { return positions_on_interf_; };
  const DoubleTab& norm_interf() const { return normal_on_interf_; };
  const DoubleTab& dist_interf() const { return dist_to_interf_; };

  virtual void set_pas_a_jour() { projected_on_interface_flag_ = false; };
  virtual bool is_a_jour() const { return projected_on_interface_flag_; };

  // On projete le point bary_vap qui est localisé sur une face sur le plan
  // passant par bary, selon la normale qui est donnee. Ce point est appelé ici
  // bary_vap, mais ça peut très bien être bary_liqu ou n'importe quoi.
  static void projete_interface(const Vecteur3& normale,
                                const Vecteur3& point_plan, const Vecteur3& coo,
                                Vecteur3& posit_proj);
  static void distance_point_point(
    const Vecteur3& point_1,
    const Vecteur3& point_2,
    double& distance);
  static void distance_point_point_signe(
    const Vecteur3& point_1,
    const Vecteur3& point_2,
    const Vecteur3& normal_interf,
    double& distance);

protected:
  // Intersection_Interface_ijk_face n'est pas propriétaire de ces objets, pas
  // de gestion de la mémoire.
  const IJK_Interfaces *interfaces_;
  /*
   * TODO: create pointers that point to object
   * shared with IJK_Interfaces ?
   */
  const IJK_Splitting *splitting_;

  // Le nombre de cellules diphasiques.
  int n_diph_;
  // La projection de chaque bary de face mouille sur l'interface
  DoubleTab positions_on_interf_;
  // La distance a l'interface
  DoubleTab dist_to_interf_;
  // La normale utilisée pour faire la projection
  DoubleTab normal_on_interf_;
  // Flag qui est remis a 0 a chaque maj des projections des barys des
  // faces mouillées sur l'interface
  bool projected_on_interface_flag_;

  // Ici on recupere l'interface moyenne sur une cellule telle quelle etait au debut du
  // pas de temps.
  void get_mean_interface_cell(const int elem, Vecteur3& normale, Vecteur3& bary) const;
};

class Intersection_Interface_ijk_face : public Intersection_Interface_ijk
{
  Declare_instanciable_sans_constructeur(Intersection_Interface_ijk_face);
public:
  Intersection_Interface_ijk_face() : Intersection_Interface_ijk() {};
  int initialize(const IJK_Splitting& splitting,
                 const IJK_Interfaces& interfaces) override;

  // const int& operator()(int i_diph) const { return ijkf_interfaces_(i_diph); }
  const int& operator()(int i, int j, int k, int dir) const
  {
    return idiph_ijk_[dir](i, j, k);
  }

  void maj_interpolation_coo_on_interfaces();

protected:
  // Le tableau qui donne pour chaque face coupée par l'interface :
  // - les coordonnees ijk dans le field
  // - la face correspondant dans cette maille (FACE I/J/K)
  IntTab ijkf_interfaces_;
  // Le tableau qui donne pour un triplet ijk le numéro de la cellule diphasique
  // correspondante (-1. si la cellule n'est pas diphasique).
  FixedVector<IJK_Field_int, 3> idiph_ijk_;

  // Calcul les positions des points sur l'interface.
  // Les centres de grav des faces mouillées sont projecté
  // sur l'interface.
  // L'interface est construite par moyenne sur les deux cellules qui
  // correspondent à la face.
  void calcul_projection_bary_face_mouillee_interface_moy(
    DoubleTab& positions, IntTab& indices, DoubleTab& normales_de_la_proj, DoubleTab& distance_barys_interface);

  // La on moyenne les deux interfaces des cellules mitoyennes de la face.
  // Cette face doit être traversée par l'interface.
  void compute_mean_interface_face(const int i, const int j, const int k,
                                   const int dir, Vecteur3& normale, Vecteur3& bary) const;

};

class Intersection_Interface_ijk_cell : public Intersection_Interface_ijk
{
  Declare_instanciable_sans_constructeur(Intersection_Interface_ijk_cell);
public:
  Intersection_Interface_ijk_cell() : Intersection_Interface_ijk() {};
  int initialize(const IJK_Splitting& splitting,
                 const IJK_Interfaces& interfaces) override;

  // int& operator()(const int i,const int j,const int k);
  // int& operator()(const int i_diph);
  // const int& operator()(int i_diph) const { return ijk_interfaces_(i_diph); }
  const int& operator()(int i_diph, int dir) const { return ijk_interfaces_(i_diph, dir); }
  const int& operator()(int i, int j, int k) const { return idiph_ijk_(i, j, k); }
  const DoubleTab& pos_pure_faces_interf() const { return positions_pure_faces_on_interf_; };
  const DoubleTab& dist_pure_faces_interf() const { return dist_pure_faces_to_interf_; };
  const int& get_ijk_pure_face_neighbours(int i_diph, int neighbour) const { return ijk_pure_face_neighbours_(i_diph, neighbour); };
  FixedVector<DoubleVect, 3>& get_set_ijk_pure_face_to_correct() { return ijk_pure_face_to_correct_; };

  int get_nb_diph()
  {
    return ijk_interfaces_.dimension(0);
  }
  void update_interpolations_cell_centres_on_interface();
  void update_interpolations_cell_faces_on_interface();

  void update_interpolations_cell_centres_on_interface_new();
  void update_interpolations_cell_faces_on_interface_new();

  void update_interpolations_cell_centres_on_interface(const IJK_Field_double& interfaces);
  void update_interpolations_cell_faces_on_interface(const IJK_Field_double& interfaces);

  void set_pas_a_jour() override { projected_on_interface_flag_ = false; face_centres_projected_on_interface_flag_ = false;};
  bool is_a_jour() const override { return (projected_on_interface_flag_ && face_centres_projected_on_interface_flag_); };


protected:
  // Le tableau qui donne pour chaque face coupée par l'interface :
  // - les coordonnees ijk dans le field
  IntTab ijk_interfaces_;
  /*
   * Store the pure neighbouring cells
   * -1 if the neighbouring cells in not pure
   */
  IntTab ijk_pure_face_neighbours_;
  FixedVector<DoubleVect, 3> ijk_pure_face_to_correct_;
  // Le tableau qui donne pour un triplet ijk le numéro de la cellule diphasique
  // correspondante (-1. si la cellule n'est pas diphasique).
  IJK_Field_int idiph_ijk_;

  bool face_centres_projected_on_interface_flag_ = false;

  DoubleTab positions_pure_faces_on_interf_;
  DoubleTab dist_pure_faces_to_interf_;

  // Calcul les positions des points sur l'interface.
  // Les centres de cellules diphasiques sont projecté
  // sur l'interface.
  // L'interface est construite par moyenne sur la cellule
  void calcul_projection_centre_sur_interface_moy(const IJK_Field_double& indicatrice,
                                                  DoubleTab& positions, IntTab& indices,
                                                  DoubleTab& normales_de_la_proj,
                                                  DoubleTab& distance_centre_interface);

  void calcul_projection_centre_faces_sur_interface_moy(const IJK_Field_double& indicatrice,
                                                        const IntTab& indices,
                                                        const DoubleTab& normale_interf,
                                                        DoubleTab& positions,
                                                        IntTab& indices_voisins,
                                                        FixedVector<DoubleVect, 3>& indices_faces_corrections,
                                                        DoubleTab& distance_centre_faces_interface) const;

  void compute_face_to_correct();
};

#endif
