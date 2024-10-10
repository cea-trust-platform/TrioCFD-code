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
// File      : SurfaceVapeurIJKComputation.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef SurfaceVapeurIJKComputation_included
#define SurfaceVapeurIJKComputation_included

#include <FixedVector.h>
#include <IJK_Field.h> // est-ce que j'en ai vraiment besoin ?
#include <Linear_algebra_tools_impl.h>
#include <Maillage_FT_IJK.h>
#include <medcoupling++.h>
#ifdef MEDCOUPLING_
#include <MEDCouplingCMesh.hxx>
#include <MEDCouplingFieldDouble.hxx>
#include <MEDCouplingUMesh.hxx>
#include <MEDLoader.hxx>
//#include <DataArrayDouble.hxx>
//#include <DataArrayInt.hxx>
#include <MCAuto.hxx>
using MEDCoupling::DataArrayDouble;
using MEDCoupling::DataArrayIdType;
using MEDCoupling::DataArrayInt;
using MEDCoupling::MCAuto;
using MEDCoupling::MEDCouplingCMesh;
using MEDCoupling::MEDCouplingFieldDouble;
using MEDCoupling::MEDCouplingUMesh;
using DAI = MCAuto<DataArrayInt>;
using DAD = MCAuto<DataArrayDouble>;
using MCT = MCAuto<DataArrayIdType>;
using MCC = MCAuto<MEDCouplingCMesh>;
using MCU = MCAuto<MEDCouplingUMesh>;
using MCF = MCAuto<MEDCouplingFieldDouble>;
#endif


static const int max_authorized_nb_of_groups_ = 3;
static const int max_authorized_nb_of_components_ = 3;
static const double EPS_ = 1.e-12;

void get_coo_to_keep(const int d, std::vector<int>& COO2KEEP);
void get_coo_inv_to_keep(const int d, std::array<int, 3>& COOINV);
void print_int_med(const DataArrayIdType *data);
void print_double_med(const DataArrayDouble *data);
void print_umesh_conn(const DataArrayIdType *data);

class SurfaceVapeurIJKComputation
{
public:
  SurfaceVapeurIJKComputation() {}
  ~SurfaceVapeurIJKComputation() {}
  void initialize(const IJK_Splitting& splitting);
  // Met à jour les valeurs de surface_vapeur_par_face_ et
  // barycentre_vapeur_par_face_.
  int compute_surf_and_barys(
    const Maillage_FT_IJK& maillage_ft_ijk,
    const IJK_Field_double& indicatrice_ft,
    const FixedVector<IJK_Field_double, 3>& normale_of_interf,
    FixedVector<IJK_Field_double, 3>& surface_vapeur_par_face,
    FixedVector<FixedVector<IJK_Field_double, 3>, 3>& barycentre_vapeur_par_face
  );
  // Permet de recuperer un mcu qui est une vue de maillage_bulles_ft_ijk. Il
  // n'y a pas de copie memoire, seulement des passages de case memoire.
  // splitting.get_grid_geometry().get_origin, get_constant_delta(DIRECTION_X),
  // ...
  static void get_maillage_MED_from_IJK_FT(MEDCouplingUMesh *maillage_bulles_mcu,
                                           const Maillage_FT_IJK& maillage_bulles_ft_ijk);
  void set_compute_surfaces_mouillees() { compute_surf_mouillees_ = true; }

protected:
  // Interdit le constructeur par copie (car constructeurs par copie interdits
  // pour parcours_ et autres
  // SurfaceVapeurIJKComputation(const SurfaceVapeurIJKComputation& x) :
  //   maillage_ft_ijk_(x.maillage_ft_ijk_)
  // {
  //   Cerr << "Erreur IJK_Interfaces(const IJK_Interfaces&)" << finl;
  //   Process::exit();
  // }
  const SurfaceVapeurIJKComputation& operator=(const SurfaceVapeurIJKComputation&)
  {
    Cerr << "Erreur IJK_Interfaces& operator=" << finl;
    Process::exit();
    return *this;
  }
  int rempli_surface_vapeur_par_face_interieur_bulles(FixedVector<IJK_Field_double, 3>& surface_vapeur_par_face, const IJK_Field_double& indicatrice_ft);
  // Cette methode appelle la methode statique get_maillage_MED_from_IJK_FT sur ses propres membres. Elle met donc a jour le maillage maillage_bulles_med_.
  void set_maillage_MED(const Maillage_FT_IJK& maillage_ft_ijk);

  // TODO: ecrire dans un fichier pour voir ce que ça donne. Utiliser nom_du_cas
  // ou toto.med Methode qui calcule les faces mouillees du maillage IJK par la
  // phase vapeur à l'aide de MEDCoupling. Les faces mouillees sont ensuite
  // utilisees pour faire des operations d'operateur discret de FT dans la
  // thermique principalement. surfaces: une vecteur comme celui de la vitesse
  // qui porte non pas les composantes x,y,z de la vitesse mais les surfaces des
  // faces I,J,K mouillees par la phase vapeur. barycentres: un vecteur qui
  // porte les valeurs sur les 3 faces I,J,K un vecteur de 3 composantes
  // rapporte à une cellule. Ces 3 composantes sont les coordonnees du
  // barycentre de la phase vapeur.
  void calculer_surfaces_et_barys_faces_mouillees_vapeur(
    const Maillage_FT_IJK& maillage_ft_ijk,
    const FixedVector<IJK_Field_double, 3>& normale_of_interf,
    FixedVector<IJK_Field_double, 3>& surfaces,
    FixedVector<FixedVector<IJK_Field_double, 3>, 3>& barycentres);

  // Cette methode calcule le vecteur que va d'un barycentre à l'autre entre
  // deux sous-cellules d'une cellule du maillage IJ. Params:
  //- barycentre: le tableau des barycentre pour le maillage merge
  //    (indice par les elements du maillage merge)
  //- ids_diph: le tableau des indices des sub-volumes pour chaque volume
  // diphasique
  //    (indice de la meme manière que le tableau qui donne la correspondance
  //    avec l'indice de l'element dans le maillage IJ).
  // Returns:
  //- le vecteur dont il est question, un tableau de taille egale au nombre de
  // cellules
  //    du mesh2d.
  void get_vect_from_sub_cells_tuple(const int dim, const DataArrayDouble *bary0,
                                     const DataArrayIdType *cIcellsIdinMesh0, const DataArrayIdType *cellsIdinMesh0,
                                     DataArrayDouble *vect) const;

  // Slice the bubble.
  // Params:
  //     intersect_pt : la coordonnee d'intersection
  //     dim: la direction normale à laquelle on fait la coupe
  // MEDCouplingUMesh* slice_bubble(const double intersect_pt, const int dim,
  // DataArrayIdType* cutcellsid, bool& plan_cut_some_bubble) const;
  void slice_bubble(const double intersect_pt, const int dim, DataArrayIdType *cutcellsid, bool& plan_cut_some_bubble,
                    MCU& mesh1dfil) const;

  // Cette methode trouve les doublons dans mesh_merge et leurs id.
  // Params:
  //    mesh_merge:
  //        le tableau de conversion qui donne l'indice old pour les indices new
  //        (donc new to old).
  //    n_tot_mesh2d:
  //        le nombre de maille du maillage IJ.
  // Returns:
  //    un tableau avec les indices des 2 cellules du maillage diphasique en
  //    couple et un tableau avec les cellules du maillage non decoupe
  //    correspondantes.
  void findCommonTuples(const DataArrayIdType *mesh_merge, const int n_tot_mesh2d, DataArrayIdType *tab_id_subcells,
                        DataArrayIdType *tab_id_cut_cells) const;

  // Cette methode ordonne les noeuds par ordre croissant cellule par cellule
  static void order_elem_mesh_filaire(MEDCouplingUMesh *mesh1D);

  // Cette methode permet de recuperer les indices IJK en partant des
  // informations concernant la geometrie IJK et le plan de coupe que l'on est
  // en train de parcourir. Elle renvoie les coordonnees IJK correspondant, dans
  // le bon ordre.
  void get_IJK_ind_from_ind2d(const int dim, const int i_plan, const int i_2d, const int nx,
                              std::array<int, 3>& ijk_coo) const;

  // Cette methode verifie si les subcells sont dans le bon sens (vapeur ->
  // liquide) et les echange sinon. Params:
  //- vector: TODO
  //- orthoprojected: TODO
  //- ids_diph: TODO
  //- ids_IJ_cell_from_diph: TODO.
  // Returns:
  //- TODO
  void check_if_vect_is_from_liquid2vapor(
    const FixedVector<IJK_Field_double, 3>& normale_of_interf,
    const DataArrayDouble *vector,
    const int dim, const int i_plan, const int nx,
    const DataArrayIdType *ids_diph,
    DataArrayIdType *ids_IJ_cell_from_diph) const;

  // Maillage au format MEDCouplingUMesh
  MCU maillage_bulles_med_;
  // Ref to maillage_ft_ijk
  // const FixedVector<IJK_Field_double, 3> normale_of_interf_;
  OBS_PTR(IJK_Splitting) ref_splitting_;
  bool desactive_med_;
  bool compute_surf_mouillees_;
};

#endif
