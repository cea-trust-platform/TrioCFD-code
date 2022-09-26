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
// File      : SurfaceVapeurIJKComputation.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <SurfaceVapeurIJKComputation.h>


void SurfaceVapeurIJKComputation::initialize(const IJK_Splitting& splitting)
{
  maillage_bulles_med_ = MEDCouplingUMesh::New("bubbles_surf_mesh", 2);
  desactive_med_ = false;
  compute_surf_mouillees_ = false;
  ref_splitting_ = splitting;
}


void SurfaceVapeurIJKComputation::get_maillage_MED_from_IJK_FT(MEDCouplingUMesh *maillage_bulles_mcu,
                                                               const Maillage_FT_IJK& maillage_bulles_ft_ijk)
{
  const int nbOfNodes = maillage_bulles_ft_ijk.nb_sommets();
  // Cerr << "nbre de sommet : " << nbOfNodes << finl;

  // Setting nodes coordinates of maillage_bulles_mcu
  auto coo_ptr = maillage_bulles_ft_ijk.sommets().addr();
  DAD coordsArr = DataArrayDouble::New();
  coordsArr->useArray(coo_ptr, false, MEDCoupling::DeallocType::CPP_DEALLOC, nbOfNodes, 3);
  maillage_bulles_mcu->setCoords(coordsArr);

  // Setting connectivity of maillage_bulles_mcu
  // getting connectivity from facettes
  const auto nbFacettes = maillage_bulles_ft_ijk.nb_facettes();
  const auto connectivity = maillage_bulles_ft_ijk.facettes();
  MCT mcu_connectivity = DataArrayIdType::New();
  mcu_connectivity->alloc(4 * nbFacettes, 1);
  auto mcu_connectivity_ptr = mcu_connectivity->getPointer();
  MCT mc_conn_idx = DataArrayIdType::New();
  mc_conn_idx->alloc(nbFacettes + 1, 1);
  auto mc_conn_idx_ptr = mc_conn_idx->getPointer();
  for (int ifa7 = 0; ifa7 < nbFacettes; ifa7++)
    {
      mcu_connectivity_ptr[4 * ifa7] = 3;
      mcu_connectivity_ptr[4 * ifa7 + 1] = connectivity(ifa7, 0);
      mcu_connectivity_ptr[4 * ifa7 + 2] = connectivity(ifa7, 1);
      mcu_connectivity_ptr[4 * ifa7 + 3] = connectivity(ifa7, 2);
      mc_conn_idx_ptr[ifa7] = 4 * ifa7;
    }
  mc_conn_idx_ptr[nbFacettes] = 4 * nbFacettes;
  maillage_bulles_mcu->setConnectivity(mcu_connectivity, mc_conn_idx, true);
  MCT no_use3 = maillage_bulles_mcu->zipConnectivityTraducer(0);
  MCT do_not_use5 = maillage_bulles_mcu->zipCoordsTraducer();
}

void SurfaceVapeurIJKComputation::set_maillage_MED(const Maillage_FT_IJK& maillage_ft_ijk)
{
  get_maillage_MED_from_IJK_FT(maillage_bulles_med_, maillage_ft_ijk);
  // const char fileName[] = "testExample_MEDCouplingFieldDouble_WriteVTK.vtk";
  // test sur le maillage
  Cerr << "Le maillage MED est de taille " << maillage_bulles_med_->getNumberOfNodes() << finl;
  const auto coords = maillage_bulles_med_->getCoords();
  for (int x = 0; x < 3; x++)
    {
      double mini = 100000000;
      double maxi = -100000000;
      for (int c = 0; c < coords->getNumberOfTuples(); c++)
        {
          const auto val = coords->getIJ(c, x);
          if (val < mini)
            mini = val;
          if (val > maxi)
            maxi = val;
        }
      Cout << "coord min dans la direction " << x << " : " << mini << finl;
      Cout << "coord max dans la direction " << x << " : " << maxi << finl;
    }
  MCT do_not_use = maillage_bulles_med_->convertDegeneratedCellsAndRemoveFlatOnes();
  auto connectivity = maillage_bulles_med_->getNodalConnectivity();
  const int nb_compo = static_cast<int>(connectivity->getNumberOfComponents());
  Cerr << "Connectivity mesh bulles : " << finl;
  Cerr << "Nb nodes : " << connectivity->getNumberOfTuples() << finl;
  Cerr << "Nb compo : " << nb_compo << finl;
  print_umesh_conn(connectivity);
  // Cerr << "Connectivity bis" << finl;
  // for (unsigned int comp = 0; comp < static_cast<int>(connectivity->getNumberOfComponents());
  //      comp++) {
  //   Cerr  << " [ ";
  //   for (int c = 0; c < connectivity->getNumberOfTuples(); c++) {
  //     Cerr << connectivity->getIJ(c, comp) << ", ";
  //   }
  //   Cerr << "]" << finl;
  // }
  maillage_bulles_med_->checkConsistencyLight();
  // maillage_bulles_med_->writeVTK(fileName, true);
  // Cout << "Le maillage a probablement ete ecrit qq part" << finl;
}

void SurfaceVapeurIJKComputation::slice_bubble(const double intersect_pt, const int dim, DataArrayIdType *cutcellsid,
                                               bool& plan_cut_some_bubble, MCU& mesh1dfil) const
{
  double intersect[3] = {0., 0., 0.};
  double vect_norm[3] = {0., 0., 0.};
  intersect[dim] = intersect_pt;
  vect_norm[dim] = 1.;
  std::vector<int> COO2KEEP(2);
  get_coo_to_keep(dim, COO2KEEP);

  // Création du maillage filaire coupé
  std::vector<int> selected_dims {dim};
  // MEDCouplingUMesh * mesh1dfil;
  try
    {
      mesh1dfil = maillage_bulles_med_->buildSlice3DSurf(intersect, vect_norm, 10.e-8, cutcellsid);
      mesh1dfil->checkConsistencyLight();
      const DataArrayDouble *coords = mesh1dfil->getCoords();
      for (int x = 0; x < 3; x++)
        {
          double mini = 10000;
          double maxi = -10000;
          // Cout << "Coords fil compo : " << x << finl;
          // Cout << "[ ";
          for (int c = 0; c < coords->getNumberOfTuples(); c++)
            {
              const auto val = coords->getIJ(c, x);
              // Cout << val << ", ";
              if (val < mini)
                mini = val;
              if (val > maxi)
                maxi = val;
            }
          // Cout << "]" << finl;
          Cout << "coord fil min dans la direction " << x << " : " << mini << finl;
          Cout << "coord fil max dans la direction " << x << " : " << maxi << finl;
        }
      const auto connectivity = mesh1dfil->getNodalConnectivity();
      const int nb_compo = static_cast<int>(connectivity->getNumberOfComponents());
      Cout << "Connectivity fil : " << finl;
      Cout << "Nb nodes : " << connectivity->getNumberOfTuples() << finl;
      Cout << "Nb compo : " << nb_compo << finl;
      print_umesh_conn(connectivity);
      Cout << "Youpi la coupe marche" << finl;
    }
  catch (INTERP_KERNEL::Exception& e)
    {
      mesh1dfil = MEDCouplingUMesh::New();
      plan_cut_some_bubble = false;
    }
  if (plan_cut_some_bubble)
    {
      mcIdType nbOfNodes;
      bool areNodesMerged;
      MCT oldNodes = mesh1dfil->mergeNodes(EPS_, areNodesMerged, nbOfNodes);
      MCT do_not_use = mesh1dfil->zipCoordsTraducer();
      MCT do_not_use2 = mesh1dfil->zipConnectivityTraducer(0);
      DataArrayDouble *coord = mesh1dfil->getCoords();
      std::vector<long unsigned int> COO2KEEPLONG(begin(COO2KEEP), end(COO2KEEP));
      DAD coords = coord->keepSelectedComponents(COO2KEEPLONG);
      mesh1dfil->setCoords(coords);
    }
  // return mesh1dfil;
}

// DONE: plusieurs bulles (plus que 1) peuvent être présente dans la cellule,
// auquel cas il faut le prendre en compte. Donc on doit retenir les n+1 phases
// dans une cellule qui a n compo connexes différentes.
void SurfaceVapeurIJKComputation::findCommonTuples(const DataArrayIdType *new_to_old_id, const int n_tot_mesh2d,
                                                   DataArrayIdType *tab_id_subcells, DataArrayIdType *tab_id_cut_cell) const
{
  // temp est un tableau temporaire qui fait la taille du nbre de cellules old
  // et qui est initialement rempli de -1. Il sert à identifier les cellules
  // diphasiques (une cellule qui a déjà été renseignée dans old et sur laquelle
  // on retombe en parcourant new_to_old est une cellule diphasique qu'il faut
  // donc renseigner dans tab_id_cut_cell et tab_id_subcells Son deuxieme indice
  // permet de savoir en quelle position était le tuple diphasique
  const auto n_new_to_old = new_to_old_id->getNumberOfTuples();
  const auto new_to_old_id_ptr = new_to_old_id->getConstPointer();

  DAI temp = DataArrayInt::New();
  temp->alloc(n_tot_mesh2d, 2);
  temp->fillWithValue(-1.);
  int *temp_ptr = temp->getPointer();
  // n_cut n'est pas le bon nombre ! Mais comme je n'ai aucune;
  // idée de sa valeur je prend la différence entre new et old qui me
  // donne le bon résultat si toutes les surfaces ne sont coupées que par une
  // bulle; n_cut est censé correspondre au nbre de cellules du maillage IJ;
  // traversée par l'interface de la bulle. Je reduis la taille du;
  // tableau ensuite.;
  auto n_cut = n_new_to_old - n_tot_mesh2d;
  assert(n_cut > 0);
  // tab_id_cut_cell est un tableau qui renvoie la correspondance entre l'indice
  // de la maille diphasique dans tab_id_subcells et son id deja defini DAI
  // tab_id_cut_cell = DataArrayInt::New();
  tab_id_cut_cell->alloc(n_cut, 1);
  tab_id_cut_cell->fillWithValue(-1);
  int *tab_id_cut_cell_ptr = tab_id_cut_cell->getPointer();
  // tab_id_subcells est un tableau de tuple des ids des bouts de cellules qui
  // sont en fait dans la meme cellule. La taille de ce tuple est celle du nbre
  // max de compo connexe plus 1 (la phase liquide). Le tuple est rempli au fur
  // et à mesure et complété par des -1 deja defini mais pas encore alloué
  tab_id_subcells->alloc(n_cut, max_authorized_nb_of_components_ + 1);
  tab_id_subcells->fillWithValue(-1);
  int *tab_id_subcells_ptr = tab_id_subcells->getPointer();
  int ind = 0;
  // ind est l'indice du dernier doublon qui avance au fur et à mesure de la
  // découverte d'un nouveau doublon
  for (int i_new = 0; i_new < n_new_to_old; i_new++)
    {
      const auto i_old = new_to_old_id_ptr[i_new];
      assert(i_old < n_tot_mesh2d);
      if (temp_ptr[2 * i_old] == -1)
        temp_ptr[2 * i_old] = i_new;
      else
        {
          // TODO: est-ce que c'est ind ou i_old ?
          if (tab_id_subcells_ptr[(max_authorized_nb_of_components_ + 1) * ind] == -1)
            {
              // la cellule n'est pas encore marquee comme diphasique, on ajoute les
              // deux subcells identifies pour le moment et on augmente le nbre de
              // cellule diph
              tab_id_cut_cell_ptr[ind] = i_old;
              tab_id_subcells_ptr[(max_authorized_nb_of_components_ + 1) * ind] = temp_ptr[2 * i_old];
              tab_id_subcells_ptr[(max_authorized_nb_of_components_ + 1) * ind + 1] = i_new;
              // On garde une trace du fait que ce tuple est a la ind eme position
              // dans tab_id
              temp_ptr[2 * i_old + 1] = ind;
              // on a effectivement un doublon de plus, si on vient de remplir le
              // tuple 0, maintenant ind doit valeur 1, soit le nbre de cellules
              // diphasiques trouvées
              ind += 1;
            }
          else
            {
              // la cellule avait deja été marquee comme diphasique, on va juste
              // modifier le bon tuple des subcells
              const int ind_tuple = temp_ptr[2 * i_old + 1];
              assert(ind_tuple >= 0);
              int last_compo_remplie = 1;
              // ca plante si on a plus de compo dans une cellule que nb_max_compo + 1
              for (int i_compo = 1; tab_id_subcells_ptr[(max_authorized_nb_of_components_ + 1) * ind_tuple + i_compo] > -1;
                   i_compo++)
                {
                  last_compo_remplie = i_compo;
                }
              assert(last_compo_remplie > 0);
              tab_id_subcells_ptr[(max_authorized_nb_of_components_ + 1) * ind_tuple + last_compo_remplie + 1] = i_new;
            }
        }
    }
  // {

  tab_id_cut_cell->reAlloc(ind);
  tab_id_subcells->reAlloc(ind);
}

void get_coo_to_keep(int d, std::vector<int>& COO2KEEP)
{
  if (d == 0)
    {
      COO2KEEP[0] = 1;
      COO2KEEP[1] = 2;
    }
  else if (d == 1)
    {
      COO2KEEP[0] = 2;
      COO2KEEP[1] = 0;
    }
  else if (d == 2)
    {
      COO2KEEP[0] = 0;
      COO2KEEP[1] = 1;
    }
  else
    {
      Cerr << "La dimension renseignée est trop grande" << finl;
    }
}

void get_coo_inv_to_keep(int d, std::array<int, 3>& COOINV)
{
  assert((d >= 0) && (d < 3));
  if (d == 0)
    {
      COOINV[0] = 2;
      COOINV[1] = 0;
      COOINV[2] = 1;
    }
  else if (d == 1)
    {
      COOINV[0] = 1;
      COOINV[1] = 2;
      COOINV[2] = 0;
    }
  else
    {
      COOINV[0] = 0;
      COOINV[1] = 1;
      COOINV[2] = 2;
    }
}

void SurfaceVapeurIJKComputation::get_IJK_ind_from_ind2d(const int i_dim, const int i_plan, const int i_2d, const int nx,
                                                         std::array<int, 3>& ijk_coo) const
{
  std::vector<int> COO2KEEP(2);
  std::array<int, 3> COOINV;
  get_coo_to_keep(i_dim, COO2KEEP);
  get_coo_inv_to_keep(i_dim, COOINV);
  // renvoie les coordonnees x et y correspondant au numero de la
  // cellule du maillage cartesien du plan XY
  // Pour rappel mon tableau est de taille un peu supérieur car
  // il représente un chmap aux faces et non au centre des cellules.
  const int y = i_2d / (nx - 1);
  const int x = i_2d - (nx - 1) * y;

  // On remet dans l'ordre ijk les indices selon la direction normale
  // au plan de coupe.
  // TODO: vérifier avec std::array<int, 3> coo {x, y, i_plan};
  std::array<int, 3> coo {x, y, i_plan};
  for (int c = 0; c < 3; c++)
    {
      ijk_coo[c] = coo[COOINV[c]];
    }
}

void SurfaceVapeurIJKComputation::check_if_vect_is_from_liquid2vapor(
  const FixedVector<IJK_Field_double, 3>& normale_of_interf,
  const DataArrayDouble *vector,
  const int dim, const int i_plan,
  const int nx, const DataArrayIdType *ids_old,
  DataArrayIdType *ids_new) const
{
  const int *ids_old_ptr = ids_old->getConstPointer();
  int *ids_new_ptr = ids_new->getPointer();
  const double *vector_ptr = vector->getConstPointer();
  const int n_compo_ids_IJ = static_cast<int>(ids_new->getNumberOfComponents());
  const int n_compo_v = static_cast<int>(vector->getNumberOfComponents());
  const int n_tup_v = vector->getNumberOfTuples();
  std::vector<int> coo2keep(2);
  get_coo_to_keep(dim, coo2keep);

  // {
  //   // debug part
  //   Cerr << "vector : " << finl;
  //   print_double_med(vector);
  // }

  for (int i_diph = 0; i_diph < n_tup_v; i_diph++)
    {
      std::array<int, 3> coo_ijk = {0, 0, 0};
      get_IJK_ind_from_ind2d(dim, i_plan, ids_old_ptr[i_diph], nx, coo_ijk);
      // debut debug
      // Cerr << "coo_ijk : { " << static_cast<int>(coo_ijk[0]) << ", " <<
      // static_cast<int>(coo_ijk[1]) << ", " << static_cast<int>(coo_ijk[2]) << "
      // }" << finl; fin debug
      double scal_prod = 0.;
      for (int c = 0; c < 2; c++)
        {
          // Je fais la moyenne de deux vecteurs normaux des cellules droites et
          // gauches pour avoir le vecteur de la normale a la face. Dans le cas ou
          // il ya plusieurs compos connexes je ne sais pas trop comment gérer la
          // situation.
          // TODO: PB!!! je vais chercher dans le field normale_of_interf l'element
          // 0,0,0 qui visiblement n'existe pas. bizarrement il y a un talbeau non
          // initialisé dans les tableaux suivants lors du tt premier passage à
          // chaque pas de temps. Je vais chercher la valeur pour la premiere compo
          // connexe dans l'element et j'espere qu'il s'agit bien de la mm dans
          // l'element d'à coté
          const double normale_of_interf_cell = normale_of_interf[coo2keep[c]](coo_ijk[0], coo_ijk[1], coo_ijk[2]);
          double normale_of_interf_next_cell;
          if (dim == 0)
            normale_of_interf_next_cell = normale_of_interf[coo2keep[c]](coo_ijk[0] + 1, coo_ijk[1], coo_ijk[2]);
          else if (dim == 1)
            normale_of_interf_next_cell = normale_of_interf[coo2keep[c]](coo_ijk[0], coo_ijk[1] + 1, coo_ijk[2]);
          else
            normale_of_interf_next_cell = normale_of_interf[coo2keep[c]](coo_ijk[0], coo_ijk[1], coo_ijk[2] + 1);

          const double normale_of_interf_moy = (normale_of_interf_cell + normale_of_interf_next_cell) / 2.;
          const double scal_prod_compo = vector_ptr[n_compo_v * i_diph + c] * normale_of_interf_moy;
          scal_prod += scal_prod_compo;
        }
      if (scal_prod < 0.)
        {
          // En fait on n'a plus besoin de la valeur de la premiere cellule si le
          // vecteur n'est pas dans le bon sens c'est que c'est la valeur de la
          // phase liquide ou que plusieurs compos connexes sont dans la cellule,
          // auquel cas on laisse tomber la rigueur du calcul dans la mesure ou on
          // ne prendra pas en compte cette complexite dans un operateur. int interm
          // = ids_IJ_cell_from_diph->getIJ(i_diph, 0);
          ids_new_ptr[n_compo_ids_IJ * i_diph] = ids_new_ptr[n_compo_ids_IJ * i_diph + 1];
        }
    }
}

void SurfaceVapeurIJKComputation::get_vect_from_sub_cells_tuple(const int dim, const DataArrayDouble *bary0,
                                                                const DataArrayIdType *cIcellsIdinMesh0,
                                                                const DataArrayIdType *cellsIdinMesh0,
                                                                DataArrayDouble *vect) const
{
  // vect est de la taille du nombre de cellules diphasiques, on ne regarde vect
  // que quand la face est coupée.
  const int n_vect = cIcellsIdinMesh0->getNumberOfTuples();
  vect->alloc(n_vect, 2);
  double *vect_ptr = vect->getPointer();
  // const int n_diph = bary0->getNumberOfTuples();
  // const int dim_bary0 = static_cast<int>(bary0->getNumberOfComponents());
  std::vector<int> coo_to_keep(2);
  get_coo_to_keep(dim, coo_to_keep);
  const int *cellsIdinMesh0_ptr = cellsIdinMesh0->getConstPointer();
  // const int* cIcellsIdinMesh0_ptr = cIcellsIdinMesh0->getConstPointer();
  const int n_cIcellsIdinMesh0 = cIcellsIdinMesh0->getNumberOfTuples();

  // {
  //   // debug part
  //   Cerr << "cellsIdinMesh0 : " << finl;
  //   print_int_med(cellsIdinMesh0);
  //   Cerr << "cIcellsIdinMesh0 : " << finl;
  //   print_int_med(cIcellsIdinMesh0);
  //   // end of the debug
  // }

  for (int c = 0; c < n_cIcellsIdinMesh0; c++)
    {
      // on considere arbitrairement les 2 premiers, je ne vois pas comment on
      // pourrait faire dans le cas ou il y en a plus.
      const int ind_phase0 = cellsIdinMesh0_ptr[(max_authorized_nb_of_components_ + 1) * c];
      assert(ind_phase0 >= 0);
      assert(ind_phase0 < bary0->getNumberOfTuples());
      const int ind_phase1 = cellsIdinMesh0_ptr[(max_authorized_nb_of_components_ + 1) * c + 1];
      assert(ind_phase1 >= 0);
      assert(ind_phase1 < bary0->getNumberOfTuples());
      // const int ind_maillage_2d_cart = cIcellsIdinMesh0_ptr[c];
      // Cerr << "Ind old : " << ind_maillage_2d_cart << finl;
      // Cerr << "Ind phase : " << ind_phase0 << ", " << ind_phase1 << finl;
      for (int dir = 0; dir < 2; dir++)
        {
          // const int coo_kept = coo_to_keep[dir];
          // assert(coo_kept <= dim_bary0);
          // Cerr << "Coo to keep : " << coo_kept << finl;
          // const double compo = bary0->getIJ(ind_phase1, coo_kept) -
          // bary0->getIJ(ind_phase0, coo_kept);
          const double compo = bary0->getIJ(ind_phase1, dir) - bary0->getIJ(ind_phase0, dir);
          // c plutot que ind_maillage_2d_cart
          vect_ptr[2 * c + dir] = compo;
        }
    }
}

void SurfaceVapeurIJKComputation::calculer_surfaces_et_barys_faces_mouillees_vapeur(
  const Maillage_FT_IJK& maillage_ft_ijk,
  const FixedVector<IJK_Field_double, 3>& normale_of_interf,
  FixedVector<IJK_Field_double, 3>& surfaces,
  FixedVector<FixedVector<IJK_Field_double, 3>, 3>& barycentres)
{
  // Lors de l'appel de cette méthode (une fois à chaque passage dans la boucle
  // principale en temps), on convertit le maillage au format MED.
  set_maillage_MED(maillage_ft_ijk);

  // Pass through 3 dimensions to cut the mesh of the bubble with multiple plans
  const auto& geom = surfaces[0].get_splitting().get_grid_geometry();
  const auto& splitting = surfaces[0].get_splitting();
  std::vector<int> N = {splitting.get_nb_elem_local(0) + 1, splitting.get_nb_elem_local(1) + 1,
                        splitting.get_nb_elem_local(2) + 1
                       };
  Cerr << "Dimensions du cas : ";
  for (int c = 0; c < 3; c++)
    Cerr << ", " << N[c];
  Cerr << finl;

  // std::vector<int> L_min = {geom.get_origin(0), geom.get_origin(1),
  // geom.get_origin(2)} std::vector<int> L_max = {geom.get_origin,
  // geom.get_origin(1), geom.get_origin(2)}

  // Dans chaque direction on va parcourir tous les plans pour découper les
  // bulles en tranches.
  std::vector<int> COO2KEEP(2);
  std::array<int, 3> COOINV;

  for (int d = 0; d < 3; d++)
    {
      get_coo_to_keep(d, COO2KEEP);
      get_coo_inv_to_keep(d, COOINV);

      // debut debug
      // Cerr << "COOINV : { " << static_cast<int>(COOINV[0]) << ", " <<
      // static_cast<int>(COOINV[1]) << ", " << static_cast<int>(COOINV[2]) << "
      // }" << finl; fin debug

      // Building plan mesh
      const auto X1 = COO2KEEP[0];
      DAD xarr = DataArrayDouble::New();
      xarr->alloc(N[X1], 1);
      xarr->iota();
      const auto dx = geom.get_constant_delta(X1);
      const auto xmin = geom.get_origin(X1) + dx * splitting.get_offset_local(X1);
      // Marche si le forecasting marche
      // xarr = xarr * dx + xmin;
      auto xarr_ptr = xarr->getPointer();
      for (int c = 0; c < xarr->getNumberOfTuples(); c++)
        {
          xarr_ptr[c] *= dx;
          xarr_ptr[c] += xmin;
        }

      const int X2 = COO2KEEP[1];
      DAD yarr = DataArrayDouble::New();
      yarr->alloc(N[X2], 1);
      yarr->iota();
      const auto dy = geom.get_constant_delta(X2);
      const auto ymin = geom.get_origin(X2) + dy * splitting.get_offset_local(X2);
      // yarr = yarr * dy + ymin;
      auto yarr_ptr = yarr->getPointer();
      for (int c = 0; c < yarr->getNumberOfTuples(); c++)
        {
          yarr_ptr[c] *= dy;
          yarr_ptr[c] += ymin;
        }

      MCC cmesh = MEDCouplingCMesh::New();
      cmesh->setCoords(xarr, yarr);
      MCU mesh = cmesh->buildUnstructured();
      mesh->setName("Plan_IJ_IK_JK");

      // Orient correctly grid
      mesh->changeSpaceDimension(3);
      const double vect[3] = {0., 0., -1.};
      mesh->orientCorrectly2DCells(vect, false);
      mesh->changeSpaceDimension(2);
      MCT do_not_use3 = mesh->zipCoordsTraducer();
      const int n_cell_tot = mesh->getNumberOfCells();
      if (Process::je_suis_maitre())
        MEDCoupling::WriteUMesh("mesh.med", mesh, true);

      DAD zarr = DataArrayDouble::New();
      zarr->alloc(N[d], 1);
      zarr->iota();
      const auto dz = geom.get_constant_delta(d);
      const auto zmin = geom.get_origin(d) + dz * splitting.get_offset_local(d);
      // Marche si le forecasting marche, sinon il faut parcourir en boucle
      auto zarr_ptr = zarr->getPointer();
      const int n_tup_zarr = zarr->getNumberOfTuples();
      for (int c = 0; c < n_tup_zarr; c++)
        {
          zarr_ptr[c] *= dz;
          zarr_ptr[c] += zmin;
        }

      // For each plan cut through the bubble mesh and compute the surfaces
      // and barycenters of the cut faces
      if (Process::je_suis_maitre())
        MEDCoupling::WriteUMesh("mesh_bulles.med", maillage_bulles_med_, true);
      for (int i_plan = 0; i_plan < N[d]; i_plan++)
        {
          MCT cellId = DataArrayIdType::New();
          bool plan_cut_some_bubble = true;
          MCU mesh1D = MEDCouplingUMesh::New();
          slice_bubble(zarr_ptr[i_plan], d, cellId, plan_cut_some_bubble, mesh1D);

          // Si on ne coupe pas de bulles, pas la peine de faire tout ça, la surface
          // vapeur reste nulle il n'y a que du liquide sur toutes ces faces. En
          // fait ça fait gagner en temps de calcul.
          if (plan_cut_some_bubble)
            {
              Cerr << finl;
              Cerr << "Point d'intersection : " << zarr_ptr[i_plan] << finl;
              Cerr << "Direction : " << d << finl;
              Cerr << finl;

              Cerr << "Mesh 1D coordinates : " << finl;
              print_double_med(mesh1D->getCoords());
              Cerr << "Mesh 1D connectivity" << finl;
              print_umesh_conn(mesh1D->getNodalConnectivity());
              // Cerr << "Mesh 1D connectivity bis" << finl;
              // for (unsigned int comp = 0; comp <
              // mesh1D->static_cast<int>(getNodalConnectivity()->getNumberOfComponents());
              //      comp++) {
              //   Cerr  << " [ ";
              //   for (int c = 0; c <
              //   mesh1D->getNodalConnectivity()->getNumberOfTuples(); c++) {
              //     Cerr << mesh1D->getNodalConnectivity()->getIJ(c, comp) << ", ";
              //   }
              //   Cerr << "]" << finl;
              // }
              MCT do_not_use4 = mesh1D->zipCoordsTraducer();
              MCT no_use1 = mesh1D->zipConnectivityTraducer(2);
              if (Process::je_suis_maitre())
                MEDCoupling::WriteUMesh("mesh1d.med", mesh1D, true);
              // récupérer le bon champ
              // auto cellId_ptr = cellId->getPointer();
              // const auto n_cell_diph = cellId->getNumberOfTuples();
              // MCU orthoProjected =
              // maillage_bulles_med_->buildPartOfMySelf(cellId_ptr,
              // cellId_ptr+n_cell_diph, true);

              MCT cellOrder = mesh1D->orderConsecutiveCells1D();

              // debut debug
              // Cerr << "Reordering mesh1D : " <<finl;
              // print_int_med(cellOrder);
              // fin debug

              MCU mesh1D_ordered = mesh1D->buildPartOfMySelf(cellOrder->begin(), cellOrder->end(), false);
              MCT no_use2 = mesh1D_ordered->zipConnectivityTraducer(2);
              MCT do_not_use6 = mesh1D_ordered->zipCoordsTraducer();
              Cerr << "Mesh 1D connectivity cell ordered" << finl;
              print_umesh_conn(mesh1D_ordered->getNodalConnectivity());
              // Cerr << "Mesh 1D connectivity index ordered" << finl;
              // print_int_med(mesh1D_ordered->getNodalConnectivityIndex());
              order_elem_mesh_filaire(mesh1D_ordered);
              Cerr << "Mesh 1D connectivity nodes ordered" << finl;
              print_umesh_conn(mesh1D_ordered->getNodalConnectivity());
              mesh1D_ordered->checkConsistencyLight();

              // MCU mesh1D_f = mesh1D_ordered->sortCellsInMEDFileFrmt();
              // DONE with false as 3rd arg: mesh1D_ordered->zipCoords();

              // debut debug
              // Cerr << "Coordonnees mesh1D ordered : " <<finl;
              // const auto coo = mesh1D_ordered->getCoords();
              // print_double_med(coo);
              // fin debug

              if (Process::je_suis_maitre())
                MEDCoupling::WriteUMesh("mesh1d_ordered.med", mesh1D_ordered, true);

              // {
              //   const auto connectivity = mesh1D_ordered->getNodalConnectivity();
              //   const int nb_compo = static_cast<int>(connectivity->getNumberOfComponents());
              //   Cerr << "Connectivity fil : " << finl;
              //   Cerr << "Nb nodes : " << connectivity->getNumberOfTuples() << finl;
              //   Cerr << "Nb compo : " << nb_compo << finl;
              //   print_umesh_conn(connectivity);
              // }

              // Compute the merge mesh and stock it
              MCU merge0, mesh_filaire_decoupe0;
              MCT cellIdInMesh2d0, mesh1dmergeId0;
              {
                MEDCouplingUMesh *merge0Ptr(0), *meshFilDecoupePtr(0);
                DataArrayIdType *cellIdInMesh2d0Ptr(0), *mesh1dmergeId0Ptr(0);

                MEDCouplingUMesh::Intersect2DMeshWith1DLine(mesh, mesh1D_ordered, EPS_, merge0Ptr, meshFilDecoupePtr,
                                                            cellIdInMesh2d0Ptr, mesh1dmergeId0Ptr);
                merge0 = merge0Ptr;
                mesh_filaire_decoupe0 = meshFilDecoupePtr;
                cellIdInMesh2d0 = cellIdInMesh2d0Ptr;
                mesh1dmergeId0 = mesh1dmergeId0Ptr;
              }
              if (Process::je_suis_maitre())
                {
                  // MCT otn = merge0->sortCellsInMEDFileFrmt();
                  merge0->setName("Merge0");
                  MEDCoupling::WriteUMesh("merge0.med", merge0, true);
                }

              // Compute surf and bary
              DAD bary0 = merge0->computeCellCenterOfMass();
              DAD surf0 = merge0->getMeasureField(false)->getArray(); // ->getValueOnMulti(bary0->getConstPointer(),
              // bary0->getNumberOfTuples());

              // Compute the table of cut cells in mesh2d
              MCT cellsIdinMesh0 = DataArrayIdType::New();
              MCT cIcellsIdinMesh0 = DataArrayIdType::New();
              findCommonTuples(cellIdInMesh2d0, n_cell_tot, cellsIdinMesh0, cIcellsIdinMesh0);
              // cellsIdinMesh1, cIcellsIdinMesh1 = findCommonTuples(cellIdInMesh2d1,;
              //     n_cell_tot);

              // Compute the vect from sub_cell0 to sub_cell1 in diph cells and check
              // if it goes from vapour to liquid-> If not switches sub_cell0 and
              // sub_cell1 ids->
              DAD vect0 = DataArrayDouble::New();
              // TODO la methode rempli le vecteur a partir des barycentres des deux
              // premieres sous cellules de chaque cellule diphasique.
              get_vect_from_sub_cells_tuple(d, bary0, cIcellsIdinMesh0, cellsIdinMesh0, vect0);

              // vect1 = get_vect_from_sub_cells_tuple(bary1, cIcellsIdinMesh1,
              //         cellsIdinMesh1, n_cell_tot);

              // Get ortho from 2d mesh of the bubble
              // get_ortho_per_cell(orthoProjected, N[X1]-1, N[Y1]-1, xmin, ymin,
              //         dx, dy, normale_of_interf_[old()]);

              // Check if orientation is Ok or switches it
              check_if_vect_is_from_liquid2vapor(normale_of_interf, vect0, d, i_plan, N[X1], cIcellsIdinMesh0, cellsIdinMesh0);
              // cIcellsIdinMesh1 = check_if_vect_is_from_liquid2vapor(vect1,
              // orthoPerCell,
              //                                                       cIcellsIdinMesh1,
              //                                                       cellsIdinMesh1);

              // Select values from the barys of the vapor phase and set it in the res
              // array
              // auto surf_ptr = surf0->getPointer();
              // const int len_surf =
              //     surf0->getNumberOfTuples() * static_cast<int>(surf0->getNumberOfComponents());
              // for (int c = 0; c < len_surf; c++)
              //   surf_ptr[c] /= (dx * dy);

              // // debug surf0
              // Cerr << "Surf0 : " << finl;
              // print_double_med(surf0);
              // //
              // // debug bary0
              // Cerr << "bary0 : " << finl;
              // print_double_med(bary0);

              // Cette boucle est un parcours des cellules diphasiques seulement.
              // Elle permet de remplir de manière efficace les tableaux surfaces et
              // barycentres
              const int n_diph = cellsIdinMesh0->getNumberOfTuples();
              const int n_compo_subcell = static_cast<int>(cellsIdinMesh0->getNumberOfComponents());
              const int *subcellPtr = cellsIdinMesh0->getConstPointer();
              const int *cIPtr = cIcellsIdinMesh0->getConstPointer();
              const double *surf0Ptr = surf0->getConstPointer();
              const double *bary0Ptr = bary0->getConstPointer();
              const int n_compo_bary = static_cast<int>(bary0->getNumberOfComponents());

              for (int c = 0; c < n_diph; c++)
                {
                  // TODO renvoie les coordonnees x et y correspondant au numero de la
                  // cellule du maillage cartesien du plan XY
                  const int i_cell_cart = cIPtr[c];
                  // On prend le premier de la ligne parce que pour rappel on les a
                  // réordonné afin d'avoir le premier qui correspond à la subcell
                  // vapeur
                  const int i_subcell = subcellPtr[c * n_compo_subcell];
                  std::array<int, 3> coo_inv {0, 0, 0};
                  get_IJK_ind_from_ind2d(d, i_plan, i_cell_cart, N[X1], coo_inv);

                  // J'attribue à la face d de la maille i,j,k
                  // la premiere surface de la sous cellule X1 X2
                  const double surf_vap = surf0Ptr[i_subcell];
                  surfaces[d](coo_inv[0], coo_inv[1], coo_inv[2]) = surf_vap;
                  // je dois remettre dans le bon ordre les deux autres composantes
                  const double bary_x1 = bary0Ptr[i_subcell * n_compo_bary];
                  barycentres[X1][d](coo_inv[0], coo_inv[1], coo_inv[2]) = bary_x1;
                  const double bary_x2 = bary0Ptr[i_subcell * n_compo_bary + 1];
                  barycentres[X2][d](coo_inv[0], coo_inv[1], coo_inv[2]) = bary_x2;
                  barycentres[d][d](coo_inv[0], coo_inv[1], coo_inv[2]) = zarr_ptr[i_plan];
                }
            }
        }
    }
  for (int c = 0; c < 3; c++)
    {
      for (int cc = 0; cc < 3; cc++)
        barycentres[c][cc].echange_espace_virtuel(barycentres[0][0].ghost());
      surfaces[c].echange_espace_virtuel(surfaces[0].ghost());
    }
}

int SurfaceVapeurIJKComputation::compute_surf_and_barys(
  const Maillage_FT_IJK& maillage_ft_ijk,
  const IJK_Field_double& indicatrice_ft,
  const FixedVector<IJK_Field_double, 3>& normale_of_interf,
  FixedVector<IJK_Field_double, 3>& surface_vapeur_par_face,
  FixedVector<FixedVector<IJK_Field_double, 3>, 3>& barycentre_vapeur_par_face) //next()
{
  for (int dir = 0; dir < 3; dir++)
    surface_vapeur_par_face[dir].data() = 0.;

  if ((!desactive_med_) and compute_surf_mouillees_)
    {
      Cout << "Calcul surfaces et barys" << finl;
      calculer_surfaces_et_barys_faces_mouillees_vapeur(
        maillage_ft_ijk,
        normale_of_interf,
        surface_vapeur_par_face,
        barycentre_vapeur_par_face);
      Cerr << "Le calcul des surfaces et barys est fini" << finl;
    }
  const int nb_surface_mouillees = rempli_surface_vapeur_par_face_interieur_bulles(surface_vapeur_par_face, indicatrice_ft);
  return nb_surface_mouillees;
}

int SurfaceVapeurIJKComputation::rempli_surface_vapeur_par_face_interieur_bulles(
  FixedVector<IJK_Field_double, 3>& surface_vapeur_par_face,
  const IJK_Field_double& indicatrice_ft) // next()
{
  const int ni = surface_vapeur_par_face[0].ni();
  const int nj = surface_vapeur_par_face[0].nj();
  const int nk = surface_vapeur_par_face[0].nk();
  const double eps_diph = 10.e-6;
  int n_faces_mouilles = 0;
  double surf_cell = 1.;
  for (int dir = 0; dir < 3; dir++)
    {
      surf_cell = 1.;
      for (int c = 0; c < 3; c++)
        {
          if (c != dir)
            surf_cell *= ref_splitting_->get_grid_geometry().get_constant_delta(c);
        }
      for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++)
          for (int k = 0; k < nk; k++)
            {
              double& surf_vap = surface_vapeur_par_face[dir](i, j, k);
              if (    (surf_vap / surf_cell > eps_diph)
                      and (surf_vap / surf_cell < 1. - eps_diph))
                {
                  n_faces_mouilles++;
                  // Cerr << "Cette face est diphasique ; dir " << dir << " i " << i
                  // << ", j " << j << ", k " << k << finl;
                }
              else
                {
                  double mean_indic = 0.;
                  if (dir == 0)
                    mean_indic = (indicatrice_ft(i - 1, j, k) + indicatrice_ft(i, j, k)) / 2.;
                  if (dir == 1)
                    mean_indic = (indicatrice_ft(i, j - 1, k) + indicatrice_ft(i, j, k)) / 2.;
                  if (dir == 2)
                    mean_indic = (indicatrice_ft(i, j, k - 1) + indicatrice_ft(i, j, k)) / 2.;
                  bool test = (mean_indic < 0.5);
                  // Alors on met la surface de la face comme vapeur pure
                  surf_vap = (test ? surf_cell : 0.);
                }
            }
    }
  return n_faces_mouilles;
}

void SurfaceVapeurIJKComputation::order_elem_mesh_filaire(MEDCouplingUMesh *mesh1D)
{
  DataArrayIdType *n = mesh1D->getNodalConnectivity();
  int *n_ptr = n->getPointer();
  DataArrayIdType *nI = mesh1D->getNodalConnectivityIndex();
  const int *nI_ptr = nI->getConstPointer();
  int n_cells = mesh1D->getNumberOfCells();
  Cerr << "Coordonnees de mesh1D : " << finl;
  print_double_med(mesh1D->getCoords());
  if (n_cells > 1)
    {
      int last_node, second_node;
      last_node = -1;
      int node0, node1, node2, node3;

      for (int i_cell = 0; i_cell < n_cells; i_cell++)
        {
          // Dans ce cas là on traite un debut de boucle par exemple.
          // Pour changer de boucle il faut atteindre à nouveau le
          // premier noeud (qui est donc le last_node).
          if (last_node == -1)
            {
              assert(i_cell + 1 < n_cells);
              node0 = n_ptr[nI_ptr[i_cell] + 1];
              node1 = n_ptr[nI_ptr[i_cell] + 2];
              node2 = n_ptr[nI_ptr[i_cell + 1] + 1];
              node3 = n_ptr[nI_ptr[i_cell + 1] + 2];
              if ((node0 == node2) || (node0 == node3))
                {
                  last_node = node1;
                  second_node = node0;
                  n_ptr[nI_ptr[i_cell] + 1] = last_node;
                  n_ptr[nI_ptr[i_cell] + 2] = second_node;
                }
              else if ((node1 == node2) || (node1 == node3))
                {
                  last_node = node0;
                  second_node = node1;
                }
              else
                throw "Erreur au debut d'une boucle deux element ne se suivent pas";
            }
          // Dans ce cas là on est dans un boucle, et on traite
          // l'élément suivant qu'on ordonne dans le mm sens que
          // le précédent.
          else
            {
              node0 = n_ptr[nI_ptr[i_cell] + 1];
              node1 = n_ptr[nI_ptr[i_cell] + 2];
              // Si c'est node1 qui est égal au second noeud de l'elem precedent
              // on inverse les noeuds 0 et 1.
              if (node1 == second_node)
                {
                  n_ptr[nI_ptr[i_cell] + 1] = node1;
                  n_ptr[nI_ptr[i_cell] + 2] = node0;
                  second_node = node0;
                }
              else if (node0 == second_node)
                {
                  second_node = node1;
                }
              // Si aucun des 2 noeud n'est dans l'element precedent c'est qu on a
              // fini un arc brise et qu on commence une nouvelle boucle.
              // On va retraiter cet element comme un debut de boucle.
              else
                {
                  // throw "Erreur les cellules ne sont pas faites avec des noeuds "
                  //       "consecutifs";
                  i_cell--;
                  last_node = -1;
                }
              // On a fini la boucle sur cet element, le prochain est dans une
              // nouvelle boucle.
              if (second_node == last_node)
                last_node = -1;
            }
        }
    }

  // double x1 = coords[n_compo_coo * node1];
  // double y1 = coords[n_compo_coo * node1 + 1];
  // double cellNormX = norms2D[n_compo_coo * cellNumber];
  // double cellNormY = norms2D[n_compo_coo * cellNumber + 1];
  // if ((cellNormX * (x1 - centerX) + cellNormY * ( y1-centerY ) ) < 0)
  // {
  //   n[nI[cellNumber] + 1] = node2;
  //   n[nI[cellNumber] + 2] = node1;
  // }
}

void print_int_med(const DataArrayIdType *data)
{
  const int n_tup = data->getNumberOfTuples();
  const int n_compo = static_cast<int>(data->getNumberOfComponents());
  const int *dataPtr = data->getConstPointer();
  if (n_compo == 1)
    {
      Cerr << "[ ";
      for (int i = 0; i < n_tup; i++)
        Cerr << dataPtr[i] << ", ";
      Cerr << "]" << finl;
    }
  else
    {
      for (int i = 0; i < n_tup; i++)
        {
          for (int j = 0; j < n_compo; j++)
            Cerr << dataPtr[i * n_compo + j] << ", ";
          Cerr << finl;
        }
    }
}

void print_double_med(const DataArrayDouble *data)
{
  const int n_tup = data->getNumberOfTuples();
  const int n_compo = static_cast<int>(data->getNumberOfComponents());
  const double *dataPtr = data->getConstPointer();
  if (n_compo == 1)
    {
      Cerr << "[ ";
      for (int i = 0; i < n_tup; i++)
        Cerr << dataPtr[i] << ", ";
      Cerr << "]" << finl;
    }
  else
    {
      for (int i = 0; i < n_tup; i++)
        {
          for (int j = 0; j < n_compo; j++)
            Cerr << dataPtr[i * n_compo + j] << ", ";
          Cerr << finl;
        }
    }
}

void print_umesh_conn(const DataArrayIdType *data)
{
  const int n_tup = data->getNumberOfTuples();
  const int n_compo = static_cast<int>(data->getNumberOfComponents());
  const int *dataPtr = data->getConstPointer();
  if (n_compo == 1)
    {
      // Cerr << "[ ";
      int i = 0;
      while (i < n_tup)
        {
          int elem_type = dataPtr[i];
          if (elem_type == 1)
            elem_type++;
          for (int j = 0; j < elem_type; j++)
            {
              Cerr << dataPtr[i + 1] << " ";
              i++;
            }
          i++;
          Cerr << ", ";
        }
      Cerr << finl;
      // Cerr << "]" << finl;
    }
  else
    {
      Cerr << "Attention on ne plot pas le tableau des connectivites" << finl;
    }
}

