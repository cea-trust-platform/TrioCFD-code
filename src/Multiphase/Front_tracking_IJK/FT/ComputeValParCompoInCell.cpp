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
// File      : ComputeValParCompoInCell.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <ComputeValParCompoInCell.h>
#ifndef NUM_COMPO_INVALID
#define NUM_COMPO_INVALID (-2000000000)
#endif

// Pour remplacer calculer_normale_et_bary_element_par_compo
void ComputeValParCompoInCell::calculer_moyennes_interface_element_pour_compo(
  const int num_compo,
  const int elem,
  double& surface_tot,
  Vecteur3& normale,
  Vecteur3& bary
) const
{
  const Intersections_Elem_Facettes& intersections = mesh_->intersections_elem_facettes();
  const IntTab& facettes = mesh_->facettes();
  const DoubleTab& sommets = mesh_->sommets();
  const ArrOfDouble& surface_facettes = mesh_->get_update_surface_facettes();
  const DoubleTab& normale_facettes = mesh_->get_update_normale_facettes();
  const ArrOfInt& compo_connexe = mesh_->compo_connexe_facettes();

  int index = intersections.index_elem()[elem];
  surface_tot = 0.;
  normale = 0.;
  bary = 0.;

  if (index < 0)
    return; // Aucune facette dans cet element.

  // Boucle sur les facettes qui traversent l'element elem :
  int count = 0;
  int num_som = 0;
  int fa7 = 0;
  while (index >= 0)
    {
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
      fa7 = data.numero_facette_;
      const int icompo = compo_connexe[fa7];
      if (icompo == num_compo)
        {
          const double surface_facette = surface_facettes[fa7];
          const double surf = data.fraction_surface_intersection_ * surface_facette;
          // Les coordonnees du barycentre de la fraction de facette :
          Vecteur3 coord_barycentre_fraction(0., 0., 0.);
          for (int dir = 0; dir < 3; dir++)
            {
              const double nx = normale_facettes(fa7, dir);
              normale[dir] += nx * surf;
            }
          for (int isom = 0; isom < 3; isom++)
            {
              // numero du sommet dans le tableau sommets
              num_som = facettes(fa7, isom);
              // Coordonnees barycentriques du centre de gravite de l'intersection
              // par rapport aux trois sommets de la facette.
              const double bary_som = data.barycentre_[isom];
              // pour avoir la courbure au centre du triangle
              for (int dir = 0; dir < 3; dir++)
                coord_barycentre_fraction[dir] += bary_som * sommets(num_som, dir);
            }
          coord_barycentre_fraction *= surf;
          bary += coord_barycentre_fraction;
          surface_tot += surf;
        }
      index = data.index_facette_suivante_;
      count++;
    }

  if (surface_tot > 0.)
    {
      normale *= 1. / surface_tot;
      bary *= 1. / surface_tot;
    }
  else
    {
      Cerr << " Trouble in "
           "ComputeValParCompoInCell::calculer_normale_et_bary_element_pour_compo."
           << finl;
      Cerr << "L'element " << elem << " contient des facettes de surface totale "
           << surface_tot << " nulle!" << finl;
      for (int dir = 0; dir < 3; dir++)
        {
          bary[dir] = sommets(num_som, dir);
          normale[dir] = normale_facettes(fa7, dir);
        }
      Cerr << "bary: " << bary[0] << " " << bary[1] << " " << bary[2] << finl;
      Cerr << "normale: " << normale[0] << " " << normale[1] << " " << normale[2] << finl;
      Process::exit();

    }

  const double norm = normale[0] * normale[0] + normale[1] * normale[1] + normale[2] * normale[2];
  if (norm < 0.9)
    {
      Cerr << " Dans calculer_normale_et_bary_element_pour_compo." << finl;
      Cerr << "Petite : " << count << " facettes dans l'element " << elem << ". surface_tot = " << surface_tot
           << "Norm**2 = " << norm << finl;
    }
}

void ComputeValParCompoInCell::calculer_moy_par_compo(
#ifdef SMOOTHING_RHO
  const double delta_rho,
#endif
  IJK_Field_int& nb_compo_traversante,
  FixedVector<IJK_Field_int, max_authorized_nb_of_components_>& compos_traversantes,
  FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& normale_par_compo,
  FixedVector<IJK_Field_double, 3 * max_authorized_nb_of_components_>& bary_par_compo,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& indic_par_compo,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& surface_par_compo
) const
{

  Cout << "Calcul moyennes de l'interface dans chaque cellule par compo" << finl;

  // Initialisation a NUM_COMPO_INVALID signifie aucune compo presente.
  nb_compo_traversante.data() = NUM_COMPO_INVALID;
  // Peut-etre pas necessaire
  for (int c = 0; c < max_authorized_nb_of_components_; c++)
    {
      // TODO: AYM verifier l'init
      compos_traversantes[c].data() = NUM_COMPO_INVALID;
      surface_par_compo[c].data() = 0.;
      for (int compo = 0; compo < 3; compo++)
        {
          normale_par_compo[3 * c + compo].data() = 0.;
          bary_par_compo[3 * c + compo].data() = 0.;
        }
    }

  // Boucle sur les elements:
  const int ni = ref_splitting_->get_nb_elem_local(DIRECTION_I);
  const int nj = ref_splitting_->get_nb_elem_local(DIRECTION_J);
  const int nk = ref_splitting_->get_nb_elem_local(DIRECTION_K);

  ArrOfInt liste_composantes_connexes_dans_element;
  liste_composantes_connexes_dans_element.set_smart_resize(1);
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          assert(mesh_->ref_splitting() == ref_splitting_);
          // A present, elle est dans le splitting :
          const int elem = ref_splitting_->convert_ijk_cell_to_packed(i, j, k);
          // Pour chaque element, est-il traverse par une ou plusieurs interface ?
          // (seules les surfaces non-nulles sont comptees)
          const int nb_compo_traversantes =
            compute_list_compo_connex_in_element(elem, liste_composantes_connexes_dans_element);
          if (nb_compo_traversantes > max_authorized_nb_of_components_)
            {
              Cerr << "ComputeValParCompoInCell::ajouter_terme_source_interfaces. Trop de "
                   "compo connexes dans "
                   << "l'element " << elem << "[ " << i << " " << j << " " << k << " "
                   << " ]" << finl;
              Cerr << "Augmenter la taille max_authorized_nb_of_components_." << finl;
              assert(0);
            }

          nb_compo_traversante(i, j, k) = nb_compo_traversantes;
          // On boucle sur les composantes connexes trouvees dans l'element
          for (int i_compo = 0; i_compo < nb_compo_traversantes; i_compo++)
            {
              const int num_compo = liste_composantes_connexes_dans_element[i_compo];
              compos_traversantes[i_compo](i, j, k) = num_compo;

              double indic;
              calculer_indic_elem_pour_compo(num_compo, elem, indic);

              double surface;
              Vecteur3 normale;
              Vecteur3 bary_facettes_dans_elem;
              calculer_moyennes_interface_element_pour_compo(num_compo, elem, surface, normale, bary_facettes_dans_elem);

              indic_par_compo[i_compo](i, j, k) = indic;
              surface_par_compo[i_compo](i, j, k) = surface;
              for (int dir = 0; dir < 3; dir++)
                {
                  const int idx = i_compo * 3 + dir;
                  normale_par_compo[idx](i, j, k) = normale[dir];
                  bary_par_compo[idx](i, j, k) = bary_facettes_dans_elem[dir];
                }
            }
        }

#ifdef SMOOTHING_RHO
  if (smooth_density)
    {
      for (int icompo = 0; icompo < max_authorized_nb_of_components_; icompo++)
        {
          const int nx = indic_par_compo[icompo].ni();
          const int ny = indic_par_compo[icompo].nj();
          const int nz = indic_par_compo[icompo].nk();
          // Attention, rho n'est pas par compo, donc ca ne marche pas en
          // multi-bulles.
          for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
              for (int i = 0; i < nx; i++)
                {
                  //      double rho = smooth_rho(i,j,k);
                  double rho = rho_field(i, j, k);
                  indic_par_compo[icompo](i, j, k) = (rho - rho_v) / delta_rho;
                }
        }
    }
#endif

  // Mise a jour des espaces virtuels :
  nb_compo_traversante.echange_espace_virtuel(nb_compo_traversante.ghost());
  compos_traversantes.echange_espace_virtuel();
  surface_par_compo.echange_espace_virtuel();
  indic_par_compo.echange_espace_virtuel();
  normale_par_compo.echange_espace_virtuel();
  bary_par_compo.echange_espace_virtuel();
}


void ComputeValParCompoInCell::calculer_valeur_par_compo(
#ifdef SMOOTHING_RHO
  const double delta_rho,
#endif
  const double time,
  const int itstep,
  IJK_Field_int& nb_compo_trav,
  FixedVector<IJK_Field_int, max_authorized_nb_of_components_>& compos_trav,
  FixedVector<IJK_Field_double, 3*max_authorized_nb_of_components_>& normale_par_compo,
  FixedVector<IJK_Field_double, 3*max_authorized_nb_of_components_>& bary_par_compo,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& indicatrice_par_compo,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& surface_par_compo,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& courbure_par_compo
) // next()
{

  calculer_moy_par_compo(
#ifdef SMOOTHING_RHO
    const double delta_rho,
#endif
    nb_compo_trav,
    compos_trav,
    normale_par_compo,
    bary_par_compo,
    indicatrice_par_compo,
    surface_par_compo
  );

  // calcul courbure
  calculer_moy_field_sommet_par_compo(
    mesh_->get_update_courbure_sommets(),
    courbure_par_compo);
}

void ComputeValParCompoInCell::calculer_moy_field_sommet_par_compo(
  const ArrOfDouble& val_on_sommet,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& field_par_compo) const
{

  const ArrOfInt& compo_connexe = mesh_->compo_connexe_facettes();
  const ArrOfDouble& surface_facettes = mesh_->get_update_surface_facettes();
  const Intersections_Elem_Facettes& intersections = mesh_->intersections_elem_facettes();
  const IntTab& facettes = mesh_->facettes();

  for (int c = 0; c < max_authorized_nb_of_components_; c++)
    field_par_compo[c].data() = 0.;

  // Boucle sur les elements:
  const int ni = field_par_compo[0].ni();
  const int nj = field_par_compo[0].nj();
  const int nk = field_par_compo[0].nk();

  ArrOfInt liste_composantes_connexes_dans_element;
  liste_composantes_connexes_dans_element.set_smart_resize(1);
  liste_composantes_connexes_dans_element = 0;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          assert(mesh_->ref_splitting() == ref_splitting_);
          // A present, elle est dans le splitting :
          const int elem = ref_splitting_->convert_ijk_cell_to_packed(i, j, k);
          const int nb_compo_traversantes =
            compute_list_compo_connex_in_element(elem, liste_composantes_connexes_dans_element);
          // Pour chaque element, est-il traverse par une ou plusieurs interface ?
          assert(nb_compo_traversantes < max_authorized_nb_of_components_);

          // On boucle sur les composantes connexes trouvees dans l'element
          for (int i_compo = 0; i_compo < nb_compo_traversantes; i_compo++)
            {
              const int num_compo = liste_composantes_connexes_dans_element[i_compo];
              // Calcul de la moyenne de field dans l'element :
              double moy_field = 0.;
              int index = intersections.index_elem()[elem];
              assert(index >= 0); // Aucune facette dans cet element.

              double surface = 0.;
              double moy_field_fa7 = 0.;
              // Boucle sur les facettes qui traversent l'element elem :
              while (index >= 0)
                {
                  const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                  const int fa7 = data.numero_facette_;
                  const int icompo = compo_connexe[fa7];
                  moy_field_fa7 = 0.;
                  if (icompo == num_compo)
                    {
                      const double surface_facette = surface_facettes[fa7];
                      const double surf = data.fraction_surface_intersection_ * surface_facette;
                      // Les coordonnees du barycentre de la fraction de facette :
                      //Vecteur3 coord_barycentre_fraction(0., 0., 0.);
                      for (int isom = 0; isom < 3; isom++)
                        {
                          // numero du sommet dans le tableau sommets
                          const int num_som = facettes(fa7, isom);
                          const double val = val_on_sommet[num_som];
                          // Coordonnees barycentriques du centre de gravite de
                          // l'intersection par rapport aux trois sommets de la facette.
                          const double bary_som = data.barycentre_[isom];
                          // pour avoir la moyenne du champ au centre du triangle
                          moy_field_fa7 += bary_som * val;
                        }
                      surface += surf;
                      moy_field += moy_field_fa7 * surf;
                    }
                  index = data.index_facette_suivante_;
                }

              if (surface > 0.)
                moy_field *= 1. / surface;
              else
                moy_field = moy_field_fa7;
              field_par_compo[i_compo](i, j, k) = moy_field;
            }
        }

  // Mise a jour des espaces virtuels :
  field_par_compo.echange_espace_virtuel();
}

void ComputeValParCompoInCell::calculer_moy_field_fa7_par_compo(
  const ArrOfDouble& val_on_fa7,
  FixedVector<IJK_Field_double, max_authorized_nb_of_components_>& field_par_compo) const
{

  const ArrOfInt& compo_connexe = mesh_->compo_connexe_facettes();
  const ArrOfDouble& surface_facettes = mesh_->get_update_surface_facettes();
  const Intersections_Elem_Facettes& intersections = mesh_->intersections_elem_facettes();

  for (int c = 0; c < max_authorized_nb_of_components_; c++)
    field_par_compo[c].data() = 0.;

  // Boucle sur les elements:
  const int ni = field_par_compo[0].ni();
  const int nj = field_par_compo[0].nj();
  const int nk = field_par_compo[0].nk();

  ArrOfInt liste_composantes_connexes_dans_element;
  liste_composantes_connexes_dans_element.set_smart_resize(1);
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          assert(mesh_->ref_splitting() == ref_splitting_);
          // A present, elle est dans le splitting :
          const int elem = ref_splitting_->convert_ijk_cell_to_packed(i, j, k);
          const int nb_compo_traversantes =
            compute_list_compo_connex_in_element(elem, liste_composantes_connexes_dans_element);
          // Pour chaque element, est-il traverse par une ou plusieurs interface ?
          assert(nb_compo_traversantes > max_authorized_nb_of_components_);

          // On boucle sur les composantes connexes trouvees dans l'element
          for (int i_compo = 0; i_compo < nb_compo_traversantes; i_compo++)
            {
              const int num_compo = liste_composantes_connexes_dans_element[i_compo];
              // Calcul de la moyenne de field dans l'element :
              double moy_field = 0.;
              int index = intersections.index_elem()[elem];
              assert(index >= 0); // Aucune facette dans cet element.

              double surface = 0.;

              // Boucle sur les facettes qui traversent l'element elem :
              while (index >= 0)
                {
                  const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                  const int fa7 = data.numero_facette_;
                  const int icompo = compo_connexe[fa7];
                  if (icompo == num_compo)
                    {
                      const double surface_facette = surface_facettes[fa7];
                      const double surf = data.fraction_surface_intersection_ * surface_facette;
                      const double moy_field_fa7 = val_on_fa7[fa7];
                      surface += surf;
                      moy_field += moy_field_fa7 * surf;
                    }
                  index = data.index_facette_suivante_;
                }

              assert(surface > 0.);
              moy_field *= 1. / surface;
              field_par_compo[i_compo](i, j, k) = moy_field;
            }
        }
  field_par_compo.echange_espace_virtuel();
}
//
// Attention : elem ne doit pas etre un element virtuel :
int ComputeValParCompoInCell::compute_list_compo_connex_in_element(
  const int elem,
  ArrOfInt& liste_composantes_connexes_dans_element
) const
{
  const Intersections_Elem_Facettes& intersections = mesh_->intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  const ArrOfInt& compo_connexe_facettes = mesh_->compo_connexe_facettes();

  liste_composantes_connexes_dans_element.set_smart_resize(1);
  liste_composantes_connexes_dans_element.resize_array(0);

  // L'element, est-il traverse par une interface ?
  int index_next_facette = index_elem[elem];
  while (index_next_facette >= 0)
    {
      // Oui. Par quelles composantes connexes est-il traverse ?
      // On boucle sur les faces qui traversent l'interface, on marque
      // les composantes deja vues au fur et a mesure:
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index_next_facette);
      if (data.fraction_surface_intersection_ == 0.)
        {
#ifdef GB_VERBOSE
          Cerr << "compute_list_compo_connex_in_element : elem= " << elem
               << " Intersection issue du parcours : numero_facette_= " << data.numero_facette_
               << " numero_element_= " << data.numero_element_ << " fraction= " << data.fraction_surface_intersection_
               << finl;
#endif
          // Si on a une facette de fraction_surface_intersection nulle,
          // c'est qu'en realite elle n'intersecte pas.
          // On ne la traite pas.
        }
      else
        {
          int compo = compo_connexe_facettes[data.numero_facette_];
          // Cherche si cette composante est deja dans la liste
          int ii;
          const int nmax = liste_composantes_connexes_dans_element.size_array();
          for (ii = 0; ii < nmax; ii++)
            {
              if (liste_composantes_connexes_dans_element[ii] == compo)
                break; // sort de la boucle for, on aura ii < nmax
            }
          if (ii == nmax) // vrai si on n'a pas trouve compo dans
            // liste_composantes_connexes_dans_element
            liste_composantes_connexes_dans_element.append_array(compo);
        }
      index_next_facette = data.index_facette_suivante_;
    }

  const int nb_compo_traversantes = liste_composantes_connexes_dans_element.size_array();
  return nb_compo_traversantes;
}

// static int calculer_phi_et_indic_element_pour_compo(
int ComputeValParCompoInCell::calculer_indic_elem_pour_compo(
  const int icompo,
  const int elem,
  double& indic) const
{
  const Intersections_Elem_Facettes& intersections = mesh_->intersections_elem_facettes();
  const IntTab& facettes = mesh_->facettes();
  const DoubleTab& sommets = mesh_->sommets();
  const ArrOfInt& compo_connexe = mesh_->compo_connexe_facettes();

  indic = -1.; // valeur invalide.
  double somme_contrib = 0.;
  int index = intersections.index_elem()[elem];
  if (index < 0)
    return 0; // Aucune facette dans cet element.

  // Boucle sur les facettes qui traversent l'element elem :
  while (index >= 0)
    {
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
      const int fa7 = data.numero_facette_;
      const int num_compo = compo_connexe[fa7];
      if (icompo == num_compo)
        {
          // Calcul du potentiel au centre de l'intersection
          // Les coordonnees du barycentre de la fraction de facette :
          Vecteur3 coord_barycentre_fraction(0., 0., 0.);
          for (int isom = 0; isom < 3; isom++)
            {
              const int num_som = facettes(fa7, isom); // numero du sommet dans le tableau sommets
              // Coordonnees barycentriques du centre de gravite de l'intersection
              // par rapport aux trois sommets de la facette.
              const double bary_som = data.barycentre_[isom];
              for (int dir = 0; dir < 3; dir++)
                coord_barycentre_fraction[dir] += bary_som * sommets(num_som, dir);
            }
          // On ne somme la contribution que de la composante icompo :
          somme_contrib += data.contrib_volume_phase1_;
        }
      index = data.index_facette_suivante_;
    }

  // retour entre 0 et 1 :
  while (somme_contrib > 1.)
    somme_contrib -= 1.;
  while (somme_contrib < 0.)
    somme_contrib += 1.;

  if (somme_contrib == 0.)
    {
      Cerr << "calculer_indicatrice ne fait rien dans la boucle while... "
           "pourquoi? "
           << finl;
      // Process::exit();
    }

  indic = somme_contrib;
  return 1;
}
