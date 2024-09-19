/****************************************************************************
* Copyright (c) 2024, CEA
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
// File:        Modele_turbulence_hyd_LES_Fst_sel_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_LES_Fst_sel_VEF.h>
#include <VEF_discretisation.h>
#include <Domaine_VEF.h>
#include <Champ_P1NC.h>
#include <math.h>

Implemente_instanciable_sans_constructeur(Modele_turbulence_hyd_LES_Fst_sel_VEF, "Modele_turbulence_hyd_sous_maille_fst_selectif_VEF", Modele_turbulence_hyd_LES_Fst_VEF);

Modele_turbulence_hyd_LES_Fst_sel_VEF::Modele_turbulence_hyd_LES_Fst_sel_VEF()
{
  C1_ = 0.777 * 0.18247 * 0.18247 * 1.356;  // C_Fst_sel=C_Fst*1.356 car modele selectif
}

Sortie& Modele_turbulence_hyd_LES_Fst_sel_VEF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_LES_Fst_sel_VEF::readOn(Entree& s)
{
  return Modele_turbulence_hyd_LES_Fst_VEF::readOn(s);
}

void Modele_turbulence_hyd_LES_Fst_sel_VEF::discretiser()
{
  if (dimension != 3)
    {
      Cerr << " The model sous_maille_fst_selectif has sense only for dimension 3." << finl;
      exit();

    }
  Modele_turbulence_hyd_LES_base::discretiser();
  const VEF_discretisation& dis = ref_cast(VEF_discretisation, mon_equation_->discretisation());
  dis.vorticite(mon_equation_->domaine_dis(), mon_equation_->inconnue(), la_vorticite_);
}

int Modele_turbulence_hyd_LES_Fst_sel_VEF::a_pour_Champ_Fonc(const Motcle& mot,
                                                             REF(Champ_base) &ch_ref) const
{
  Motcles les_motcles(3);
  {
    les_motcles[0] = "viscosite_turbulente";
    les_motcles[1] = "k";
    les_motcles[2] = "vorticite";

  }
  int rang = les_motcles.search(mot);
  switch(rang)
    {
    case 0:
      {
        ch_ref = la_viscosite_turbulente_.valeur();
        return 1;
      }
    case 1:
      {
        ch_ref = energie_cinetique_turb_.valeur();
        return 1;
      }
    case 2:
      {
        ch_ref = la_vorticite_.valeur();
        return 1;
      }
    default:
      return 0;
    }
}

void Modele_turbulence_hyd_LES_Fst_sel_VEF::calculer_racine()
{
  //  Modele_turbulence_hyd_LES_VEF::calculer_fonction_structure(FS,delta_c);
  Modele_turbulence_hyd_LES_Fst_VEF::calculer_racine();
  cutoff();
}

// Fonction qui permet d'appliquer un filtre sur la fonction de structure
// La fonction de structure d'un element est mise a zero si il existe une
// deviation inferieure a N degres entre son vecteur vorticite et le
// vecteur moyen des vorticites des 6 elements les plus proches
// l angle de coupure varie en fonction du delta c (regression log)

void Modele_turbulence_hyd_LES_Fst_sel_VEF::cutoff()
{
  double Sin2Angl;
  const Champ_P1NC& vitesse = ref_cast(Champ_P1NC, mon_equation_->inconnue());
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_VF_.valeur());
  const int nb_elem = domaine_VEF.nb_elem();
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  //  const Domaine& domaine = domaine_VEF.domaine();
  //  int nfac = domaine.nb_faces_elem();
  //  int nfac = 4; // en 3D 4 faces!!!
  DoubleTab& vorticite = la_vorticite_->valeurs();
  const DoubleTab& xp = domaine_VEF.xp();

  la_vorticite_->mettre_a_jour(vitesse.temps());
  vorticite.echange_espace_virtuel();

  //  int el0,el1,el2,el3;
  double norme, norme_moyen, prod, angle, delta;
  DoubleVect vorti_moyen(3);
  IntVect elem_nn(4);
  IntVect elem_autour(4);
  DoubleVect dist(4);
  int nb_elem_nn, compteur, elem;
  double d;

  for (int num_elem = 0; num_elem < nb_elem; num_elem++)
    {
      delta = l_[num_elem];
      angle = 0.494 * pow(delta, 0.4075);
      Sin2Angl = sin(angle);
      Sin2Angl *= Sin2Angl;

      elem_autour[0] = face_voisins(elem_faces(num_elem, 0), 0);
      if (elem_autour[0] == num_elem)
        elem_autour[0] = face_voisins(elem_faces(num_elem, 0), 1);
      elem_autour[1] = face_voisins(elem_faces(num_elem, 1), 0);
      if (elem_autour[1] == num_elem)
        elem_autour[1] = face_voisins(elem_faces(num_elem, 1), 1);
      elem_autour[2] = face_voisins(elem_faces(num_elem, 2), 0);
      if (elem_autour[2] == num_elem)
        elem_autour[2] = face_voisins(elem_faces(num_elem, 2), 1);
      elem_autour[3] = face_voisins(elem_faces(num_elem, 3), 0);
      if (elem_autour[3] == num_elem)
        elem_autour[3] = face_voisins(elem_faces(num_elem, 3), 1);

      //      double d0,d1,d2,d3,d;
      double x, y, z;
      x = xp(num_elem, 0);
      y = xp(num_elem, 1);
      z = xp(num_elem, 2);

      nb_elem_nn = -1;
      elem_nn = 0;
      compteur = 0;
      d = 0.;

      while (compteur < 4)
        {
          if (elem_autour[compteur] != -1)
            {
              nb_elem_nn++;
              elem_nn[nb_elem_nn] = elem_autour[compteur];
            }
          compteur++;
        }

      for (elem = 0; elem < nb_elem_nn; elem++)
        {
          dist[elem] = (x - xp(elem_nn[elem], 0)) * (x - xp(elem_nn[elem], 0)) + (y - xp(elem_nn[elem], 1)) * (y - xp(elem_nn[elem], 1)) + (z - xp(elem_nn[elem], 2)) * (z - xp(elem_nn[elem], 2));
          dist[elem] = 1. / dist[elem];
          d += dist[elem];
        }

      vorti_moyen = 0.;
      for (elem = 0; elem < nb_elem_nn; elem++)
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) += dist[elem] * vorticite(elem_nn[elem], k) / d;
        }

      // On rend nul le vecteur vorti_moyen(k) ce qui provoquera la mise a zero de FS
      if (nb_elem_nn == -1)  // Cas d'un element coin ; on met FS a zero
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = 0;
        }

      // Calcul du produit vectoriel entre la vorticite dans l'element
      // et le vecteur des vorticites des elements voisins

      norme = 0;
      int k;
      for (k = 0; k < 3; k++)
        norme += carre(vorticite(num_elem, k));

      norme_moyen = 0;
      for (k = 0; k < 3; k++)
        norme_moyen += carre(vorti_moyen(k));

      if ((norme > 1.e-10) && (norme_moyen > 1.e-10))
        {
          prod = carre(vorti_moyen(1) * vorticite(num_elem, 2) - vorti_moyen(2) * vorticite(num_elem, 1)) + carre(vorti_moyen(2) * vorticite(num_elem, 0) - vorti_moyen(0) * vorticite(num_elem, 2))
                 + carre(vorti_moyen(0) * vorticite(num_elem, 1) - vorti_moyen(1) * vorticite(num_elem, 0));
          prod /= (norme * norme_moyen);

          if (prod <= Sin2Angl)
            Racine_(num_elem) = 0;
        }
      else
        // bruit numerique ou element de coin
        Racine_(num_elem) = 0;
      Racine_.echange_espace_virtuel();

    }
}
