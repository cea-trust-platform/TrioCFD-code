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
// File:        Modele_turbulence_hyd_LES_selectif_mod_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_LES_selectif_mod_VDF.h>
#include <VDF_discretisation.h>
#include <Champ_Face_VDF.h>
#include <Domaine_VDF.h>
#include <Equation_base.h>
#include <TRUSTTrav.h>
#include <Param.h>
#include <math.h>

#define K_DEF 0.415

Implemente_instanciable_sans_constructeur(Modele_turbulence_hyd_LES_selectif_mod_VDF, "Modele_turbulence_hyd_sous_maille_selectif_mod_VDF", Modele_turbulence_hyd_LES_VDF);

Sortie& Modele_turbulence_hyd_LES_selectif_mod_VDF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_LES_selectif_mod_VDF::readOn(Entree& s)
{
  return Modele_turbulence_hyd_LES_VDF::readOn(s);
}

void Modele_turbulence_hyd_LES_selectif_mod_VDF::set_param(Param& param)
{
  Modele_turbulence_hyd_LES_VDF::set_param(param);
  param.ajouter_non_std("Canal", (this));
  param.ajouter_non_std("THI", (this));
}

int Modele_turbulence_hyd_LES_selectif_mod_VDF::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Canal")
    {
      is >> demi_h_;
      is >> dir_par_;
      canal_ = 1;
      Cerr << "The half height of the chanel is :" << demi_h_ << " m" << finl;
      Cerr << "Orientation of the walls is : " << dir_par_ << finl;
      return 1;
    }
  else if (mot == "THI")
    {
      is >> ki_;
      is >> kc_;
      thi_ = 1;
      Cerr << "At initial time, the spectrum peak is located at : " << ki_ << " m-1" << finl;
      Cerr << "The spectrum cut-off is located at : " << kc_ << " m-1" << finl;
      return 1;
    }
  else
    return Modele_turbulence_hyd_LES_VDF::lire_motcle_non_standard(mot, is);
}


void Modele_turbulence_hyd_LES_selectif_mod_VDF::discretiser()
{
  Modele_turbulence_hyd_LES_base::discretiser();
  const VDF_discretisation& dis = ref_cast(VDF_discretisation, mon_equation_->discretisation());
  dis.vorticite(mon_equation_->domaine_dis(), mon_equation_->inconnue(), la_vorticite_);
}

int Modele_turbulence_hyd_LES_selectif_mod_VDF::a_pour_Champ_Fonc(const Motcle& mot,
                                                                  REF(Champ_base) &ch_ref) const
{
  Motcles les_motcles(7);
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

void Modele_turbulence_hyd_LES_selectif_mod_VDF::calculer_fonction_structure()
{

  Modele_turbulence_hyd_LES_VDF::calculer_fonction_structure();
  cutoff();
}

// Fonction qui permet d'appliquer un filtre sur la fonction de structure
// La fonction de structure d'un element est mise a zero si il existe une
// deviation inferieure a N degres entre son vecteur vorticite et le
// vecteur moyen des vorticites des 6 elements les plus proches
// l angle de coupure varie en fonction du delta c (2 expressions possibles :
// lineaire ou log)

void Modele_turbulence_hyd_LES_selectif_mod_VDF::cutoff()
{
  double Sin2Angl;
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF, mon_equation_->inconnue());
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  const DoubleTab& xp = domaine_VDF.xp();
  int nb_poly = domaine_VDF.domaine().nb_elem();
  DoubleTab& vorticite = la_vorticite_->valeurs();

  la_vorticite_->mettre_a_jour(vitesse.temps());
  vorticite.echange_espace_virtuel();

  int elx0, elx1, ely0, ely1, elz0, elz1;
  double norme, norme_moyen, prod, angle;
  DoubleTrav vorti_moyen(3);
  int i;
  int num_elem;

  double lm, rapport;
  double dy, y_elem;
  static double kappa = K_DEF;
  static double a = 23. * M_PI / 180.;
  static double b = -0.4;
  static double borne_y;
  if (canal_ != 0)
    borne_y = 0.2 * demi_h_;

  if ((thi_ == 0) && (canal_ == 0))
    {
      Cerr << "Problem with the input data of the sous_maille_selectif_mod model." << finl;
      Cerr << "A necessary parameter : THI or Canal has not been specified." << finl;
      exit();
    }

  if ((nb_points_ == 4) && (dir_par_ != dir3_))
    {
      Cerr << "WARNING!!! WARNING!!! Modele_turbulence_hyd_LES_selectif_mod_VDF" << finl;
      Cerr << "A 4 points formulation is used." << finl;
      Cerr << "Your planes are orthogonal to the direction : " << dir3_ << finl;
      Cerr << "while the given direction for the walls is :" << dir_par_ << finl;
      Cerr << "Check that you really wish these directions to be different ..." << finl;
      Cerr << "WARNING!!! WARNING!!! dans Modele_turbulence_hyd_LES_selectif_mod_VDF" << finl;
    }

  for (num_elem = 0; num_elem < nb_poly; num_elem++)
    {

      if (thi_ != 0)
        {
          rapport = 2. * kc_ / ki_;
        }
      else
        {
          // y_elem : distance entre le centre de l element et la paroi!!
          y_elem = xp(num_elem, dir_par_);
          if (y_elem >= demi_h_)
            y_elem = 2. * demi_h_ - y_elem;

          // lm : longueur de melange dans le canal
          if (y_elem < borne_y)
            lm = kappa * y_elem;
          else
            lm = 0.2 * kappa * demi_h_;

          // dy : pas de maillage local en y (largeur de filtre local)
          dy = domaine_VDF.dim_elem(num_elem, dir_par_);

          // rapport pour la dependance de l angle
          rapport = lm / (2. * dy);
        }

      if (rapport <= 1.)
        angle = a;
      else if ((rapport >= 1.) && (rapport < 10.))
        angle = a * pow(rapport, b);
      else
        angle = 9. * M_PI / 180.;

      Sin2Angl = sin(angle);
      Sin2Angl *= Sin2Angl;

      elx0 = face_voisins(elem_faces(num_elem, 0), 0);
      elx1 = face_voisins(elem_faces(num_elem, 3), 1);
      ely0 = face_voisins(elem_faces(num_elem, 1), 0);
      ely1 = face_voisins(elem_faces(num_elem, 4), 1);
      elz0 = face_voisins(elem_faces(num_elem, 2), 0);
      elz1 = face_voisins(elem_faces(num_elem, 5), 1);
      //dist_elem_period(int n1, int n2, int k)
      double dx0 = 0., dx1 = 0., dy0 = 0., dy1 = 0., dz0 = 0., dz1 = 0.;

      for (i = 0; i < 3; i++)
        {
          if (elx0 != -1)
            dx0 += domaine_VDF.dist_elem_period(num_elem, elx0, i) * domaine_VDF.dist_elem_period(num_elem, elx0, i);
          if (elx1 != -1)
            dx1 += domaine_VDF.dist_elem_period(num_elem, elx1, i) * domaine_VDF.dist_elem_period(num_elem, elx1, i);
          if (ely0 != -1)
            dy0 += domaine_VDF.dist_elem_period(num_elem, ely0, i) * domaine_VDF.dist_elem_period(num_elem, ely0, i);
          if (ely1 != -1)
            dy1 += domaine_VDF.dist_elem_period(num_elem, ely1, i) * domaine_VDF.dist_elem_period(num_elem, ely1, i);
          if (elz0 != -1)
            dz0 += domaine_VDF.dist_elem_period(num_elem, elz0, i) * domaine_VDF.dist_elem_period(num_elem, elz0, i);
          if (elz1 != -1)
            dz1 += domaine_VDF.dist_elem_period(num_elem, elz1, i) * domaine_VDF.dist_elem_period(num_elem, elz1, i);
        }

      if (std::fabs(dx0) > DMINFLOAT)
        dx0 = 1. / sqrt(dx0);
      if (std::fabs(dx1) > DMINFLOAT)
        dx1 = 1. / sqrt(dx1);
      if (std::fabs(dy0) > DMINFLOAT)
        dy0 = 1. / sqrt(dy0);
      if (std::fabs(dy1) > DMINFLOAT)
        dy1 = 1. / sqrt(dy1);
      if (std::fabs(dz0) > DMINFLOAT)
        dz0 = 1. / sqrt(dz0);
      if (std::fabs(dz1) > DMINFLOAT)
        dz1 = 1. / sqrt(dz1);

      if ((elx0 != -1) && (elx1 != -1) && (ely0 != -1) && (ely1 != -1) && (elz0 != -1) && (elz1 != -1))
        // Cas d'un element interne
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = (dx0 * vorticite(elx0, k) + dx1 * vorticite(elx1, k) + dy0 * vorticite(ely0, k) + dy1 * vorticite(ely1, k) + dz0 * vorticite(elz0, k) + dz1 * vorticite(elz1, k))
                             / (dx0 + dx1 + dy0 + dy1 + dz0 + dz1);
        }
      else if ((elx0 != -1) && (elx1 != -1) && (ely0 != -1) && (ely1 != -1))
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = (dx0 * vorticite(elx0, k) + dx1 * vorticite(elx1, k) + dy0 * vorticite(ely0, k) + dy1 * vorticite(ely1, k)) / (dx0 + dx1 + dy0 + dy1);
        }
      else if ((elx0 != -1) && (elx1 != -1) && (elz0 != -1) && (elz1 != -1))
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = (dx0 * vorticite(elx0, k) + dx1 * vorticite(elx1, k) + dz0 * vorticite(elz0, k) + dz1 * vorticite(elz1, k)) / (dx0 + dx1 + dz0 + dz1);
        }
      else if ((ely0 != -1) && (ely1 != -1) && (elz0 != -1) && (elz1 != -1))
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = (dy0 * vorticite(ely0, k) + dy1 * vorticite(ely1, k) + dz0 * vorticite(elz0, k) + dz1 * vorticite(elz1, k)) / (dy0 + dy1 + dz0 + dz1);
        }
      else if ((elx0 != -1) && (elx1 != -1))
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = (dx0 * vorticite(elx0, k) + dx1 * vorticite(elx1, k)) / (dx0 + dx1);
        }
      else if ((ely0 != -1) && (ely1 != -1))
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = (dy0 * vorticite(ely0, k) + dy1 * vorticite(ely1, k)) / (dy0 + dy1);
        }
      else if ((elz0 != -1) && (elz1 != -1))
        {
          for (int k = 0; k < 3; k++)
            vorti_moyen(k) = (dz0 * vorticite(elz0, k) + dz1 * vorticite(elz1, k)) / (dz0 + dz1);
        }
      else  // Cas d'un element coin ; on met FS a zero
        // On rend nul le vecteur vorti_moyen(k) ce qui provoquera la mise a zero de FS
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
            F2_(num_elem) = 0;
        }
      else
        // bruit numerique ou element de coin
        F2_(num_elem) = 0;

    }
  F2_.echange_espace_virtuel();
}
