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
// File:        Modele_turbulence_hyd_LES_DSGS_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_LES_DSGS_VDF.h>
#include <Champ_Fonc_P0_VDF.h>
#include <Schema_Temps_base.h>
#include <Champ_Face_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Equation_base.h>
#include <Domaine_VDF.h>
#include <TRUSTTrav.h>
#include <Debog.h>

Implemente_instanciable(Modele_turbulence_hyd_LES_DSGS_VDF, "Modele_turbulence_hyd_sous_maille_DSGS_VDF", Modele_turbulence_hyd_LES_Smago_VDF);

Sortie& Modele_turbulence_hyd_LES_DSGS_VDF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_LES_DSGS_VDF::readOn(Entree& s)
{
  return Modele_turbulence_hyd_LES_Smago_VDF::readOn(s);
}

void Modele_turbulence_hyd_LES_DSGS_VDF::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  Modele_turbulence_hyd_LES_Smago_VDF::associer(domaine_dis, domaine_Cl_dis);

  Cerr << "Discretisation de coeff_field" << finl;
  coeff_field.typer("Champ_Fonc_P0_VDF");
  Champ_Fonc_P0_VDF& coeff = ref_cast(Champ_Fonc_P0_VDF, coeff_field.valeur());
  coeff.associer_domaine_dis_base(le_dom_VF_.valeur());
  coeff.nommer("dynamic_coefficient");
  coeff.fixer_nb_comp(1);
  coeff.fixer_nb_valeurs_nodales(le_dom_VF_->nb_elem());
  coeff.fixer_unite("adim");
  coeff.changer_temps(0.);

  model_coeff.ref(coeff_field->valeurs());
  champs_compris_.ajoute_champ(coeff_field);
  // model_coeff.resize(le_dom_VDF->nb_elem_tot());
}

Champ_Fonc& Modele_turbulence_hyd_LES_DSGS_VDF::calculer_viscosite_turbulente()
{
  double temps = mon_equation_->inconnue().temps();
  DoubleTab& visco_turb = la_viscosite_turbulente_->valeurs();
  int nb_elem = ref_cast(Domaine_VDF, le_dom_VF_.valeur()).nb_elem();
  int nb_elem_tot = ref_cast(Domaine_VDF, le_dom_VF_.valeur()).nb_elem_tot();

  DoubleTrav Sij_test_scale(nb_elem_tot, dimension, dimension);
  DoubleTrav Sij_grid_scale(nb_elem_tot, dimension, dimension);
  DoubleTrav cell_cent_vel(nb_elem_tot, dimension);
  DoubleTrav filt_vel(nb_elem_tot, dimension);
  DoubleTrav Lij(nb_elem_tot, dimension, dimension);
  DoubleTrav Mij(nb_elem_tot, dimension, dimension);

  SMA_barre_.resize(nb_elem_tot);

  calculer_cell_cent_vel(cell_cent_vel);
  calculer_filter_field(cell_cent_vel, filt_vel);
  calculer_Lij(cell_cent_vel, filt_vel, Lij);
  calculer_Sij(cell_cent_vel, Sij_grid_scale);
  calculer_Sij(filt_vel, Sij_test_scale);
  calculer_Mij(Sij_grid_scale, Sij_test_scale, Mij);
  calculer_model_coefficient(Lij, Mij);

  calculer_S_barre();

  if (visco_turb.size() != nb_elem)
    {
      Cerr << "erreur dans la taille du DoubleTab valeurs de la viscosite" << finl;
      exit();
    }

  Debog::verifier("Modele_turbulence_hyd_LES_DSGS_VDF::calculer_viscosite_turbulente visco_turb 0", visco_turb);

  for (int elem = 0; elem < nb_elem; elem++)
    visco_turb[elem] = model_coeff[elem] * l_[elem] * l_[elem] * sqrt(SMA_barre_[elem]);

  Debog::verifier("Modele_turbulence_hyd_LES_DSGS_VDF::calculer_viscosite_turbulente visco_turb 1", visco_turb);

  coeff_field->changer_temps(temps);

  la_viscosite_turbulente_->changer_temps(temps);
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_cell_cent_vel(DoubleTab& cell_cent_vel)
{
  const DoubleTab& vitesse = mon_equation_->inconnue().valeurs();
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  int element_number;
  int num0, num1, num2, num3, num4 = -1, num5 = -1;

  // This is to calculate the cell centered velocity values

  for (element_number = 0; element_number < nb_elem_tot; element_number++)
    {
      num0 = elem_faces(element_number, 0);
      num1 = elem_faces(element_number, 1);
      num2 = elem_faces(element_number, 2);
      num3 = elem_faces(element_number, 3);
      if (dimension == 3)
        {
          num4 = elem_faces(element_number, 4);
          num5 = elem_faces(element_number, 5);
        }

      if (dimension == 2)
        {
          cell_cent_vel(element_number, 0) = 0.5 * (vitesse[num0] + vitesse[num2]);
          cell_cent_vel(element_number, 1) = 0.5 * (vitesse[num1] + vitesse[num3]);
        }
      else
        {
          cell_cent_vel(element_number, 0) = 0.5 * (vitesse[num0] + vitesse[num3]);
          cell_cent_vel(element_number, 1) = 0.5 * (vitesse[num1] + vitesse[num4]);
          cell_cent_vel(element_number, 2) = 0.5 * (vitesse[num2] + vitesse[num5]);
        }
    }
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_filter_field(const DoubleTab& in_vel, DoubleTab& out_vel)
{
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  int element_number;
  int num0, num1, num2, num3, num4, num5;
  int f0, f1, f2, f3, f4, f5;
  int i;

  DoubleTrav temp1(nb_elem_tot, dimension);
  DoubleTrav temp2(nb_elem_tot, dimension);

  // This is to calculate the filtered velocity field

  if (dimension == 2)
    {
      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f0 = elem_faces(element_number, 0);
          num0 = face_voisins(f0, 0);
          f2 = elem_faces(element_number, 2);
          num2 = face_voisins(f2, 1);
          if ((num0 == -1) || (num2 == -1))
            for (i = 0; i < dimension; i++)
              temp1(element_number, i) = in_vel(element_number, i);
          else
            {
              for (i = 0; i < dimension; i++)
                {
                  temp1(element_number, i) = 0.25 * (in_vel(num0, i) + 2.0 * in_vel(element_number, i) + in_vel(num2, i));
                }
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f1 = elem_faces(element_number, 1);
          num1 = face_voisins(f1, 0);
          f3 = elem_faces(element_number, 3);
          num3 = face_voisins(f3, 1);
          if ((num1 == -1) || (num3 == -1))
            for (i = 0; i < dimension; i++)
              out_vel(element_number, i) = temp1(element_number, i);
          else
            {
              for (i = 0; i < dimension; i++)
                {
                  out_vel(element_number, i) = 0.25 * (temp1(num1, i) + 2.0 * temp1(element_number, i) + temp1(num3, i));
                }
            }
        }
    }
  else
    {
      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f0 = elem_faces(element_number, 0);
          num0 = face_voisins(f0, 0);
          f3 = elem_faces(element_number, 3);
          num3 = face_voisins(f3, 1);
          if ((num0 == -1) || (num3 == -1))
            for (i = 0; i < dimension; i++)
              temp1(element_number, i) = in_vel(element_number, i);
          else
            {
              for (i = 0; i < dimension; i++)
                {
                  temp1(element_number, i) = 0.25 * (in_vel(num0, i) + 2.0 * in_vel(element_number, i) + in_vel(num3, i));
                }
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f1 = elem_faces(element_number, 1);
          num1 = face_voisins(f1, 0);
          f4 = elem_faces(element_number, 4);
          num4 = face_voisins(f4, 1);
          if ((num1 == -1) || (num4 == -1))
            for (i = 0; i < dimension; i++)
              temp2(element_number, i) = temp1(element_number, i);
          else
            {
              for (i = 0; i < dimension; i++)
                {
                  temp2(element_number, i) = 0.25 * (temp1(num1, i) + 2.0 * temp1(element_number, i) + temp1(num4, i));
                }
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f2 = elem_faces(element_number, 2);
          num2 = face_voisins(f2, 0);
          f5 = elem_faces(element_number, 5);
          num5 = face_voisins(f5, 1);
          if ((num2 == -1) || (num5 == -1))
            for (i = 0; i < dimension; i++)
              out_vel(element_number, i) = temp2(element_number, i);
          else
            {
              for (i = 0; i < dimension; i++)
                {
                  out_vel(element_number, i) = 0.25 * (temp2(num2, i) + 2.0 * temp2(element_number, i) + temp2(num5, i));
                }
            }
        }
    }
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_filter_tensor(DoubleTab& in_vel)
{
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  int element_number;
  int num0, num1, num2, num3, num4, num5;
  int f0, f1, f2, f3, f4, f5;
  int i, j;

  DoubleTrav temp1(nb_elem_tot, dimension, dimension);
  DoubleTrav temp2(nb_elem_tot, dimension, dimension);

  // This is to calculate the filtered tensor field

  if (dimension == 2)
    {
      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f0 = elem_faces(element_number, 0);
          num0 = face_voisins(f0, 0);
          f2 = elem_faces(element_number, 2);
          num2 = face_voisins(f2, 1);
          if ((num0 == -1) || (num2 == -1))
            for (i = 0; i < dimension; i++)
              for (j = 0; j < dimension; j++)
                temp1(element_number, i, j) = in_vel(element_number, i, j);
          else
            {
              for (i = 0; i < dimension; i++)
                for (j = 0; j < dimension; j++)
                  {
                    temp1(element_number, i, j) = 0.25 * (in_vel(num0, i, j) + 2.0 * in_vel(element_number, i, j) + in_vel(num2, i, j));
                  }
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f1 = elem_faces(element_number, 1);
          num1 = face_voisins(f1, 0);
          f3 = elem_faces(element_number, 3);
          num3 = face_voisins(f3, 1);
          if ((num1 == -1) || (num3 == -1))
            for (i = 0; i < dimension; i++)
              for (j = 0; j < dimension; j++)
                in_vel(element_number, i, j) = temp1(element_number, i, j);
          else
            {
              for (i = 0; i < dimension; i++)
                for (j = 0; j < dimension; j++)
                  {
                    in_vel(element_number, i, j) = 0.25 * (temp1(num1, i, j) + 2.0 * temp1(element_number, i, j) + temp1(num3, i, j));
                  }
            }
        }
    }
  else
    {
      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f0 = elem_faces(element_number, 0);
          num0 = face_voisins(f0, 0);
          f3 = elem_faces(element_number, 3);
          num3 = face_voisins(f3, 1);
          if ((num0 == -1) || (num3 == -1))
            for (i = 0; i < dimension; i++)
              for (j = 0; j < dimension; j++)
                temp1(element_number, i, j) = in_vel(element_number, i, j);
          else
            {
              for (i = 0; i < dimension; i++)
                for (j = 0; j < dimension; j++)
                  {
                    temp1(element_number, i, j) = 0.25 * (in_vel(num0, i, j) + 2.0 * in_vel(element_number, i, j) + in_vel(num3, i, j));
                  }
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f1 = elem_faces(element_number, 1);
          num1 = face_voisins(f1, 0);
          f4 = elem_faces(element_number, 4);
          num4 = face_voisins(f4, 1);
          if ((num1 == -1) || (num4 == -1))
            for (i = 0; i < dimension; i++)
              for (j = 0; j < dimension; j++)
                temp2(element_number, i, j) = temp1(element_number, i, j);
          else
            {
              for (i = 0; i < dimension; i++)
                for (j = 0; j < dimension; j++)
                  {
                    temp2(element_number, i, j) = 0.25 * (temp1(num1, i, j) + 2.0 * temp1(element_number, i, j) + temp1(num4, i, j));
                  }
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f2 = elem_faces(element_number, 2);
          num2 = face_voisins(f2, 0);
          f5 = elem_faces(element_number, 5);
          num5 = face_voisins(f5, 1);
          if ((num2 == -1) || (num5 == -1))
            for (i = 0; i < dimension; i++)
              for (j = 0; j < dimension; j++)
                in_vel(element_number, i, j) = temp2(element_number, i, j);
          else
            {
              for (i = 0; i < dimension; i++)
                for (j = 0; j < dimension; j++)
                  {
                    in_vel(element_number, i, j) = 0.25 * (temp2(num2, i, j) + 2.0 * temp2(element_number, i, j) + temp2(num5, i, j));
                  }
            }
        }
    }
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_filter_coeff(DoubleVect& in_vel)
{
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  int element_number;
  int num0, num1, num2, num3, num4, num5;
  int f0, f1, f2, f3, f4, f5;

  DoubleTrav temp1(nb_elem_tot);
  DoubleTrav temp2(nb_elem_tot);

  // This is to calculate the filtered coefficient field

  if (dimension == 2)
    {
      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f0 = elem_faces(element_number, 0);
          num0 = face_voisins(f0, 0);
          f2 = elem_faces(element_number, 2);
          num2 = face_voisins(f2, 1);
          if ((num0 == -1) || (num2 == -1))
            temp1(element_number) = 3.0 * in_vel(element_number);
          else
            {
              temp1(element_number) = in_vel(num0) + in_vel(element_number) + in_vel(num2);
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f1 = elem_faces(element_number, 1);
          num1 = face_voisins(f1, 0);
          f3 = elem_faces(element_number, 3);
          num3 = face_voisins(f3, 1);
          if ((num1 == -1) || (num3 == -1))
            in_vel(element_number) = temp1(element_number) / 3.0;
          else
            {
              in_vel(element_number) = (temp1(num1) + temp1(element_number) + temp1(num3)) / 9.0;
            }
        }
    }
  else
    {
      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f0 = elem_faces(element_number, 0);
          num0 = face_voisins(f0, 0);
          f3 = elem_faces(element_number, 3);
          num3 = face_voisins(f3, 1);
          if ((num0 == -1) || (num3 == -1))
            temp1(element_number) = 3.0 * in_vel(element_number);
          else
            {
              temp1(element_number) = in_vel(num0) + in_vel(element_number) + in_vel(num3);
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f1 = elem_faces(element_number, 1);
          num1 = face_voisins(f1, 0);
          f4 = elem_faces(element_number, 4);
          num4 = face_voisins(f4, 1);
          if ((num1 == -1) || (num4 == -1))
            temp2(element_number) = 3.0 * temp1(element_number);
          else
            {
              temp2(element_number) = temp1(num1) + temp1(element_number) + temp1(num4);
            }
        }

      for (element_number = 0; element_number < nb_elem_tot; element_number++)
        {
          f2 = elem_faces(element_number, 2);
          num2 = face_voisins(f2, 0);
          f5 = elem_faces(element_number, 5);
          num5 = face_voisins(f5, 1);
          if ((num2 == -1) || (num5 == -1))
            in_vel(element_number) = temp2(element_number) / 9.0;
          else
            {
              in_vel(element_number) = (temp2(num2) + temp2(element_number) + temp2(num5)) / 27.0;
            }
        }
    }
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_Lij(const DoubleTab& cell_cent_vel, const DoubleTab& filt_vel, DoubleTab& Lij)
{
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot();
  int element_number;

  DoubleTrav uij_filt(nb_elem_tot, dimension, dimension);

  // This is to calculate the Lij term for the C coefficient

  for (element_number = 0; element_number < nb_elem_tot; element_number++)
    {
      for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
          uij_filt(element_number, i, j) = cell_cent_vel(element_number, i) * cell_cent_vel(element_number, j);
    }

  calculer_filter_tensor(uij_filt);

  for (element_number = 0; element_number < nb_elem_tot; element_number++)
    {
      for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
          Lij(element_number, i, j) = uij_filt(element_number, i, j) - filt_vel(element_number, i) * filt_vel(element_number, j);
    }
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_Mij(const DoubleTab& Sij_grid_scale, const DoubleTab& Sij_test_scale, DoubleTab& Mij)
{
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot();
  int element_number;

  DoubleTrav sij_filt(nb_elem_tot, dimension, dimension);
  double temp, coeff = 0;
  const double alpha = sqrt(6.0);
  int i, j;

  sij_filt = Sij_grid_scale;

  // This is to calculate the Mij term for the C coefficient

  for (element_number = 0; element_number < nb_elem_tot; element_number++)
    {
      temp = 0.;
      for (i = 0; i < dimension; i++)
        for (j = 0; j < dimension; j++)
          temp += Sij_grid_scale(element_number, i, j) * Sij_grid_scale(element_number, i, j);

      coeff = sqrt(2.0 * temp);
    }
  //GF
  Cerr << "dans " << que_suis_je() << " le calcul de coeff est debile!!!!" << finl;
  Cerr << " On fait planter pour corriger plus tard " << finl;
  // en sortie de la boucle coeff=sqrt(2*tmp) tmp calculer qu avec le dernier
  // elt
  assert(0);
  exit();
  sij_filt *= coeff;

  calculer_filter_tensor(sij_filt);

  for (element_number = 0; element_number < nb_elem_tot; element_number++)
    {
      temp = 0.;
      for (i = 0; i < dimension; i++)
        for (j = 0; j < dimension; j++)
          temp += Sij_test_scale(element_number, i, j) * Sij_test_scale(element_number, i, j);

      coeff = sqrt(2.0 * temp);
      for (i = 0; i < dimension; i++)
        for (j = 0; j < dimension; j++)
          Mij(element_number, i, j) = l_[element_number] * l_[element_number] * (alpha * alpha * coeff * Sij_test_scale(element_number, i, j) - sij_filt(element_number, i, j));
    }
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_model_coefficient(const DoubleTab& Lij, const DoubleTab& Mij)
{
  int nb_elem_tot = ref_cast(Domaine_VDF, le_dom_VF_.valeur()).nb_elem_tot();
  double temp1;
  double temp2;

  // Evaluate the dynamic model coeficient C
  int elem;
  for (elem = 0; elem < nb_elem_tot; elem++)
    {
      temp1 = 0.;
      temp2 = 0.;
      for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
          {
            temp1 += Lij(elem, i, j) * Mij(elem, i, j);
            temp2 += Mij(elem, i, j) * Mij(elem, i, j);
          }
      if (std::fabs(temp2) < 1.e-12)
        model_coeff[elem] = 0.;
      else
        model_coeff[elem] = -0.5 * temp1 / temp2;
    }

  calculer_filter_coeff(model_coeff);
  //  double coeff_moy = 0.;
  for (elem = 0; elem < nb_elem_tot; elem++)
    {
      if (model_coeff[elem] < 0.0)
        model_coeff[elem] = 0.0;
      if (model_coeff[elem] > 0.5)
        model_coeff[elem] = 0.5;
    }
}

void Modele_turbulence_hyd_LES_DSGS_VDF::calculer_Sij(const DoubleTab& in_vel, DoubleTab& out_vel)
{
  Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, mon_equation_->inconnue());
  const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, le_dom_VF_.valeur());

  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot();
  //Nouvaeu calcul de Sij
  DoubleTab duidxj(nb_elem_tot, dimension, dimension);
  vit.calcul_duidxj(in_vel, duidxj);

  for (int elem = 0; elem < nb_elem_tot; elem++)
    for (int i = 0; i < dimension; i++)
      for (int j = 0; j < dimension; j++)
        out_vel(elem, i, j) = 0.5 * (duidxj(elem, i, j) + duidxj(elem, j, i));
}
