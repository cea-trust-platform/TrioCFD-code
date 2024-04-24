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
// File:        Modele_turbulence_hyd_LES_Smago_filtre_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_LES_Smago_filtre_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Domaine_VEF.h>
#include <Champ_P1NC.h>

Implemente_instanciable(Modele_turbulence_hyd_LES_Smago_filtre_VEF, "Modele_turbulence_hyd_sous_maille_Smago_filtre_VEF", Modele_turbulence_hyd_LES_Smago_VEF);

Sortie& Modele_turbulence_hyd_LES_Smago_filtre_VEF::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Modele_turbulence_hyd_LES_Smago_filtre_VEF::readOn(Entree& s) { return Modele_turbulence_hyd_LES_Smago_VEF::readOn(s); }

void Modele_turbulence_hyd_LES_Smago_filtre_VEF::calculer_S_barre()
{
  const DoubleTab& la_vitesse = mon_equation_->inconnue().valeurs();
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF, le_dom_Cl_.valeur());
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_VF_.valeur());
  const int nb_elem = domaine_VEF.nb_elem();

  const DoubleVect& vol = domaine_VEF.volumes();
  const Domaine& domaine = domaine_VEF.domaine();
  int nb_faces_elem = domaine.nb_faces_elem();

  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const int nb_face = domaine_VEF.nb_faces();

  int i, elem;
  int fac = 0;

  //////////////////////////////
  //Filtrage du champ de vitesse
  //////////////////////////////
  DoubleTab vitesse(la_vitesse);

  for (; fac < nb_face; fac++)
    {
      int num1;
      num1 = face_voisins(fac, 0);
      int num2;
      num2 = face_voisins(fac, 1);

      int fac1, fac2, facel;

      vitesse = 0.;

      for (facel = 0; facel < nb_faces_elem; facel++)
        {
          fac1 = elem_faces(num1, facel);
          //Correction pour avoir le champ de vitesse par face
          for (i = 0; i < dimension; i++)
            vitesse(fac, i) += la_vitesse(fac1, i) / double(2 * nb_faces_elem);
          fac2 = elem_faces(num2, facel);
          for (i = 0; i < dimension; i++)
            vitesse(fac, i) += la_vitesse(fac2, i) / double(2 * nb_faces_elem);
        }

    }

  Champ_P1NC::calcul_S_barre(vitesse, SMA_barre_, domaine_Cl_VEF);

  // On recalcule la longueur caracteristique de l'element

  for (elem = 0; elem < nb_elem; elem++)
    {
      double voltot = vol(elem);
      double eldif = 1.;
      int num1, facel;

      for (facel = 0; facel < nb_faces_elem; facel++)
        {
          fac = elem_faces(elem, facel);
          num1 = face_voisins(fac, 0);
          if (num1 == elem)
            num1 = face_voisins(fac, 1);

          if (num1 != -1)
            {
              voltot += vol(num1);
              eldif += 1.;
            }
        }

      voltot /= eldif;
      l_(elem) = 2.0 * pow(6. * voltot, 1. / double(dimension));
    }
}
