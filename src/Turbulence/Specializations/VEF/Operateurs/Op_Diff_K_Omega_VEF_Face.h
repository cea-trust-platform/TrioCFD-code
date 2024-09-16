/****************************************************************************
* Copyright (c) 2023, CEA
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
/*! @brief class Op_Diff_K_Omega_VEF_Face Cette classe represente l'operateur de diffusion turbulente
 *
 *    La discretisation est VEF
 *   Les methodes pour l'implicite sont codees.
 *
 */
#ifndef Op_Diff_K_Omega_VEF_Face_included
#define Op_Diff_K_Omega_VEF_Face_included

#include <Op_Diff_K_Omega_VEF_base.h>
#include <Op_VEF_Face.h>
#include <Domaine_VEF.h>
#include <TRUST_Ref.h>
#include <Champ_Fonc.h>
#include <Champ_Inc.h>



class Domaine_Cl_VEF;
class Champ_P1NC;
class Champ_Don_base;

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Op_Diff_K_Omega_VEF_Face
//
//////////////////////////////////////////////////////////////////////////////

class Op_Diff_K_Omega_VEF_Face : public Op_Diff_K_Omega_VEF_base, public Op_VEF_Face
{

  Declare_instanciable(Op_Diff_K_Omega_VEF_Face);

public:

  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&,
                const Champ_Inc&) override;
  void associer_diffusivite_turbulente() override;
  const Champ_Fonc& diffusivite_turbulente() const;
  DoubleTab& ajouter(const DoubleTab&, DoubleTab&) const override;
  DoubleTab& calculer(const DoubleTab&, DoubleTab&) const override;
  inline double viscA(int, int, int, double) const;
  void calc_visc(ArrOfDouble& diffu_tot, const Domaine_VEF& le_dom,
                 int num_face, int num2, int dimension, int num_elem,
                 double diffu_turb, const DoubleTab& diffu, int is_mu_unif,
                 const ArrOfDouble& inv_Prdt) const;

  double blender(double const val1, double const val2, int const face) const;
  // Methodes pour l implicite.

  inline void dimensionner(Matrice_Morse&) const override;
  void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;
  inline void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  inline void contribuer_au_second_membre(DoubleTab&) const override;
  void contribue_au_second_membre(DoubleTab&) const;
  void ajouter_contribution(const DoubleTab&, Matrice_Morse&) const;

protected :
  REF(Domaine_VEF) le_dom_vef;
  REF(Domaine_Cl_VEF) la_zcl_vef;
  REF(Champ_P1NC) inconnue_;
};

// ATTENTION le diffu intervenant dans les fonctions n'est que LOCAL (on appelle d_mu apres)
// Fonction utile visc
// mu <Si, Sj> / |K|
inline double Op_Diff_K_Omega_VEF_Face::viscA(int num_face, int num2,
                                              int num_elem, double diffu) const
{
  const Domaine_VEF& domaine = le_dom_vef.valeur();
  const IntTab& face_voisins = domaine.face_voisins();
  const DoubleTab& face_normales = domaine.face_normales();
  const DoubleVect& inverse_volumes = domaine.inverse_volumes();
  double pscal = face_normales(num_face, 0)*face_normales(num2, 0)
                 + face_normales(num_face, 1)*face_normales(num2, 1);
  if (Objet_U::dimension == 3)
    pscal += face_normales(num_face, 2)*face_normales(num2, 2);

  // *equation().milieu().porosite_elem(num_elem));
  double signe = 1;
  if ((face_voisins(num_face, 0) == face_voisins(num2, 0)) ||
      (face_voisins(num_face, 1) == face_voisins(num2, 1)))
    signe = -1;
  return signe*(pscal*diffu)*inverse_volumes(num_elem);
}

/*! @brief on dimensionne notre matrice au moyen de la methode dimensionner de la classe Op_VEF_Face.
 *
 */
inline void Op_Diff_K_Omega_VEF_Face::dimensionner(Matrice_Morse& matrice) const
{
  Op_VEF_Face::dimensionner(le_dom_vef.valeur(), la_zcl_vef.valeur(), matrice);
}

/*! @brief on assemble la matrice des inconnues implicite.
 *
 */
inline void Op_Diff_K_Omega_VEF_Face::contribuer_a_avec(const DoubleTab& inco,
                                                        Matrice_Morse& matrice) const
{
  ajouter_contribution(inco, matrice);
}

/*! @brief on ajoute la contribution du second membre.
 *
 */
inline void Op_Diff_K_Omega_VEF_Face::contribuer_au_second_membre(DoubleTab& resu) const
{
  contribue_au_second_membre(resu);
}
#endif
