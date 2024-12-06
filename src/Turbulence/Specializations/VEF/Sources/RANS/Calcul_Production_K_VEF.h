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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Calcul_Production_K_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Calcul_Production_K_VEF_included
#define Calcul_Production_K_VEF_included

#include <TRUSTTabs_forward.h>

class Champ_Don_base;
class Domaine_Cl_VEF;
class Domaine_VEF;

//  CLASS Calcul_Production_K_VEF
// Classe qui porte les fonctions de calcul des termes de production de l'energie cinetique turbulente et de destruction
// de cette energie. Cette classe ne derive d'aucune autre car nous voulons l'utiliser pour faire de l'heritage multiple.
class Calcul_Production_K_VEF
{
protected:
  Calcul_Production_K_VEF() { }

  // Standard TKE production
  DoubleTab& calculer_terme_production_K(const Domaine_VEF&, const Domaine_Cl_VEF&, DoubleTab&, const DoubleTab&, const DoubleTab&, const DoubleTab&, const int& interpol_visco, const double& limiteur) const;
  void loop_for_internal_or_periodic_faces(DoubleTab& prodK, const DoubleTab& gradient_elem, const DoubleTab& visco_turb, const DoubleVect& volumes, const IntTab& face_voisins, const int nfaceinit, const int nfaceend, const int interpol_visco, const double limiteur) const;
  void loop_for_non_periodic_boundaries(DoubleTab& prodK, const DoubleTab& gradient_elem, const DoubleTab& visco_turb, const DoubleVect& volumes, const IntTab& face_voisins, const int nfaceinit, const int nfaceend, const int interpol_visco, const double limiteur) const;


  DoubleTab& calculer_terme_production_K_BiK(const Domaine_VEF&, const Domaine_Cl_VEF&, DoubleTab&, const DoubleTab&, const DoubleTab&, const DoubleTab&, const DoubleTab&, const int& interpol_visco, const double& limiteur) const;

  // EASM
  DoubleTab& calculer_terme_production_K_EASM(const Domaine_VEF&, const Domaine_Cl_VEF&, DoubleTab&, const DoubleTab&, const DoubleTab&, const DoubleTab&, const DoubleTab&, const int& interpol_visco, const double& limiteur) const;
  void compute_production_term_EASM(const int face, const double visco_face, const DoubleTab& Re_face, const DoubleTab& gradient_face, DoubleTab& P) const;
  DoubleTab& calcul_tenseur_face(DoubleTab&, const DoubleTab&, const Domaine_VEF&, const Domaine_Cl_VEF&) const;

  // Commons
  double get_turbulent_viscosity(const DoubleTab& visco_turb, const DoubleVect& volumes, const int type_interpo, const int poly1, const int poly2, const double limiteur) const;

  // TKE destruction
  DoubleTab& calculer_terme_destruction_K_gen(const Domaine_VEF&, const Domaine_Cl_VEF&, DoubleTab&, const DoubleTab&, const DoubleTab&, const Champ_Don_base&, const DoubleVect&, int) const;
  void compute_utheta_nbConsti_le_1_nbCompo_eq_0(const Domaine_VEF& domaine_VEF, const Domaine_Cl_VEF& zcl_VEF, const IntTab& face_voisins, const DoubleVect& volumes, const DoubleTab& tab_beta, const DoubleTab& alpha_turb, const DoubleTrav& gradient_elem, DoubleTrav& u_theta) const;
  void compute_utheta_nbConsti_le_1_nbCompo_eq_1(const Domaine_VEF& domaine_VEF, const Domaine_Cl_VEF& zcl_VEF, const IntTab& face_voisins, const DoubleVect& volumes, const DoubleTab& tab_beta, const DoubleTab& alpha_turb, const DoubleTrav& gradient_elem, DoubleTrav& u_theta) const;
  void compute_utheta_nbConsti_le_1_nbCompo_gt_1(const Domaine_VEF& domaine_VEF, const Domaine_Cl_VEF& zcl_VEF, const IntTab& face_voisins, const DoubleVect& volumes, const DoubleTab& tab_beta, const DoubleTab& alpha_turb, const DoubleTrav& gradient_elem, DoubleTrav& u_theta) const;
  void compute_utheta_nbConsti_gt_1_nbCompo_eq_0(const Domaine_VEF& domaine_VEF, const Domaine_Cl_VEF& zcl_VEF, const IntTab& face_voisins, const DoubleVect& volumes, const DoubleTab& tab_beta, const DoubleTab& alpha_turb, const DoubleTrav& gradient_elem, const int nb_consti, DoubleTrav& u_theta) const;
  void compute_utheta_nbConsti_gt_1_nbCompo_eq_1(const Domaine_VEF& domaine_VEF, const Domaine_Cl_VEF& zcl_VEF, const IntTab& face_voisins, const DoubleVect& volumes, const DoubleTab& tab_beta, const DoubleTab& alpha_turb, const DoubleTrav& gradient_elem, const int nb_consti, DoubleTrav& u_theta) const;
  void compute_utheta_nbConsti_gt_1_nbCompo_gt_1(const Domaine_VEF& domaine_VEF, const Domaine_Cl_VEF& zcl_VEF, const IntTab& face_voisins, const DoubleVect& volumes, const DoubleTab& tab_beta, const DoubleTab& alpha_turb, const DoubleTrav& gradient_elem, const int nb_consti, DoubleTrav& u_theta) const;

  void mettre_a_jour(double temps) { }
};

#endif /* Calcul_Production_K_VEF_included */
