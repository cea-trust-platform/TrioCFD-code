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
// File:        Modele_turbulence_hyd_LES_SMAGO_DYN_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_LES_SMAGO_DYN_VDF_included
#define Modele_turbulence_hyd_LES_SMAGO_DYN_VDF_included

#include <Modele_turbulence_hyd_LES_Smago_VDF.h>

class Modele_turbulence_hyd_LES_SMAGO_DYN_VDF: public Modele_turbulence_hyd_LES_Smago_VDF
{
  Declare_instanciable(Modele_turbulence_hyd_LES_SMAGO_DYN_VDF);
public:
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis) override;
  void mettre_a_jour(double) override;
  int preparer_calcul() override;

  static void calculer_length_scale(DoubleVect&, const Domaine_VDF&);
  static void calculer_cell_cent_vel(DoubleTab&, const Domaine_VDF&, Champ_Inc_base&);
  static void calculer_filter_field(const DoubleTab&, DoubleTab&, const Domaine_VDF&);
  static void calculer_Sij(DoubleTab&, const Domaine_VDF&, const Domaine_Cl_VDF&, Champ_Inc_base&);
  static void calculer_Sij_vel_filt(const DoubleTab&, DoubleTab&, const Domaine_VDF&);
  static void calculer_S_norme(const DoubleTab&, DoubleVect&, int);
  static void interpole(const IntVect&, const DoubleVect&, const DoubleVect&, double&);
protected:
  OWN_PTR(Champ_Fonc_base)  coeff_field_;
  Motcle methode_stabilise_;
  int N_c_ = -123;
  IntVect compt_c_;
  IntVect corresp_c_;
  IntTab elem_elem_;
  DoubleTab cell_cent_vel_;

  void calculer_filter_tensor(DoubleTab&);
  void calculer_Lij(const DoubleTab&, const DoubleTab&, DoubleTab&);
  void calculer_Mij(const DoubleTab&, const DoubleTab&, const DoubleVect&, DoubleTab&);
  void calculer_model_coefficient(const DoubleTab&, const DoubleTab&);
  Champ_Fonc_base& calculer_viscosite_turbulente(const DoubleVect&, const DoubleVect&);
  Champ_Fonc_base& calculer_viscosite_turbulente() override;
  void calculer_energie_cinetique_turb() override;
  Champ_Fonc_base& calculer_energie_cinetique_turb(const DoubleVect&, const DoubleVect&);
  void controler_grandeurs_turbulentes();

  void stabilise_moyenne(const DoubleTab&, const DoubleTab&);
  void stabilise_moyenne_6_points(const DoubleTab&, const DoubleTab&);
  void stabilise_moyenne_plans_paralleles(const DoubleTab&, const DoubleTab&);
  void stabilise_moyenne_euler_lagrange(const DoubleTab&, const DoubleTab&);
  void calcul_voisins(const int, IntVect&, DoubleVect&);
  void calc_elem_elem();
  void calcul_tableaux_correspondance(int&, IntVect&, IntVect&);
};

#endif /* Modele_turbulence_hyd_LES_SMAGO_DYN_VDF_included */
