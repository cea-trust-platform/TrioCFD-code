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

#include <Transport_Interfaces_FT_Disc.h>
#include <Dirichlet_paroi_fixe.h>
#include <Fluide_Diphasique.h>
#include <Champ_Face_VDF.h>
#include <Champ_Face_VDF.h>
#include <Champ_Uniforme.h>
#include <Probleme_base.h>
#include <Zone_Cl_VDF.h>

void moy_2D_vit(const DoubleVect& vit, int elem, int iori, const Zone_VDF& zone, double& u)
{
  int num1, num2;
  if (iori == 0)
    {
      num1 = zone.elem_faces(elem, 1);
      num2 = zone.elem_faces(elem, 3);
    }
  else if (iori == 1)
    {
      num1 = zone.elem_faces(elem, 0);
      num2 = zone.elem_faces(elem, 2);
    }
  else
    {
      Cerr << "valeur de iori " << iori << " impossible en 2D" << finl;
      Process::exit();
      num1 = num2 = -1;
    }
  u = 0.5 * (vit(num1) + vit(num2));
}

double norm_2D_vit(const DoubleVect& vit, int elem, int iori, const Zone_VDF& zone, double& u)
{
  double v;
  moy_2D_vit(vit, elem, iori, zone, v);
  v = std::fabs(v);
  if (v == 0)
    u = 0;
  else
    u = 1;
  return v;
}

void moy_3D_vit(const DoubleVect& vit, int elem, int iori, const Zone_VDF& zone, double& val1, double& val2)
{
  int num1, num2, num3, num4;
  if (iori == 0)
    {
      num1 = zone.elem_faces(elem, 1);
      num2 = zone.elem_faces(elem, 4);
      num3 = zone.elem_faces(elem, 2);
      num4 = zone.elem_faces(elem, 5);
    }
  else if (iori == 1)
    {
      num1 = zone.elem_faces(elem, 0);
      num2 = zone.elem_faces(elem, 3);
      num3 = zone.elem_faces(elem, 2);
      num4 = zone.elem_faces(elem, 5);
    }
  else if (iori == 2)
    {
      num1 = zone.elem_faces(elem, 0);
      num2 = zone.elem_faces(elem, 3);
      num3 = zone.elem_faces(elem, 1);
      num4 = zone.elem_faces(elem, 4);
    }
  else
    {
      Cerr << "valeur de iori " << iori << " impossible en 3D" << finl;
      Process::exit();
      num1 = num2 = num3 = num4 = -1;
    }
  val1 = 0.5 * (vit(num1) + vit(num2));
  val2 = 0.5 * (vit(num3) + vit(num4));
}

double norm_3D_vit(const DoubleVect& vit, int elem, int iori, const Zone_VDF& zone, double& val1, double& val2)
{
  moy_3D_vit(vit, elem, iori, zone, val1, val2);
  double v1 = std::fabs(val1);
  double v2 = std::fabs(val2);
  double norm_vit = sqrt(v1 * v1 + v2 * v2);
  val1 = v1 / (norm_vit + DMINFLOAT);
  val2 = v2 / (norm_vit + DMINFLOAT);
  return norm_vit;
}

void Champ_Face_VDF::calcul_y_plus_diphasique(DoubleTab& y_plus, const Zone_Cl_VDF& zone_Cl_VDF)
{
  // On initialise le champ y_plus avec une valeur negative,
  // comme ca lorsqu'on veut visualiser le champ pres de la paroi,
  // on n'a qu'a supprimer les valeurs negatives et n'apparaissent
  // que les valeurs aux parois.

  int ndeb, nfin, elem, ori, l_unif;
  double norm_tau, u_etoile, norm_v = 0, dist, val0, val1, val2, d_visco = 0, visco_ph0 = 1.;
  y_plus = -1.;

  const Champ_Face_VDF& vit = *this;
  const Zone_VDF& zone_VDF = zone_vdf();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const IntVect& orientation = zone_VDF.orientation();
  const Equation_base& eqn_hydr = equation();

  // Physical properties of both phases
  const Fluide_Diphasique& le_fluide = ref_cast(Fluide_Diphasique, eqn_hydr.milieu());
  const Fluide_base& phase_1 = le_fluide.fluide_phase(1);
  const Fluide_base& phase_0 = le_fluide.fluide_phase(0);
  const Champ_Don& ch_visco_cin_ph1 = phase_1.viscosite_cinematique();
  const Champ_Don& ch_visco_cin_ph0 = phase_0.viscosite_cinematique();
  const DoubleTab& tab_visco_ph1 = phase_1.viscosite_cinematique().valeur().valeurs();
  const DoubleTab& tab_visco_ph0 = phase_0.viscosite_cinematique().valeur().valeurs();
  const double delta_nu = tab_visco_ph1(0, 0) - tab_visco_ph0(0, 0);

  // One way to get the Transport equation to pass the indicator DoubleTab
  const Zone_Cl_dis_base& zone_Cl_dis_base = eqn_hydr.zone_Cl_dis().valeur();
  const Equation_base& eqn_trans = zone_Cl_dis_base.equation().probleme().equation("Transport_Interfaces_FT_Disc");
  const Transport_Interfaces_FT_Disc& eqn_interf = ref_cast(Transport_Interfaces_FT_Disc, eqn_trans);
  const DoubleTab& indic = eqn_interf.inconnue().valeurs();

  if (sub_type(Champ_Uniforme,ch_visco_cin_ph1.valeur()) && sub_type(Champ_Uniforme, ch_visco_cin_ph0.valeur()))
    {
      visco_ph0 = std::max(tab_visco_ph0(0, 0), DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;

  if ((!l_unif) && ((tab_visco_ph1.local_min_vect() < DMINFLOAT) || (tab_visco_ph0.local_min_vect() < DMINFLOAT)))
    {
      Cerr << "Negative viscosity !!!" << finl;
      Process::exit(-1);
    }

  DoubleTab yplus_faces(1, 1); // will contain yplus values if available
  int yplus_already_computed = 0; // flag

  const RefObjU& modele_turbulence = eqn_hydr.get_modele(TURBULENCE);
  if (modele_turbulence.non_nul() && sub_type(Mod_turb_hyd_base, modele_turbulence.valeur()))
    {
      const Mod_turb_hyd_base& mod_turb = ref_cast(Mod_turb_hyd_base, modele_turbulence.valeur());
      const Turbulence_paroi_base& loipar = mod_turb.loi_paroi();
      yplus_faces.resize(zone_VDF.nb_faces_tot());
      yplus_faces.ref(loipar.tab_d_plus());
      yplus_already_computed = 1;
    }

  for (int n_bord = 0; n_bord < zone_VDF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VDF.les_conditions_limites(n_bord);

      if (sub_type(Dirichlet_paroi_fixe, la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (int num_face = ndeb; num_face < nfin; num_face++)
            {

              if (face_voisins(num_face, 0) != -1)
                elem = face_voisins(num_face, 0);
              else
                elem = face_voisins(num_face, 1);

              if (yplus_already_computed)
                {
                  // y+ is only defined on faces so we take the face value to put in the element
                  y_plus(elem) = yplus_faces(num_face);
                }
              else
                {
                  if (dimension == 2)
                    {
                      ori = orientation(num_face);
                      norm_v = norm_2D_vit(vit.valeurs(), elem, ori, zone_VDF, val0);
                    }
                  else if (dimension == 3)
                    {
                      ori = orientation(num_face);
                      norm_v = norm_3D_vit(vit.valeurs(), elem, ori, zone_VDF, val1, val2);
                    } // dim 3

                  if (axi)
                    dist = zone_VDF.dist_norm_bord_axi(num_face);
                  else
                    dist = zone_VDF.dist_norm_bord(num_face);

                  if (l_unif)
                    d_visco = visco_ph0 + indic(elem) * delta_nu;
                  else
                    d_visco = (tab_visco_ph0.nb_dim() == 1 ? (tab_visco_ph0(elem) + indic(elem) * delta_nu) : (tab_visco_ph0(elem, 0) + indic(elem) * delta_nu));

                  // PQ : 01/10/03 : corrections par rapport a la version premiere
                  norm_tau = d_visco * norm_v / dist;

                  u_etoile = sqrt(norm_tau);
                  y_plus(elem) = dist * u_etoile / d_visco;

                } // else yplus already computed
            } // loop on faces
        } // Fin paroi fixe
    } // Fin boucle sur les bords
}
