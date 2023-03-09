/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Source_Robin_Scalaire.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Robin_Scalaire.h>
#include <Domaine_dis.h>
#include <Domaine_Cl_dis.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Equation_base.h>
#include <Fluide_base.h>
#include <distances_VEF.h>
#include <Navier_Stokes_Turbulent.h>
#include <Source_Robin.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Paroi_decalee_Robin.h>
#include <Paroi_std_scal_hyd_VEF.h>

Implemente_instanciable(Source_Robin_Scalaire,"Source_Robin_Scalaire_VEF_P1NC",Source_base);


// printOn
Sortie& Source_Robin_Scalaire::printOn(Sortie& s) const
{
  return s << que_suis_je();
}


// readOn
Entree& Source_Robin_Scalaire::readOn(Entree& s)
{
  int nb;
  s >> nb;
  noms_parois.dimensionner(nb);
  T_parois.resize(nb);
  for (int i=0; i<nb; i++) s>> noms_parois[i] >> T_parois(i);
//  s >> dt_post;
  return s;
}


// associer_pb
void Source_Robin_Scalaire::associer_pb(const Probleme_base& pb)
{
}


// completer
void Source_Robin_Scalaire::completer()
{
  Source_base::completer();
  const Equation_base& eqn_NS = equation().probleme().equation(0);
  const Sources& srcs = eqn_NS.sources();
  for (const auto& itr : srcs)
    {
      if (sub_type(Source_Robin,itr.valeur()))
        {
          const Source_Robin& ts_Robin = ref_cast(Source_Robin,itr.valeur());
          tab_u_star = ts_Robin.tab_u_star();
          tab_d_plus = ts_Robin.tab_d_plus();
        }
    }
}


// ajouter
DoubleTab& Source_Robin_Scalaire::ajouter(DoubleTab& resu) const
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const Domaine_Cl_VEF& domaine_Cl_VEF = le_dom_Cl_VEF.valeur();
  const DoubleTab& temperature = equation().inconnue().valeurs();
  const Fluide_base& fluide = ref_cast(Fluide_base,equation().milieu());
  const Champ_Don& lambda = fluide.conductivite();
  int lambda_uniforme = (sub_type(Champ_Uniforme,lambda.valeur()) ? 1 : 0);
  const Convection_Diffusion_Temperature& eq_th = ref_cast(Convection_Diffusion_Temperature,equation());
  const Modele_turbulence_scal_base& le_modele_scalaire = ref_cast(Modele_turbulence_scal_base,eq_th.get_modele(TURBULENCE).valeur());
  const DoubleTab& lambda_t = le_modele_scalaire.conductivite_turbulente().valeurs();
  const Paroi_scal_hyd_base_VEF& loi_de_paroi = ref_cast(Paroi_scal_hyd_base_VEF,le_modele_scalaire.loi_paroi().valeur());
  const DoubleVect& surfaces_face = ref_cast_non_const(Domaine_VEF,domaine_VEF).face_surfaces();

  for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      int face;
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      int pos_Paroi = noms_parois.search(la_cl->frontiere_dis().le_nom());
      const DoubleVect& dist_equiv = loi_de_paroi.equivalent_distance(n_bord);

      if (sub_type(Paroi_decalee_Robin,la_cl.valeur()))
        {
          // double acc_loc_tot = 0;
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          double Tw = T_parois(pos_Paroi);

          for (face=ndeb; face<nfin; face++)
            {
              int elem = domaine_VEF.face_voisins(face,0);
              if (elem==-1) elem = domaine_VEF.face_voisins(face,1);
              double d_lambda = (lambda_uniforme ? lambda(0,0) : lambda(elem,0));
              double acc_loc = - (d_lambda + lambda_t(elem)) * (temperature(face) - Tw) / dist_equiv[face-ndeb] * surfaces_face(face);
              // acc_loc_tot += acc_loc;
              resu(face) += acc_loc;
            }
        }
    }
  return resu;
}


// calculer
DoubleTab& Source_Robin_Scalaire::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}


// associer_domaines
void Source_Robin_Scalaire::associer_domaines(const Domaine_dis& z, const Domaine_Cl_dis& zcl)
{
  le_dom_VEF = ref_cast(Domaine_VEF,z.valeur());
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF,zcl.valeur());
}
