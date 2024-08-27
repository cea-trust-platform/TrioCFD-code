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
// File:        Source_Transport_K_Eps_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <TRUSTTabs_forward.h>
#include <Source_Transport_K_Eps_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Transport_K_Eps.h>
#include <Milieu_base.h>
#include <Domaine_VEF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_VEF_Face,
                                          "Source_Transport_K_Eps_VEF_P1NC",
                                          Source_Transport_VEF_Face_base);

Sortie& Source_Transport_K_Eps_VEF_Face::printOn(Sortie& s) const { return s << que_suis_je() ; }

Entree& Source_Transport_K_Eps_VEF_Face::readOn(Entree& is)
{
  Source_Transport_VEF_Face_base::verifier_pb_keps(mon_equation->probleme(),que_suis_je());
  return Source_Transport_VEF_Face_base::readOn(is);
}

void Source_Transport_K_Eps_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_VEF_Face_base::associer_pb(pb);
  mon_eq_transport_K_Eps = ref_cast(Transport_K_Eps, equation());
}
const DoubleTab& Source_Transport_K_Eps_VEF_Face::get_visc_turb() const
{
  return mon_eq_transport_K_Eps->modele_turbulence().viscosite_turbulente()->valeurs();
}

const DoubleTab& Source_Transport_K_Eps_VEF_Face::get_cisaillement_paroi() const
{
  const Modele_turbulence_hyd_K_Eps& mod  = ref_cast(Modele_turbulence_hyd_K_Eps,mon_eq_transport_K_Eps->modele_turbulence());
  return mod.loi_paroi()->Cisaillement_paroi();
}

const DoubleTab& Source_Transport_K_Eps_VEF_Face::get_K_pour_production() const
{
  return mon_eq_transport_K_Eps->inconnue()->valeurs();
}

const OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)& Source_Transport_K_Eps_VEF_Face::get_modele_fonc_bas_reyn() const
{
  return ref_cast(Modele_turbulence_hyd_K_Eps,mon_eq_transport_K_Eps->modele_turbulence()).associe_modele_fonction();
}

void Source_Transport_K_Eps_VEF_Face::calcul_tabs_bas_reyn(const DoubleTrav& P, const DoubleTab& vit, const DoubleTab& visco_turb, const Champ_Don& ch_visco_cin,
                                                           const Champ_base& ch_visco_cin_ou_dyn, DoubleTab& D, DoubleTab& E, DoubleTab& F1, DoubleTab& F2) const
{
  const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue()->valeurs();
  get_modele_fonc_bas_reyn()->Calcul_D(D,mon_eq_transport_K_Eps->domaine_dis(),mon_eq_transport_K_Eps->domaine_Cl_dis(),vit,K_eps,ch_visco_cin);
  get_modele_fonc_bas_reyn()->Calcul_E(E,mon_eq_transport_K_Eps->domaine_dis(),mon_eq_transport_K_Eps->domaine_Cl_dis(),vit,K_eps,ch_visco_cin,visco_turb);
  D.echange_espace_virtuel();
  E.echange_espace_virtuel();
  get_modele_fonc_bas_reyn()->Calcul_F1(F1,mon_eq_transport_K_Eps->domaine_dis(),mon_eq_transport_K_Eps->domaine_Cl_dis(), P, K_eps,ch_visco_cin_ou_dyn);
  get_modele_fonc_bas_reyn()->Calcul_F2(F2,D,mon_eq_transport_K_Eps->domaine_dis(),K_eps, ch_visco_cin_ou_dyn  );
}

const Nom Source_Transport_K_Eps_VEF_Face::get_type_paroi() const
{
  const Modele_turbulence_hyd_K_Eps& mod  = ref_cast(Modele_turbulence_hyd_K_Eps,mon_eq_transport_K_Eps->modele_turbulence());
  return mod.loi_paroi()->que_suis_je();
}

void Source_Transport_K_Eps_VEF_Face::calcul_tenseur_reyn(const DoubleTab& visco_turb, const DoubleTab& gradient_elem, DoubleTab& Re) const
{
  get_modele_fonc_bas_reyn()->calcul_tenseur_Re(visco_turb, gradient_elem, Re);
}

void Source_Transport_K_Eps_VEF_Face::fill_resu_bas_rey(const DoubleVect& volumes_entrelaces, const DoubleTrav& P, const DoubleTab& D, const DoubleTab& E, const DoubleTab& F1, const DoubleTab& F2, DoubleTab& resu) const
{
  const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue()->valeurs();
  const double LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_K_MIN();
  for (int fac = 0; fac < le_dom_VEF->nb_faces(); fac++)
    {
      resu(fac,0) += (P(fac)-K_eps(fac,1)-D(fac))*volumes_entrelaces(fac);
      if (K_eps(fac,0) >= LeK_MIN)
        resu(fac,1) += ((C1*P(fac)*F1(fac)-C2*K_eps(fac,1)*F2(fac))*K_eps(fac,1)/(K_eps(fac,0))+E(fac)) *volumes_entrelaces(fac);
    }
}

void Source_Transport_K_Eps_VEF_Face::fill_resu(const DoubleVect& volumes_entrelaces, const DoubleTrav& P, DoubleTab& resu) const
{
  const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue()->valeurs();
  const double LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_K_MIN();
  for (int fac = 0; fac < le_dom_VEF->nb_faces(); fac++)
    {
      resu(fac,0) += (P(fac)-K_eps(fac,1))*volumes_entrelaces(fac);
      if (K_eps(fac,0) >= LeK_MIN)
        resu(fac,1) += (C1*P(fac)-C2*K_eps(fac,1))*volumes_entrelaces(fac)*K_eps(fac,1)/K_eps(fac,0);
    }
}

DoubleTab& Source_Transport_K_Eps_VEF_Face::ajouter(DoubleTab& resu) const
{
  return Source_Transport_VEF_Face_base::ajouter_keps(resu);
}

void Source_Transport_K_Eps_VEF_Face::contribuer_a_avec(const DoubleTab& a, Matrice_Morse& matrice) const
{
  const DoubleTab& K_eps = equation().inconnue()->valeurs();
  const double LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_K_MIN();
  const DoubleVect& porosite_face = mon_eq_transport_K_Eps->milieu().porosite_face();
  const DoubleVect& volumes_entrelaces = le_dom_VEF->volumes_entrelaces();
  // on implicite le -eps et le -eps^2/k
  for (int face = 0; face < K_eps.dimension(0); face++)
    if (K_eps(face, 0) >= LeK_MIN) // -eps*vol  donne +vol dans la bonne case
      {
        double coef_k = porosite_face(face) * volumes_entrelaces(face) * K_eps(face, 1) / K_eps(face, 0);
        matrice(face * 2, face * 2) += coef_k;
        double coef_eps = C2 * K_eps(face, 1) / K_eps(face, 0) * volumes_entrelaces(face) * porosite_face(face);
        matrice(face * 2 + 1, face * 2 + 1) += coef_eps;
      }
}
