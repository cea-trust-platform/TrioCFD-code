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
// File:        Source_Transport_K_VDF_Elem.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <Source_Transport_K_VDF_Elem.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Source_Transport_K_VDF_Elem,"Source_Transport_K_VDF_P0_VDF",Source_Transport_VDF_Elem_base);

Sortie& Source_Transport_K_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_VDF_Elem::readOn(Entree& is)
{
  Source_Transport_VDF_Elem_base::verifier_pb_keps(mon_equation->probleme(),que_suis_je());
  return Source_Transport_VDF_Elem_base::readOn_nothing(is);
}

void Source_Transport_K_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_VDF_Elem_base::associer_pb(pb);
  mon_eq_transport_K = ref_cast(Transport_K_ou_Eps,equation());
  mon_eq_transport_Eps = ref_cast(Transport_K_ou_Eps, mon_eq_transport_K.valeur().modele_turbulence().eqn_transp_Eps());
}

const DoubleTab& Source_Transport_K_VDF_Elem::get_visc_turb() const
{
  return mon_eq_transport_K->modele_turbulence().viscosite_turbulente().valeurs();
}

void Source_Transport_K_VDF_Elem::calculer_terme_production(const Champ_Face& vitesse, const DoubleTab& visco_turb, const DoubleTab& vit, DoubleVect& P) const
{
  const DoubleTab& K   = mon_eq_transport_K->inconnue().valeurs();

  if (axi) calculer_terme_production_K_BiK_Axi(la_zone_VDF.valeur(),vitesse,P,K,visco_turb);
  else calculer_terme_production_K_BiK(la_zone_VDF.valeur(),la_zone_Cl_VDF.valeur(),P,K,vit,vitesse,visco_turb);
}

const Modele_Fonc_Bas_Reynolds& Source_Transport_K_VDF_Elem::get_modele_fonc_bas_reyn() const
{
  return ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,mon_eq_transport_K->modele_turbulence()).associe_modele_fonction();
}

void Source_Transport_K_VDF_Elem::calcul_D_E(const DoubleTab& vit, const DoubleTab& visco_turb, const Champ_Don& ch_visco_cin, DoubleTab& D, DoubleTab& E) const
{
  const DoubleTab& K = mon_eq_transport_K->inconnue().valeurs(), &Eps = mon_eq_transport_Eps->inconnue().valeurs();
  get_modele_fonc_bas_reyn().Calcul_D_BiK(D,mon_eq_transport_K->zone_dis(),mon_eq_transport_K->zone_Cl_dis(),vit,K, Eps,ch_visco_cin);
  get_modele_fonc_bas_reyn().Calcul_E_BiK(E,mon_eq_transport_K->zone_dis(),mon_eq_transport_K->zone_Cl_dis(),vit,K, Eps,ch_visco_cin,visco_turb);
}

void Source_Transport_K_VDF_Elem::calcul_F1_F2(const Champ_base& ch_visco_cin_ou_dyn, DoubleTab& P_tab, DoubleTab& D, DoubleTab& F1, DoubleTab& F2) const
{
  const DoubleTab& K = mon_eq_transport_K->inconnue().valeurs(), &Eps = mon_eq_transport_Eps->inconnue().valeurs();
  get_modele_fonc_bas_reyn().Calcul_F1_BiK(F1,mon_eq_transport_K->zone_dis(),mon_eq_transport_K->zone_Cl_dis(), P_tab, K, Eps,ch_visco_cin_ou_dyn);
  get_modele_fonc_bas_reyn().Calcul_F2_BiK(F2,D,mon_eq_transport_K->zone_dis(),K, Eps, ch_visco_cin_ou_dyn  );
}

void Source_Transport_K_VDF_Elem::fill_resu_bas_rey(const DoubleVect& P, const DoubleTab& D, const DoubleTab& E, const DoubleTab& F1, const DoubleTab& F2, DoubleTab& resu) const
{
  const DoubleVect& volumes = la_zone_VDF->volumes(), &porosite_vol = la_zone_VDF->porosite_elem();
  const DoubleTab& Eps = mon_eq_transport_Eps->inconnue().valeurs();
  for (int elem = 0; elem < la_zone_VDF->nb_elem(); elem++)
    resu(elem) += (P(elem)-Eps(elem)-D(elem))*volumes(elem)*porosite_vol(elem);
}

void Source_Transport_K_VDF_Elem::fill_resu(const DoubleVect& P, DoubleTab& resu) const
{
  const DoubleVect& volumes = la_zone_VDF->volumes(), &porosite_vol = la_zone_VDF->porosite_elem();
  const DoubleTab& Eps = mon_eq_transport_Eps->inconnue().valeurs();
  for (int elem = 0; elem < la_zone_VDF->nb_elem(); elem++)
    resu(elem) += (P(elem)-Eps(elem))*volumes(elem)*porosite_vol(elem);
}

void Source_Transport_K_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Source_Transport_VDF_Elem_base::ajouter_keps(secmem);

  const std::string& nom_inco = equation().inconnue().le_nom().getString();
  Matrice_Morse* mat = matrices.count(nom_inco) ? matrices.at(nom_inco) : NULL;
  if(!mat) return;

  const DoubleTab& K   = mon_eq_transport_K->inconnue().valeurs(), &Eps = mon_eq_transport_Eps->inconnue().valeurs();
  const DoubleVect& porosite = la_zone_VDF->porosite_elem(), &volumes = la_zone_VDF->volumes();
  const int size = K.dimension(0);
  // on implicite le -eps et le -eps^2/k
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc=ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,mon_eq_transport_K->modele_turbulence()).associe_modele_fonction();
  const int is_modele_fonc=(mon_modele_fonc.non_nul());
  DoubleTab F2;
  if (is_modele_fonc)
    {
      DoubleTrav D(0);
      F2.resize(K.dimension_tot(0));
      const Zone_dis& zone_dis_k =mon_eq_transport_K->zone_dis();
      const Champ_base& ch_visco_cin_ou_dyn =ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();
      mon_modele_fonc.Calcul_F2_BiK(F2,D,zone_dis_k,K,Eps, ch_visco_cin_ou_dyn  );
    }

  for (int c=0; c<size; c++)
    {
      // -eps*vol  donne +vol dans la bonne case
      if (K(c)>DMINFLOAT)
        {
          double coef_k=porosite(c)*volumes(c)*Eps(c)/K(c);
          (*mat)(c,c)+=coef_k;
        }
    }
}
