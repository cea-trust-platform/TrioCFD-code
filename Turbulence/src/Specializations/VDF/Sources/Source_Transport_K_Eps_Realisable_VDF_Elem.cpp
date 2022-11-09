/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Source_Transport_K_Eps_Realisable_VDF_Elem.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_Realisable_VDF_Elem.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Milieu_base.h>
#include <TRUSTTrav.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_Realisable_VDF_Elem,"Source_Transport_K_Eps_Realisable_VDF_P0_VDF",Source_Transport_Realisable_VDF_Elem_base);

Sortie& Source_Transport_K_Eps_Realisable_VDF_Elem::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_Eps_Realisable_VDF_Elem::readOn(Entree& is ) { return Source_Transport_Realisable_VDF_Elem_base::readOn(is); }

void Source_Transport_K_Eps_Realisable_VDF_Elem::associer_pb(const Probleme_base& pb )
{
  Source_Transport_Realisable_VDF_Elem_base::associer_pb(pb);
  eqn_keps_Rea = ref_cast(Transport_K_Eps_Realisable,equation());
}

const DoubleTab& Source_Transport_K_Eps_Realisable_VDF_Elem::get_visc_turb() const
{
  return eqn_keps_Rea->modele_turbulence().viscosite_turbulente().valeurs();
}

const Modele_Fonc_Realisable_base& Source_Transport_K_Eps_Realisable_VDF_Elem::get_modele_fonc() const
{
  const Modele_turbulence_hyd_K_Eps_Realisable& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence());
  return mod_turb.associe_modele_fonction();
}

void Source_Transport_K_Eps_Realisable_VDF_Elem::calculer_terme_production_real(const Champ_Face& vitesse, const DoubleTab& visco_turb, const DoubleTab& vit, DoubleTrav& P) const
{
  const DoubleTab& K_eps_Rea = eqn_keps_Rea->inconnue().valeurs();

  if (axi) calculer_terme_production_K_Axi(la_zone_VDF.valeur(),vitesse,P,K_eps_Rea,visco_turb);
  else calculer_terme_production_K(la_zone_VDF.valeur(),la_zone_Cl_VDF.valeur(),P,K_eps_Rea,vit,vitesse,visco_turb);
}

void Source_Transport_K_Eps_Realisable_VDF_Elem::fill_resu_real(const int is_visco_const, const DoubleTab& tab_visco, const DoubleTrav& P, const DoubleTrav& CC1, const DoubleTrav& S, double& visco, DoubleTab& resu) const
{
  const DoubleTab& K_eps_Rea = eqn_keps_Rea->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Realisable& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence());
  const DoubleVect& volumes = la_zone_VDF->volumes(), &porosite_vol = la_zone_Cl_VDF->equation().milieu().porosite_elem();
  const double LeK_MIN = mod_turb.get_LeK_MIN(), LeEPS_MIN = mod_turb.get_LeEPS_MIN();

  for (int elem = 0; elem < la_zone_VDF->nb_elem(); elem++)
    {
      if (!is_visco_const) visco =  tab_visco(elem);
      assert(visco > 0.);
      resu(elem,0) += ( P(elem)-K_eps_Rea(elem,1) )*volumes(elem)*porosite_vol(elem);
      if ( ( K_eps_Rea(elem,0) >= LeK_MIN ) and ( K_eps_Rea(elem,1) >= LeEPS_MIN ) )
        resu(elem,1) += K_eps_Rea(elem,1)*( CC1(elem)*S(elem) - ( C2*K_eps_Rea(elem,1)/( K_eps_Rea(elem,0) + sqrt( visco*K_eps_Rea(elem,1) ) ) ) )*volumes(elem)*porosite_vol(elem);
    }
}

void Source_Transport_K_Eps_Realisable_VDF_Elem::mettre_a_jour(double temps)
{
  const DoubleTab& K_eps_Rea = eqn_keps_Rea->inconnue().valeurs(), &vit = eq_hydraulique->inconnue().valeurs();
  const double epsilon_minimum = eqn_keps_Rea->modele_turbulence().get_LeEPS_MIN();
  Modele_Fonc_Realisable_base& mon_modele_fonc = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_keps_Rea->modele_turbulence()).associe_modele_fonction();
  mon_modele_fonc.Contributions_Sources(eqn_keps_Rea ->zone_dis(),eqn_keps_Rea->zone_Cl_dis(),vit,K_eps_Rea,epsilon_minimum);
  Source_Transport_Realisable_VDF_Elem_base::mettre_a_jour(temps);
}

void Source_Transport_K_Eps_Realisable_VDF_Elem::fill_coeff_matrice(const int is_visco_const, const DoubleTab& tab_visco, const DoubleVect& volumes, const DoubleVect& porosite, double& visco, Matrice_Morse& matrice) const
{
  const DoubleTab& val = equation().inconnue().valeurs();
  for (int c = 0; c < val.dimension(0); c++)
    {
      if (!is_visco_const) visco = tab_visco(c);
      assert(visco > 0.);
      if (val(c,0) > DMINFLOAT)
        {
          double coef_k = porosite(c)*volumes(c)*val(c,1)/( val(c,0) + sqrt( visco*val(c,1) ) );
          matrice(c*2,c*2) += coef_k;
          double coef_eps = C2*coef_k;
          matrice(c*2+1,c*2+1) += coef_eps;
        }
    }
}

void Source_Transport_K_Eps_Realisable_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Source_Transport_Realisable_VDF_Elem_base::ajouter_blocs(matrices, secmem, semi_impl);
  Source_Transport_Realisable_VDF_Elem_base::ajouter_keps_real(secmem);

}

