/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Source_Transport_K_Realisable_VDF_Elem.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.h>
#include <Source_Transport_K_Realisable_VDF_Elem.h>
#include <Champ_Face.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Source_Transport_K_Realisable_VDF_Elem,"Source_Transport_K_Realisable_VDF_P0_VDF",Source_Transport_Realisable_VDF_Elem_base);

Sortie& Source_Transport_K_Realisable_VDF_Elem::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_Realisable_VDF_Elem::readOn(Entree& is ) { return Source_Transport_Realisable_VDF_Elem_base::readOn_nothing(is,que_suis_je()); }

void Source_Transport_K_Realisable_VDF_Elem::associer_pb(const Probleme_base& pb )
{
  Source_Transport_Realisable_VDF_Elem_base::associer_pb(pb);
  eqn_k_Rea = ref_cast(Transport_K_ou_Eps_Realisable,equation());
  eqn_eps_Rea = ref_cast(Transport_K_ou_Eps_Realisable,eqn_eps_Rea.valeur().modele_turbulence().eqn_transp_Eps());
}

void Source_Transport_K_Realisable_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Source_Transport_Realisable_VDF_Elem_base::ajouter_blocs(matrices, secmem, semi_impl);
  const DoubleTab& K_Rea = eqn_k_Rea->inconnue().valeurs(), &eps_Rea = eqn_eps_Rea->inconnue().valeurs(), &vit = eq_hydraulique->inconnue().valeurs();
  Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());
  const DoubleTab& visco_turb = eqn_k_Rea->modele_turbulence().viscosite_turbulente().valeurs();
  const DoubleVect& volumes = la_zone_VDF->volumes(), &porosite_vol = la_zone_VDF->porosite_elem();
  DoubleTrav P(visco_turb);

  if (axi) calculer_terme_production_K_BiK_Axi(la_zone_VDF.valeur(),vitesse,P,K_Rea,visco_turb);
  else calculer_terme_production_K_BiK(la_zone_VDF.valeur(),la_zone_Cl_VDF.valeur(),P,K_Rea,vit,vitesse,visco_turb);
  P.echange_espace_virtuel();

  for (int elem = 0; elem < la_zone_VDF->nb_elem(); elem++) secmem(elem) += ( P(elem)-eps_Rea(elem) )*volumes(elem)*porosite_vol(elem);

}

void Source_Transport_K_Realisable_VDF_Elem::mettre_a_jour(double temps)
{
  const DoubleTab& K_Rea  = eqn_k_Rea->inconnue().valeurs(), &eps_Rea  = eqn_eps_Rea->inconnue().valeurs(), &vit  = eq_hydraulique->inconnue().valeurs();
  const double epsilon_minimum = eqn_k_Rea->modele_turbulence().get_LeEPS_MIN();
  Modele_Fonc_Realisable_base& mon_modele_fonc = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale,eqn_k_Rea->modele_turbulence()).associe_modele_fonction();
  mon_modele_fonc.Contributions_Sources_BiK(eqn_k_Rea ->zone_dis(),eqn_k_Rea->zone_Cl_dis(),vit,K_Rea,eps_Rea,epsilon_minimum);
  Source_Transport_Realisable_VDF_Elem_base::mettre_a_jour(temps);
}

void Source_Transport_K_Realisable_VDF_Elem::fill_coeff_matrice(const int is_visco_const, const DoubleTab& tab_visco, const DoubleVect& volumes, const DoubleVect& porosite, double& visco, Matrice_Morse& matrice) const
{
  const DoubleTab& K_Rea  = eqn_k_Rea->inconnue().valeurs(), &eps_Rea  = eqn_eps_Rea->inconnue().valeurs();
  for (int c = 0; c < K_Rea.dimension(0); c++)
    {
      if (!is_visco_const) visco = tab_visco(c);
      assert(visco > 0.);
      if (K_Rea(c)>DMINFLOAT)
        matrice(c,c)+=porosite(c)*volumes(c)*eps_Rea(c)/( K_Rea(c) + sqrt( visco*eps_Rea(c) ) );
    }
}
