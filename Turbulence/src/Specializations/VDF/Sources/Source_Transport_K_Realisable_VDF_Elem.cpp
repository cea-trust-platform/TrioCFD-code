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
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////


#include <Source_Transport_K_Realisable_VDF_Elem.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Fluide_base.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Champ_Face.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable_Bicephale.h>
#include <DoubleTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Convection_Diffusion_Temperature_Turbulent.h>
#include <Ref_Transport_K_ou_Eps_Realisable.h>
#include <Constituant.h>
#include <Source_base.h>
#include <Param.h>

Implemente_instanciable(Source_Transport_K_Realisable_VDF_Elem,"Source_Transport_K_Realisable_VDF_P0_VDF",Source_base);


Sortie& Source_Transport_K_Realisable_VDF_Elem::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


Entree& Source_Transport_K_Realisable_VDF_Elem::readOn(Entree& is )
{
  Param param(que_suis_je());
  param.lire_avec_accolades(is);
  return is ;
}



void Source_Transport_K_Realisable_VDF_Elem::associer_pb(const Probleme_base& pb )
{
  eq_hydraulique = pb.equation(0);
  eqn_k_Rea      = ref_cast(Transport_K_ou_Eps_Realisable,equation());
  eqn_eps_Rea    = ref_cast(Transport_K_ou_Eps_Realisable,eqn_eps_Rea.valeur().modele_turbulence().eqn_transp_Eps());
}

void Source_Transport_K_Realisable_VDF_Elem::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF, zone_Cl_dis.valeur());
}

DoubleTab& Source_Transport_K_Realisable_VDF_Elem::ajouter(DoubleTab& resu) const
{
  const Zone_VDF&                                         zone_VDF = ref_cast(Zone_VDF,eq_hydraulique->zone_dis().valeur());
  const Zone_Cl_VDF&                                       zcl_VDF = ref_cast(Zone_Cl_VDF,eq_hydraulique->zone_Cl_dis().valeur());
  const DoubleTab&                                           K_Rea = eqn_k_Rea->inconnue().valeurs();
  const DoubleTab&                                         eps_Rea = eqn_eps_Rea->inconnue().valeurs();
  const Modele_turbulence_hyd_K_Eps_Realisable_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale,eqn_k_Rea->modele_turbulence());
  const DoubleTab&                                      visco_turb = mod_turb.viscosite_turbulente().valeurs();
  const DoubleTab&                                            vit  = eq_hydraulique->inconnue().valeurs();

  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();
  int nb_elem = zone_VDF.nb_elem();
  //  int nb_elem_tot = zone_VDF.nb_elem_tot();
  DoubleTrav P(visco_turb);

  // Rq : la distinction entre zone_cl est importante pour les deux equations pour l imposition des conditions aux limites!!!!!!
  if (axi)
    {
      Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());

      calculer_terme_production_K_BiK_Axi(zone_VDF,vitesse,P,K_Rea,visco_turb);
    }
  else
    {
      Champ_Face& vitesse = ref_cast_non_const(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K_BiK(zone_VDF,zcl_VDF,P,K_Rea,vit,vitesse,visco_turb);
    }

  P.echange_espace_virtuel();

  for (int elem=0; elem<nb_elem; elem++)
    {
      resu(elem) += ( P(elem)-eps_Rea(elem) )*volumes(elem)*porosite_vol(elem);
    }

  return resu;
}

DoubleTab& Source_Transport_K_Realisable_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Source_Transport_K_Realisable_VDF_Elem::mettre_a_jour(double temps)
{
  const Zone_Cl_dis&   zcl_keps  = eqn_k_Rea->zone_Cl_dis();
  const Zone_dis& zone_dis_keps  = eqn_k_Rea ->zone_dis();
  const DoubleTab&        K_Rea  = eqn_k_Rea->inconnue().valeurs();
  const DoubleTab&      eps_Rea  = eqn_eps_Rea->inconnue().valeurs();
  const DoubleTab&          vit  = eq_hydraulique->inconnue().valeurs();

  Modele_turbulence_hyd_K_Eps_Realisable_Bicephale&     mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable_Bicephale,eqn_k_Rea->modele_turbulence());
  Modele_Fonc_Realisable_base&                   mon_modele_fonc = mod_turb.associe_modele_fonction();
  double                                         epsilon_minimum = eqn_k_Rea.valeur().modele_turbulence().get_LeEPS_MIN();

  mon_modele_fonc.Contributions_Sources_BiK(zone_dis_keps,zcl_keps,vit,K_Rea,eps_Rea,epsilon_minimum);

  Calcul_Production_K_VDF::mettre_a_jour(temps);
}

void  Source_Transport_K_Realisable_VDF_Elem::contribuer_a_avec(const DoubleTab& a,  Matrice_Morse& matrice) const
{
  const DoubleTab&        K_Rea  = eqn_k_Rea->inconnue().valeurs();
  const DoubleTab&      eps_Rea  = eqn_eps_Rea->inconnue().valeurs();
  int size=K_Rea.dimension(0);

  const Zone_VDF&   zone_VDF = la_zone_VDF.valeur();
  const DoubleVect& porosite = zone_VDF.porosite_elem();

  const Fluide_base&  fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const Champ_Don&        ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab&           tab_visco = ch_visco_cin->valeurs();

  int is_visco_const=sub_type(Champ_Uniforme,ch_visco_cin.valeur());

  double visco=-1;
  if (is_visco_const)
    {
      visco = std::max(tab_visco(0,0),DMINFLOAT);
    }

  // on implicite le -eps et le -eps^2/k

  const DoubleVect& volumes=zone_VDF.volumes();
  for (int c=0; c<size; c++)
    {
      if (!is_visco_const)
        {
          visco =  tab_visco(c);
        }

      assert(visco>0.);

      // -eps*vol  donne +vol dans la bonne case
      if (K_Rea(c)>DMINFLOAT)

        {
          matrice(c,c)+=porosite(c)*volumes(c)*eps_Rea(c)/( K_Rea(c) + sqrt( visco*eps_Rea(c) ) );
        }
    }

}

