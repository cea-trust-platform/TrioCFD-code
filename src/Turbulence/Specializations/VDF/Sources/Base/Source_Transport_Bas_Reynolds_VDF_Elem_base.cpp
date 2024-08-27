/****************************************************************************
* Copyright (c) 2022, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_Bas_Reynolds_VDF_Elem_base.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VDF/Sources/new
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Bas_Reynolds_VDF_Elem_base.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>
#include <Domaine_Cl_VDF.h>
#include <Fluide_base.h>
#include <TRUSTTrav.h>
#include <Champ_Face_VDF.h>

Implemente_base_sans_constructeur( Source_Transport_Bas_Reynolds_VDF_Elem_base, "Source_Transport_Bas_Reynolds_VDF_Elem_base", Source_Transport_VDF_Elem_base ) ;

Sortie& Source_Transport_Bas_Reynolds_VDF_Elem_base::printOn( Sortie& os ) const { return os << que_suis_je(); }
Entree& Source_Transport_Bas_Reynolds_VDF_Elem_base::readOn( Entree& is ) { return Source_Transport_VDF_Elem_base::readOn(is); }

void Source_Transport_Bas_Reynolds_VDF_Elem_base::associer_pb(const Probleme_base& pb )
{
  Source_Transport_VDF_Elem_base::associer_pb(pb);
  eqn_keps_bas_re = ref_cast(Transport_K_Eps_Bas_Reynolds,equation());
}

// XXX : Elie Saikali
// a revoir si ca donne des soucis pour Source_Transport_K_Eps_Bas_Reynolds_W_VDF_Elem
// pour le moment j'ai fait comme Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem
void Source_Transport_Bas_Reynolds_VDF_Elem_base::ajouter_blocs(matrices_t matrices, DoubleTab& resu, const tabs_t& semi_impl) const
{
  const DoubleTab& K_eps_Bas_Re = eqn_keps_bas_re->inconnue()->valeurs();
  const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds,eqn_keps_bas_re->modele_turbulence());
  const DoubleTab& visco_turb = mod_turb.viscosite_turbulente()->valeurs();
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc = mod_turb.associe_modele_fonction();
  const Fluide_base& fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& vit = eq_hydraulique->inconnue()->valeurs();
  const Domaine_Cl_dis& zcl_keps=eqn_keps_bas_re->domaine_Cl_dis();
  const Domaine_dis& domaine_dis_keps =eqn_keps_bas_re ->domaine_dis();
  Champ_Face_VDF& vitesse = ref_cast_non_const(Champ_Face_VDF,eq_hydraulique->inconnue().valeur());
  const int nb_elem = le_dom_VDF->nb_elem();

  DoubleTrav P(visco_turb), D(visco_turb), E(visco_turb), F1(nb_elem), F2(nb_elem);
  mon_modele_fonc.Calcul_D(D,domaine_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin);
  D.echange_espace_virtuel();
  mon_modele_fonc.Calcul_E(E,domaine_dis_keps,zcl_keps,vit,K_eps_Bas_Re,ch_visco_cin,visco_turb);
  E.echange_espace_virtuel();
  mon_modele_fonc.Calcul_F2(F2,D,domaine_dis_keps,K_eps_Bas_Re, ch_visco_cin);

  // Rq : la distinction entre domaine_cl est importante pour les deux equations pour l imposition des conditions aux limites!!!!!!
  if (axi) calculer_terme_production_K_Axi(le_dom_VDF.valeur(),vitesse,P,K_eps_Bas_Re,visco_turb);
  else calculer_terme_production_K(le_dom_VDF.valeur(),le_dom_Cl_VDF.valeur(),P,K_eps_Bas_Re,vit,vitesse,visco_turb);
  P.echange_espace_virtuel();

  fill_resu_bas_reyn(P,D,E,F1,F2,resu);
  resu.echange_espace_virtuel();
}
