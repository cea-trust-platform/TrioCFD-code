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
// File:        Source_Transport_K_Eps_VDF_Elem.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/41
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_Eps_VDF_Elem.h>
#include <Transport_K_Eps.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Concentration.h>
#include <Modele_turbulence_scal_base.h>
#include <Probleme_base.h>
#include <IntTrav.h>
#include <Champ_Uniforme.h>
#include <Zone_VDF.h>
#include <Champ_Face.h>
#include <Zone_Cl_VDF.h>
#include <Fluide_Quasi_Compressible.h>
#include <Debog.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <DoubleTrav.h>
#include <Pb_Hydraulique_Turbulent.h>
#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Pb_Thermohydraulique_Turbulent_QC.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Param.h>
#include <Constituant.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_Eps_VDF_Elem,"Source_Transport_K_Eps_VDF_P0_VDF",Source_base);


Sortie& Source_Transport_K_Eps_VDF_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je() ;
}


Entree& Source_Transport_K_Eps_VDF_Elem::readOn(Entree& is)
{
  const Probleme_base& problem = mon_equation->probleme();
  if (!sub_type(Pb_Hydraulique_Turbulent,problem) && !sub_type(Pb_Thermohydraulique_Turbulent_QC,problem)) error(que_suis_je(),problem.que_suis_je());
  Param param(que_suis_je());
  param.ajouter("C1_eps", &C1);
  param.ajouter("C2_eps", &C2);
  param.lire_avec_accolades(is);
  Cerr << "C1_eps = " << C1 << finl;
  Cerr << "C2_eps = " << C2 << finl;
  return is ;
}



void Source_Transport_K_Eps_VDF_Elem::associer_zones(const Zone_dis& zone_dis,
                                                     const Zone_Cl_dis&  )
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
}

// remplit les references eq_hydraulique et mon_eq_transport_K_Eps
void Source_Transport_K_Eps_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  eq_hydraulique = pb.equation(0);
  mon_eq_transport_K_Eps = ref_cast(Transport_K_Eps,equation());
}

DoubleTab& Source_Transport_K_Eps_VDF_Elem::ajouter(DoubleTab& resu) const
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_Cl_VDF& zcl_VDF = ref_cast(Zone_Cl_VDF,eq_hydraulique->zone_Cl_dis().valeur());
  const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue().valeurs();
  const DoubleTab& visco_turb = mon_eq_transport_K_Eps->modele_turbulence().viscosite_turbulente().valeurs();
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = zone_VDF.porosite_elem();

  int nb_elem = zone_VDF.nb_elem();

  // Ajout d'un espace virtuel au tableu P
  DoubleVect P;

  zone_VDF.zone().creer_tableau_elements(P);

  if (axi)
    {
      const Champ_Face& vitesse = ref_cast(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K_Axi(zone_VDF,vitesse,P,K_eps,visco_turb);
    }
  else
    {
      const Champ_Face& vitesse = ref_cast(Champ_Face,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K(zone_VDF,zcl_VDF,P,K_eps,vit,vitesse,visco_turb);
    }
  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc=ref_cast(Modele_turbulence_hyd_K_Eps,mon_eq_transport_K_Eps->modele_turbulence()).associe_modele_fonction();
  int is_modele_fonc=(mon_modele_fonc.non_nul());
  //is_modele_fonc=0;
  if (is_modele_fonc)
    {

      DoubleTab& D=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("D").valeurs());
      DoubleTab& E=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("E").valeurs());
      DoubleTab& F1=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("F1").valeurs());
      DoubleTab& F2=ref_cast_non_const(DoubleTab,mon_modele_fonc.valeur().get_champ("F2").valeurs());
      const Fluide_base& fluide=ref_cast(Fluide_base,eq_hydraulique.valeur().milieu());
      const Champ_Don& ch_visco_cin = fluide.viscosite_cinematique();
      // const DoubleTab& tab_visco = ch_visco_cin->valeurs();
      const Zone_Cl_dis& zcl_keps=mon_eq_transport_K_Eps->zone_Cl_dis();
      const Zone_dis& zone_dis_keps =mon_eq_transport_K_Eps->zone_dis();
      mon_modele_fonc.Calcul_D(D,zone_dis_keps,zcl_keps,vit,K_eps,ch_visco_cin);
      mon_modele_fonc.Calcul_E(E,zone_dis_keps,zcl_keps,vit,K_eps,ch_visco_cin,visco_turb);

      D.echange_espace_virtuel();
      E.echange_espace_virtuel();
      const Champ_base& ch_visco_cin_ou_dyn =ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();
      DoubleTab P_tab;
      P_tab.ref(P);
      mon_modele_fonc.Calcul_F1(F1,zone_dis_keps,zcl_keps, P_tab, K_eps,ch_visco_cin_ou_dyn);


      mon_modele_fonc.Calcul_F2(F2,D,zone_dis_keps,K_eps, ch_visco_cin_ou_dyn  );
      Debog::verifier("D",D);
      Debog::verifier("E",E);
      Debog::verifier("F2",F2);
      Debog::verifier("F1",F1);
      Debog::verifier("avt",resu);
      for (int elem=0; elem<nb_elem; elem++)
        {
          //if (K_eps(elem,0) > 1.e-10 && K_eps(elem,1) > 1.e-10)
          {
            // resu(elem,0) += (P(elem)-eps_sur_k(elem)*K_eps(elem,0))*volumes(elem)*porosite_vol(elem)-D(elem);
            resu(elem,0) += (P(elem)-K_eps(elem,1)-D(elem))*volumes(elem)*porosite_vol(elem);
            // on en profite pour diviser D  par le volume

            //assert(eps_sur_k(elem)==K_eps(elem,1)/K_eps(elem,0));
            //   if (K_eps(elem,0) >= 10.e-10)
            {
              double coef=1.;
              //if (K_eps(elem,0)<2e-2)  coef=0;
              //  D(elem)/=volumes(elem);

              //ving directory `/users/fauchet/VUES/168/Baltik/BR/build/src/exec_opt'
              //resu(elem,1) += ((C1*P(elem)*F1(elem)- C2*F2(elem)*coef*(K_eps(elem,1)-0*D(elem)))*eps_sur_k(elem)+E(elem))*volumes(elem)*porosite_vol(elem);
              resu(elem,1) += ((C1*P(elem)*F1(elem)- C2*F2(elem)*coef*(K_eps(elem,1)))*K_eps(elem,1)/(K_eps(elem,0)+DMINFLOAT)+E(elem))*volumes(elem)*porosite_vol(elem);
            }
          }
        }
      Debog::verifier("ap",resu);

      // int elem=0;
      //Cerr<<" ici "<<D(elem)<< " "<<K_eps(elem,1)<<finl;
    }
  else
    {
      double LeK_MIN = mon_eq_transport_K_Eps->modele_turbulence().get_LeK_MIN();

      for (int elem=0; elem<nb_elem; elem++)
        {
          resu(elem,0) += (P(elem)-K_eps(elem,1))*volumes(elem)*porosite_vol(elem);

          if (K_eps(elem,0) >= LeK_MIN)
            resu(elem,1) += (C1*P(elem)- C2*K_eps(elem,1))*volumes(elem)*porosite_vol(elem)*K_eps(elem,1)/K_eps(elem,0);

        }
    }

  //Debog::verifier("Source_Transport_K_Eps_VDF_Elem::ajouter resu",resu);
  resu.echange_espace_virtuel();
  return resu;
}

DoubleTab& Source_Transport_K_Eps_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu =0.;
  return ajouter(resu);
}

void  Source_Transport_K_Eps_VDF_Elem::contribuer_a_avec(const DoubleTab& a,  Matrice_Morse& matrice) const
{

  const DoubleTab& val=equation().inconnue().valeurs();
  int size=val.dimension(0);

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const DoubleVect& porosite = zone_VDF.porosite_elem();
  // on implicite le -eps et le -eps^2/k

  const Modele_Fonc_Bas_Reynolds& mon_modele_fonc=ref_cast(Modele_turbulence_hyd_K_Eps,mon_eq_transport_K_Eps->modele_turbulence()).associe_modele_fonction();
  int is_modele_fonc=(mon_modele_fonc.non_nul());
  //is_modele_fonc=0;
  DoubleTab F2;
  if (is_modele_fonc)
    {

      DoubleTrav D(0);
      F2.resize(val.dimension_tot(0));
      const Zone_dis& zone_dis_keps =mon_eq_transport_K_Eps->zone_dis();

      const Champ_base& ch_visco_cin_ou_dyn =ref_cast(Op_Diff_K_Eps_base, equation().operateur(0).l_op_base()).diffusivite();

      mon_modele_fonc.Calcul_F2(F2,D,zone_dis_keps,val, ch_visco_cin_ou_dyn  );

    }

  {
    const DoubleVect& volumes=zone_VDF.volumes();
    for (int c=0; c<size; c++)
      {
        // -eps*vol  donne +vol dans la bonne case
        if (val(c,0)>DMINFLOAT)

          {

            double coef_k=porosite(c)*volumes(c)*val(c,1)/val(c,0);
            matrice(c*2,c*2)+=coef_k;

            double coef_eps=C2*coef_k;
            if (is_modele_fonc) coef_eps*=F2(c);
            matrice(c*2+1,c*2+1)+=coef_eps;
          }
      }
  }
}


void Source_Transport_K_Eps_VDF_Elem::mettre_a_jour(double temps)
{

  Calcul_Production_K_VDF::mettre_a_jour(temps);

}
