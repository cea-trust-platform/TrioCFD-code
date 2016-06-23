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
// File:        Source_rayo_semi_transp_QC_VEF_P1NC.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VEF
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_rayo_semi_transp_QC_VEF_P1NC.h>
#include <Zone_VEF.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Champ_Uniforme.h>
#include <Debog.h>

Implemente_instanciable(Source_rayo_semi_transp_QC_VEF_P1NC,"Source_rayo_semi_transp_QC_VEF_P1NC",Source_rayo_semi_transp_base);



Sortie& Source_rayo_semi_transp_QC_VEF_P1NC::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Source_rayo_semi_transp_QC_VEF_P1NC::readOn(Entree& s )
{
  return s;
}

void Source_rayo_semi_transp_QC_VEF_P1NC::associer_modele_rayo(Modele_rayo_semi_transp& mod)
{
  le_modele_ = mod;
  Source_rayo_semi_transp_base& source=ref_cast(Source_rayo_semi_transp_base,le_source_rayo.valeur());
  source.associer_modele_rayo(mod);
}

DoubleTab& Source_rayo_semi_transp_QC_VEF_P1NC::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}


DoubleTab& Source_rayo_semi_transp_QC_VEF_P1NC::ajouter(DoubleTab& resu) const
{
  Debog::verifier("Source_rayo_semi_transp_QC_VEF_P1NC::ajouter resu debut ",resu);
  DoubleTab resu_tmp(resu);
  le_source_rayo.calculer(resu_tmp);
  // C'est la que l'on multiplie le terme source par rho*cp

  const Zone_VEF& zvef = ref_cast(Zone_VEF,equation().zone_dis().valeur());
  int nb_faces = zvef.nb_faces();
  const Milieu_base& fluide = equation().milieu();
  const DoubleTab& rho = fluide.masse_volumique().valeurs();
  const DoubleTab& cp = fluide.capacite_calorifique().valeurs();

  int face;
  for (face=0; face<nb_faces; face++)
    {
      double val_rho;
      assert(fluide.masse_volumique().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide.masse_volumique().valeur()))
        val_rho = rho(0,0);
      else if (rho.nb_dim() == 1)
        val_rho = rho(face);
      else
        val_rho = rho(face,0);

      double val_cp;
      assert(fluide.capacite_calorifique().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide.capacite_calorifique().valeur()))
        val_cp = cp(0,0);
      else
        ////val_cp = cp(face,0);
        if (cp.nb_dim() == 1)
          val_cp = cp(face);
        else
          val_cp = cp(face,0);

      resu_tmp(face) *= val_rho*val_cp;
    }

  resu+=resu_tmp;
  resu.echange_espace_virtuel();
  Debog::verifier("Source_rayo_semi_transp_QC_VEF_P1NC::ajouter resu fin ",resu);
  return resu;
}


void Source_rayo_semi_transp_QC_VEF_P1NC::associer_pb(const Probleme_base& pb)
{
  ;
}


void Source_rayo_semi_transp_QC_VEF_P1NC::associer_zones(const Zone_dis& zone_dis,
                                                         const Zone_Cl_dis& zone_Cl_dis)
{
  ;
}

void Source_rayo_semi_transp_QC_VEF_P1NC::mettre_a_jour(double temps)
{
  le_source_rayo.mettre_a_jour(temps);
}

void Source_rayo_semi_transp_QC_VEF_P1NC::completer()
{
  Source_base::completer();
  le_source_rayo.typer_direct("Source_rayo_semi_transp_VEF_P1NC");
  le_source_rayo->associer_eqn(equation());
  le_source_rayo.completer();
}
