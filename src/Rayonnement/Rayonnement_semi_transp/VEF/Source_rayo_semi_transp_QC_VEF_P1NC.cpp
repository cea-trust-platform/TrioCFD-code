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
#include <Domaine_VEF.h>
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
  return le_source_rayo.ajouter(resu);
}


void Source_rayo_semi_transp_QC_VEF_P1NC::associer_pb(const Probleme_base& pb)
{
  ;
}


void Source_rayo_semi_transp_QC_VEF_P1NC::associer_domaines(const Domaine_dis_base& domaine_dis,
                                                            const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  ;
}

void Source_rayo_semi_transp_QC_VEF_P1NC::mettre_a_jour(double temps)
{
  le_source_rayo->mettre_a_jour(temps);
}

void Source_rayo_semi_transp_QC_VEF_P1NC::completer()
{
  Source_base::completer();
  le_source_rayo.typer_direct("Source_rayo_semi_transp_VEF_P1NC");
  le_source_rayo->associer_eqn(equation());
  le_source_rayo->completer();
}
