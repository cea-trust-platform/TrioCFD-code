/****************************************************************************
* Copyright (c) 2020, CEA
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
// File      : Operateur_Conv_sensibility.cpp
// Directory : $Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_Conv_sensibility.h>

Implemente_base( Operateur_Conv_sensibility, "Operateur_Conv_sensibility", Operateur_Conv_base ) ;
// XD convection_sensibility convection_deriv sensibility 0 A convective scheme for the sensibility problem.
// XD  attr opconv bloc_convection opconv 0 Choice between: amont and muscl   NL2 Example: convection {  Sensibility { amont } }

Sortie& Operateur_Conv_sensibility::printOn( Sortie& os ) const
{
  return os;
}

Entree& Operateur_Conv_sensibility::readOn( Entree& is )
{
  return is;
}
DoubleTab& Operateur_Conv_sensibility::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  op_conv.ajouter(inco, resu);
  return resu;
}
void Operateur_Conv_sensibility::associer_vitesse(const Champ_base& vit)
{
  la_vitesse = ref_cast(Champ_Inc_base,vit);
}
void Operateur_Conv_sensibility::associer(const Domaine_dis_base& zdis,
                                          const Domaine_Cl_dis_base& zcl_dis,
                                          const Champ_Inc& inco)
{
  dom=inco->domaine();
  op_conv.l_op_base().associer(zdis, zcl_dis, inco);
}
DoubleTab& Operateur_Conv_sensibility::calculer(const DoubleTab& inco, DoubleTab& resu) const
{
  op_conv.calculer(inco, resu);
  return resu;
}


