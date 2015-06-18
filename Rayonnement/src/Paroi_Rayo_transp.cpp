/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Paroi_Rayo_transp.cpp
// Directory:   $TRUST_ROOT/src/Rayonnement
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_Rayo_transp.h>
#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Milieu_base.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Front_VF.h>
#include <Fluide_Quasi_Compressible.h>
Implemente_base(Paroi_Rayo_transp,"Paroi_Rayo_transp",Neumann_paroi);

//
// printOn et readOn

Sortie& Paroi_Rayo_transp::printOn(Sortie& s ) const
{
  return s;
}

Entree& Paroi_Rayo_transp::readOn(Entree& is )
{
  return is ;
}

double Paroi_Rayo_transp::flux_impose(int i) const
{
  const Milieu_base& mil=ma_zone_cl_dis->equation().milieu();
  double d_Cp = 1;
  double d_rho= 1;
  // si on est en QC on n'a pas a diviser le flux par rho*cp
  if (!sub_type(Fluide_Quasi_Compressible,mil))
    {
      const Champ_Uniforme& rho=ref_cast(Champ_Uniforme,mil.masse_volumique().valeur());
      const Champ_Uniforme& Cp=ref_cast(Champ_Uniforme,mil.capacite_calorifique().valeur());
      d_Cp = Cp(0,0);
      d_rho= rho(0,0);
    }
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF,frontiere_dis());
  int ndeb = la_frontiere_VF.num_premiere_face();
  double flux_radia=le_modele_rayo->flux_radiatif(i+ndeb);
  if (le_champ_front.valeurs().size()==1)
    return (le_champ_front(0,0)-flux_radia)/(d_rho*d_Cp);
  else if (le_champ_front.valeurs().dimension(1)==1)
    return (le_champ_front(i,0)-flux_radia)/(d_rho*d_Cp);
  else
    Cerr << "Paroi_Rayo_transp::flux_impose erreur" << finl;
  exit();
  return 0.;
}

double Paroi_Rayo_transp::flux_impose(int i,int j) const
{
  Cerr << "Paroi_Rayo_transp::flux_impose(i,j)  erreur" << finl;
  exit();
  return 0;

}
