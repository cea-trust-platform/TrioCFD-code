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
// File:        Modele_F1F2FMU_unitaire_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_F1F2FMU_unitaire_VDF.h>
#include <Domaine_VDF.h>
#include <Champ_Uniforme.h>

Implemente_instanciable(Modele_F1F2FMU_unitaire_VDF,"Modele_F1F2FMU_unitaire_VDF",Modele_Jones_Launder_VDF);



///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////
// printOn et readOn

Sortie& Modele_F1F2FMU_unitaire_VDF::printOn(Sortie& s ) const
{
  Modele_Jones_Launder_VDF::printOn(s);
  return s;
}

Entree& Modele_F1F2FMU_unitaire_VDF::readOn(Entree& is )
{
  Modele_Jones_Launder_VDF::readOn(is);
  return is;
}

Entree& Modele_F1F2FMU_unitaire_VDF::lire(const Motcle& , Entree& is)
{
  return is;
}

void  Modele_F1F2FMU_unitaire_VDF::associer(const Domaine_dis_base& domaine_dis,
                                            const Domaine_Cl_dis_base& domaine_Cl_dis)
{
}

DoubleTab&  Modele_F1F2FMU_unitaire_VDF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  Fmu=1;
//  Cerr<<Fmu.mp_min_vect()<<" Fmu "<<Fmu.mp_max_vect()<<finl;
  return Fmu;
}
DoubleTab& Modele_F1F2FMU_unitaire_VDF::Calcul_F1( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{
  F1= 1.;
  return F1;
}

DoubleTab& Modele_F1F2FMU_unitaire_VDF::Calcul_F2( DoubleTab& F2, DoubleTab& D, const Domaine_dis_base& domaine_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  F2=1;
  return F2;
}

DoubleTab&  Modele_F1F2FMU_unitaire_VDF::Calcul_Fmu_BiK( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  Fmu=1;
  Cerr<<Fmu.mp_min_vect()<<" Fmu "<<Fmu.mp_max_vect()<<finl;
  return Fmu;
}

DoubleTab& Modele_F1F2FMU_unitaire_VDF::Calcul_F2_BiK( DoubleTab& F2, DoubleTab& D, const Domaine_dis_base& domaine_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco ) const
{
  F2=1;
  return F2;
}

DoubleTab& Modele_F1F2FMU_unitaire_VDF::Calcul_F1_BiK( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const
{
  F1= 1.;
  return F1;
}



