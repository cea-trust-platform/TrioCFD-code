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
// File:        Modele_Launder_Sharma_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Launder_Sharma_VEF_included
#define Modele_Launder_Sharma_VEF_included

#include <Modele_Lam_Bremhorst_VEF.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>
#include <Zone_Cl_dis.h>
#include <Param.h>

class Zone_dis;
class Zone_Cl_dis;
#include <TRUSTTabs_forward.h>
#include <TRUSTTabs_forward.h>
class Zone_Cl_VEF;

class Modele_Launder_Sharma_VEF : public Modele_Lam_Bremhorst_VEF
{

  Declare_instanciable(Modele_Launder_Sharma_VEF);

public :

  DoubleTab& Calcul_D(DoubleTab&, const Zone_dis&, const Zone_Cl_dis&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const override;
  DoubleTab& Calcul_E(DoubleTab&,const Zone_dis&,const Zone_Cl_dis&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const override ;

  //  virtual DoubleTab& Calcul_F1(DoubleTab&, const Zone_dis& ) const ;
  DoubleTab& Calcul_F1( DoubleTab& F1, const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& P,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const override;
  DoubleTab& Calcul_F2(DoubleTab&, DoubleTab&,const Zone_dis&,const DoubleTab&,const Champ_base&) const override ;
  DoubleTab& Calcul_Fmu ( DoubleTab&,const Zone_dis&,const Zone_Cl_dis&,const DoubleTab&,const Champ_Don& )const override ;
  void lire_distance_paroi( ) override;

  DoubleTab& Calcul_D_BiK(DoubleTab&, const Zone_dis&, const Zone_Cl_dis&,const DoubleTab&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const override;
  DoubleTab& Calcul_E_BiK(DoubleTab&,const Zone_dis&,const Zone_Cl_dis&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const override ;
  DoubleTab& Calcul_F1_BiK( DoubleTab& F1, const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& P,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const override;
  DoubleTab& Calcul_F2_BiK(DoubleTab&, DoubleTab&,const Zone_dis&,const DoubleTab&,const DoubleTab&,const Champ_base&) const override ;
  DoubleTab& Calcul_Fmu_BiK ( DoubleTab&,const Zone_dis&,const Zone_Cl_dis&,const DoubleTab&,const DoubleTab&,const Champ_Don& )const override ;

protected:

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zone_Cl_VEF
  ;
};

#endif

