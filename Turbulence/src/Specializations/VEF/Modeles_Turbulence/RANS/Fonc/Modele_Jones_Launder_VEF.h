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
// File:        Modele_Jones_Launder_VEF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Jones_Launder_VEF_included
#define Modele_Jones_Launder_VEF_included

#include <Modele_Fonc_Bas_Reynolds_Base.h>
#include <Domaine_Cl_dis.h>
#include <TRUST_Ref.h>

class Domaine_dis;
class Domaine_Cl_dis;
#include <TRUSTTabs_forward.h>
class Domaine_Cl_VEF;
class Domaine_VEF;

class Modele_Jones_Launder_VEF : public Modele_Fonc_Bas_Reynolds_Base
{

  Declare_instanciable(Modele_Jones_Launder_VEF);

public :


  DoubleTab& Calcul_D(DoubleTab&, const Domaine_dis&, const Domaine_Cl_dis&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const override;
  DoubleTab& Calcul_E(DoubleTab&,const Domaine_dis&,const Domaine_Cl_dis&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const override ;

//  virtual DoubleTab& Calcul_F1(DoubleTab&, const Domaine_dis& ) const ;
  DoubleTab& Calcul_F1( DoubleTab& F1, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const override;
  DoubleTab& Calcul_F2(DoubleTab&, DoubleTab&,const Domaine_dis&,const DoubleTab&,const Champ_base&) const override ;
  DoubleTab& Calcul_Fmu ( DoubleTab&,const Domaine_dis&,const Domaine_Cl_dis&,const DoubleTab&,const Champ_Don& )const override ;
  Entree& lire(const Motcle&, Entree&);
  // void associer_pb(const Probleme_base& );
  void associer(const Domaine_dis& , const Domaine_Cl_dis& ) override;
  void mettre_a_jour(double) override;

  DoubleTab& Calcul_D_BiK(DoubleTab&, const Domaine_dis&, const Domaine_Cl_dis&,const DoubleTab&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const override;
  DoubleTab& Calcul_E_BiK(DoubleTab&,const Domaine_dis&,const Domaine_Cl_dis&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const override ;
//
  DoubleTab& Calcul_F1_BiK( DoubleTab& F1, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const override;
  DoubleTab& Calcul_F2_BiK(DoubleTab&, DoubleTab&,const Domaine_dis&,const DoubleTab&,const DoubleTab&,const Champ_base&) const override ;
  DoubleTab& Calcul_Fmu_BiK ( DoubleTab&,const Domaine_dis&,const Domaine_Cl_dis&,const DoubleTab&,const DoubleTab&,const Champ_Don& )const override ;

protected:

  REF(Domaine_VEF) le_dom_VEF;
  REF(Domaine_Cl_VEF) le_dom_Cl_VEF;
  /*   REF(Fluide_Incompressible) le_fluide; */
  /*   REF(Champ_Inc) la_vitesse_transportante; */
  /*   REF(Transport_K_Eps_Bas_Reynolds) eq_transport_K_Eps_Bas_Re; */
  /*   REF(Champ_Don) visco; */
};

#endif



