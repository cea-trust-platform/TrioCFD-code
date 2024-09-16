/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Modele_Jones_Launder_Thermique_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Jones_Launder_Thermique_VDF_included
#define Modele_Jones_Launder_Thermique_VDF_included

#include <Modele_Fonc_Bas_Reynolds_Thermique_Base.h>
#include <Transport_Fluctuation_Temperature_W_Bas_Re.h>
#include <TRUSTTabs_forward.h>
#include <TRUST_Ref.h>

#include <Domaine_Cl_dis.h>
class Domaine_Cl_VDF;
class Domaine_VDF;
class Fluide_base;
class Transport_K_Eps_Bas_Reynolds;
class Champ_Don_base;
class Champ_Face_VDF;

class Modele_Jones_Launder_Thermique_VDF : public Modele_Fonc_Bas_Reynolds_Thermique_Base
{

  Declare_instanciable(Modele_Jones_Launder_Thermique_VDF);

public :

  DoubleTab& Calcul_D(DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis&,const DoubleTab&,const DoubleTab&, double) const override;
  DoubleTab& Calcul_E(DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis&,const DoubleTab&,const DoubleTab&,double,const DoubleTab& ) const override;
  DoubleTab& Calcul_F1(DoubleTab&, const Domaine_dis_base&,const DoubleTab&,const DoubleTab&,double,double) const override ;
  DoubleTab& Calcul_F2(DoubleTab&, const Domaine_dis_base&,const DoubleTab&,const DoubleTab&,double,double) const override ;
  DoubleTab& Calcul_F3(DoubleTab&, const Domaine_dis_base&,const DoubleTab&,const DoubleTab&,double,double) const override ;
  DoubleTab& Calcul_F4(DoubleTab&, const Domaine_dis_base&,const DoubleTab&,const DoubleTab&,double,double) const override ;
  DoubleTab& Calcul_Flambda ( DoubleTab&,const Domaine_dis_base&,const DoubleTab&,const DoubleTab&,double,double) const override ;
  Entree& lire(const Motcle&, Entree&);
  void associer_pb(const Probleme_base& ) override;
  void associer(const Domaine_dis_base& , const Domaine_Cl_dis& ) override;
  void mettre_a_jour(double) override;

protected:

  REF(Domaine_VDF) le_dom_VDF;
  REF(Domaine_Cl_VDF) le_dom_Cl_VDF;
  REF(Fluide_base) le_fluide;
  REF(Champ_Inc) la_vitesse_transportante;
  REF(Transport_Fluctuation_Temperature_W_Bas_Re) eq_transport_Fluctu_Temp_Bas_Re;

};

#endif



