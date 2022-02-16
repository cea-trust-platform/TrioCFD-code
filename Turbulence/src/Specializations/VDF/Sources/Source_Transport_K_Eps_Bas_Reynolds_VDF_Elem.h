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
// File:        Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem_included
#define Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem_included


#define C11_DEFAULT   1.55 // 1.44   // Valeurs par defaut des constantes qui interviennent
#define C21_DEFAULT  2. // 1.92   // dans l'equation de k_esp

//
// .DESCRIPTION class Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem
//
// .SECTION voir aussi

#include <Source_Transport_K_Eps_VDF_Elem.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Ref_Transport_Flux_Chaleur_Turbulente.h>
#include <Modele_Fonc_Bas_Reynolds.h>
#include <Ref_Zone_dis.h>
#include <Zone_Cl_dis.h>
#include <Ref_Transport_K_Eps_Bas_Reynolds.h>

class Probleme_base;
class Champ_Don_base;
class DoubleVect;
class DoubleTab;
class Champ_Face;

class Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem : public Source_base,public Calcul_Production_K_VDF


{
  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem);

public :

  inline Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem(double cte1 = C11_DEFAULT,
                                                      double cte2 = C21_DEFAULT );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  inline Modele_Fonc_Bas_Reynolds&  associe_modele_fonc();
  inline const Modele_Fonc_Bas_Reynolds&  associe_modele_fonc() const;
  virtual void associer_pb(const Probleme_base& );
  void mettre_a_jour(double temps)
  {
    Calcul_Production_K_VDF::mettre_a_jour(temps);
  }

protected :

  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  double C1;
  double C2;

  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;
  REF(Transport_K_Eps_Bas_Reynolds) eqn_keps_bas_re;
  REF(Transport_Flux_Chaleur_Turbulente) eqn_flux_chaleur;
  REF(Equation_base) eq_hydraulique;

};




//////////////////////////////////////////////////////////////////////////////
//
//   Fonctions inline de la classe Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem::
Source_Transport_K_Eps_Bas_Reynolds_VDF_Elem(double cte1,double cte2)

  : C1(cte1), C2(cte2) {}




#endif
