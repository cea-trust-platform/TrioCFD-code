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
// File:        Source_Transport_Fluctuation_Temperature_VDF_Elem.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources/UNUSED
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Transport_Fluctuation_Temperature_VDF_Elem_included
#define Source_Transport_Fluctuation_Temperature_VDF_Elem_included


/*! @brief class Source_Transport_Fluctuation_Temperature_VDF_Elem
 *
 */
#define Ca_DEFAULT 0.8   // Valeurs par defaut des constantes qui interviennent
#define Cb_DEFAULT 2.0   // dans le calcul des termes sources des equations
#define Cc_DEFAULT 1.96  // de transport des fluctuations thermiques.
#define Cd_DEFAULT 0.8

#include <Source_Transport_K_Eps_VDF_Elem.h>
#include <Transport_K_Eps_Bas_Reynolds.h>
#include <Transport_Fluctuation_Temperature.h>
#include <Transport_Flux_Chaleur_Turbulente.h>
#include <TRUSTTabs_forward.h>
#include <TRUST_Ref.h>



class Probleme_base;
class Champ_Don_base;
class Domaine_Cl_VDF;
class Champ_Face_VDF;

// La classe derive de Source_base et peut etre d'un terme source
class Source_Transport_Fluctuation_Temperature_VDF_Elem : public Source_base,
  public Calcul_Production_K_VDF
{
  Declare_instanciable_sans_constructeur(Source_Transport_Fluctuation_Temperature_VDF_Elem);

public :

  void associer_pb(const Probleme_base& ) override;
  inline Source_Transport_Fluctuation_Temperature_VDF_Elem(double ctea = Ca_DEFAULT, double cteb = Cb_DEFAULT,  double ctec = Cc_DEFAULT,  double cted = Cd_DEFAULT);
  DoubleTab& calculer_Prod_uteta_T(const Domaine_VDF&,const Domaine_Cl_VDF&, const DoubleTab&,const  DoubleTab&, DoubleTab&) const;
  DoubleTab& ajouter(DoubleTab& ) const override;
  DoubleTab& calculer(DoubleTab& ) const override;
  void mettre_a_jour(double temps) override { Calcul_Production_K_VDF::mettre_a_jour(temps); }

protected :

  double Ca, Cb, Cc, Cd;
  OBS_PTR(Domaine_VDF) le_dom_VDF;
  OBS_PTR(Domaine_Cl_VDF) le_dom_Cl_VDF;
  OBS_PTR(Equation_base) eq_hydraulique;
  OBS_PTR(Transport_K_Eps_Bas_Reynolds)  mon_eq_transport_K_Eps_Bas_Re_;
  OBS_PTR(Transport_Fluctuation_Temperature) mon_eq_transport_Fluctu_Temp;
  OBS_PTR(Transport_Flux_Chaleur_Turbulente) mon_eq_transport_Flux_Chaleur_Turb_;
  OBS_PTR(Convection_Diffusion_Temperature) eq_thermique;
  OBS_PTR(Champ_Don_base) beta_t;
  OBS_PTR(Champ_Don_base) gravite;
  void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base& ) override;

};

inline Source_Transport_Fluctuation_Temperature_VDF_Elem::
Source_Transport_Fluctuation_Temperature_VDF_Elem(double ctea, double cteb,  double ctec,  double cted)
  : Ca(ctea), Cb(cteb), Cc(ctec), Cd(cted) {}

#endif
