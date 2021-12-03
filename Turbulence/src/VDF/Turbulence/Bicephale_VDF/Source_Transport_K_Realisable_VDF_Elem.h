/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Source_Transport_K_Realisable_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Transport_K_Realisable_VDF_Elem_included
#define Source_Transport_K_Realisable_VDF_Elem_included


#define C21_DEFAULT_K_EPS_REALISABLE  1.9   // dans l'equation de k_eps
#define C3_DEFAULT_K_EPS_REALISABLE 1.0    // de transport de K et Eps source: Chabard et N3S

//
// .DESCRIPTION class Source_Transport_K_Realisable_VDF_Elem
//
// .SECTION voir aussi

#include <Source_Transport_K_Eps_VDF_Elem.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Ref_Transport_Flux_Chaleur_Turbulente.h>
#include <Modele_Fonc_Realisable_base.h>
#include <Ref_Zone_dis.h>
#include <Zone_Cl_dis.h>
#include <Ref_Transport_K_ou_Eps_Realisable.h>
#include <Ref_Convection_Diffusion_Temperature.h>
#include <Ref_Convection_Diffusion_Concentration.h>
#include <Calcul_Production_K_VDF.h>
#include <Source_base.h>

class Probleme_base;
class Champ_Don_base;
class DoubleVect;
class DoubleTab;
class Champ_Face;
class Zone_dis;
class Zone_Cl_dis;
class Zone_Cl_VDF;
class Champ_Face;

class Source_Transport_K_Realisable_VDF_Elem : public Source_base,public Calcul_Production_K_VDF


{
  Declare_instanciable(Source_Transport_K_Realisable_VDF_Elem);

public :

  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  inline Modele_Fonc_Realisable_base&  associe_modele_fonc();
  inline const Modele_Fonc_Realisable_base&  associe_modele_fonc() const;
  virtual void associer_pb(const Probleme_base& );
  void mettre_a_jour(double temps);
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const ;

protected :

  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;
  REF(Transport_K_ou_Eps_Realisable) eqn_k_Rea;
  REF(Transport_K_ou_Eps_Realisable) eqn_eps_Rea;
  REF(Transport_Flux_Chaleur_Turbulente) eqn_flux_chaleur;
  REF(Equation_base) eq_hydraulique;

};


//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Realisable_anisotherme_VDF_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de Navier_Stokes
// sont couplees a l'equation de la thermique
// On suppose que le coefficient de variation de la masse volumique
// du fluide en fonction de ce scalaire est un coefficient uniforme.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Realisable_anisotherme_VDF_Elem :
  public Source_Transport_K_Realisable_VDF_Elem
{

  Declare_instanciable(Source_Transport_K_Realisable_anisotherme_VDF_Elem);

public:

  virtual void associer_pb(const Probleme_base& );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;

protected:
  REF(Convection_Diffusion_Temperature) eq_thermique;
  REF(Champ_Don) beta_t;
  REF(Champ_Don_base) gravite;

};

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Realisable_aniso_concen_VDF_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de Navier_Stokes
// sont couplees a l'equation d'un transport d'un constituant
// On suppose que le coefficient de variation de la masse volumique
// du fluide en fonction de ce scalaire est un coefficient uniforme.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Realisable_aniso_concen_VDF_Elem :
  public Source_Transport_K_Realisable_VDF_Elem
{

  Declare_instanciable(Source_Transport_K_Realisable_aniso_concen_VDF_Elem);

public:

  virtual void associer_pb(const Probleme_base& );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;

protected:
  REF(Convection_Diffusion_Concentration) eq_concentration;
  REF(Champ_Don) beta_c;
  REF(Champ_Don_base) gravite;

};

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Realisable_aniso_therm_concen_VDF_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de
// Navier_Stokes sont couplees a l'equation de convection diffusion
// d'une concentration et a l'equation de la thermique
// Les champs beta_t et beta_c sont uniformes
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Realisable_aniso_therm_concen_VDF_Elem :
  public Source_Transport_K_Realisable_VDF_Elem
{

  Declare_instanciable(Source_Transport_K_Realisable_aniso_therm_concen_VDF_Elem);

public:

  virtual void associer_pb(const Probleme_base& );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;

protected:
  REF(Convection_Diffusion_Temperature) eq_thermique;
  REF(Convection_Diffusion_Concentration) eq_concentration;
  REF(Champ_Don) beta_t;
  REF(Champ_Don) beta_c;
  REF(Champ_Don_base) gravite;

};

#endif
