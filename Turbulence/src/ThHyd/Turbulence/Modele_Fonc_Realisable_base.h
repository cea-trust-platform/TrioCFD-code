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
// File:        Modele_Fonc_Realisable_base.h
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Fonc_Realisable_base_included
#define Modele_Fonc_Realisable_base_included


#include <Ref_Champ_base.h>
#include <Champ_Fonc.h>
#include <Ref_Fluide_base.h>
#include <Ref_Champ_Inc.h>
#include <Ref_Equation_base.h>
#include <Ref_Champ_Don.h>
#include <Champs_compris.h>
#include <Champs_compris_interface.h>

class Motcle;
class Zone_dis;
class Zone_Cl_dis;
/* class Champ_Don_base; */
#include <Champ_Don.h>

class Equation_base;
class Probleme_base;
class Discretisation_base;
class Champ_base;

class Modele_Fonc_Realisable_base : public Champs_compris_interface, public Objet_U
{

  Declare_base(Modele_Fonc_Realisable_base);

public:

  inline const Equation_base& equation() const;
  inline  Equation_base& equation();
  inline const DoubleTab& get_S() const;
  inline  DoubleTab& get_S();
  inline const DoubleTab& get_Cmu( void ) const;
  inline DoubleTab& get_Cmu( void );
  inline const DoubleTab& get_C1( void ) const;
  inline DoubleTab& get_C1( void );


  virtual int preparer_calcul();
  virtual void mettre_a_jour(double ) =0;
  virtual void discretiser();
  virtual void completer();
  virtual void associer_pb(const Probleme_base& ) ;
  virtual void associer(const Zone_dis& , const Zone_Cl_dis& )= 0;
  virtual int sauvegarder(Sortie& ) const;
  virtual int reprendre(Entree& );
  virtual int Calcul_is_Reynolds_stress_isotrope() const;
  virtual int Calcul_is_Cmu_constant() const;

  virtual void Calcul_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse) =0;
  virtual void Calcul_C1(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) =0;
  virtual void Calcul_Cmu_et_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN)  =0;
  virtual void Contributions_Sources(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) =0;
  virtual void Contributions_Sources_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                           const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) =0;

  Entree& lire(const Motcle&, Entree&);

  //Methodes de l interface des champs postraitables
  /////////////////////////////////////////////////////
  virtual void creer_champ(const Motcle& motlu);
  virtual const Champ_base& get_champ(const Motcle& nom) const;
  virtual void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const;
  /////////////////////////////////////////////////////

  void associer_eqn(const Equation_base& );

protected :

  DoubleTab S_;
  DoubleTab Cmu_;
  DoubleTab C1_;

  REF(Equation_base) mon_equation;

  REF(Fluide_base) le_fluide;
  REF(Champ_Inc) la_vitesse_transportante;
  REF(Equation_base) eq_hydraulique;
  REF(Champ_Don) visco_;

  Nom nom_fic;
  Champ_Fonc BR_wall_length_;
  int is_Cmu_constant_;
  int is_Reynolds_stress_isotrope_;

private :

  Champs_compris champs_compris_;

};

//
// fonctions inline
//



inline const Equation_base& Modele_Fonc_Realisable_base::equation() const
{
  return mon_equation.valeur();
}

inline Equation_base& Modele_Fonc_Realisable_base::equation()
{
  return mon_equation.valeur();
}

inline const DoubleTab& Modele_Fonc_Realisable_base::get_S() const
{
  return S_;
}

inline DoubleTab& Modele_Fonc_Realisable_base::get_S()
{
  return S_;
}

inline DoubleTab& Modele_Fonc_Realisable_base::get_Cmu( void )
{
  return Cmu_;
}

inline const DoubleTab& Modele_Fonc_Realisable_base::get_Cmu( void ) const
{
  return Cmu_;
}

inline DoubleTab& Modele_Fonc_Realisable_base::get_C1( void )
{
  return C1_;
}

inline const DoubleTab& Modele_Fonc_Realisable_base::get_C1( void ) const
{
  return C1_;
}


#endif



