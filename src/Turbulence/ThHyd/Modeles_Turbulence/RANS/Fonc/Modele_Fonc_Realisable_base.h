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
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Fonc_Realisable_base_included
#define Modele_Fonc_Realisable_base_included


#include <Champ_Fonc.h>
#include <Champs_compris.h>
#include <Champs_compris_interface.h>
#include <Champ_Don.h>
#include <TRUST_Ref.h>
#include <Champ_Inc.h>



class Fluide_base;
class Motcle;
class Equation_base;
class Probleme_base;
class Discretisation_base;
class Champ_base;

class Modele_Fonc_Realisable_base : public Champs_compris_interface, public Objet_U
{
  Declare_base(Modele_Fonc_Realisable_base);
public:

  static void typer_lire_Modele_Fonc_Realisable(OWN_PTR(Modele_Fonc_Realisable_base)&, const Equation_base&, Entree& is );

  inline const Equation_base& equation() const;
  inline  Equation_base& equation();
  inline const Equation_base& seconde_equation() const;
  inline  Equation_base& seconde_equation();
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
  virtual void associer(const Domaine_dis_base& , const Domaine_Cl_dis_base& )= 0;
  int sauvegarder(Sortie& ) const override;
  int reprendre(Entree& ) override;
  virtual int Calcul_is_Reynolds_stress_isotrope() const;
  virtual int Calcul_is_Cmu_constant() const;

  virtual void Calcul_S(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse) =0;
  virtual void Calcul_C1(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) =0;
  virtual void Calcul_Cmu_et_S(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN)  =0;
  virtual void Contributions_Sources(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) =0;
  virtual void Contributions_Sources_Paroi(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                           const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) =0;

  virtual void Calcul_C1_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) =0;
  virtual void Calcul_Cmu_et_S_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN)  =0;
  virtual void Contributions_Sources_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) =0;
  virtual void Contributions_Sources_Paroi_BiK(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN,
                                               const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt) =0;

  Entree& lire(const Motcle&, Entree&);

  //Methodes de l interface des champs postraitables
  /////////////////////////////////////////////////////
  void creer_champ(const Motcle& motlu) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  /////////////////////////////////////////////////////

  void associer_eqn(const Equation_base& );

  virtual void associer_eqn_2(const Equation_base& );

  bool has_seconde_equation() const { return ma_seconde_equation.non_nul(); }

public:
  REF(Equation_base) ma_seconde_equation;

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

inline const Equation_base& Modele_Fonc_Realisable_base::seconde_equation() const
{
  if (ma_seconde_equation.non_nul()==0)
    {
      Cerr << "\nError in Modele_Fonc_Realisable_base::seconde_equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return ma_seconde_equation.valeur();
}

inline Equation_base& Modele_Fonc_Realisable_base::seconde_equation()
{
  if (ma_seconde_equation.non_nul()==0)
    {
      Cerr << "\nError in Modele_Fonc_Realisable_base::seconde_equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return ma_seconde_equation.valeur();
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



