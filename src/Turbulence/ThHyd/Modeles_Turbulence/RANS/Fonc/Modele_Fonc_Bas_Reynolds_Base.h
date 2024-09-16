/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Modele_Fonc_Bas_Reynolds_Base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Fonc_Bas_Reynolds_Base_included
#define Modele_Fonc_Bas_Reynolds_Base_included


#include <Champs_compris_interface.h>
#include <Champs_compris.h>


#include <Champ_Fonc.h>
#include <Champ_Don.h>
#include <Champ_Inc.h>
#include <TRUST_Ref.h>

class Fluide_base;
class Motcle;
class Equation_base;
class Probleme_base;
class Discretisation_base;
class Champ_base;

class Modele_Fonc_Bas_Reynolds_Base : public Champs_compris_interface, public Objet_U
{

  Declare_base(Modele_Fonc_Bas_Reynolds_Base);

public:

  static void typer_lire_Modele_Fonc_Bas_Reynolds(OWN_PTR(Modele_Fonc_Bas_Reynolds_Base)&, const Equation_base&, Entree& is );

  inline const Equation_base& equation() const;
  inline  Equation_base& equation();
  inline const Equation_base& seconde_equation() const;
  inline  Equation_base& seconde_equation();
  virtual int preparer_calcul();
  virtual void mettre_a_jour(double ) =0;
  virtual void discretiser();
  virtual void completer();
  virtual void associer_pb(const Probleme_base& ) ;
  virtual void associer_eqn(const Equation_base& );
  virtual void associer_eqn_2(const Equation_base& );
  virtual void associer(const Domaine_dis_base& , const Domaine_Cl_dis_base& )= 0;
  int sauvegarder(Sortie& ) const override;
  int reprendre(Entree& ) override;
  virtual DoubleTab& Calcul_D(DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const=0;
  virtual int Calcul_is_Reynolds_stress_isotrope() const;
  virtual int Calcul_is_Cmu_constant() const;
  virtual DoubleTab& Calcul_E(DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const =0 ;

//  virtual DoubleTab& Calcul_F1(DoubleTab&, const Domaine_dis_base& ) const =0 ;
  virtual DoubleTab& Calcul_F1( DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&, const DoubleTab&,const DoubleTab&,const Champ_base&) const=0;
  virtual DoubleTab& Calcul_F2(DoubleTab&, DoubleTab&,const Domaine_dis_base&,const DoubleTab&,const Champ_base&) const =0 ;
  virtual DoubleTab& Calcul_Fmu ( DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const Champ_Don& )const =0 ;
  virtual DoubleTab& Calcul_Cmu(DoubleTab&,const Domaine_dis_base&, const Domaine_Cl_dis_base&, const DoubleTab&, const DoubleTab&, const double) const;
  virtual DoubleTab& Calcul_Cmu_Paroi(DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&,
                                      const DoubleTab& , const DoubleTab& ,
                                      const DoubleTab& ,const int,
                                      const DoubleTab&, const DoubleTab&,
                                      const double) const;

  virtual bool calcul_tenseur_Re(const DoubleTab&, const DoubleTab&, DoubleTab&) const;

  virtual DoubleTab& Calcul_D_BiK(DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const DoubleTab&, const Champ_Don&) const=0;
  virtual DoubleTab& Calcul_E_BiK(DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const Champ_Don&, const DoubleTab& ) const =0 ;

  virtual DoubleTab& Calcul_F1_BiK( DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&, const DoubleTab&,const DoubleTab&,const DoubleTab&,const Champ_base&) const=0;
  virtual DoubleTab& Calcul_F2_BiK(DoubleTab&, DoubleTab&,const Domaine_dis_base&,const DoubleTab&,const DoubleTab&,const Champ_base&) const =0 ;
  virtual DoubleTab& Calcul_Fmu_BiK ( DoubleTab&,const Domaine_dis_base&,const Domaine_Cl_dis_base&,const DoubleTab&,const DoubleTab&,const Champ_Don& )const =0 ;
  virtual DoubleTab& Calcul_Cmu_BiK(DoubleTab&,const Domaine_dis_base&, const Domaine_Cl_dis_base&, const DoubleTab&, const DoubleTab&, const DoubleTab&, const double) const;
  virtual DoubleTab& Calcul_Cmu_Paroi_BiK(DoubleTab&, const Domaine_dis_base&, const Domaine_Cl_dis_base&,
                                          const DoubleTab& , const DoubleTab& ,
                                          const DoubleTab& ,const int,
                                          const DoubleTab&, const DoubleTab&, const DoubleTab&,
                                          const double) const;
  virtual bool calcul_tenseur_Re_BiK(const DoubleTab&, const DoubleTab&, DoubleTab&) const;

  //Methodes de l interface des champs postraitables
  /////////////////////////////////////////////////////
  void creer_champ(const Motcle& motlu) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  /////////////////////////////////////////////////////
  virtual void lire_distance_paroi( );

  bool has_seconde_equation() const { return ma_seconde_equation.non_nul(); }

public:
  REF(Equation_base) ma_seconde_equation;

protected :
  REF(Equation_base) mon_equation;
  REF(Fluide_base) le_fluide;
  REF(Champ_Inc) la_vitesse_transportante;
  REF(Equation_base) eq_hydraulique;
  REF(Champ_Don) visco_;

  Nom nom_fic;
  Champ_Fonc BR_wall_length_;
  int is_Cmu_constant_;
  int is_Reynolds_stress_isotrope_;
  Champ_Don D_,E_,F1_,F2_;
private :

  Champs_compris champs_compris_;

};

//
// fonction inline
//



inline const Equation_base& Modele_Fonc_Bas_Reynolds_Base::equation() const
{
  if (mon_equation.non_nul()==0)
    {
      Cerr << "\nError in Modele_Fonc_Bas_Reynolds_Base::equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return mon_equation.valeur();
}

inline Equation_base& Modele_Fonc_Bas_Reynolds_Base::equation()
{
  if (mon_equation.non_nul()==0)
    {
      Cerr << "\nError in Modele_Fonc_Bas_Reynolds_Base::equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return mon_equation.valeur();
}

inline const Equation_base& Modele_Fonc_Bas_Reynolds_Base::seconde_equation() const
{
  if (ma_seconde_equation.non_nul()==0)
    {
      Cerr << "\nError in Modele_Fonc_Bas_Reynolds_Base::seconde_equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return ma_seconde_equation.valeur();
}

inline Equation_base& Modele_Fonc_Bas_Reynolds_Base::seconde_equation()
{
  if (ma_seconde_equation.non_nul()==0)
    {
      Cerr << "\nError in Modele_Fonc_Bas_Reynolds_Base::seconde_equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return ma_seconde_equation.valeur();
}

#endif



