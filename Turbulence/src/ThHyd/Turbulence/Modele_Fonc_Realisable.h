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
// File:        Modele_Fonc_Realisable.h
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_Fonc_Realisable_included
#define Modele_Fonc_Realisable_included

#include <MorEqn.h>
#include <Modele_Fonc_Realisable_base.h>

Declare_deriv(Modele_Fonc_Realisable_base);



class Modele_Fonc_Realisable : public MorEqn, public DERIV(Modele_Fonc_Realisable_base)
{
  Declare_instanciable(Modele_Fonc_Realisable);

public:

  inline const DoubleTab& get_S() const;
  inline  DoubleTab& get_S();
  inline const DoubleTab& get_Cmu( void ) const;
  inline DoubleTab& get_Cmu( void );
  inline const DoubleTab& get_C1( void ) const;
  inline DoubleTab& get_C1( void );

  inline int preparer_calcul();
  inline void mettre_a_jour(double );
  inline void discretiser();
  inline void completer();
  inline int sauvegarder(Sortie& os) const;
  inline int reprendre(Entree& is);
  inline int Calcul_is_Reynolds_stress_isotrope() const;
  inline int Calcul_is_Cmu_constant() const;

  inline void Calcul_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse) ;
  inline void Calcul_C1(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) ;
  inline void Calcul_Cmu_et_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN)  ;
  inline void Contributions_Sources(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN) ;
  inline void Contributions_Sources_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                          const DoubleTab& visco, const DoubleTab& visco_turb,const DoubleTab& loi_paroi,const int idt) ;

};

///
//    implementation des fonctions inline
///

inline const DoubleTab& Modele_Fonc_Realisable::get_S() const
{
  return valeur().get_S();
}

inline  DoubleTab& Modele_Fonc_Realisable::get_S()
{
  return valeur().get_S();
}

inline const DoubleTab& Modele_Fonc_Realisable::get_Cmu( void ) const
{
  return valeur().get_Cmu();
}

inline DoubleTab& Modele_Fonc_Realisable::get_Cmu( void )
{
  return valeur().get_Cmu();
}

inline const DoubleTab& Modele_Fonc_Realisable::get_C1( void ) const
{
  return valeur().get_C1();
}

inline DoubleTab& Modele_Fonc_Realisable::get_C1( void )
{
  return valeur().get_C1();
}

inline int Modele_Fonc_Realisable::preparer_calcul()
{
  return valeur().preparer_calcul();
}

inline void Modele_Fonc_Realisable::mettre_a_jour(double temps)
{
  valeur().mettre_a_jour(temps);
}

inline void Modele_Fonc_Realisable::discretiser()
{
  valeur().discretiser();
}

inline void Modele_Fonc_Realisable::completer()
{
  valeur().completer();
}

inline int Modele_Fonc_Realisable::sauvegarder(Sortie& os) const
{
  return valeur().sauvegarder(os);
}

inline int Modele_Fonc_Realisable::reprendre(Entree& is)
{
  return valeur().reprendre(is);
}

inline int Modele_Fonc_Realisable::Calcul_is_Cmu_constant() const
{
  return valeur().Calcul_is_Cmu_constant();
}

inline int Modele_Fonc_Realisable::Calcul_is_Reynolds_stress_isotrope() const
{
  return valeur().Calcul_is_Reynolds_stress_isotrope();
}

inline void Modele_Fonc_Realisable::Calcul_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse)
{
  valeur().Calcul_S( zone_dis, zone_Cl_dis, vitesse );
}

inline void Modele_Fonc_Realisable::Calcul_C1(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN)
{
  valeur().Calcul_C1( zone_dis, zone_Cl_dis, vitesse, K_Eps, EPS_MIN );
}

inline void Modele_Fonc_Realisable::Calcul_Cmu_et_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN)
{
  valeur().Calcul_Cmu_et_S( zone_dis, zone_Cl_dis, vitesse, K_Eps, EPS_MIN );
}

inline void Modele_Fonc_Realisable::Contributions_Sources(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN)
{
  valeur().Contributions_Sources( zone_dis, zone_Cl_dis, vitesse, K_Eps, EPS_MIN );
}

inline void Modele_Fonc_Realisable::Contributions_Sources_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,const DoubleTab& visco, const DoubleTab& visco_turb,const DoubleTab& loi_paroi,const int idt)
{
  valeur().Contributions_Sources_Paroi( zone_dis, zone_Cl_dis, vitesse, K_Eps, EPS_MIN,visco,visco_turb,loi_paroi,idt);
}

#endif
