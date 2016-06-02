//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_F1F2FMU_unitaire_VDF.cpp
// Directory:   $TRUST_ROOT/VDF/Turbulence
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_F1F2FMU_unitaire_VDF.h>
#include <Zone_VDF.h>
#include <Champ_Uniforme.h>

Implemente_instanciable(Modele_F1F2FMU_unitaire_VDF,"Modele_F1F2FMU_unitaire_VDF",Modele_Jones_Launder_VDF);



///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////
// printOn et readOn

Sortie& Modele_F1F2FMU_unitaire_VDF::printOn(Sortie& s ) const
{
  Modele_Jones_Launder_VDF::printOn(s);
  return s;
}

Entree& Modele_F1F2FMU_unitaire_VDF::readOn(Entree& is )
{
  Modele_Jones_Launder_VDF::readOn(is);
  return is;
}

Entree& Modele_F1F2FMU_unitaire_VDF::lire(const Motcle& , Entree& is)
{
  return is;
}

void  Modele_F1F2FMU_unitaire_VDF::associer(const Zone_dis& zone_dis,
                                            const Zone_Cl_dis& zone_Cl_dis)
{
}

DoubleTab&  Modele_F1F2FMU_unitaire_VDF::Calcul_Fmu( DoubleTab& Fmu,const Zone_dis& zone_dis,const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  Fmu=1;
  Cerr<<Fmu.mp_min_vect()<<" Fmu "<<Fmu.mp_max_vect()<<finl;
  return Fmu;
}
DoubleTab& Modele_F1F2FMU_unitaire_VDF::Calcul_F1( DoubleTab& F1, const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{
  F1= 1.;
  return F1;
}

DoubleTab& Modele_F1F2FMU_unitaire_VDF::Calcul_F2( DoubleTab& F2, DoubleTab& D, const Zone_dis& zone_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  F2=1;
  return F2;
}
