//////////////////////////////////////////////////////////////////////////////
//
// File:        K_Eps_residu_P1NC.h
//
//////////////////////////////////////////////////////////////////////////////

#ifndef K_Eps_residu_P1NC_included
#define K_Eps_residu_P1NC_included


#include <Champ_Fonc_P1NC.h>
#include <Ref_Champ_P1NC.h>
#include <DoubleVect.h>
#include <Zone_dis.h>
#include <Nom.h>

class K_Eps_residu_P1NC : public Champ_Fonc_P1NC
{

  Declare_instanciable( K_Eps_residu_P1NC );

public:

  inline const Champ_P1NC& mon_champ( ) const;
  inline void mettre_a_jour( double );
  void associer_champ( const Champ_P1NC& );
  void associer_zone( const Zone_dis& zone_dis );
  void me_calculer( double );
  /* void init_zone( void ); */
protected:

  REF( Champ_P1NC ) K_Eps_;
  DoubleVect volumes_ ;
  Zone_dis zone_dis_;
};

inline void K_Eps_residu_P1NC::mettre_a_jour( double tps )
{
  me_calculer( tps );
  changer_temps( tps );
  Champ_Fonc_base::mettre_a_jour( tps );
}




#endif
