//////////////////////////////////////////////////////////////////////////////
//
// File:        ISC_Qch_Exponential_Linearized_P1NC.cpp
//
//////////////////////////////////////////////////////////////////////////////

#include <K_Eps_residu_P1NC.h>
#include <Champ_P1NC.h>
#include <DoubleTab.h>
#include <Zone_VEF.h>
#include <Zone_VF.h>
#include <Champ_Generique_base.h>
#include <Champ.h>
#include <Schema_Temps_base.h>

Implemente_instanciable( K_Eps_residu_P1NC, "K_Eps_residu_P1NC", Champ_Fonc_P1NC );

Sortie& K_Eps_residu_P1NC::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& K_Eps_residu_P1NC::readOn(Entree& s)
{
  return s ;
}

// void K_Eps_residu_P1NC::init_zone( void )
// {
//   volumes_.ref( ref_cast( Zone_VF,zone_dis_.valeur( ) ).volumes_entrelaces( ) );
// }

void K_Eps_residu_P1NC::associer_champ( const Champ_P1NC& K_Eps )
{
  K_Eps_= K_Eps ;
}

void K_Eps_residu_P1NC::associer_zone( const Zone_dis& zone_dis )
{
  zone_dis_ = zone_dis ;
}

void K_Eps_residu_P1NC::me_calculer( double tps )
{
  int nb_faces = zone_vef( ).nb_faces( );
  const DoubleTab K_Eps = K_Eps_.valeur( ).passe( );
  DoubleTab& champ = le_champ( ).valeurs( );

  if( nb_faces != K_Eps.dimension( 0 ) )
    {
      Cerr << "Error in K_Eps_residu_P1NC::me_calculer "<<finl;
      Cerr << "There are "<<nb_faces<<" faces and "<<K_Eps.dimension( 0 )<<" values for K_Eps field "<<finl;
      Process::abort( );
    }

  if( tps > 0.0 )
	  champ = K_Eps_.valeur( ).equation( ).get_residuals();

  else
    {
      champ = -10000.0 ;
      Cerr << "[Information] K_Eps_residu_P1NC::me_calculer : le residu est mis a -10000.0 au temps initial"<<finl;
    }

}
