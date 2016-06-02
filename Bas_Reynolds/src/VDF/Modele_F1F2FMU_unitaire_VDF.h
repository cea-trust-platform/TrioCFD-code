//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_F1F2FMU_unitaire_VDF.h
// Directory:   $TRUST_ROOT/VDF/Turbulence
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_F1F2FMU_unitaire_VDF_included
#define Modele_F1F2FMU_unitaire_VDF_included

#include <Modele_Jones_Launder_VDF.h>

class Zone_dis;
class Zone_Cl_dis;
class DoubleVect;
class DoubleTab;
class Zone_Cl_VDF;
class Champ_Face;

class Modele_F1F2FMU_unitaire_VDF : public Modele_Jones_Launder_VDF
{

  Declare_instanciable(Modele_F1F2FMU_unitaire_VDF);

public :

  DoubleTab& Calcul_Fmu (DoubleTab&,const Zone_dis&,const Zone_Cl_dis&,const DoubleTab&,const Champ_Don&) const;
  DoubleTab& Calcul_F1(DoubleTab&, const Zone_dis&, const Zone_Cl_dis&, const DoubleTab&,const DoubleTab&,const Champ_base& ) const ;
  DoubleTab& Calcul_F2(DoubleTab&, DoubleTab&,const Zone_dis&,const DoubleTab&,const Champ_base&) const ;

  void associer(const Zone_dis& , const Zone_Cl_dis& );
  Entree& lire(const Motcle&, Entree&);

protected:

  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;
};

#endif

