//////////////////////////////////////////////////////////////////////////////
//
// File:        Trait_part_NS_canal_VDF.h
// Directory:   $TRIO_U_ROOT/VDF
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_canal_VDF_inclus
#define Traitement_particulier_NS_canal_VDF_inclus

#include <Trait_part_NS_canal.h>

class Navier_Stokes_Turbulent;

//////////////////////////////////////////////////////////////////////////////
// .NOM  Traitement_particulier_NS_canal_VDF
// .ENTETE  Trio_U ThHyd/Turbulence
// .LIBRAIRIE  libvdf
// .FILE  Trait_part_NS_canal_VDF.h
// .FILE  Trait_part_NS_canal_VDF.cpp
// 
// .DESCRIPTION 
//     classe Traitement_particulier_NS_canal_VDF
//     Cette classe permet de faire les traitements particuliers
//     pour le calcul d'un canal plan :
//         * conservation du debit
//         * calculs de moyennes
//     
// .SECTION voir aussi 
//      Navier_Stokes_Turbulent, Traitement_particulier_base,
//      Traitement_particulier_VDF
// .CONTRAINTES 
// .INVARIANTS 
// .HTML 
// .EPS 
//////////////////////////////////////////////////////////////////////////////
class Traitement_particulier_NS_canal_VDF : public Traitement_particulier_NS_canal
{ 
  Declare_instanciable(Traitement_particulier_NS_canal_VDF);
      
    public :

  Entree& lire(const Motcle& , Entree& );
  Entree& lire(Entree& );

  protected :
  void remplir_Tab_recap(DoubleTab& ); // Ajour F.A 15/02/11
  void remplir_Y(DoubleVect&, DoubleVect&, entier& ) const;
  void calculer_moyenne_spatiale_vitesse_rho_mu(DoubleTab&) const;
  void calculer_moyenne_spatiale_nut(DoubleTab&) const;
  void calculer_moyenne_spatiale_Temp(DoubleTab&) const;
  void calculer_Temp(DoubleTab& val_moy,const DoubleTab& diffusivite_turb) const;
};
#endif
