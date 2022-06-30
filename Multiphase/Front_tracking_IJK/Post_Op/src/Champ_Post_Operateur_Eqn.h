//////////////////////////////////////////////////////////////////////////////
//
// File:        Champ_Post_Operateur_Eqn.h
// Directory:   $TRIO_U_ROOT/Kernel/Champs
// Version:     /main/1
//
//////////////////////////////////////////////////////////////////////////////

// .NOM Champ_Post_Operateur_Eqn 
// .ENTETE Trio_U Kernel/Champs
// .LIBRAIRIE 
// .FILE Champ_Post_Operateur_Eqn.h
// .FILE Champ_Post_Operateur_Eqn.cpp

#ifndef Champ_Post_Operateur_Eqn_inclus
#define Champ_Post_Operateur_Eqn_inclus

#include <Champ_Generique_Operateur_base.h>
#include <Operateur.h>
#include <Ref_Eqn_base.h>
//
// .DESCRIPTION class Champ_Post_Operateur_Eqn
// Champ destine à post-traiter le gradient d un champ generique
// La classe porte un operateur statistique "gradient"

//// Syntaxe à respecter pour jdd
//
// "nom_champ" Champ_Post_Statistiques_Eqn { 
//		source Champ_Post...{ ...source Champ_Post_ref_Champ { Pb_champ "nom_pb" "nom_champ_discret" } }
//	       }
// "nom_champ" fixe par utilisateur sera le nom du champ generique
//Ce type de champ implique que le champ source possede des conditions limites
//Son application est restreinte à certains champs discrets (pression VDF et VEF ou temperature en VEF)

class Champ_Post_Operateur_Eqn : public Champ_Generique_Operateur_base
{
   
  Declare_instanciable(Champ_Post_Operateur_Eqn);
      
    public:
    
  virtual const Noms get_property(const Motcle & query) const;
  virtual Entity  get_localisation(const entier index = -1) const;
  virtual const Champ_base&  get_champ(Champ& espace_stockage) const;
   
  const Operateur_base& Operateur() const;
  Operateur_base& Operateur();
  void completer(const Postraitement_base& post);
  void nommer_source();
  //virtual Entree &   lire(const Motcle & motcle, Entree & is);
  void set_param(Param& param);
  
  
 protected:
  int numero_source_,numero_op_;
  REF(Equation_base) ref_eq_;    
  int sans_solveur_masse_; 
  Entity localisation_inco_;
};

#endif
