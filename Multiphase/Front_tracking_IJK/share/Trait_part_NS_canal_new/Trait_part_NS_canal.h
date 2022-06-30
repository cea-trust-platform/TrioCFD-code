//////////////////////////////////////////////////////////////////////////////
//
// File:        Trait_part_NS_canal.h
// Directory:   $TRIO_U_ROOT/ThHyd/Turbulence
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_canal_inclus
#define Traitement_particulier_NS_canal_inclus

#include <Trait_part_NS_base.h>

class Navier_Stokes_Turbulent;

//////////////////////////////////////////////////////////////////////////////
// .NOM  Traitement_particulier_NS_canal
// .ENTETE  Trio_U ThHyd/Turbulence
// .LIBRAIRIE  libhydturb
// .FILE  Trait_part_NS_canal.h
// .FILE  Trait_part_NS_canal.cpp
// 
// .DESCRIPTION 
//     classe Traitement_particulier_NS_canal
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
class Traitement_particulier_NS_canal : public Traitement_particulier_NS_base
{ 
  Declare_base(Traitement_particulier_NS_canal);
      
    public :

  Entree& lire(Entree& );
  void preparer_calcul_particulier() ;
  virtual void reprendre_stat();
  virtual void sauver_stat() const;
  void post_traitement_particulier();
  inline void en_cours_de_resolution(int , DoubleTab&, DoubleTab& ,double) ;
  inline entier a_pour_Champ_Fonc(const Motcle& mot, REF(Champ_base)& ch_ref) const  ; 
  inline entier comprend_champ(const Motcle& mot) const  ; 

  protected :
  
  //void init_calcul_stats(void);
  void reprendre_stat_canal(DoubleTab&, const Nom& );
  void sauver_stat_canal(const DoubleTab&, const Nom& ) const;
  virtual void remplir_Y(DoubleVect&,  DoubleVect&, entier& ) const = 0;
  
  virtual void remplir_Tab_recap(DoubleTab& ); // Ajour F.A 15/02/11
  
  void remplir_reordonne_Y_tot(const DoubleVect&,  DoubleVect&) const;
  virtual void calculer_moyenne_spatiale_vitesse_rho_mu(DoubleTab&) const = 0;
  virtual void calculer_moyenne_spatiale_nut(DoubleTab&) const = 0;
  virtual void calculer_moyenne_spatiale_Temp(DoubleTab&) const = 0;
  void calcul_reynolds_tau();
  void ecriture_fichiers_moy_vitesse_rho_mu(const DoubleTab&, const Nom&, const double&, const int&) const;
  void ecriture_fichiers_moy_nut(const DoubleTab&, const Nom&, const double&, const int&) const;
  void ecriture_fichiers_moy_Temp(const DoubleTab&, const Nom&, const double&, const int&) const;

  void ecriture_fichiers_moy_vitesse_rho_mu_old(const DoubleTab&, const Nom&, const double&, const int&) const;
  void ecriture_fichiers_moy_nut_old(const DoubleTab&, const Nom&, const double&, const int&) const; // modified AT 5/06/09
  void ecriture_fichiers_moy_Temp_old(const DoubleTab&, const Nom&, const double&, const int&) const;
  
  entier     Ny,Nval,Nphase;
  DoubleTab  val_moy_tot,val_moy_temp,val_moy_phase,Tab_recap; // ajout de Tab_recap F.A 15/02/11
  DoubleVect Y,Y_tot;
  DoubleVect compt,compt_tot;
  IntVect    Nb_ech_phase;     // nombre d'echantillons par phase
  double     w,freq;
  int	     ind_phase;

  int    oui_profil_nu_t,oui_profil_Temp,oui_repr,oui_pulse;
  double temps_deb,temps_fin,debut_phase;
  double dt_impr_moy_spat,dt_impr_moy_temp;
  Nom    fich_repr;
  
  REF(Champ_base) Temp;

  DoubleVect uv_moy_temp,wv_moy_temp,uw_moy_temp,u2_moy_temp,v2_moy_temp,w2_moy_temp,nu_t_temp;
};
#endif


inline entier Traitement_particulier_NS_canal::a_pour_Champ_Fonc(const Motcle& mot, REF(Champ_base)& ch_ref) const
{
  return 0 ;
}

inline entier Traitement_particulier_NS_canal::comprend_champ(const Motcle& mot) const
{
  return 0 ;
}

inline void Traitement_particulier_NS_canal::en_cours_de_resolution(int i, DoubleTab& a, DoubleTab& b,double c) 
{
  ;
}
