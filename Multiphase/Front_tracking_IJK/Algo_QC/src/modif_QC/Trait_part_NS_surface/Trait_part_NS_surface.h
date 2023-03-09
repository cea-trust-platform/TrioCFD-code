/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Trait_part_NS_surface.h
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_surface
//
/////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_surface_inclus
#define Traitement_particulier_NS_surface_inclus

#include <Traitement_particulier_NS_base.h>
#include <IntTab.h>
class Navier_Stokes_Turbulent;

/*! @brief classe Traitement_particulier_NS_surface Cette classe permet de faire les traitements particuliers
 *
 *      pour le calcul d'un canal surface :
 *          * conservation du debit
 *          * calculs de moyennes
 *
 *
 * @sa Navier_Stokes_Turbulent, Traitement_particulier_base,, Traitement_particulier_VDF
 */
class Traitement_particulier_NS_surface : public Traitement_particulier_NS_base
{
  Declare_base(Traitement_particulier_NS_surface);

public :

  Entree& lire(Entree& );
  void preparer_calcul_particulier() ;
  virtual void reprendre_stat();
  virtual void sauver_stat() const;
  void post_traitement_particulier();
  inline void en_cours_de_resolution(int , DoubleTab&, DoubleTab& ,double) ;
  inline int a_pour_Champ_Fonc(const Motcle& mot, REF(Champ_base)& ch_ref) const  ;
  inline int comprend_champ(const Motcle& mot) const  ;

protected :

  virtual void remplir_XYZ(DoubleVect& ,DoubleVect& ,DoubleVect&  ,int&, int& ,int&,IntTab&) const = 0;

  void remplir_reordonne_Z_tot(const DoubleVect&,  DoubleVect&) const;
  void remplir_reordonne_Y_tot(const DoubleVect&,  DoubleVect&) const;
  void remplir_reordonne_X_tot(const DoubleVect&,  DoubleVect&) const;
  virtual void recuperation_grandeurs(DoubleTab&) const = 0;

  void calcul_reynolds_tau();
  void ecriture_fichiers_moy_vitesse_rho_mu(const DoubleTab&, const Nom&, int&, int&) const;
  void ecriture_fichiers_moy_vitesse_rho_mu_gnuplot(const DoubleTab&, const Nom&, int&, int&) const;
  void ecriture_fichiers_moy_nut(const DoubleTab&, const Nom&, const double&, const int&) const;
  void ecriture_fichiers_moy_Temp(const DoubleTab&, const Nom&, const double&, const int&) const;

  void ecriture_fichiers_moy_vitesse_rho_mu_old(const DoubleTab&, const Nom&, const double&, const int&) const;
  void ecriture_fichiers_moy_Temp_old(const DoubleTab&, const Nom&, const double&, const int&) const;

  int     Ny,Nx,Np,Nz,Nval,Nphase;
  DoubleTab  val_moy_tot,val_moy_temp,val_moy_phase;
  DoubleVect Z,X,Y,Z_tot,X_tot,Y_tot,Yp;
  DoubleTab  compt,compt_tot;
  IntTab     Tab_post,Tab_recap,yt,PosX,PosY; // Tab_post contient le tableau range par surface des elements a post traite
  IntVect    Nb_ech_phase;     // nombre d'echantillons par phase
  double     w,freq;
  int	     ind_phase;
  int choix_fichier;
  int q; // Variable de decalage d'indice

  int    oui_profil_nu_t,oui_profil_Temp,oui_repr,oui_pulse;
  double temps_deb,temps_fin,debut_phase;
  double dt_impr_moy_spat,dt_impr_moy_temp;
  Nom    fich_repr;

  REF(Champ_base) Temp;

  DoubleVect uv_moy_temp,wv_moy_temp,uw_moy_temp,u2_moy_temp,v2_moy_temp,w2_moy_temp,nu_t_temp;
};
#endif


inline int Traitement_particulier_NS_surface::a_pour_Champ_Fonc(const Motcle& mot, REF(Champ_base)& ch_ref) const
{
  return 0 ;
}

inline int Traitement_particulier_NS_surface::comprend_champ(const Motcle& mot) const
{
  return 0 ;
}

inline void Traitement_particulier_NS_surface::en_cours_de_resolution(int i, DoubleTab& a, DoubleTab& b,double c)
{
  ;
}
