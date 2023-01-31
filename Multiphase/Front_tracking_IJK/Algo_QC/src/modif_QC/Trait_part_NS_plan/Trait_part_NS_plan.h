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
// File      : Trait_part_NS_plan.h
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_plan
//
/////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_plan_inclus
#define Traitement_particulier_NS_plan_inclus

#include <Traitement_particulier_NS_base.h>
#include <IntTab.h>
#include <IJK_Splitting.h>
#include <VDF_to_IJK.h>

class Navier_Stokes_Turbulent;

/*! @brief classe Traitement_particulier_NS_plan Cette classe permet de faire les traitements particuliers
 *
 *      pour le calcul d'un canal plan :
 *          * conservation du debit
 *          * calculs de moyennes
 *
 *
 * @sa Navier_Stokes_Turbulent, Traitement_particulier_base,, Traitement_particulier_VDF
 */
class Traitement_particulier_NS_plan : public Traitement_particulier_NS_base
{
  Declare_base(Traitement_particulier_NS_plan);

public :

  Entree& lire(Entree& );
  void preparer_calcul_particulier() ;
  void post_traitement_particulier();
  void sauver_stat() const {};
  void reprendre_stat() {};
  inline void en_cours_de_resolution(int , DoubleTab&, DoubleTab& ,double) ;
  inline int a_pour_Champ_Fonc(const Motcle& mot, REF(Champ_base)& ch_ref) const  ;
  inline int comprend_champ(const Motcle& mot) const  ;

protected :
  virtual void calculer_valeur_spatiale_vitesse_rho_mu(DoubleTab&,DoubleTab&,double) = 0;
  int calculer_critere_post() const;
  void calculer_moyennes(DoubleTab& Moyennes, DoubleTab& rrhouf, DoubleTab& Racinederho) const;
  void ecriture_val_post(const DoubleTab& val_post, double current_time) const;

#ifdef OLD_CODE
  int     Ny,Nx,Np,Nz,Nval;
  int oui_repr, oui_calcul ; // Indique si il y a une reprise et si on imprime.

  double Pth_passe; // stock la valeur de Pth au pas de temps precedent.
  double temps_deb,temps_fin, dt_impr_moy , temps_de_moyenne; // variables regissant les moyennes.

  DoubleTab  Moyennes_temporelle,Moyennes_finale,sauv_moyennes,Rrhouf,racinederho;

  DoubleVect Z,X,Y,Z_tot,X_tot,Y_tot,Yp;

  IntTab     Tab_post,Tab_recap,yt,Ref_Y,Envoyer_val,PtoTX,PtoTZ; // tableau qui va chercher a savoir si il faut envoyer les donnees; // Tab_post contient le tableau range par plan des elements a post traite
  IntTab    for_each_plane_recv_flag_; // tableau de stokage des informations des plans qui post-traites
#else
  // a revoir pour ijk:
  IntVect Ref_Y; // pour chaque element de la domaine VDF sur ce processeur, indice global de l'element dans la direction k
  DoubleTab  Rrhouf, racinederho;
  IntTab Tab_recap;
#endif

  IJK_Splitting ijk_splitting_;
  VDF_to_IJK vdf_to_ijk_; // converter to ijk field

  // Index of each plane to process (in y direction)
  ArrOfInt indices_planes_;
  // This flag tells if we shall
  // Physical time at which we start integration and postprocessing
  double start_time_;
  // Time interval between writing spatial average values
  double dt_write_spatial_avg_;

  int oui_calcul_;
  int oui_repr_;
  Nom    fich_repr_; // Nom du fichier pour la reprise des stats (Moyennes_temporelles_)
  int nb_moyennes_temporelles() const
  {
    return 43;
  }
  int nb_valeurs_spatiales_x() const
  {
    return 29;
  }
  // Pour chaque plan de maillage: dimensin(0) = nb_elem_z, dimension(1) = nb_moyennes_temporelles()
  DoubleTab Moyennes_temporelles_; // alloue uniquement sur le maitre
  DoubleTab  Moyennes_finale_;
  DoubleTab sauv_moyennes_; // alloue uniquement sur le maitre
  double temps_de_moyenne_;

  REF(Champ_base) ref_champ_temperature_;
  double Pth_passe_; // stock la valeur de Pth au pas de temps precedent.

};
#endif


inline int Traitement_particulier_NS_plan::a_pour_Champ_Fonc(const Motcle& mot, REF(Champ_base)& ch_ref) const
{
  return 0 ;
}

inline int Traitement_particulier_NS_plan::comprend_champ(const Motcle& mot) const
{
  return 0 ;
}

inline void Traitement_particulier_NS_plan::en_cours_de_resolution(int i, DoubleTab& a, DoubleTab& b,double c)
{
  ;
}
