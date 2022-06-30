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
// File      : Couplage_Tubes_IBC.h
// Directory : $IJK_ROOT/src/IBC
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Couplage_Tubes_IBC_H
#define Couplage_Tubes_IBC_H
#include <Objet_U.h>
#include <TRUSTTab.h>
#include <IJK_Field.h>
#include <Tube_base.h>

#define MOUVEMENT_LIBRE 0
#define MOUVEMENT_IMPOSE 1
#define ETALEMENT 0
#define PAS_ETALEMENT 1
#define MIROIR 0
#define PAS_MIROIR 1
#define SYMETRIE_PLANE 2
#define IBC0 0
#define IBC_DIFFUSE 1
#define IBC_LOCALISEE 2
#define IBC_LOCALISEE_QDM 3
#define CYLINDRE 0
#define CUBE 1

class EFichier;
class SFichier;

class Couplage_Tubes_IBC : public Objet_U
{
  Declare_instanciable(Couplage_Tubes_IBC);
public:
  void initialize(const IJK_Splitting&);

  void force_ibc(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                 const IJK_Field_double& rho_field,
                 double timestep, const IJK_Field_double& pressure, double current_time);
  void calcul_force_post_projection(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                    const IJK_Field_double& rho_field,
                                    double timestep, double current_time);
  void champ_miroir(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                    const IJK_Field_double& rho_field,
                    double timestep, const IJK_Field_double& pressure);
  void forcage_anticipe(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                        const IJK_Field_double& rho_field,
                        double timestep, const IJK_Field_double& pressure);
  void force_ibc_velocity_anticipe_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                        const IJK_Field_double& rho_field,
                                        DoubleTab& masse_fluide_cylindres,
                                        DoubleTab& volume_cylindres,
                                        DoubleTab& integrale_force,
                                        const Faisceau_Tubes& ,
                                        const double timestep) const;

  void force_ibc_velocity(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                          const IJK_Field_double& rho_field,
                          DoubleTab& masse_fluide_cylindres,
                          DoubleTab& volume_cylindres,
                          DoubleTab& integrale_force,
                          const Faisceau_Tubes& ) const;
  void force_ibc_velocity_frac_vol(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                   const IJK_Field_double& rho_field,
                                   DoubleTab& masse_fluide_cylindres,
                                   DoubleTab& volume_cylindres,
                                   DoubleTab& integrale_force,
                                   const Faisceau_Tubes& ) const;
  void ibc0_velocity_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                          const IJK_Field_double& rho_field,
                          DoubleTab& masse_fluide_cylindres,
                          DoubleTab& volume_cylindres,
                          DoubleTab& integrale_force,
                          const Faisceau_Tubes& ) const;
  void ibc_diffuse_velocity_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                 const IJK_Field_double& rho_field,
                                 DoubleTab& masse_fluide_cylindres,
                                 DoubleTab& volume_cylindres,
                                 DoubleTab& integrale_force,
                                 const Faisceau_Tubes& ) const;
  void ibc_localisee_velocity_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                   const IJK_Field_double& rho_field,
                                   DoubleTab& masse_fluide_cylindres,
                                   DoubleTab& volume_cylindres,
                                   DoubleTab& integrale_force,
                                   const Faisceau_Tubes& , double current_time) const;
  void ibc_localisee_velocity_cube_qdm(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                       const IJK_Field_double& rho_field,
                                       DoubleTab& masse_fluide_cylindres,
                                       DoubleTab& volume_cylindres,
                                       DoubleTab& integrale_force,
                                       const Faisceau_Tubes& , double current_time) const;
  void force_ibc_velocity_miroir(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                 const IJK_Field_double& rho_field,
                                 DoubleTab& masse_fluide_cylindres,
                                 DoubleTab& volume_cylindres,
                                 DoubleTab& integrale_force,
                                 const Faisceau_Tubes& ) const;
  void force_ibc_velocity_symetrie_plane(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                         const IJK_Field_double& rho_field,
                                         DoubleTab& masse_fluide_cylindres,
                                         DoubleTab& volume_cylindres,
                                         DoubleTab& integrale_force,
                                         const Faisceau_Tubes& ) const;
  void calcul_F_pression(const IJK_Field_double& pressure, const IJK_Field_double& vx,
                         DoubleTab& integrale_force_pression, DoubleTab& pression_teta,
                         const Faisceau_Tubes& ) const;
  void calcul_F_pression2(const IJK_Field_double& pressure, const IJK_Field_double& vx,
                          DoubleTab& integrale_force_pression, DoubleTab& pression_teta,
                          const Faisceau_Tubes& ) const;
  void ibc0_force_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                       const IJK_Field_double& rho_field,
                       DoubleTab& masse_fluide_cylindres,
                       DoubleTab& volume_cylindres,
                       DoubleTab& integrale_force,
                       const Faisceau_Tubes& ) const;
  void ibc_diffuse_force_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                              const IJK_Field_double& rho_field,
                              DoubleTab& masse_fluide_cylindres,
                              DoubleTab& volume_cylindres,
                              DoubleTab& integrale_force,
                              const Faisceau_Tubes& ) const;
  void ibc_localisee_force_cube(IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz,
                                const IJK_Field_double& rho_field,
                                DoubleTab& masse_fluide_cylindres,
                                DoubleTab& volume_cylindres,
                                DoubleTab& integrale_force,
                                const Faisceau_Tubes& , double current_time) const;
  void update(double current_time,
              double timestep);

  void reprendre_probleme(Entree& fichier);

  void sauvegarder_probleme(SFichier& f);
  void sauvegarder_pression(SFichier& f);

protected:
#if 0
  double oscillation_cylindre(double& az,
                              double& vz,
                              int tstep,
                              double position_initiale_cylindre_oscillant,
                              const double timestep,
                              const double masse_fluide_dans_cylindre,
                              const double vitesse_entree) const;
  void newmark_implicit_update(double& vz,
                               double& az,
                               DoubleTab& donnees_tubes_,
                               int tstep,
                               double position_initiale_cylindre_oscillant,
                               const double timestep,
                               const double vitesse_entree) const;
#endif
protected:
  REF(IJK_Splitting) ref_splitting_; // cette variable est initialisee dans la methode initialize()

  int lissage_;
  double epaisseur_lissage_;
  int n_P_;
  int champ_miroir_;
  int methode_IBC_;
  int solide_;
  double L_cube_;

  // force IBC : tableau a trois colonnes et N lignes contenant les forces selon x, y et z sur les N cylindres
  DoubleTab integrale_force_;
  DoubleTab integrale_force_post_proj_;
  DoubleTab integrale_force_pression_;
  DoubleTab pression_teta_;
  DoubleTab integrale_force_N_moins_1_;
  DoubleTab integrale_force_N_moins_1_post_proj_;
  DoubleTab masse_fluide_cylindres_;
  DoubleTab d_integrale_force_;
  DoubleTab d_integrale_force_post_proj_;
  DoubleTab volume_cylindres_;

  double rho_fluide_pour_adim_;
  double vitesse_pour_adim_;

  // L'ensemble des tubes...
  Faisceau_Tubes faisceau_;
};
#endif
