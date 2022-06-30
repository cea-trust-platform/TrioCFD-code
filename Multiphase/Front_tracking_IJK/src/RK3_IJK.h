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
// File      : RK3_IJK.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#ifndef RK3_IJK_H
#define RK3_IJK_H
#include <FixedVector.h>
#include <IJK_Field.h>
#include <IJK_Splitting.h>
#include <OpDiffTurbIJK.h>
#include <OpConvIJKQuickSharp.h>
#include <Multigrille_Adrien.h>
#include <Interprete.h>
#include <IJK_Lata_writer.h>
#include <Linear_algebra_tools.h>
#include <Couplage_Tubes_IBC.h>
#include <OpConvIJKQuickScalar.h>
#include <Boundary_Conditions.h>


#define ACTIVE 0
#define DESACTIVE 1
#define DESPRES_SAMUEL 0
#define DESPRES_DELLACHERIE 1
#define AMONT 2


class IJK_problem_double : public Interprete
{
  Declare_instanciable(IJK_problem_double);
public:
  Entree& interpreter(Entree&) override;
  void run();
  void euler_time_step();
protected:
  void reprendre_probleme(const char *fichier_reprise);
  void initialise();
  void terme_source_gravite(IJK_Field_double& dv, int k_index, int dir) const;
  void recalculer_rho_de_levelset(const IJK_Field_double& levelset,
                                  IJK_Field_double& rho,
                                  int epaisseur_joint) const;
  void recalculer_mu_de_rho(const IJK_Field_double& rho,
                            IJK_Field_double& molecular_mu,
                            int epaisseur_joint) const;
  void euler_explicit_update(const IJK_Field_double& dv, IJK_Field_double& v,
                             const int k_layer) const;
  const IJK_Grid_Geometry& get_geometry() const
  {
    return splitting_.get_grid_geometry();
  }
  double intermediate_dt(int step) const
  {
    const double intermediate_tstep[3] = { 1./3., 5./12., 1./4. };
    assert(step>=0 && step<3);
    return intermediate_tstep[step] * timestep_;
  }

  void sauvegarder_probleme(const char *fichier_sauvegarde);
  void sauvegarder_pression();
  void posttraiter_champs_instantanes(const char * lata_name, double time);
  void scalar_convection_op_Triton(const IJK_Field_double& field,
                                   const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                   IJK_Field_double& result);
  void scalar_convection_op_Despres_Samuel(IJK_Field_double& field,
                                           const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz);
  void scalar_convection_op_Despres_Dellacherie(IJK_Field_double& field,
                                                const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz);
  void scalar_convection_op_DESPRES_alpha_AMONT_rho(IJK_Field_double& alpha_field, IJK_Field_double& rho_field,
                                                    const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz);
  void scalar_convection_op_Amont_rho(IJK_Field_double& field,
                                      const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz);
  void calculer_rho2(const IJK_Field_double& alpha_field_, IJK_Field_double& rho_field_bis_);
  void calcul_integrale_rho( IJK_Field_double& rho_field_) const;

  int dt_post_;
  Motcles liste_post_instantanes_; // liste des champs instantanes a postraiter
  // Pour numeroter les fichiers .lata il faut compter combien on en a ecrit:
  int compteur_post_instantanes_;
  Nom check_stop_file_; // Nom du fichier stop
  int dt_sauvegarde_;
  Nom nom_sauvegarde_;
  Nom nom_reprise_;
  // Le jeu de donnees doit fournir soit des fichiers de reprise: ..
  Nom fichier_reprise_rho_;
  int timestep_reprise_rho_;
  Nom fichier_reprise_alpha_;
  int timestep_reprise_alpha_;
  Nom fichier_reprise_vitesse_;
  int timestep_reprise_vitesse_;
  // ... soit des expressions f(x,y,z)
  Nom expression_rho_initiale_;
  Nom expression_alpha_initiale_;
  Noms expression_vitesse_initiale_; // on attend trois expressions

  IJK_Splitting splitting_;
  ArrOfDouble_with_ghost delta_z_local_;

  Boundary_Conditions boundary_conditions_;

  // Inconnues du probleme (a sauvegarder et a reprendre)
  // Velocity field:
  FixedVector<IJK_Field_double, 3> velocity_;
  FixedVector<IJK_Field_double, 3> velocity_tmp_;

  // Fonction level set
  //IJK_Field_double fonction_levelset_;

  // Masse volumique:
  IJK_Field_double rho_field_;

  // Epaisseur de joint du champ de vitesse pour le calcul du champ ibc mirroir
  int ghost_size_pour_ibc_;
  IJK_Field_double rho_field_bis_;

  // Taux de vide
  IJK_Field_double alpha_field_;



  // Temporary storage for the derivative
  FixedVector<IJK_Field_double, 3> d_velocity_;

  IJK_Field_double d_rho_field_;
  //IJK_Field_double d_fonction_levelset_;

  // Pressure field
  IJK_Field_double pressure_;
  // Molecular diffusivity (see diffusion operator)
  IJK_Field_double molecular_mu_;
  // right hand side for pressure solver
  IJK_Field_double pressure_rhs_;
  // Operators and pressure solver
  OpDiffIJK_double velocity_diffusion_op_;
  OpConvQuickSharpIJK_double velocity_convection_op_;
  OpConvIJKQuickScalar_double scalar_convection_op_;

  Multigrille_Adrien poisson_solver_;
  // Simulation parameters
  int nb_timesteps_;
  double timestep_;
  double source_pressure_gradient_;
  double mu_liquide_, mu_gaz_;
  double timestep_facsec_;

  double vitesse_entree_;
  double rho_liquide_;
  double rho_gaz_;
  ArrOfDouble gravite_; // vecteur de taille 3 a lire dans le jeu de donnees
  int check_divergence_;

  Couplage_Tubes_IBC couplage_tubes_ibc_;

  double current_time_;

  int champ_miroir_pour_RHS_;
  int choix_schema_transport_rho_;
};
#endif
