/****************************************************************************
* Copyright (c) 2022, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Domaine_ALE.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Domaine_ALE_included
#define Domaine_ALE_included

#include <Domaine.h>
#include <TRUSTLists.h>
#include <Champs_front.h>
#include <Champ_P1NC.h>
#include <Beam_model.h>
#include <Structural_dynamic_mesh_model.h>
#include <Champs_front_ALE_projection.h>
#include <TRUST_Ref.h>

class Equation_base;
class Domaine_dis;
class Beam_model;
class Structural_dynamic_mesh_model;
//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Domaine_ALE
//    Cette classe est un interprete qui sert a lire l'attribut axi.
//    Directive:
//        Domaine_ALE
//    Cette directive optionelle permets de faire les calculs en
//    coordonnees cylindriques. En l'absence de cette directive les calculs
//    se font en coordonnees cartesiennes.
// .SECTION voir aussi
//    Interprete Objet_U
//////////////////////////////////////////////////////////////////////////////
class Domaine_ALE : public Domaine
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Domaine_ALE);

public :
  Domaine_ALE();
  ~Domaine_ALE();
  inline const double& get_dt() const;
  void set_dt(double& dt) override;
  inline const DoubleTab& vitesse() const;
  inline DoubleTab& vitesse_faces();
  inline const DoubleTab& vitesse_faces() const;
  void mettre_a_jour (double temps, Domaine_dis&, Probleme_base&) override;
  void update_after_post(double temps) override;
  void initialiser (double temps, Domaine_dis&, Probleme_base&) override;
  DoubleTab calculer_vitesse(double temps,Domaine_dis&, Probleme_base&, bool&);
  DoubleTab& calculer_vitesse_faces(DoubleTab&, int, int, IntTab&);
  void reading_vit_bords_ALE(Entree& is);
  void reading_solver_moving_mesh_ALE(Entree& is);
  void reading_beam_model(Entree& is);
  void reading_projection_ALE_boundary(Entree& is);
  void reading_structural_dynamic_mesh_model(Entree& is);
  void  update_ALE_projection(double, Nom&, Champ_front_ALE_projection& , int);
  void  update_ALE_projection(const double);
  DoubleTab& laplacien(Domaine_dis&, Probleme_base&, const DoubleTab&, DoubleTab&);
  int update_or_not_matrix_coeffs() const;
  void update_ALEjacobians(DoubleTab&, DoubleTab&, int);
  void resumptionJacobian(DoubleTab&, DoubleTab&);
  inline const DoubleTab& getOldJacobian();
  inline const DoubleTab& getOldJacobian() const;
  inline const DoubleTab& getNewJacobian();
  inline const DoubleTab& getNewJacobian() const;
  inline int getMeshMotionModel() const ;


  DoubleVect interpolationOnThe3DSurface(const double& x, const double& y, const double& z, const DoubleTab& u, const DoubleTab& R) const;
  void initializationBeam (double velocity) ;
  double computeDtBeam(Domaine_dis&);
  const DoubleTab& getBeamDisplacement(int i) const;
  const DoubleTab& getBeamRotation(int i) const;
  const int& getBeamDirection() const;
  DoubleVect& getBeamVelocity(const double& tps, const double& dt);
  const int& getBeamNbModes();
  void computeFluidForceOnBeam();
  const DoubleVect& getFluidForceOnBeam();
  Equation_base& getEquation() ;
  inline void associer_equation(const Equation_base& une_eq);
  const DoubleVect& getMeshPbPressure() const ;
  const DoubleVect& getMeshPbVonMises() const ;
  const DoubleTab& getMeshPbForceFace() const ;
protected:

  double dt_;
  DoubleTab ALE_mesh_velocity;
  DoubleTab vf; //faces velocity
  IntTab som_faces_bords;
  SolveurSys solv;
  Matrice_Morse_Sym mat;
  Champs_front les_champs_front;
  int nb_bords_ALE;
  Bords les_bords_ALE;
  int update_or_not_matrix_coeffs_; //=1 in case of zero ALE boundary/mesh velocity, =0 otherwise (see Domaine_ALE::calculer_vitesse).
  DoubleTab ALEjacobian_old; // n
  DoubleTab ALEjacobian_new; // n+1
  int resumption; //1 if resumption of calculation else 0
  Beam_model *beam; // Mechanical model: a beam model
  Structural_dynamic_mesh_model *str_mesh_model; // Fictitious structural model for mesh motion
  REF(Equation_base) eq;
  DoubleVect fluidForceOnBeam; //Fluid force acting on the IFS boundary
  Champs_front_ALE_projection field_ALE_projection_; // Definition of the modes of vibration in view of projection of the IFS force
  Noms name_ALE_boundary_projection_; // Names of the ALE boundary where the projection is computed
  bool associate_eq;
  double tempsComputeForceOnBeam; // Time at which the fluid force acting on the Beam is computed.
  int meshMotionModel_ = 0 ; // Model for ALE mesh motion: 0 = Laplacien, 1 = Structural_dynamics
  void solveDynamicMeshProblem_(const double temps, const DoubleTab& imposedVelocity, const IntVect& imposedVelocityTag,
                                DoubleTab& outputMeshVelocity, const int nbSom, const int nbElem, const int nbSomElem,
                                const IntTab& sommets, const int nbFace, const int nbSomFace, const IntTab& face_sommets) ;

};


inline const DoubleTab& Domaine_ALE::vitesse() const
{
  return ALE_mesh_velocity;
}

inline const double& Domaine_ALE::get_dt() const
{
  return dt_;
}

inline DoubleTab& Domaine_ALE::vitesse_faces()
{
  return vf;
}
inline const DoubleTab& Domaine_ALE::vitesse_faces() const
{
  return vf;
}
inline int Domaine_ALE::update_or_not_matrix_coeffs() const
{
  return update_or_not_matrix_coeffs_;
}
inline const DoubleTab& Domaine_ALE::getOldJacobian()
{
  return ALEjacobian_old;
}
inline const DoubleTab& Domaine_ALE::getOldJacobian() const
{
  return ALEjacobian_old;
}

inline const DoubleTab& Domaine_ALE::getNewJacobian()
{
  return ALEjacobian_new;
}
inline const DoubleTab& Domaine_ALE::getNewJacobian() const
{
  return ALEjacobian_new;
}

inline void Domaine_ALE::associer_equation(const Equation_base& une_eq)
{
  eq = une_eq;
}

inline int Domaine_ALE::getMeshMotionModel() const
{
  return meshMotionModel_ ;
}

#endif
