/****************************************************************************
* Copyright (c) 2021, CEA
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
// File      : Beam_model.h
// Directory : $BEAM_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Beam_model_included
#define Beam_model_included

#include <TRUST_List.h>
#include <Interprete_geometrique_base.h>
#include <Nom.h>
#include <Bords.h>
#include <Motcle.h>
#include <SFichier.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Beam_model
//
// <Description of class Beam_model>
// Fluid-structure interaction model.
// Reduced mechanical model: a beam model (1d).
// Resolution based on a modal analysis.
// Temporal discretization: Newmark finite differences, Newmark mean acceleration,  or Hilber-Hughes-Taylor (HHT)
// ALE coupling with the reduced beam model
/////////////////////////////////////////////////////////////////////////////

class Beam_model : public Interprete_geometrique_base
{

  Declare_instanciable_sans_constructeur_ni_destructeur( Beam_model ) ;

public :
  Beam_model();
  ~Beam_model();
  Entree& interpreter_(Entree&) override;
  inline const int& getNbModes() const;
  inline double getTime() const;
  inline void setNbModes(const int&) ;
  inline const int& getDirection() const;
  inline void setDirection(const int&) ;
  inline const double& getYoung() const;
  inline void setYoung(const double&);
  inline const double& getRhoBeam() const;
  inline void setRhoBeam(const double&) ;
  inline const Nom& getBeamName() const;
  inline void setBeamName(const Nom&) ;
  inline const Nom& getFileName() const;
  inline void setFileName(const Nom&) ;
  inline const Nom& getTimeScheme() const;
  inline void setTimeScheme(const Nom&, double&) ;
  inline void setOutputPosition1D(const DoubleVect&) ;
  inline void setOutputPosition3D(const DoubleTab&) ;
  void readInputMassStiffnessFiles (Nom& masse_and_stiffness_file_name);
  void readInputAbscFiles (Nom& absc_file_name);
  void readInputModalDeformation(Noms& modal_deformation_file_name);
  void readInputCIFile(Nom& CI_file_name);
  void readRestartFile(Nom& Restart_file_name);
  //void interpolationOnThe3DSurface(const Bords& les_bords_ALE);
  DoubleVect interpolationOnThe3DSurface(const double& x, const double& y, const double& z, const DoubleTab& u, const DoubleTab& R) const;
  void initialization(double displacement);
  void initialization();
  DoubleVect& NewmarkSchemeFD (const double& dt);
  DoubleVect& NewmarkSchemeMA (const double& dt);
  DoubleVect& TimeSchemeHHT (const double& dt);
  DoubleVect& getVelocity(const double& tps, const double& dt);
  inline double soundSpeed();
  inline  double getMass(int i);
  inline  double getStiffness(int i);
  inline  const DoubleTab& getDisplacement(int i) const;
  inline  const DoubleTab& getRotation(int i) const;
  void saveBeamForRestart() const;
  void printOutputBeam1D(bool first_writing=false) const;
  void printOutputBeam3D(bool first_writing=false) const;
  void printOutputFluidForceOnBeam(bool first_writing=false) const;
  inline void setFluidForceOnBeam(const DoubleVect&);
  inline void setTempsComputeForceOnBeam(const double&);
  inline const double& getTempsComputeForceOnBeam() const;
  void setCenterCoordinates(const double&,const double&, const double&);
protected :
  int nbModes_;
  int direction_; //x=0, y=1, z=2
  double young_; // Young module
  double rho_; // solid density
  DoubleVect mass_; //diagonal modal mass matrix
  DoubleVect stiffness_; //diagonal modal stiffness_ matrix
  DoubleVect damping_; //diagonal modal damping_ matrix
  DoubleVect abscissa_; // coordinates of the Beam
  LIST(DoubleTab) u_;  // Beam modal displacement constituting the modal shape
  LIST(DoubleTab) R_; // Beam modal rotation constituting the modal shape
  DoubleVect qSpeed_; // Beam 1d speed
  DoubleVect qHalfSpeed_; // Beam 1d speed at dt/2
  DoubleVect qAcceleration_; //Beam 1d acceleration
  DoubleVect qDisplacement_; //Beam 1d displacement
  Nom timeScheme_; // Time discretization scheme
  //DoubleTab phi3D_;
  double temps_;
  DoubleVect output_position_1D_; //post-processing of the 1d position of the points (points on the Beam)
  DoubleTab output_position_3D_; // post-processing of the 3d position of the points (points on the 3d surface)
  Nom beamName_;
  double alpha_;
  DoubleVect fluidForceOnBeam_; //Fluid force acting on the IFS boundary
  double tempsComputeForceOnBeam_;
  double x0_; //x-coordinate of the center of the Beam base
  double y0_; //y-coordinate of the center of the Beam base
  double z0_; //z-coordinate of the center of the Beam base
  mutable SFichier displacement_out_1d_; //output files of the displacement
  mutable SFichier speed_out_1d_; //output files of the speed
  mutable SFichier acceleration_out_1d_; //output files of the acceleration
  mutable SFichier displacement_out_3d_; //output files of the displacement
  mutable SFichier speed_out_3d_; //output files of the speed
  mutable SFichier acceleration_out_3d_; //output files of the acceleration
  mutable SFichier fluidForceOnBeam_out_; //output files of the fluid force acting on the IFS boundary
};
inline const int& Beam_model::getNbModes() const
{
  return nbModes_;
}
inline double Beam_model::getTime() const
{
  return temps_;
}
inline void Beam_model::setNbModes(const int& modes)
{
  nbModes_=modes;
}
inline const int& Beam_model::getDirection() const
{
  return direction_;
}
inline void Beam_model::setDirection(const int& direction)
{
  direction_=direction;
}
inline const double& Beam_model::getYoung() const
{
  return young_;
}
inline void Beam_model::setYoung(const double& young)
{
  young_ = young;
}
inline const double& Beam_model::getRhoBeam() const
{
  return rho_;
}
inline void Beam_model::setRhoBeam(const double& rho)
{
  rho_=rho;
}

inline const Nom& Beam_model::getTimeScheme() const
{
  return timeScheme_;
}
inline void Beam_model::setTimeScheme(const Nom& timeScheme, double& alpha)
{
  timeScheme_=timeScheme;
  alpha_=alpha;
}
inline double Beam_model::soundSpeed()
{
  return sqrt(young_/rho_);
}

inline double  Beam_model::getMass(int i)
{
  assert(i<nbModes_);
  return mass_[i];
}
inline double  Beam_model::getStiffness(int i)
{
  assert(i<nbModes_);
  return stiffness_[i];
}

inline const DoubleTab& Beam_model::getDisplacement(int i) const
{
  assert(i<nbModes_);
  return u_(i);
}
inline const DoubleTab& Beam_model::getRotation(int i) const
{
  assert(i<nbModes_);
  return R_(i);

}

inline void Beam_model::setOutputPosition1D(const DoubleVect& pos)
{
  output_position_1D_=pos;

}
inline void Beam_model::setOutputPosition3D(const DoubleTab& pos)
{
  output_position_3D_=pos;

}

inline const Nom& Beam_model::getBeamName() const
{
  return beamName_;
}
inline void Beam_model::setBeamName(const Nom& value)
{
  beamName_=value;
}
inline void Beam_model::setFluidForceOnBeam(const DoubleVect& force)
{

  fluidForceOnBeam_=force;
}

inline void Beam_model::setTempsComputeForceOnBeam(const double& tps)
{
  tempsComputeForceOnBeam_=tps;
}

inline const double& Beam_model::getTempsComputeForceOnBeam() const
{
  return tempsComputeForceOnBeam_;
}

#endif /* Beam_model_included */
