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

#include <List.h>
#include <Interprete_geometrique_base.h>
#include <Nom.h>
#include <Bords.h>
#include <Motcle.h>

Declare_liste(DoubleTab);
/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Beam_model
//
// <Description of class Beam_model>
//
/////////////////////////////////////////////////////////////////////////////

class Beam_model : public Interprete_geometrique_base
{

  Declare_instanciable( Beam_model ) ;

public :
  Entree& interpreter_(Entree&);
  inline const int& getNbModes() const;
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
  inline const bool& getActivate() const;
  inline void setActivate(const bool&) ;
  inline const bool& getTimeScheme() const;
  inline void setTimeScheme(const bool&) ;
  void readInputMassStiffnessFiles (Nom& masse_and_stiffness_file_name);
  void readInputAbscFiles (Nom& absc_file_name);
  void readInputModalDeformation(Noms& modal_deformation_file_name);
  void readInputCIFile(Nom& CI_file_name);
  //void interpolationOnThe3DSurface(const Bords& les_bords_ALE);
  DoubleVect interpolationOnThe3DSurface(const double& x, const double& y, const double& z, const DoubleTab& u, const DoubleTab& R) const;
  void initialization(double velocity);
  DoubleVect& NewmarkSchemeFD (const double& dt, const DoubleVect& fluidForce);
  DoubleVect& NewmarkSchemeMA (const double& dt, const DoubleVect& fluidForce);
  DoubleVect& getVelocity(const double& dt, const DoubleVect& fluidForce);
  inline double soundSpeed();
  inline  double getMass(int i);
  inline  double getStiffness(int i);
  inline  const DoubleTab& getDisplacement(int i) const;
  inline  const DoubleTab& getRotation(int i) const;

protected :
  int nbModes_;
  int direction_; //x=0, y=1, z=2
  double young_; // Young module
  double rho_; // solid density
  bool activate_=false;
  DoubleVect mass_;
  DoubleVect stiffness_;
  DoubleVect damping_;
  DoubleVect abscissa_;
  LIST(DoubleTab) u_;
  LIST(DoubleTab) R_;
  DoubleVect qSpeed_;
  DoubleVect qHalfSpeed_;
  DoubleVect qAcceleration_;
  DoubleVect qDisplacement_;
  bool timeScheme_=true;
  //DoubleTab phi3D_;

};
inline const int& Beam_model::getNbModes() const
{
  return nbModes_;
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

inline const bool& Beam_model::getActivate() const
{
  return activate_;
}
inline void Beam_model::setActivate(const bool& activate)
{
  activate_=activate;
}
inline const bool& Beam_model::getTimeScheme() const
{
  return timeScheme_;
}
inline void Beam_model::setTimeScheme(const bool& timeScheme)
{
  timeScheme_=timeScheme;
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

const DoubleTab& Beam_model::getDisplacement(int i) const
{
  assert(i<nbModes_);
  return u_(i);
}
const DoubleTab& Beam_model::getRotation(int i) const
{
  assert(i<nbModes_);
  return R_(i);

}


#endif /* Beam_model_included */
