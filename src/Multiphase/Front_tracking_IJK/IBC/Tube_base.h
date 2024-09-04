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
// File      : Tube_base.h
// Directory : $IJK_ROOT/src/IBC
//
/////////////////////////////////////////////////////////////////////////////
#include <Objet_U.h>
#include <Linear_algebra_tools.h>
#include <Noms.h>
#include <IJK_Splitting.h>
#include <TRUST_Vector.h>
#include <TRUST_Deriv.h>
#include <TRUST_Ref.h>

class Tube_base : public Objet_U
{
  Declare_base(Tube_base);
public:
  void initialize(const IJK_Splitting&);
  const Vecteur3& get_current_pos() const
  {
    return pos_;
  }
  const Vecteur3& get_current_velocity() const
  {
    return v_;
  }
  const Vecteur3& get_current_acceleration() const
  {
    return acceleration_;
  }
  double get_omega() const
  {
    return omega_;
  }
  double get_rayon() const
  {
    return rayon_;
  }
  // current_time: temps physique actuel du calcul
  // force appliquee est la force appliquee par le fluide sur le tube
  // dt intervalle de temps depuis le current_time precedent
  virtual void update_vitesse_position(double current_time, double dt, const Vecteur3& force_appliquee) = 0;
  const Nom& nom_du_tube() const
  {
    return nom_;
  }
protected:
  REF(IJK_Splitting) ref_splitting_;
  double omega_; // vitesse de rotation du tube en rad/s
  double rayon_; // rayon du tube
  Vecteur3 pos_;
  Vecteur3 v_;
  Vecteur3 acceleration_;
  Nom nom_;
};


class Tube_impose : public Tube_base
{
  Declare_instanciable(Tube_impose);
public:
  void update_vitesse_position(double current_time, double dt, const Vecteur3& force_appliquee) override;
protected:
  // Troix chaines de caracteres (fonctions position(t))
  Noms expressions_position_;
};


class Tube_libre : public Tube_base
{
  Declare_instanciable(Tube_libre);
public:
  void update_vitesse_position(double current_time, double dt, const Vecteur3& force_appliquee) override;
protected:
  Vecteur3 position_equilibre_;
  ArrOfInt blocage_;
  double rho_cylindre_, amortissement_, raideur_;

};
#if 0
class Tube_couple : public Tube_base
{
};
#endif

class Faisceau_Tubes : public VECT(OWN_PTR(Tube_base))
{
  Declare_instanciable(Faisceau_Tubes);
};
