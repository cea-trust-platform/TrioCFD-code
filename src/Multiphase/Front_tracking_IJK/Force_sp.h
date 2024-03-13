/****************************************************************************
* Copyright (c) 2019, CEA
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
// File      : Force_sp.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef force_sp_included
#define force_sp_included

// #include <iostream>
#include <string>
#include <vector>

#include <FixedVector.h>
#include <IJK_Splitting.h>
#include <Objet_U.h>
#include <string>
#include <iostream>
#include <math.h>
// #include <fftw3.h>
#include <Force_ph.h>
#include <Random_process.h>

#include <IJK_Field.h>
// #include <init_forcage_THI.h>

class Force_sp : public Objet_U
{

  Declare_instanciable_sans_constructeur_ni_destructeur( Force_sp ) ;

public:

  Force_sp();
  ~Force_sp();
  void initialise(int nl, int nn, int nm, int momin, int momax, double kmin, double kmax,std::string nom_fichier);
  void initialise(int nl, int nn, int nm, int momin, int momax, double kmin, double kmax, double amplitude, std::string nom_fichier);
  // void initialise(int nl, int nn, int nm, double kmin, double kmax, double amplitude);
  // void initialise(int nl, int nn, int nm, double kmin, double kmax);
  // void compute_step2(const std:: vector <double >& process);
  void compute_step2(ArrOfDouble& process_flt);
  void field_advection(const ArrOfDouble& advection_length, const double dt, const int it);
  void set_zero();
  void compute_dirac_board();
  void compute_dirac_point_div_nulle();
  void compute_dirac_point_uniZ();
  void compute_dirac_point_uniY();
  void compute_dirac_point_uniX_alongX();
  void compute_dirac_point_uniX_alongY();
  void compute_dirac_point();
  void compute_door_rope();
  void compute_door_cube();
  void compute_force_kappa();
  void compute_diracs_for_cos_squarred();
  void compute_diracs_for_t_times_cos_squarred(double time);
  void write(std::string nom_fichier_sortie, double time);
  void write_separate(std::string nom_fichier_sortie, double t);
  void compute_energie();
  double get_energie();
  double get_force(int cpx, int dir, int ind);
  std::vector< std::vector< std:: vector <double >>> get_coeff();
  // std:: vector <double > get_coeff_flt();
  ArrOfDouble& get_coeff_flt();
private:

  int nl,nm,nn,n_lmn;
  int momin,momax;
  double amplitude;
  double kmin,kmax;
  double energie;

  std::vector< std::vector< std:: vector <double >>> force;
  ArrOfDouble force_flt;
};



#endif
