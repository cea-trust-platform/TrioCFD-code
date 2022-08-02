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
// File      : Init_spectral.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Chouippe_included
#define Chouippe_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <fftw3.h>
#pragma GCC diagnostic push
#if __GNUC__ > 5
#pragma GCC diagnostic ignored "-Wsuggest-override"
#endif
#include <fftw3-mpi.h>
#pragma GCC diagnostic pop
#include <vector>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Chouippe
//
// <Description of class Chouippe>
//
/////////////////////////////////////////////////////////////////////////////


int unidentified_function(int argc, char **argv);
void Create_file(std::string const& fichname, std::string const& name);
void Add_in_file(std::string const& fichname, std::string const& information);
void Add_in_file(std::string const& fichname, int const information);
void Add_in_file(std::string const& fichname, ArrOfInt const& information);
void Add_in_file(std::string const& fichname, std::string const& info, int const information);

// void fftw_org_real_data_MPI_transform(int argc, char **argv);
std::string fftw_org_real_data_MPI_transform();
// int fftw_org_unidentified_function(int argc, char **argv);
std::string fftw_org_unidentified_function();


class forcage_spectral : public Objet_U
{

  Declare_instanciable( forcage_spectral ) ;

public :
  // Constructeurs
  // forcage_spectral();
  // Destructeur;
  // ~forcage_spectral();

  void fftw_org_multi_D_MPI_DFT_real_data();
  double get_compteur0();
  double get_compteur1();
  double get_compteur2();
  void set_nk_kmin_kmax(const int number_k, const double kmin, const double kmax);
  void set_spectral_domain();
  void set_a_force();

protected :
  //void compute_inital_velocity_spectral(FixedVector<IJK_Field_double, 3>);
private:

  double k_max,k_min;
  int nk;
  int nk_tot;
  // std::vector<double> kx(), ky(),kz();
  ArrOfDouble kx, ky,kz;
  // std::vector<double> fx(),fy(),fz();
  ArrOfDouble fx,fy,fz;

  double compteur[3];
};





#endif /* Init_spectral_included */
