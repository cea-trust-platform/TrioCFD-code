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
// File      : Force_ph.h
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef force_ph_included
#define force_ph_included


#include <iostream>
#include <string>
#include <vector>
#include <Force_sp.h>
#include <FixedVector.h>
// #include <fftw3.h>
#include <IJK_Splitting.h>
#include <communications.h>
#include <Objet_U.h>

#include <IJK_Field.h>

class Force_ph : public Objet_U
{

  Declare_instanciable_sans_constructeur_ni_destructeur( Force_ph ) ;

public:

  Force_ph();
  ~Force_ph();

  // void initialise(int nproc_tot, int ni, int nj, int nk, int nl,int nm,int nn,
  // 		double Lx, double Ly, double Lz, double Ox,double Oy,double Oz, int momin, int momax, double kmin, double kmax,
  // 	  std::string nom_fichier, const IJK_Splitting& splitting
  // 	);
  void initialise(int nproc_tot, int ni, int nj, int nk, int nl,int nm,int nn,
                  double Lx, double Ly, double Lz, double Ox,double Oy,double Oz, int momin, int momax, double kmin, double kmax,
                  std::string nom_fichier, const IJK_Splitting& splitting,
                  int a_i_offset, int a_j_offset, int a_k_offset
                 );
  void from_spect_to_phys2(const std:: vector <double >& coeff_force);
  void from_spect_to_phys_opti2(ArrOfDouble& coeff_force);
  void from_spect_to_phys_opti2_advection(ArrOfDouble& coeff_force, const ArrOfDouble& advection_length);
  void from_spect_to_phys_bis(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
  void set_zero();
  void cheat_function();
  // void from_spect_to_phys_fftw(fftw_complex* in);
  void write(std::string nom_fichier_sortie, double t);
  void write_separate(std::string nom_fichier_sortie, double t);
  void write_offset_index_position( const IJK_Splitting& my_splitting);
  void compute_energie();
  double get_energie();
  FixedVector<IJK_Field_double, 3> get_force_attribute();
  FixedVector<IJK_Field_double, 3>& get_force_attribute2();

  void gbz_gather(FixedVector<IJK_Field_double, 3> force_);

private:
  int nproc_tot;

  int ni,nj,nk,n_ijk;
  int nl,nm,nn,n_lmn;
  double Lx, Ly, Lz;
  double Ox, Oy, Oz;
  double kmin,kmax;
  int momin,momax;
  FixedVector<IJK_Field_double, 3> force_;
  std::vector<std::vector< std:: vector < double > > > force;
  double energie;

  int nproc_i,nproc_j,nproc_k;
  int i_offset;
  int j_offset;
  int k_offset;


};

std::vector< std::vector< std:: vector <double >>> set_dimensions(std::vector< std::vector< std:: vector <double >>> the_vector, int dim_one, int dim_two, int dim_three);

// FixedVector<IJK_Field_double, 3> set_to_zero(FixedVector<IJK_Field_double, 3> vector)
// {
// 		for (int dir=0; dir<3; dir++)
// 			for (int i=0; i<vector[0].ni(); i++)
// 				for (int j=0; j<vector[1].nj(); j++)
// 					for (int k=0; k<vector[2].nk(); k++)
// 						vector[dir](i,j,k) = 0.;
// 		return vector;
// }

#endif
