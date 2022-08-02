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
// File      : Init_spectral.cpp
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <forcage_spectral.h>
#include <string>
#include <fstream>
#include <algorithm>    // std::random_shuffle

// GAB pour avoir une distribution Gaussienne, il nous faut la librairie :
#include <random>

Implemente_instanciable( forcage_spectral, "forcage_spectral", Objet_U ) ;

Sortie& forcage_spectral::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& forcage_spectral::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}



double my_function_d_2(const int i, const int j, const int k)
{
  double b(2.9);
  return b;
  // return i/(j*k);
}

// 6.5 Multi-dimensional MPI DFTs of Real Data
void forcage_spectral::fftw_org_multi_D_MPI_DFT_real_data()
{
  const ptrdiff_t L = 7, M = 8, N = 9;
  fftw_plan plan;
  double *rin;
  fftw_complex *c_out;
  ptrdiff_t alloc_local, local_n0, local_0_start;
  int i, j, k;

  // fftw_mpi_init();

  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_3d(L, M, N/2+1, MPI_COMM_WORLD,
                                       &local_n0, &local_0_start);
  rin = fftw_alloc_real(2 * alloc_local);
// double[2]
  c_out = fftw_alloc_complex(alloc_local);

  // fftw_complex *ux_f = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex)*kmax*kmax*kmax );

  /* create plan for out-of-place r2c DFT */
  plan = fftw_mpi_plan_dft_r2c_3d(L, M, N, rin, c_out, MPI_COMM_WORLD,
                                  FFTW_MEASURE);
  // int rank(3);
  // int n_pointeur[3];
  // n_pointeur[0]=5,n_pointeur[1]=10,n_pointeur[2]=15;
  // fftw_complex in_pointer;
  // fftw_complex out_pointer;
  // plan2 = fftw_mpi_plan_dft_r2c(rank, n_pointer, in_pointer, out_pointer, flags)

  /* initialize rin to some function my_func(x,y,z) */
  for (i = 0; i < local_n0; ++i)
    {
      compteur[0]+=1;
      for (j = 0; j < M; ++j)
        {
          compteur[1]+=1;
          for (k = 0; k < N; ++k)
            {
              compteur[2]+=1;
              // ux_f[][0] =
              // ux_f[][1] =
              rin[(i*M + j) * (2*(N/2+1)) + k]= my_function_d_2((int)local_0_start+i, j, k);
              // rin[(i*M + j) * (2*(N/2+1)) + k][1] = my_function_d_2(local_0_start+i, j, k);
            }
        }
    }
  /* compute transforms as many times as desired */
  fftw_free(rin); // Dans la doc, ils disent que c'est que si on utilise malloc, mais bon, la ca plante s'il n'y est pas...
  // fftw_free(c_out);
  fftw_execute(plan);
  // fftw_free(rin); // Dans la doc, ils disent que c'est que si on utilise malloc, mais bon, la ca plante s'il n'y est pas...

  fftw_destroy_plan(plan);
}

void forcage_spectral::set_nk_kmin_kmax(const int number_k,
                                        const double kmin, const double kmax)
{
  nk = number_k;
  k_max = kmax;
  k_min = kmin;
}

static int myPow(int x, unsigned int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = myPow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}

void forcage_spectral::set_spectral_domain()
{
  int i,j,k,ind;
  nk_tot = myPow(nk,3)-1;
  kx.resize_array(nk_tot);
  ky.resize_array(nk_tot);
  kz.resize_array(nk_tot);
  for (i=0; i<nk; i++)
    {
      for (j=0; j<nk; j++)
        {
          for (k=0; k<nk; k++)
            {
              ind = i*nk*nk+j*nk+k;
              kx[ind] = k_min + i*(k_max/nk);
              ky[ind] = k_min + j*(k_max/nk);
              kz[ind] = k_min + k*(k_max/nk);
            }
        }
    }
}

void forcage_spectral::set_a_force()
{
  const double pi(3.14159265358979);
  int i,j,k,ind;
  nk_tot = myPow(nk,3)-1;
  fx.resize_array(nk_tot);
  fy.resize_array(nk_tot);
  fz.resize_array(nk_tot);
  for (i=0; i<nk; i++)
    {
      for (j=0; j<nk; j++)
        {
          for (k=0; k<nk; k++)
            {
              ind = i*nk*nk+j*nk+k;
              fx[ind] = cos(2*pi*kx[ind]);
              fy[ind] = cos(2*pi*ky[ind]);
              fz[ind] = cos(2*pi*kz[ind]);
            }
        }
    }

}
////////////////////////////////////////////////////////////////////////////////
///////////// RELATED FUNCTIONS ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void Create_file(std::string const& fichname, std::string const& name)
{

  std::ofstream myFlux(fichname.c_str());
  if(myFlux)
    {
      myFlux << "----------- " << name << " --------" << endl;
    }
}

void Add_in_file(std::string const& fichname, std::string const& information)
{

  std::ofstream myFlux(fichname.c_str(), ios::app);
  if(myFlux)
    {
      myFlux << information << endl;
    }
}
void Add_in_file(std::string const& fichname, int const information)
{

  std::ofstream myFlux(fichname.c_str(), ios::app);
  if(myFlux)
    {
      myFlux << information << endl;
    }
}
void Add_in_file(std::string const& fichname, std::string const& info, int const information)
{

  std::ofstream myFlux(fichname.c_str(), ios::app);
  if(myFlux)
    {
      myFlux << info << " : "<<information << endl;
    }
}




// int fftw_org_unidentified_function(int argc, char **argv)
std::string fftw_org_unidentified_function()
{
  const ptrdiff_t N0 = 2, N1 = 5;
  fftw_plan plan;
  fftw_complex *data;
  ptrdiff_t alloc_local, local_n0, local_0_start;

  // MPI_Init(&argc, &argv);
  // fftw_mpi_init();

  /* get local data size and allocate */
  // you call fftw_mpi_local_size_2d to find out what portion of the array resides
  // on each processor, and how much space to allocate. Here, the portion of the
  // array on each process is a local_n0 by N1 slice of the total array, starting
  // at index local_0_start. The total number of fftw_complex numbers to
  // allocate is given by the alloc_local return value, which may be greater than
  // local_n0 * N1
  // //////////////////////////////
  alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                       &local_n0, &local_0_start);
  data = fftw_alloc_complex(alloc_local);

  /* create plan for in-place forward DFT */
  plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                              FFTW_FORWARD, FFTW_ESTIMATE);

  /* initialize data to some function my_function(x,y) */
  // for (i = 0; i < local_n0; ++i) for (j = 0; j < N1; ++j)
  //     {}

  /* compute transforms, in-place, as many times as desired */
  fftw_execute(plan);


  // One good way to output a large multi-dimensional distributed array in MPI
  // to a portable binary file is to use the free HDF5 library; see the HDF home
  // page.
  fftw_free(data);
  fftw_destroy_plan(plan);

  std::string msg("MPI unditentified");
  return msg;
}

double my_function_d(const int i, const int j)
{
  return i/j;
}

// void fftw_org_real_data_MPI_transform(int argc, char **argv)
std::string fftw_org_real_data_MPI_transform()
{
  const ptrdiff_t L = 6, M = 9;
  fftw_plan plan;
  double *data;
  ptrdiff_t alloc_local, local_n0, local_0_start;
  int i, j;

  fftw_mpi_init();
  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_2d(L, M, MPI_COMM_WORLD,
                                       &local_n0, &local_0_start);
  data = fftw_alloc_real(alloc_local);

  /* create plan for in-place REDFT10 x RODFT10 */
  plan = fftw_mpi_plan_r2r_2d(L, M, data, data, MPI_COMM_WORLD,
                              FFTW_REDFT10, FFTW_RODFT10, FFTW_MEASURE);

  /* initialize data to some function my_function(x,y) */
  for (i = 0; i < local_n0; ++i)
    for (j = 0; j < M; ++j)
      data[i*M + j] = my_function_d((int)local_0_start + i, j);

  /* compute transforms, in-place, as many times as desired */
  fftw_execute(plan);
  fftw_free(data);
  fftw_destroy_plan(plan);

  std::string msg("MPI r2r_2d");
  return msg;
}


double forcage_spectral::get_compteur0()
{
  return compteur[0];
}
double forcage_spectral::get_compteur1()
{
  return compteur[1];
}
double forcage_spectral::get_compteur2()
{
  return compteur[2];
}
