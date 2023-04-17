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

#include <Init_spectral.h>
#include <string>

Implemente_instanciable( Init_spectral, "Init_spectral", Objet_U ) ;

Sortie& Init_spectral::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& Init_spectral::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

double randInRange(double min, double max)
{
  return min + ((double)(rand())/(double)(RAND_MAX)) * (max-min);
}

void compute_inital_velocity_spectral(FixedVector<IJK_Field_double, 3>& velocity)
{
  cout << endl;
  cout << " = Initialisation du champ de vitesse =" << endl;
  cout << endl;

  fftw_plan plan_ux;
  fftw_plan plan_uy;
  fftw_plan plan_uz;

  unsigned int seed = 1;
  int kmax = velocity[0].ni();
  // double kpic = 10;
  // double s = 2;
  // DoubleTab coord_velocity_fourier(kmax,2);
  DoubleTab angles(kmax,3); // (theta_1, theta_2 et phi)

  // Tirage aleatoire
  srand(seed);
  for(int i=0; i<kmax; i++)
    {
      angles(i,0) = randInRange(0.0, 2.*M_PI);
      angles(i,1) = randInRange(0.0, 2.*M_PI);
      angles(i,2) = randInRange(0.0, 2.*M_PI);
    }

  // On définit le tableau des nombres d'ondes de 0 à Kmax/2-1 puis de -Kmax/2 à -1
  DoubleTab kx(kmax, 1);
  for (int i=0 ; i<kmax/2; i++)
    {
      kx[i] = i;
    }
  for (int i=-kmax/2; i<0; i++)
    {
      kx[i+kmax] = i;
    }
  kx[0] = 1; // pour la division par zéro


  fftw_complex *ux_f = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex) * kmax*kmax*kmax );
  fftw_complex *uy_f = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex) * kmax*kmax*kmax );
  fftw_complex *uz_f = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex) * kmax*kmax*kmax );
  // On va de -16 à 16 (on ne passe pas par zéro)
  // Les données sont stockés de 1 à 16 puis de -16 à -1
  for(int i = 0; i<kmax; i++)
    {
      for(int j = 0; j<kmax; j++)
        {
          // double Sk1k2 = sqrt(kx[i]*kx[i] + kx[j]*kx[j]);
          for(int k = 0; k<kmax; k++)
            {
              // double normK = sqrt(kx[i]*kx[i] + kx[j]*kx[j] + kx[k]*kx[k]);
              // double Ek = 1/kpic * pow(normK/kpic, s) * exp(-s/2. * pow(normK/kpic, 2.));
              // double sqrt_Ek_4piK2 = sqrt(Ek/(4.0*M_PI*normK*normK));
              // const double theta1_k =  angles((int)floor(normK), 0);
              // const double theta2_k =  angles((int)floor(normK), 1);
              // const double phi_k =  angles((int)floor(normK), 2);
              // double cTh1 = cos(theta1_k);
              // double sTh1 = sin(theta1_k);
              // double cTh2 = cos(theta2_k);
              // double sTh2 = sin(theta2_k);
              // double cPhi = cos(phi_k);
              // double sPhi = sin(phi_k);
              // double alphaKr = sqrt_Ek_4piK2 * cTh1 * cPhi;
              // double alphaKi = sqrt_Ek_4piK2 * sTh1 * cPhi;
              // double betaKr = sqrt_Ek_4piK2 * cTh2 * sPhi;
              // double betaKi = sqrt_Ek_4piK2 * sTh2 * sPhi;
              ux_f[(i*kmax+j)*kmax+k][0] = 0; // (alphaKr*normK*j + betaKr*i*k)/(normK*Sk1k2);
              ux_f[(i*kmax+j)*kmax+k][1] = 0; // (alphaKi*normK*j + betaKi*i*k)/(normK*Sk1k2);
              uy_f[(i*kmax+j)*kmax+k][0] = 0; // (betaKr*j*k - alphaKr*normK*i)/(normK*Sk1k2);
              uy_f[(i*kmax+j)*kmax+k][1] = 0; // (betaKi*j*k - alphaKi*normK*i)/(normK*Sk1k2);
              uz_f[(i*kmax+j)*kmax+k][0] = 0; // betaKr*Sk1k2/normK;
              uz_f[(i*kmax+j)*kmax+k][1] = 0; // betaKi*Sk1k2/normK;
            }
        }
    }
  // ux_f[(0*kmax+1)*kmax+0][0] = 1;
  // ux_f[(0*kmax+1)*kmax+0][1] = 1;
  ux_f[((kmax/2+1)*kmax+kmax/2+1)*kmax+kmax/2+1][0] = 1;
  uy_f[((kmax/2+1)*kmax+kmax/2+1)*kmax+kmax/2+1][0] = 1;

  // Create output fields in real space
  fftw_complex *ux_r = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex) * kmax*kmax*kmax );
  fftw_complex *uy_r = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex) * kmax*kmax*kmax );
  fftw_complex *uz_r = (fftw_complex*) fftw_malloc ( sizeof(fftw_complex) * kmax*kmax*kmax );

  // And fill them with (inverse?) Fourier Transform 3D
  plan_ux = fftw_plan_dft_3d(kmax, kmax, kmax, ux_f, ux_r, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_uy = fftw_plan_dft_3d(kmax, kmax, kmax, uy_f, uy_r, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_uz = fftw_plan_dft_3d(kmax, kmax, kmax, uz_f, uz_r, FFTW_BACKWARD, FFTW_ESTIMATE);

  // And execute Fourier 3D
  fftw_execute(plan_ux);
  fftw_execute(plan_uy);
  fftw_execute(plan_uz);

  //for (int dir = 0; dir < 3; dir++)
  {
    //const IJK_Field_double& v = velocity[dir];
    IJK_Field_double& vx = velocity[0];
    IJK_Field_double& vy = velocity[1];
    IJK_Field_double& vz = velocity[2];

    const int ni = vx.ni();
    const int nj = vx.nj();
    const int nk = vx.nk();
    if ((ni != nj) || (nj != nk))
      Process::exit();

    for (int k = 0; k < nk; k++)
      {
        for (int j = 0; j < nj; j++)
          {
            for (int i = 0; i < ni; i++)
              {
                // On ne prend que la partie réelle [0]
                vx(i,j,k) =ux_r[(i*kmax +j)*kmax+k][0]/sqrt(2);
                vy(i,j,k) =uy_r[(i*kmax +j)*kmax+k][0]/sqrt(2);
                vz(i,j,k) =uz_r[(i*kmax +j)*kmax+k][0]/sqrt(2);
              }
          }
      }

    // Vérification
    // double *Ek_x = new double [kmax];
    // double *Ek_y = new double [kmax];
    // double *Ek_z = new double [kmax];
    // double *nn = new double [kmax];
    // std::fill_n(Ek_x, kmax, static_cast<double>(0));
    // std::fill_n(Ek_y, kmax, static_cast<double>(0));
    // std::fill_n(Ek_z, kmax, static_cast<double>(0));
    // std::fill_n(nk, kmax, static_cast<double>(0));
    // for (int i = 0; i < nk; i++)
    //   {
    //     for (int j = 0; j < nj; j++)
    //       {
    //         for (int k = 0; k < ni; k++)
    //           {
    //             // Test: calcul des spectres
    //             int k2 = (int)(i*i+j*j+k*k);
    //             int indice = (int)(sqrt(k2));
    //             if (indice <= kmax)
    //               {
    //                 Ek_x[indice] += 4.*M_PI*k2*(pow(ux_f[((i*kmax)+j)*kmax+k][0], 2.)+pow(ux_f[((i*kmax)+j)*kmax+k][1], 2.));
    //                 Ek_y[indice] += 4.*M_PI*k2*(pow(uy_f[((i*kmax)+j)*kmax+k][0], 2.)+pow(uy_f[((i*kmax)+j)*kmax+k][1], 2.));
    //                 Ek_z[indice] += 4.*M_PI*k2*(pow(uz_f[((i*kmax)+j)*kmax+k][0], 2.)+pow(uz_f[((i*kmax)+j)*kmax+k][1], 2.));
    //                 nn[indice] += 1;
    //               }
    //           }
    //       }
    //   }

    // Écriture vérif (ne marche pas :( )
    // std::string nomfic("spectres.dat");
    // std::ofstream sortie(nomfic.c_str(), ios::out | ios::trunc);
    // for (int i = 0; i < kmax; i++)
    //   {
    //     Ek_x[i] /= nn[i];
    //     Ek_y[i] /= nn[i];
    //     Ek_z[i] /= nn[i];
    //     sortie << i << ", " << Ek_x[i] << ", " << Ek_y[i] << ", " << Ek_z[i] << std::endl;
    //   }
    // sortie.close();
  }
}
