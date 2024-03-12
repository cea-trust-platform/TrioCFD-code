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
// File      : Turbulent_viscosity.cpp
// Directory : $NEW_ALGO_QC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <Turbulent_viscosity.h>
#include <DNS_QC.h>

static inline double calculer_somme_mimi(double m_i,
                                         double m_j,
                                         double m_k)
{
  return   m_i * m_i
           + m_j * m_j
           + m_k * m_k;
}

static inline double calculer_somme_mijmij(double m_ii,
                                           double m_ij,
                                           double m_ik,
                                           double m_ji,
                                           double m_jj,
                                           double m_jk,
                                           double m_ki,
                                           double m_kj,
                                           double m_kk)
{
  return   m_ii * m_ii
           + m_ij * m_ij
           + m_ik * m_ik
           + m_ji * m_ji
           + m_jj * m_jj
           + m_jk * m_jk
           + m_ki * m_ki
           + m_kj * m_kj
           + m_kk * m_kk;
}

static inline double calculer_somme_mijmji(double m_ii,
                                           double m_ij,
                                           double m_ik,
                                           double m_ji,
                                           double m_jj,
                                           double m_jk,
                                           double m_ki,
                                           double m_kj,
                                           double m_kk)
{
  return   m_ii * m_ii
           + m_ij * m_ji
           + m_ik * m_ki
           + m_ji * m_ij
           + m_jj * m_jj
           + m_jk * m_kj
           + m_ki * m_ik
           + m_kj * m_jk
           + m_kk * m_kk;
}

double Turbulent_viscosity_constant::operator()(double constant,
                                                double dx,
                                                double dy,
                                                double dz,
                                                double rho,
                                                FixedVector<FixedVector<double, 3>, 3>& g,
                                                FixedVector<double, 3>& q)
{
  return constant;
}

double Turbulent_viscosity_unsrho::operator()(double constant,
                                              double dx,
                                              double dy,
                                              double dz,
                                              double rho,
                                              FixedVector<FixedVector<double, 3>, 3>& g,
                                              FixedVector<double, 3>& q)
{
  return constant/rho;
}

double Turbulent_viscosity_smagorinsky::operator()(double smagorinsky_constant,
                                                   double dx,
                                                   double dy,
                                                   double dz,
                                                   double rho,
                                                   FixedVector<FixedVector<double, 3>, 3>& g,
                                                   FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele de smagorinsky
  // ----------------
  // tau_ij = -2 (c_s delta)^2 |S| S_ij
  //   |S| = sqrt(2sum_SijSij) = sqrt(sum_gijgij + sum_gijgji)
  const double somme_gijgij = calculer_somme_mijmij(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double somme_gijgji = calculer_somme_mijmji(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double norme_s = sqrt(somme_gijgij + somme_gijgji);

  const double filter_volume = dx * dy * dz;
  const double filter_lenght = cbrt(filter_volume);

  return smagorinsky_constant*smagorinsky_constant * filter_lenght*filter_lenght * norme_s;
}

double Turbulent_viscosity_vreman::operator()(double vreman_constant,
                                              double dx,
                                              double dy,
                                              double dz,
                                              double rho,
                                              FixedVector<FixedVector<double, 3>, 3>& g,
                                              FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele de vreman
  // ----------------
  // tau_ij = -2 nu_sgs S_ij
  //   nu_sgs = C_v PI^g
  //   PI^g = sqrt(B_b^g/(a_kl a_kl))
  //   B_b^g = b_11^g b_22^g - b_12^g b_12^g + b_11^g b_33^g - b_13^g b_13^g + b_22^g b_33^g -b_23^g b_23^g
  //   b_ij^g = delta_m^2 a_mi a_mj
  //   a_ij = duj/dxi = gji
  // New variables to use Vreman's notations
  const double a_ii = g_ii;
  const double a_ij = g_ji;
  const double a_ik = g_ki;
  const double a_ji = g_ij;
  const double a_jj = g_jj;
  const double a_jk = g_kj;
  const double a_ki = g_ik;
  const double a_kj = g_jk;
  const double a_kk = g_kk;

  const double b_ii = dx*dx*a_ii*a_ii + dy*dy*a_ji*a_ji + dz*dz*a_ki*a_ki;
  const double b_ij = dx*dx*a_ii*a_ij + dy*dy*a_ji*a_jj + dz*dz*a_ki*a_kj;
  const double b_ik = dx*dx*a_ii*a_ik + dy*dy*a_ji*a_jk + dz*dz*a_ki*a_kk;
  const double b_jj = dx*dx*a_ij*a_ij + dy*dy*a_jj*a_jj + dz*dz*a_kj*a_kj;
  const double b_jk = dx*dx*a_ij*a_ik + dy*dy*a_jj*a_jk + dz*dz*a_kj*a_kk;
  const double b_kk = dx*dx*a_ik*a_ik + dy*dy*a_jk*a_jk + dz*dz*a_kk*a_kk;

  const double b_b =  b_ii*b_jj - b_ij*b_ij
                      + b_ii*b_kk - b_ik*b_ik
                      + b_jj*b_kk - b_jk*b_jk;

  const double somme_aijaij = calculer_somme_mijmji(a_ii, a_ij, a_ik, a_ji, a_jj, a_jk, a_ki, a_kj, a_kk);

  double vreman_pi;
  if (somme_aijaij == 0.)
    {
      vreman_pi = 0;
    }
  else
    {
      if (b_b/(somme_aijaij) < 0)
        {
          vreman_pi = 0;
          //Cerr << "Warning: In the vreman's model,"
          //     << " the value of b_b/(somme_aijaij)=" << b_b/(somme_aijaij) << " is negative, which should not be possible."
          //     << " If the value is very close to zero, this might also be a numerical error." << finl;
        }
      else
        {
          vreman_pi = (somme_aijaij == 0) ? 0. : sqrt(b_b/(somme_aijaij));
        }
    }
  return vreman_constant * vreman_pi;
}

double Turbulent_viscosity_sigma::operator()(double sigma_constant,
                                             double dx,
                                             double dy,
                                             double dz,
                                             double rho,
                                             FixedVector<FixedVector<double, 3>, 3>& g,
                                             FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele sigma
  // ----------------
  // tau_ij = -2 (c_sigma delta)^2 d_sigma S_ij
  //   d_sigma = sigma_3 (sigma_1 - sigma_2) (sigma_2 - sigma_3) / sigma_1^2
  //   sigma_[123] : the three singular value of the matrix g
  //                 ordered such that sigma_1 >= sigma_2 >= sigma_3
  // method b of: nicoud, baya toda, nicoud, cabrit, bose, lee (2013) sigma
  //   1.   G = g^t g      ... gtg_ij = sum_m g_mi g_mj
  //   2.1. I1 = tr(G)     ... i1 = gtg_ii + gtg_jj + gtg_kk
  //   2.2. I2 = 0.5 (tr(G)^2 - tr(G^2))
  //                       ... gtg_carre_ij = sum_m gtg_im gtg_mj
  //                       ... i2 = &c.
  //   2.3. I3 = det(G)    ... i3 =  gtg_ii * gtg_jj * gtg_kk
  //                               - gtg_ii * gtg_kj * gtg_jk
  //                               - gtg_ji * gtg_ij * gtg_kk
  //                               + gtg_ji * gtg_kj * gtg_ik
  //                               + gtg_ki * gtg_ij * gtg_jk
  //                               - gtg_ki * gtg_jj * gtg_ik
  //   3.1. alpha_1 = i1*i1/9. - i2/3.
  //   3.2. alpha_2 = i1*i1*i1/27. - i1*i2/6. + i3/2.
  //   3.3. alpha_3 = 1/3. * arccos(alpha_2/alpha_1^(3./2.))
  //   4.1. sigma_1 = sqrt(i1/3. + 2sqrt(alpha_1)cos(alpha_3))
  //   4.2. sigma_2 = sqrt(i1/3. - 2sqrt(alpha_1)cos(pi/3 + alpha_3))
  //   4.3. sigma_3 = sqrt(i1/3. - 2sqrt(alpha_1)cos(pi/3 - alpha_3))
  const double gtg_ii = g_ii*g_ii + g_ji*g_ji + g_ki*g_ki;
  const double gtg_ij = g_ii*g_ij + g_ji*g_jj + g_ki*g_kj;
  const double gtg_ik = g_ii*g_ik + g_ji*g_jk + g_ki*g_kk;
  const double gtg_ji = g_ij*g_ii + g_jj*g_ji + g_kj*g_ki;
  const double gtg_jj = g_ij*g_ij + g_jj*g_jj + g_kj*g_kj;
  const double gtg_jk = g_ij*g_ik + g_jj*g_jk + g_kj*g_kk;
  const double gtg_ki = g_ik*g_ii + g_jk*g_ji + g_kk*g_ki;
  const double gtg_kj = g_ik*g_ij + g_jk*g_jj + g_kk*g_kj;
  const double gtg_kk = g_ik*g_ik + g_jk*g_jk + g_kk*g_kk;

  const double gtg_carre_ii = gtg_ii*gtg_ii + gtg_ij*gtg_ji + gtg_ik*gtg_ki;
  const double gtg_carre_jj = gtg_ji*gtg_ij + gtg_jj*gtg_jj + gtg_jk*gtg_kj;
  const double gtg_carre_kk = gtg_ki*gtg_ik + gtg_kj*gtg_jk + gtg_kk*gtg_kk;

  const double trace_gtg = gtg_ii + gtg_jj + gtg_kk;
  const double trace_gtg_carre = gtg_carre_ii + gtg_carre_jj + gtg_carre_kk;
  const double det_gtg =  gtg_ii * gtg_jj * gtg_kk
                          - gtg_ii * gtg_kj * gtg_jk
                          - gtg_ji * gtg_ij * gtg_kk
                          + gtg_ji * gtg_kj * gtg_ik
                          + gtg_ki * gtg_ij * gtg_jk
                          - gtg_ki * gtg_jj * gtg_ik;

  const double i1 = trace_gtg;
  const double i2 = 0.5*(trace_gtg*trace_gtg - trace_gtg_carre);
  const double i3 = det_gtg;

  const double alpha_1 = i1*i1/9. - i2/3.;
  const double alpha_2 = i1*i1*i1/27. - i1*i2/6. + i3/2;
  double alpha_3;
  if (alpha_2/(sqrt(alpha_1*alpha_1*alpha_1)) > 1)
    {
      alpha_3 = 0;
      //Cerr << "Warning: In the sigma model,"
      //     << " the value of alpha_2/(sqrt(alpha_1*alpha_1*alpha_1))=" << alpha_2/(sqrt(alpha_1*alpha_1*alpha_1)) << " is greater than one, which should not be possible."
      //     << " If the value is very close to one, this might also be a numerical error." << finl;
    }
  else if (alpha_2/(sqrt(alpha_1*alpha_1*alpha_1)) < -1)
    {
      alpha_3 = M_PI;
      //Cerr << "Warning: In the sigma model,"
      //     << " the value of alpha_2/(sqrt(alpha_1*alpha_1*alpha_1))=" << alpha_2/(sqrt(alpha_1*alpha_1*alpha_1)) << " is less than minus one, which should not be possible."
      //     << " If the value is very close to minus one, this might also be a numerical error." << finl;
    }
  else
    {
      alpha_3 = (1./3.)*acos(alpha_2/(sqrt(alpha_1*alpha_1*alpha_1)));
    }

  double sigma_1;
  if (i1/3. + 2.*sqrt(alpha_1)*cos(alpha_3) < 0)
    {
      sigma_1 = 0;
      //Cerr << "Warning: In the sigma model,"
      //     << " the value of i1/3. + 2.*sqrt(alpha_1)*cos(alpha_3)=" << i1/3. + 2.*sqrt(alpha_1)*cos(alpha_3) << " is negative, which should not be possible."
      //     << " If the value is very close to zero, this might also be a numerical error." << finl;
    }
  else
    {
      sigma_1 = sqrt(i1/3. + 2.*sqrt(alpha_1)*cos(alpha_3));
    }

  double sigma_2;
  if (i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. + alpha_3) < 0)
    {
      sigma_2 = 0;
      //Cerr << "Warning: In the sigma model,"
      //     << " the value of i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. + alpha_3)=" << i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. + alpha_3) << " is negative, which should not be possible."
      //     << " If the value is very close to zero, this might also be a numerical error." << finl;
    }
  else
    {
      sigma_2 = sqrt(i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. + alpha_3));
    }

  double sigma_3;
  if (i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. - alpha_3) < 0)
    {
      sigma_3 = 0;
      //Cerr << "Warning: In the sigma model,"
      //     << " the value of i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. - alpha_3)=" << i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. - alpha_3) << " is negative, which should not be possible."
      //     << " If the value is very close to zero, this might also be a numerical error." << finl;
    }
  else
    {
      sigma_3 = sqrt(i1/3. - 2.*sqrt(alpha_1)*cos(M_PI/3. - alpha_3));
    }

  const double d_sigma = sigma_3*(sigma_1 - sigma_2)*(sigma_2 - sigma_3)/(sigma_1*sigma_1);

  const double filter_volume = dx * dy * dz;
  const double filter_lenght = cbrt(filter_volume);

  return sigma_constant*sigma_constant * filter_lenght*filter_lenght * d_sigma;
}

double Turbulent_viscosity_wale::operator()(double wale_constant,
                                            double dx,
                                            double dy,
                                            double dz,
                                            double rho,
                                            FixedVector<FixedVector<double, 3>, 3>& g,
                                            FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele wale
  // ----------------
  // tau_ij = -2 (c_wale delta)^2 f_wale S_ij
  //
  //                      (SSdij SSdij)^(3./2.)
  //   f_wale = -----------------------------------------
  //            (Sij Sij)^(5./2.) + (SSdij SSdij)^(5./4.)
  //
  //   Sij = 0.5 * (gij + gji)
  //   SSdij = 0.5 * (gimgmj + gjmgmi) - delta_ij (1./3.) gkmgmk
  //                    ^        ^
  //               gij_carre   gji_carre
  const double g_carre_ii = g_ii*g_ii + g_ij*g_ji + g_ik*g_ki;
  const double g_carre_ij = g_ii*g_ij + g_ij*g_jj + g_ik*g_kj;
  const double g_carre_ik = g_ii*g_ik + g_ij*g_jk + g_ik*g_kk;
  const double g_carre_ji = g_ji*g_ii + g_jj*g_ji + g_jk*g_ki;
  const double g_carre_jj = g_ji*g_ij + g_jj*g_jj + g_jk*g_kj;
  const double g_carre_jk = g_ji*g_ik + g_jj*g_jk + g_jk*g_kk;
  const double g_carre_ki = g_ki*g_ii + g_kj*g_ji + g_kk*g_ki;
  const double g_carre_kj = g_ki*g_ij + g_kj*g_jj + g_kk*g_kj;
  const double g_carre_kk = g_ki*g_ik + g_kj*g_jk + g_kk*g_kk;

  const double trace_g_carre = g_carre_ii + g_carre_jj + g_carre_kk;

  const double ssd_ii = 0.5 * (g_carre_ii + g_carre_ii) - (1./3.)*trace_g_carre;
  const double ssd_ij = 0.5 * (g_carre_ij + g_carre_ji);
  const double ssd_ik = 0.5 * (g_carre_ik + g_carre_ki);
  const double ssd_ji = 0.5 * (g_carre_ji + g_carre_ij);
  const double ssd_jj = 0.5 * (g_carre_jj + g_carre_jj) - (1./3.)*trace_g_carre;
  const double ssd_jk = 0.5 * (g_carre_jk + g_carre_kj);
  const double ssd_ki = 0.5 * (g_carre_ki + g_carre_ik);
  const double ssd_kj = 0.5 * (g_carre_kj + g_carre_jk);
  const double ssd_kk = 0.5 * (g_carre_kk + g_carre_kk) - (1./3.)*trace_g_carre;

  const double somme_ssdijssdij = calculer_somme_mijmji(ssd_ii, ssd_ij, ssd_ik, ssd_ji, ssd_jj, ssd_jk, ssd_ki, ssd_kj, ssd_kk);

  const double s_ii = 0.5 * (g_ii + g_ii);
  const double s_ij = 0.5 * (g_ij + g_ji);
  const double s_ik = 0.5 * (g_ik + g_ki);
  const double s_ji = 0.5 * (g_ji + g_ij);
  const double s_jj = 0.5 * (g_jj + g_jj);
  const double s_jk = 0.5 * (g_jk + g_kj);
  const double s_ki = 0.5 * (g_ki + g_ik);
  const double s_kj = 0.5 * (g_kj + g_jk);
  const double s_kk = 0.5 * (g_kk + g_kk);

  const double somme_sijsij = calculer_somme_mijmji(s_ii, s_ij, s_ik, s_ji, s_jj, s_jk, s_ki, s_kj, s_kk);

  const double f_wale = ((pow(somme_sijsij,5./2.)+pow(somme_ssdijssdij,5./4.))==0) ? 0. : pow(somme_ssdijssdij,3./2.)/(pow(somme_sijsij,5./2.)+pow(somme_ssdijssdij,5./4.));

  const double filter_volume = dx * dy * dz;
  const double filter_lenght = cbrt(filter_volume);

  return wale_constant*wale_constant * filter_lenght*filter_lenght * f_wale;
}

double Turbulent_viscosity_amd::operator()(double amd_constant,
                                           double dx,
                                           double dy,
                                           double dz,
                                           double rho,
                                           FixedVector<FixedVector<double, 3>, 3>& g,
                                           FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele amd
  // ----------------
  // tau_ij = -2 c_amd f_amd S_ij
  //   f_amd = max(-cijsij,0)/gklgkl
  //   c_ij = delta_m^2 g_im g_jm
  const double s_ii = 0.5 * (g_ii + g_ii);
  const double s_ij = 0.5 * (g_ij + g_ji);
  const double s_ik = 0.5 * (g_ik + g_ki);
  const double s_ji = 0.5 * (g_ji + g_ij);
  const double s_jj = 0.5 * (g_jj + g_jj);
  const double s_jk = 0.5 * (g_jk + g_kj);
  const double s_ki = 0.5 * (g_ki + g_ik);
  const double s_kj = 0.5 * (g_kj + g_jk);
  const double s_kk = 0.5 * (g_kk + g_kk);

  const double c_ii = dx*dx*g_ii*g_ii + dy*dy*g_ij*g_ij + dz*dz*g_ik*g_ik;
  const double c_ij = dx*dx*g_ii*g_ji + dy*dy*g_ij*g_jj + dz*dz*g_ik*g_jk;
  const double c_ik = dx*dx*g_ii*g_ki + dy*dy*g_ij*g_kj + dz*dz*g_ik*g_kk;
  const double c_ji = dx*dx*g_ji*g_ii + dy*dy*g_jj*g_ij + dz*dz*g_jk*g_ik;
  const double c_jj = dx*dx*g_ji*g_ji + dy*dy*g_jj*g_jj + dz*dz*g_jk*g_jk;
  const double c_jk = dx*dx*g_ji*g_ki + dy*dy*g_jj*g_kj + dz*dz*g_jk*g_kk;
  const double c_ki = dx*dx*g_ki*g_ii + dy*dy*g_kj*g_ij + dz*dz*g_kk*g_ik;
  const double c_kj = dx*dx*g_ki*g_ji + dy*dy*g_kj*g_jj + dz*dz*g_kk*g_jk;
  const double c_kk = dx*dx*g_ki*g_ki + dy*dy*g_kj*g_kj + dz*dz*g_kk*g_kk;

  const double somme_gijgij = calculer_somme_mijmij(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double somme_cijsij = c_ii * s_ii
                              + c_ij * s_ij
                              + c_ik * s_ik
                              + c_ji * s_ji
                              + c_jj * s_jj
                              + c_jk * s_jk
                              + c_ki * s_ki
                              + c_kj * s_kj
                              + c_kk * s_kk;
  const double f_amd = (somme_gijgij==0) ? 0. : max(0., -somme_cijsij)/somme_gijgij;

  return amd_constant * f_amd;
}

double Turbulent_viscosity_amdnoclip::operator()(double amdnoclip_constant,
                                                 double dx,
                                                 double dy,
                                                 double dz,
                                                 double rho,
                                                 FixedVector<FixedVector<double, 3>, 3>& g,
                                                 FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele amdnoclip
  // ----------------
  // tau_ij = -2 c_amdnoclip f_amdnoclip S_ij
  //   f_amdnoclip = -cijsij/gklgkl
  //   c_ij = delta_m^2 g_im g_jm
  const double s_ii = 0.5 * (g_ii + g_ii);
  const double s_ij = 0.5 * (g_ij + g_ji);
  const double s_ik = 0.5 * (g_ik + g_ki);
  const double s_ji = 0.5 * (g_ji + g_ij);
  const double s_jj = 0.5 * (g_jj + g_jj);
  const double s_jk = 0.5 * (g_jk + g_kj);
  const double s_ki = 0.5 * (g_ki + g_ik);
  const double s_kj = 0.5 * (g_kj + g_jk);
  const double s_kk = 0.5 * (g_kk + g_kk);

  const double c_ii = dx*dx*g_ii*g_ii + dy*dy*g_ij*g_ij + dz*dz*g_ik*g_ik;
  const double c_ij = dx*dx*g_ii*g_ji + dy*dy*g_ij*g_jj + dz*dz*g_ik*g_jk;
  const double c_ik = dx*dx*g_ii*g_ki + dy*dy*g_ij*g_kj + dz*dz*g_ik*g_kk;
  const double c_ji = dx*dx*g_ji*g_ii + dy*dy*g_jj*g_ij + dz*dz*g_jk*g_ik;
  const double c_jj = dx*dx*g_ji*g_ji + dy*dy*g_jj*g_jj + dz*dz*g_jk*g_jk;
  const double c_jk = dx*dx*g_ji*g_ki + dy*dy*g_jj*g_kj + dz*dz*g_jk*g_kk;
  const double c_ki = dx*dx*g_ki*g_ii + dy*dy*g_kj*g_ij + dz*dz*g_kk*g_ik;
  const double c_kj = dx*dx*g_ki*g_ji + dy*dy*g_kj*g_jj + dz*dz*g_kk*g_jk;
  const double c_kk = dx*dx*g_ki*g_ki + dy*dy*g_kj*g_kj + dz*dz*g_kk*g_kk;

  const double somme_gijgij = calculer_somme_mijmij(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double somme_cijsij = c_ii * s_ii
                              + c_ij * s_ij
                              + c_ik * s_ik
                              + c_ji * s_ji
                              + c_jj * s_jj
                              + c_jk * s_jk
                              + c_ki * s_ki
                              + c_kj * s_kj
                              + c_kk * s_kk;
  const double f_amdnoclip = (somme_gijgij==0) ? 0. : -somme_cijsij/somme_gijgij;

  return amdnoclip_constant * f_amdnoclip;
}

double Turbulent_viscosity_amdscalar::operator()(double amdscalar_constant,
                                                 double dx,
                                                 double dy,
                                                 double dz,
                                                 double rho,
                                                 FixedVector<FixedVector<double, 3>, 3>& g,
                                                 FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  const double& q_i = q[0];
  const double& q_j = q[1];
  const double& q_k = q[2];

  // modele amdscalar
  // ----------------
  // tau_ij = -2 c_amdscalar f_amdscalar S_ij
  //   f_amdscalar = -ciqi/qkqk
  //   c_i = delta_m^2 g_im q_m
  const double c_i = dx*dx*g_ii*q_i + dy*dy*g_ij*q_j + dz*dz*g_ik*q_k;
  const double c_j = dx*dx*g_ji*q_i + dy*dy*g_jj*q_j + dz*dz*g_jk*q_k;
  const double c_k = dx*dx*g_ki*q_i + dy*dy*g_kj*q_j + dz*dz*g_kk*q_k;

  const double somme_qiqi = calculer_somme_mimi(q_i, q_j, q_k);
  const double somme_ciqi = c_i * q_i
                            + c_j * q_j
                            + c_k * q_k;
  const double f_amdscalar = (somme_qiqi==0) ? 0. : max(0., -somme_ciqi)/somme_qiqi;

  return amdscalar_constant * f_amdscalar;
}

double Turbulent_viscosity_amdscalarnoclip::operator()(double amdscalarnoclip_constant,
                                                       double dx,
                                                       double dy,
                                                       double dz,
                                                       double rho,
                                                       FixedVector<FixedVector<double, 3>, 3>& g,
                                                       FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  const double& q_i = q[0];
  const double& q_j = q[1];
  const double& q_k = q[2];

  // modele amdscalarnoclip
  // ----------------
  // tau_ij = -2 c_amdscalarnoclip f_amdscalarnoclip S_ij
  //   f_amdscalarnoclip = -ciqi/qkqk
  //   c_i = delta_m^2 g_im q_m
  const double c_i = dx*dx*g_ii*q_i + dy*dy*g_ij*q_j + dz*dz*g_ik*q_k;
  const double c_j = dx*dx*g_ji*q_i + dy*dy*g_jj*q_j + dz*dz*g_jk*q_k;
  const double c_k = dx*dx*g_ki*q_i + dy*dy*g_kj*q_j + dz*dz*g_kk*q_k;

  const double somme_qiqi = calculer_somme_mimi(q_i, q_j, q_k);
  const double somme_ciqi = c_i * q_i
                            + c_j * q_j
                            + c_k * q_k;
  const double f_amdscalarnoclip = (somme_qiqi==0) ? 0. : -somme_ciqi/somme_qiqi;

  return amdscalarnoclip_constant * f_amdscalarnoclip;
}

double Turbulent_viscosity_rds::operator()(double rds_constant,
                                           double dx,
                                           double dy,
                                           double dz,
                                           double rho,
                                           FixedVector<FixedVector<double, 3>, 3>& g,
                                           FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele rds (reversed dynamic structure)
  // ----------------
  // tau_ij = -1/12 c_rds (c_kk/|S|) S_ij
  //   |S| = sqrt(2sum_SijSij) = sqrt(sum_gijgij + sum_gijgji)
  //   c_ij = delta_m^2 g_im g_jm
  const double somme_gijgij = calculer_somme_mijmij(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double somme_gijgji = calculer_somme_mijmji(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double norme_s = sqrt(somme_gijgij + somme_gijgji);

  const double c_ii = dx*dx*g_ii*g_ii + dy*dy*g_ij*g_ij + dz*dz*g_ik*g_ik;
  const double c_jj = dx*dx*g_ji*g_ji + dy*dy*g_jj*g_jj + dz*dz*g_jk*g_jk;
  const double c_kk = dx*dx*g_ki*g_ki + dy*dy*g_kj*g_kj + dz*dz*g_kk*g_kk;

  const double trace_c = c_ii + c_jj + c_kk;

  const double f_rds = (norme_s==0) ? 0. : trace_c/norme_s;

  return rds_constant * (1./12.) * f_rds;
}

double Turbulent_viscosity_vss::operator()(double vss_constant,
                                           double dx,
                                           double dy,
                                           double dz,
                                           double rho,
                                           FixedVector<FixedVector<double, 3>, 3>& g,
                                           FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele vss (reversed dynamic structure)
  // ----------------
  // tau_ij = -2 (d_vss delta)^2 S_ij
  //   d_vss = c_r (Rij Rij)^(3./2.)/(Sij Sij)^{5./2.}
  //   R_ij = beta_i duj/dxj (pas de somme sur j)
  //   [beta_i] = (S_23, S_13, S_12)
  //   c_r = 1.3 (valeur recommandee)
  const double s_ii = 0.5 * (g_ii + g_ii);
  const double s_ij = 0.5 * (g_ij + g_ji);
  const double s_ik = 0.5 * (g_ik + g_ki);
  const double s_ji = 0.5 * (g_ji + g_ij);
  const double s_jj = 0.5 * (g_jj + g_jj);
  const double s_jk = 0.5 * (g_jk + g_kj);
  const double s_ki = 0.5 * (g_ki + g_ik);
  const double s_kj = 0.5 * (g_kj + g_jk);
  const double s_kk = 0.5 * (g_kk + g_kk);

  const double somme_sijsij = calculer_somme_mijmji(s_ii, s_ij, s_ik, s_ji, s_jj, s_jk, s_ki, s_kj, s_kk);

  const double r_ii = s_jk * g_ii;
  const double r_ij = s_jk * g_jj;
  const double r_ik = s_jk * g_kk;
  const double r_ji = s_ik * g_ii;
  const double r_jj = s_ik * g_jj;
  const double r_jk = s_ik * g_kk;
  const double r_ki = s_ij * g_ii;
  const double r_kj = s_ij * g_jj;
  const double r_kk = s_ij * g_kk;

  const double somme_rijrij = calculer_somme_mijmji(r_ii, r_ij, r_ik, r_ji, r_jj, r_jk, r_ki, r_kj, r_kk);

  const double d_vss = (pow(somme_sijsij,5./2.)==0) ? 0. : (pow(somme_rijrij,3./2.))/(pow(somme_sijsij,5./2.));

  const double filter_volume = dx * dy * dz;
  const double filter_lenght = cbrt(filter_volume);

  return vss_constant*vss_constant * filter_lenght*filter_lenght * d_vss;
}

double Turbulent_viscosity_kobayashi::operator()(double kobayashi_constant,
                                                 double dx,
                                                 double dy,
                                                 double dz,
                                                 double rho,
                                                 FixedVector<FixedVector<double, 3>, 3>& g,
                                                 FixedVector<double, 3>& q)
{
  const double& g_ii = g[0][0];
  const double& g_ij = g[0][1];
  const double& g_ik = g[0][2];
  const double& g_ji = g[1][0];
  const double& g_jj = g[1][1];
  const double& g_jk = g[1][2];
  const double& g_ki = g[2][0];
  const double& g_kj = g[2][1];
  const double& g_kk = g[2][2];

  // modele de kobayashi
  // ----------------
  // tau_ij = -2 c_k (delta)^2 (1-F)|F|^(3./2.) |S| S_ij
  //   |S| = sqrt(2sum_SijSij) = sqrt(sum_gijgij + sum_gijgji)
  //   F = Q/E = - gijgji/gklgkl
  //   Q = 0.5(WijWij - SijSij)
  //   E = 0.5(WijWij + SijSij)
  //   Sij = 0.5(gij + gji)
  //   Wij = 0.5(gij - gji)
  const double somme_gijgij = calculer_somme_mijmij(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double somme_gijgji = calculer_somme_mijmji(g_ii, g_ij, g_ik, g_ji, g_jj, g_jk, g_ki, g_kj, g_kk);
  const double norme_s = sqrt(somme_gijgij + somme_gijgji);

  const double fcs = - somme_gijgji/somme_gijgij;

  const double f_koba = pow(fabs(fcs),3./2.) * (1 - fcs) * norme_s;

  const double filter_volume = dx * dy * dz;
  const double filter_lenght = cbrt(filter_volume);

  return kobayashi_constant * filter_lenght*filter_lenght * f_koba;
}

