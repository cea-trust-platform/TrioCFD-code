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
// File      : Filter_kernel.cpp
// Directory : $NEW_ALGO_QC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <Filter_kernel.h>
#include <IJK_Field.h>

FixedVector<double, 21> Filter_kernel_box::inhomogeneous(const bool elem,
                                                         const int k,
                                                         const int kg,
                                                         const int nktot,
                                                         const double delta,
                                                         const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1./2.); //
      filter_kernel[11] = (1./2.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./2.);
      filter_kernel[10] = (1./2.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./3.);
      filter_kernel[10] = (1./3.); //
      filter_kernel[11] = (1./3.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_box::uniform(double delta,
                                                   double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = 0.;
  filter_kernel[9] = (1./3.);
  filter_kernel[10] = (1./3.); //
  filter_kernel[11] = (1./3.);
  filter_kernel[12] = 0.;
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_13_13_base::inhomogeneous(const bool elem,
                                                                       const int k,
                                                                       const int kg,
                                                                       const int nktot,
                                                                       const double delta,
                                                                       const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==1)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (2./3.); //
      filter_kernel[11] = (1./3.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-2))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./3.);
      filter_kernel[10] = (2./3.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./3.);
      filter_kernel[10] = (1./3.); //
      filter_kernel[11] = (1./3.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_13_13_base::uniform(double delta,
                                                                 double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = 0.;
  filter_kernel[9] = (1./3.);
  filter_kernel[10] = (1./3.); //
  filter_kernel[11] = (1./3.);
  filter_kernel[12] = 0.;
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_12_14_base::inhomogeneous(const bool elem,
                                                                       const int k,
                                                                       const int kg,
                                                                       const int nktot,
                                                                       const double delta,
                                                                       const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (3./4.); //
      filter_kernel[11] = (1./4.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-2))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./4.);
      filter_kernel[10] = (3./4.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./4.);
      filter_kernel[10] = (1./2.); //
      filter_kernel[11] = (1./4.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_12_14_base::uniform(double delta,
                                                                 double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = 0.;
  filter_kernel[9] = (1./4.);
  filter_kernel[10] = (1./2.); //
  filter_kernel[11] = (1./4.);
  filter_kernel[12] = 0.;
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_23_16_base::inhomogeneous(const bool elem,
                                                                       const int k,
                                                                       const int kg,
                                                                       const int nktot,
                                                                       const double delta,
                                                                       const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==1)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (5./6.); //
      filter_kernel[11] = (1./6.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-2))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./6.);
      filter_kernel[10] = (5./6.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (1./6.);
      filter_kernel[10] = (2./3.); //
      filter_kernel[11] = (1./6.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_23_16_base::uniform(double delta,
                                                                 double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = 0.;
  filter_kernel[9] = (1./6.);
  filter_kernel[10] = (2./3.); //
  filter_kernel[11] = (1./6.);
  filter_kernel[12] = 0.;
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_14_14_18_base::inhomogeneous(const bool elem,
                                                                          const int k,
                                                                          const int kg,
                                                                          const int nktot,
                                                                          const double delta,
                                                                          const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==1)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1./2.); //
      filter_kernel[11] = (3./8.);
      filter_kernel[12] = (1./8.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==2)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (3./8.);
      filter_kernel[10] = (1./4.); //
      filter_kernel[11] = (1./4.);
      filter_kernel[12] = (1./8.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-3))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (1./8.);
      filter_kernel[9] = (1./4.);
      filter_kernel[10] = (1./4.); //
      filter_kernel[11] = (3./8.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-2))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (1./8.);
      filter_kernel[9] = (3./8.);
      filter_kernel[10] = (1./2.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (1./8.);
      filter_kernel[9] = (1./4.);
      filter_kernel[10] = (1./4.); //
      filter_kernel[11] = (1./4.);
      filter_kernel[12] = (1./8.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_14_14_18_base::uniform(double delta,
                                                                    double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = (1./8.);
  filter_kernel[9] = (1./4.);
  filter_kernel[10] = (1./4.); //
  filter_kernel[11] = (1./4.);
  filter_kernel[12] = (1./8.);
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_15_15_15_base::inhomogeneous(const bool elem,
                                                                          const int k,
                                                                          const int kg,
                                                                          const int nktot,
                                                                          const double delta,
                                                                          const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==1)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (3./5.); //
      filter_kernel[11] = (2./5.);
      filter_kernel[12] = (1./5.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==2)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (2./5.);
      filter_kernel[10] = (1./5.); //
      filter_kernel[11] = (1./5.);
      filter_kernel[12] = (1./5.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-3))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (1./5.);
      filter_kernel[9] = (1./5.);
      filter_kernel[10] = (1./5.); //
      filter_kernel[11] = (2./5.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-2))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (1./5.);
      filter_kernel[9] = (1./5.);
      filter_kernel[10] = (3./5.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = (1./5.);
      filter_kernel[9] = (1./5.);
      filter_kernel[10] = (1./5.); //
      filter_kernel[11] = (1./5.);
      filter_kernel[12] = (1./5.);
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_15_15_15_base::uniform(double delta,
                                                                    double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = (1./5.);
  filter_kernel[9] = (1./5.);
  filter_kernel[10] = (1./5.); //
  filter_kernel[11] = (1./5.);
  filter_kernel[12] = (1./5.);
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_16_16_16_112_base::inhomogeneous(const bool elem,
                                                                              const int k,
                                                                              const int kg,
                                                                              const int nktot,
                                                                              const double delta,
                                                                              const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (0.);
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = (0.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==1)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (0.);
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (3./12.); //
      filter_kernel[11] = (5./12.);
      filter_kernel[12] = (3./12.);
      filter_kernel[13] = (1./12.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==2)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (0.);
      filter_kernel[8] = (0.);
      filter_kernel[9] = (5./12.);
      filter_kernel[10] = (1./6.); //
      filter_kernel[11] = (1./6.);
      filter_kernel[12] = (1./6.);
      filter_kernel[13] = (1./12.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==3)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (0.);
      filter_kernel[8] = (3./12.);
      filter_kernel[9] = (1./6.);
      filter_kernel[10] = (1./6.); //
      filter_kernel[11] = (1./6.);
      filter_kernel[12] = (1./6.);
      filter_kernel[13] = (1./12.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-4))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (1./12.);
      filter_kernel[8] = (1./6.);
      filter_kernel[9] = (1./6.);
      filter_kernel[10] = (1./6.); //
      filter_kernel[11] = (1./6.);
      filter_kernel[12] = (3./12.);
      filter_kernel[13] = (0.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-3))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (1./12.);
      filter_kernel[8] = (1./6.);
      filter_kernel[9] = (1./6.);
      filter_kernel[10] = (1./6.); //
      filter_kernel[11] = (5./12.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = (0.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-2))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (1./12.);
      filter_kernel[8] = (3./12.);
      filter_kernel[9] = (5./12.);
      filter_kernel[10] = (3./12.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = (0.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (0.);
      filter_kernel[8] = (0.);
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = (0.);
      filter_kernel[13] = (0.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = (1./12.);
      filter_kernel[8] = (1./6.);
      filter_kernel[9] = (1./6.);
      filter_kernel[10] = (1./6.); //
      filter_kernel[11] = (1./6.);
      filter_kernel[12] = (1./6.);
      filter_kernel[13] = (1./12.);
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_16_16_16_112_base::uniform(double delta,
                                                                        double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = (1./12.);
  filter_kernel[8] = (1./6.);
  filter_kernel[9] = (1./6.);
  filter_kernel[10] = (1./6.); //
  filter_kernel[11] = (1./6.);
  filter_kernel[12] = (1./6.);
  filter_kernel[13] = (1./12.);
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_14_38_base::inhomogeneous(const bool elem,
                                                                       const int k,
                                                                       const int kg,
                                                                       const int nktot,
                                                                       const double delta,
                                                                       const ArrOfDouble_with_ghost& delta_z)
{
  FixedVector<double, 21> filter_kernel;
  if (kg==0)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==1)
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (5./8.); //
      filter_kernel[11] = (3./8.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-2))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (3./8.);
      filter_kernel[10] = (5./8.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else if (kg==(nktot-1))
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (0.);
      filter_kernel[10] = (1.); //
      filter_kernel[11] = (0.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  else
    {
      filter_kernel[0] = 0.;
      filter_kernel[1] = 0.;
      filter_kernel[2] = 0.;
      filter_kernel[3] = 0.;
      filter_kernel[4] = 0.;
      filter_kernel[5] = 0.;
      filter_kernel[6] = 0.;
      filter_kernel[7] = 0.;
      filter_kernel[8] = 0.;
      filter_kernel[9] = (3./8.);
      filter_kernel[10] = (1./4.); //
      filter_kernel[11] = (3./8.);
      filter_kernel[12] = 0.;
      filter_kernel[13] = 0.;
      filter_kernel[14] = 0.;
      filter_kernel[15] = 0.;
      filter_kernel[16] = 0.;
      filter_kernel[17] = 0.;
      filter_kernel[18] = 0.;
      filter_kernel[19] = 0.;
      filter_kernel[20] = 0.;
    }
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_weight_14_38_base::uniform(double delta,
                                                                 double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = 0.;
  filter_kernel[9] = (3./8.);
  filter_kernel[10] = (1./4.); //
  filter_kernel[11] = (3./8.);
  filter_kernel[12] = 0.;
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_laplacian::inhomogeneous(const bool elem,
                                                               const int k,
                                                               const int kg,
                                                               const int nktot,
                                                               const double delta,
                                                               const ArrOfDouble_with_ghost& delta_z)
{
  double dz;
  double dm;
  double dp;

  if (elem)
    {
      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_m1_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z[k-1];
      const double dz_p1_glo = (kg+1<0 || kg+1>(nktot-1)) ? 0. : delta_z[k+1];
      const double delta_m_glo = kg==0 ? 0.5*dz_glo : 0.5*(dz_glo + dz_m1_glo);
      const double delta_p_glo = kg==(nktot-1) ? 0.5*dz_glo : 0.5*(dz_glo + dz_p1_glo);

      dz = dz_glo;
      dm = delta_m_glo;
      dp = delta_p_glo;
    }
  else
    {
      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_m1_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z[k-1];
      const double delta_m_glo = kg==0 ? 0.5*dz_glo : 0.5*(dz_glo + dz_m1_glo);

      dz = delta_m_glo;
      dm = dz_m1_glo;
      dp = dz_glo;
    }

  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = 0.;
  filter_kernel[9] = delta==0. ? 0. : (delta*delta*dm)/(24.*dz*dp*dm);
  filter_kernel[10] = delta==0. ? 0. : 1 - (delta*delta*(dp+dm))/(24.*dz*dp*dm); //
  filter_kernel[11] = delta==0. ? 0. : (delta*delta*dp)/(24.*dz*dp*dm);
  filter_kernel[12] = 0.;
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

FixedVector<double, 21> Filter_kernel_laplacian::uniform(double delta,
                                                         double dz)
{
  FixedVector<double, 21> filter_kernel;
  filter_kernel[0] = 0.;
  filter_kernel[1] = 0.;
  filter_kernel[2] = 0.;
  filter_kernel[3] = 0.;
  filter_kernel[4] = 0.;
  filter_kernel[5] = 0.;
  filter_kernel[6] = 0.;
  filter_kernel[7] = 0.;
  filter_kernel[8] = 0.;
  filter_kernel[9] = delta==0. ? 0. : (delta*delta)/(24.*dz*dz);
  filter_kernel[10] = delta==0. ? 0. : 1 - (delta*delta*2)/(24.*dz*dz); //
  filter_kernel[11] = delta==0. ? 0. : (delta*delta)/(24.*dz*dz);
  filter_kernel[12] = 0.;
  filter_kernel[13] = 0.;
  filter_kernel[14] = 0.;
  filter_kernel[15] = 0.;
  filter_kernel[16] = 0.;
  filter_kernel[17] = 0.;
  filter_kernel[18] = 0.;
  filter_kernel[19] = 0.;
  filter_kernel[20] = 0.;
  return filter_kernel;
}

