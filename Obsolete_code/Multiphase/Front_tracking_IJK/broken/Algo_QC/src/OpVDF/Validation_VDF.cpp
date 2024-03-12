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
// File      : Validation_VDF.cpp
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#if 0
#include <Validation_VDF.h>

void Validation_VDF::compare(const IJK_Field_double& ref, const IJK_Field_double& x, const char *msg)
{
  int ni = x.ni();
  int nj = x.nj();
  int nk = x.nk();
  int i0 = x.get_splitting().get_offset_local(DIRECTION_I);
  int j0 = x.get_splitting().get_offset_local(DIRECTION_J);
  int k0 = x.get_splitting().get_offset_local(DIRECTION_K);
  double maxd = 0.;
  for(int k=0; k<nk; k++)
    for(int j=0; j<nj; j++)
      for(int i=0; i<ni; i++)
        {
          double d = fabs(x(i,j,k)-ref(i,j,k));
          maxd = max(d, maxd);
          if (d > 1e-15)
            {
              char s[1000];
              snprintf(s, 1000, " local_ijk=( %3d %3d %3d )  global_ijk=( %3d %3d %3d ) ref= %10g  current= %10g  delta=%5g",
                      i, j, k, i+i0, j+j0, k+k0, ref(i,j,k), x(i,j,k), d);
              Process::Journal() << msg << s << finl;
            }
        }
  maxd = Process::mp_max(maxd);
  Process::Journal() << msg << " max_error= " << maxd << finl;
  Cerr << msg << " max_error= " << maxd << finl;
}

// Called by completer()
void Validation_VDF::init(const char * message)
{
  msg_ = message;

  // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
  const char * ijkdis_name = IJK_discretization::get_conventional_name();
  const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));
  const IJK_Splitting& split = ijkdis.get_IJK_splitting();
  vx_.allocate(split, IJK_Splitting::FACES_I, 1);
  vy_.allocate(split, IJK_Splitting::FACES_J, 1);
  vz_.allocate(split, IJK_Splitting::FACES_K, 1);
  inputx_.allocate(split, IJK_Splitting::FACES_I, 2);
  inputy_.allocate(split, IJK_Splitting::FACES_J, 2);
  inputz_.allocate(split, IJK_Splitting::FACES_K, 2);
  dvx_.allocate(split, IJK_Splitting::FACES_I, 0);
  dvy_.allocate(split, IJK_Splitting::FACES_J, 0);
  dvz_.allocate(split, IJK_Splitting::FACES_K, 0);

  Nom tmp(message);
  counter = statistiques().new_counter(0, message + "operateur ijk");
  counter2 = statistiques().new_counter(0,message + "operateur vdf");
}

// Called before computing operators
void Validation_VDF::prepare(const DoubleTab& inco, const DoubleTab& convecting_velocity) const
{
  // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
  const char * ijkdis_name = IJK_discretization::get_conventional_name();
  const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));

  const VDF_to_IJK& vdf_to_ijk_i = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_I);
  const VDF_to_IJK& vdf_to_ijk_j = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_J);
  const VDF_to_IJK& vdf_to_ijk_k = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_K);

  vdf_to_ijk_i.convert_to_ijk(inco, inputx_);
  vdf_to_ijk_j.convert_to_ijk(inco, inputy_);
  vdf_to_ijk_k.convert_to_ijk(inco, inputz_);
  vdf_to_ijk_i.convert_to_ijk(convecting_velocity, vx_);
  vdf_to_ijk_j.convert_to_ijk(convecting_velocity, vy_);
  vdf_to_ijk_k.convert_to_ijk(convecting_velocity, vz_);
  vx_.echange_espace_virtuel(vx_.ghost());
  vy_.echange_espace_virtuel(vy_.ghost());
  vz_.echange_espace_virtuel(vz_.ghost());
  inputx_.echange_espace_virtuel(inputx_.ghost());
  inputy_.echange_espace_virtuel(inputy_.ghost());
  inputz_.echange_espace_virtuel(inputz_.ghost());

}
void Validation_VDF::check(const DoubleTab& vdf_result)
{
  vdf_to_ijk_i.convert_to_ijk(vdf_result, vx_);
  vdf_to_ijk_j.convert_to_ijk(vdf_result, vy_);
  vdf_to_ijk_k.convert_to_ijk(vdf_result, vz_);

  compare(vx_, dvx_, msg_ + "vx");
  compare(vy_, dvy_, msg_ + "vy");
  compare(vz_, dvz_, msg_ + "vz");
}
#endif
