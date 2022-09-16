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
// File      : MonofluidVar.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <MonofluidVar.h>

void IJK_MonofluidVar::initialize(const double liqu_val, const double vap_val)
{
  liqu_val_ = liqu_val;
  vap_val_ = vap_val;
}

///////////////
// Accessors //
///////////////

double IJK_MonofluidVar::operator() (const double indic_l) const
{
  assert ((indic_l > -1.e-12) && (indic_l < 1. + 1.e-12));
  const double res = liqu_val_ * indic_l + vap_val_ * (1. - indic_l);
  return res;
}
double IJK_MonofluidVar::operator() (const IJK_Field_double& indic_l, const int i, const int j, const int k) const
{
  const double res = operator()(indic_l(i,j,k));
  return res;
}
double IJK_MonofluidVar::harmo(const double indic_l) const
{
  assert ((indic_l > -1.e-12) && (indic_l < 1. + 1.e-12));
  const double res = 1. / ( indic_l/liqu_val_ + (1.-indic_l)/vap_val_);
  return res;
}
double IJK_MonofluidVar::harmo(const IJK_Field_double indic_l, const int i, const int j, const int k) const
{
  const double res = harmo(indic_l(i,j,k));
  return res;
}
void IJK_MonofluidVar::harmo(const IJK_Field_double& indic_l, IJK_Field_double& res) const
{
  const int ni = indic_l.ni();
  const int nj = indic_l.nj();
  const int nk = indic_l.nk();
  // res.allocate(indic_l.get_splitting(), IJK_Splitting::ELEM, 2);
  for (int i=0; i < ni; i++)
    for (int j=0; j < nj; j++)
      for (int k=0; k < nk; k++)
        {
          res(i,j,k) = 1. / ( indic_l(i,j,k) / liqu_val_ + (1. - indic_l(i,j,k)) / vap_val_ ) ;
        }
  res.echange_espace_virtuel(res.ghost());
}

void IJK_MonofluidVar::operator()(const IJK_Field_double& indic_l, IJK_Field_double& res) const
{
  const int ni = indic_l.ni();
  const int nj = indic_l.nj();
  const int nk = indic_l.nk();
  // res.allocate(indic_l.get_splitting(), IJK_Splitting::ELEM, 2);
  for (int i=0; i < ni; i++)
    for (int j=0; j < nj; j++)
      for (int k=0; k < nk; k++)
        {
          res(i,j,k) = liqu_val_ * indic_l(i,j,k) + vap_val_ * (1. - indic_l(i,j,k));
        }
  res.echange_espace_virtuel(res.ghost());
}




///////////////
// Modifiers //
///////////////

IJK_MonofluidVar& IJK_MonofluidVar::operator*= (const IJK_MonofluidVar& mono2)
{
  liqu_val_ *= mono2.liqu();
  vap_val_ *= mono2.vap();
  return *this;
}

IJK_MonofluidVar& IJK_MonofluidVar::operator/= (const IJK_MonofluidVar& mono2)
{
  liqu_val_ /= mono2.liqu();
  vap_val_ /= mono2.vap();
  return *this;
}

IJK_MonofluidVar& IJK_MonofluidVar::operator+= (const IJK_MonofluidVar& mono2)
{
  liqu_val_ += mono2.liqu();
  vap_val_ += mono2.vap();
  return *this;
}

/////////////
// getters //
/////////////

double IJK_MonofluidVar::liqu() const { return liqu_val_; }
double IJK_MonofluidVar::vap() const { return vap_val_; }


// Non member functions


IJK_MonofluidVar operator* (const IJK_MonofluidVar& mono1, const IJK_MonofluidVar& mono2)
{
  IJK_MonofluidVar res;
  res.initialize(mono1.liqu() * mono2.liqu(), mono2.vap() * mono2.vap());
  return res;
}

IJK_MonofluidVar operator/ (const IJK_MonofluidVar& mono1, const IJK_MonofluidVar& mono2)
{
  IJK_MonofluidVar res;
  res.initialize(mono1.liqu() / mono2.liqu(), mono1.vap() / mono2.vap());
  return res;
}

IJK_MonofluidVar operator+ (const IJK_MonofluidVar& mono1, const IJK_MonofluidVar& mono2)
{
  IJK_MonofluidVar res;
  res.initialize(mono1.liqu() + mono2.liqu(), mono1.vap() + mono2.vap());
  return res;
}
