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
// File      : Filter_kernel.h
// Directory : $NEW_ALGO_QC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Filter_kernel_H
#define Filter_kernel_H
#include <FixedVector.h>
#include <IJK_Field.h>

class Filter_kernel_base
{
public:
  Filter_kernel_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : ghost_size_(ghost_size), n_mailles_(n_mailles), ponderation_(ponderation), normalisation_(normalisation) {}
  virtual ~Filter_kernel_base() {};

  virtual FixedVector<double, 21> inhomogeneous(const bool elem,
                                                const int k,
                                                const int kg,
                                                const int nktot,
                                                const double delta,
                                                const ArrOfDouble_with_ghost& delta_z) = 0;
  virtual FixedVector<double, 21> uniform(double delta,
                                          double dz) = 0;
  inline int ghost_size()
  {
    return ghost_size_;
  }
  inline double n_mailles()
  {
    return n_mailles_;
  }
  inline bool ponderation()
  {
    return ponderation_;
  }
  inline bool normalisation()
  {
    return normalisation_;
  }

  inline int size_uniform()
  {
    const int s = this->ghost_size();
    return 2*s+1;
  }
  inline int size_k_elem(const int kg, const int nktot)
  {
    const int s = this->ghost_size();
    return min(2*s+1, min(s+kg+2, s+nktot-kg+1));
  }
  inline int size_k_face(const int kg, const int nktot)
  {
    const int s = this->ghost_size();
    return min(2*s+1, min(s+kg+1, s+nktot-kg+1));
  }
  inline int shift_uniform()
  {
    const int s = this->ghost_size();
    return s;
  }
  inline int shift_k_elem(const int kg)
  {
    const int s = this->ghost_size();
    return min(s, kg+1);
  }
  inline int shift_k_face(const int kg)
  {
    const int s = this->ghost_size();
    return min(s, kg);
  }
  inline bool is_at_wall_elem(const int kg, const int nktot)
  {
    return (kg==-1 || kg==nktot);
  }
  inline bool is_at_wall_face(const int kg, const int nktot)
  {
    return (kg==0 || kg==nktot);
  }

private:
  int ghost_size_;
  double n_mailles_;
  bool ponderation_;
  bool normalisation_;
};

class Filter_kernel_box : public Filter_kernel_base
{
public:
  Filter_kernel_box(int ghost_size) : Filter_kernel_base(1, 3., true, true) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_13_13_base : public Filter_kernel_base
{
public:
  Filter_kernel_weight_13_13_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : Filter_kernel_base(ghost_size, n_mailles, ponderation, normalisation) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_13_13_pondere : public Filter_kernel_weight_13_13_base
{
public:
  Filter_kernel_weight_13_13_pondere(int ghost_size) : Filter_kernel_weight_13_13_base(1, 3., true, true) {}
};

class Filter_kernel_weight_13_13_sansponderation : public Filter_kernel_weight_13_13_base
{
public:
  Filter_kernel_weight_13_13_sansponderation(int ghost_size) : Filter_kernel_weight_13_13_base(1, 3., false, false) {}
};

class Filter_kernel_weight_13_13_conservatif : public Filter_kernel_weight_13_13_base
{
public:
  Filter_kernel_weight_13_13_conservatif(int ghost_size) : Filter_kernel_weight_13_13_base(1, 3., true, false) {}
};

class Filter_kernel_weight_12_14_base : public Filter_kernel_base
{
public:
  Filter_kernel_weight_12_14_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : Filter_kernel_base(ghost_size, n_mailles, ponderation, normalisation) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_12_14_pondere : public Filter_kernel_weight_12_14_base
{
public:
  Filter_kernel_weight_12_14_pondere(int ghost_size) : Filter_kernel_weight_12_14_base(1, 2., true, true) {}
};

class Filter_kernel_weight_12_14_sansponderation : public Filter_kernel_weight_12_14_base
{
public:
  Filter_kernel_weight_12_14_sansponderation(int ghost_size) : Filter_kernel_weight_12_14_base(1, 2., false, false) {}
};

class Filter_kernel_weight_12_14_conservatif : public Filter_kernel_weight_12_14_base
{
public:
  Filter_kernel_weight_12_14_conservatif(int ghost_size) : Filter_kernel_weight_12_14_base(1, 2., true, false) {}
};

class Filter_kernel_weight_23_16_base : public Filter_kernel_base
{
public:
  Filter_kernel_weight_23_16_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : Filter_kernel_base(ghost_size, n_mailles, ponderation, normalisation) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_23_16_pondere : public Filter_kernel_weight_23_16_base
{
public:
  Filter_kernel_weight_23_16_pondere(int ghost_size) : Filter_kernel_weight_23_16_base(1, 2., true, true) {}
};

class Filter_kernel_weight_23_16_sansponderation : public Filter_kernel_weight_23_16_base
{
public:
  Filter_kernel_weight_23_16_sansponderation(int ghost_size) : Filter_kernel_weight_23_16_base(1, 2., false, false) {}
};

class Filter_kernel_weight_23_16_conservatif : public Filter_kernel_weight_23_16_base
{
public:
  Filter_kernel_weight_23_16_conservatif(int ghost_size) : Filter_kernel_weight_23_16_base(1, 2., true, false) {}
};

class Filter_kernel_weight_14_14_18_base : public Filter_kernel_base
{
public:
  Filter_kernel_weight_14_14_18_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : Filter_kernel_base(ghost_size, n_mailles, ponderation, normalisation) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_14_14_18_pondere : public Filter_kernel_weight_14_14_18_base
{
public:
  Filter_kernel_weight_14_14_18_pondere(int ghost_size) : Filter_kernel_weight_14_14_18_base(2, 4., true, true) {}
};

class Filter_kernel_weight_14_14_18_sansponderation : public Filter_kernel_weight_14_14_18_base
{
public:
  Filter_kernel_weight_14_14_18_sansponderation(int ghost_size) : Filter_kernel_weight_14_14_18_base(2, 4., false, false) {}
};

class Filter_kernel_weight_14_14_18_conservatif : public Filter_kernel_weight_14_14_18_base
{
public:
  Filter_kernel_weight_14_14_18_conservatif(int ghost_size) : Filter_kernel_weight_14_14_18_base(2, 4., true, false) {}
};

class Filter_kernel_weight_15_15_15_base : public Filter_kernel_base
{
public:
  Filter_kernel_weight_15_15_15_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : Filter_kernel_base(ghost_size, n_mailles, ponderation, normalisation) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_15_15_15_pondere : public Filter_kernel_weight_15_15_15_base
{
public:
  Filter_kernel_weight_15_15_15_pondere(int ghost_size) : Filter_kernel_weight_15_15_15_base(2, 5., true, true) {}
};

class Filter_kernel_weight_15_15_15_sansponderation : public Filter_kernel_weight_15_15_15_base
{
public:
  Filter_kernel_weight_15_15_15_sansponderation(int ghost_size) : Filter_kernel_weight_15_15_15_base(2, 5., false, false) {}
};

class Filter_kernel_weight_15_15_15_conservatif : public Filter_kernel_weight_15_15_15_base
{
public:
  Filter_kernel_weight_15_15_15_conservatif(int ghost_size) : Filter_kernel_weight_15_15_15_base(2, 5., true, false) {}
};

class Filter_kernel_weight_16_16_16_112_base : public Filter_kernel_base
{
public:
  Filter_kernel_weight_16_16_16_112_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : Filter_kernel_base(ghost_size, n_mailles, ponderation, normalisation) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_16_16_16_112_pondere : public Filter_kernel_weight_16_16_16_112_base
{
public:
  Filter_kernel_weight_16_16_16_112_pondere(int ghost_size) : Filter_kernel_weight_16_16_16_112_base(3, 6., true, true) {}
};

class Filter_kernel_weight_16_16_16_112_sansponderation : public Filter_kernel_weight_16_16_16_112_base
{
public:
  Filter_kernel_weight_16_16_16_112_sansponderation(int ghost_size) : Filter_kernel_weight_16_16_16_112_base(3, 6., false, false) {}
};

class Filter_kernel_weight_16_16_16_112_conservatif : public Filter_kernel_weight_16_16_16_112_base
{
public:
  Filter_kernel_weight_16_16_16_112_conservatif(int ghost_size) : Filter_kernel_weight_16_16_16_112_base(3, 6., true, false) {}
};

class Filter_kernel_weight_14_38_base : public Filter_kernel_base
{
public:
  Filter_kernel_weight_14_38_base(int ghost_size, double n_mailles, bool ponderation, bool normalisation) : Filter_kernel_base(ghost_size, n_mailles, ponderation, normalisation) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

class Filter_kernel_weight_14_38_pondere : public Filter_kernel_weight_14_38_base
{
public:
  Filter_kernel_weight_14_38_pondere(int ghost_size) : Filter_kernel_weight_14_38_base(1, 3., true, true) {}
};

class Filter_kernel_weight_14_38_sansponderation : public Filter_kernel_weight_14_38_base
{
public:
  Filter_kernel_weight_14_38_sansponderation(int ghost_size) : Filter_kernel_weight_14_38_base(1, 3., false, false) {}
};

class Filter_kernel_weight_14_38_conservatif : public Filter_kernel_weight_14_38_base
{
public:
  Filter_kernel_weight_14_38_conservatif(int ghost_size) : Filter_kernel_weight_14_38_base(1, 3., true, false) {}
};

class Filter_kernel_laplacian : public Filter_kernel_base
{
public:
  Filter_kernel_laplacian(int ghost_size) : Filter_kernel_base(1, -9999., false, false) {}
  FixedVector<double, 21> inhomogeneous(const bool elem,
                                        const int k,
                                        const int kg,
                                        const int nktot,
                                        const double delta,
                                        const ArrOfDouble_with_ghost& delta_z);
  FixedVector<double, 21> uniform(double delta,
                                  double dz);
};

#endif
