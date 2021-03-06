//TRUST_NO_INDENT
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
// File      : simd_tools_nosse.h
// Directory : $IJK_ROOT/src/IJK/solveur_mg
//
/////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file simd_tools_nosse.h.P
//
// This is a non optimized and portable implementation of the SimdFloat SimdDouble etc... classes
// (eg, without any kind of simd vectorisation)
// Simd vector types are of size 1.
#ifndef SIMD_TOOLS_NOSSE_H
#define SIMD_TOOLS_NOSSE_H

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

// Wrapper functions to allocate simd aligned blocs.
// Description: allocates a memory bloc of given size with proper alignment for SIMD.
//  nosimd=>fall back to malloc()
template<typename T>
inline T* simd_malloc (std::size_t size)
{
	T* res_ptr = (T*)malloc(size * sizeof(T));
	return res_ptr;
}

// Description: frees a memory bloc previously allocated with simd_malloc()
inline void simd_free(void * ptr)
{
  free(ptr);
}

// uintptr_t should be defined in stdint.h
//  (this type is the result of pointer operations like ptr1-ptr2)
typedef uintptr_t uintptr_type;
// Description: returns 1 if pointer is aligned on size bytes, 0 otherwise
//  Warn: size must be a power of 2.
inline int aligned(const void *ptr, int size)
{
  return ((uintptr_type)ptr & (uintptr_type)(size-1)) == 0;
}

// Empty macro, no alignment required for simd data structures
#define _SimdAligned_

// Implementation for single precision type
#define Simd_floatSIZE 1
// .DESCRIPTION
// This class provides a generic access to simd operations on IA32 and Intel 64 architecture.
// Functionalities provided by the class are designed to match those provided by common
// processor architectures (Altivec, SSE, etc):
//  - load vector size aligned data from memory (SimdGet)
//  - getting x[i-1] and x[i+1] efficiently for finite difference algorithms 
//    (SimdGetAtLeft, SimdGetAtRight, etc)
//  - arithmetic operations (+ - * /)
//  - conditional affectation (SimdSelect)
// See simd_malloc() and simd_free() to allocate aligned blocs of memory.
class Simd_float
{
public:
  typedef float value_type;
 
  Simd_float() {};
  // Size of the vector
  static int size() {
    return 1;
  }

  void operator+=(Simd_float a) {
    data_ += a.data_;
  }
  void operator*=(Simd_float a) {
    data_ *= a.data_;
  }
  void operator-=(Simd_float a) {
    data_ -= a.data_;
  }

  // The type below is architecture specific.
  // Code using it will be non portable.
  float data_;
  // Commodity default constructor (provides implicit conversion)
  Simd_float(float x) : data_(x) {};
};

// Description: Returns the vector found at address data. 
//  data must be aligned for the architecture (see simd_malloc())
inline Simd_float SimdGet(const float *data)
{
  return *data;
}

// Description: Stores vector x at address data. 
//  data must be aligned for the architecture (see simd_malloc())
inline void SimdPut(float *data, Simd_float x)
{
  *data = x.data_;
}

// Description: Returns the vector x starting at adress data+1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation.
inline Simd_float SimdGetAtRight(const float *data)
{
  return data[1];
}

// Description: Returns the vector x starting at adress data-1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation.
inline Simd_float SimdGetAtLeft(const float *data)
{
  return data[-1];
}

// Description: Returns the vector left and center starting at adress data-1 and data
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation
inline void SimdGetLeftCenter(const float *data, Simd_float &left, Simd_float &center)
{
  left = data[-1];
  center = data[0];
}

// Description: Returns the vector center and right starting at adress data and data+1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation
inline void SimdGetCenterRight(const float *data, 
				   Simd_float &center,
				   Simd_float &right)
{
  center = data[0];
  right = data[1];
}

// Description: Returns the vectors left, center and right starting at adress data-1, data and data+1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs three vector loads and two shift operations
inline void SimdGetLeftCenterRight(const float *data, 
				   Simd_float &left,
				   Simd_float &center,
				   Simd_float &right)
{
  left = data[-1];
  center = data[0];
  right = data[1];
}

inline void SimdGetLeftleftLeftCenterRight(const float *data, 
					   Simd_float &leftleft,
					   Simd_float &left,
					   Simd_float &center,
					   Simd_float &right)
{
  leftleft = data[-2];
  left = data[-1];
  center = data[0];
  right = data[1];
}

// Description: returns a+b
inline Simd_float operator+(Simd_float a, Simd_float b)
{
  return a.data_ + b.data_;
}

// Description: returns a-b
inline Simd_float operator-(Simd_float a, Simd_float b)
{
  return a.data_ - b.data_;
}

// Description: returns a*b
inline Simd_float operator*(Simd_float a, Simd_float b)
{
  return a.data_ * b.data_;
}

inline Simd_float Simd_absolute_value(Simd_float a)
{
	return abs(a.data_);
}

// Description: This function performs the following operation on the vectors
// for (i=0; i<size())
//   if (x1[i] < x2[i])
//     result[i] = value_if_x1_lower_than_x2[i]
//   else
//     result[i] = value_otherwise[i]
inline Simd_float SimdSelect(Simd_float x1,
			       Simd_float x2,
			       Simd_float value_if_x1_lower_than_x2,
			       Simd_float value_otherwise)
{
  return (x1.data_ < x2.data_) ? value_if_x1_lower_than_x2 : value_otherwise;
}

// Returns a vector built with min(a[i],b[i]) (element wise)
inline Simd_float SimdMin(const Simd_float & a, const Simd_float  & b)
{
  return (a.data_ < b.data_) ? a : b;
}

// Returns a vector built with max(a[i],b[i]) (element wise)
inline Simd_float SimdMax(const Simd_float & a, const Simd_float  & b)
{
  return (a.data_ > b.data_) ? a : b;
}

// Returns a 12 bits accurate result of a/b
inline Simd_float SimdDivideLow(const Simd_float & a, const Simd_float & b)
{
  return a.data_ / b.data_;
}

// Returns a 22 bits accurate result of a/b
inline Simd_float SimdDivideMed(const Simd_float & a, const Simd_float & b)
{
  return a.data_ / b.data_;
}

// Returns a 22 bits accurate result of 1/b
inline Simd_float SimdReciprocalMed(const Simd_float & b)
{
  return 1. / b.data_;
}

// Implementation for double precision type
#define Simd_doubleSIZE 1
// .DESCRIPTION
// This class provides a generic access to simd operations on IA32 and Intel 64 architecture.
// Functionalities provided by the class are designed to match those provided by common
// processor architectures (Altivec, SSE, etc):
//  - load vector size aligned data from memory (SimdGet)
//  - getting x[i-1] and x[i+1] efficiently for finite difference algorithms 
//    (SimdGetAtLeft, SimdGetAtRight, etc)
//  - arithmetic operations (+ - * /)
//  - conditional affectation (SimdSelect)
// See simd_malloc() and simd_free() to allocate aligned blocs of memory.
class Simd_double
{
public:
  typedef double value_type;
 
  Simd_double() {};
  // Size of the vector
  static int size() {
    return 1;
  }

  void operator+=(Simd_double a) {
    data_ += a.data_;
  }
  void operator*=(Simd_double a) {
    data_ *= a.data_;
  }
  void operator-=(Simd_double a) {
    data_ -= a.data_;
  }

  // The type below is architecture specific.
  // Code using it will be non portable.
  double data_;
  // Commodity default constructor (provides implicit conversion)
  Simd_double(double x) : data_(x) {};
};

// Description: Returns the vector found at address data. 
//  data must be aligned for the architecture (see simd_malloc())
inline Simd_double SimdGet(const double *data)
{
  return *data;
}

// Description: Stores vector x at address data. 
//  data must be aligned for the architecture (see simd_malloc())
inline void SimdPut(double *data, Simd_double x)
{
  *data = x.data_;
}

// Description: Returns the vector x starting at adress data+1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation.
inline Simd_double SimdGetAtRight(const double *data)
{
  return data[1];
}

// Description: Returns the vector x starting at adress data-1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation.
inline Simd_double SimdGetAtLeft(const double *data)
{
  return data[-1];
}

// Description: Returns the vector left and center starting at adress data-1 and data
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation
inline void SimdGetLeftCenter(const double *data, Simd_double &left, Simd_double &center)
{
  left = data[-1];
  center = data[0];
}

// Description: Returns the vector center and right starting at adress data and data+1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs two vector loads and a shift operation
inline void SimdGetCenterRight(const double *data, 
				   Simd_double &center,
				   Simd_double &right)
{
  center = data[0];
  right = data[1];
}

// Description: Returns the vectors left, center and right starting at adress data-1, data and data+1
//  data must be aligned for the architecture (see simd_malloc())
//  The implementation usually needs three vector loads and two shift operations
inline void SimdGetLeftCenterRight(const double *data, 
				   Simd_double &left,
				   Simd_double &center,
				   Simd_double &right)
{
  left = data[-1];
  center = data[0];
  right = data[1];
}

inline void SimdGetLeftleftLeftCenterRight(const double *data, 
					   Simd_double &leftleft,
					   Simd_double &left,
					   Simd_double &center,
					   Simd_double &right)
{
  leftleft = data[-2];
  left = data[-1];
  center = data[0];
  right = data[1];
}

// Description: returns a+b
inline Simd_double operator+(Simd_double a, Simd_double b)
{
  return a.data_ + b.data_;
}

// Description: returns a-b
inline Simd_double operator-(Simd_double a, Simd_double b)
{
  return a.data_ - b.data_;
}

// Description: returns a*b
inline Simd_double operator*(Simd_double a, Simd_double b)
{
  return a.data_ * b.data_;
}

inline Simd_double Simd_absolute_value(Simd_double a)
{
	return abs(a.data_);
}

// Description: This function performs the following operation on the vectors
// for (i=0; i<size())
//   if (x1[i] < x2[i])
//     result[i] = value_if_x1_lower_than_x2[i]
//   else
//     result[i] = value_otherwise[i]
inline Simd_double SimdSelect(Simd_double x1,
			       Simd_double x2,
			       Simd_double value_if_x1_lower_than_x2,
			       Simd_double value_otherwise)
{
  return (x1.data_ < x2.data_) ? value_if_x1_lower_than_x2 : value_otherwise;
}

// Returns a vector built with min(a[i],b[i]) (element wise)
inline Simd_double SimdMin(const Simd_double & a, const Simd_double  & b)
{
  return (a.data_ < b.data_) ? a : b;
}

// Returns a vector built with max(a[i],b[i]) (element wise)
inline Simd_double SimdMax(const Simd_double & a, const Simd_double  & b)
{
  return (a.data_ > b.data_) ? a : b;
}

// Returns a 22 bits accurate result of a/b
inline Simd_double SimdDivideMed(const Simd_double & a, const Simd_double & b)
{
  return a.data_ / b.data_;
}

inline Simd_double SimdReciprocalMed(const Simd_double & b)
{
  return 1. / b.data_;
}


class Simd_int
{
public:
  typedef int value_type;
 
  Simd_int() {};
  static int size() {
    return 1;
  }
  void operator&=(Simd_int a) {
    data_ &= a.data_;
  }
  void operator|=(Simd_int a) {
    data_ |= a.data_;
  }

  // The type below is architecture specific.
  // Code using it will be non portable.
  int data_;

  // Commodity default constructor (provides implicit conversion)
  Simd_int(int x) { data_ = x; }
};

// Description: Returns the vector found at address data. 
//  data must be aligned for the architecture (see simd_malloc())
inline Simd_int SimdGet(const int *data)
{
  return *data;
}

// Description: Stores vector x at address data. 
//  data must be aligned for the architecture (see simd_malloc())
inline void SimdPut(int *data, Simd_int x)
{
  *data = x.data_;
}

// Returns 0 if x1!=x2 and value_if_equal if x1==x2
inline Simd_int SimdTestEqual(Simd_float x1, Simd_float x2, Simd_int value_if_equal)
{
  return (x1.data_ == x2.data_) ? value_if_equal.data_ : 0;
}

inline Simd_int SimdTestEqual(Simd_float x1, Simd_float x2, 
			      Simd_int value_if_equal, Simd_int value_if_not_equal)
{
  return (x1.data_ == x2.data_) ? value_if_equal.data_ : value_if_not_equal.data_;
}


inline Simd_int SimdSelect(Simd_float x1,
			   Simd_float x2,
			   Simd_int value_if_x1_lower_than_x2,
			   Simd_int value_otherwise)
{
  return (x1.data_ < x2.data_) ? value_if_x1_lower_than_x2.data_ : value_otherwise.data_;
}

// Look the code... sorry :)
inline void SimdCompareAndSetIfLower(const Simd_float & x_new, Simd_float & x, 
				     const Simd_int & i_new, Simd_int & i)
{
  if (x_new.data_ < x.data_) {
    x = x_new;
    i = i_new;
  }
}



#endif


