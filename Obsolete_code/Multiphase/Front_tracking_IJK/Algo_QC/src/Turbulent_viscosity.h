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
// File      : Turbulent_viscosity.h
// Directory : $NEW_ALGO_QC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Turbulent_viscosity_H
#define Turbulent_viscosity_H
#include <FixedVector.h>

class Turbulent_viscosity_base
{
public:
  virtual ~Turbulent_viscosity_base() {};
  virtual double operator()(double base_constant,
                            double dx,
                            double dy,
                            double dz,
                            double rho,
                            FixedVector<FixedVector<double, 3>, 3>& g,
                            FixedVector<double, 3>& q) = 0;
};

class Turbulent_viscosity_constant : public Turbulent_viscosity_base
{
public:
  double operator()(double constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_unsrho : public Turbulent_viscosity_base
{
public:
  double operator()(double constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_smagorinsky : public Turbulent_viscosity_base
{
public:
  double operator()(double smagorinsky_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_vreman : public Turbulent_viscosity_base
{
public:
  double operator()(double vreman_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_sigma : public Turbulent_viscosity_base
{
public:
  double operator()(double sigma_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_wale : public Turbulent_viscosity_base
{
public:
  double operator()(double wale_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_amd : public Turbulent_viscosity_base
{
public:
  double operator()(double amd_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_amdnoclip : public Turbulent_viscosity_base
{
public:
  double operator()(double amdnoclip_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_amdscalar : public Turbulent_viscosity_base
{
public:
  double operator()(double amdscalar_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_amdscalarnoclip : public Turbulent_viscosity_base
{
public:
  double operator()(double amdscalarnoclip_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_rds : public Turbulent_viscosity_base
{
public:
  double operator()(double rds_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_vss : public Turbulent_viscosity_base
{
public:
  double operator()(double vss_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

class Turbulent_viscosity_kobayashi : public Turbulent_viscosity_base
{
public:
  double operator()(double kobayashi_constant,
                    double dx,
                    double dy,
                    double dz,
                    double rho,
                    FixedVector<FixedVector<double, 3>, 3>& g,
                    FixedVector<double, 3>& q);
};

#endif
