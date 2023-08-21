/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : IJK_Bubble_tools.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Bubble_tools_included
#define IJK_Bubble_tools_included

#include <IJK_Field.h>
#include <IJK_Interfaces.h>
#include <Maillage_FT_IJK.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Bubble_tools
//
// <Description of class IJK_Bubble_tools>
//
/////////////////////////////////////////////////////////////////////////////
void compute_bounding_box_fill_compo(const IJK_Interfaces& interfaces, DoubleTab& bounding_box, IJK_Field_double& eulerian_compo_connex);
void compute_rising_velocity(const FixedVector<IJK_Field_double, 3>& velocity, const IJK_Interfaces& interfaces,
                             const IJK_Field_double& eulerian_compo_connex_ns, const int& gravity_dir,
                             DoubleTab& rising_velocities, DoubleTab& rising_vectors);
void fill_rising_velocity(const IJK_Field_double& eulerian_compo_connex_ns, const DoubleTab& rising_velocities,
                          IJK_Field_double& eulerian_rising_velocity);

#endif /* IJK_Bubble_tools_included */
