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
// File      : ParcoursIJKDir.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#ifndef ParcoursIJKDir_included
#define ParcoursIJKDir_included

#include <Objet_U.h>
#include <MonofluidVar.h>
#include <Boundary_Conditions_Thermique.h>
#include <IJK_Splitting.h>
#include <Parser.h>
#include <IJK_Interfaces.h>
#include <IJK_Lata_writer.h>
#include <Intersection_Interface_ijk.h>
#include <Ouvrir_fichier.h>
#include <TRUST_Ref.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class ParcoursIJKDir
//
// <Description of class ParcoursIJKDir>
//
/////////////////////////////////////////////////////////////////////////////

class ParcoursIJKDir : Objet_U
{

  Declare_instanciable_sans_constructeur( ParcoursIJKDir );

public:
  ParcoursIJKDir() {next_elem_.resize(3); indices_to_keep_.resize(2); };
  // void initialize();
  FixedVector<int, 3> elem( const int q ) const
  {
    FixedVector<int, 3> res;
    res[0] = i_ + q * next_elem_[0];
    res[1] = j_ + q * next_elem_[1];
    res[2] = k_ + q * next_elem_[2];
    return res;
  }
  const int& i() const {return i_;}
  const int& j() const {return j_;}
  const int& k() const {return k_;}
  const int& face() const {return dir_;}
  void set_dir(const int dir)
  {
    dir_ = dir;
    set_next_elem();
    set_indices_to_keep();
  }
  void set_elem(const int i1, const int j1, const int k1)
  {
    i_ = i1;
    j_ = j1;
    k_ = k1;
  }
  const ArrOfInt& get_indices_to_keep() const { return indices_to_keep_;}
  const ArrOfInt& get_normale_vec() const {return next_elem_;}

  double calculer_surface_face(const IJK_Grid_Geometry& geom) const
  {
    double surface = 1.;
    for (int c = 0; c < 2; c++)
      surface *= geom.get_constant_delta(indices_to_keep_[c]);
    return surface;
  }

protected:
  void set_next_elem();
  void set_indices_to_keep();
  int dir_;
  int i_;
  int j_;
  int k_;
  ArrOfInt next_elem_;
  ArrOfInt indices_to_keep_;
};

#endif /* ParcoursIJKDir_included */
