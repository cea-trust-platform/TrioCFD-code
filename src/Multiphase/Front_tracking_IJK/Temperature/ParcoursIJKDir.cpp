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
// File      : ParcoursIJKDir.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <ParcoursIJKDir.h>

Implemente_instanciable_sans_constructeur( ParcoursIJKDir, "ParcoursIJKDir", Objet_U ) ;

Sortie& ParcoursIJKDir::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& ParcoursIJKDir::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void ParcoursIJKDir::set_indices_to_keep()
{
  if (dir_ == 0)
    {
      indices_to_keep_[0] = 1;
      indices_to_keep_[1] = 2;
    }
  else if (dir_ == 1)
    {
      indices_to_keep_[0] = 2;
      indices_to_keep_[1] = 0;
    }
  else
    {
      indices_to_keep_[0] = 0;
      indices_to_keep_[1] = 1;
    }
}

void ParcoursIJKDir::set_next_elem()
{
  if (dir_ == 0)
    {
      next_elem_[0] = 1;
      next_elem_[1] = 0;
      next_elem_[2] = 0;
    }
  else if (dir_ == 1)
    {
      next_elem_[0] = 0;
      next_elem_[1] = 1;
      next_elem_[2] = 0;
    }
  else
    {
      next_elem_[0] = 0;
      next_elem_[1] = 0;
      next_elem_[2] = 1;
    }
}
