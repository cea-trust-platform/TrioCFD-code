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
// File      : IJK_Composantes_Connex.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <IJK_Composantes_Connex.h>
#include <IJK_FT.h>
#include <IJK_Interfaces.h>
#include <IJK_Bubble_tools.h>

Implemente_instanciable( IJK_Composantes_Connex, "IJK_Composantes_Connex", Objet_U ) ;

Sortie& IJK_Composantes_Connex::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& IJK_Composantes_Connex::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

int IJK_Composantes_Connex::initialize(const IJK_Splitting& splitting, const IJK_Interfaces& interfaces)
{
  int nalloc = 0;
  interfaces_ = &interfaces;

  eulerian_compo_connex_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
  nalloc += 1;
  eulerian_compo_connex_ft_.data() = -1.;
  eulerian_compo_connex_ft_.echange_espace_virtuel(eulerian_compo_connex_ft_.ghost());

  eulerian_compo_connex_ghost_ft_.allocate(ref_ijk_ft_->get_splitting_ft(), IJK_Splitting::ELEM, 2);
  nalloc += 1;
  eulerian_compo_connex_ghost_ft_.data() = -1.;
  eulerian_compo_connex_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_ghost_ft_.ghost());

  eulerian_compo_connex_ns_.allocate(splitting, IJK_Splitting::ELEM, 0);
  nalloc += 1;
  eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());

  eulerian_compo_connex_ghost_ns_.allocate(splitting, IJK_Splitting::ELEM, 0);
  nalloc += 1;
  eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());

  return nalloc;
}

void IJK_Composantes_Connex::associer(const IJK_FT_double& ijk_ft)
{
  ref_ijk_ft_ = ijk_ft;
}

void IJK_Composantes_Connex::compute_bounding_box_fill_compo_connex()
{
  compute_bounding_box_fill_compo(ref_ijk_ft_->itfce(),
                                  bounding_box_,
                                  min_max_larger_box_,
                                  eulerian_compo_connex_ft_,
                                  eulerian_compo_connex_ghost_ft_,
                                  bubbles_barycentre_);
  eulerian_compo_connex_ft_.echange_espace_virtuel(eulerian_compo_connex_ft_.ghost());
  eulerian_compo_connex_ghost_ft_.echange_espace_virtuel(eulerian_compo_connex_ghost_ft_.ghost());

  eulerian_compo_connex_ns_.data() = -1;
  eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());
  ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_compo_connex_ft_, eulerian_compo_connex_ns_);
  eulerian_compo_connex_ns_.echange_espace_virtuel(eulerian_compo_connex_ns_.ghost());

  eulerian_compo_connex_ghost_ns_.data() = -1;
  eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
  ref_ijk_ft_->redistribute_from_splitting_ft_elem(eulerian_compo_connex_ghost_ft_, eulerian_compo_connex_ghost_ns_);
  eulerian_compo_connex_ghost_ns_.echange_espace_virtuel(eulerian_compo_connex_ghost_ns_.ghost());
}
