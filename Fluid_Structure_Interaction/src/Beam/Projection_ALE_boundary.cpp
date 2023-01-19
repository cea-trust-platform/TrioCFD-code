/****************************************************************************
* Copyright (c) 2021, CEA
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
// File      : Projection_ALE_boundary.cpp
// Directory : $BEAM_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Projection_ALE_boundary.h>
#include <Zone_ALE.h>

Implemente_instanciable(Projection_ALE_boundary, "Projection_ALE_boundary", Interprete_geometrique_base ) ;
//XD  Projection_ALE_boundary interprete Projection_ALE_boundary 0 block to compute the projection of a modal function on a mobile boundary. Use to compute modal added coefficients in FSI.
//XD  attr dom ref_domaine dom 0 Name of domain.
//XD attr bloc bloc_lecture bloc 0 between the braces, you must specify the numbers of the mobile borders then list these mobile borders and indicate the modal function which must be projected on these boundaries.  NL2 Example:  Projection_ALE_boundary dom_name  { 1  boundary_name   3  0.sin(pi*x)*1.e-4 0. }

Sortie& Projection_ALE_boundary::printOn( Sortie& os ) const
{
  return Interprete::printOn( os );
}

Entree& Projection_ALE_boundary::readOn( Entree& is )
{
  return Interprete::readOn( is );
}

Entree& Projection_ALE_boundary::interpreter_(Entree& is)
{
  associer_domaine(is);
  Zone_ALE& dom=ref_cast(Zone_ALE, domaine());
  dom.reading_projection_ALE_boundary(is);
  return is;
}
