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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Echange_contact_VDF_Zoom_base.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Cond_Lim
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Echange_contact_VDF_Zoom_base_included
#define Echange_contact_VDF_Zoom_base_included


#include <Echange_externe_impose.h>
class Front_VF;
class Zone_VDF;
class Faces;

/*! @brief classe : Echange_contact_VDF_Zoom_base
 *
 *
 *
 */

////////////////////////////////////////////////////////////////

class Echange_contact_VDF_Zoom_base  : public Echange_externe_impose
{

  Declare_base(Echange_contact_VDF_Zoom_base);
public :
  void mettre_a_jour(double ) override =0;
protected :
  double h_paroi;

};

int l_elem_bord(Faces& les_faces,int face);
DoubleVect& trace_face_raccord_distant(const Front_VF& fr_vf,const DoubleVect& y,DoubleVect& x);

DoubleTab& trace_raccord_distant(const Front_VF& fr_vf,const Zone_VDF& zvdf,const DoubleTab& y, DoubleTab& x);


#endif
