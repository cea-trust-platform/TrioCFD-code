/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Cond_Lim_Rayo.cpp
// Directory:   $TRUST_ROOT/src/Rayonnement
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_Lim_Rayo.h>
#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Zone_VF.h>
#include <Zone_Cl_dis_base.h>

void Cond_Lim_Rayo::associer_modele_rayo(Modele_Rayonnement_base& mod)
{
  //  Cerr<<"Cond_Lim_Rayo::associer_modele_rayo"<<finl;
  le_modele_rayo = ref_cast(Modele_Rayonnement_Milieu_Transparent, mod);
}
void Cond_Lim_Rayo::completer()
{
  Cerr<<"Cond_Lim_Rayo::doit etre surchargee"<<finl;
  assert(0);
  Process::exit();
}
void  Cond_Lim_Rayo::preparer_surface(const Frontiere_dis_base& fr ,const Zone_Cl_dis_base& zcl)
{
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF,fr);
  int ndeb = la_frontiere_VF.num_premiere_face();
  int nb_faces_bord = la_frontiere_VF.nb_faces();
  //dimensionner Teta_i et Surf_i a nb_faces_bord.
  Surf_i.resize(nb_faces_bord);
  Teta_i.resize(nb_faces_bord);
  // recuperation des surfaces de bords.
  const Zone_VF& zone=ref_cast(Zone_VF,zcl.zone_dis().valeur());
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    Surf_i[numfa]= zone.face_surfaces(numfa+ndeb);
}
