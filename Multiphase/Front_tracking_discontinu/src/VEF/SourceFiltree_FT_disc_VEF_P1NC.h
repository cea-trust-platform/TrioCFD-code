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
// File:        SourceFiltree_FT_disc_VEF_P1NC.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src/VEF
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#ifndef SourceFiltree_FT_disc_VEF_P1NC_included
#define SourceFiltree_FT_disc_VEF_P1NC_included

#include <Source_base.h>
#include <SourceFiltree_FT_disc_base.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
// class SourceFiltree_FT_disc_VEF_P1NC
// la Classe SourceFiltree permet d'ajouter un terme source agissant uniquement
// sur une phase donnee
// .SECTION voir aussi
//     Probleme_FT_Disc_gen
//
//////////////////////////////////////////////////////////////////////////////

class SourceFiltree_FT_disc_VEF_P1NC : public Source_base, public SourceFiltree_FT_disc_base
{

  Declare_instanciable(SourceFiltree_FT_disc_VEF_P1NC);
public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  void completer();

protected :
  Entree& lire(Entree& is);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zone_Cl_VEF;
  REF(Champ_Inc) Indic_;
};

#endif
