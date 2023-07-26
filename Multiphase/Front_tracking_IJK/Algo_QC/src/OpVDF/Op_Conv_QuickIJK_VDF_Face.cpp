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
// File      : Op_Conv_QuickIJK_VDF_Face.cpp
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#if 0
#include <Op_Conv_QuickIJK_VDF_Face.h>
#include <IJK_discretization.h>
#include <Interprete_bloc.h>

// This class is only here to validate the ijk implementation against the old vdf one


Implemente_instanciable(Op_Conv_QuickIJK_VDF_Face,"Op_Conv_QuickIJK_VDF_Face",Op_Conv_Quick_VDF_Face);

Sortie& Op_Conv_QuickIJK_VDF_Face::printOn(Sortie& s) const
{
  return s;
}

Entree& Op_Conv_QuickIJK_VDF_Face::readOn(Entree& s)
{
  return s;
}

void Op_Conv_QuickIJK_VDF_Face::associer_vitesse(const Champ_base& vconv)
{
  vconv_ = vconv;
  Op_Conv_Quick_VDF_Face::associer_vitesse(vconv);
}

DoubleTab& Op_Conv_QuickIJK_VDF_Face::calculer(const DoubleTab& inco, DoubleTab& resu) const
{
  valid_.prepare(inco, vconv_.valeur().valeurs());

  statistiques().begin_count(counter_ijk);
  op_ijk_.calculer(inputx_, inputy_, inputz_, vx_, vy_, vz_, dvx_, dvy_, dvz_);
  statistiques().end_count(counter_ijk);
  Cerr << " Temps pour operateur ijk : " << statistiques().last_time(counter_ijk) << finl;
  statistiques().begin_count(counter_vdf);
  Op_Conv_centre4_VDF_Face::calculer(inco, resu);
  statistiques().end_count(counter_vdf);
  Cerr << " Temps pour operateur vdf: " << statistiques().last_time(counter_vdf) << finl;

  valid_.check(resu);

  return resu;
}

DoubleTab& Op_Conv_QuickIJK_VDF_Face::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  DoubleTab tmp(resu);
  tmp = 0.;
  calculer(inco, tmp);
  resu += tmp;
  return resu;
}

void Op_Conv_QuickIJK_VDF_Face::completer()
{
  Cerr << " Op_Conv_QuickIJK_VDF_Face::completer()" << endl;

  valid_.init(que_suis_je());

  // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
  const char * ijkdis_name = IJK_discretization::get_conventional_name();
  const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));
  const IJK_Splitting& split = ijkdis.get_IJK_splitting();
  op_ijk_.initialize(split);

  const Domaine_VDF& domaine_vf = ref_cast(Domaine_VDF, equation().domaine_dis().valeur());
  flux_bords_.resize(domaine_vf.nb_faces_bord(),dimension);

  Op_Conv_Quick_VDF_Face::completer();
}
#endif
