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
// File      : Op_Dift_IJK_VDF_var_Face.cpp
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#include <Op_Dift_IJK_VDF_Face.h>
#include <Statistiques.h>
#include <Mod_turb_hyd_ss_maille_VDF.h>
#include <IJK_discretization.h>
#include <Interprete_bloc.h>
#include <SChaine.h>
#include <EChaine.h>

Implemente_instanciable(Op_Dift_IJK_VDF_Face, "Op_Dift_VDF_var_FaceIJK", Op_Dift_VDF_var_Face);

Sortie& Op_Dift_IJK_VDF_Face::printOn(Sortie& s) const
{
  return s;
}

Entree& Op_Dift_IJK_VDF_Face::readOn(Entree& s)
{

  return s;
}

void Op_Dift_IJK_VDF_Face::completer()
{
  Op_Dift_VDF_var_Face::completer();
  // get reference to k from the turbulence model
  // (copied from Op_Dift_VDF_var_Face::completer())
  const RefObjU& ref_model = equation().get_modele(TURBULENCE);
  if (!sub_type(Mod_turb_hyd_ss_maille, ref_model.valeur()))
    {
      Cerr << "Error in Op_Dift_IJK_VDF_Face::completer(): turbulence model is not hyd_ss_maille" << finl;
      exit();
    }
  const Mod_turb_hyd_ss_maille& turbulence_model = ref_cast(Mod_turb_hyd_ss_maille, ref_model.valeur());
  vdf_nut_ = turbulence_model.viscosite_turbulente();
  vdf_kenergy_ = turbulence_model.energie_cinetique_turbulente();

  // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
  const char * ijkdis_name = IJK_discretization::get_conventional_name();
  const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));
  const IJK_Splitting& split = ijkdis.get_IJK_splitting();

#if 1
  SChaine instructions;
  instructions << "{" << finl;
  instructions << "  bctype_kmin Symetrie" << finl;
  instructions << "  bctype_kmax Symetrie" << finl;
  instructions << "}" << finl;
  Cerr << "Interpretation de la chaine suivante:" << finl << instructions.get_str();

  EChaine is(instructions.get_str());

  is >> bc_;
#endif

  op_ijk_.initialize(split, bc_);
  vx_.allocate(split, IJK_Splitting::FACES_I, 1);
  vy_.allocate(split, IJK_Splitting::FACES_J, 1);
  vz_.allocate(split, IJK_Splitting::FACES_K, 1);
  nu_.allocate(split, IJK_Splitting::ELEM, 1);
  nut_.allocate(split, IJK_Splitting::ELEM, 1);
  kenergy_.allocate(split, IJK_Splitting::ELEM, 1);
  dvx_.allocate(split, IJK_Splitting::FACES_I, 0);
  dvy_.allocate(split, IJK_Splitting::FACES_J, 0);
  dvz_.allocate(split, IJK_Splitting::FACES_K, 0);

  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, equation().domaine_dis().valeur());
  flux_bords_.resize(domaine_vdf.nb_faces_bord(),dimension);
}

void Op_Dift_IJK_VDF_Face::associer_diffusivite(const Champ_base& nu)
{
  vdf_nu_ = nu;
  Op_Dift_VDF_var_Face::associer_diffusivite(nu);
}

#define CHECK_VAL
#ifdef CHECK_VAL

static void compare(const IJK_Field_double& ref, const IJK_Field_double& x, const char *msg)
{
  int ni = x.ni();
  int nj = x.nj();
  int nk = x.nk();
  double maxd = 0.;
  for(int k=0; k<nk; k++)
    for(int j=0; j<nj; j++)
      for(int i=0; i<ni; i++)
        {
          double d = fabs(x(i,j,k)-ref(i,j,k));
          maxd = max(d,maxd);
          if (/*fabs(x(i,j,k)) > 1e-15 || fabs(ref(i,j,k))>1e-15*/ d > 1e-15)
            {
              char s[100];
              snprintf(s, 100, "%10g  %10g      %5g", ref(i,j,k), x(i,j,k), d);
              Process::Journal() << "check diffu " << msg << " (i,j,k)="
                                 << i << "," << j << "," << k
                                 << " ref val= " << s << finl;
            }
        }
  maxd=Process::mp_max(maxd);
  Process::Journal() << "check difft faces max error = " << maxd << finl;
}
#endif

DoubleTab& Op_Dift_IJK_VDF_Face::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  static Stat_Counter_Id counter = statistiques().new_counter(0, "Op_Dift_IJK_VDF_Face::ajouter_ijk");
  static Stat_Counter_Id counter2 = statistiques().new_counter(0, "Op_Dift_IJK_VDF_Face::ajouter");
  static Stat_Counter_Id counterconvert = statistiques().new_counter(0, "Op_Dift_IJK_VDF_Face ijk conversion");

  // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
  const char * ijkdis_name = IJK_discretization::get_conventional_name();
  const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));

  const VDF_to_IJK& vdf_to_ijk_i = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_I);
  const VDF_to_IJK& vdf_to_ijk_j = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_J);
  const VDF_to_IJK& vdf_to_ijk_k = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_K);
  const VDF_to_IJK& vdf_to_ijk_elem = ijkdis.get_vdf_to_ijk(IJK_Splitting::ELEM);


  //#define SINGLEVAL
#ifdef SINGLEVAL
  vx_.data() = 0.;
  vy_.data() = 0.;
  vz_.data() = 0.;
  vx_(3,3,3) = 1.;
  vy_(3,3,10) = 1.;
  vz_(3,3,18) = 1.;
  vdf_to_ijk_i.convert_from_ijk(vx_, ref_cast(DoubleTab&,inco));
  vdf_to_ijk_j.convert_from_ijk(vy_, ref_cast(DoubleTab,inco));
  vdf_to_ijk_k.convert_from_ijk(vz_, ref_cast(DoubleTab,inco));
  ref_csat(DoubleTab,inco)).echange_espace_virtuel();
  //(DoubleTab&)vdf_kenergy_.valeur().valeurs() = 0.;
  //(DoubleTab&)vdf_nut_.valeur().valeurs() = 1.;
  //(DoubleTab&)vdf_nu_.valeur().valeurs() = 0.;

#endif
  Cerr << "Op_Dift_IJK_VDF_Face::ajouter" << endl;
  statistiques().begin_count(counter2);
  Op_Dift_VDF_var_Face::ajouter(inco, resu);
  statistiques().end_count(counter2);
  Cerr << " Temps pour operateur diffusion vdf: " << statistiques().last_time(counter2) << finl;

  statistiques().begin_count(counterconvert);
  vdf_to_ijk_i.convert_to_ijk(inco, vx_);
  vdf_to_ijk_j.convert_to_ijk(inco, vy_);
  vdf_to_ijk_k.convert_to_ijk(inco, vz_);
  vdf_to_ijk_elem.convert_to_ijk(vdf_nu_.valeur().valeurs(), nu_);
  vdf_to_ijk_elem.convert_to_ijk(vdf_nut_.valeur().valeurs(), nut_);
  vdf_to_ijk_elem.convert_to_ijk(vdf_kenergy_.valeur().valeurs(), kenergy_);

  vx_.echange_espace_virtuel(1);
  vy_.echange_espace_virtuel(1);
  vz_.echange_espace_virtuel(1);
  nu_.echange_espace_virtuel(1);
  nut_.echange_espace_virtuel(1);
  kenergy_.echange_espace_virtuel(1);
  statistiques().end_count(counterconvert);
  Cerr << " Temps pour operateur diffusion vdf IJK conversion: " << statistiques().last_time(counterconvert) << finl;

  statistiques().begin_count(counter);
  op_ijk_.calculer(vx_, vy_, vz_,
                   nu_, nut_, kenergy_,
                   dvx_, dvy_, dvz_);
  statistiques().end_count(counter2);
  Cerr << " Temps pour operateur diffusion vdf IJK: " << statistiques().last_time(counter) << finl;

#ifdef CHECK_VAL
  DoubleTab toto(resu);
  toto = 0.;
  Op_Dift_VDF_var_Face::ajouter(inco, toto);

  vdf_to_ijk_i.convert_to_ijk(toto, vx_);
  vdf_to_ijk_j.convert_to_ijk(toto, vy_);
  vdf_to_ijk_k.convert_to_ijk(toto, vz_);
  compare(vx_,dvx_,"vx");
  compare(vy_,dvy_,"vy");
  compare(vz_,dvz_,"vz");

#endif
  return resu;
}
