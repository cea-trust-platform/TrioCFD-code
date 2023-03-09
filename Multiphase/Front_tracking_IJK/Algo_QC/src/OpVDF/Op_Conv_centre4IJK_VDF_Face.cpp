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
// File      : Op_Conv_centre4IJK_VDF_Face.cpp
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
#include <Op_Conv_centre4IJK_VDF_Face.h>
#include <Statistiques.h>
#include <IJK_discretization.h>
#include <Interprete_bloc.h>

Implemente_instanciable(Op_Conv_centre4IJK_VDF_Face,"Op_Conv_centre4IJK_VDF_Face",Op_Conv_centre4_VDF_Face);

Sortie& Op_Conv_centre4IJK_VDF_Face::printOn(Sortie& s) const
{
  return s;
}

Entree& Op_Conv_centre4IJK_VDF_Face::readOn(Entree& s)
{
  return s;
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
              Process::Journal() << "check centre4 " << msg << " (i,j,k)="
                                 << i << "," << j << "," << k
                                 << " ref val=" << s << finl;
            }
        }
  maxd=Process::mp_max(maxd);
  Process::Journal() << "check centre4 max error = " << maxd << finl;
}
#endif
void Op_Conv_centre4IJK_VDF_Face::associer_vitesse(const Champ_base& vconv)
{
  vconv_ = vconv;
  Op_Conv_centre4_VDF_Face::associer_vitesse(vconv);
}

DoubleTab& Op_Conv_centre4IJK_VDF_Face::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr << " Op_Conv_centre4IJK_VDF_Face::ajouter !!!" << endl;

  // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
  const char * ijkdis_name = IJK_discretization::get_conventional_name();
  const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));

  const VDF_to_IJK& vdf_to_ijk_i = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_I);
  const VDF_to_IJK& vdf_to_ijk_j = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_J);
  const VDF_to_IJK& vdf_to_ijk_k = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_K);


  static Stat_Counter_Id counter = statistiques().new_counter(0, "Op_Conv_centre4IJK_VDF_Face::ajouter_ijk");
  static Stat_Counter_Id counter2 = statistiques().new_counter(0, "Op_Conv_centre4IJK_VDF_Face::ajouter");
  static Stat_Counter_Id counterconvert = statistiques().new_counter(0, "Op_Conv_centre4IJK_VDF_Face ijk conversion");
  statistiques().begin_count(counterconvert);
  vdf_to_ijk_i.convert_to_ijk(inco, inputx_);
  vdf_to_ijk_j.convert_to_ijk(inco, inputy_);
  vdf_to_ijk_k.convert_to_ijk(inco, inputz_);
  vdf_to_ijk_i.convert_to_ijk(vconv_.valeur().valeurs(), vx_);
  vdf_to_ijk_j.convert_to_ijk(vconv_.valeur().valeurs(), vy_);
  vdf_to_ijk_k.convert_to_ijk(vconv_.valeur().valeurs(), vz_);
  vx_.echange_espace_virtuel(1);
  vy_.echange_espace_virtuel(1);
  vz_.echange_espace_virtuel(1);
  inputx_.echange_espace_virtuel(2);
  inputy_.echange_espace_virtuel(2);
  inputz_.echange_espace_virtuel(2);
  statistiques().end_count(counterconvert);
#ifdef CHECK_VAL

  // Check
  DoubleVect toto(resu);
  {
    toto=-1e10;
    vdf_to_ijk_i.convert_from_ijk(inputx_, toto);
    vdf_to_ijk_j.convert_from_ijk(inputy_, toto);
    vdf_to_ijk_k.convert_from_ijk(inputz_, toto);
    const int n = toto.size_reelle();
    double m = 0, m2 = 0;
    for (int i = 0; i < n; i++)
      {
        double d = fabs(inco[i] - toto[i]);
        m = max(m, d);
        m2 = max(fabs(inco[i]), m2);
        if (d > 1e-8)
          {
            Journal() << "check ijk convert error, i=" << i << " src="<<inco[i]<<" dest=" << toto[i] << finl;
          }
      }
    Journal() << "check ijk convert error Maxval= " << m2 << " Max error = " << m << finl;

#if 0
    // Single value:
    vx_.data() = 0.;
    vy_.data() = 0.;
    vz_.data() = 0.;
    //vx_(3,3,3) = 1.;
    inputx_.data() = 0.;
    inputy_.data() = 0.;
    inputz_.data() = 0.;
    //vx_(3,3,3) = 1.;
    //inputx_(3,3,3) = 1.;
    //vz_(5,5,5) = 1.;
    //inputz_(5,5,5) = 1.;

    /*    vdf_to_ijk_i.convert_from_ijk(inputx_, ref_cast(DoubleTab;inco));
    vdf_to_ijk_j.convert_from_ijk(inputy_, ref_cst(DoubleTab,inco));
    vdf_to_ijk_k.convert_from_ijk(inputz_, ref_cst(DoubleTab,inco));
    vdf_to_ijk_i.convert_from_ijk(vx_, ref_cst(DoubleTab,vconv_.valeur().valeurs()));
    vdf_to_ijk_j.convert_from_ijk(vy_, ref_cst(DoubleTab,vconv_.valeur().valeurs()));
    vdf_to_ijk_k.convert_from_ijk(vz_, ref_cst(DoubleTab,vconv_.valeur().valeurs()));
    */
#endif
  }
#endif

  statistiques().begin_count(counter);
  op_ijk_.calculer(inputx_, inputy_, inputz_, vx_, vy_, vz_, dvx_, dvy_, dvz_);
  statistiques().end_count(counter);
  /*  statistiques().begin_count(counterconvert);

  vdf_to_ijk_i.convert_from_ijk(dvx_, resu, true);
     vdf_to_ijk_j.convert_from_ijk(dvy_, resu, true);
     vdf_to_ijk_k.convert_from_ijk(dvz_, resu, true);

  statistiques().end_count(counterconvert); */
  Cerr << " Temps pour operateur convection ijk (conversion): " << statistiques().last_time(counterconvert) << finl;
  Cerr << " Temps pour operateur convection ijk (calcul)    : " << statistiques().last_time(counter) << finl;

  statistiques().begin_count(counter2);
  toto = resu;
  Op_Conv_centre4_VDF_Face::calculer(inco, resu);
  statistiques().end_count(counter2);
  Cerr << " Temps pour operateur convection vdf: " << statistiques().last_time(counter2) << finl;
#ifdef CHECK_VAL
  {
    vdf_to_ijk_i.convert_to_ijk(resu, vx_);
    vdf_to_ijk_j.convert_to_ijk(resu, vy_);
    vdf_to_ijk_k.convert_to_ijk(resu, vz_);

    compare(vx_,dvx_,"vx");
    compare(vy_,dvy_,"vy");
    compare(vz_,dvz_,"vz");
  }
#endif
  resu += toto;
  return resu;
}

void Op_Conv_centre4IJK_VDF_Face::completer()
{
  Cerr << " Op_Conv_centre4IJK_VDF_Face::completer()" << endl;

  // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
  const char * ijkdis_name = IJK_discretization::get_conventional_name();
  const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));
  const IJK_Splitting& split = ijkdis.get_IJK_splitting();
  op_ijk_.initialize(split);


  const Domaine_VDF& domaine_vf = ref_cast(Domaine_VDF, equation().domaine_dis().valeur());

  op_ijk_.initialize(split);
  vx_.allocate(split, IJK_Splitting::FACES_I, 1);
  vy_.allocate(split, IJK_Splitting::FACES_J, 1);
  vz_.allocate(split, IJK_Splitting::FACES_K, 1);
  inputx_.allocate(split, IJK_Splitting::FACES_I, 2);
  inputy_.allocate(split, IJK_Splitting::FACES_J, 2);
  inputz_.allocate(split, IJK_Splitting::FACES_K, 2);
  dvx_.allocate(split, IJK_Splitting::FACES_I, 0);
  dvy_.allocate(split, IJK_Splitting::FACES_J, 0);
  dvz_.allocate(split, IJK_Splitting::FACES_K, 0);

  flux_bords_.resize(domaine_vf.nb_faces_bord(),dimension);

  Op_Conv_centre4_VDF_Face::completer();
}
