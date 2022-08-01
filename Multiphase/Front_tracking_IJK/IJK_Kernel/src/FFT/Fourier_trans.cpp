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
// File      : Fourier_trans.cpp
// Directory : $NEW_ALGO_QC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <Fourier_trans.h>
#include <IJK_Grid_Geometry.h>
#include <TRUSTTab.h>
#include <communications.h>
#include <Param.h>
#include <IJK_Navier_Stokes_tools.h> // pour initialiser les IJK_double
#include <Statistiques.h>
#include <IJK_Lata_writer.h>
#include <SFichier.h>
#include <EFichier.h>
#include <LecFicDistribue_sansnum.h>
Implemente_instanciable_sans_constructeur(Fourier_trans,"Fourier_trans",Objet_U);

Fourier_trans::Fourier_trans()
{
  N_echantillons_ = 0;
}

// ici on va # define tous les valeur intermediaire a calculer !!!
#define U_FOURIER 0
#define V_FOURIER 1
#define W_FOURIER 2

#define TRIADIQUE_INPLANE_I_FOURIER 3
#define TRIADIQUE_INPLANE_J_FOURIER 4
#define TRIADIQUE_INPLANE_K_FOURIER 5

#define TRIADIQUE_INTERPLANE_I_FOURIER 6
#define TRIADIQUE_INTERPLANE_J_FOURIER 7
#define TRIADIQUE_INTERPLANE_K_FOURIER 8

#define DEFORMATION_TURBULENTE_I_FOURIER 9
#define DEFORMATION_TURBULENTE_J_FOURIER 10
#define DEFORMATION_TURBULENTE_K_FOURIER 11

#define P_FOURIER 12
#define PSURRHO_FOURIER 13
#define DIFFUSION_PRESSION_RHOFLUC_ADERIVE_FOURIER 14

#define DIFFUSION_PRESSION_DERIVRHO_I_FOURIER 15
#define DIFFUSION_PRESSION_DERIVRHO_J_FOURIER 16
#define DIFFUSION_PRESSION_DERIVRHO_K_FOURIER 17

#define DISSIPATION_COMPRESSIBLE_UN_II_FOURIER 18
#define DISSIPATION_COMPRESSIBLE_UN_IJ_FOURIER 19
#define DISSIPATION_COMPRESSIBLE_UN_IK_FOURIER 20
#define DISSIPATION_COMPRESSIBLE_UN_JI_FOURIER 21
#define DISSIPATION_COMPRESSIBLE_UN_JJ_FOURIER 22
#define DISSIPATION_COMPRESSIBLE_UN_JK_FOURIER 23
#define DISSIPATION_COMPRESSIBLE_UN_KI_FOURIER 24
#define DISSIPATION_COMPRESSIBLE_UN_KJ_FOURIER 25
#define DISSIPATION_COMPRESSIBLE_UN_KK_FOURIER 26

#define DISSIPATION_COMPRESSIBLE_DEUX_II_FOURIER 27
#define DISSIPATION_COMPRESSIBLE_DEUX_IJ_FOURIER 28
#define DISSIPATION_COMPRESSIBLE_DEUX_IK_FOURIER 29
#define DISSIPATION_COMPRESSIBLE_DEUX_JI_FOURIER 30
#define DISSIPATION_COMPRESSIBLE_DEUX_JJ_FOURIER 31
#define DISSIPATION_COMPRESSIBLE_DEUX_JK_FOURIER 32
#define DISSIPATION_COMPRESSIBLE_DEUX_KI_FOURIER 33
#define DISSIPATION_COMPRESSIBLE_DEUX_KJ_FOURIER 34
#define DISSIPATION_COMPRESSIBLE_DEUX_KK_FOURIER 35

#define DISSIPATION_COMPRESSIBLE_TROIS_I_FOURIER 36
#define DISSIPATION_COMPRESSIBLE_TROIS_J_FOURIER 37
#define DISSIPATION_COMPRESSIBLE_TROIS_K_FOURIER 38

#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_I_FOURIER 39
#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_J_FOURIER 40
#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_K_FOURIER 41

#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_I_FOURIER 42
#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_J_FOURIER 43
#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_K_FOURIER 44

#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS_FOURIER 45

#define DIFFUSION_VISQUEUSE_DERIVRHO_UN_I_FOURIER 46
#define DIFFUSION_VISQUEUSE_DERIVRHO_UN_J_FOURIER 47
#define DIFFUSION_VISQUEUSE_DERIVRHO_UN_K_FOURIER 48

#define DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_I_FOURIER 49
#define DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_J_FOURIER 50
#define DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_K_FOURIER 51

#define DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_I_FOURIER 52
#define DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_J_FOURIER 53
#define DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_K_FOURIER 54

#define TERME_NON_CLASSE_I_FOURIER 55
#define TERME_NON_CLASSE_J_FOURIER 56
#define TERME_NON_CLASSE_K_FOURIER 57

#define DUIDXJ_II_FOURIER 58
#define DUIDXJ_IJ_FOURIER 59
#define DUIDXJ_IK_FOURIER 60
#define DUIDXJ_JI_FOURIER 61
#define DUIDXJ_JJ_FOURIER 62
#define DUIDXJ_JK_FOURIER 63
#define DUIDXJ_KI_FOURIER 64
#define DUIDXJ_KJ_FOURIER 65
#define DUIDXJ_KK_FOURIER 66

#define DIVU_FOURIER 67

#define DUIDXJ_II_NONCENTRE_FOURIER 68
#define DUIDXJ_IJ_NONCENTRE_FOURIER 69
#define DUIDXJ_IK_NONCENTRE_FOURIER 70
#define DUIDXJ_JI_NONCENTRE_FOURIER 71
#define DUIDXJ_JJ_NONCENTRE_FOURIER 72
#define DUIDXJ_JK_NONCENTRE_FOURIER 73
#define DUIDXJ_KI_NONCENTRE_FOURIER 74
#define DUIDXJ_KJ_NONCENTRE_FOURIER 75
#define DUIDXJ_KK_NONCENTRE_FOURIER 76

#define UIUJ_II_FOURIER 77
#define UIUJ_IJ_FOURIER 78
#define UIUJ_IK_FOURIER 79
//#define UIUJ_JI_FOURIER
#define UIUJ_JJ_FOURIER 80
#define UIUJ_JK_FOURIER 81
//#define UIUJ_KI_FOURIER
//#define UIUJ_KJ_FOURIER
#define UIUJ_KK_FOURIER 82

#define UJDUIDXJ_I_FOURIER 83
#define UJDUIDXJ_J_FOURIER 84
#define UJDUIDXJ_K_FOURIER 85

#define NVAL_INTER 86

// ICI ON DEFINIT LES TERMES DE SORTIE  !!!
#define RES_UU 0
#define RES_VV 1
#define RES_WW 2
#define RES_PROD_TURB 3
#define RES_TRIADIQUE_INPLANE 4
#define RES_TRIADIQUE_INTERPLANE 5
#define RES_DEFORMATION_TURBULENTE 6
#define RES_CORRELATION_PRESSION_DIVERGENCE 7
#define RES_P_WFLUC 8
#define RES_DIFFUSION_PRESSION_RHOFLUC_ADERIVE 9
#define RES_DIFFUSION_PRESSION_DERIVRHO 10
#define RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INPLANE 11
#define RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INTERPLANE 12
#define RES_DISSIPATION_INCOMPRESSIBLE_DEUX_AFOISNU 13
#define RES_DISSIPATION_COMPRESSIBLE_UN 14
#define RES_DISSIPATION_COMPRESSIBLE_DEUX 15
#define RES_DISSIPATION_COMPRESSIBLE_TROIS 16
#define RES_DIFFUSION_VISQUEUSE_INCOMPRESSIBLE_DEUX 17
#define RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN 18
#define RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX 19
#define RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS 20
#define RES_DIFFUSION_VISQUEUSE_DERIVRHO_UN 21
#define RES_DIFFUSION_VISQUEUSE_DERIVRHO_DEUX 22
#define RES_DIFFUSION_VISQUEUSE_DERIVRHO_TROIS 23
#define RES_TERME_NON_CLASSE 24
#define RES_CONVECTIONTURBULENTE 25
#define RES_TERMEPUREMENTSPECTRAL 26

#define N_RES 27

// Pour generer cette liste, envoyer la liste des #define dans
//  sed 's/#define /"/g;s/_MOY.*/",/g' >resultat
// et copier ci-dessous:
static const char *noms_moyennes[] =
{
  "RES_UU",
  "RES_VV",
  "RES_WW",
  "RES_PROD_TURB",
  "RES_TRIADIQUE_INPLANE",
  "RES_TRIADIQUE_INTERPLANE",
  "RES_DEFORMATION_TURBULENTE",
  "RES_CORRELATION_PRESSION_DIVERGENCE",
  "RES_P_WFLUC",
  "RES_DIFFUSION_PRESSION_RHOFLUC_ADERIVE",
  "RES_DIFFUSION_PRESSION_DERIVRHO",
  "RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INPLANE",
  "RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INTERPLANE",
  "RES_DISSIPATION_INCOMPRESSIBLE_DEUX_AFOISNU",
  "RES_DISSIPATION_COMPRESSIBLE_UN",
  "RES_DISSIPATION_COMPRESSIBLE_DEUX",
  "RES_DISSIPATION_COMPRESSIBLE_TROIS",
  "RES_DIFFUSION_VISQUEUSE_INCOMPRESSIBLE_DEUX",
  "RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN",
  "RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX",
  "RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS",
  "RES_DIFFUSION_VISQUEUSE_DERIVRHO_UN",
  "RES_DIFFUSION_VISQUEUSE_DERIVRHO_DEUX",
  "RES_DIFFUSION_VISQUEUSE_DERIVRHO_TROIS",
  "RES_TERME_NON_CLASSE",
  "RES_CONVECTIONTURBULENTE",
  "RES_TERMEPUREMENTSPECTRAL"
};

/*
static inline double calculer_lambda_air(double temperature)
{
  const double fac_a = -5.05628e-18 ;
  const double fac_b = 2.469e-14 ;
  const double fac_c = -4.98344e-11 ;
  const double fac_d = 7.06714e-08 ;
  const double fac_e = 1.0894e-06 ;
  const double facteur = 1.93198026315789000e-3 / 1.461e-6;

  double val = temperature;
  double calc = val * fac_a + fac_b;
  calc = val * calc  + fac_c;
  calc = val * calc  + fac_d;
  calc = val * calc  + fac_e;
  return calc * facteur;
} */

/*
static inline double calculer_mu_air(double temperature)
{
  const double fac_a = -5.05628e-18 ;
  const double fac_b = 2.469e-14 ;
  const double fac_c = -4.98344e-11 ;
  const double fac_d = 7.06714e-08 ;
  const double fac_e = 1.0894e-06 ;

  double val = temperature;
  double calc = val * fac_a + fac_b;
  calc = val * calc  + fac_c;
  calc = val * calc  + fac_d;
  calc = val * calc  + fac_e;
  return calc;
} */

// calcule la derivee premiere en maillage anisotrope !!
static inline double derivee_aniso(  double alpha /* rapport des distances a la maille sup et inf */
                                     , double un_sur_alpha /* inverse de alpha */
                                     , double un_sur_dp_plus_dm /* inverse de la somme des distances aux elems inf et sup */
                                     , double val_moins /* valeur dans la maille inf */
                                     , double val_centre /* valeur dans la maille */
                                     , double val_plus /* valeur dans la maille sup */ )
{
  double derivee  = (val_plus   - val_centre)*un_sur_alpha ; // contribution sup pondere
  derivee += (val_centre - val_moins )*alpha ; // contribution moins pondere
  derivee *= un_sur_dp_plus_dm ; // division
  return derivee;
}

// static inline double derivee_2_aniso( double un_sur_delta_moins /* inverse da la distance a la maille inf */
//                                      , double un_sur_delta_plus /* inverse da la distance a la maille sup */
//                                      , double un_sur_dp_plus_dm /* inverse de la somme des distances aux elems inf et sup */
//                                      , double val_moins /* valeur dans la maille inf */
//                                      , double val_centre/* valeur dans la maille */
//                                      , double val_plus  /* valeur dans la maille sup */)
//{
//  double derivee_2  = (val_plus   - val_centre)*un_sur_delta_plus ; // contribution sup pondere
//  derivee_2 -= (val_centre - val_moins )*un_sur_delta_moins ; // contribution moins pondere
//  derivee_2 *= 2. * un_sur_dp_plus_dm ; // division (le 2. disparait naturellement lorsque dp = dm => 2* 1/ (dp + dm ) = 2*1/ 2d  = 1/d
//  return derivee_2;
//}

static void triABulle(ArrOfDouble& tab)
{
  int longueur = tab.size_array();

  bool permutation;

  do
    {
      permutation = false;
      for(int i=0; i<(longueur-1); i++)
        {
          if(tab[i]>tab[i+1])
            {
              double x = tab[i+1];
              double y = tab[i];
              tab[i] =x ;
              tab[i+1]=y;
              permutation = true;
            }
        }
      longueur--;
    }
  while(permutation);
}

static int retirer_doublons(ArrOfDouble& tab)
{
  int doublons = 0;
  const int taille = tab.size_array();

  /* ne marche que pour des tableaux trie */
  for (int i =0 ; i < (taille-1) ; i++)
    {
      const double x = tab[i];
      for ( int j = i+1 ; j < taille ; j++)
        {
          const double y = tab[j];
          if ( fabs(y-x) <= (y*1.e-3 ))
            {
              tab[j] = tab[taille-1]+i; // on va set a une valeur tres grande
              doublons++;
              i++;
            }
          else
            {
              break; // si on a trouver toute ls differences on s'arrette
            }
        }
    }
  return doublons;

}

static void remplir_ref_ki(ArrOfDouble& tab_k_tri ,DoubleTab& tab_k,IntTab& ref_k)
{
  const int Nvect = tab_k_tri.size_array();
  const int Ni = ref_k.dimension(0);
  const int Nj = ref_k.dimension(1);
  for (int i = 0 ; i < Ni; i++ )
    for (int j = 0 ; j < Nj ; j++)
      {
        double x = tab_k(i,j) ; // on cherche cette valeur
        bool trouver = false;
        for (int ou = 0; ou < Nvect ; ou++)
          {
            double y = tab_k_tri[ou];
            double compa = x > y ? x : y;
            if ( fabs(y-x) <= (compa*1.e-3) )
              {
                ref_k(i,j) = ou;
                trouver = true;
                break; // on a trouver donc on sort
              }
          }

        if ( !trouver )
          {
            Cerr << " Impossible de remplir correctement le tableau des references aux vecteurs, verifier la precision " << finl;
            Process::exit();
          }
      }
}

/* on met le signal au centre du domaine  pour pouvoir en faire une TF !!! */
static void recentrer_spectres(IJK_Field_double& field,IJK_Field_double& tmp )
{
  /* on prefere passer un tableau pour ne pas faire d'allocation (tres tres long) */
  int kmax= tmp.nk();
  int jmax= tmp.nj();
  int imax= tmp.ni();
  for(int k =0 ; k < kmax ; k++)
    for(int j = 0 ; j < jmax ; j++)
      for(int i =0 ; i < imax ; i++)
        tmp(i,j,k)=field(i,j,k);
  // pour l'instant on test, si ca marche on mettra ca directement dans la TF !!!!
  for(int k =0 ; k < kmax ; k++)
    for(int j = 0 ; j < jmax ; j++)
      for(int i =0 ; i < imax ; i++)
        {
          int i_cor = (i+imax/2)%imax;
          int j_cor = (j+jmax/2)%jmax;
          // if( i>=(imax/2)) { i_cor -= imax ; i_cor = abs(i_cor); }
          // if( j>=(jmax/2)) { j_cor -= jmax ; j_cor = abs(j_cor); }
          field(i_cor,j_cor,k)= tmp(i,j,k);
        }
}
// Attention, cette methode est appelee apres readOn(),
// il ne faut pas casser les donnees lues
void Fourier_trans::initialize(const IJK_Splitting& split_origine,const IJK_Splitting& split,VECT(ArrOfDouble) vmoy, ArrOfDouble rhomoy, ArrOfDouble numoy, double T_KMAX , double T_KMIN, int avec_reprise)
{

  static Stat_Counter_Id cnt_TFprepa = statistiques().new_counter(1, "TF_prepa");
  /* d'abord on initialise la redistribution */

  Post_splitting_ = split ;
  Cerr << "Initialisation du post-traitement spectral" << finl;
  redistribute_to_post_splitting_elem_.initialize(split_origine,Post_splitting_,IJK_Splitting::ELEM);
  redistribute_to_post_splitting_faces_[0].initialize(split_origine,Post_splitting_,IJK_Splitting::FACES_I);
  redistribute_to_post_splitting_faces_[1].initialize(split_origine,Post_splitting_,IJK_Splitting::FACES_J);
  redistribute_to_post_splitting_faces_[2].initialize(split_origine,Post_splitting_,IJK_Splitting::FACES_K);

  /* on recupere des infos sur le maillage */
  const IJK_Grid_Geometry& geom = Post_splitting_.get_grid_geometry();
  const int Nk_tot = geom.get_nb_elem_tot(DIRECTION_K);

  dx_=geom.get_constant_delta(0); //modif AT 20/06/2013
  dy_=geom.get_constant_delta(1); //modif AT 20/06/2013
  tab_dz_=geom.get_delta(2); //modif AT 20/06/2013
  elem_coord_.resize_array(Nk_tot);

  const ArrOfDouble& coord_z = geom.get_node_coordinates(DIRECTION_K);
  for (int i = 0; i < Nk_tot; i++)
    elem_coord_[i] = (coord_z[i] + coord_z[i+1]) * 0.5;

  /* DD 16/10/2015: on recopie le tableau des masses volumiques moyennes */
  /* il a deja ete adimensionne par statistiques_dns_jk_ */
  rho_moy_.resize_array(Nk_tot);
  rho_moy_ = rhomoy;
  /* DD 16/10/2015: on recopie le tableau des viscosites cinematiques moyennes */
  /* il a deja ete adimensionne par statistiques_dns_jk_ */
  nu_moy_.resize_array(Nk_tot);
  nu_moy_ = numoy;

  /* on recopie le tableau des vitesses moyennes */
  /* il a deja ete adimensionne par statistiques_dns_jk_ */
  vit_moy_.dimensionner(3);
  for (int i = 0; i < 3; i++)
    vit_moy_[i].resize_array(Nk_tot);

  vit_moy_[0] = vmoy[0];
  vit_moy_[1] = vmoy[1];
  vit_moy_[2] = vmoy[2];

  // recuperation des temperatures des CLs
  TCL_kmax_ = T_KMAX;
  TCL_kmin_ = T_KMIN;

  /* on dimensionne les tableaux de resultats */

  const int N = N_RES ;
  for (int i = 0 ; i < N ; i++)
    resultat_[i].allocate(Post_splitting_, IJK_Splitting::ELEM, 0);

  const int imax  = Post_splitting_.get_nb_elem_local(DIRECTION_I);
  const int jmax  = Post_splitting_.get_nb_elem_local(DIRECTION_J);
  const int kmax  = Post_splitting_.get_nb_elem_local(DIRECTION_K);
  statistiques().begin_count(cnt_TFprepa);
  //  allou les tableaux des entrees/sortie de la FFT
  in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(imax*jmax));
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(imax*jmax));

  plan = fftw_plan_dft_2d(jmax,imax, in, out,FFTW_FORWARD, FFTW_MEASURE );
  inverseplan = fftw_plan_dft_2d(jmax,imax, in, out,FFTW_BACKWARD, FFTW_MEASURE );

  const int Nval_inter = NVAL_INTER ; // chiffre au hasard pour le moment !!
  /*
    Div_U.allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
    Div_U.data()=0;


    DUiDx.allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
    DUiDx.data()=0;


    DUjDy.allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
    DUjDy.data()=0;


    DUkDz.allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
    DUkDz.data()=0;
    */
  allocate_velocity(vitesse, Post_splitting_, 1);
  champ_mu.allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
  pression.allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
  masse_vol.allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
  for ( int type = 0 ; type < Nval_inter ; type ++ )
    {
      Avant_TF[type].allocate(Post_splitting_, IJK_Splitting::ELEM, 1);
      Reel_TF[type].allocate(Post_splitting_,  IJK_Splitting::ELEM, 1);
      Imag_TF[type].allocate(Post_splitting_,  IJK_Splitting::ELEM, 1);
    }

  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   * ici on calcule les dimensions des tableaux liee a l'integrale
   * on commence par calculer la longueur des vecteurs d'ondes
   * Puis on tri les vecteurs
   * Puis on retire les doublons
   * Puis on retri le tableau (pour mettre les doublons a la fin
   * Enfin on recopie le tableau */


  if (avec_reprise)
    {
      reprise();
    }
  else
    {
      DoubleTab Taille_ki(imax,jmax);
      ArrOfDouble tem_k_tri(imax*jmax);

      const double coeff_x  = 1. / ( dx_ * (double)(imax) * dx_ * (double)(imax) );
      const double coeff_y  = 1. / ( dy_ * (double)(jmax) * dy_ * (double)(jmax) );

      for(int j = 0 ; j < jmax ; j++)
        for(int i =0 ; i < imax ; i++)
          {
            const int Indice = i + imax *j;
            int i_cor = i;
            int j_cor = j;
            if( i>=(imax/2))
              {
                i_cor -= imax ;
              }
            if( j>=(jmax/2))
              {
                j_cor -= jmax ;
              }
            double x = sqrt( double(i_cor*i_cor)*coeff_x + double(j_cor*j_cor)*coeff_y );
            Taille_ki(i,j) = x;
            tem_k_tri[Indice] = x;
          }
      // on copie le tableau pour apres

      // on tri le tableau et on supprime les doublons !
      triABulle(tem_k_tri);
      int doublons = retirer_doublons(tem_k_tri); // on set a un tres grand nombre les doublons;
      triABulle(tem_k_tri);
      const int taille = tem_k_tri.size_array()-doublons;
      tri_ki.resize_array(taille); // on dimensionne le tableau

      Cerr << " On part de " << imax*jmax << " vecteurs " << finl;
      Cerr << " On a retire " << doublons << " doublons  " << "il reste " <<  taille  << " vecteurs " << finl;

      for( int i = 0 ; i < tri_ki.size_array() ; i++)
        tri_ki[i]=tem_k_tri[i];

      // ici on definie un tableau qui indique comment les donnees vont etre utilisees pour l'integrale !!!
      ref_ki.resize(imax, jmax);
      ref_ki= -1;
      remplir_ref_ki(tri_ki,Taille_ki,ref_ki);

      TFx_.dimensionner(N);
      TFy_.dimensionner(N);
      TFi_.dimensionner(N);

      for ( int n = 0 ; n < N ; n++ )
        {
          TFx_[n].resize_array(imax/2*kmax);
          TFy_[n].resize_array(jmax/2*kmax);
          TFi_[n].resize_array(taille*kmax);
        }
    } // test pour reprise !!
  statistiques().end_count(cnt_TFprepa);
}


void Fourier_trans::update(const FixedVector<IJK_Field_double, 3>& vitesse_ori,
                           const IJK_Field_double& pression_ori,
                           //const IJK_Field_double &temperature,
                           const IJK_Field_double& masse_vol_ori,
                           const IJK_Field_double& champ_mu_ori)
{
  /* on dimensionne les tableaux d'origine */

  // tableau globaux qui vont contenir les nombres avant et apres TF
  // attention pas d'echange espace virtuel possible, dans le cas general ca crash !!!!!
  // ATTENTION ici on declare les tableaux pour stocker les variables temporaires necessaires aux calculs !!

  static Stat_Counter_Id cnt_redistribution = statistiques().new_counter(2, "TF Redistribution");
  statistiques().begin_count(cnt_redistribution);
  /* on recoit les tableaux dans le spliting de resolution */
  /* -------------- conversion -------------------*/
  /* ----------ET mise a jour des ghost ----------*/
  for (int i = 0; i < 3; i++)
    {
      redistribute_to_post_splitting_faces_[i].redistribute(vitesse_ori[i], vitesse[i]);
      vitesse[i].echange_espace_virtuel(vitesse[i].ghost()); // Ghost
    }
  redistribute_to_post_splitting_elem_.redistribute(masse_vol_ori,masse_vol);
  redistribute_to_post_splitting_elem_.redistribute(pression_ori,pression);
  redistribute_to_post_splitting_elem_.redistribute(champ_mu_ori,champ_mu);
  // Ghost
  masse_vol.echange_espace_virtuel(masse_vol.ghost());
  pression.echange_espace_virtuel(pression.ghost());
  champ_mu.echange_espace_virtuel(champ_mu.ghost());

  statistiques().end_count(cnt_redistribution);


  static Stat_Counter_Id cnt_tableaux = statistiques().new_counter(2, "TF tableau");
  statistiques().begin_count(cnt_tableaux);

  const int Nval = NVAL_INTER ; // chiffre au hasard pour le moment !!

  // copie des champs moyen existant
  const ArrOfDouble& vit_moy_i = vit_moy_[0];
  const ArrOfDouble& vit_moy_j = vit_moy_[1];
  const ArrOfDouble& vit_moy_k = vit_moy_[2];  // Ici c'est bien la vitesse aux faces !!!
  const ArrOfDouble& rho_moy = rho_moy_;  // DD 16/10/2015
  const ArrOfDouble& nu_moy = nu_moy_;  // DD 16/10/2015

  // Nombre total de mailles en K
  const int nktot = Post_splitting_.get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  // Nombre local de mailles en K
  const int imax = pression.ni();
  const int jmax = pression.nj();
  const int kmax = pression.nk();
  const int offsetk = Post_splitting_.get_offset_local(DIRECTION_K);

  const double unsurdx=1./dx_;
  const double unsurdy=1./dy_;

  // Copie des champs utiles des vitesses
  const IJK_Field_double& vitesse_i = vitesse[0];
  const IJK_Field_double& vitesse_j = vitesse[1];
  const IJK_Field_double& vitesse_k = vitesse[2];

  statistiques().end_count(cnt_tableaux);

  // on calcule les proprietes du fluides sur la Cl :
  /*
  double mu_kmin = calculer_mu_air(TCL_kmin_);
  double mu_kmax = calculer_mu_air(TCL_kmax_);

  double lambda_kmin = calculer_lambda_air(TCL_kmin);
  double lambda_kmax = calculer_lambda_air(TCL_kmax);
  */

  static Stat_Counter_Id cnt_calcul_termes = statistiques().new_counter(2, "TF_calcul_termes");
  statistiques().begin_count(cnt_calcul_termes);

  /* ******** BOUCLE PRINCIPALE ********* */
  /* ATTENTION VOLONTAIREMENT ON NE CALCULE PAS LES BORD  */

  for (int k = 0; k < kmax; k++)
    {
      // Calcul des moyennes spatiales sur le plan ij
      const int kg=k+offsetk; // position en k sur la grille globale
      const int paroi_basse = kg == 0 ? 1 : 0 ;
      const int paroi_haute = ( kg == (nktot-1) ) ? 1 : 0;
      const int plan_de_paroi = paroi_haute + paroi_basse ;

      if ( !plan_de_paroi )
        {
          const double unsurdz=1./tab_dz_[kg]; // taille de la maille dans la direction z
          double delta_m  = ( tab_dz_[kg] + tab_dz_[kg-1] ) * 0.5 ; // distance entre les centre de gravite avec gestion de la paroi !!!
          double delta_p  = ( tab_dz_[kg] + tab_dz_[kg+1] ) * 0.5 ; // distance entre les centre de gravite


          double alpha    = delta_p/delta_m ; // rapport des distance
          double unsalpha = delta_m/delta_p ; //
          double unsdpdm  = 1. / (  delta_p + delta_m );
          double unsdp    = 1./ delta_p;
          double unsdm    = 1./ delta_m;
          // On y stocke la somme de toutes les valeurs sur le plan ij, on divisera apres
          /* for (int j = 0; j < jmax; j++)
             for (int i = 0; i < imax; i++) */
          for (int j = 0; j < jmax; j++)
            for (int i = 0; i < imax; i++)
              {

                // calcul des fluctuations de vitesses aux faces.
                const double uf_i   = vitesse_i( i , j , k ) - vit_moy_i[kg];
                const double uf_ip1 = vitesse_i(i+1, j , k ) - vit_moy_i[kg];

                const double vf_j   = vitesse_j( i , j , k ) - vit_moy_j[kg];
                const double vf_jp1 = vitesse_j( i ,j+1, k ) - vit_moy_j[kg];

                const double wf_k   = vitesse_k( i , j , k ) - vit_moy_k[kg];
                const double wf_kp1 = vitesse_k( i , j ,k+1) - vit_moy_k[kg+1];

                // Fluctuations de vitesses aux elems
                const double ue = 0.5*(uf_ip1+uf_i);
                const double ve = 0.5*(vf_jp1+vf_j);
                const double we = 0.5*(wf_kp1+wf_k);

                // vitesses aux elems
                const double uetot = 0.5*(vitesse_i(i+1,j,k) + vitesse_i(i,j,k));
                const double vetot = 0.5*(vitesse_j(i,j+1,k) + vitesse_j(i,j,k));
                const double wetot = 0.5*(vitesse_k(i,j,k+1) + vitesse_k(i,j,k));

                // calcul des grandeurs composites a partir de ces vitesse !!!
                const double duidx=(uf_ip1-uf_i)*unsurdx;
                const double dujdy=(vf_jp1-vf_j)*unsurdy;
                const double dukdz=(wf_kp1-wf_k)*unsurdz;
                const double divu= duidx+dujdy+dukdz;

                // calcul des grandeurs scalaire !!
                const double unsurrho=1./masse_vol(i,j,k);
                const double rho = masse_vol(i,j,k);
                const double mu = champ_mu(i,j,k);
                const double nu = mu * unsurrho;
                const double p = pression(i,j,k);
                const double unsurrho_moy = 1./ rho_moy[kg];
                const double rho_fluc = (rho - rho_moy[kg]);
                const double nu_fluc = (nu - nu_moy[kg]);

                const double drhodx = ( masse_vol(i+1,j,k) - masse_vol(i-1,j,k) ) * 0.5 * unsurdx;
                const double drhody = ( masse_vol(i,j+1,k) - masse_vol(i,j-1,k) ) * 0.5 * unsurdx;
                const double drhodz = derivee_aniso(alpha, unsalpha, unsdpdm, masse_vol(i,j,k-1), rho, masse_vol(i,j,k+1));

                // derivation direction z
                double Uf_mk_0 = kg == 0 ? 0.         : vitesse_i(i,j,k-1) - vit_moy_i[kg-1];
                double Uf_pk_0 = kg == (nktot-1) ? 0. : vitesse_i(i,j,k+1) - vit_moy_i[kg+1];
                double duidz_0 = derivee_aniso(alpha, unsalpha, unsdpdm, Uf_mk_0, uf_i, Uf_pk_0);

                double Uf_mk_1 = kg == 0 ? 0.         : vitesse_i(i+1,j,k-1) - vit_moy_i[kg-1];
                double Uf_pk_1 = kg == (nktot-1) ? 0. : vitesse_i(i+1,j,k+1) - vit_moy_i[kg+1];
                double duidz_1 = derivee_aniso(alpha, unsalpha, unsdpdm, Uf_mk_1, uf_i, Uf_pk_1);

                double Uf_mk_0_tot = kg == 0 ? 0.           : vitesse_i(i,j,k-1);
                double Uf_pk_0_tot = kg == (nktot-1) ? 0.   : vitesse_i(i,j,k+1);
                double duidz_0_tot = derivee_aniso(alpha, unsalpha, unsdpdm, Uf_mk_0_tot, vitesse_i(i,j,k), Uf_pk_0_tot);

                double Uf_mk_1_tot = kg == 0 ? 0.         : vitesse_i(i+1,j,k-1);
                double Uf_pk_1_tot = kg == (nktot-1) ? 0. : vitesse_i(i+1,j,k+1);
                double duidz_1_tot = derivee_aniso(alpha, unsalpha, unsdpdm, Uf_mk_1_tot, vitesse_i(i+1,j,k), Uf_pk_1_tot);

                double Vf_mk= kg == 0 ? 0.            : vitesse_j(i,j,k-1)-vit_moy_j[kg-1];
                double Vf_pk= kg == (nktot-1) ? 0.    : vitesse_j(i,j,k+1)-vit_moy_j[kg+1];
                double dujdz_0 = derivee_aniso(alpha, unsalpha, unsdpdm, Vf_mk, vf_j, Vf_pk);

                double Vf_mk_1= kg == 0 ? 0.          : vitesse_j(i,j+1,k-1)-vit_moy_j[kg-1]; // on gere la paroi
                double Vf_pk_1= kg == (nktot-1) ? 0.  : vitesse_j(i,j+1,k+1)-vit_moy_j[kg+1]; // on gere la paroi
                double dujdz_1 = derivee_aniso(alpha , unsalpha ,unsdpdm ,Vf_mk_1, vf_jp1 , Vf_pk_1);

                // Transfert Triadique
                // on commence par calculer les produits colineaires dans l'axe.
                double duiidx = (uf_ip1*uf_ip1 - uf_i*uf_i) * unsurdx;
                double dujjdy = (vf_jp1*vf_jp1 - vf_j*vf_j) * unsurdy;
                double dukkdz = (wf_kp1*wf_kp1 - wf_k*wf_k) * unsurdz;

                // puis les produit mixtes  on decompose la somme pour faire apparaitre la composante de div de u !!!!
                // pour u_i
                // on va calculer la valeur de duidxj sur les deux face qui portent uj puis faire la moyenne !!!
                double vjduidy    = vf_j   * 0.5 * unsurdy * ( (  uf_ip1 - vitesse_i(i+1,j-1, k ) + vit_moy_i[ kg ] )
                                                               + (  uf_i   - vitesse_i( i ,j-1, k ) + vit_moy_i[ kg ] ) );
                double vj_p1duidy = vf_jp1 * 0.5 * unsurdy * ( ( -uf_ip1 + vitesse_i(i+1,j+1, k ) - vit_moy_i[ kg ] )
                                                               + ( -uf_i   + vitesse_i( i ,j+1, k ) - vit_moy_i[ kg ] ) );

                double wkduidz    = wf_k   * 0.5 * unsdm   * ( (  uf_ip1 - vitesse_i(i+1, j ,k-1) + vit_moy_i[kg-1] )
                                                               + (  uf_i   - vitesse_i( i , j ,k-1) + vit_moy_i[kg-1] ) );
                double wk_p1duidz = wf_kp1 * 0.5 * unsdp   * ( ( -uf_ip1 + vitesse_i(i+1, j ,k+1) - vit_moy_i[kg+1] )
                                                               + ( -uf_i   + vitesse_i( i , j ,k+1) - vit_moy_i[kg+1] ) );

                double duijdy = ue * dujdy + 0.5 * ( vjduidy + vj_p1duidy );
                double duikdz = ue * dukdz + 0.5 * ( wkduidz + wk_p1duidz );

                // pout u_j
                double uidujdx    = uf_i   * 0.5 * unsurdx * ( (  vf_jp1 - vitesse_j(i-1,j+1, k ) + vit_moy_j[ kg ] )
                                                               + (  vf_j   - vitesse_j(i-1, j , k ) + vit_moy_j[ kg ] ) );
                double ui_p1dujdx = uf_ip1 * 0.5 * unsurdx * ( ( -vf_jp1 + vitesse_j(i+1,j+1, k ) - vit_moy_j[ kg ] )
                                                               + ( -vf_j   + vitesse_j(i+1, j , k ) - vit_moy_j[ kg ] ) );

                double wkdujdz    = wf_k   * 0.5 * unsdm   * ( (  vf_jp1 - vitesse_j( i ,j+1,k-1) + vit_moy_j[kg-1] )
                                                               + (  vf_j   - vitesse_j( i , j ,k-1) + vit_moy_j[kg-1] ) );
                double wk_p1dujdz = wf_kp1 * 0.5 * unsdp   * ( ( -vf_jp1 + vitesse_j( i ,j+1,k+1) - vit_moy_j[kg+1] )
                                                               + ( -vf_j   + vitesse_j( i , j ,k+1) - vit_moy_j[kg+1] ) );

                double dujidx = ve * duidx + 0.5 * ( uidujdx + ui_p1dujdx );
                double dujkdz = ve * dukdz + 0.5 * ( wkdujdz + wk_p1dujdz );

                // pour u_k
                double uidukdx    = uf_i   * 0.5 * unsurdx * ( (  wf_kp1 - vitesse_k(i-1, j ,k+1) + vit_moy_k[kg+1] )
                                                               + (  wf_k   - vitesse_k(i-1, j , k ) + vit_moy_k[ kg ] ) );
                double ui_p1dukdx = uf_ip1 * 0.5 * unsurdx * ( ( -wf_kp1 + vitesse_k(i+1, j ,k+1) - vit_moy_k[kg+1] )
                                                               + ( -wf_k   + vitesse_k(i+1, j , k ) - vit_moy_k[ kg ] ) );

                double vjdukdy    = vf_j   * 0.5 * unsurdy * ( (  wf_kp1 - vitesse_k( i ,j-1,k+1) + vit_moy_k[kg+1] )
                                                               + (  wf_k   - vitesse_k( i ,j-1, k ) + vit_moy_k[ kg ] ) );
                double vj_p1dukdy = vf_jp1 * 0.5 * unsurdy * ( ( -wf_kp1 + vitesse_k( i ,j+1,k+1) - vit_moy_k[kg+1] )
                                                               + ( -wf_k   + vitesse_k( i ,j+1, k ) - vit_moy_k[ kg ] ) );

                double dukidx = we * duidx + 0.5 * ( uidukdx + ui_p1dukdx );
                double dukjdy = we * dujdy + 0.5 * ( vjdukdy + vj_p1dukdy );

                // ---
                double Uf_mj = vitesse_i(i,j-1,k)-vit_moy_i[kg];
                double duidy_0 = (uf_i  - Uf_mj ) * unsurdy;

                double Vf_im =  vitesse_j(i-1,j,k)-vit_moy_j[kg];
                double dujdx_0 = (vf_j - Vf_im ) * unsurdx;

                double Wf_mi =  vitesse_k(i-1,j,k)-vit_moy_k[kg];
                double dukdx_0 = (wf_k - Wf_mi ) * unsurdx;

                double Wf_mipk = vitesse_k(i-1,j,k+1)-vit_moy_k[kg+1];
                double dukdx_1 = (wf_kp1 - Wf_mipk ) * unsurdx;

                double Wf_mj = vitesse_k(i,j-1,k)-vit_moy_k[kg];
                double dukdy_0 = (wf_k  - Wf_mj ) * unsurdy;

                double Wf_mj_pk = vitesse_k(i,j-1,k+1)-vit_moy_k[kg+1];
                double dukdy_1 = (wf_kp1  - Wf_mj_pk ) * unsurdy;

                // ---

                double nu_arrette_ij = 0.25*( ( champ_mu(i,j,k)      /masse_vol(i,j,k))
                                              + ( champ_mu(i-1,j,k)  /masse_vol(i-1,j,k))
                                              + ( champ_mu(i,j-1,k)  /masse_vol(i,j-1,k))
                                              + ( champ_mu(i-1,j-1,k)/masse_vol(i-1,j-1,k)) );

                double nu_face_i = 0.5 *( ( champ_mu(i,j,k)    /masse_vol(i,j,k))
                                          + ( champ_mu(i-1,j,k)/masse_vol(i-1,j,k)) );

                double nu_face_j =  0.5 *( ( champ_mu(i,j,k)    /masse_vol(i,j,k))
                                           + ( champ_mu(i,j-1,k)/masse_vol(i,j-1,k)) );

                // ------------------------------------------------------------
                // ============================================================
                const double duidz = 0.5*(duidz_0 + duidz_1);
                const double dujdz = 0.5*(dujdz_0 + dujdz_1);
                const double duidz_tot = 0.5*(duidz_0_tot + duidz_1_tot);
                const double dukdz_tot = (vitesse_k(i,j,k+1) - vitesse_k(i,j,k)) * unsurdz;
                const double divtotu = duidx + dujdy + dukdz_tot;

                const double Ve_mi = 0.5 * ( (vitesse_j(i-1,j,k) - vit_moy_j[kg]) + (vitesse_j(i-1,j+1,k) - vit_moy_j[kg]) );
                const double Ve_pi = 0.5 * ( (vitesse_j(i+1,j,k) - vit_moy_j[kg]) + (vitesse_j(i+1,j+1,k) - vit_moy_j[kg]) );
                const double dujdx = (Ve_pi  - Ve_mi) * 0.5 * unsurdx;

                const double We_mi = 0.5 * ( (vitesse_k(i-1,j,k) - vit_moy_k[kg]) + (vitesse_k(i-1,j,k+1) - vit_moy_k[kg+1]) );
                const double We_pi = 0.5 * ( (vitesse_k(i+1,j,k) - vit_moy_k[kg]) + (vitesse_k(i+1,j,k+1) - vit_moy_k[kg+1]) );
                const double dukdx = (We_pi  - We_mi) * 0.5 * unsurdx;

                const double Ue_mj = 0.5*( (vitesse_i(i,j-1,k) - vit_moy_i[kg]) + (vitesse_i(i+1,j-1,k) - vit_moy_i[kg]) );
                const double Ue_pj = 0.5*( (vitesse_i(i,j+1,k) - vit_moy_i[kg]) + (vitesse_i(i+1,j+1,k) - vit_moy_i[kg]) );
                const double duidy = (Ue_pj  - Ue_mj) * 0.5 * unsurdy;

                const double We_mj = 0.5*( (vitesse_k(i,j-1,k) - vit_moy_k[kg]) + (vitesse_k(i,j-1,k+1) - vit_moy_k[kg+1]) );
                const double We_pj = 0.5*( (vitesse_k(i,j+1,k) - vit_moy_k[kg]) + (vitesse_k(i,j+1,k+1) - vit_moy_k[kg+1]) );
                const double dukdy = (We_pj  - We_mj ) * 0.5 * unsurdy;

                // ------------------------------------------------------------
                // PRODUCTION
                // ------------------------------------------------------------
                // production incompressible
                // <\w{ux}*\ \w{uy}\> [d<Ux>/dy]

                // production thermique
                // <\w{uy}*\ \w{uy}\> [d<Uy>/dy]

                // ------------------------------------------------------------
                // ADVECTION
                // ------------------------------------------------------------
                // [<Uy>] [d<Ec>/dy] ; Ec=w{ui}* w{ui} /2

                // ------------------------------------------------------------
                // TRANSFERT TRIADIQUE
                // ------------------------------------------------------------
                // transfert triadique dans le plan
                // <w{ui}* w{d(ui uj)/dx_j}>, j=x,z
                const double duijdxj_i = duiidx + duijdy;
                const double duijdxj_j = dujidx + dujjdy;
                const double duijdxj_k = dukidx + dukjdy;

                // transfert triadique interplan
                // <w{ui}* w{d(ui uj)/dxj}>, j=y
                // const double duikdz_i = duikdz;
                // const double duikdz_j = dujkdz;
                // const double duikdz_k = dukkdz;

                // ------------------------------------------------------------
                // DEFORMATION TURBULENTE
                // ------------------------------------------------------------
                // <w{ui}* w{ui dujdx_j}>
                double uidujdxj_i = ue * divu;
                double uidujdxj_j = ve * divu;
                double uidujdxj_k = we * divu;

                // ------------------------------------------------------------
                // CORRELATION PRESSION-DILATATION
                // ------------------------------------------------------------
                // <\w{duidx_i}*\ w{P/rho}>
                const double p_unsrho = p * unsurrho;

                // ------------------------------------------------------------
                // DIFFUSION PAR LA PRESSION
                // ------------------------------------------------------------
                // diffusion par la pression incompressible
                // [1/<rho>] <[d/dy](\w{uy}*\ \w{P}\)>

                // diffusion par la pression compressible, en d<rho>/dy
                // [<w{uy}* w{P}>/<rho>2] [<d<rho>/dy]
                // meme terme que pour "diffusion par la pression incompressible"

                // diffusion par la pression compressible, en rho' et rho
                // <[d/dy]((w{uy}* w{P rho')/(<rho> rho)})>
                const double rhofluc_unsrho_unsrhomoy_P = p * rho_fluc * unsurrho_moy * unsurrho;

                // diffusion par la pression compressible, en drho
                // <w{ui}* w{P/rho2 drho/dxi}>
                const double unsrho_unsrho_P_drhodxi_i = p * unsurrho * unsurrho * drhodx;
                const double unsrho_unsrho_P_drhodxi_j = p * unsurrho * unsurrho * drhody;
                const double unsrho_unsrho_P_drhodxi_k = p * unsurrho * unsurrho * drhodz;

                // ------------------------------------------------------------
                // DISSIPATION
                // ------------------------------------------------------------
                // dissipation incompressible, premiere partie
                // [<nu>] <\w{dui/dx_j}*\ \w{dui/dx_j}\> = 2 [<nu>] [k2 <Ec>] ; Ec=w{ui}* w{ui} /2

                // dissipation incompressible, seconde partie
                // [<nu>] <\w{dui/dx_j}*\ \w{duj/dx_i}\>

                // dissipation compressible, premiere partie
                // <\w{dui/dx_j}*\ w{nu' dUi/dx_j}>
                const double somme_nufluc_dUidxj_ii = nu_fluc * duidx;
                const double somme_nufluc_dUidxj_ij = (nu_arrette_ij - nu_moy[kg]) * duidy_0;
                const double somme_nufluc_dUidxj_ik = (nu_face_i - nu_moy[kg]) * duidz_0_tot;
                const double somme_nufluc_dUidxj_ji = (nu_arrette_ij - nu_moy[kg]) * dujdx_0;
                const double somme_nufluc_dUidxj_jj = nu_fluc * dujdy;
                const double somme_nufluc_dUidxj_jk = (nu_face_j - nu_moy[kg]) * dujdz_0;
                const double somme_nufluc_dUidxj_ki = (nu_face_i - nu_moy[kg]) * ((dukdx_0 + dukdx_1)*0.5);
                const double somme_nufluc_dUidxj_kj = (nu_face_j - nu_moy[kg]) * ((dukdy_0 + dukdy_1)*0.5);
                const double somme_nufluc_dUidxj_kk = nu_fluc * dukdz_tot;

                // dissipation compressible, seconde partie
                // <\w{dui/dx_j}*\ w{nu' dUj/dx_i}>
                // meme terme que pour "dissipation compressible, premiere partie"
                const double somme_nufluc_dUjdxi_ii = nu_fluc * duidx;
                const double somme_nufluc_dUjdxi_ij = nu_fluc * dujdx;
                const double somme_nufluc_dUjdxi_ik = nu_fluc * dukdx;
                const double somme_nufluc_dUjdxi_ji = nu_fluc * duidy;
                const double somme_nufluc_dUjdxi_jj = nu_fluc * dujdy;
                const double somme_nufluc_dUjdxi_jk = nu_fluc * dukdy;
                const double somme_nufluc_dUjdxi_ki = nu_fluc * duidz_tot;
                const double somme_nufluc_dUjdxi_kj = nu_fluc * dujdz;
                const double somme_nufluc_dUjdxi_kk = nu_fluc * dukdz_tot;

                // dissipation compressible, troisieme partie
                // <\w{dui/dx_i}*\ w{2/3 nu' dUj/dx_j}>
                const double deuxtiers_nu_divU = (2./3.) * nu * divtotu;

                // ------------------------------------------------------------
                // DIFFUSION VISQUEUSE
                // ------------------------------------------------------------
                // diffusion visqueuse incompressible, premiere partie
                // [<nu>] [d2<Ec>/dy2] ; Ec=w{ui}* w{ui} /2

                // diffusion visqueuse incompressible, seconde partie
                // [<nu>] <[d/dy](\w{ui}*\ w{duy/dx_i})>
                // const double dukdx_i = dukdx;
                // const double dukdx_j = dukdy;
                // const double dukdx_k = dukdz;

                // diffusion visqueuse compressible, en d<nu>/dy, premiere partie
                // [d<nu>dy] [d<Ec>/dy] ; Ec=w{ui}* w{ui} /2

                // diffusion visqueuse compressible, en d<nu>/dy, seconde partie
                // [d<nu>dy] [d(w{ui}* w{duy/dx_i})/dy]
                // meme terme que pour "diffusion visqueuse incompressible, seconde partie"

                // diffusion visqueuse compressible, en nu' et nu, premiere partie
                // <[d/dy](\w{ui}*\ w{nu' dUidy})>
                const double somme_nufluc_dUidxk_i = nu_fluc * duidz_tot;
                const double somme_nufluc_dUidxk_j = nu_fluc * dujdz;
                const double somme_nufluc_dUidxk_k = nu_fluc * dukdz_tot;

                // diffusion visqueuse compressible, en nu' et nu, seconde partie
                // <[d/dy](\w{ui}*\ w{nu' dUydx_i})>
                const double somme_nufluc_dUkdxi_i = nu_fluc * dukdx;
                const double somme_nufluc_dUkdxi_j = nu_fluc * dukdy;
                const double somme_nufluc_dUkdxi_k = nu_fluc * dukdz_tot;

                // diffusion visqueuse compressible, en nu' et nu, troisieme partie
                // <2/3 [d/dy](\w{uy}*\ w{nu dUydy})>
                const double somme_deuxtiers_nu_divU = (2./3.) * nu * divtotu;

                // diffusion visqueuse compressible, en drho, premiere partie
                // <\w{ui}*\ w{nu/rho dUi/dx_j drho/dx_j}>
                const double somme_nusrho_dUidxj_drhodxj_i = nu * unsurrho * (drhodx*duidx + drhody*duidy + drhodz*duidz_tot);
                const double somme_nusrho_dUidxj_drhodxj_j = nu * unsurrho * (drhodx*dujdx + drhody*dujdy + drhodz*dujdz);
                const double somme_nusrho_dUidxj_drhodxj_k = nu * unsurrho * (drhodx*dukdx + drhody*dukdy + drhodz*dukdz_tot);

                // diffusion visqueuse compressible, en drho, seconde partie
                // <\w{ui}*\ w{nu/rho dUj/dx_i drho/dx_j}>
                const double somme_nusrho_dUjdxi_drhodxj_i = nu * unsurrho * (drhodx*duidx     + drhody*dujdx + drhodz*dukdx);
                const double somme_nusrho_dUjdxi_drhodxj_j = nu * unsurrho * (drhodx*duidy     + drhody*dujdy + drhodz*dukdy);
                const double somme_nusrho_dUjdxi_drhodxj_k = nu * unsurrho * (drhodx*duidz_tot + drhody*dujdz + drhodz*dukdz_tot);

                // diffusion visqueuse compressible, en drho, troisieme partie
                // <\w{ui}*\ w{2/3 nu/rho dUj/dx_j drho/dx_i}>
                const double somme_deuxtiers_nusrho_divU_drhodxi_i = (2./3.) * nu * unsurrho * drhodx * divtotu;
                const double somme_deuxtiers_nusrho_divU_drhodxi_j = (2./3.) * nu * unsurrho * drhody * divtotu;
                const double somme_deuxtiers_nusrho_divU_drhodxi_k = (2./3.) * nu * unsurrho * drhodz * divtotu;

                // ------------------------------------------------------------
                // TERME NON CLASSE
                // ------------------------------------------------------------
                // <\w{ui}*\ w{1/rho Uj ui drho/dx_j}>
                const double somme_unsrho_ui_Uj_drhodxj_i = unsurrho * ue * (uetot*drhodx + vetot*drhody + wetot*drhodz);
                const double somme_unsrho_ui_Uj_drhodxj_j = unsurrho * ve * (uetot*drhodx + vetot*drhody + wetot*drhodz);
                const double somme_unsrho_ui_Uj_drhodxj_k = unsurrho * we * (uetot*drhodx + vetot*drhody + wetot*drhodz);

                // ------------------------------------------------------------
                // ============================================================

#define REMPLIR(quoi,avec) Avant_TF[quoi](i,j,k) = avec

                REMPLIR(U_FOURIER, ue);
                REMPLIR(V_FOURIER, ve);
                REMPLIR(W_FOURIER, wf_k);

                REMPLIR(TRIADIQUE_INPLANE_I_FOURIER, - duijdxj_i);
                REMPLIR(TRIADIQUE_INPLANE_J_FOURIER, - duijdxj_j);
                REMPLIR(TRIADIQUE_INPLANE_K_FOURIER, - duijdxj_k);

                REMPLIR(TRIADIQUE_INTERPLANE_I_FOURIER, - duikdz);
                REMPLIR(TRIADIQUE_INTERPLANE_J_FOURIER, - dujkdz);
                REMPLIR(TRIADIQUE_INTERPLANE_K_FOURIER, - dukkdz);

                REMPLIR(DEFORMATION_TURBULENTE_I_FOURIER, uidujdxj_i);
                REMPLIR(DEFORMATION_TURBULENTE_J_FOURIER, uidujdxj_j);
                REMPLIR(DEFORMATION_TURBULENTE_K_FOURIER, uidujdxj_k);

                REMPLIR(P_FOURIER, p);
                REMPLIR(PSURRHO_FOURIER, p_unsrho);
                REMPLIR(DIFFUSION_PRESSION_RHOFLUC_ADERIVE_FOURIER, rhofluc_unsrho_unsrhomoy_P);

                REMPLIR(DIFFUSION_PRESSION_DERIVRHO_I_FOURIER, - unsrho_unsrho_P_drhodxi_i);
                REMPLIR(DIFFUSION_PRESSION_DERIVRHO_J_FOURIER, - unsrho_unsrho_P_drhodxi_j);
                REMPLIR(DIFFUSION_PRESSION_DERIVRHO_K_FOURIER, - unsrho_unsrho_P_drhodxi_k);

                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_II_FOURIER, - somme_nufluc_dUidxj_ii);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_IJ_FOURIER, - somme_nufluc_dUidxj_ij);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_IK_FOURIER, - somme_nufluc_dUidxj_ik);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_JI_FOURIER, - somme_nufluc_dUidxj_ji);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_JJ_FOURIER, - somme_nufluc_dUidxj_jj);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_JK_FOURIER, - somme_nufluc_dUidxj_jk);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_KI_FOURIER, - somme_nufluc_dUidxj_ki);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_KJ_FOURIER, - somme_nufluc_dUidxj_kj);
                REMPLIR(DISSIPATION_COMPRESSIBLE_UN_KK_FOURIER, - somme_nufluc_dUidxj_kk);

                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_II_FOURIER, - somme_nufluc_dUjdxi_ii);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_IJ_FOURIER, - somme_nufluc_dUjdxi_ij);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_IK_FOURIER, - somme_nufluc_dUjdxi_ik);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_JI_FOURIER, - somme_nufluc_dUjdxi_ji);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_JJ_FOURIER, - somme_nufluc_dUjdxi_jj);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_JK_FOURIER, - somme_nufluc_dUjdxi_jk);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_KI_FOURIER, - somme_nufluc_dUjdxi_ki);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_KJ_FOURIER, - somme_nufluc_dUjdxi_kj);
                REMPLIR(DISSIPATION_COMPRESSIBLE_DEUX_KK_FOURIER, - somme_nufluc_dUjdxi_kk);

                REMPLIR(DISSIPATION_COMPRESSIBLE_TROIS_I_FOURIER, deuxtiers_nu_divU);
                REMPLIR(DISSIPATION_COMPRESSIBLE_TROIS_J_FOURIER, deuxtiers_nu_divU);
                REMPLIR(DISSIPATION_COMPRESSIBLE_TROIS_K_FOURIER, deuxtiers_nu_divU);

                REMPLIR(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_I_FOURIER, somme_nufluc_dUidxk_i);
                REMPLIR(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_J_FOURIER, somme_nufluc_dUidxk_j);
                REMPLIR(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_K_FOURIER, somme_nufluc_dUidxk_k);

                REMPLIR(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_I_FOURIER, somme_nufluc_dUkdxi_i);
                REMPLIR(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_J_FOURIER, somme_nufluc_dUkdxi_j);
                REMPLIR(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_K_FOURIER, somme_nufluc_dUkdxi_k);

                REMPLIR(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS_FOURIER,  - somme_deuxtiers_nu_divU);

                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_UN_I_FOURIER, somme_nusrho_dUidxj_drhodxj_i);
                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_UN_J_FOURIER, somme_nusrho_dUidxj_drhodxj_j);
                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_UN_K_FOURIER, somme_nusrho_dUidxj_drhodxj_k);

                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_I_FOURIER, somme_nusrho_dUjdxi_drhodxj_i);
                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_J_FOURIER, somme_nusrho_dUjdxi_drhodxj_j);
                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_K_FOURIER, somme_nusrho_dUjdxi_drhodxj_k);

                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_I_FOURIER, - somme_deuxtiers_nusrho_divU_drhodxi_i);
                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_J_FOURIER, - somme_deuxtiers_nusrho_divU_drhodxi_j);
                REMPLIR(DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_K_FOURIER, - somme_deuxtiers_nusrho_divU_drhodxi_k);

                REMPLIR(TERME_NON_CLASSE_I_FOURIER, somme_unsrho_ui_Uj_drhodxj_i);
                REMPLIR(TERME_NON_CLASSE_J_FOURIER, somme_unsrho_ui_Uj_drhodxj_j);
                REMPLIR(TERME_NON_CLASSE_K_FOURIER, somme_unsrho_ui_Uj_drhodxj_k);

                REMPLIR(DUIDXJ_II_FOURIER, duidx);
                REMPLIR(DUIDXJ_IJ_FOURIER, duidy);
                REMPLIR(DUIDXJ_IK_FOURIER, duidz);
                REMPLIR(DUIDXJ_JI_FOURIER, dujdx);
                REMPLIR(DUIDXJ_JJ_FOURIER, dujdy);
                REMPLIR(DUIDXJ_JK_FOURIER, dujdz);
                REMPLIR(DUIDXJ_KI_FOURIER, dukdx);
                REMPLIR(DUIDXJ_KJ_FOURIER, dukdy);
                REMPLIR(DUIDXJ_KK_FOURIER, dukdz);

                REMPLIR(DIVU_FOURIER, divu);

                // pour la dissipation incompressible_un ou compressible_un
                // on a pas besoin et l'on ne veut pas de derivees centrees
                REMPLIR(DUIDXJ_II_NONCENTRE_FOURIER, duidx);
                REMPLIR(DUIDXJ_IJ_NONCENTRE_FOURIER, duidy_0);
                REMPLIR(DUIDXJ_IK_NONCENTRE_FOURIER, duidz_0);
                REMPLIR(DUIDXJ_JI_NONCENTRE_FOURIER, dujdx_0);
                REMPLIR(DUIDXJ_JJ_NONCENTRE_FOURIER, dujdy);
                REMPLIR(DUIDXJ_JK_NONCENTRE_FOURIER, dujdz_0);
                REMPLIR(DUIDXJ_KI_NONCENTRE_FOURIER, (dukdx_0 + dukdx_1) * 0.5);
                REMPLIR(DUIDXJ_KJ_NONCENTRE_FOURIER, (dukdy_0 + dukdy_1) * 0.5);
                REMPLIR(DUIDXJ_KK_NONCENTRE_FOURIER, dukdz);


                /* -----------------------------------------------------------
                 * ajout de nouveaux terms suite a l'ecriture de part1 et
                 * l'identifications de decompositions plus pertinentes.
                 * [+] turbulent convection, purely spectral term.
                 * [-] triadic transfer.
                 */

                REMPLIR(UIUJ_II_FOURIER, ue * ue);
                REMPLIR(UIUJ_IJ_FOURIER, ue * ve);
                REMPLIR(UIUJ_IK_FOURIER, ue * we);
                //REMPLIR(UIUJ_JI_FOURIER, ve * ue);
                REMPLIR(UIUJ_JJ_FOURIER, ve * ve);
                REMPLIR(UIUJ_JK_FOURIER, ve * we);
                //REMPLIR(UIUJ_KI_FOURIER, we * ue);
                //REMPLIR(UIUJ_KJ_FOURIER, we * ve);
                REMPLIR(UIUJ_KK_FOURIER, we * we);

                REMPLIR(UJDUIDXJ_I_FOURIER, ue * duidx  +  0.5 * ( vjduidy + vj_p1duidy )  +  0.5 * ( wkduidz + wk_p1duidz ));
                REMPLIR(UJDUIDXJ_J_FOURIER, 0.5 * ( uidujdx + ui_p1dujdx )  +  ve * dujdy  +  0.5 * ( wkdujdz + wk_p1dujdz ));
                REMPLIR(UJDUIDXJ_K_FOURIER, 0.5 * ( uidukdx + ui_p1dukdx )  +  0.5 * ( vjdukdy + vj_p1dukdy )  +  we * dukdz);

#undef REMPLIR
              }
        }
      else // si on est a la paroi on set a zero !!
        {
          for (int j = 0; j < jmax; j++)
            for (int i = 0; i < imax; i++)
              for (int val = 0 ; val < Nval ; val++ )
                {
                  Avant_TF[val](i,j,k) =0;
                  Reel_TF[val](i,j,k)  =0;
                  Imag_TF[val](i,j,k)  =0 ;
                }
        }
    }

  // pour eviter de faire des echanges espaces virtuel (que je pense plus couteux que le calcul de quelques TF de plus, on va remplir les bord
  int KMAX = kmax;
  if ((kmax+offsetk)!=nktot)
    KMAX++;
  int KMIN = 0;
  if (offsetk!=0)// si je ne suis pas a la paroi je remplit le ghost
    KMIN--;


  {
    int k = KMAX-1;
    int kg = k + offsetk;
    for (int j = 0; j < jmax; j++)
      for (int i = 0; i < imax; i++)
        {
          double uf_i = vitesse_i(i,j,k) - vit_moy_i[kg];
          double vf_j = vitesse_j(i,j,k) - vit_moy_j[kg];
          double wf_k = vitesse_k(i,j,k) - vit_moy_k[kg];
          const double uf_ip1 = vitesse_i(i+1, j , k ) - vit_moy_i[kg];
          const double vf_jp1 = vitesse_j( i ,j+1, k ) - vit_moy_j[kg];
          const double ue = 0.5*(uf_ip1+uf_i);
          const double ve = 0.5*(vf_jp1+vf_j);
          double p = pression(i,j,k);
#define REMPLIR(quoi,avec) Avant_TF[quoi](i,j,k) = avec
          // vitesse
          REMPLIR(U_FOURIER,ue);
          REMPLIR(V_FOURIER,ve);
          REMPLIR(W_FOURIER,wf_k);
          REMPLIR(P_FOURIER,p);
#undef REMPLIR
        }

    k = KMIN;
    kg = k + offsetk;
    for (int j = 0; j < jmax; j++)
      for (int i = 0; i < imax; i++)
        {
          double uf_i = vitesse_i(i,j,k) - vit_moy_i[kg];
          double vf_j = vitesse_j(i,j,k) - vit_moy_j[kg];
          double wf_k = vitesse_k(i,j,k) - vit_moy_k[kg];
          const double uf_ip1 = vitesse_i(i+1, j , k ) - vit_moy_i[kg];
          const double vf_jp1 = vitesse_j( i ,j+1, k ) - vit_moy_j[kg];
          const double ue = 0.5*(uf_ip1+uf_i);
          const double ve = 0.5*(vf_jp1+vf_j);
          double p = pression(i,j,k);
#define REMPLIR(quoi,avec) Avant_TF[quoi](i,j,k) = avec
          // vitesse
          REMPLIR(U_FOURIER,ue);
          REMPLIR(V_FOURIER,ve);
          REMPLIR(W_FOURIER,wf_k);
          REMPLIR(P_FOURIER,p);
#undef REMPLIR
        }
  }
  statistiques().end_count(cnt_calcul_termes);

  static Stat_Counter_Id cnt_TFpure = statistiques().new_counter(2, "TF_PURE");
  statistiques().begin_count(cnt_TFpure);
  // Ici on a calculer les valeur maintenant il faut les envoyer aux autres processeur !!!
  // calcule la TF 2D de chaque grandeur et la stocke
  /* EXPLICATION : le retour de la TF est un plan 2D de Ni,Nj points.
   * Ainsi, il se stock naturellement dans un tableau ijk !!!
   * pour se balader le long des k1 on se deplace sur les i pour les k2 on se deplacera sur j */
  for ( int type = 0 ; type < Nval ; type ++ )
    Traitement_spectral(Avant_TF[type],Reel_TF[type],Imag_TF[type],1);


  {
    /* ghost et bord des vitesses et de la pression */
    Traitement_spectral_klayer_direct(Avant_TF[U_FOURIER],Reel_TF[U_FOURIER],Imag_TF[U_FOURIER],KMIN);
    Traitement_spectral_klayer_direct(Avant_TF[U_FOURIER],Reel_TF[U_FOURIER],Imag_TF[U_FOURIER],KMAX-1);

    Traitement_spectral_klayer_direct(Avant_TF[V_FOURIER],Reel_TF[V_FOURIER],Imag_TF[V_FOURIER],KMIN);
    Traitement_spectral_klayer_direct(Avant_TF[V_FOURIER],Reel_TF[V_FOURIER],Imag_TF[V_FOURIER],KMAX-1);

    Traitement_spectral_klayer_direct(Avant_TF[W_FOURIER],Reel_TF[W_FOURIER],Imag_TF[W_FOURIER],KMIN);
    Traitement_spectral_klayer_direct(Avant_TF[W_FOURIER],Reel_TF[W_FOURIER],Imag_TF[W_FOURIER],KMAX-1);

    Traitement_spectral_klayer_direct(Avant_TF[P_FOURIER],Reel_TF[P_FOURIER],Imag_TF[P_FOURIER],KMIN);
    Traitement_spectral_klayer_direct(Avant_TF[P_FOURIER],Reel_TF[P_FOURIER],Imag_TF[P_FOURIER],KMAX-1);
  }
  static Stat_Counter_Id TF_Combinaisons = statistiques().new_counter(2, "TF Combinaisons");
  statistiques().begin_count(TF_Combinaisons);
  /* Combinaisons des termes spectraux */
  /* Pour le moment, je ne fait pas de fonction car c'est trop compliquer a implemnter a cause du FixedVector
  Combiner_spectres(Reel_TF,Imag_TF); */
  statistiques().end_count(cnt_TFpure);

  for (int k = 0; k < kmax; k++)
    {
      const int kg=k+offsetk; // position en k sur la grille globale
      const int paroi_basse = kg == 0 ? 1 : 0 ;
      const int paroi_haute = ( kg == (nktot-1) ) ? 1 : 0;
      const int plan_de_paroi = paroi_haute + paroi_basse ;
      /*
       * const double unsurdz = 1./tab_dz_[kg];
       * // AUTRE FACON DE CALCULER Lx
       * // const IJK_Grid_Geometry & geom = Post_splitting_.get_grid_geometry();
       * // const int Nx_tot = geom.get_nb_elem_tot(0) + 1;
       * // const int Ny_tot = geom.get_nb_elem_tot(1) + 1;
       * // const double Lx = geom.get_node_coordinates(0)[Nx_tot - 1] - geom.get_origin(0);
       * // const double Ly = geom.get_node_coordinates(1)[Ny_tot - 1] - geom.get_origin(1);
       * const double Lx = dx_*imax;
       * const double Ly = dy_*jmax;
       */
      if ( !plan_de_paroi )
        {
          /*
           * double delta_m  = ( tab_dz_[kg] + tab_dz_[kg-1] ) * 0.5 ; // distance entre les centre de gravite avec gestion de la paroi !!!
           * double delta_p  = ( tab_dz_[kg] + tab_dz_[kg+1] ) * 0.5 ; // distance entre les centre de gravite
           * double alpha    = delta_p/delta_m;
           * double unsalpha = delta_m/delta_p;
           * double unsdpdm  = 1. / (  delta_p + delta_m );
           */
          for (int j = 0; j < jmax; j++)
            for (int i = 0; i < imax; i++)
              {
                /*
                      * const int i_cor = (i < imax/2) ? i : abs(i - imax);
                      * const int j_cor = (j < jmax/2) ? j : abs(j - jmax);
                      * double veconde_kx = ((double)(i_cor)/Lx)*2*3.14159265358979323846;
                      * double veconde_ky = ((double)(j_cor)/Ly)*2*3.14159265358979323846;
                      * //double veconde_k = sqrt(veconde_kx*veconde_kx + veconde_ky*veconde_ky);
                */

                double uu  = Reel_TF[U_FOURIER](i,j,k)*Reel_TF[U_FOURIER](i,j,k) + Imag_TF[U_FOURIER](i,j,k)*Imag_TF[U_FOURIER](i,j,k);
                double vv  = Reel_TF[V_FOURIER](i,j,k)*Reel_TF[V_FOURIER](i,j,k) + Imag_TF[V_FOURIER](i,j,k)*Imag_TF[V_FOURIER](i,j,k);
                double ww  = Reel_TF[W_FOURIER](i,j,k)*Reel_TF[W_FOURIER](i,j,k) + Imag_TF[W_FOURIER](i,j,k)*Imag_TF[W_FOURIER](i,j,k);
                ww += Reel_TF[W_FOURIER](i,j,k+1)*Reel_TF[W_FOURIER](i,j,k+1) + Imag_TF[W_FOURIER](i,j,k+1)*Imag_TF[W_FOURIER](i,j,k+1);
                ww *= 0.5;
                resultat_[RES_UU](i,j,k) = uu ;
                resultat_[RES_VV](i,j,k) = vv ; // fait aussi la prod_thermique !!!
                resultat_[RES_WW](i,j,k) = ww ;
                const double RUE = Reel_TF[U_FOURIER](i,j,k);
                const double IUE = Imag_TF[U_FOURIER](i,j,k);
                const double RVE = Reel_TF[V_FOURIER](i,j,k);
                const double IVE = Imag_TF[V_FOURIER](i,j,k);
                // la TF est au face alors on revient au elems
                const double RWE = 0.5 * (Reel_TF[W_FOURIER](i,j,k) + Reel_TF[W_FOURIER](i,j,k+1));
                const double IWE = 0.5 * (Imag_TF[W_FOURIER](i,j,k) + Imag_TF[W_FOURIER](i,j,k+1));

                // derivee de chaque composante de u dans la direction aniso
                /*
                       * const double deriv_U_re = derivee_aniso(alpha, unsalpha, unsdpdm, Reel_TF[U_FOURIER](i,j,k-1), Reel_TF[U_FOURIER](i,j,k), Reel_TF[U_FOURIER](i,j,k+1));
                       * const double deriv_U_im = derivee_aniso(alpha, unsalpha, unsdpdm, Imag_TF[U_FOURIER](i,j,k-1), Imag_TF[U_FOURIER](i,j,k), Imag_TF[U_FOURIER](i,j,k+1));
                       * const double deriv_V_re = derivee_aniso(alpha, unsalpha, unsdpdm, Reel_TF[V_FOURIER](i,j,k-1), Reel_TF[V_FOURIER](i,j,k), Reel_TF[V_FOURIER](i,j,k+1));
                 *
                       * const double deriv_V_im = derivee_aniso(alpha, unsalpha, unsdpdm, Imag_TF[V_FOURIER](i,j,k-1), Imag_TF[V_FOURIER](i,j,k), Imag_TF[V_FOURIER](i,j,k+1));
                       * const double deriv_W_re = (Reel_TF[W_FOURIER](i,j,k+1) - Reel_TF[W_FOURIER](i,j,k)) * unsurdz;
                 *
                       * const double deriv_W_im = (Imag_TF[W_FOURIER](i,j,k+1) - Imag_TF[W_FOURIER](i,j,k)) * unsurdz;
                 */


                double prod = RUE*RWE
                              + IUE*IWE;
                resultat_[RES_PROD_TURB](i,j,k) = prod;

                double Tin_i = RUE*Reel_TF[TRIADIQUE_INPLANE_I_FOURIER](i,j,k)
                               + IUE*Imag_TF[TRIADIQUE_INPLANE_I_FOURIER](i,j,k);
                double Tin_j = RVE*Reel_TF[TRIADIQUE_INPLANE_J_FOURIER](i,j,k)
                               + IVE*Imag_TF[TRIADIQUE_INPLANE_J_FOURIER](i,j,k);
                double Tin_k = RWE*Reel_TF[TRIADIQUE_INPLANE_K_FOURIER](i,j,k)
                               + IWE*Imag_TF[TRIADIQUE_INPLANE_K_FOURIER](i,j,k);
                resultat_[RES_TRIADIQUE_INPLANE](i,j,k) = Tin_i + Tin_j + Tin_k;

                double Tit_i = RUE*Reel_TF[TRIADIQUE_INTERPLANE_I_FOURIER](i,j,k)
                               + IUE*Imag_TF[TRIADIQUE_INTERPLANE_I_FOURIER](i,j,k);
                double Tit_j = RVE*Reel_TF[TRIADIQUE_INTERPLANE_J_FOURIER](i,j,k)
                               + IVE*Imag_TF[TRIADIQUE_INTERPLANE_J_FOURIER](i,j,k);
                double Tit_k = RWE*Reel_TF[TRIADIQUE_INTERPLANE_K_FOURIER](i,j,k)
                               + IWE*Imag_TF[TRIADIQUE_INTERPLANE_K_FOURIER](i,j,k);
                resultat_[RES_TRIADIQUE_INTERPLANE](i,j,k) = Tit_i + Tit_j + Tit_k;

                double Tdef_i = RUE*Reel_TF[DEFORMATION_TURBULENTE_I_FOURIER](i,j,k)
                                + IUE*Imag_TF[DEFORMATION_TURBULENTE_I_FOURIER](i,j,k);
                double Tdef_j = RVE*Reel_TF[DEFORMATION_TURBULENTE_J_FOURIER](i,j,k)
                                + IVE*Imag_TF[DEFORMATION_TURBULENTE_J_FOURIER](i,j,k);
                double Tdef_k = RWE*Reel_TF[DEFORMATION_TURBULENTE_K_FOURIER](i,j,k)
                                + IWE*Imag_TF[DEFORMATION_TURBULENTE_K_FOURIER](i,j,k);
                resultat_[RES_DEFORMATION_TURBULENTE](i,j,k) = Tdef_i + Tdef_j + Tdef_k;

                double Pdiv = Reel_TF[DIVU_FOURIER](i,j,k)*Reel_TF[PSURRHO_FOURIER](i,j,k)
                              + Imag_TF[DIVU_FOURIER](i,j,k)*Imag_TF[PSURRHO_FOURIER](i,j,k);
                resultat_[RES_CORRELATION_PRESSION_DIVERGENCE](i,j,k) = Pdiv;

                /*
                 * double Pdiv_ver2 = ( - veconde_kx*IUE - veconde_ky*IVE + deriv_W_re)*Reel_TF[PSURRHO_FOURIER](i,j,k)
                *                  + (veconde_kx*RUE + veconde_ky*RVE + deriv_W_im)*Imag_TF[PSURRHO_FOURIER](i,j,k);
                 * resultat_[RES_CORRELATION_PRESSION_DIVERGENCE_VER2](i,j,k) = Pdiv_ver2;
                */

                double wp = RWE*Reel_TF[P_FOURIER](i,j,k) + IWE*Imag_TF[P_FOURIER](i,j,k);
                resultat_[RES_P_WFLUC](i,j,k) = wp;

                double Pflu = RWE*Reel_TF[DIFFUSION_PRESSION_RHOFLUC_ADERIVE_FOURIER](i,j,k)
                              + IWE*Imag_TF[DIFFUSION_PRESSION_RHOFLUC_ADERIVE_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_PRESSION_RHOFLUC_ADERIVE](i,j,k) = Pflu;

                double Pdrho_i = RUE*Reel_TF[DIFFUSION_PRESSION_DERIVRHO_I_FOURIER](i,j,k)
                                 + IUE*Imag_TF[DIFFUSION_PRESSION_DERIVRHO_I_FOURIER](i,j,k);
                double Pdrho_j = RVE*Reel_TF[DIFFUSION_PRESSION_DERIVRHO_J_FOURIER](i,j,k)
                                 + IVE*Imag_TF[DIFFUSION_PRESSION_DERIVRHO_J_FOURIER](i,j,k);
                double Pdrho_k = RWE*Reel_TF[DIFFUSION_PRESSION_DERIVRHO_K_FOURIER](i,j,k)
                                 + IWE*Imag_TF[DIFFUSION_PRESSION_DERIVRHO_K_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_PRESSION_DERIVRHO](i,j,k) = Pdrho_i + Pdrho_j + Pdrho_k;

                double Dincun_ii = Reel_TF[DUIDXJ_II_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_II_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_II_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_II_NONCENTRE_FOURIER](i,j,k);
                double Dincun_ij = Reel_TF[DUIDXJ_IJ_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_IJ_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_IJ_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_IJ_NONCENTRE_FOURIER](i,j,k);
                double Dincun_ik = Reel_TF[DUIDXJ_IK_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_IK_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_IK_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_IK_NONCENTRE_FOURIER](i,j,k);
                double Dincun_ji = Reel_TF[DUIDXJ_JI_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_JI_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_JI_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_JI_NONCENTRE_FOURIER](i,j,k);
                double Dincun_jj = Reel_TF[DUIDXJ_JJ_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_JJ_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_JJ_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_JJ_NONCENTRE_FOURIER](i,j,k);
                double Dincun_jk = Reel_TF[DUIDXJ_JK_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_JK_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_JK_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_JK_NONCENTRE_FOURIER](i,j,k);
                double Dincun_ki = Reel_TF[DUIDXJ_KI_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_KI_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_KI_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_KI_NONCENTRE_FOURIER](i,j,k);
                double Dincun_kj = Reel_TF[DUIDXJ_KJ_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_KJ_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_KJ_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_KJ_NONCENTRE_FOURIER](i,j,k);
                double Dincun_kk = Reel_TF[DUIDXJ_KK_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DUIDXJ_KK_NONCENTRE_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_KK_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DUIDXJ_KK_NONCENTRE_FOURIER](i,j,k);
                resultat_[RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INPLANE](i,j,k) = Dincun_ii + Dincun_ij + Dincun_ji + Dincun_jj + Dincun_ki + Dincun_kj;
                resultat_[RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INTERPLANE](i,j,k) = Dincun_ik + Dincun_jk + Dincun_kk;

                /*
                 * double Dincun_ver2_ii = veconde_kx*veconde_kx*(uu);
                 * double Dincun_ver2_ij = veconde_ky*veconde_ky*(uu);
                 * double Dincun_ver2_ik = deriv_U_re*deriv_U_re + deriv_U_im*deriv_U_im;
                 * double Dincun_ver2_ji = veconde_kx*veconde_kx*(vv);
                 * double Dincun_ver2_jj = veconde_ky*veconde_ky*(vv);
                 * double Dincun_ver2_jk = deriv_V_re*deriv_V_re + deriv_V_im*deriv_V_im;
                 * double Dincun_ver2_ki = veconde_kx*veconde_kx*(ww);
                 * double Dincun_ver2_kj = veconde_ky*veconde_ky*(ww);
                 * double Dincun_ver2_kk = deriv_W_re*deriv_W_re + deriv_W_im*deriv_W_im;
                 * resultat_[RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INPLANE_VER2](i,j,k) = Dincun_ver2_ii + Dincun_ver2_ij + Dincun_ver2_ji + Dincun_ver2_jj + Dincun_ver2_ki + Dincun_ver2_kj;
                 * resultat_[RES_DISSIPATION_INCOMPRESSIBLE_UN_AFOISNU_INTERPLANE_VER2](i,j,k) = Dincun_ver2_ik + Dincun_ver2_jk + Dincun_ver2_kk;
                */

                double Dincdeux_ii = Reel_TF[DUIDXJ_II_FOURIER](i,j,k)*Reel_TF[DUIDXJ_II_FOURIER](i,j,k)
                                     + Imag_TF[DUIDXJ_II_FOURIER](i,j,k)*Imag_TF[DUIDXJ_II_FOURIER](i,j,k);
                double Dincdeux_ij = Reel_TF[DUIDXJ_IJ_FOURIER](i,j,k)*Reel_TF[DUIDXJ_JI_FOURIER](i,j,k)
                                     + Imag_TF[DUIDXJ_IJ_FOURIER](i,j,k)*Imag_TF[DUIDXJ_JI_FOURIER](i,j,k);
                double Dincdeux_ik = Reel_TF[DUIDXJ_IK_FOURIER](i,j,k)*Reel_TF[DUIDXJ_KI_FOURIER](i,j,k)
                                     + Imag_TF[DUIDXJ_IK_FOURIER](i,j,k)*Imag_TF[DUIDXJ_KI_FOURIER](i,j,k);
                double Dincdeux_jj = Reel_TF[DUIDXJ_JJ_FOURIER](i,j,k)*Reel_TF[DUIDXJ_JJ_FOURIER](i,j,k)
                                     + Imag_TF[DUIDXJ_JJ_FOURIER](i,j,k)*Imag_TF[DUIDXJ_JJ_FOURIER](i,j,k);
                double Dincdeux_jk = Reel_TF[DUIDXJ_JK_FOURIER](i,j,k)*Reel_TF[DUIDXJ_KJ_FOURIER](i,j,k)
                                     + Imag_TF[DUIDXJ_JK_FOURIER](i,j,k)*Imag_TF[DUIDXJ_KJ_FOURIER](i,j,k);
                double Dincdeux_kk = Reel_TF[DUIDXJ_KK_FOURIER](i,j,k)*Reel_TF[DUIDXJ_KK_FOURIER](i,j,k)
                                     + Imag_TF[DUIDXJ_KK_FOURIER](i,j,k)*Imag_TF[DUIDXJ_KK_FOURIER](i,j,k);
                resultat_[RES_DISSIPATION_INCOMPRESSIBLE_DEUX_AFOISNU](i,j,k) = Dincdeux_ii + 2*Dincdeux_ij + 2*Dincdeux_ik + Dincdeux_jj + 2*Dincdeux_jk + Dincdeux_kk;

                /*
                 * double Dincdeux_ver2_ii = veconde_kx*veconde_kx*(uu);
                 * double Dincdeux_ver2_ij = veconde_kx*veconde_ky*(uu);
                 * double Dincdeux_ver2_ik = veconde_kx*(deriv_U_im*RWE
                 *                                     - deriv_U_re*IWE);
                 * double Dincdeux_ver2_jj = veconde_ky*veconde_ky*(vv);
                 * double Dincdeux_ver2_jk = veconde_ky*(deriv_V_im*RWE
                 *                                     - deriv_V_re*IWE);
                 * double Dincdeux_ver2_kk = deriv_W_re*deriv_W_re + deriv_W_im*deriv_W_im;
                 * resultat_[RES_DISSIPATION_INCOMPRESSIBLE_DEUX_AFOISNU_VER2](i,j,k) = Dincdeux_ver2_ii + 2*Dincdeux_ver2_ij + 2*Dincdeux_ver2_ik + Dincdeux_ver2_jj + 2*Dincdeux_ver2_jk + Dincdeux_ver2_kk;
                */

                double Dcom_un_ii = Reel_TF[DUIDXJ_II_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_II_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_II_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_II_FOURIER](i,j,k);
                double Dcom_un_ij = Reel_TF[DUIDXJ_IJ_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_IJ_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_IJ_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_IJ_FOURIER](i,j,k);
                double Dcom_un_ik = Reel_TF[DUIDXJ_IK_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_IK_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_IK_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_IK_FOURIER](i,j,k);
                double Dcom_un_ji = Reel_TF[DUIDXJ_JI_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_JI_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_JI_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_JI_FOURIER](i,j,k);
                double Dcom_un_jj = Reel_TF[DUIDXJ_JJ_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_JJ_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_JJ_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_JJ_FOURIER](i,j,k);
                double Dcom_un_jk = Reel_TF[DUIDXJ_JK_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_JK_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_JK_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_JK_FOURIER](i,j,k);
                double Dcom_un_ki = Reel_TF[DUIDXJ_KI_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_KI_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_KI_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_KI_FOURIER](i,j,k);
                double Dcom_un_kj = Reel_TF[DUIDXJ_KJ_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_KJ_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_KJ_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_KJ_FOURIER](i,j,k);
                double Dcom_un_kk = Reel_TF[DUIDXJ_KK_NONCENTRE_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_KK_FOURIER](i,j,k)
                                    + Imag_TF[DUIDXJ_KK_NONCENTRE_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_KK_FOURIER](i,j,k);
                resultat_[RES_DISSIPATION_COMPRESSIBLE_UN](i,j,k) = Dcom_un_ii + Dcom_un_ij + Dcom_un_ik + Dcom_un_ji + Dcom_un_jj + Dcom_un_jk + Dcom_un_ki + Dcom_un_kj + Dcom_un_kk;

                /*
                       * double Dcom_un_ver2_ii = veconde_kx*(RUE*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_II_FOURIER](i,j,k)
                       *                                    - IUE*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_II_FOURIER](i,j,k));
                       * double Dcom_un_ver2_ij = veconde_ky*(RUE*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_IJ_FOURIER](i,j,k)
                       *                                    - IUE*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_IJ_FOURIER](i,j,k));
                       * double Dcom_un_ver2_ik = deriv_U_re*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_IK_FOURIER](i,j,k)
                       *                        + deriv_U_im*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_IK_FOURIER](i,j,k);
                       * double Dcom_un_ver2_ji = veconde_kx*(RVE*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_JI_FOURIER](i,j,k)
                       *                                    - IVE*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_JI_FOURIER](i,j,k));
                       * double Dcom_un_ver2_jj = veconde_ky*(RVE*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_JJ_FOURIER](i,j,k)
                       *                                    - IVE*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_JJ_FOURIER](i,j,k));
                       * double Dcom_un_ver2_jk = deriv_V_re*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_JK_FOURIER](i,j,k)
                       *                        + deriv_V_im*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_JK_FOURIER](i,j,k);
                       * double Dcom_un_ver2_ki = veconde_kx*(RWE*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_KI_FOURIER](i,j,k)
                       *                                    - IWE*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_KI_FOURIER](i,j,k));
                       * double Dcom_un_ver2_kj = veconde_ky*(RWE*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_KJ_FOURIER](i,j,k)
                       *                                    - IWE*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_KJ_FOURIER](i,j,k));
                       * double Dcom_un_ver2_kk = deriv_W_re*Reel_TF[DISSIPATION_COMPRESSIBLE_UN_KK_FOURIER](i,j,k)
                       *                        + deriv_W_im*Imag_TF[DISSIPATION_COMPRESSIBLE_UN_KK_FOURIER](i,j,k);
                       * resultat_[RES_DISSIPATION_COMPRESSIBLE_UN_VER2](i,j,k) = Dcom_un_ver2_ii + Dcom_un_ver2_ij + Dcom_un_ver2_ik + Dcom_un_ver2_ji + Dcom_un_ver2_jj + Dcom_un_ver2_jk + Dcom_un_ver2_ki + Dcom_un_ver2_kj + Dcom_un_ver2_kk;
                 */

                double Dcom_deux_ii = Reel_TF[DUIDXJ_II_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_II_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_II_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_II_FOURIER](i,j,k);
                double Dcom_deux_ij = Reel_TF[DUIDXJ_IJ_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_IJ_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_IJ_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_IJ_FOURIER](i,j,k);
                double Dcom_deux_ik = Reel_TF[DUIDXJ_IK_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_IK_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_IK_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_IK_FOURIER](i,j,k);
                double Dcom_deux_ji = Reel_TF[DUIDXJ_JI_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_JI_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_JI_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_JI_FOURIER](i,j,k);
                double Dcom_deux_jj = Reel_TF[DUIDXJ_JJ_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_JJ_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_JJ_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_JJ_FOURIER](i,j,k);
                double Dcom_deux_jk = Reel_TF[DUIDXJ_JK_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_JK_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_JK_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_JK_FOURIER](i,j,k);
                double Dcom_deux_ki = Reel_TF[DUIDXJ_KI_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_KI_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_KI_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_KI_FOURIER](i,j,k);
                double Dcom_deux_kj = Reel_TF[DUIDXJ_KJ_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_KJ_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_KJ_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_KJ_FOURIER](i,j,k);
                double Dcom_deux_kk = Reel_TF[DUIDXJ_KK_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_KK_FOURIER](i,j,k)
                                      + Imag_TF[DUIDXJ_KK_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_KK_FOURIER](i,j,k);
                resultat_[RES_DISSIPATION_COMPRESSIBLE_DEUX](i,j,k) = Dcom_deux_ii + Dcom_deux_ij + Dcom_deux_ik + Dcom_deux_ji + Dcom_deux_jj + Dcom_deux_jk + Dcom_deux_ki + Dcom_deux_kj + Dcom_deux_kk;

                /*
                       * double Dcom_deux_ver2_ii = veconde_kx*(RUE*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_II_FOURIER](i,j,k)
                       *                                      - IUE*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_II_FOURIER](i,j,k));
                       * double Dcom_deux_ver2_ij = veconde_ky*(RUE*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_IJ_FOURIER](i,j,k)
                       *                                      - IUE*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_IJ_FOURIER](i,j,k));
                       * double Dcom_deux_ver2_ik = deriv_U_re*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_IK_FOURIER](i,j,k)
                       *                          + deriv_U_im*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_IK_FOURIER](i,j,k);
                       * double Dcom_deux_ver2_ji = veconde_kx*(RVE*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_JI_FOURIER](i,j,k)
                       *                                      - IVE*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_JI_FOURIER](i,j,k));
                       * double Dcom_deux_ver2_jj = veconde_ky*(RVE*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_JJ_FOURIER](i,j,k)
                       *                                      - IVE*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_JJ_FOURIER](i,j,k));
                       * double Dcom_deux_ver2_jk = deriv_V_re*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_JK_FOURIER](i,j,k)
                       *                          + deriv_V_im*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_JK_FOURIER](i,j,k);
                       * double Dcom_deux_ver2_ki = veconde_kx*(RWE*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_KI_FOURIER](i,j,k)
                       *                                      - IWE*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_KI_FOURIER](i,j,k));
                       * double Dcom_deux_ver2_kj = veconde_ky*(RWE*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_KJ_FOURIER](i,j,k)
                       *                                      - IWE*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_KJ_FOURIER](i,j,k));
                       * double Dcom_deux_ver2_kk = deriv_W_re*Reel_TF[DISSIPATION_COMPRESSIBLE_DEUX_KK_FOURIER](i,j,k)
                       *                          + deriv_W_im*Imag_TF[DISSIPATION_COMPRESSIBLE_DEUX_KK_FOURIER](i,j,k);
                       * resultat_[RES_DISSIPATION_COMPRESSIBLE_DEUX_VER2](i,j,k) = Dcom_deux_ver2_ii + Dcom_deux_ver2_ij + Dcom_deux_ver2_ik + Dcom_deux_ver2_ji + Dcom_deux_ver2_jj + Dcom_deux_ver2_jk + Dcom_deux_ver2_ki + Dcom_deux_ver2_kj + Dcom_deux_ver2_kk;
                 */

                double Dcom_trois = Reel_TF[DIVU_FOURIER](i,j,k)*Reel_TF[DISSIPATION_COMPRESSIBLE_TROIS_I_FOURIER](i,j,k)
                                    + Imag_TF[DIVU_FOURIER](i,j,k)*Imag_TF[DISSIPATION_COMPRESSIBLE_TROIS_I_FOURIER](i,j,k);
                resultat_[RES_DISSIPATION_COMPRESSIBLE_TROIS](i,j,k) = Dcom_trois;

                /*
                       * double Dcom_trois_ver2 = ( - veconde_kx*IUE - veconde_ky*IVE + deriv_W_re)*Reel_TF[DISSIPATION_COMPRESSIBLE_TROIS_I_FOURIER](i,j,k)
                       *                        + (veconde_kx*RUE + veconde_ky*RVE + deriv_W_im)*Imag_TF[DISSIPATION_COMPRESSIBLE_TROIS_I_FOURIER](i,j,k);
                       * resultat_[RES_DISSIPATION_COMPRESSIBLE_TROIS_VER2](i,j,k) = Dcom_trois_ver2;
                 */

                double Vincdeux_i = RUE*Reel_TF[DUIDXJ_KI_FOURIER](i,j,k)
                                    + IUE*Imag_TF[DUIDXJ_KI_FOURIER](i,j,k);
                double Vincdeux_j = RVE*Reel_TF[DUIDXJ_KJ_FOURIER](i,j,k)
                                    + IVE*Imag_TF[DUIDXJ_KJ_FOURIER](i,j,k);
                double Vincdeux_k = RWE*Reel_TF[DUIDXJ_KK_FOURIER](i,j,k)
                                    + IWE*Imag_TF[DUIDXJ_KK_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_VISQUEUSE_INCOMPRESSIBLE_DEUX](i,j,k) = Vincdeux_i + Vincdeux_j + Vincdeux_k;

                double Vflu_un_i = RUE*Reel_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_I_FOURIER](i,j,k)
                                   + IUE*Imag_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_I_FOURIER](i,j,k);
                double Vflu_un_j = RVE*Reel_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_J_FOURIER](i,j,k)
                                   + IVE*Imag_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_J_FOURIER](i,j,k);
                double Vflu_un_k = RWE*Reel_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_K_FOURIER](i,j,k)
                                   + IWE*Imag_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_K_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN](i,j,k) = Vflu_un_i + Vflu_un_j + Vflu_un_k;

                double Vflu_deux_i = RUE*Reel_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_I_FOURIER](i,j,k)
                                     + IUE*Imag_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_I_FOURIER](i,j,k);
                double Vflu_deux_j = RVE*Reel_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_J_FOURIER](i,j,k)
                                     + IVE*Imag_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_J_FOURIER](i,j,k);
                double Vflu_deux_k = RWE*Reel_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_K_FOURIER](i,j,k)
                                     + IWE*Imag_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_K_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX](i,j,k) = Vflu_deux_i + Vflu_deux_j + Vflu_deux_k;

                double Vflu_trois = RWE*Reel_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS_FOURIER](i,j,k)
                                    + IWE*Imag_TF[DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS](i,j,k) = Vflu_trois;

                double Vdrho_un_i = RUE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_UN_I_FOURIER](i,j,k)
                                    + IUE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_UN_I_FOURIER](i,j,k);
                double Vdrho_un_j = RVE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_UN_J_FOURIER](i,j,k)
                                    + IVE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_UN_J_FOURIER](i,j,k);
                double Vdrho_un_k = RWE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_UN_K_FOURIER](i,j,k)
                                    + IWE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_UN_K_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_VISQUEUSE_DERIVRHO_UN](i,j,k) = Vdrho_un_i + Vdrho_un_j + Vdrho_un_k;

                double Vdrho_deux_i = RUE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_I_FOURIER](i,j,k)
                                      + IUE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_I_FOURIER](i,j,k);
                double Vdrho_deux_j = RVE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_J_FOURIER](i,j,k)
                                      + IVE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_J_FOURIER](i,j,k);
                double Vdrho_deux_k = RWE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_K_FOURIER](i,j,k)
                                      + IWE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_K_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_VISQUEUSE_DERIVRHO_DEUX](i,j,k) = Vdrho_deux_i + Vdrho_deux_j + Vdrho_deux_k;

                double Vdrho_trois_i = RUE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_I_FOURIER](i,j,k)
                                       + IUE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_I_FOURIER](i,j,k);
                double Vdrho_trois_j = RVE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_J_FOURIER](i,j,k)
                                       + IVE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_J_FOURIER](i,j,k);
                double Vdrho_trois_k = RWE*Reel_TF[DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_K_FOURIER](i,j,k)
                                       + IWE*Imag_TF[DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_K_FOURIER](i,j,k);
                resultat_[RES_DIFFUSION_VISQUEUSE_DERIVRHO_TROIS](i,j,k) = Vdrho_trois_i + Vdrho_trois_j + Vdrho_trois_k;

                double TNC_i = RUE*Reel_TF[TERME_NON_CLASSE_I_FOURIER](i,j,k)
                               + IUE*Imag_TF[TERME_NON_CLASSE_I_FOURIER](i,j,k);
                double TNC_j = RVE*Reel_TF[TERME_NON_CLASSE_J_FOURIER](i,j,k)
                               + IVE*Imag_TF[TERME_NON_CLASSE_J_FOURIER](i,j,k);
                double TNC_k = RWE*Reel_TF[TERME_NON_CLASSE_K_FOURIER](i,j,k)
                               + IWE*Imag_TF[TERME_NON_CLASSE_K_FOURIER](i,j,k);
                resultat_[RES_TERME_NON_CLASSE](i,j,k) = TNC_i + TNC_j + TNC_k;

                double CT_i = RUE*Reel_TF[UIUJ_IK_FOURIER](i,j,k)
                              + IUE*Imag_TF[UIUJ_IK_FOURIER](i,j,k);
                double CT_j = RVE*Reel_TF[UIUJ_JK_FOURIER](i,j,k)
                              + IVE*Imag_TF[UIUJ_JK_FOURIER](i,j,k);
                double CT_k = RWE*Reel_TF[UIUJ_KK_FOURIER](i,j,k)
                              + IWE*Imag_TF[UIUJ_KK_FOURIER](i,j,k);
                resultat_[RES_CONVECTIONTURBULENTE](i,j,k) = - 0.5 * (CT_i + CT_j + CT_k);

                double TPS_un_ii = Reel_TF[DUIDXJ_II_FOURIER](i,j,k)*Reel_TF[UIUJ_II_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_II_FOURIER](i,j,k)*Imag_TF[UIUJ_II_FOURIER](i,j,k);
                double TPS_un_ij = Reel_TF[DUIDXJ_IJ_FOURIER](i,j,k)*Reel_TF[UIUJ_IJ_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_IJ_FOURIER](i,j,k)*Imag_TF[UIUJ_IJ_FOURIER](i,j,k);
                double TPS_un_ik = Reel_TF[DUIDXJ_IK_FOURIER](i,j,k)*Reel_TF[UIUJ_IK_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_IK_FOURIER](i,j,k)*Imag_TF[UIUJ_IK_FOURIER](i,j,k);
                double TPS_un_ji = Reel_TF[DUIDXJ_JI_FOURIER](i,j,k)*Reel_TF[UIUJ_IJ_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_JI_FOURIER](i,j,k)*Imag_TF[UIUJ_IJ_FOURIER](i,j,k);
                double TPS_un_jj = Reel_TF[DUIDXJ_JJ_FOURIER](i,j,k)*Reel_TF[UIUJ_JJ_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_JJ_FOURIER](i,j,k)*Imag_TF[UIUJ_JJ_FOURIER](i,j,k);
                double TPS_un_jk = Reel_TF[DUIDXJ_JK_FOURIER](i,j,k)*Reel_TF[UIUJ_JK_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_JK_FOURIER](i,j,k)*Imag_TF[UIUJ_JK_FOURIER](i,j,k);
                double TPS_un_ki = Reel_TF[DUIDXJ_KI_FOURIER](i,j,k)*Reel_TF[UIUJ_IK_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_KI_FOURIER](i,j,k)*Imag_TF[UIUJ_IK_FOURIER](i,j,k);
                double TPS_un_kj = Reel_TF[DUIDXJ_KJ_FOURIER](i,j,k)*Reel_TF[UIUJ_JK_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_KJ_FOURIER](i,j,k)*Imag_TF[UIUJ_JK_FOURIER](i,j,k);
                double TPS_un_kk = Reel_TF[DUIDXJ_KK_FOURIER](i,j,k)*Reel_TF[UIUJ_KK_FOURIER](i,j,k)
                                   + Imag_TF[DUIDXJ_KK_FOURIER](i,j,k)*Imag_TF[UIUJ_KK_FOURIER](i,j,k);
                double TPS_deux_i = RUE*Reel_TF[UJDUIDXJ_I_FOURIER](i,j,k)
                                    + IUE*Imag_TF[UJDUIDXJ_I_FOURIER](i,j,k);
                double TPS_deux_j = RVE*Reel_TF[UJDUIDXJ_J_FOURIER](i,j,k)
                                    + IVE*Imag_TF[UJDUIDXJ_J_FOURIER](i,j,k);
                double TPS_deux_k = RWE*Reel_TF[UJDUIDXJ_K_FOURIER](i,j,k)
                                    + IWE*Imag_TF[UJDUIDXJ_K_FOURIER](i,j,k);
                double TPS_un = 0.5 * (TPS_un_ii + TPS_un_ij + TPS_un_ik + TPS_un_ji + TPS_un_jj + TPS_un_jk + TPS_un_ki + TPS_un_kj + TPS_un_kk);
                double TPS_deux = - 0.5 * (TPS_deux_i + TPS_deux_j + TPS_deux_k);
                resultat_[RES_TERMEPUREMENTSPECTRAL](i,j,k) = TPS_un + TPS_deux;
              }
        }
    }
  const int N = N_RES;
  const double normalisation = (2./ (imax*jmax))*(2./ (imax*jmax));
  for ( int val = 0 ; val < N ; val ++ )
    for (int k = 0; k < kmax; k++)
      for (int j = 0; j < jmax; j++)
        for (int i = 0; i < imax; i++)
          resultat_[val](i,j,k) *= normalisation ;
  statistiques().end_count(TF_Combinaisons);

#define DEBUG 0
  if ( DEBUG )
    {
      static Stat_Counter_Id TF_lata = statistiques().new_counter(2, "TF lata");
      statistiques().begin_count(TF_lata);
      Cerr << " Attention mode debug dans le post-traitement spectral dump de " << N_RES+Nval*2+1 << " lata par dt " << finl;
      /* ici au besoin pour peut mettre un dump de lata pour avoir l'info complete en cas de debug !!!! */
      Nom lata_name("DEBUG_spectre.lata");
      dumplata_header(lata_name, resultat_[0] /* on passe un champ pour ecrire la geometrie */);
      dumplata_newtime(lata_name,N_echantillons_); // c'est moche ca

      for ( int val = 0 ; val < N ; val ++ )
        {
          recentrer_spectres(resultat_[val],Avant_TF[0]);
          dumplata_scalar(lata_name,Nom(val), resultat_[val], 0);
        }
      for (int val = 0 ; val < Nval ; val ++)
        {
          Nom NOMREEL = "PARTIE_REELLE_";
          NOMREEL+= Nom(val);
          Nom NOMIMAG = "PARTIE_IMAGINAIRE_";
          NOMIMAG+= Nom(val);
          recentrer_spectres(Reel_TF[val],Avant_TF[0]);
          recentrer_spectres(Imag_TF[val],Avant_TF[0]);
          dumplata_scalar(lata_name,NOMREEL, Reel_TF[val], 0);
          dumplata_scalar(lata_name,NOMIMAG, Imag_TF[val], 0);
        }
      dumplata_vector(lata_name,"REVIT", vitesse[0],vitesse[1],vitesse[2], 0);
      statistiques().end_count(TF_lata);
    }
  /* toto */

  /* On a calculer les termes, on va maintenant stocker + post-traiter les TFx et TFz
    * Pour ceci on utilise deux tableaux 2D contenant toute l'information pour un dump direct;
    */

  static Stat_Counter_Id TF_tfx_tfy = statistiques().new_counter(2, "TF partie TFx, TFy");
  statistiques().begin_count(TF_tfx_tfy);
  for ( int val = 0 ; val < N ; val ++ )
    for (int k = 0; k < kmax; k++)
      {
        for (int i = 0; i < imax/2; i++)
          {
            // Bon on fait un cas de test ;;
            const int Indice = i + imax/2*(k) ;
            const double x = resultat_[val](i,0,k);
            TFx_[val][Indice]+= x;
          }
        for (int j = 0; j < jmax/2; j++)
          {
            // Bon on fait un cas de test ;;
            const int Indice = j + jmax/2*(k);
            const double y = resultat_[val](0,j,k);
            TFy_[val][Indice]+= y;
          }
      }
  statistiques().end_count(TF_tfx_tfy);

  static Stat_Counter_Id TF_integrale = statistiques().new_counter(2, "TF partie integrale");
  statistiques().begin_count(TF_integrale);

  const int taille = tri_ki.size_array(); // donne la taille du tableaux
  for ( int val = 0 ; val < N ; val ++ )
    for (int k = 0 ; k < kmax ; k++)
      {
        const int decalage = taille*(k);
        for (int j = 0 ; j < jmax ; j++)
          for (int i = 0 ; i < imax ; i++)
            {
              const int Indice = ref_ki(i,j) + decalage;
              const double x = resultat_[val](i,j,k);
              TFi_[val][Indice]+= x/2.;
            }
      }
  statistiques().end_count(TF_integrale);

  // a ne pas oublier
  N_echantillons_++ ;
}


/* converti un tableau dans l'espace spectral */
void Fourier_trans::Traitement_spectral(const IJK_Field_double& entree, IJK_Field_double& reel, IJK_Field_double& imag, int sens)
{
  const int Nk  = Post_splitting_.get_nb_elem_local(DIRECTION_K);
  if ( sens == 1 )
    {
      for (int k = 0 ; k < Nk ; k++)
        Traitement_spectral_klayer_direct(entree,reel,imag,k);
    }
  else if (sens == -1 )
    {
      for (int k = 0 ; k < Nk ; k++)
        Traitement_spectral_klayer_inverse(reel,imag,k); //* attention ici les donnees font etre ecrasee
    }
  else
    {
      Cerr << "parametre incompatible dans Traitement_spectral " << finl;
      Process::exit();
    }
}

void Fourier_trans::Traitement_spectral_klayer_direct(const IJK_Field_double& entree, IJK_Field_double& reel, IJK_Field_double& imag, int k)
{
  const int Ni  = Post_splitting_.get_nb_elem_local(DIRECTION_I);
  const int Nj  = Post_splitting_.get_nb_elem_local(DIRECTION_J);

  for (int j=0 ; j < Nj ; j++)
    for (int i=0; i < Ni ; i++)
      {
        const double x = entree(i,j,k);
        in[i+j*Ni][0]= x;
        in[i+j*Ni][1] = 0;
        out[i+j*Ni][0]= 0;
        out[i+j*Ni][1]= 0;
      }

  fftw_execute(plan); // Calcule la TF
  /* EXPLICATION : le retour de la TF est un plan 2D de Ni,Nj points ainsi, il se stock naturellement dans un tableau ijk !!!
   * pour se balader le long des k1 on se deplace sur les i pour les k2 on se deplacera sur j */

  for (int j=0 ; j < Nj ; j++)
    for (int i=0; i < Ni ; i++)
      {
        const double x = out[i+j*Ni][0]; // on fait ceci pour aider le compilateur
        const double y = out[i+j*Ni][1];
        reel(i,j,k)= x;
        imag(i,j,k)= y;
      }
}

void Fourier_trans::Traitement_spectral_klayer_inverse(IJK_Field_double& reel, IJK_Field_double& imag, int k)
{
  const int Ni  = Post_splitting_.get_nb_elem_local(DIRECTION_I);
  const int Nj  = Post_splitting_.get_nb_elem_local(DIRECTION_J);

  for (int j=0 ; j < Nj ; j++)
    for (int i=0; i < Ni ; i++)
      {
        const double x = reel(i,j,k);
        const double y = imag(i,j,k);
        in[i+j*Ni][0] = x;
        in[i+j*Ni][1] = y;
        out[i+j*Ni][0]= 0;
        out[i+j*Ni][1]= 0;
      }

  fftw_execute(inverseplan); // Calcule la TF
  /* EXPLICATION : le retour de la TF est un plan 2D de Ni,Nj points ainsi, il se stock naturellement dans un tableau ijk !!!
   * pour se balader le long des k1 on se deplace sur les i pour les k2 on se deplacera sur j */

  for (int j=0 ; j < Nj ; j++)
    for (int i=0; i < Ni ; i++)
      {
        const double x = out[i+j*Ni][0]; // on fait ceci pour aider le compilateur
        const double y = out[i+j*Ni][1];
        reel(i,j,k)= x;
        imag(i,j,k)= y;
      }
}

void Fourier_trans::postraiter(Nom& base) const
{
  Nom Nom_TFi = base + "_TFi_";
  Nom_TFi += Nom(me());
  Nom_TFi += ".txt" ;

  SFichier fi(Nom_TFi);
  // F.A modification de la precision pour allez chercher les 4 ordres
  fi.setf(ios::scientific);
  fi.precision(15);
  postraiter_TFi(fi,0);

  Nom Nom_TFx  = base + "_TFx_";
  Nom_TFx += Nom(me());
  Nom_TFx += ".txt" ; ;

  SFichier fx(Nom_TFx);
  // F.A modification de la precision pour allez chercher les 4 ordres
  fx.setf(ios::scientific);
  fx.precision(15);
  postraiter_TFx(fx,0);

  Nom Nom_TFy  = base + "_TFy_";
  Nom_TFy += Nom(me());
  Nom_TFy += ".txt" ;

  SFichier fy(Nom_TFy);
  // F.A modification de la precision pour allez chercher les 4 ordres
  fy.setf(ios::scientific);
  fy.precision(15);
  postraiter_TFy(fy,0);

  Nom Entete = base + "_Entete.txt";
  if (je_suis_maitre())
    {
      SFichier fe(Entete);
      Ecrire_entete(fe,0);
    }
}

void Fourier_trans::Ecrire_entete(Sortie& os, int flag_valeur_instantanee) const
{
  const int N = N_RES;
  if (flag_valeur_instantanee == 0)
    {
      os << "# N_echantillons " <<  N_echantillons_ << finl;
      os << "# Impression des moyennes temporelles" << finl;
    }
  else if (flag_valeur_instantanee == 1)
    {
      os << "# Impression des moyennes spatiales instantanee" << finl;
    }
  else
    {
      Cerr << "Erreur dans Fourier_trans::postraiter: flag inconnu" << finl;
      Process::exit();
    }
  os << "# colonne 1 : coordonnee z" << finl;
  os << "# colonne 2 : vecteur d'onde k" << finl;
  for (int i = 0; i < N; i++)
    {
      os << "# colonne " << i+3 << " : " << noms_moyennes[i] << finl;
    }
}

// If flag_valeur_instantanee==1 : write the last computed instantaneous value,
// if flag_valeur_instantanee==0 : write the time averaged value
void Fourier_trans::postraiter_TFi(Sortie& os, int flag_valeur_instantanee) const
{
  const int N = N_RES;
  // tous les processeurs ecrivent en meme temps mais chaun son fichier !!!

  const int kmax  = Post_splitting_.get_nb_elem_local(DIRECTION_K);
  const int offsetk = Post_splitting_.get_offset_local(DIRECTION_K);
  const int Taille = tri_ki.size_array();

  for (int z = 0; z < kmax; z++)
    {
      for (int k =0; k < Taille ; k++ )
        {
          char s[100];
          sprintf(s, "%16.16e ", elem_coord_[z+offsetk]);
          os << s;
          sprintf(s, "%16.16e ", tri_ki[k]);
          os << s;
          for (int val = 0; val < N; val++)
            {
              double x;
              const int Indice = k +z* Taille;
              if (flag_valeur_instantanee == 0)
                x = TFi_[val][Indice] / N_echantillons_;
              else
                x = TFi_[val][Indice];
              sprintf(s, "%16.16e ", x);
              os << s;
            }
          os << finl;
        }
      os << " " << finl;
    }
}

void Fourier_trans::postraiter_TFx(Sortie& os, int flag_valeur_instantanee) const
{
  const int N = N_RES;
  // tous les processeurs ecrivent en meme temps mais chaun son fichier !!!

  const int kmax  = Post_splitting_.get_nb_elem_local(DIRECTION_K);
  const int imax  = Post_splitting_.get_nb_elem_local(DIRECTION_I);
  const int offsetk = Post_splitting_.get_offset_local(DIRECTION_K);
  const IJK_Grid_Geometry& geom = Post_splitting_.get_grid_geometry();
  const int Nx_tot = geom.get_nb_elem_tot(0) + 1;
  const double Lx = geom.get_node_coordinates(0)[Nx_tot - 1] - geom.get_origin(0);
  const double fac = 1./Lx;
  for (int z = 0; z < kmax; z++)
    {
      for (int k =0; k < imax/2 ; k++ )
        {
          char s[100];
          sprintf(s, "%16.16e ", elem_coord_[z+offsetk]);
          os << s;
          sprintf(s, "%16.16e ", k*fac);
          os << s;
          for (int val = 0; val < N; val++)
            {
              double x;
              const int Indice = k +z*imax/2;
              if (flag_valeur_instantanee == 0)
                x = TFx_[val][Indice] / N_echantillons_;
              else
                x = TFx_[val][Indice];
              sprintf(s, "%16.16e ", x);
              os << s;
            }
          os << finl;
        }
      os << " " << finl;
    }
}

void Fourier_trans::postraiter_TFy(Sortie& os, int flag_valeur_instantanee) const
{
  const int N = N_RES;
  // tous les processeurs ecrivent en meme temps mais chaun son fichier !!!

  const int kmax  = Post_splitting_.get_nb_elem_local(DIRECTION_K);
  const int jmax  = Post_splitting_.get_nb_elem_local(DIRECTION_J);
  const int offsetk = Post_splitting_.get_offset_local(DIRECTION_K);
  const IJK_Grid_Geometry& geom = Post_splitting_.get_grid_geometry();
  const int Ny_tot = geom.get_nb_elem_tot(1) + 1;
  const double Ly = geom.get_node_coordinates(1)[Ny_tot - 1] - geom.get_origin(1);
  const double fac = 1./Ly;

  for (int z = 0; z < kmax; z++)
    {
      for (int k =0; k < jmax/2 ; k++ )
        {
          char s[100];
          sprintf(s, "%16.16e ", elem_coord_[z+offsetk]);
          os << s;
          sprintf(s, "%16.16e ", k*fac);
          os << s;
          for (int val = 0; val < N; val++)
            {
              double x;
              const int Indice = k +z*jmax/2;
              if (flag_valeur_instantanee == 0)
                x = TFy_[val][Indice] / N_echantillons_;
              else
                x = TFy_[val][Indice];
              sprintf(s, "%16.16e ", x);
              os << s;
            }
          os << finl;
        }
      os << " " << finl;
    }
}

void Fourier_trans::sauvegarde() const
{
  Cerr << "Sauvegarde parrallele des stats spectrales " << finl;
  Nom nom_sauv = "Sauvegarde_spectrale_";
  // on va indiquer la position dans le mapping, on ne sais jamais des fois que cela bouge !!!!
  int pos = Post_splitting_.get_local_slice_index(2);
  nom_sauv += Nom(pos) + ".sauv";
  SFichier fichier(nom_sauv);
  fichier.precision(17);
  fichier.setf(std::ios_base::scientific);

  const int N = N_RES;
  fichier << "{\n";
  fichier  << " N_ECHANTILLONS " << N_echantillons_ << "\n";
  fichier  << " VECT_K " << tri_ki << "\n";
  fichier  << " REF_KI " << ref_ki << "\n";

  for (True_int i = 0; i < N; i++)
    {
      char nom_champ[100];
      sprintf(nom_champ, "TFI_%d",i );
      fichier << nom_champ << " " << TFi_[i] << finl;
      sprintf(nom_champ, "TFX_%d",i );
      fichier << nom_champ << " " << TFx_[i] << finl;
      sprintf(nom_champ, "TFY_%d",i );
      fichier << nom_champ << " " << TFy_[i] << finl;
    }

  fichier << "}" << finl;
  Cout << "Ecriture des donnees fourier pour reprise: N_echantillons=" << N_echantillons_ << finl;
}

void Fourier_trans::reprise()
{
  // on a deja lu tous le  fichier de reprise de base on va lire les specifiques pour ici !!!
  Nom nom_sauv = "Sauvegarde_spectrale_";
  // on va indiquer la position dans le mapping, on ne sais jamais des fois que cela bouge !!!!
  int pos = Post_splitting_.get_local_slice_index(2);
  nom_sauv += Nom(pos) + ".sauv";
  Cerr << " Reprise des stats spectrales dans "<< nom_sauv << finl;

  Entree*  fichier = new LecFicDistribue_sansnum(nom_sauv);
  //LecFicDiffuse_JDD fichier(nom_sauv);
  const int N = N_RES;
  TFi_.dimensionner(N);
  TFx_.dimensionner(N);
  TFy_.dimensionner(N);
  N_echantillons_ = 0;

  Param param("Post_Fourier");
  param.ajouter("N_ECHANTILLONS", &N_echantillons_);
  param.ajouter("VECT_K", &tri_ki);
  param.ajouter("REF_KI", &ref_ki);
  for (True_int i = 0; i < N; i++)
    {
      char nom_champ[100];
      sprintf(nom_champ, "TFI_%d",i );
      param.ajouter(nom_champ, &TFi_[i]);
      sprintf(nom_champ, "TFX_%d",i );
      param.ajouter(nom_champ, &TFx_[i]);
      sprintf(nom_champ, "TFY_%d",i );
      param.ajouter(nom_champ, &TFy_[i]);
    }

  param.lire_avec_accolades(*fichier);

  Cout << "Reprise des donnees statistiques: N_echantillons=" << N_echantillons_ << finl;
  Cout << "Reprise des donnees statistiques: N vecteur =" << tri_ki.size_array() << finl;
}
// Impression des donnees pour reprise
Sortie& Fourier_trans::printOn(Sortie& os) const
{
  return os;
}

// Reprise des donnees stat dans un fichier reprise
// Attention, cette methode peut etre appelee avant initialize() !
Entree& Fourier_trans::readOn(Entree& is)
{
  return is;
}
