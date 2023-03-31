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

#include <Statistiques_dns_ijk.h>
#include <IJK_Grid_Geometry.h>
#include <TRUSTTab.h>
#include <communications.h>
#include <IJK_Navier_Stokes_tools.h> // pour initialiser les IJK_double
#include <Statistiques.h>

Implemente_instanciable_sans_constructeur(Statistiques_dns_ijk,"Statistiques_dns_ijk",Objet_U);

#define U_MOY 0
#define V_MOY 1
#define W_MOY 2
#define UU_MOY 3
#define VV_MOY 4
#define WW_MOY 5
#define UV_MOY 6
#define VW_MOY 7
#define UW_MOY 8
#define UP_MOY 9
#define VP_MOY 10
#define WP_MOY 11
#define UUU_MOY 12
#define VVV_MOY 13
#define WWW_MOY 14
#define RHO_MOY 15
#define T_MOY 16
#define UT_MOY 17
#define VT_MOY 18
#define WT_MOY 19
#define URHO_MOY 20
#define VRHO_MOY 21
#define WRHO_MOY 22
#define RHOP_MOY 23
#define TP_MOY 24
#define P_MOY 25
#define T2_MOY 26
#define T3_MOY 27
#define T4_MOY 28
#define RHO2_MOY 29
#define RHO3_MOY 30
#define RHO4_MOY 31
#define UUW_MOY 32
#define VVW_MOY 33
#define U4_MOY 34
#define V4_MOY 35
#define W4_MOY 36
#define MU_MOY 37
#define LAMBDA_MOY 38
#define PP_MOY 39
#define WFACE_MOY 40
#define NU_MOY 41
#define UN_SUR_RHO_MOY 42
#define NU2_MOY 43
#define MU2_MOY 44
#define NUTURB_MOY 45
#define LAMBDADTDZ_MOY 46
#define KAPPATURB_MOY 47
#define RHOWFACE_MOY 48
#define RHOTWFACE_MOY 49
#define RHOUU_MOY 50
#define RHOVV_MOY 51
#define RHOWW_MOY 52
#define RHOUV_MOY 53
#define RHOVW_MOY 54
#define RHOUW_MOY 55
#define RHOUT_MOY 56
#define RHOVT_MOY 57
#define RHOWT_MOY 58
#define NUTURB_XX_MOY 59
#define NUTURB_XY_MOY 60
#define NUTURB_XZ_MOY 61
#define NUTURB_YY_MOY 62
#define NUTURB_YZ_MOY 63
#define NUTURB_ZZ_MOY 64
#define KAPPATURB_X_MOY 65
#define KAPPATURB_Y_MOY 66
#define KAPPATURB_Z_MOY 67
#define NUTURB_XX_DUDX_MOY 68
#define NUTURB_XY_DUDY_MOY 69
#define NUTURB_XZ_DUDZ_MOY 70
#define NUTURB_XY_DVDX_MOY 71
#define NUTURB_YY_DVDY_MOY 72
#define NUTURB_YZ_DVDZ_MOY 73
#define NUTURB_XZ_DWDX_MOY 74
#define NUTURB_YZ_DWDY_MOY 75
#define NUTURB_ZZ_DWDZ_MOY 76
#define KAPPATURB_X_DSCALARDX_MOY 77
#define KAPPATURB_Y_DSCALARDY_MOY 78
#define KAPPATURB_Z_DSCALARDZ_MOY 79
#define STRUCTURAL_UU_MOY 80
#define STRUCTURAL_UV_MOY 81
#define STRUCTURAL_UW_MOY 82
#define STRUCTURAL_VV_MOY 83
#define STRUCTURAL_VW_MOY 84
#define STRUCTURAL_WW_MOY 85
#define STRUCTURAL_USCALAR_MOY 86
#define STRUCTURAL_VSCALAR_MOY 87
#define STRUCTURAL_WSCALAR_MOY 88
#define NUTURB_XX_DUDX_DUDX_MOY 89
#define NUTURB_XY_DUDY_DUDY_MOY 90
#define NUTURB_XZ_DUDZ_DUDZ_MOY 91
#define NUTURB_XY_DVDX_DVDX_MOY 92
#define NUTURB_YY_DVDY_DVDY_MOY 93
#define NUTURB_YZ_DVDZ_DVDZ_MOY 94
#define NUTURB_XZ_DWDX_DWDX_MOY 95
#define NUTURB_YZ_DWDY_DWDY_MOY 96
#define NUTURB_ZZ_DWDZ_DWDZ_MOY 97
#define NUTURB_XY_DVDX_DUDY_MOY 98
#define NUTURB_XZ_DWDX_DUDZ_MOY 99
#define NUTURB_YZ_DWDY_DVDZ_MOY 100
#define KAPPATURB_X_DSCALARDX_DSCALARDX_MOY 101
#define KAPPATURB_Y_DSCALARDY_DSCALARDY_MOY 102
#define KAPPATURB_Z_DSCALARDZ_DSCALARDZ_MOY 103
#define STRUCTURAL_UU_DUDX_MOY 104
#define STRUCTURAL_UV_DUDY_MOY 105
#define STRUCTURAL_UW_DUDZ_MOY 106
#define STRUCTURAL_UV_DVDX_MOY 107
#define STRUCTURAL_VV_DVDY_MOY 108
#define STRUCTURAL_VW_DVDZ_MOY 109
#define STRUCTURAL_UW_DWDX_MOY 110
#define STRUCTURAL_VW_DWDY_MOY 111
#define STRUCTURAL_WW_DWDZ_MOY 112
#define STRUCTURAL_USCALAR_DSCALARDX_MOY 113
#define STRUCTURAL_VSCALAR_DSCALARDY_MOY 114
#define STRUCTURAL_WSCALAR_DSCALARDZ_MOY 115
#define LAMBDADTDZ2_MOY 116

// termes de l'equation d'evolution de la demi-trace du tenseur des correlations de fluctuation de vitesse
#define KW_MOY 0

#define CORRELATION_EC_DIVERGENCE_MOY 1
#define CORRELATION_PRESSION_DIVERGENCE_MOY 2

#define P_WFLUC_MOY 3
#define DIFFUSION_PRESSION_RHOFLUC_ADERIVE_MOY 4
#define DIFFUSION_PRESSION_DERIVRHO_MOY 5

#define DISSIPATION_INCOMPRESSIBLE_UN_AFOISNUMOY_MOY 6
#define DISSIPATION_INCOMPRESSIBLE_DEUX_AFOISNUMOY_MOY 7
#define DISSIPATION_COMPESSIBLE_UN_MOY 8
#define DISSIPATION_COMPESSIBLE_DEUX_MOY 9
#define DISSIPATION_COMPESSIBLE_TROIS_MOY 10

#define DIFFUSION_VISQUEUSE_INCOMPRESSIBLE_DEUX_MOY 11
#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_MOY 12
#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_MOY 13
#define DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS_MOY 14
#define DIFFUSION_VISQUEUSE_DERIVRHO_UN_MOY 15
#define DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_MOY 16
#define DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_MOY 17

#define TERME_NON_CLASSE_MOY 18

#define FLUC_U_MOY 19
#define FLUC_V_MOY 20
#define FLUC_W_MOY 21
#define FLUC_RHO_MOY 22
#define FLUC_NU_MOY 23
#define FLUC_DIV_U_MOY 24
#define TOTAL_DIV_U_MOY 25
#define MEAN_DIV_U_MOY 26
#define TERME_SOURCE_MOY 27

Statistiques_dns_ijk::Statistiques_dns_ijk()
{

  check_converge_ = 0;
  t_integration_k_ = 0.;
  t_integration_ = 0.;

  // Pour generer cette liste, envoyer la liste des #define dans
  //  sed 's/#define /"/g;s/_MOY.*/",/g' >resultat
  // et copier ci-dessous:
  const char *noms_moyennes_prov[] =
  {
    "U",
    "V",
    "W",
    "UU",
    "VV",
    "WW",
    "UV",
    "VW",
    "UW",
    "UP",
    "VP",
    "WP",
    "UUU",
    "VVV",
    "WWW",
    "RHO",
    "T",
    "UT",
    "VT",
    "WT",
    "URHO",
    "VRHO",
    "WRHO",
    "RHOP",
    "TP",
    "P",
    "T2",
    "T3",
    "T4",
    "RHO2",
    "RHO3",
    "RHO4",
    "UUW",
    "VVW",
    "U4",
    "V4",
    "W4",
    "MU",
    "LAMBDA",
    "PP",
    "WFACE",
    "NU",
    "UN_SUR_RHO",
    "NU2",
    "MU2",
    "NUTURB",
    "LAMBDADTDZ",
    "KAPPATURB",
    "RHOWFACE",
    "RHOTWFACE",
    "RHOUU_MOY",
    "RHOVV_MOY",
    "RHOWW_MOY",
    "RHOUV_MOY",
    "RHOVW_MOY",
    "RHOUW_MOY",
    "RHOUT_MOY",
    "RHOVT_MOY",
    "RHOWT_MOY",
    "NUTURB_XX",
    "NUTURB_XY",
    "NUTURB_XZ",
    "NUTURB_YY",
    "NUTURB_YZ",
    "NUTURB_ZZ",
    "KAPPATURB_X",
    "KAPPATURB_Y",
    "KAPPATURB_Z",
    "NUTURB_XX_DUDX",
    "NUTURB_XY_DUDY",
    "NUTURB_XZ_DUDZ",
    "NUTURB_XY_DVDX",
    "NUTURB_YY_DVDY",
    "NUTURB_YZ_DVDZ",
    "NUTURB_XZ_DWDX",
    "NUTURB_YZ_DWDY",
    "NUTURB_ZZ_DWDZ",
    "KAPPATURB_X_DSCALARDX",
    "KAPPATURB_Y_DSCALARDY",
    "KAPPATURB_Z_DSCALARDZ",
    "STRUCTURAL_UU",
    "STRUCTURAL_UV",
    "STRUCTURAL_UW",
    "STRUCTURAL_VV",
    "STRUCTURAL_VW",
    "STRUCTURAL_WW",
    "STRUCTURAL_USCALAR",
    "STRUCTURAL_VSCALAR",
    "STRUCTURAL_WSCALAR",
    "NUTURB_XX_DUDX_DUDX_MOY",
    "NUTURB_XY_DUDY_DUDY_MOY",
    "NUTURB_XZ_DUDZ_DUDZ_MOY",
    "NUTURB_XY_DVDX_DVDX_MOY",
    "NUTURB_YY_DVDY_DVDY_MOY",
    "NUTURB_YZ_DVDZ_DVDZ_MOY",
    "NUTURB_XZ_DWDX_DWDX_MOY",
    "NUTURB_YZ_DWDY_DWDY_MOY",
    "NUTURB_ZZ_DWDZ_DWDZ_MOY",
    "NUTURB_XY_DVDX_DUDY_MOY",
    "NUTURB_XZ_DWDX_DUDZ_MOY",
    "NUTURB_YZ_DWDY_DVDZ_MOY",
    "KAPPATURB_X_DSCALARDX_DSCALARDX_MOY",
    "KAPPATURB_Y_DSCALARDY_DSCALARDY_MOY",
    "KAPPATURB_Z_DSCALARDZ_DSCALARDZ_MOY",
    "STRUCTURAL_UU_DUDX_MOY",
    "STRUCTURAL_UV_DUDY_MOY",
    "STRUCTURAL_UW_DUDZ_MOY",
    "STRUCTURAL_UV_DVDX_MOY",
    "STRUCTURAL_VV_DVDY_MOY",
    "STRUCTURAL_VW_DVDZ_MOY",
    "STRUCTURAL_UW_DWDX_MOY",
    "STRUCTURAL_VW_DWDY_MOY",
    "STRUCTURAL_WW_DWDZ_MOY",
    "STRUCTURAL_USCALAR_DSCALARDX_MOY",
    "STRUCTURAL_VSCALAR_DSCALARDY_MOY",
    "STRUCTURAL_WSCALAR_DSCALARDZ_MOY",
    "LAMBDADTDZ2"
  };
  nval_=117;
  noms_moyennes_.dimensionner(nval_);
  for (int i=0; i<nval_; i++)
    noms_moyennes_[i]=noms_moyennes_prov[i];

  static const char *noms_k_prov[] =
  {
    "KW_MOY",
    "CORRELATION_EC_DIVERGENCE_MOY",
    "CORRELATION_PRESSION_DIVERGENCE_MOY",
    "P_WFLUC_MOY",
    "DIFFUSION_PRESSION_RHOFLUC_ADERIVE_MOY",
    "DIFFUSION_PRESSION_DERIVRHO_MOY",
    "DISSIPATION_INCOMPRESSIBLE_UN_AFOISNUMOY_MOY",
    "DISSIPATION_INCOMPRESSIBLE_DEUX_AFOISNUMOY_MOY",
    "DISSIPATION_COMPESSIBLE_UN_MOY",
    "DISSIPATION_COMPESSIBLE_DEUX_MOY",
    "DISSIPATION_COMPESSIBLE_TROIS_MOY",
    "DIFFUSION_VISQUEUSE_INCOMPRESSIBLE_DEUX_MOY",
    "DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_MOY",
    "DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_MOY",
    "DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS_MOY",
    "DIFFUSION_VISQUEUSE_DERIVRHO_UN_MOY",
    "DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_MOY",
    "DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_MOY",
    "TERME_NON_CLASSE_MOY",
    "FLUC_U_MOY",
    "FLUC_V_MOY",
    "FLUC_W_MOY",
    "FLUC_RHO_MOY",
    "FLUC_NU_MOY",
    "FLUC_DIV_U_MOY",
    "TOTAL_DIV_U_MOY",
    "MEAN_DIV_U_MOY",
    "TERME_SOURCE_MOY",
  };
  kval_=28;
  noms_k_.dimensionner(kval_);
  for (int i=0; i<kval_; i++)
    noms_k_[i]=noms_k_prov[i];
}

/*
static inline double calculer_lambda_air(double temperature)
{
  const double fac_a = -5.05628e-18;
  const double fac_b = 2.469e-14;
  const double fac_c = -4.98344e-11;
  const double fac_d = 7.06714e-08;
  const double fac_e = 1.0894e-06;
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
  const double fac_a = -5.05628e-18;
  const double fac_b = 2.469e-14;
  const double fac_c = -4.98344e-11;
  const double fac_d = 7.06714e-08;
  const double fac_e = 1.0894e-06;

  double val = temperature;
  double calc = val * fac_a + fac_b;
  calc = val * calc  + fac_c;
  calc = val * calc  + fac_d;
  calc = val * calc  + fac_e;
  return calc;
} */

static inline double calculer_rho_air(double temperature, double constante_specifique_gaz, double pression_thermodynamique)
{
  return ( pression_thermodynamique / ( constante_specifique_gaz * temperature ) );
}

// calcule la derivee premiere en maillage anisotrope !!
static inline double derivee_aniso( double alpha /* rapport des distances a la maille sup et inf */
                                    , double un_sur_alpha /* inverse de alpha */
                                    , double un_sur_dp_plus_dm /* inverse de la somme des distances aux elems inf et sup */
                                    , double val_moins /* valeur dans la maille inf */
                                    , double val_centre /* valeur dans la maille */
                                    , double val_plus /* valeur dans la maille sup */ )
{
  double derivee  = (val_plus - val_centre) * un_sur_alpha; // contribution sup pondere
  derivee += (val_centre - val_moins ) * alpha; // contribution moins pondere
  derivee *= un_sur_dp_plus_dm; // division
  return derivee;
}
//static inline double derivee_2_aniso( double un_sur_delta_moins /* inverse da la distance a la maille inf */
//                                     , double un_sur_delta_plus /* inverse da la distance a la maille sup */
//                                      , double un_sur_dp_plus_dm /* inverse de la somme des distances aux elems inf et sup */
//                                      , double val_moins /* valeur dans la maille inf */
//                                      , double val_centre/* valeur dans la maille */
//                                      , double val_plus  /* valeur dans la maille sup */)
//{
//  double derivee_2  = (val_plus - val_centre) * un_sur_delta_plus; // contribution sup pondere
//  derivee_2 -= (val_centre - val_moins ) * un_sur_delta_moins; // contribution moins pondere
//  derivee_2 *= 2. * un_sur_dp_plus_dm; // division (le 2. disparait naturellement lorsque dp = dm => 2* 1/ (dp + dm ) = 2*1/ 2d  = 1/d
//  return derivee_2;
//}

// Attention, cette methode est appelee apres readOn(),
// il ne faut pas casser les donnees lues
void Statistiques_dns_ijk::initialize(const IJK_Grid_Geometry& geom,double T_KMAX , double T_KMIN, double constante_specifique_gaz)
{
  // F.A ajout recuperation des temperatures des CLs
  TCL_kmax_ = T_KMAX;
  TCL_kmin_ = T_KMIN;
  constante_specifique_gaz_ = constante_specifique_gaz;
  // check_converge_ =1; // !!! uniquement pour tester les stats k !!!
  dx_=geom.get_constant_delta(0); //modif AT 20/06/2013
  dy_=geom.get_constant_delta(1); //modif AT 20/06/2013
  tab_dz_=geom.get_delta(2); //modif AT 20/06/2013
  const int n = geom.get_nb_elem_tot(DIRECTION_K);
  elem_coord_.resize_array(n);
  int i;
  const ArrOfDouble& coord_z = geom.get_node_coordinates(DIRECTION_K);
  for (i = 0; i < n; i++)
    elem_coord_[i] = (coord_z[i] + coord_z[i+1]) * 0.5;

  if (Process::je_suis_maitre())
    {
      moyenne_spatiale_instantanee_.dimensionner(nval_);
      for (i = 0; i < nval_; i++)
        {
          moyenne_spatiale_instantanee_[i].resize_array(n);
        }
      moyenne_spatiale_instantanee_[WFACE_MOY].resize_array(n+1);//En effet, il y a une face de plus que d'elem et ce tableau est au face
    }

  int flag = 0;
  if (integrale_temporelle_.size() == 0)
    {
      integrale_temporelle_.dimensionner(nval_);
      t_integration_ = 0.;
    }
  else
    flag = 1;

  for (i = 0; i < nval_; i++)
    {
      if (flag)
        {
          if ( ((integrale_temporelle_[i].size_array() != n) && (i!=WFACE_MOY)) ||
               ((integrale_temporelle_[WFACE_MOY].size_array() != n+1) && (i==WFACE_MOY)) )
            {
              Cerr << "Erreur dans Statistiques_dns_ijk::initialize: reprise avec mauvais nombre de mailles en z" << finl;
              Process::exit();
            }
        }
      else
        {
          if (i!=WFACE_MOY)
            integrale_temporelle_[i].resize_array(n);
          else
            integrale_temporelle_[i].resize_array(n+1);
        }
    }
  //modif AT 20/06/2013
  if (check_converge_)
    {
      // DD 16/10/2015: on recupere la masse_volumique moyenne pour calculer les fluctuations de masse_volumique
      rho_moy_.resize_array(n);
      rho_moy_=integrale_temporelle_[15];
      rho_moy_/=t_integration_; // Il ne faut pas oublier de diviser par le temps d'integration !!!
      // DD 16/10/2015: on recupere la viscosite_cinematique moyenne pour calculer les fluctuations de viscosite_cinematique
      nu_moy_.resize_array(n);
      nu_moy_=integrale_temporelle_[41];
      nu_moy_/=t_integration_; // Il ne faut pas oublier de diviser par le temps d'integration !!!

      //on recupere les vitesses moyennes aux faces pour calculer les fluctuations de vitesse.
      vit_moy_.dimensionner(3);
      for (i = 0; i < 3; i++)
        vit_moy_[i].resize_array(n);

      vit_moy_[0]=integrale_temporelle_[0];
      vit_moy_[0]/=t_integration_; // Il ne faut pas oublier de diviser par le temps d'integration !!!
      vit_moy_[1]=integrale_temporelle_[1];
      vit_moy_[1]/=t_integration_;
      vit_moy_[2]=integrale_temporelle_[40];//car on veut w aux faces !
      vit_moy_[2]/=t_integration_;
      int flag2 = 0;
      if ((integrale_k_.size() == 0 )|| (integrale_k_[0].size_array() ==0))
        {
          integrale_k_.dimensionner(kval_);
          t_integration_k_ = 0.;
        }
      else
        flag2 = 1;

      for (i = 0; i < kval_; i++)
        {
          if (flag2)
            {
              if (integrale_k_[i].size_array() != n)
                {
                  Cerr << "Erreur dans Statistiques_dns_ijk::initialize: reprise stat k avec mauvais nombre de mailles en z" << finl;
                  Process::exit();
                }
            }
          else
            {
              integrale_k_[i].resize_array(n);
            }
        }
      if (Process::je_suis_maitre())
        {
          moyenne_spatiale_ec_.dimensionner(kval_);
          for (int ii = 0; ii < kval_; ii++)
            {
              moyenne_spatiale_ec_[ii].resize_array(n);
            }
        }
    }
}

void Statistiques_dns_ijk::update_stat(const FixedVector<IJK_Field_double, 3>& vitesse,
                                       const IJK_Field_double& pression,
                                       const IJK_Field_double& temperature,
                                       const IJK_Field_double& masse_vol,
                                       const IJK_Field_double& champ_mu,
                                       const IJK_Field_double& champ_lambda,
                                       const ArrOfDouble_with_ghost& delta_z_local_pour_delta,
                                       const bool flag_nu_anisotropic,
                                       const int flag_turbulent_viscosity,
                                       const IJK_Field_double& champ_turbulent_mu_xx,
                                       const IJK_Field_double& champ_turbulent_mu_xy,
                                       const IJK_Field_double& champ_turbulent_mu_xz,
                                       const IJK_Field_double& champ_turbulent_mu_yy,
                                       const IJK_Field_double& champ_turbulent_mu_yz,
                                       const IJK_Field_double& champ_turbulent_mu_zz,
                                       const bool flag_kappa_anisotropic,
                                       const int flag_turbulent_diffusivity,
                                       const IJK_Field_double& champ_turbulent_kappa_x,
                                       const IJK_Field_double& champ_turbulent_kappa_y,
                                       const IJK_Field_double& champ_turbulent_kappa_z,
                                       const int flag_structural_uu,
                                       const FixedVector<IJK_Field_double, 6>& structural_uu_tensor,
                                       const int flag_structural_uscalar,
                                       const FixedVector<IJK_Field_double, 3>& structural_uscalar_vector,
                                       const int flag_formulation_favre,
                                       const int flag_formulation_velocity,
                                       const double cp_gaz,
                                       const double pression_thermodynamique,
                                       double dt)
{
  if (elem_coord_.size_array() == 0)
    {
      Cerr << "Erreur dans Statistiques_dns_ijk::update_stat: non initialise" << finl;
      Process::exit();
    }

  const IJK_Field_double& vitesse_i = vitesse[0];
  const IJK_Field_double& vitesse_j = vitesse[1];
  const IJK_Field_double& vitesse_k = vitesse[2];

  const IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  const IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  const IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  const IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  const IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  const IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const IJK_Field_double& structural_uscalar_x = structural_uscalar_vector[0];
  const IJK_Field_double& structural_uscalar_y = structural_uscalar_vector[1];
  const IJK_Field_double& structural_uscalar_z = structural_uscalar_vector[2];

  // Nombre total de mailles en K
  const int nktot = pression.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);

  // Nombre local de mailles en K
  const int kmax = pression.nk();
  const int offset = pression.get_splitting().get_offset_local(DIRECTION_K);

  DoubleTab tmp(nktot, nval_);
  const int imax = pression.ni();
  const int jmax = pression.nj();

  const double rho_kmin = calculer_rho_air(TCL_kmin_, constante_specifique_gaz_, pression_thermodynamique);
  const double rho_kmax = calculer_rho_air(TCL_kmax_, constante_specifique_gaz_, pression_thermodynamique);

  const double unsurdx = 1./dx_;
  const double unsurdy = 1./dy_;
  const double nu_dx_delta = flag_nu_anisotropic ? dx_ : 1.;
  const double nu_dy_delta = flag_nu_anisotropic ? dy_ : 1.;
  const double kappa_dx_delta = flag_kappa_anisotropic ? dx_ : 1.;
  const double kappa_dy_delta = flag_kappa_anisotropic ? dy_ : 1.;

  for (int k = 0; k < kmax; k++)
    {
      const int kg=k+offset;
      double unsurdz = 1./tab_dz_[kg];
      double delta_m = kg== 0 ? tab_dz_[kg] * 0.5 : ( tab_dz_[kg] + tab_dz_[kg-1] ) * 0.5;
      double delta_p = kg==(nktot-1) ? tab_dz_[kg] * 0.5 : ( tab_dz_[kg] + tab_dz_[kg+1] ) * 0.5;
      double alpha = delta_p/delta_m;
      double unsalpha = delta_m/delta_p;
      double unsdpdm = 1./(delta_p + delta_m);
      double unsurdm = 1./delta_m;
      double unsurdp = 1./delta_p;
      double nu_dz_delta = flag_nu_anisotropic ? delta_z_local_pour_delta[k] : 1.;
      double nu_delta_m_delta = flag_nu_anisotropic ? (kg==0 ? 0. : (delta_z_local_pour_delta[k] + delta_z_local_pour_delta[k-1])*0.5) : 1.;
      double nu_delta_p_delta = flag_nu_anisotropic ? (kg==(nktot-1) ? 0. : (delta_z_local_pour_delta[k] + delta_z_local_pour_delta[k+1])*0.5) : 1.;
      double kappa_delta_m_delta = flag_kappa_anisotropic ? (kg==0 ? 0. : (delta_z_local_pour_delta[k] + delta_z_local_pour_delta[k-1])*0.5) : 1.;
      double kappa_delta_p_delta = flag_kappa_anisotropic ? (kg==(nktot-1) ? 0. : (delta_z_local_pour_delta[k] + delta_z_local_pour_delta[k+1])*0.5) : 1.;

      // Calcul des moyennes spatiales sur le plan ij
      // On y stocke la somme de toutes les valeurs sur le plan ij, on divisera apres
      ArrOfDouble moy(nval_);
      for (int i = 0; i < nval_; i++)
        {
          moy[i] = 0.;
        }
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              // Vitesses au centre de la maille i,j,k
              // On peut ici interpoler plus finement si on veut:
              double uf = vitesse_i(i,j,k);
              double vf = vitesse_j(i,j,k);
              double wf_k = vitesse_k(i,j,k);
              double wf_kp1 = vitesse_k(i,j,k+1);
              double ue = (vitesse_i(i,j,k) + vitesse_i(i+1, j, k)) * 0.5;
              double ve = (vitesse_j(i,j,k) + vitesse_j(i, j+1, k)) * 0.5;
              double we = (vitesse_k(i,j,k) + vitesse_k(i, j, k+1)) * 0.5;
              double p = pression(i,j,k);

              double t = temperature(i,j,k);
              double t_plus_k  = kg==(nktot-1) ? TCL_kmax_ : temperature(i,j,k+1);
              double t_moins_k = kg==0 ? TCL_kmin_ : temperature(i,j,k-1);
              double dtdz = derivee_aniso(alpha, unsalpha, unsdpdm, t_moins_k, t, t_plus_k);
              double tf_k = (t + t_moins_k) * 0.5;
              //double tf_kp1 = (t + t_plus_k) * 0.5;

              double rho = masse_vol(i,j,k);
              double rho_moins_k = kg==0 ? rho_kmin : masse_vol(i,j,k-1);
              double rhof_k = (rho + rho_moins_k) * 0.5;
              double rho_plus_k  = kg==(nktot-1) ? rho_kmax : masse_vol(i,j,k+1);
              //double rhof_kp1 = (rho + rho_plus_k) * 0.5;
              //double drhodz = derivee_aniso(alpha, unsalpha, unsdpdm, rho_moins_k, rho, rho_plus_k);

              double mu = champ_mu(i,j,k);
              double lambda = champ_lambda(i,j,k);
              double nu = mu/rho;
              double unsrho = 1/rho;

              // derivees de vitesse pour viscosite turbulente
              const double nu_deltaunsurdx = nu_dx_delta * unsurdx;
              const double nu_deltaunsurdy = nu_dy_delta * unsurdy;

              const double nu_deltaunsurdz = nu_dz_delta * unsurdz;
              const double nu_deltaunsurdelta_m = nu_delta_m_delta * unsurdm;
              const double nu_deltaunsurdelta_p = nu_delta_p_delta * unsurdp;

              const double uf_i = vitesse_i(i,j,k);
              const double uf_ip1 = vitesse_i(i+1,j,k);
              const double uf_i_jm1 = vitesse_i(i,j-1,k);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);
              const double uf_i_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i,j,k+1);

              const double vf_j = vitesse_j(i,j,k);
              const double vf_jp1 = vitesse_j(i,j+1,k);
              const double vf_j_im1 = vitesse_j(i-1,j,k);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);
              const double vf_j_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j,k+1);

              const double wf_k_im1 = vitesse_k(i-1,j,k);
              const double wf_k_jm1 = vitesse_k(i,j-1,k);
              const double wf_kp1_im1 = kg==(nktot-1) ? 0. : vitesse_k(i-1,j,k+1);
              const double wf_kp1_jm1 = kg==(nktot-1) ? 0. : vitesse_k(i,j-1,k+1);

              const double duidx = nu_deltaunsurdx * (uf_ip1 - uf_i);
              const double duidy_ij = nu_deltaunsurdy * (uf_i - uf_i_jm1);
              const double duidz_ik = nu_deltaunsurdelta_m * (uf_i - uf_i_km1);
              const double duidz_ikp1 = nu_deltaunsurdelta_p * (uf_i_kp1 - uf_i);
              const double dujdx_ij = nu_deltaunsurdx * (vf_j - vf_j_im1);
              const double dujdy = nu_deltaunsurdy * (vf_jp1 - vf_j);
              const double dujdz_jk = nu_deltaunsurdelta_m * (vf_j - vf_j_km1);
              const double dujdz_jkp1 = nu_deltaunsurdelta_p * (vf_j_kp1 - vf_j);
              const double dukdx_ik = nu_deltaunsurdx * (wf_k - wf_k_im1);
              const double dukdx_ikp1 = nu_deltaunsurdx * (wf_kp1 - wf_kp1_im1);
              const double dukdy_jk = nu_deltaunsurdy * (wf_k - wf_k_jm1);
              const double dukdy_jkp1 = nu_deltaunsurdy * (wf_kp1 - wf_kp1_jm1);
              const double dukdz = nu_deltaunsurdz * (wf_kp1 - wf_k);

              // scalaire pour diffusivite turbulente
              const double kappa_deltaunsurdx = kappa_dx_delta * unsurdx;
              const double kappa_deltaunsurdy = kappa_dy_delta * unsurdy;

              const double kappa_deltaunsurdelta_m = kappa_delta_m_delta * unsurdm;
              const double kappa_deltaunsurdelta_p = kappa_delta_p_delta * unsurdp;

              double scalar = 0.;
              double scalar_im1 = 0.;
              double scalar_jm1 = 0.;
              double scalar_km1 = 0.;
              double scalar_kp1 = 0.;

              if (flag_formulation_favre)
                {
                  scalar = 1/masse_vol(i,j,k);
                  scalar_im1 = 1/masse_vol(i-1,j,k);
                  scalar_jm1 = 1/masse_vol(i,j-1,k);
                  scalar_km1 = kg==0 ? 1/rho_kmin : 1/masse_vol(i,j,k-1);
                  scalar_kp1 = kg==(nktot-1) ? 1/rho_kmax : 1/masse_vol(i,j,k+1);
                }
              if (flag_formulation_velocity)
                {
                  scalar = masse_vol(i,j,k);
                  scalar_im1 = masse_vol(i-1,j,k);
                  scalar_jm1 = masse_vol(i,j-1,k);
                  scalar_km1 = kg==0 ? rho_kmin : masse_vol(i,j,k-1);
                  scalar_kp1 = kg==(nktot-1) ? rho_kmax : masse_vol(i,j,k+1);
                }
              const double dscalardxf_i = kappa_deltaunsurdx * (scalar - scalar_im1);
              const double dscalardyf_j = kappa_deltaunsurdy * (scalar - scalar_jm1);
              const double dscalardzf_k = kappa_deltaunsurdelta_m * (scalar - scalar_km1);
              const double dscalardzf_kp1 = kappa_deltaunsurdelta_p * (scalar_kp1 - scalar);

              // viscosite ou diffusivite turbulente
              double nuturb_xx = 0.;
              double nuturb_xy = 0.;
              double nuturb_xz = 0.;
              double nuturb_yy = 0.;
              double nuturb_yz = 0.;
              double nuturb_zz = 0.;

              double nuturb_xx_dudx = 0.;
              double nuturb_xy_dudy = 0.;
              double nuturb_xz_dudz = 0.;
              double nuturb_xy_dvdx = 0.;
              double nuturb_yy_dvdy = 0.;
              double nuturb_yz_dvdz = 0.;
              double nuturb_xz_dwdx = 0.;
              double nuturb_yz_dwdy = 0.;
              double nuturb_zz_dwdz = 0.;

              double nuturb_xx_dudx_dudx = 0.;
              double nuturb_xy_dudy_dudy = 0.;
              double nuturb_xz_dudz_dudz = 0.;
              double nuturb_xy_dvdx_dvdx = 0.;
              double nuturb_yy_dvdy_dvdy = 0.;
              double nuturb_yz_dvdz_dvdz = 0.;
              double nuturb_xz_dwdx_dwdx = 0.;
              double nuturb_yz_dwdy_dwdy = 0.;
              double nuturb_zz_dwdz_dwdz = 0.;
              double nuturb_xy_dvdx_dudy = 0.;
              double nuturb_xz_dwdx_dudz = 0.;
              double nuturb_yz_dwdy_dvdz = 0.;

              if (flag_turbulent_viscosity)
                {
                  double nuturb_xy_im1 = 0.;
                  double nuturb_xy_jm1 = 0.;
                  double nuturb_xy_im1_jm1 = 0.;

                  double nuturb_xz_im1 = 0.;
                  double nuturb_xz_km1 = 0.;
                  double nuturb_xz_im1_km1 = 0.;
                  double nuturb_xz_kp1 = 0.;
                  double nuturb_xz_im1_kp1 = 0.;

                  double nuturb_yz_jm1 = 0.;
                  double nuturb_yz_km1 = 0.;
                  double nuturb_yz_jm1_km1 = 0.;
                  double nuturb_yz_kp1 = 0.;
                  double nuturb_yz_jm1_kp1 = 0.;

                  if (flag_formulation_favre)
                    {
                      const double rho_im1 = masse_vol(i-1,j,k);
                      const double rho_jm1 = masse_vol(i,j-1,k);
                      const double rho_im1_jm1 = masse_vol(i-1,j-1,k);
                      const double rho_im1_km1 = kg==0 ? rho_kmin : masse_vol(i-1,j,k-1);
                      const double rho_im1_kp1 = kg==(nktot-1) ? rho_kmax : masse_vol(i-1,j,k+1);
                      const double rho_jm1_km1 = kg==0 ? rho_kmin : masse_vol(i,j-1,k-1);
                      const double rho_jm1_kp1 = kg==(nktot-1) ? rho_kmax : masse_vol(i,j-1,k+1);

                      nuturb_xx = champ_turbulent_mu_xx(i,j,k)/rho;

                      nuturb_xy = champ_turbulent_mu_xy(i,j,k)/rho;
                      nuturb_xy_im1 = champ_turbulent_mu_xy(i-1,j,k)/rho_im1;
                      nuturb_xy_jm1 = champ_turbulent_mu_xy(i,j-1,k)/rho_jm1;
                      nuturb_xy_im1_jm1 = champ_turbulent_mu_xy(i-1,j-1,k)/rho_im1_jm1;

                      nuturb_xz = champ_turbulent_mu_xz(i,j,k)/rho;
                      nuturb_xz_im1 = champ_turbulent_mu_xz(i-1,j,k)/rho_im1;
                      nuturb_xz_km1 = kg==0 ? 0. : champ_turbulent_mu_xz(i,j,k-1)/rho_moins_k;
                      nuturb_xz_im1_km1 = kg==0 ? 0. : champ_turbulent_mu_xz(i-1,j,k-1)/rho_im1_km1;
                      nuturb_xz_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_xz(i,j,k+1)/rho_plus_k;
                      nuturb_xz_im1_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_xz(i-1,j,k+1)/rho_im1_kp1;

                      nuturb_yy = champ_turbulent_mu_yy(i,j,k)/rho;

                      nuturb_yz = champ_turbulent_mu_yz(i,j,k)/rho;
                      nuturb_yz_jm1 = champ_turbulent_mu_yz(i,j-1,k)/rho_jm1;
                      nuturb_yz_km1 = kg==0 ? 0. : champ_turbulent_mu_yz(i,j,k-1)/rho_moins_k;
                      nuturb_yz_jm1_km1 = kg==0 ? 0. : champ_turbulent_mu_yz(i,j-1,k-1)/rho_jm1_km1;
                      nuturb_yz_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_yz(i,j,k+1)/rho_plus_k;
                      nuturb_yz_jm1_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_yz(i,j-1,k+1)/rho_jm1_kp1;

                      nuturb_zz = champ_turbulent_mu_zz(i,j,k)/rho;
                    }
                  if (flag_formulation_velocity)
                    {
                      nuturb_xx = champ_turbulent_mu_xx(i,j,k);

                      nuturb_xy = champ_turbulent_mu_xy(i,j,k);
                      nuturb_xy_im1 = champ_turbulent_mu_xy(i-1,j,k);
                      nuturb_xy_jm1 = champ_turbulent_mu_xy(i,j-1,k);
                      nuturb_xy_im1_jm1 = champ_turbulent_mu_xy(i-1,j-1,k);

                      nuturb_xz = champ_turbulent_mu_xz(i,j,k);
                      nuturb_xz_im1 = champ_turbulent_mu_xz(i-1,j,k);
                      nuturb_xz_km1 = kg==0 ? 0. : champ_turbulent_mu_xz(i,j,k-1);
                      nuturb_xz_im1_km1 = kg==0 ? 0. : champ_turbulent_mu_xz(i-1,j,k-1);
                      nuturb_xz_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_xz(i,j,k+1);
                      nuturb_xz_im1_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_xz(i-1,j,k+1);

                      nuturb_yy = champ_turbulent_mu_yy(i,j,k);

                      nuturb_yz = champ_turbulent_mu_yz(i,j,k);
                      nuturb_yz_jm1 = champ_turbulent_mu_yz(i,j-1,k);
                      nuturb_yz_km1 = kg==0 ? 0. : champ_turbulent_mu_yz(i,j,k-1);
                      nuturb_yz_jm1_km1 = kg==0 ? 0. : champ_turbulent_mu_yz(i,j-1,k-1);
                      nuturb_yz_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_yz(i,j,k+1);
                      nuturb_yz_jm1_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_mu_yz(i,j-1,k+1);

                      nuturb_zz = champ_turbulent_mu_zz(i,j,k);
                    }

                  const double nuturb_xy_ij = 0.25 * (nuturb_xy + nuturb_xy_im1 + nuturb_xy_jm1 + nuturb_xy_im1_jm1);
                  const double nuturb_xz_ik = 0.25 * (nuturb_xz + nuturb_xz_im1 + nuturb_xz_km1 + nuturb_xz_im1_km1);
                  const double nuturb_xz_ikp1 = 0.25 * (nuturb_xz_kp1 + nuturb_xz_im1_kp1 + nuturb_xz + nuturb_xz_im1);
                  const double nuturb_yz_jk = 0.25 * (nuturb_yz + nuturb_yz_jm1 + nuturb_yz_km1 + nuturb_yz_jm1_km1);
                  const double nuturb_yz_jkp1 = 0.25 * (nuturb_yz_kp1 + nuturb_yz_jm1_kp1 + nuturb_yz + nuturb_yz_jm1);

                  nuturb_xx_dudx = nuturb_xx*duidx;
                  nuturb_xy_dudy = nuturb_xy_ij*duidy_ij;
                  nuturb_xz_dudz = 0.5 * (nuturb_xz_ikp1*duidz_ikp1 + nuturb_xz_ik*duidz_ik);
                  nuturb_xy_dvdx = nuturb_xy_ij*dujdx_ij;
                  nuturb_yy_dvdy = nuturb_yy*dujdy;
                  nuturb_yz_dvdz = 0.5 * (nuturb_yz_jkp1*dujdz_jkp1 + nuturb_yz_jk*dujdz_jk);
                  nuturb_xz_dwdx = 0.5 * (nuturb_xz_ikp1*dukdx_ikp1 + nuturb_xz_ik*dukdx_ik);
                  nuturb_yz_dwdy = 0.5 * (nuturb_yz_jkp1*dukdy_jkp1 + nuturb_yz_jk*dukdy_jk);
                  nuturb_zz_dwdz = nuturb_zz*dukdz;

                  nuturb_xx_dudx_dudx = nuturb_xx*duidx*duidx;
                  nuturb_xy_dudy_dudy = nuturb_xy_ij*duidy_ij*duidy_ij;
                  nuturb_xz_dudz_dudz = 0.5 * (nuturb_xz_ikp1*duidz_ikp1*duidz_ikp1 + nuturb_xz_ik*duidz_ik*duidz_ik);
                  nuturb_xy_dvdx_dvdx = nuturb_xy_ij*dujdx_ij*dujdx_ij;
                  nuturb_yy_dvdy_dvdy = nuturb_yy*dujdy*dujdy;
                  nuturb_yz_dvdz_dvdz = 0.5 * (nuturb_yz_jkp1*dujdz_jkp1*dujdz_jkp1 + nuturb_yz_jk*dujdz_jk*dujdz_jk);
                  nuturb_xz_dwdx_dwdx = 0.5 * (nuturb_xz_ikp1*dukdx_ikp1*dukdx_ikp1 + nuturb_xz_ik*dukdx_ik*dukdx_ik);
                  nuturb_yz_dwdy_dwdy = 0.5 * (nuturb_yz_jkp1*dukdy_jkp1*dukdy_jkp1 + nuturb_yz_jk*dukdy_jk*dukdy_jk);
                  nuturb_zz_dwdz_dwdz = nuturb_zz*dukdz*dukdz;
                  nuturb_xy_dvdx_dudy = nuturb_xy_ij*dujdx_ij*duidy_ij;
                  nuturb_xz_dwdx_dudz = 0.5 * (nuturb_xz_ikp1*dukdx_ikp1*duidz_ikp1 + nuturb_xz_ik*dukdx_ik*duidz_ik);
                  nuturb_yz_dwdy_dvdz = 0.5 * (nuturb_yz_jkp1*dukdy_jkp1*dujdz_jkp1 + nuturb_yz_jk*dukdy_jk*dujdz_jk);
                }

              double kappaturb_x = 0.;
              double kappaturb_y = 0.;
              double kappaturb_z = 0.;

              double kappaturb_x_dscalardx = 0.;
              double kappaturb_y_dscalardy = 0.;
              double kappaturb_z_dscalardz = 0.;

              double kappaturb_x_dscalardx_dscalardx = 0.;
              double kappaturb_y_dscalardy_dscalardy = 0.;
              double kappaturb_z_dscalardz_dscalardz = 0.;

              if (flag_turbulent_diffusivity)
                {
                  double kappaturb_x_im1 = 0.;
                  double kappaturb_y_jm1 = 0.;
                  double kappaturb_z_km1 = 0.;
                  double kappaturb_z_kp1 = 0.;

                  if (flag_formulation_favre)
                    {
                      const double rho_im1 = masse_vol(i-1,j,k);
                      const double rho_jm1 = masse_vol(i,j-1,k);

                      kappaturb_x = champ_turbulent_kappa_x(i,j,k)/(rho*cp_gaz);
                      kappaturb_x_im1 = champ_turbulent_kappa_x(i-1,j,k)/(rho_im1*cp_gaz);

                      kappaturb_y = champ_turbulent_kappa_y(i,j,k)/(rho*cp_gaz);
                      kappaturb_y_jm1 = champ_turbulent_kappa_y(i,j-1,k)/(rho_jm1*cp_gaz);

                      kappaturb_z = champ_turbulent_kappa_z(i,j,k)/(rho*cp_gaz);
                      kappaturb_z_km1 = kg==0 ? 0. : champ_turbulent_kappa_z(i,j,k-1)/(rho_moins_k*cp_gaz);
                      kappaturb_z_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_kappa_z(i,j,k+1)/(rho_plus_k*cp_gaz);
                    }
                  if (flag_formulation_velocity)
                    {
                      kappaturb_x = champ_turbulent_kappa_x(i,j,k);
                      kappaturb_x_im1 = champ_turbulent_kappa_x(i-1,j,k);

                      kappaturb_y = champ_turbulent_kappa_y(i,j,k);
                      kappaturb_y_jm1 = champ_turbulent_kappa_y(i,j-1,k);

                      kappaturb_z = champ_turbulent_kappa_z(i,j,k);
                      kappaturb_z_km1 = kg==0 ? 0. : champ_turbulent_kappa_z(i,j,k-1);
                      kappaturb_z_kp1 = kg==(nktot-1) ? 0. : champ_turbulent_kappa_z(i,j,k+1);
                    }

                  const double kappaturb_xf_i = 0.5 * (kappaturb_x + kappaturb_x_im1);
                  const double kappaturb_yf_j = 0.5 * (kappaturb_y + kappaturb_y_jm1);
                  const double kappaturb_zf_k = 0.5 * (kappaturb_z + kappaturb_z_km1);
                  const double kappaturb_zf_kp1 = 0.5 * (kappaturb_z + kappaturb_z_kp1);

                  kappaturb_x_dscalardx = kappaturb_xf_i*dscalardxf_i;
                  kappaturb_y_dscalardy = kappaturb_yf_j*dscalardyf_j;
                  kappaturb_z_dscalardz = 0.5 * (kappaturb_zf_k*dscalardzf_k + kappaturb_zf_kp1*dscalardzf_kp1);

                  kappaturb_x_dscalardx_dscalardx = kappaturb_xf_i*dscalardxf_i*dscalardxf_i;
                  kappaturb_y_dscalardy_dscalardy = kappaturb_yf_j*dscalardyf_j*dscalardyf_j;
                  kappaturb_z_dscalardz_dscalardz = 0.5 * (kappaturb_zf_k*dscalardzf_k*dscalardzf_k + kappaturb_zf_kp1*dscalardzf_kp1*dscalardzf_kp1);
                }

              // terme sous maille structurel
              double structural_uu = 0.;
              double structural_uv = 0.;
              double structural_uw = 0.;
              double structural_vv = 0.;
              double structural_vw = 0.;
              double structural_ww = 0.;

              double structural_uu_dudx = 0.;
              double structural_uv_dudy = 0.;
              double structural_uw_dudz = 0.;
              double structural_uv_dvdx = 0.;
              double structural_vv_dvdy = 0.;
              double structural_vw_dvdz = 0.;
              double structural_uw_dwdx = 0.;
              double structural_vw_dwdy = 0.;
              double structural_ww_dwdz = 0.;

              if (flag_structural_uu)
                {
                  const double structural_uu_e = structural_uu_xx(i,j,k);

                  const double structural_uv_ij = structural_uu_xy(i,j,k);

                  const double structural_uw_ik = structural_uu_xz(i,j,k);
                  const double structural_uw_ikp1 = kg==(nktot-1) ? 0. : structural_uu_xz(i,j,k+1);

                  const double structural_vv_e = structural_uu_yy(i,j,k);

                  const double structural_vw_jk = structural_uu_yz(i,j,k);
                  const double structural_vw_jkp1 = kg==(nktot-1) ? 0. : structural_uu_yz(i,j,k+1);

                  const double structural_ww_e = structural_uu_zz(i,j,k);

                  structural_uu = structural_uu_e;
                  structural_uv = structural_uv_ij;
                  structural_uw = 0.5 * (structural_uw_ikp1 + structural_uw_ik);
                  structural_vv = structural_vv_e;
                  structural_vw = 0.5 * (structural_vw_jkp1 + structural_vw_jk);
                  structural_ww = structural_ww_e;

                  structural_uu_dudx = structural_uu_e*duidx;
                  structural_uv_dudy = structural_uv_ij*duidy_ij;
                  structural_uw_dudz = 0.5 * (structural_uw_ikp1*duidz_ikp1 + structural_uw_ik*duidz_ik);
                  structural_uv_dvdx = structural_uv_ij*dujdx_ij;
                  structural_vv_dvdy = structural_vv_e*dujdy;
                  structural_vw_dvdz = 0.5 * (structural_vw_jkp1*dujdz_jkp1 + structural_vw_jk*dujdz_jk);
                  structural_uw_dwdx = 0.5 * (structural_uw_ikp1*dukdx_ikp1 + structural_uw_ik*dukdx_ik);
                  structural_vw_dwdy = 0.5 * (structural_vw_jkp1*dukdy_jkp1 + structural_vw_jk*dukdy_jk);
                  structural_ww_dwdz = structural_ww_e*dukdz;
                }

              double structural_uscalar = 0.;
              double structural_vscalar = 0.;
              double structural_wscalar = 0.;

              double structural_uscalar_dscalardx = 0.;
              double structural_vscalar_dscalardy = 0.;
              double structural_wscalar_dscalardz = 0.;

              if (flag_structural_uscalar)
                {
                  const double structural_uscalarf_i = structural_uscalar_x(i,j,k);

                  const double structural_vscalarf_j = structural_uscalar_y(i,j,k);

                  const double structural_wscalarf_k = structural_uscalar_z(i,j,k);
                  const double structural_wscalarf_kp1 = kg==(nktot-1) ? 0. : structural_uscalar_z(i,j,k+1);

                  structural_uscalar = structural_uscalarf_i;
                  structural_vscalar = structural_vscalarf_j;
                  structural_wscalar = 0.5 * (structural_wscalarf_kp1 + structural_wscalarf_k);

                  structural_uscalar_dscalardx = structural_uscalarf_i*dscalardxf_i;
                  structural_vscalar_dscalardy = structural_vscalarf_j*dscalardyf_j;
                  structural_wscalar_dscalardz = 0.5 * (structural_wscalarf_kp1*dscalardzf_kp1 + structural_wscalarf_k*dscalardzf_k);
                }

#define AJOUT(somme,val) moy[somme] += val
              // moyennes
              AJOUT(U_MOY,uf);
              AJOUT(V_MOY,vf);
              AJOUT(W_MOY,we);
              // moyennes des carres
              AJOUT(UU_MOY,uf*uf);
              AJOUT(VV_MOY,vf*vf);
              AJOUT(WW_MOY,0.5*(wf_k*wf_k+wf_kp1*wf_kp1));
              // correlations
              AJOUT(UV_MOY,ue*ve);
              AJOUT(VW_MOY,ve*we);
              AJOUT(UW_MOY,ue*we);
              AJOUT(UP_MOY,ue*p);
              AJOUT(VP_MOY,ve*p);
              AJOUT(WP_MOY,we*p);
              // 3 eme ordre
              AJOUT(UUU_MOY,uf*uf*uf);
              AJOUT(VVV_MOY,vf*vf*vf);
              AJOUT(WWW_MOY,0.5*(wf_k*wf_k*wf_k+wf_kp1*wf_kp1*wf_kp1));
              // energie
              AJOUT(RHO_MOY,rho);
              AJOUT(T_MOY,t);
              AJOUT(UT_MOY,ue*t);
              AJOUT(VT_MOY,ve*t);
              AJOUT(WT_MOY,we*t);
              AJOUT(URHO_MOY,ue*rho);
              AJOUT(VRHO_MOY,ve*rho);
              AJOUT(WRHO_MOY,we*rho);
              AJOUT(RHOP_MOY,rho*p);
              AJOUT(TP_MOY,t*p);
              AJOUT(P_MOY,p);
              AJOUT(T2_MOY,t*t);
              AJOUT(T3_MOY,t*t*t);
              AJOUT(T4_MOY,t*t*t*t);
              AJOUT(RHO2_MOY,rho*rho);
              AJOUT(RHO3_MOY,rho*rho*rho);
              AJOUT(RHO4_MOY,rho*rho*rho*rho);
              AJOUT(UUW_MOY,ue*ue*ve);
              AJOUT(VVW_MOY,we*we*ve);
              AJOUT(U4_MOY,uf*uf*uf*uf);
              AJOUT(V4_MOY,vf*vf*vf*vf);
              AJOUT(W4_MOY,0.5*(wf_k*wf_k*wf_k*wf_k+wf_kp1*wf_kp1*wf_kp1*wf_kp1));
              AJOUT(MU_MOY,mu);
              AJOUT(LAMBDA_MOY,lambda);
              AJOUT(PP_MOY,p*p);
              AJOUT(WFACE_MOY,wf_k);
              AJOUT(NU_MOY,nu);
              AJOUT(UN_SUR_RHO_MOY,unsrho);
              AJOUT(NU2_MOY,nu*nu);
              AJOUT(MU2_MOY,mu*mu);
              AJOUT(NUTURB_MOY,nuturb_xx);
              AJOUT(LAMBDADTDZ_MOY,lambda*dtdz);
              AJOUT(KAPPATURB_MOY,kappaturb_x);
              AJOUT(RHOWFACE_MOY,rhof_k*wf_k);
              AJOUT(RHOTWFACE_MOY,rhof_k*tf_k*wf_k);
              AJOUT(RHOUU_MOY,rho*uf*uf);
              AJOUT(RHOVV_MOY,rho*vf*vf);
              AJOUT(RHOWW_MOY,rho*0.5*(wf_k*wf_k+wf_kp1*wf_kp1));
              AJOUT(RHOUV_MOY,rho*ue*ve);
              AJOUT(RHOVW_MOY,rho*ve*we);
              AJOUT(RHOUW_MOY,rho*ue*we);
              AJOUT(RHOUT_MOY,rho*ue*t);
              AJOUT(RHOVT_MOY,rho*ve*t);
              AJOUT(RHOWT_MOY,rho*we*t);
              AJOUT(NUTURB_XX_MOY,nuturb_xx);
              AJOUT(NUTURB_XY_MOY,nuturb_xy);
              AJOUT(NUTURB_XZ_MOY,nuturb_xz);
              AJOUT(NUTURB_YY_MOY,nuturb_yy);
              AJOUT(NUTURB_YZ_MOY,nuturb_yz);
              AJOUT(NUTURB_ZZ_MOY,nuturb_zz);
              AJOUT(KAPPATURB_X_MOY,kappaturb_x);
              AJOUT(KAPPATURB_Y_MOY,kappaturb_y);
              AJOUT(KAPPATURB_Z_MOY,kappaturb_z);
              AJOUT(NUTURB_XX_DUDX_MOY,nuturb_xx_dudx);
              AJOUT(NUTURB_XY_DUDY_MOY,nuturb_xy_dudy);
              AJOUT(NUTURB_XZ_DUDZ_MOY,nuturb_xz_dudz);
              AJOUT(NUTURB_XY_DVDX_MOY,nuturb_xy_dvdx);
              AJOUT(NUTURB_YY_DVDY_MOY,nuturb_yy_dvdy);
              AJOUT(NUTURB_YZ_DVDZ_MOY,nuturb_yz_dvdz);
              AJOUT(NUTURB_XZ_DWDX_MOY,nuturb_xz_dwdx);
              AJOUT(NUTURB_YZ_DWDY_MOY,nuturb_yz_dwdy);
              AJOUT(NUTURB_ZZ_DWDZ_MOY,nuturb_zz_dwdz);
              AJOUT(KAPPATURB_X_DSCALARDX_MOY,kappaturb_x_dscalardx);
              AJOUT(KAPPATURB_Y_DSCALARDY_MOY,kappaturb_y_dscalardy);
              AJOUT(KAPPATURB_Z_DSCALARDZ_MOY,kappaturb_z_dscalardz);
              AJOUT(STRUCTURAL_UU_MOY,structural_uu);
              AJOUT(STRUCTURAL_UV_MOY,structural_uv);
              AJOUT(STRUCTURAL_UW_MOY,structural_uw);
              AJOUT(STRUCTURAL_VV_MOY,structural_vv);
              AJOUT(STRUCTURAL_VW_MOY,structural_vw);
              AJOUT(STRUCTURAL_WW_MOY,structural_ww);
              AJOUT(STRUCTURAL_USCALAR_MOY,structural_uscalar);
              AJOUT(STRUCTURAL_VSCALAR_MOY,structural_vscalar);
              AJOUT(STRUCTURAL_WSCALAR_MOY,structural_wscalar);
              AJOUT(NUTURB_XX_DUDX_DUDX_MOY,nuturb_xx_dudx_dudx);
              AJOUT(NUTURB_XY_DUDY_DUDY_MOY,nuturb_xy_dudy_dudy);
              AJOUT(NUTURB_XZ_DUDZ_DUDZ_MOY,nuturb_xz_dudz_dudz);
              AJOUT(NUTURB_XY_DVDX_DVDX_MOY,nuturb_xy_dvdx_dvdx);
              AJOUT(NUTURB_YY_DVDY_DVDY_MOY,nuturb_yy_dvdy_dvdy);
              AJOUT(NUTURB_YZ_DVDZ_DVDZ_MOY,nuturb_yz_dvdz_dvdz);
              AJOUT(NUTURB_XZ_DWDX_DWDX_MOY,nuturb_xz_dwdx_dwdx);
              AJOUT(NUTURB_YZ_DWDY_DWDY_MOY,nuturb_yz_dwdy_dwdy);
              AJOUT(NUTURB_ZZ_DWDZ_DWDZ_MOY,nuturb_zz_dwdz_dwdz);
              AJOUT(NUTURB_XY_DVDX_DUDY_MOY,nuturb_xy_dvdx_dudy);
              AJOUT(NUTURB_XZ_DWDX_DUDZ_MOY,nuturb_xz_dwdx_dudz);
              AJOUT(NUTURB_YZ_DWDY_DVDZ_MOY,nuturb_yz_dwdy_dvdz);
              AJOUT(KAPPATURB_X_DSCALARDX_DSCALARDX_MOY,kappaturb_x_dscalardx_dscalardx);
              AJOUT(KAPPATURB_Y_DSCALARDY_DSCALARDY_MOY,kappaturb_y_dscalardy_dscalardy);
              AJOUT(KAPPATURB_Z_DSCALARDZ_DSCALARDZ_MOY,kappaturb_z_dscalardz_dscalardz);
              AJOUT(STRUCTURAL_UU_DUDX_MOY,structural_uu_dudx);
              AJOUT(STRUCTURAL_UV_DUDY_MOY,structural_uv_dudy);
              AJOUT(STRUCTURAL_UW_DUDZ_MOY,structural_uw_dudz);
              AJOUT(STRUCTURAL_UV_DVDX_MOY,structural_uv_dvdx);
              AJOUT(STRUCTURAL_VV_DVDY_MOY,structural_vv_dvdy);
              AJOUT(STRUCTURAL_VW_DVDZ_MOY,structural_vw_dvdz);
              AJOUT(STRUCTURAL_UW_DWDX_MOY,structural_uw_dwdx);
              AJOUT(STRUCTURAL_VW_DWDY_MOY,structural_vw_dwdy);
              AJOUT(STRUCTURAL_WW_DWDZ_MOY,structural_ww_dwdz);
              AJOUT(STRUCTURAL_USCALAR_DSCALARDX_MOY,structural_uscalar_dscalardx);
              AJOUT(STRUCTURAL_VSCALAR_DSCALARDY_MOY,structural_vscalar_dscalardy);
              AJOUT(STRUCTURAL_WSCALAR_DSCALARDZ_MOY,structural_wscalar_dscalardz);
              AJOUT(LAMBDADTDZ2_MOY,lambda*dtdz*lambda*dtdz);
#undef AJOUT
            }
        }
      // facteur 1./(ni*nj) car sommation de ni*nj valeurs sur des mailles de meme taille
      // facteur delta_z / taille_totale_en_z  car mailles non uniformes en z
      const int ni_tot = pression.get_splitting().get_grid_geometry().get_nb_elem_tot(DIRECTION_I);
      const int nj_tot = pression.get_splitting().get_grid_geometry().get_nb_elem_tot(DIRECTION_J);

      double facteur = 1./(double)(ni_tot * nj_tot);

      for (int i = 0; i < nval_; i++)
        tmp(k + offset, i) = moy[i] * facteur;
    }
  // Somme sur tous les processeurs:
  mp_sum_for_each_item(tmp);

  // Sur processeur 0, ajout de la contribution a l'integrale temporelle:
  if (Process::je_suis_maitre())
    {
      for (int i = 0; i < nval_; i++)
        {
          for (int k = 0; k < nktot; k++)
            {
              integrale_temporelle_[i][k] += tmp(k, i) * dt;
              moyenne_spatiale_instantanee_[i][k] = tmp(k, i);
            }
        }
      t_integration_ += dt;
    }
}

//modif AT 20/06/2013
// ce post traitement calcule les correlations en 1 points sur lenergie cinetique turbulente.

void Statistiques_dns_ijk::update_stat_k(const FixedVector<IJK_Field_double, 3>& vitesse,
                                         const IJK_Field_double& pression,
                                         //const IJK_Field_double &temperature,
                                         const IJK_Field_double& masse_vol,
                                         const IJK_Field_double& champ_mu,
                                         //const IJK_Field_double &champ_lambda,
                                         const double pression_thermodynamique,
                                         const double terme_source_acceleration,
                                         double dt)
{
  if (elem_coord_.size_array() == 0)
    {
      Cerr << "Erreur dans Statistiques_dns_ijk::update_stat: non initialise" << finl;
      Process::exit();
    }
  Cerr.setf(ios::scientific);
  Cerr.precision(20);
  // Copie des champs utiles des vitesses
  const IJK_Field_double& vitesse_i = vitesse[0];
  const IJK_Field_double& vitesse_j = vitesse[1];
  const IJK_Field_double& vitesse_k = vitesse[2];

  // copie des champs moyen existant
  const ArrOfDouble& vit_moy_i = vit_moy_[0];
  const ArrOfDouble& vit_moy_j = vit_moy_[1];
  const ArrOfDouble& vit_moy_k = vit_moy_[2];  // Ici c'est bien la vitesse aux face !!!
  const ArrOfDouble& rho_moy = rho_moy_;  // DD 16/10/2015
  const ArrOfDouble& nu_moy = nu_moy_;  // DD 16/10/2015

  // Nombre total de mailles en K
  const int nktot = pression.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);

  // Nombre local de mailles en K
  const int kmax = pression.nk();
  const int offsetk = pression.get_splitting().get_offset_local(DIRECTION_K);

  DoubleTab tmp(nktot, kval_);
  const int imax = pression.ni();
  const int jmax = pression.nj();

  const double unsurdx=1./dx_;
  const double unsurdy=1./dy_;

  // on calcule les proprietes du fluides sur la Cl :
  const double rho_kmin = calculer_rho_air(TCL_kmin_, constante_specifique_gaz_, pression_thermodynamique);
  const double rho_kmax = calculer_rho_air(TCL_kmax_, constante_specifique_gaz_, pression_thermodynamique);
  /*
  double lambda_kmin = calculer_lambda_air(TCL_kmin);
  double lambda_kmax = calculer_lambda_air(TCL_kmax);
  */

  const IJK_Splitting& split = pression.get_splitting();
  IJK_Field_double Ener_cin;


  Ener_cin.allocate(split, IJK_Splitting::ELEM, 1);
  Ener_cin.data()=0;

  IJK_Field_double Div_fluc_U;
  Div_fluc_U.allocate(split, IJK_Splitting::ELEM, 1);
  Div_fluc_U.data()=0;

  IJK_Field_double Div_tot_U;
  Div_tot_U.allocate(split, IJK_Splitting::ELEM, 1);
  Div_tot_U.data()=0;
  // on precalcule l'ener cinetique fluctuations de pression !!!    (le div u
  // pour la paroi
  int KMAX = kmax;
  if ((kmax+offsetk)!=nktot) // si je ne suis pas a la paroi je remplit le ghost
    KMAX++;
  int KMIN = 0;
  if (offsetk!=0)// si je ne suis pas a la paroi je remplit le ghost
    KMIN--;

  for (int k = KMIN; k < KMAX; k++)
    {
      const int kg=k+offsetk; // position en k sur la grille globale
      const double unsurdz=1./tab_dz_[kg]; // taille de la maille dans la direction z
      for (int j = -1; j < (jmax+1); j++)
        for (int i = -1; i < (imax+1); i++)
          {
            double uf_i = vitesse_i(i,j,k);
            double vf_j = vitesse_j(i,j,k);
            double wf_k = vitesse_k(i,j,k);

            double uf_ip1 = vitesse_i(i+1,j,k);
            double vf_jp1 = vitesse_j(i,j+1,k);
            double wf_kp1 = vitesse_k(i,j,k+1);

            // divergence de la valeur totale
            double duidx=(uf_ip1-uf_i)*unsurdx;
            double dujdy=(vf_jp1-vf_j)*unsurdy;
            double dukdz=(wf_kp1-wf_k)*unsurdz;
            double divu= duidx+dujdy+dukdz;
            Div_tot_U(i,j,k)= divu;

            // ici on calcule les fluctuations
            uf_i -= vit_moy_i[kg];
            vf_j -= vit_moy_j[kg];
            wf_k -= vit_moy_k[kg];

            uf_ip1 -= vit_moy_i[kg];
            vf_jp1 -= vit_moy_j[kg];
            wf_kp1 -= vit_moy_k[kg+1];

            double Ec=0.25*(uf_i*uf_i+uf_ip1*uf_ip1+vf_j*vf_j+vf_jp1*vf_jp1+wf_k*wf_k+wf_kp1*wf_kp1);
            Ener_cin(i,j,k)=Ec;

            // divergence des fluctuations
            duidx=(uf_ip1-uf_i)*unsurdx;
            dujdy=(vf_jp1-vf_j)*unsurdy;
            dukdz=(wf_kp1-wf_k)*unsurdz;
            divu= duidx+dujdy+dukdz;
            Div_fluc_U(i,j,k)=divu;
          }
    }
  for (int k = 0; k < kmax; k++)
    {
      // Calcul des moyennes spatiales sur le plan ij
      const int kg=k+offsetk; // position en k sur la grille globale
      const double unsurdz=1./tab_dz_[k+offsetk];
      double delta_m   = kg== 0 ? tab_dz_[kg] *0.5 : ( tab_dz_[kg] + tab_dz_[kg-1] ) * 0.5; // distance entre les centre de gravite avec gestion de la paroi !!!
      double delta_p   =  kg==(nktot-1) ? tab_dz_[kg] *0.5 : ( tab_dz_[kg] + tab_dz_[kg+1] ) * 0.5; // distance entre les centre de gravite

      double alpha      = delta_p/delta_m; // rapport des distance
      double unsalpha      = delta_m/delta_p; //
      double unsdpdm      = 1. / (  delta_p + delta_m );
      // On y stocke la somme de toutes les valeurs sur le plan ij, on divisera apres
      ArrOfDouble moy(kval_);
      for (int i = 0; i < kval_; i++)
        {
          moy[i] = 0.;
        }
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              // calcul des fluctuations de vitesses aux faces.
              const double uf_i   = vitesse_i(i,j,k)  -vit_moy_i[kg];
              const double uf_ip1 = vitesse_i(i+1,j,k)-vit_moy_i[kg];
              const double vf_j   = vitesse_j(i,j,k)  -vit_moy_j[kg];
              const double vf_jp1 = vitesse_j(i,j+1,k)-vit_moy_j[kg];
              const double wf_k   = vitesse_k(i,j,k)  -vit_moy_k[kg];
              const double wf_kp1 = vitesse_k(i,j,k+1)-vit_moy_k[kg+1];

              // Fluctuations de vitesses aux elems
              const double ue = (uf_ip1+uf_i)*0.5;
              const double ve = (vf_jp1+vf_j)*0.5;
              const double we = (wf_kp1+wf_k)*0.5;

              // vitesses aux elems
              const double uetot = 0.5*(vitesse_i(i+1,j,k) + vitesse_i(i,j,k));
              const double vetot = 0.5*(vitesse_j(i,j+1,k) + vitesse_j(i,j,k));
              const double wetot = 0.5*(vitesse_k(i,j,k+1) + vitesse_k(i,j,k));


              // calcul des grandeurs composites a partir de ces vitesse !!!
              const double duidx=(uf_ip1-uf_i)*unsurdx;
              const double dujdy=(vf_jp1-vf_j)*unsurdy;
              const double dukdz=(wf_kp1-wf_k)*unsurdz;
              const double divu= Div_fluc_U(i,j,k);

              // calcul des grandeurs scalaire !!
              const double unsurrho = 1./masse_vol(i,j,k);
              const double rho = masse_vol(i,j,k);
              const double mu = champ_mu(i,j,k);
              const double nu = mu * unsurrho;
              const double p = pression(i,j,k);
              const double unsurrho_moy = 1./rho_moy[kg];
              const double rho_fluc = (rho - rho_moy[kg]);
              const double nu_fluc = (nu - nu_moy[kg]);

              // derivee de Ec
              const double Ec         = Ener_cin(i,j,k);
              const double Ec_plus_i  = Ener_cin(i+1,j,k);
              const double Ec_plus_j  = Ener_cin(i,j+1,k);
              const double Ec_plus_k  = kg==(nktot-1)? 0. :Ener_cin(i,j,k+1);

              const double Ec_moins_i = Ener_cin(i-1,j,k);
              const double Ec_moins_j = Ener_cin(i,j-1,k);
              const double Ec_moins_k = kg==0 ? 0. : Ener_cin(i,j,k-1);

              double dEcdx = (Ec_plus_i  - Ec_moins_i) * 0.5 * unsurdx;
              double dEcdy = (Ec_plus_j  - Ec_moins_j) * 0.5 * unsurdy;
              double dEcdz = derivee_aniso(alpha, unsalpha, unsdpdm, Ec_moins_k, Ec, Ec_plus_k);

              // derivee de Rho
              const double rho_plus_i  = masse_vol(i+1,j,k);
              const double rho_plus_j  = masse_vol(i,j+1,k);
              const double rho_plus_k  = kg==(nktot-1) ? rho_kmax : masse_vol(i,j,k+1);

              const double rho_moins_i = masse_vol(i-1,j,k);
              const double rho_moins_j = masse_vol(i,j-1,k);
              const double rho_moins_k = kg==0 ? rho_kmin : masse_vol(i,j,k-1);

              const double drhodx = (rho_plus_i  - rho_moins_i ) * 0.5 * unsurdx;
              const double drhody = (rho_plus_j  - rho_moins_j ) * 0.5 * unsurdy;
              const double drhodz = derivee_aniso(alpha , unsalpha ,unsdpdm ,rho_moins_k, rho , rho_plus_k);

              // derivation direction x de Uy
              double Vf_im =  vitesse_j(i-1,j,k)-vit_moy_j[kg];
              double dujdx_0 = (vf_j - Vf_im ) * unsurdx;

              double nu_arrette_ij = 0.25*( ( champ_mu(i,j,k)      /masse_vol(i,j,k))
                                            + ( champ_mu(i-1,j,k)  /masse_vol(i-1,j,k))
                                            + ( champ_mu(i,j-1,k)  /masse_vol(i,j-1,k))
                                            + ( champ_mu(i-1,j-1,k)/masse_vol(i-1,j-1,k)) );
              double contrib_dujdx = dujdx_0*dujdx_0;

              // derivation direction x de de Uz
              double Wf_mi =  vitesse_k(i-1,j,k)-vit_moy_k[kg];
              double dukdx_0 = (wf_k - Wf_mi ) * unsurdx;

              double Wf_mipk = vitesse_k(i-1,j,k+1)-vit_moy_k[kg+1];
              double dukdx_1 = (wf_kp1 - Wf_mipk ) * unsurdx;

              double nu_face_i = 0.5 *( ( champ_mu(i,j,k)    /masse_vol(i,j,k))
                                        + ( champ_mu(i-1,j,k)/masse_vol(i-1,j,k)) );
              double contrib_dukdx = (dukdx_0* dukdx_0 + dukdx_1 * dukdx_1) *0.5;

              // derivation direction y de Ux
              double Uf_mj = vitesse_i(i,j-1,k)-vit_moy_i[kg];
              double duidy_0 = (uf_i  - Uf_mj ) * unsurdy;
              double contrib_duidy = duidy_0*duidy_0;

              //  derivation direction y de Uz
              double Wf_mj = vitesse_k(i,j-1,k)-vit_moy_k[kg];
              double dukdy_0 = (wf_k  - Wf_mj ) * unsurdy;

              double Wf_mj_pk = vitesse_k(i,j-1,k+1)-vit_moy_k[kg+1];
              double dukdy_1 = (wf_kp1  - Wf_mj_pk ) * unsurdy;

              double nu_face_j =  0.5 *( ( champ_mu(i,j,k)    /masse_vol(i,j,k))
                                         + ( champ_mu(i,j-1,k)/masse_vol(i,j-1,k)) );
              double contrib_dukdy = (dukdy_0* dukdy_0 + dukdy_1 * dukdy_1) *0.5;

              // derivation direction z de Ux
              double Uf_mk_0 = kg == 0 ? 0.         : vitesse_i(i,j,k-1)-vit_moy_i[kg-1];
              double Uf_pk_0 = kg == (nktot-1) ? 0. : vitesse_i(i,j,k+1)-vit_moy_i[kg+1];
              double duidz_0 = derivee_aniso(alpha, unsalpha, unsdpdm, Uf_mk_0, uf_i, Uf_pk_0);

              double Uf_mk_1 = kg == 0 ? 0.         : vitesse_i(i+1,j,k-1)-vit_moy_i[kg-1];
              double Uf_pk_1 = kg == (nktot-1) ? 0. : vitesse_i(i+1,j,k+1)-vit_moy_i[kg+1];
              double duidz_1 = derivee_aniso(alpha, unsalpha, unsdpdm, Uf_mk_1, uf_i, Uf_pk_1);

              double Uf_mk_0_tot = kg == 0 ? 0.           : vitesse_i(i,j,k-1);
              double Uf_pk_0_tot = kg == (nktot-1) ? 0.   : vitesse_i(i,j,k+1);
              double duidz_0_tot = derivee_aniso(alpha , unsalpha ,unsdpdm ,Uf_mk_0_tot, vitesse_i(i,j,k) , Uf_pk_0_tot);

              double Uf_mk_1_tot = kg == 0 ? 0.         : vitesse_i(i+1,j,k-1);
              double Uf_pk_1_tot = kg == (nktot-1) ? 0. : vitesse_i(i+1,j,k+1);
              double duidz_1_tot = derivee_aniso(alpha , unsalpha ,unsdpdm ,Uf_mk_1_tot, vitesse_i(i+1,j,k) , Uf_pk_1_tot);

              double contrib_duidz = duidz_0 * duidz_0;
              double contrib_duidz_tot = duidz_0 * duidz_0_tot;

              //  derivation direction z de Uy
              double Vf_mk= kg == 0 ? 0.            : vitesse_j(i,j,k-1)-vit_moy_j[kg-1];
              double Vf_pk= kg == (nktot-1) ? 0.    : vitesse_j(i,j,k+1)-vit_moy_j[kg+1];
              double dujdz_0 = derivee_aniso(alpha , unsalpha ,unsdpdm ,Vf_mk, vf_j , Vf_pk);

              double Vf_mk_1= kg == 0 ? 0.          : vitesse_j(i,j+1,k-1)-vit_moy_j[kg-1];
              double Vf_pk_1= kg == (nktot-1) ? 0.  : vitesse_j(i,j+1,k+1)-vit_moy_j[kg+1];
              double dujdz_1 = derivee_aniso(alpha , unsalpha ,unsdpdm ,Vf_mk_1, vf_jp1 , Vf_pk_1);

              double contrib_dujdz = dujdz_0 * dujdz_0;

              double contrib_div = (duidx*duidx + dujdy*dujdy + dukdz*dukdz);

              //  derivation direction z de Uz
              double dukdz_tot = ( vitesse_k(i,j,k+1) - vitesse_k(i,j,k) ) * unsurdz;

              // c'est parti pour les vitesses moyennes
              // dUx/dz
              double Umoyf_0  = vit_moy_i[kg];
              double Umoyf_1  = kg == (nktot-1) ? 0.  : vit_moy_i[kg+1];
              double Umoyf_00 = kg == 0         ? 0.  : vit_moy_i[kg-1];
              double dUmoydz  = derivee_aniso(alpha, unsalpha, unsdpdm, Umoyf_00, Umoyf_0, Umoyf_1);

              // dUz/dz
              double Wmoyf_0 = vit_moy_k[kg];
              double Wmoyf_1 = vit_moy_k[kg+1];
              double dWmoydz = (Wmoyf_1-Wmoyf_0)*unsurdz;

              // derivees centrees
              // de Uy ; ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j+1,k) aux mailles i-1 et i+1
              double Ve_mi = 0.5*( (vitesse_j(i-1,j,k) - vit_moy_j[kg]) + (vitesse_j(i-1,j+1,k) - vit_moy_j[kg]) );
              double Ve_pi = 0.5*( (vitesse_j(i+1,j,k) - vit_moy_j[kg]) + (vitesse_j(i+1,j+1,k) - vit_moy_j[kg]) );
              double dujdx = (Ve_pi  - Ve_mi) * 0.5 * unsurdx;

              // de Uz ; ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j,k+1) aux mailles i-1 et i+1
              double We_mi = 0.5*( (vitesse_k(i-1,j,k) - vit_moy_k[kg] ) + (vitesse_k(i-1,j,k+1) - vit_moy_k[kg+1]) );
              double We_pi = 0.5*( (vitesse_k(i+1,j,k) - vit_moy_k[kg] ) + (vitesse_k(i+1,j,k+1) - vit_moy_k[kg+1]) );
              double dukdx = (We_pi  - We_mi) * 0.5 * unsurdx;


              // derivation direction y
              // de Ux ; ATTENTION on veux calculer la moyenne entre (i,j,k) et (i+1,j,k) aux mailles j-1 et j+1
              double Ue_mj = 0.5*( (vitesse_i(i,j-1,k) - vit_moy_i[kg]) + (vitesse_i(i+1,j-1,k) - vit_moy_i[kg]) );
              double Ue_pj = 0.5*( (vitesse_i(i,j+1,k) - vit_moy_i[kg]) + (vitesse_i(i+1,j+1,k) - vit_moy_i[kg]) );
              double duidy = (Ue_pj  - Ue_mj) * 0.5 * unsurdy;

              // de Uz ; ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j,k+1) aux mailles j-1 et j+1
              double We_mj = 0.5*( (vitesse_k(i,j-1,k) - vit_moy_k[kg]) + (vitesse_k(i,j-1,k+1) - vit_moy_k[kg+1]) );
              double We_pj = 0.5*( (vitesse_k(i,j+1,k) - vit_moy_k[kg]) + (vitesse_k(i,j+1,k+1) - vit_moy_k[kg+1]) );
              double dukdy = (We_pj  - We_mj ) * 0.5 * unsurdy;

              /* *
                  *
                 *
                  * */

              double divtotu = Div_tot_U(i,j,k);
              double divergence_moy = (vit_moy_k[kg+1] - vit_moy_k[kg]) * unsurdz;

              // ------------------------------------------------------------
              // ============================================================

              const double contrib_divtot = (duidx*duidx + dujdy*dujdy + dukdz*dukdz_tot);
              const double duidz = 0.5*(duidz_0 + duidz_1);
              const double dujdz = 0.5*(dujdz_0 + dujdz_1);
              const double duidz_tot = 0.5*(duidz_0_tot + duidz_1_tot);

              // ------------------------------------------------------------
              // PRODUCTION
              // ------------------------------------------------------------
              // production incompressible
              // [<ux uy>] [d<Ux>/dy]

              // production thermique
              // [<uy uy>] [d<Uy>/dy]

              // ------------------------------------------------------------
              // ADVECTION
              // ------------------------------------------------------------
              // [<Uy>] [d<e>/dy] ; Ec=uiui/2

              // ------------------------------------------------------------
              // DIFFUSION TURBULENTE
              // ------------------------------------------------------------
              // d<Ec uy>/dy ; Ec=uiui/2
              const double Ec_we = Ec * we;

              // ------------------------------------------------------------
              // CORRELATION EC-DILATATION
              // ------------------------------------------------------------
              // <Ec duidx_i> ; Ec=uiui/2
              const double Ec_divu = Ec * divu;

              // ------------------------------------------------------------
              // CORRELATION PRESSION-DILATATION
              // ------------------------------------------------------------
              // <P/rho duidx_i>
              const double p_unsrho_divu = p * unsurrho * divu;

              // ------------------------------------------------------------
              // DIFFUSION PAR LA PRESSION
              // ------------------------------------------------------------
              // diffusion par la pression incompressible
              // [1/<rho>] <[d/dy](uy P)>
              const double p_we = p * we;

              // diffusion par la pression compressible, en d<rho>/dy
              // [<uy P>/<rho>2] [<d<rho>/dy]
              // meme terme que pour "diffusion par la pression incompressible"

              // diffusion par la pression compressible, en rho' et rho
              // <[d/dy]((uy P rho')/(<rho> rho))>
              const double rhofluc_unsrho_unsrhomoy_P_we = unsurrho_moy * unsurrho * we * p * rho_fluc;

              // diffusion par la pression compressible, en drho
              // <ui P/rho2 drho/dxi>
              const double somme_unsrho_unsrho_P_ui_drhodxi = unsurrho * unsurrho * p * (ue*drhodx + ve*drhody + we*drhodz);

              // ------------------------------------------------------------
              // DISSIPATION
              // ------------------------------------------------------------
              // dissipation incompressible, premiere partie
              // [<nu>] <dui/dx_j dui/dx_j>
              const double somme_duidxj_duidxj = contrib_div
                                                 + contrib_dujdz + contrib_duidz
                                                 + contrib_dukdy + contrib_duidy
                                                 + contrib_dujdx + contrib_dukdx;

              // dissipation incompressible, seconde partie
              // [<nu>] <dui/dx_j duj/dx_i>
              const double somme_duidxj_dujdxi = contrib_div
                                                 + 2*duidy*dujdx + 2*duidz*dukdx + 2*dujdz*dukdy;

              // dissipation compressible, premiere partie
              // <nu' dui/dx_j dUi/dx_j>
              const double somme_nufluc_dujdxi_dUidxj = contrib_divtot*nu_fluc
                                                        + contrib_dujdz*(nu_face_j - nu_moy[kg])     + contrib_duidz_tot*(nu_face_i - nu_moy[kg])
                                                        + contrib_dukdy*(nu_face_j - nu_moy[kg])     + contrib_duidy*(nu_arrette_ij - nu_moy[kg])
                                                        + contrib_dujdx*(nu_arrette_ij - nu_moy[kg]) + contrib_dukdx*(nu_face_i - nu_moy[kg]);

              // dissipation compressible, seconde partie
              // <nu' dui/dx_j dUj/dx_i>
              const double somme_nufluc_duidxj_dUjdxi = nu_fluc * ( contrib_divtot
                                                                    + 2*duidy*dujdx + duidz*dukdx + duidz_tot*dukdx + 2*dujdz*dukdy );

              // dissipation compressible, troisieme partie
              // <2/3 nu' dui/dx_i dUj/dx_j>
              const double deuxtiers_nu_divu_divU = (2./3.) * nu * divu * divtotu;

              // ------------------------------------------------------------
              // DIFFUSION VISQUEUSE
              // ------------------------------------------------------------
              // diffusion visqueuse incompressible, premiere partie
              // [<nu>] [d2<Ec>/dy2] ; Ec=uiui/2

              // diffusion visqueuse incompressible, seconde partie
              // [<nu>] <[d/dy](ui duy/dx_i)>
              const double somme_ui_dukdxi = ue*dukdx + ve*dukdy + we*dukdz;

              // diffusion visqueuse compressible, en d<nu>/dy, premiere partie
              // [d<nu>dy] [d<Ec>/dy] ; Ec=uiui/2

              // diffusion visqueuse compressible, en d<nu>/dy, seconde partie
              // [d<nu>dy] [d(ui duy/dx_i)/dy]
              // meme terme que pour "diffusion visqueuse incompressible, seconde partie"

              // diffusion visqueuse compressible, en nu' et nu, premiere partie
              // <[d/dy](nu' ui dUidy)>
              const double somme_nufluc_dEcdxk = nu_fluc * dEcdz;
              const double somme_nufluc_ui_dUimoydxk = nu_fluc * (ue*dUmoydz + we*dWmoydz);
              const double somme_nufluc_ui_dUidxk = somme_nufluc_dEcdxk + somme_nufluc_ui_dUimoydxk;

              // diffusion visqueuse compressible, en nu' et nu, seconde partie
              // <[d/dy](nu' ui dUydx_i)>
              const double somme_nufluc_ui_dUkdxi = nu_fluc * (ue*dukdx + ve*dukdy + we*dukdz_tot);

              // diffusion visqueuse compressible, en nu' et nu, troisieme partie
              // <2/3 [d/dy](nu uy dUydy)>
              const double somme_deuxtiers_nu_uk_divU = (2./3.) * nu * we * divtotu;

              // diffusion visqueuse compressible, en drho, premiere partie
              // <nu/rho ui dUi/dx_j drho/dx_j>
              const double somme_nusrho_dEcdxj_drhodxj = nu * unsurrho * (dEcdx*drhodx + dEcdy*drhody + dEcdz*drhodz);
              const double somme_nusrho_ui_dUimoydxk_drhodxk = nu * unsurrho * (ue*dUmoydz + we*dWmoydz) * drhodz;
              const double somme_nusrho_ui_dUidxj_drhodxj = somme_nusrho_dEcdxj_drhodxj + somme_nusrho_ui_dUimoydxk_drhodxk;

              // diffusion visqueuse compressible, en drho, seconde partie
              // <nu/rho ui dUj/dx_i drho/dx_j>
              const double ui_dUjdxi_drhodxj = ue * (drhodx*duidx     + drhody*dujdx + drhodz*dukdx);
              const double uj_dUjdxj_drhodxj = ve * (drhodx*duidy     + drhody*dujdy + drhodz*dukdy);
              const double uk_dUjdxk_drhodxj = we * (drhodx*duidz_tot + drhody*dujdz + drhodz*dukdz_tot);
              const double somme_nusrho_ui_dUjdxi_drhodxj = nu * unsurrho * (ui_dUjdxi_drhodxj + uj_dUjdxj_drhodxj + uk_dUjdxk_drhodxj);

              // diffusion visqueuse compressible, en drho, troisieme partie
              // <2/3 nu/rho ui dUj/dx_j drho/dx_i>
              const double somme_ui_divU_drhodxi = (ue * drhodx + ve * drhody + we * drhodz) * divtotu;
              const double somme_deuxtiers_nusrho_ui_divU_drhodxi = (2./3.) * nu * unsurrho * somme_ui_divU_drhodxi;

              // ------------------------------------------------------------
              // TERME NON CLASSE
              // ------------------------------------------------------------
              // <Uj Ec/rho drho/dx_j> ; Ec=uiui/2
              const double somme_unsrho_Ec_Uj_drhodxj = unsurrho * Ec * (uetot*drhodx + vetot*drhody + wetot*drhodz);

              // ------------------------------------------------------------
              // TERME SOURCE
              // ------------------------------------------------------------
              // <ux F/rho> ; M/O terme_source_acceleration = F/rho
              const double ux_terme_source_acceleration = ue * terme_source_acceleration;

              // ------------------------------------------------------------
              // ============================================================

#define AJOUT(somme,val) moy[somme] += val

              AJOUT(KW_MOY, Ec_we);
              AJOUT(CORRELATION_EC_DIVERGENCE_MOY, Ec_divu);
              AJOUT(CORRELATION_PRESSION_DIVERGENCE_MOY, p_unsrho_divu);
              AJOUT(P_WFLUC_MOY, p_we);
              AJOUT(DIFFUSION_PRESSION_RHOFLUC_ADERIVE_MOY, rhofluc_unsrho_unsrhomoy_P_we);
              AJOUT(DIFFUSION_PRESSION_DERIVRHO_MOY, - somme_unsrho_unsrho_P_ui_drhodxi);
              AJOUT(DISSIPATION_INCOMPRESSIBLE_UN_AFOISNUMOY_MOY, - somme_duidxj_duidxj);
              AJOUT(DISSIPATION_INCOMPRESSIBLE_DEUX_AFOISNUMOY_MOY, - somme_duidxj_dujdxi);
              AJOUT(DISSIPATION_COMPESSIBLE_UN_MOY, - somme_nufluc_dujdxi_dUidxj);
              AJOUT(DISSIPATION_COMPESSIBLE_DEUX_MOY, - somme_nufluc_duidxj_dUjdxi);
              AJOUT(DISSIPATION_COMPESSIBLE_TROIS_MOY, deuxtiers_nu_divu_divU);
              AJOUT(DIFFUSION_VISQUEUSE_INCOMPRESSIBLE_DEUX_MOY, somme_ui_dukdxi);
              AJOUT(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_UN_MOY, somme_nufluc_ui_dUidxk);
              AJOUT(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_DEUX_MOY, somme_nufluc_ui_dUkdxi);
              AJOUT(DIFFUSION_VISQUEUSE_NUFLUC_ADERIVE_TROIS_MOY, - somme_deuxtiers_nu_uk_divU);
              AJOUT(DIFFUSION_VISQUEUSE_DERIVRHO_UN_MOY, somme_nusrho_ui_dUidxj_drhodxj);
              AJOUT(DIFFUSION_VISQUEUSE_DERIVRHO_DEUX_MOY, somme_nusrho_ui_dUjdxi_drhodxj);
              AJOUT(DIFFUSION_VISQUEUSE_DERIVRHO_TROIS_MOY, - somme_deuxtiers_nusrho_ui_divU_drhodxi);
              AJOUT(TERME_NON_CLASSE_MOY, somme_unsrho_Ec_Uj_drhodxj);

              AJOUT(TERME_SOURCE_MOY, ux_terme_source_acceleration);

              // convergence des stats et stats annexes
              AJOUT(FLUC_U_MOY, ue);
              AJOUT(FLUC_V_MOY, ve);
              AJOUT(FLUC_W_MOY, we);
              AJOUT(FLUC_RHO_MOY, rho_fluc);
              AJOUT(FLUC_NU_MOY, nu_fluc);
              AJOUT(FLUC_DIV_U_MOY, divu);
              AJOUT(TOTAL_DIV_U_MOY, divtotu);
              AJOUT(MEAN_DIV_U_MOY, divergence_moy);

#undef AJOUT
            }
        }
      // facteur 1./(ni*nj) car sommation de ni*nj valeurs sur des mailles de meme taille
      // facteur delta_z / taille_totale_en_z  car mailles non uniformes en z
      const int ni_tot = pression.get_splitting().get_grid_geometry().get_nb_elem_tot(DIRECTION_I);
      const int nj_tot = pression.get_splitting().get_grid_geometry().get_nb_elem_tot(DIRECTION_J);

      double facteur = 1./(double)(ni_tot * nj_tot);

      for (int i = 0; i < kval_; i++)
        tmp(k + offsetk, i) = moy[i] * facteur;
    }
  // Somme sur tous les processeurs:
  mp_sum_for_each_item(tmp);

  // Sur processeur 0, ajout de la contribution a l'integrale temporelle:
  if (Process::je_suis_maitre())
    {
      for (int i = 0; i < kval_; i++)
        {
          for (int k = 0; k < nktot; k++)
            {
              integrale_k_[i][k] += tmp(k, i) * dt;
              moyenne_spatiale_ec_[i][k] = tmp(k, i);
            }
        }
      t_integration_k_ += dt;
    }
}


// If flag_valeur_instantanee==1 : write the last computed instantaneous value,
// if flag_valeur_instantanee==0 : write the time averaged value
void Statistiques_dns_ijk::postraiter(Sortie& os, int flag_valeur_instantanee) const
{
  if (!Process::je_suis_maitre())
    return;

  const int n = elem_coord_.size_array(); // nombre de points en K
  if (flag_valeur_instantanee == 0)
    {
      os << "# temps_integration " <<  t_integration_ << finl;
      os << "# Impression des moyennes temporelles" << finl;
    }
  else if (flag_valeur_instantanee == 1)
    {
      os << "# Impression des moyennes spatiales instantanee" << finl;
    }
  else
    {
      Cerr << "Erreur dans Statistiques_dns_ijk::postraiter: flag inconnu" << finl;
      Process::exit();
    }
  os << "# colonne 1 : coordonnee_K" << finl;
  for (int i = 0; i < nval_; i++)
    {
      os << "# colonne " << i+2 << " : " << noms_moyennes_[i] << finl;
    }
  for (int j = 0; j < n; j++)
    {
      char s[100];
      sprintf(s, "%16.16e ", elem_coord_[j]);
      os << s;
      for (int i = 0; i < nval_; i++)
        {
          double x;
          if (flag_valeur_instantanee == 0)
            x = integrale_temporelle_[i][j] / t_integration_;
          else
            x = moyenne_spatiale_instantanee_[i][j];
          sprintf(s, "%16.16e ", x);
          os << s;
        }
      os << finl;
    }
}

void Statistiques_dns_ijk::postraiter_k(Sortie& os, int flag_valeur_instantanee) const
{
  if (!Process::je_suis_maitre())
    return;

  const int n = elem_coord_.size_array(); // nombre de points en K
  if (flag_valeur_instantanee == 0)
    {
      os << "# temps_integration_k " <<  t_integration_k_ << finl;
      os << "# Impression des moyennes temporelles" << finl;
    }
  else if (flag_valeur_instantanee == 1)
    {
      os << "# Impression des moyennes spatiales instantanee" << finl;
    }
  else
    {
      Cerr << "Erreur dans Statistiques_dns_ijk::postraiter: flag inconnu" << finl;
      Process::exit();
    }

  os << "# colonne 1 : coordonnee_K" << finl;
  for (int i = 0; i < kval_; i++)
    {
      os << "# colonne " << i+2 << " : " << noms_k_[i] << finl;
    }
  for (int j = 0; j < n; j++)
    {
      char s[100];
      sprintf(s, "%16.16e ", elem_coord_[j]);
      os << s;
      for (int i = 0; i < kval_; i++)
        {
          double x;
          if ( flag_valeur_instantanee == 0)
            {
              if (t_integration_k_ !=0)
                {
                  x = integrale_k_[i][j] / t_integration_k_;
                }
              else
                {
                  x = integrale_k_[i][j];
                }
            }
          else
            {
              x = moyenne_spatiale_ec_[i][j];
            }
          sprintf(s, "%16.16e ", x);
          os << s;
        }
      os << finl;
    }
}

// Impression des donnees pour reprise
Sortie& Statistiques_dns_ijk::printOn(Sortie& os) const
{
  if (Process::je_suis_maitre())
    {
      os << "{\n";
      if (check_converge_)
        os << "  converge " << "\n";

      os  << "  temps_integration " << t_integration_ << "\n";
      os  << "  fields { " << "\n";
      for (int i = 0; i < nval_; i++)
        {
          os << noms_moyennes_[i] << " " << integrale_temporelle_[i] << finl;
        }
      os << "  }" << finl;
      completer_print(os); // Pour imprimer d'autres listes par exemple.
      //modif AT 20/06/2013
      if (check_converge_)
        {
          os << "  temps_integration_k " << t_integration_k_ << "\n";
          for (int i = 0; i < kval_; i++)
            {
              os << noms_k_[i] << " " << integrale_k_[i] << finl;
            }
        }
      os << "}" << finl;
      if (check_converge_)
        {
          Cout << "Ecriture des donnees statistiques d'energie pour reprise: t_integration_k_=" << t_integration_k_ << finl;
        }
      Cout << "Ecriture des donnees statistiques pour reprise: t_integration=" << t_integration_ << finl;
    }
  return os;
}

int Statistiques_dns_ijk::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  Cerr << "[Statistiques_dns_ijk] Motcle : " << mot << finl;
  Param param(que_suis_je()+Nom("_fields"));

  if (mot=="fields")
    {
      for (int i = 0; i < nval_; i++)
        {
          param.ajouter(noms_moyennes_[i], &integrale_temporelle_[i]);
        }
    }
  else
    {
      Cerr << "Motcle : " << mot <<  " inconnu" << finl;
      Cerr << "Keyword " << mot <<  " is not understood by Statistiques_dns_ijk::lire_motcle_non_standard(..)" << finl;
      Process::exit();
    }
  param.lire_avec_accolades(is);
  Cerr << "End of reading 'fields' in Statistiques_dns_ijk " << finl;
  return 1;
}

// Reprise des donnees stat dans un fichier reprise
// Attention, cette methode peut etre appelee avant initialize() !
Entree& Statistiques_dns_ijk::readOn(Entree& is)
{
  integrale_temporelle_.dimensionner(nval_);
  t_integration_ = 0.;
  check_converge_=0;//modif AT 20/06/2013
  Param param(que_suis_je());
  param.ajouter_flag("converge", &check_converge_);//modif AT 20/06/2013
  param.ajouter("temps_integration", &t_integration_);

  // Pour lire les nouveaux champs :
  param.ajouter_non_std("fields",(this));
  // ou pour lire quand meme les anciens directement sans accolade :
  for (int i = 0; i < nval_; i++)
    {
      param.ajouter(noms_moyennes_[i], &integrale_temporelle_[i]);
    }

  completer_read(param);
  integrale_k_.dimensionner(kval_);//modif AT 20/06/2013
  t_integration_k_=0.;
  param.ajouter("temps_integration_k", &t_integration_k_);

  for (int i = 0; i < kval_; i++)
    {
      param.ajouter(noms_k_[i], &integrale_k_[i]);
    }

  param.lire_avec_accolades(is);

  // // Seul le maitre conserve le tab d'integration temporelle lu en reprise :
  // // Les autres n'en ont pas besoin.
  // if (!Process::je_suis_maitre()) {
  //   integrale_temporelle_.reset();
  //   integrale_k_.reset();
  // }
  Cerr << "Check_converge " << check_converge_ << finl;
  Cout << "Reprise des donnees statistiques: t_integration=" << t_integration_ << finl;
  if ( check_converge_ )
    {
      Cout << "Reprise des donnees statistiques d'energie : t_integration_k_=" << t_integration_k_ << finl;
    }
  return is;
}

// bc_type : Le type de boundary_condition :
//    0    : Pour la vitesse, suppose une valeur nulle sur les bords.
//    1    : Pour le gradient_pression, on suppose que la valeur au bord
//               du gradient tangeant est la meme que celle au centre de
//               la premiere cellule.
double Statistiques_dns_ijk::face_to_cell_gradient(const IJK_Field_double& vitesse_i,
                                                   const IJK_Field_double& vitesse_j,
                                                   const IJK_Field_double& vitesse_k,
                                                   const int i, const int j, const int k,
                                                   const double dz,
                                                   double& duidx,
                                                   double& dujdx,
                                                   double& dukdx,
                                                   double& duidy,
                                                   double& dujdy,
                                                   double& dukdy,
                                                   double& duidz,
                                                   double& dujdz,
                                                   double& dukdz,
                                                   const bool on_the_first_cell,
                                                   const bool on_the_last_cell,
                                                   const int bc_type) const
{

  bool perio_k=vitesse_k.get_splitting().get_grid_geometry().get_periodic_flag(DIRECTION_K);
  // ******************************** //
  // derivation direction x
  // de Ux
  duidx = (vitesse_i(i+1, j, k) - vitesse_i(i,j,k))/dx_;

  // de Uy !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j+1,k) aux mailles i-1 et i+1
  double Ve_mi = ( vitesse_j(i-1,j,k) + vitesse_j(i-1,j+1,k) )*0.5;
  double Ve_pi = ( vitesse_j(i+1,j,k) + vitesse_j(i+1,j+1,k) )*0.5;
  dujdx = (Ve_pi  - Ve_mi ) / (2*dx_);

  // de Uz !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j,k+1) aux mailles i-1 et i+1
  double We_mi = ( vitesse_k(i-1,j,k) + vitesse_k(i-1,j,k+1) )*0.5;
  double We_pi = ( vitesse_k(i+1,j,k) + vitesse_k(i+1,j,k+1) )*0.5;
  dukdx = (We_pi  - We_mi ) / (2*dx_);

  // ******************************** //
  // derivation direction y
  // de Ux !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i+1,j,k) aux mailles j-1 et j+1
  double Ue_mj = ( vitesse_i(i,j-1,k) + vitesse_i(i+1,j-1,k) )*0.5;
  double Ue_pj = ( vitesse_i(i,j+1,k) + vitesse_i(i+1,j+1,k) )*0.5;
  duidy = (Ue_pj  - Ue_mj ) / (2*dy_);

  // de Uy
  dujdy = (vitesse_j(i, j+1, k) - vitesse_j(i,j,k))/dy_;

  // de Uz !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j,k+1) aux mailles j-1 et j+1
  double We_mj = ( vitesse_k(i,j-1,k) + vitesse_k(i,j-1,k+1) )*0.5;
  double We_pj = ( vitesse_k(i,j+1,k) + vitesse_k(i,j+1,k+1) )*0.5;
  dukdy = (We_pj  - We_mj ) / (2*dy_);

  // ******************************** //
  // derivation direction z
  // Si on est sur un bord, on utilise l'info que la vitesse y est nulle.
  // Cette info est a dz/2 donc on utilise la formule centree d'ordre 2 pour pas variable :
  // Formule centree (ordre 2) pour pas variable dans le domaine :
  // grad[1:-1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
  //
  if (on_the_first_cell)
    {
      if (!perio_k)
        {
          // de Ux
          double Ue_ck = ( vitesse_i(i,j,k  ) + vitesse_i(i+1,j,k  ) )*0.5;
          double Ue_pk = ( vitesse_i(i,j,k+1) + vitesse_i(i+1,j,k+1) )*0.5;
          double Ue_mk = (bc_type ? Ue_ck+(Ue_ck-Ue_pk)/8. : 0.);
          // Formule decentree (ordre 2) pour pas variable sur le bord gauche :
          duidz = (- 4*Ue_mk + 3*Ue_ck +Ue_pk  ) / (3*dz);

          // de Uy
          double Ve_ck = ( vitesse_j(i,j,k  ) + vitesse_j(i,j+1,k  ) )*0.5;
          double Ve_pk = ( vitesse_j(i,j,k+1) + vitesse_j(i,j+1,k+1) )*0.5;
          double Ve_mk =  (bc_type ? Ve_ck+(Ve_ck-Ve_pk)/8. : 0.);
          dujdz = (- 4*Ve_mk + 3*Ve_ck +Ve_pk  ) / (3*dz) ;
        }
    }
  else if (on_the_last_cell)
    {
      if (!perio_k)
        {
          // de Ux

          double Ue_mk = ( vitesse_i(i,j,k-1) + vitesse_i(i+1,j,k-1) )*0.5;
          double Ue_ck = ( vitesse_i(i,j,k  ) + vitesse_i(i+1,j,k  ) )*0.5;
          double Ue_pk = (bc_type ? Ue_ck+(Ue_ck-Ue_mk)/8. : 0.);
          // Formule decentree (ordre 2) pour pas variable sur le bord droit :
          duidz = (- Ue_mk - 3*Ue_ck +4*Ue_pk ) / (3*dz);

          // de Uy
          double Ve_mk = ( vitesse_j(i,j,k-1) + vitesse_j(i,j+1,k-1) )*0.5;
          double Ve_ck = ( vitesse_j(i,j,k  ) + vitesse_j(i,j+1,k  ) )*0.5;
          double Ve_pk = (bc_type ? Ve_ck+(Ve_ck-Ve_mk)/8. : 0.);
          dujdz = (- Ve_mk - 3*Ve_ck +4*Ve_pk ) / (3*dz) ;
        }
    }
  else
    {
      // For any interior cell... Formule centree pour pas cste

      // de Ux !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i+1,j,k) aux mailles k-1 et k+1
      double Ue_mk = ( vitesse_i(i,j,k-1) + vitesse_i(i+1,j,k-1) )*0.5;
      double Ue_pk = ( vitesse_i(i,j,k+1) + vitesse_i(i+1,j,k+1) )*0.5;
      duidz = (Ue_pk  - Ue_mk ) / (2*dz);

      // de Uy !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j+1,k) aux mailles k-1 et k+1
      double Ve_mk = ( vitesse_j(i,j,k-1) + vitesse_j(i,j+1,k-1) )*0.5;
      double Ve_pk = ( vitesse_j(i,j,k+1) + vitesse_j(i,j+1,k+1) )*0.5;
      dujdz = (Ve_pk  - Ve_mk ) / (2*dz) ;
    }

  // de Uz
  dukdz = (vitesse_k(i, j, k+1) - vitesse_k(i,j,k))/dz;

  // ******************************** //
  double dissip = duidx*duidx + dujdx*dujdx + dukdx*dukdx
                  + duidy*duidy + dujdy*dujdy + dukdy*dukdy
                  + duidz*duidz + dujdz*dujdz + dukdz*dukdz;

  return dissip;
}

// Calcul le gradient de U aux cellules a partir de la vitesse aux faces
void Statistiques_dns_ijk::compute_and_store_gradU_cell(const IJK_Field_double& vitesse_i,
                                                        const IJK_Field_double& vitesse_j,
                                                        const IJK_Field_double& vitesse_k,
                                                        /* Et les champs en sortie */
                                                        IJK_Field_double& dudx, IJK_Field_double& dvdy, IJK_Field_double& dwdx,
                                                        IJK_Field_double& dudz, IJK_Field_double& dvdz, IJK_Field_double& dwdz)
{
  const IJK_Splitting& splitting = vitesse_i.get_splitting();

  // Nombre total de mailles en K
  const int nktot = splitting.get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  // Nombre local de mailles :
  const int imax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 0);
  const int jmax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 1);
  const int kmax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 2);

  const int offset = splitting.get_offset_local(DIRECTION_K);

  for (int k = 0; k < kmax; k++)
    {
      const double dz = tab_dz_[k+offset];
      bool on_the_first_cell = false;
      bool on_the_last_cell = false;
      if (k + offset == 0)
        on_the_first_cell = true;
      if (k + offset == nktot - 1)
        on_the_last_cell = true;
      for (int j = 0; j < jmax; j++)
        {
          for (int i = 0; i < imax; i++)
            {
              // ******************************** //
              // derivation direction x
              // de Ux
              dudx(i,j,k) = (vitesse_i(i+1, j, k) - vitesse_i(i,j,k))/dx_;

              // de Uz !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j,k+1) aux mailles i-1 et i+1
              double We_mi = ( vitesse_k(i-1,j,k) + vitesse_k(i-1,j,k+1) )*0.5;
              double We_pi = ( vitesse_k(i+1,j,k) + vitesse_k(i+1,j,k+1) )*0.5;
              dwdx(i,j,k) = (We_pi  - We_mi ) / (2*dx_);

              // ******************************** //
              // derivation direction y
              // de Uy
              dvdy(i,j,k) = (vitesse_j(i, j+1, k) - vitesse_j(i,j,k))/dy_;

              // ******************************** //
              // derivation direction z
              // Si on est sur un bord, on utilise l'info que la vitesse y est nulle.
              // Cette info est a dz/2 donc on utilise la formule centree d'ordre 2 pour pas variable :
              // Formule centree (ordre 2) pour pas variable dans le domaine :
              // grad[1:-1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
              //
              if (on_the_first_cell)
                {
                  // de Ux
                  double Ue_mk = 0.;
                  double Ue_ck = ( vitesse_i(i,j,k  ) + vitesse_i(i+1,j,k  ) )*0.5;
                  double Ue_pk = ( vitesse_i(i,j,k+1) + vitesse_i(i+1,j,k+1) )*0.5;
                  // Formule decentree (ordre 2) pour pas variable sur le bord gauche :
                  dudz(i,j,k) = (- 4*Ue_mk + 3*Ue_ck +Ue_pk  ) / (3*dz);

                  // de Uy
                  double Ve_mk = 0.;
                  double Ve_ck = ( vitesse_j(i,j,k  ) + vitesse_j(i,j+1,k  ) )*0.5;
                  double Ve_pk = ( vitesse_j(i,j,k+1) + vitesse_j(i,j+1,k+1) )*0.5;
                  dvdz(i,j,k) = (- 4*Ve_mk + 3*Ve_ck +Ve_pk  ) / (3*dz) ;
                }
              else if (on_the_last_cell)
                {
                  // de Ux
                  double Ue_mk = ( vitesse_i(i,j,k-1) + vitesse_i(i+1,j,k-1) )*0.5;
                  double Ue_ck = ( vitesse_i(i,j,k  ) + vitesse_i(i+1,j,k  ) )*0.5;
                  double Ue_pk = 0.;
                  // Formule decentree (ordre 2) pour pas variable sur le bord droit :
                  dudz(i,j,k) = (- Ue_mk - 3*Ue_ck +4*Ue_pk ) / (3*dz);

                  // de Uy
                  double Ve_mk = ( vitesse_j(i,j,k-1) + vitesse_j(i,j+1,k-1) )*0.5;
                  double Ve_ck = ( vitesse_j(i,j,k  ) + vitesse_j(i,j+1,k  ) )*0.5;
                  double Ve_pk = 0.;
                  dvdz(i,j,k) = (- Ve_mk - 3*Ve_ck +4*Ve_pk ) / (3*dz) ;
                }
              else
                {
                  // For any interior cell... Formule centree pour pas cste

                  // de Ux !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i+1,j,k) aux mailles k-1 et k+1
                  double Ue_mk = ( vitesse_i(i,j,k-1) + vitesse_i(i+1,j,k-1) )*0.5;
                  double Ue_pk = ( vitesse_i(i,j,k+1) + vitesse_i(i+1,j,k+1) )*0.5;
                  dudz(i,j,k) = (Ue_pk  - Ue_mk ) / (2*dz);

                  // de Uy !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j+1,k) aux mailles k-1 et k+1
                  double Ve_mk = ( vitesse_j(i,j,k-1) + vitesse_j(i,j+1,k-1) )*0.5;
                  double Ve_pk = ( vitesse_j(i,j,k+1) + vitesse_j(i,j+1,k+1) )*0.5;
                  dvdz(i,j,k) = (Ve_pk  - Ve_mk ) / (2*dz) ;
                }

              // de Uz
              dwdz(i,j,k) = (vitesse_k(i, j, k+1) - vitesse_k(i,j,k))/dz;
            }
        }
    }

  // Mise a jour des espaces virtuels :
  dudx.echange_espace_virtuel(dudx.ghost());
  dvdy.echange_espace_virtuel(dvdy.ghost());
  dwdx.echange_espace_virtuel(dwdx.ghost());
  dudz.echange_espace_virtuel(dudz.ghost());
  dvdz.echange_espace_virtuel(dvdz.ghost());
  dwdz.echange_espace_virtuel(dwdz.ghost());
  return;
}

IJK_Field_double Statistiques_dns_ijk::compute_and_store_scalar_product_face_to_face(
  const IJK_Field_double& v1_i, const IJK_Field_double& v1_j, const IJK_Field_double& v1_k,
  const IJK_Field_double& v2_i,const IJK_Field_double& v2_j,const IJK_Field_double& v2_k)
{
  // champ scalaire qui accueuille le resultat
  IJK_Field_double v1v2 = v1_i;
  v1v2.data()=0;

  const IJK_Splitting& splitting = v1_i.get_splitting();
  // Nombre total de mailles en K
  //const int nktot = splitting.get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  // Nombre local de mailles :
  const int imax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 0);
  const int jmax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 1);
  const int kmax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 2);

  //const int offset = splitting.get_offset_local(DIRECTION_K);

  for (int k = 0; k < kmax; k++)
    for (int j = 0; j < jmax; j++)
      for (int i = 0; i < imax; i++)
        {
          v1v2(i,j,k) += v1_i(i,j,k)*v2_i(i,j,k);
          v1v2(i,j,k) += v1_j(i,j,k)*v2_j(i,j,k);
          v1v2(i,j,k) += v1_k(i,j,k)*v2_k(i,j,k);
        }
  return v1v2;
}

// Fonction pour calculer le gradient d'un champs aux elems
// Le champ retourne est aussi aux elems.
// On ne calcule pas de derivee selon z (pas necessaire, evite de se poser des questions de bords...)
// Le seul choix sur les inputs : au lieu de dwdx on peut prendre dwdy.
void Statistiques_dns_ijk::cell_to_cell_gradient(const int i, const int j, const int k,
                                                 const IJK_Field_double& dudx, const IJK_Field_double& dvdy, const IJK_Field_double& dwdx,
                                                 const IJK_Field_double& dudz, const IJK_Field_double& dvdz, const IJK_Field_double& dwdz,
                                                 /* Et les outputs en ref aussi!! */
                                                 double& ddudxy, double& ddudxz, double& ddudyz,
                                                 double& ddvdxy, double& ddvdxz, double& ddvdyz,
                                                 double& ddwdxy, double& ddwdxz, double& ddwdyz) const
{
  // ******************************** //
  // derivation direction x
  ddudxz = (dudz(i+1,j,k) - dudz(i-1,j,k) ) / (2*dx_);
  ddvdxy = (dvdy(i+1,j,k) - dvdy(i-1,j,k) ) / (2*dx_);
  ddvdxz = (dvdz(i+1,j,k) - dvdz(i-1,j,k) ) / (2*dx_);
  ddwdxz = (dwdz(i+1,j,k) - dwdz(i-1,j,k) ) / (2*dx_);

  // ******************************** //
  // derivation direction y
  ddudxy = (dudx(i,j+1,k) - dudx(i,j-1,k) ) / (2*dy_);
  ddudyz = (dudz(i,j+1,k) - dudz(i,j-1,k) ) / (2*dy_);
  ddvdyz = (dvdz(i,j+1,k) - dvdz(i,j-1,k) ) / (2*dy_);
  ddwdxy = (dwdx(i,j+1,k) - dwdx(i,j-1,k) ) / (2*dy_);
  ddwdyz = (dwdz(i,j+1,k) - dwdz(i,j-1,k) ) / (2*dy_);
  return;
}

// Calcul le gradient et la derivee seconde de la vitesse
double Statistiques_dns_ijk::calculer_gradients_vitesse(const IJK_Field_double& vitesse_i,
                                                        const IJK_Field_double& vitesse_j,
                                                        const IJK_Field_double& vitesse_k,
                                                        const int i, const int j, const int k,
                                                        const double dz,
                                                        double& duidx, double& dujdx, double& dukdx,
                                                        double& duidy, double& dujdy, double& dukdy,
                                                        double& duidz, double& dujdz, double& dukdz,
                                                        double& dduidxx, double& ddujdxx, double& ddukdxx,
                                                        double& dduidyy, double& ddujdyy, double& ddukdyy,
                                                        double& dduidzz, double& ddujdzz, double& ddukdzz,
                                                        const bool on_the_first_cell,
                                                        const bool on_the_last_cell) const
{
  const double dx2 = dx_*dx_;
  const double dy2 = dy_*dy_;
  const double dz2 = dz*dz;


  // ******************************** //
  // derivation direction x
  // de Ux
  duidx = (vitesse_i(i+1, j, k) - vitesse_i(i,j,k))/dx_;
  dduidxx = (vitesse_i(i+2, j, k) - vitesse_i(i+1,j,k) - vitesse_i(i, j, k) + vitesse_i(i-1,j,k))/(2*dx2);

  // de Uy
  double Ve_mi  = ( vitesse_j(i-1,j,k) + vitesse_j(i-1,j+1,k) )*0.5;
  double Ve_ci  = ( vitesse_j(i  ,j,k) + vitesse_j(i  ,j+1,k) )*0.5;
  double Ve_pi  = ( vitesse_j(i+1,j,k) + vitesse_j(i+1,j+1,k) )*0.5;
  dujdx = (Ve_pi  - Ve_mi ) / (2*dx_);
  ddujdxx = (Ve_pi - 2*Ve_ci + Ve_mi ) / dx2;

  // de Uz
  double We_mi  = ( vitesse_k(i-1,j,k) + vitesse_k(i-1,j,k+1) )*0.5;
  double We_ci  = ( vitesse_k(i  ,j,k) + vitesse_k(i  ,j,k+1) )*0.5;
  double We_pi  = ( vitesse_k(i+1,j,k) + vitesse_k(i+1,j,k+1) )*0.5;
  dukdx = (We_pi  - We_mi ) / (2*dx_);
  ddukdxx = (We_pi - 2*We_ci + We_mi ) / dx2;

  // ******************************** //
  // derivation direction y
  // de Ux
  double Ue_mj  = ( vitesse_i(i,j-1,k) + vitesse_i(i+1,j-1,k) )*0.5;
  double Ue_cj  = ( vitesse_i(i,j  ,k) + vitesse_i(i+1,j  ,k) )*0.5;
  double Ue_pj  = ( vitesse_i(i,j+1,k) + vitesse_i(i+1,j+1,k) )*0.5;
  duidy = (Ue_pj  - Ue_mj ) / (2*dy_);
  dduidyy = (Ue_pj  - 2*Ue_cj + Ue_mj ) / dy2;

  // de Uy
  dujdy = (vitesse_j(i, j+1, k) - vitesse_j(i,j,k))/dy_;
  ddujdyy = (vitesse_j(i, j+2, k) - vitesse_j(i,j+1,k) - vitesse_j(i, j, k) + vitesse_j(i,j-1,k))/(2*dy2);

  // de Uz
  double We_mj  = ( vitesse_k(i,j-1,k) + vitesse_k(i,j-1,k+1) )*0.5;
  double We_cj  = ( vitesse_k(i,j  ,k) + vitesse_k(i,j  ,k+1) )*0.5;
  double We_pj  = ( vitesse_k(i,j+1,k) + vitesse_k(i,j+1,k+1) )*0.5;
  dukdy = (We_pj  - We_mj ) / (2*dy_);
  ddukdyy = (We_pj  - 2*We_cj + We_mj ) / dy2;

  // ******************************** //
  // derivation direction z
  // Si on est sur un bord, on utilise l'info que la vitesse y est nulle.
  // Cette info est a dz/2 donc on utilise la formule centree d'ordre 2 pour pas variable :
  // Formule centree (ordre 2) pour pas variable dans le domaine :
  // grad[1:-1] = (h1/h2*u_pl - h2/h1*u_m + (h2**2-h1**2)/(h1*h2)*u_c) / (h1+h2)
  //
  if (on_the_first_cell)
    {
      // de Ux
      double Ue_mk = 0.;
      double Ue_ck = ( vitesse_i(i,j,k  ) + vitesse_i(i+1,j,k  ) )*0.5;
      double Ue_pk = ( vitesse_i(i,j,k+1) + vitesse_i(i+1,j,k+1) )*0.5;
      // Formule decentree (ordre 2) pour pas variable sur le bord gauche :
      duidz = (- 4*Ue_mk + 3*Ue_ck +Ue_pk  ) / (3*dz);
      dduidzz = (4*Ue_pk - 12*Ue_ck + 8*Ue_mk ) / (3*dz2);

      // de Uy
      double Ve_mk = 0.;
      double Ve_ck = ( vitesse_j(i,j,k  ) + vitesse_j(i,j+1,k  ) )*0.5;
      double Ve_pk = ( vitesse_j(i,j,k+1) + vitesse_j(i,j+1,k+1) )*0.5;
      dujdz = (- 4*Ve_mk + 3*Ve_ck +Ve_pk  ) / (3*dz) ;
      ddujdzz = (4*Ve_pk -12*Ve_ck + 8*Ve_mk ) / (3*dz2);
    }
  else if (on_the_last_cell)
    {
      // de Ux
      double Ue_mk = ( vitesse_i(i,j,k-1) + vitesse_i(i+1,j,k-1) )*0.5;
      double Ue_ck = ( vitesse_i(i,j,k  ) + vitesse_i(i+1,j,k  ) )*0.5;
      double Ue_pk = 0.;
      // Formule decentree (ordre 2) pour pas variable sur le bord droit :
      duidz = (- Ue_mk - 3*Ue_ck +4*Ue_pk ) / (3*dz);
      dduidzz = (8*Ue_pk - 12*Ue_ck + 4*Ue_mk ) / (3*dz2);

      // de Uy
      double Ve_mk = ( vitesse_j(i,j,k-1) + vitesse_j(i,j+1,k-1) )*0.5;
      double Ve_ck = ( vitesse_j(i,j,k  ) + vitesse_j(i,j+1,k  ) )*0.5;
      double Ve_pk = 0.;
      dujdz = (- Ve_mk - 3*Ve_ck +4*Ve_pk ) / (3*dz) ;
      ddujdzz = (8*Ve_pk - 12*Ve_ck + 4*Ve_mk ) / (3*dz2);
    }
  else
    {
      // For any interior cell... Formule centree pour pas cste

      // de Ux !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i+1,j,k) aux mailles k-1 et k+1
      double Ue_mk  = ( vitesse_i(i,j,k-1) + vitesse_i(i+1,j,k-1) )*0.5;
      double Ue_ck  = ( vitesse_i(i,j,k  ) + vitesse_i(i+1,j,k  ) )*0.5;
      double Ue_pk  = ( vitesse_i(i,j,k+1) + vitesse_i(i+1,j,k+1) )*0.5;
      duidz = (Ue_pk  - Ue_mk ) / (2*dz);
      dduidzz = (Ue_pk - 2*Ue_ck + Ue_mk ) / dz2;

      // de Uy !!!!!! ATTENTION on veux calculer la moyenne entre (i,j,k) et (i,j+1,k) aux mailles k-1 et k+1
      double Ve_mk  = ( vitesse_j(i,j,k-1) + vitesse_j(i,j+1,k-1) )*0.5;
      double Ve_ck  = ( vitesse_j(i,j,k  ) + vitesse_j(i,j+1,k  ) )*0.5;
      double Ve_pk  = ( vitesse_j(i,j,k+1) + vitesse_j(i,j+1,k+1) )*0.5;
      dujdz = (Ve_pk  - Ve_mk ) / (2*dz) ;
      ddujdzz = (Ve_pk - 2*Ve_ck + Ve_mk ) / dz2;
    }

  // de Uz
  dukdz = (vitesse_k(i, j, k+1) - vitesse_k(i,j,k))/dz;

  if (on_the_first_cell)
    {
      // left sided second derivative (first order) :
      ddukdzz = (vitesse_k(i, j, k+2) - 2*vitesse_k(i,j,k+1) + vitesse_k(i, j, k) )/(dz2);
    }
  else if (on_the_last_cell)
    {
      // right sided second derivative (first order) :
      ddukdzz = (vitesse_k(i,j,k+1) -2* vitesse_k(i, j, k) + vitesse_k(i,j,k-1) )/(dz2);
    }
  else
    {
      ddukdzz = (vitesse_k(i, j, k+2) - vitesse_k(i,j,k+1) - vitesse_k(i, j, k) + vitesse_k(i,j,k-1))/(2*dz2);
    }

  // ******************************** //
  /*
   * pseudo_dissip = nu (djui+djui)
   * ET ON MET PAS LE nu
  */
  double pseudo_dissip = duidx*duidx + dujdx*dujdx + dukdx*dukdx
                         + duidy*duidy + dujdy*dujdy + dukdy*dukdy
                         + duidz*duidz + dujdz*dujdz + dukdz*dukdz;
  return (pseudo_dissip);
}

double Statistiques_dns_ijk::calculer_vraie_dissipation(const double& pseudo_dissip,
                                                        const double& duidx, const double& duidy, const double& duidz,
                                                        const double& dujdx, const double& dujdy, const double& dujdz,
                                                        const double& dukdx, const double& dukdy, const double& dukdz) const
{
  /*
   * true_dissip = 2 nu sij sij = 2 nu ( 1/2(djui+diuj) 1/2(djui+diuj) )
   *             = pseudo_dissip + nu djui diuj
   * CONLCUSION : pas besoin de mettre de facteur 0.5 ou 0.25 ou quoi
   *              ou qu'est-ce.
  */
  double true_dissip = pseudo_dissip
                       + duidx*duidx + duidy*dujdx + duidz*dukdx
                       + dujdx*duidy + dujdy*dujdy + dujdz*dukdy
                       + dukdx*duidz + dukdy*dujdz + dukdz*dukdz;
  return(true_dissip);
}

double Statistiques_dns_ijk::calculer_produit_scalaire_faces_to_center(
  const IJK_Field_double& ui, const IJK_Field_double& uj, const IJK_Field_double& uk,
  const IJK_Field_double& vi, const IJK_Field_double& vj, const IJK_Field_double& vk,
  const int i, const int j, const int k)
{
  // Retourne le resultat u.v pour u et v des vecteurs aux faces.
  // Le resultat se situe au centre de l'element
  double produit_scalaire = (
                              (ui(i,j,k)*vi(i,j,k) + ui(i+1,j,k)*vi(i+1,j,k))*0.5 +
                              (uk(i,j,k)*vk(i,j,k) + uk(i,j,k+1)*vj(i,j,k+1))*0.5 +
                              (uk(i,j,k)*vk(i,j,k) + uk(i,j,k+1)*vj(i,j,k+1))*0.5
                            );

  return produit_scalaire;
}

void Statistiques_dns_ijk::compute_vecA_minus_vecB_in_vecA(FixedVector<IJK_Field_double, 3>& vecA, const FixedVector<IJK_Field_double, 3>& vecB)
{

  const IJK_Splitting& splitting = vecA.get_splitting();
  // Nombre total de mailles en K
  //const int nktot = splitting.get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  // Nombre local de mailles :
  const int imax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 0);
  const int jmax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 1);
  const int kmax = splitting.get_nb_items_local(IJK_Splitting::ELEM, 2);

  //const int offset = splitting.get_offset_local(DIRECTION_K);

  for (int k = 0; k < kmax; k++)
    for (int j = 0; j < jmax; j++)
      for (int i = 0; i < imax; i++)
        for (int dir = 0; dir < 3; ++dir)
          {
            vecA[dir](i,j,k) -= vecB[dir](i,j,k);
          }
}

