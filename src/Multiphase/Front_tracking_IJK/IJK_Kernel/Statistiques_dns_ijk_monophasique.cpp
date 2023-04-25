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

#include <Statistiques_dns_ijk_monophasique.h>
#include <IJK_Grid_Geometry.h>
#include <TRUSTTab.h>
#include <communications.h>
#include <Param.h>


Implemente_instanciable_sans_constructeur(Statistiques_dns_ijk_monophasique,"Statistiques_dns_ijk_monophasique",Statistiques_dns_ijk);

#define U_MOY 0
#define V_MOY 1
#define W_MOY 2
#define P_MOY 3
#define UU_MOY 4
#define VV_MOY 5
#define WW_MOY 6
#define UV_MOY 7
#define VW_MOY 8
#define UW_MOY 9
#define UP_MOY 10
#define VP_MOY 11
#define WP_MOY 12
#define UUU_MOY 13
#define VVV_MOY 14
#define WWW_MOY 15
#define UUV_MOY 16
#define WWV_MOY 17
#define U4_MOY 18
#define V4_MOY 19
#define W4_MOY 20
#define dUdx_MOY 21
#define dUdy_MOY 22
#define dUdz_MOY 23
#define dVdx_MOY 24
#define dVdy_MOY 25
#define dVdz_MOY 26
#define dWdx_MOY 27
#define dWdy_MOY 28
#define dWdz_MOY 29
#define EPS_MOY 30
#define dUdx2_MOY 31
#define dUdy2_MOY 32
#define dUdz2_MOY 33
#define dVdx2_MOY 34
#define dVdy2_MOY 35
#define dVdz2_MOY 36
#define dWdx2_MOY 37
#define dWdy2_MOY 38
#define dWdz2_MOY 39

//thermal fields
#define T_MOY 40
#define TU_MOY 41
#define TV_MOY 42
#define TW_MOY 43
#define TUU_MOY 44
#define TVV_MOY 45
#define TWW_MOY 46
#define TUV_MOY 47
#define TUW_MOY 48
#define TVW_MOY 49
#define TdPdx_MOY 50
#define TdPdy_MOY 51
#define TdPdz_MOY 52
#define TdUdx2_MOY 53
#define TdUdy2_MOY 54
#define TdUdz2_MOY 55
#define TdVdx2_MOY 56
#define TdVdy2_MOY 57
#define TdVdz2_MOY 58
#define TdWdx2_MOY 59
#define TdWdy2_MOY 60
#define TdWdz2_MOY 61
#define dTdx2_MOY 62
#define dTdy2_MOY 63
#define dTdz2_MOY 64
#define UdTdx2_MOY 65
#define VdTdy2_MOY 66
#define WdTdz2_MOY 67


Statistiques_dns_ijk_monophasique::Statistiques_dns_ijk_monophasique()
{

  const char *noms_moyennes_prov[] =
  {
    "U",
    "V",
    "W",
    "P",
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
    "UUV",
    "WWV",
    "U4",
    "V4",
    "W4",
    "dUdx",
    "dUdy",
    "dUdz",
    "dVdx",
    "dVdy",
    "dVdz",
    "dWdx",
    "dWdy",
    "dWdz",
    "eps",
    "dUdx2",
    "dUdy2",
    "dUdz2",
    "dVdx2",
    "dVdy2",
    "dVdz2",
    "dWdx2",
    "dWdy2",
    "dWdz2",
//thermal terms
    "T",
    "TU",
    "TV",
    "TW",
    "TUU",
    "TVV",
    "TWW",
    "TUV",
    "TUW",
    "TVW",
    "TdPdx",
    "TdPdy",
    "TdPdz",
    "TdUdx2",
    "TdUdy2",
    "TdUdz2",
    "TdVdx2",
    "TdVdy2",
    "TdVdz2",
    "TdWdx2",
    "TdWdy2",
    "TdWdz2",
    "dTdx2",
    "dTdy2",
    "dTdz2",
    "UdTdx2",
    "VdTdy2",
    "WdTdz2",
  };
  nval_=40;
  nval_+=28; //nb of thermal fields
  noms_moyennes_.dimensionner_force(nval_);
  for (int i=0; i<nval_; i++)
    noms_moyennes_[i]=noms_moyennes_prov[i];
}
void Statistiques_dns_ijk_monophasique::update_stat(const FixedVector<IJK_Field_double, 3>& vitesse,
                                                    const IJK_Field_double& pression,
                                                    const IJK_Field_double& temperature,
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




  // Nombre total de mailles en K
  const int nktot = pression.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);

  // Nombre local de mailles en K
  const int kmax = pression.nk();
  const int offset = pression.get_splitting().get_offset_local(DIRECTION_K);

  DoubleTab tmp(nktot, nval_);
  const int imax = pression.ni();
  const int jmax = pression.nj();

  for (int k = 0; k < kmax; k++)
    {
      const double dz = tab_dz_[k+offset];
      bool on_the_first_cell = false;
      bool on_the_last_cell = false;
      if (k + offset == 0)
        on_the_first_cell = true;
      if (k + offset == nktot)
        on_the_last_cell = true;
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
              double u = (vitesse_i(i,j,k) + vitesse_i(i+1, j, k)) * 0.5;
              double v = (vitesse_j(i,j,k) + vitesse_j(i, j+1, k)) * 0.5;
              double w = (vitesse_k(i,j,k) + vitesse_k(i, j, k+1)) * 0.5;
              double p = pression(i,j,k);
              double T = temperature(i,j,k);

              // ******************************** //
              //	double  duidx, dujdx, dukdx, duidy, dujdy, dukdy, duidz, dujdz, dukdz, dissip;
              double  duidx=0., dujdx=0., dukdx=0., duidy=0., dujdy=0., dukdy=0., duidz=0., dujdz=0., dukdz=0., dissip;
              dissip = face_to_cell_gradient(vitesse_i, vitesse_j, vitesse_k,
                                             i,j,k,
                                             dz,
                                             duidx, dujdx, dukdx,
                                             duidy, dujdy, dukdy,
                                             duidz, dujdz, dukdz,
                                             on_the_first_cell, on_the_last_cell,
                                             0 /* bc_type for velocity : 0 at the wall */);

#define AJOUT(somme,val) moy[somme] += val
              // moyennes
              AJOUT(U_MOY,u);
              AJOUT(V_MOY,v);
              AJOUT(W_MOY,w);
              AJOUT(P_MOY,p);
              // moyennes des carres
              AJOUT(UU_MOY,u*u);
              AJOUT(VV_MOY,v*v);
              AJOUT(WW_MOY,w*w);
              // correlations
              AJOUT(UV_MOY,u*v);
              AJOUT(VW_MOY,v*w);
              AJOUT(UW_MOY,u*w);
              AJOUT(UP_MOY,u*p);
              AJOUT(VP_MOY,v*p);
              AJOUT(WP_MOY,w*p);
              // 3 eme ordre
              AJOUT(UUU_MOY,u*u*u);
              AJOUT(VVV_MOY,v*v*v);
              AJOUT(WWW_MOY,w*w*w);
              AJOUT(UUV_MOY,u*u*v);
              AJOUT(WWV_MOY,w*w*v);
              // 4 eme ordre
              AJOUT(U4_MOY,u*u*u*u);
              AJOUT(V4_MOY,v*v*v*v);
              AJOUT(W4_MOY,w*w*w*w);
              // derivees
              AJOUT(dUdx_MOY, duidx);
              AJOUT(dUdy_MOY, duidy);
              AJOUT(dUdz_MOY, duidz);
              AJOUT(dVdx_MOY, dujdx);
              AJOUT(dVdy_MOY, dujdy);
              AJOUT(dVdz_MOY, dujdz);
              AJOUT(dWdx_MOY, dukdx);
              AJOUT(dWdy_MOY, dukdy);
              AJOUT(dWdz_MOY, dukdz);
              AJOUT(EPS_MOY, dissip);
              AJOUT(dUdx2_MOY, duidx*duidx);
              AJOUT(dUdy2_MOY, duidy*duidy);
              AJOUT(dUdz2_MOY, duidz*duidz);
              AJOUT(dVdx2_MOY, dujdx*dujdx);
              AJOUT(dVdy2_MOY, dujdy*dujdy);
              AJOUT(dVdz2_MOY, dujdz*dujdz);
              AJOUT(dWdx2_MOY, dukdx*dukdx);
              AJOUT(dWdy2_MOY, dukdy*dukdy);
              AJOUT(dWdz2_MOY, dukdz*dukdz);

              //thermique
              AJOUT(T_MOY,T);


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


// Impression des donnees pour reprise
Sortie& Statistiques_dns_ijk_monophasique::printOn(Sortie& os) const
{
  Statistiques_dns_ijk::printOn(os);
  return os;
}

// Reprise des donnees stat dans un fichier reprise
// Attention, cette methode peut etre appelee avant initialize() !
Entree& Statistiques_dns_ijk_monophasique::readOn(Entree& is)
{
  Statistiques_dns_ijk::readOn(is);
  return is;
}
// Attention, cette methode est appelee apres readOn(),
// il ne faut pas casser les donnees lues
void Statistiques_dns_ijk_monophasique::initialize(const IJK_Grid_Geometry& geom)
{
  dx_=geom.get_constant_delta(0); //modif AT 20/06/2013
  dy_=geom.get_constant_delta(1); //modif AT 20/06/2013
  tab_dz_=geom.get_delta(2); //modif GB 30/10/2014
  //  dz_=geom.get_constant_delta(2); //modif GB 07/05/2014
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
      //   moyenne_spatiale_instantanee_[WFACE_MOY].resize_array(n+1);//En effet, il y a une face de plus que d'elem et ce tableau est au face
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
          if (integrale_temporelle_[i].size_array() != n)
            {
              Cerr << "Erreur dans Statistiques_dns_ijk::initialize: reprise avec mauvais nombre de mailles en z" << finl;
              Process::exit();
            }
        }
      else
        {
          integrale_temporelle_[i].resize_array(n);

        }
    }
  //modif AT 20/06/2013
  if (check_converge_)
    {
      //on recupere les vitesses moyennes aux faces pour calculerles fluctuations de vitesse.
      vit_moy_.dimensionner(3);
      for (i = 0; i < 3; i++)
        vit_moy_[i].resize_array(n);

      vit_moy_[0]=integrale_temporelle_[0];
      vit_moy_[1]=integrale_temporelle_[1];
      vit_moy_[2]=integrale_temporelle_[40];//car on veut w aux faces !


    }

}
