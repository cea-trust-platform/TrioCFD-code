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
// File      : Force_ph.cpp
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Force_ph.h>
#include <fstream>
#include <math.h>
#include <IJK_Navier_Stokes_tools.h>
#include <communications.h>
#include <IJK_Splitting.h>

// #include <fftw3.h>

Implemente_instanciable_sans_constructeur_ni_destructeur( Force_ph, "Force_ph", Objet_U ) ;

Sortie& Force_ph::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& Force_ph::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}



Force_ph::Force_ph() : ni(0),nj(0),nk(0),n_ijk(ni*nj*nk),nl(0),nm(0),nn(0),n_lmn(nl*nn*nm),kmin(0),kmax(0),energie(0)
{

}

Force_ph::~Force_ph()
{

}

void Force_ph::initialise(int a_nproc_tot, int a_ni,int a_nj,int a_nk,int a_nl,int a_nm, int a_nn,
                          double a_Lx, double a_Ly, double a_Lz, double a_Ox,double a_Oy,double a_Oz, int a_momin, int a_momax, double a_kmin, double a_kmax,
                          std::string nom_fichier, const IJK_Splitting& splitting,
                          int a_i_offset, int a_j_offset, int a_k_offset
                         )

{
// int ind,i,j,k;
  nproc_tot = a_nproc_tot;

  ni = a_ni;
  nj = a_nj;
  nk = a_nk;
  n_ijk = ni*nj*nk;
  nl = a_nl;
  nm = a_nm;
  nn = a_nn;
  n_lmn = (2*nl+1)*(2*nm+1)*(2*nn+1);
  Lx = a_Lx;
  Ly = a_Ly;
  Lz = a_Lz;
  Ox = a_Ox;
  Oy = a_Oy;
  Oz = a_Oz;
  momin = a_momin;
  momax = a_momax;
  kmin = a_kmin;
  kmax = a_kmax;
  i_offset = a_i_offset;
  j_offset = a_j_offset;
  k_offset = a_k_offset;
// Cout << "a_i_offset : " << a_i_offset << finl;

// NB : Normalement force est juste reel, mais on verifie qd mm ...
  force = set_dimensions(force,2,3,n_ijk);
// TO DO : -->  si on colocalise allocate_cell_velocity au lieu de allocate_velocity  ...
//             voir IJK_FT_post.
  // allocate_velocity(force_, splitting, 2); //: 2 a la fin -> 2 cellules ghost a la fin
  // TENTER CA
  allocate_velocity(force_, splitting, 0); //: 2 a la fin -> 2 cellules ghost a la fin

///////////////////////////////////////////////////////////////
  std::ofstream Physical_flux(nom_fichier.c_str());
  if(Physical_flux)
    {
      Physical_flux << "-------- PHYSICAL_FORCE --------" << std::endl;
      Physical_flux << std::endl << "i,j,k \t : f_x, \tf_y, \tf_z,\t";
    }
///////////////////////////////////////////////////////////////
  std::ofstream Offset_flux("/volatile/FFTW/creeping_flow/creeping_flow_bis/src/offset.txt");
  if(Offset_flux)
    {
      Offset_flux << "-------- PHYSICAL_POSITION --------" << std::endl;
      Offset_flux << std::endl << "i,j,k \t : f_x, \tf_y, \tf_z,\t";
    }

  energie = 0.;

}

void Force_ph::set_zero()
{
  /*Mise a zero de tout le champ de force physique*/
  force_[0].data() = 0.;
  force_[1].data() = 0.;
  force_[2].data() = 0.;
}

void Force_ph::cheat_function()
{
  /*
  Fonction de test. Ecrit directement la fonction cos(x-z) ex dans le champ
  de force physique sans passer par le champ spectral.
  QUESTION : Si je mets += au lieu de = a partir du premier pdt
             j'obtiens une amplitude -6;+6 avec le time_scheme : RK3
  */
  double pos[3];
  int roc(momin+0);
  double t(2.*M_PI*roc/Lx);

  for (int k=0; k<nk; k++)
    {
      pos[2] = Oz + k*Lz/nk;
      for (int j=0; j<nj; j++)
        {
          pos[1] = Oy + j*Ly/nj;
          for (int i=0; i<ni; i++)
            {
              int ind_ijk = (k*nj + j)*ni + i;
              pos[0] = Ox + i*Lx/ni;

              // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              force[0][0][ind_ijk] = 2*cos(t*(pos[0]-pos[2]));
              force[0][2][ind_ijk] = 2*cos(t*(pos[0]-pos[2]));
              force[1][0][ind_ijk] = 0.;
              force[1][2][ind_ijk] = 0.;
              // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              force_[0](i,j,k) = 2*cos(t*(pos[0]-pos[2]));
              force_[2](i,j,k) = 2*cos(t*(pos[0]-pos[2]));
              // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            }
        }
    }
}


void Force_ph::from_spect_to_phys2(const std:: vector <double >& coeff_force)
{
  //////////////////////////////////////////////////////////////////////////////
  /*
  Version primitive de la transformee de Fourier inverse :
       - 6 boucles imbriquees (3 d'espace physique, 3 d'espace spectral).
       - coefficients de la force spectrale stockes dans un tableau
         a une dimension, on reconstruit l'indice a aller chercher.
       - coefficients de la force physique recopies dans un tableau
         a trois dimensions.
  */
  //////////////////////////////////////////////////////////////////////////////

  // Initialisation du champ de force
  for (int dir = 0; dir < 3; dir++)
    {
      force_[dir].data() = 0.;
    }

  const int n_ind_mni((2*nm+1)*(2*nn+1)*1*ni);
  const int n_ind_nij((2*nn+1)*1*ni*1*nj);
  const int n_ind_lmn((2*nl+1)*(2*nm+1)*(2*nn+1));

  std:: vector <double > bout_a_quart(2*3*n_ind_mni,0.);
  std:: vector <double > bout_b_quart(2*3*n_ind_nij,0.);

  double expo_a[2], expo_b[2], expo_c[2];
  double pos[3];
  double dz(Lz/nk), dy(Ly/nj), dx(Lx/ni);

  int n_dir = 3;

  //////////////////////////////////////////////////////////////
  // TODO : remplacer mes pos[X] par un resultat de build_local_coords
  //          --> Dans IJK_Switch ou IJK_switch_FT
  for (int k=0; k<nk; k++)
    {
      pos[2] = Oz + k*Lz/nk + dz/2.;
      for (int j=0; j<nj; j++)
        {
          pos[1] = Oy + j*Ly/nj + dy/2.;
          for (int i=0; i<ni; i++)
            {
              pos[0] = Ox + i*Lx/ni + dx/2.;
              for (int n=0; n<2*nn; n++)
                {
                  int ind_nij = ( (j*ni + i) * (2*nn+1) + n );

                  // Initialisation de la valeur de bout_b

                  for(int dir=0; dir<3; dir++)
                    {
                      int ind_RDInij((0*n_dir+dir) * n_ind_nij + ind_nij);
                      int ind_CDInij((1*n_dir+dir) * n_ind_nij + ind_nij);
                      bout_b_quart[ind_RDInij] = 0.0;
                      bout_b_quart[ind_CDInij] = 0.0;
                    }

                  double kappa_n = (- kmax + (n)*(2*kmax)/(2*nn) );
                  for (int m=0; m<2*nm; m++)
                    {
                      int ind_mni = ((i*(2*nn+1) + n) * (2*nm+1) + m );

                      // Initialisation de la valeur de bout_a

                      for(int dir=0; dir<3; dir++)
                        {
                          int ind_RDImni((0*n_dir+dir) * n_ind_mni + ind_mni);
                          int ind_CDImni((1*n_dir+dir) * n_ind_mni + ind_mni);
                          bout_a_quart[ind_RDImni] = 0.0;
                          bout_a_quart[ind_CDImni] = 0.0;
                        }

                      double kappa_m = (- kmax + (m)*(2*kmax)/(2*nm));

                      for (int l=0; l<2*nl; l++)
                        {
                          int ind_lmn = ( (n*(2*nm+1) + m) * (2*nl+1) +l);
                          double kappa_l = (- kmax + (l)*(2*kmax)/(2*nl));

                          expo_a[0] = cos(kappa_l*pos[0]);
                          expo_a[1] = sin(kappa_l*pos[0]);

                          for (int dir=0; dir<3; dir++)
                            {
                              int ind_RDImni((0*n_dir+dir)*n_ind_mni+ind_mni);
                              int ind_RDIlmn((0*n_dir+dir)*n_ind_lmn+ind_lmn);
                              int ind_CDImni((1*n_dir+dir)*n_ind_mni+ind_mni);
                              int ind_CDIlmn((1*n_dir+dir)*n_ind_lmn+ind_lmn);

                              bout_a_quart[ind_RDImni] += 1*(coeff_force[ind_RDIlmn]*expo_a[0] - coeff_force[ind_CDIlmn]*expo_a[1]);
                              bout_a_quart[ind_CDImni] += 1*(coeff_force[ind_RDIlmn]*expo_a[1] + coeff_force[ind_CDIlmn]*expo_a[0]);
                            }
                        }
                      expo_b[0] = cos(kappa_m*pos[1]);
                      expo_b[1] = sin(kappa_m*pos[1]);
                      for (int dir=0; dir<3; dir++)
                        {
                          int ind_RDInij((0*n_dir+dir)*n_ind_nij+ind_nij);
                          int ind_RDImni((0*n_dir+dir)*n_ind_mni+ind_mni);
                          int ind_CDInij((1*n_dir+dir)*n_ind_nij+ind_nij);
                          int ind_CDImni((1*n_dir+dir)*n_ind_mni+ind_mni);

                          bout_b_quart[ind_RDInij] += 1*(bout_a_quart[ind_RDImni]*expo_b[0] - bout_a_quart[ind_CDImni]*expo_b[1]);
                          bout_b_quart[ind_CDInij] += 1*(bout_a_quart[ind_RDImni]*expo_b[1] + bout_a_quart[ind_CDImni]*expo_b[0]);
                        }
                    }
                  expo_c[0] = cos(kappa_n*pos[2]);
                  expo_c[1] = sin(kappa_n*pos[2]);
                  for (int dir=0; dir<3; dir++)
                    {
                      // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      // ~ RDI : Real,Direction,Index
                      // ~ CDI : Complex,Direction,Index
                      int ind_RDInij((0*n_dir+dir)*n_ind_nij+ind_nij);
                      int ind_CDInij((1*n_dir+dir)*n_ind_nij+ind_nij);

                      // force[0][dir][ind_ijk] += 1*(bout_b_quart[ind_RDInij]*expo_c[0] - bout_b_quart[ind_CDInij]*expo_c[1]);
                      // force[1][dir][ind_ijk] += 1*(bout_b_quart[ind_RDInij]*expo_c[1] + bout_b_quart[ind_CDInij]*expo_c[0]);

                      force_[dir](i,j,k) += 1*(bout_b_quart[ind_RDInij]*expo_c[0] - bout_b_quart[ind_CDInij]*expo_c[1]);
                    }
                }
            }
        }
    }
}



void Force_ph::from_spect_to_phys_opti2(ArrOfDouble& coeff_force )
{
  //////////////////////////////////////////////////////////////////////////////
  /*
  Version la plus avancee de la transformee de Fourier inverse :
       - Les 6 boucles imbriquees sont eclatees en 3 blocs de complexite
         au plus ni*nj*nk*nn.
       - coefficients de la force spectrale stockes dans un tableau
         a une dimension, on reconstruit l'indice a aller chercher.
       - coefficients de la force physique recopies dans un tableau
         a trois dimensions (a changer, facile a changer), mais aussi stockes dans un
         FixedVector<IJK_Field_double, 3>
  */
  //////////////////////////////////////////////////////////////////////////////

  const int ni_o(ni);
  const int nj_o(nj);

  // Initialisation du champ de force
  for (int dir = 0; dir < 3; dir++)
    {
      force_[dir].data() = 0.;
    }
  int n_dir(3);
  int n_ind_mni((2*nm+1)*(2*nn+1)*1*ni_o);
  int n_ind_nij((2*nn+1)*1*ni_o*1*nj_o);
  int n_ind_lmn((2*nl+1)*(2*nm+1)*(2*nn+1));

  std:: vector <double > bout_a_quart(2*3*n_ind_mni,0.);
  std:: vector <double > bout_b_quart(2*3*n_ind_nij,0.);

  double expo_a[2], expo_b[2], expo_c[2];
  double pos[3];

  //////////////////////////////////////////////////////////////////////////////
  // Construction des posisitons. Comme la force est localisee (pour le moment)
  // aux faces, on est obliges de definir les coord pour les faces x, celles pour les faces y
  // et puis pour les faces z.
  // 1) Reconstruire a la main ferai economiser 6 ArrOfDouble
  // 2) localisr aux centres (sans reconstruire a la main) ferai econmiser 6 ArrOfDouble
  ArrOfDouble pos_x_fx, pos_y_fx, pos_z_fx;
  build_local_coords(force_[0], pos_x_fx, pos_y_fx, pos_z_fx);
  ArrOfDouble pos_x_fy, pos_y_fy, pos_z_fy;
  build_local_coords(force_[1], pos_x_fy, pos_y_fy, pos_z_fy);
  ArrOfDouble pos_x_fz, pos_y_fz, pos_z_fz;
  build_local_coords(force_[2], pos_x_fz, pos_y_fz, pos_z_fz);
  //////////////////////////////////////////////////////////////////////////////

  for (int dir=0; dir<3; dir++)
    {
      ///////////////////////////////////////////////////////////////
      for (int m=0; m<2*nm+1; m++)
        for (int n=0; n<2*nn+1; n++)
          for (int i=0; i<ni; i++)
            {
              int ind_mni = ((i*(2*nn+1) + n) * (2*nm+1) + m );
              // pos[0] = Ox + i*Lx/(ni) - dx/2.;  // Localise aux faces
              // pos[0] = Ox + i*Lx/(ni+i_offset) - dx/2.;  // Localise aux faces
              pos[0] = pos_x_fx[i]; // Ox + i*Lx/(ni+i_offset) - dx/2.;  // Localise aux faces
              for (int l=0; l<2*nl+1; l++)
                {
                  int ind_lmn = ( (n*(2*nm+1) + m) * (2*nl+1) +l);
                  double kappa_l = (- kmax + (l)*(2*kmax)/(2*nl));
                  expo_a[0] = cos(kappa_l*pos[0]);
                  expo_a[1] = sin(kappa_l*pos[0]);

                  int ind_RDImni((0*n_dir+dir) * n_ind_mni + ind_mni);
                  int ind_RDIlmn((0*n_dir+dir) * n_ind_lmn + ind_lmn);
                  int ind_CDImni((1*n_dir+dir) * n_ind_mni + ind_mni);
                  int ind_CDIlmn((1*n_dir+dir) * n_ind_lmn + ind_lmn);

                  bout_a_quart[ind_RDImni] += 1*(coeff_force[ind_RDIlmn]*expo_a[0] - coeff_force[ind_CDIlmn]*expo_a[1]);
                  bout_a_quart[ind_CDImni] += 1*(coeff_force[ind_RDIlmn]*expo_a[1] + coeff_force[ind_CDIlmn]*expo_a[0]);
                }
            }
      ///////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////
      for (int n=0; n<2*nn+1; n++)
        for (int i=0; i<ni; i++)
          for (int j=0; j<nj; j++)
            {
              int ind_nij = ( (j*ni_o + i) * (2*nn+1) + n );
              // pos[1] = Oy + j*Ly/(nj) - dy/2.; // Localise aux faces
              // pos[1] = Oy + j*Ly/(nj+j_offset) - dy/2.; // Localise aux faces
              pos[1] = pos_y_fy[j]; // Localise aux faces
              for (int m=0; m<2*nm; m++)
                {
                  int ind_mni = ((i*(2*nn+1) + n) * (2*nm+1) + m );
                  double kappa_m = (- kmax + (m)*(2*kmax)/(2*nm));
                  expo_b[0] = cos(kappa_m*pos[1]);
                  expo_b[1] = sin(kappa_m*pos[1]);

                  int ind_RDInij((0*n_dir+dir)*n_ind_nij+ind_nij);
                  int ind_RDImni((0*n_dir+dir)*n_ind_mni+ind_mni);
                  int ind_CDInij((1*n_dir+dir)*n_ind_nij+ind_nij);
                  int ind_CDImni((1*n_dir+dir)*n_ind_mni+ind_mni);

                  bout_b_quart[ind_RDInij] += 1*(bout_a_quart[ind_RDImni]*expo_b[0] - bout_a_quart[ind_CDImni]*expo_b[1]);
                  bout_b_quart[ind_CDInij] += 1*(bout_a_quart[ind_RDImni]*expo_b[1] + bout_a_quart[ind_CDImni]*expo_b[0]);
                }
            }
      ///////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////
      // GAB, question : pour cette boucle il vaut mieux mettre la boucle en k le plus a lexterieur, on est d'acc ?
      for (int k=0; k<nk; k++)
        for (int j=0; j<nj; j++)
          for (int i=0; i<ni; i++)
            {
              // pos[2] = Oz + k*Lz/(nk) - dz/2.;   // Localise aux faces
              // pos[2] = Oz + k*Lz/(nk+k_offset) - dz/2.;   // Localise aux faces
              pos[2] = pos_z_fz[k];   // Localise aux faces
              for (int n=0; n<2*nn+1; n++)
                {
                  int ind_nij = ( (j*ni_o + i) * (2*nn+1) + n );
                  double kappa_n = - kmax + (n)*(2*kmax)/(2*nn);
                  expo_c[0] = cos(kappa_n*pos[2]);
                  expo_c[1] = sin(kappa_n*pos[2]);

                  int ind_RDInij((0*n_dir+dir)*n_ind_nij+ind_nij);
                  int ind_CDInij((1*n_dir+dir)*n_ind_nij+ind_nij);
                  // Smetrie Hermitienne :
                  // force[0][dir][ind_ijk] += 1*(bout_b_quart[ind_RDInij]*expo_c[0] - bout_b_quart[ind_CDInij]*expo_c[1]);
                  // force[1][dir][ind_ijk] += 1*(bout_b_quart[ind_RDInij]*expo_c[1] + bout_b_quart[ind_CDInij]*expo_c[0]);

                  force_[dir](i,j,k) += 1*(bout_b_quart[ind_RDInij]*expo_c[0] - bout_b_quart[ind_CDInij]*expo_c[1]);
                }
            }
      ///////////////////////////////////////////////////////////////
      // Si on plante ici, mettre un 0 au moment du allocate
      force_[dir].echange_espace_virtuel(force_[dir].ghost());
    }
}

void Force_ph::from_spect_to_phys_opti2_advection(ArrOfDouble& coeff_force, const ArrOfDouble& advection_length )
{
  //////////////////////////////////////////////////////////////////////////////
  /*
  Version la plus avancee de la transformee de Fourier inverse, avec advection du champ de force :
       - Les 6 boucles imbriquees sont eclatees en 3 blocs de complexite
         au plus ni*nj*nk*nn.
       - coefficients de la force spectrale stockes dans un tableau
         a une dimension, on reconstruit l'indice a aller chercher.
       - coefficients de la force physique recopies dans un tableau
         a trois dimensions (a changer, facile a changer), mais aussi stockes dans un
         FixedVector<IJK_Field_double, 3>
  */
  //////////////////////////////////////////////////////////////////////////////

  const int ni_o(ni);
  const int nj_o(nj);

  // Initialisation du champ de force
  for (int dir = 0; dir < 3; dir++)
    {
      force_[dir].data() = 0.;
    }
  int n_dir(3);
  int n_ind_mni((2*nm+1)*(2*nn+1)*1*ni_o);
  int n_ind_nij((2*nn+1)*1*ni_o*1*nj_o);
  int n_ind_lmn((2*nl+1)*(2*nm+1)*(2*nn+1));

  std:: vector <double > bout_a_quart(2*3*n_ind_mni,0.);
  std:: vector <double > bout_b_quart(2*3*n_ind_nij,0.);

  double expo_a[2], expo_b[2], expo_c[2];
  double pos[3];

  //////////////////////////////////////////////////////////////////////////////
  // Construction des posisitons. Comme la force est localisee (pour le moment)
  // aux faces, on est obliges de definir les coord pour les faces x, celles pour les faces y
  // et puis pour les faces z.
  // 1) Reconstruire a la main ferai economiser 6 ArrOfDouble
  // 2) localisr aux centres (sans reconstruire a la main) ferai econmiser 6 ArrOfDouble
  ArrOfDouble pos_x_fx, pos_y_fx, pos_z_fx;
  build_local_coords(force_[0], pos_x_fx, pos_y_fx, pos_z_fx);
  ArrOfDouble pos_x_fy, pos_y_fy, pos_z_fy;
  build_local_coords(force_[1], pos_x_fy, pos_y_fy, pos_z_fy);
  ArrOfDouble pos_x_fz, pos_y_fz, pos_z_fz;
  build_local_coords(force_[2], pos_x_fz, pos_y_fz, pos_z_fz);
  //////////////////////////////////////////////////////////////////////////////

  for (int dir=0; dir<3; dir++)
    {
      ///////////////////////////////////////////////////////////////
      for (int m=0; m<2*nm+1; m++)
        for (int n=0; n<2*nn+1; n++)
          for (int i=0; i<ni; i++)
            {
              int ind_mni = ((i*(2*nn+1) + n) * (2*nm+1) + m );
              // pos[0] = Ox + i*Lx/(ni) - dx/2.;  // Localise aux faces
              // pos[0] = Ox + i*Lx/(ni+i_offset) - dx/2.;  // Localise aux faces
              pos[0] = pos_x_fx[i]; // Ox + i*Lx/(ni+i_offset) - dx/2.;  // Localise aux faces
              for (int l=0; l<2*nl+1; l++)
                {
                  int ind_lmn = ( (n*(2*nm+1) + m) * (2*nl+1) +l);
                  double kappa_l = (- kmax + (l)*(2*kmax)/(2*nl));
                  expo_a[0] = cos(kappa_l*(pos[0]-advection_length[0]));
                  expo_a[1] = sin(kappa_l*(pos[0]-advection_length[0]));

                  int ind_RDImni((0*n_dir+dir) * n_ind_mni + ind_mni);
                  int ind_RDIlmn((0*n_dir+dir) * n_ind_lmn + ind_lmn);
                  int ind_CDImni((1*n_dir+dir) * n_ind_mni + ind_mni);
                  int ind_CDIlmn((1*n_dir+dir) * n_ind_lmn + ind_lmn);

                  bout_a_quart[ind_RDImni] += 1*(coeff_force[ind_RDIlmn]*expo_a[0] - coeff_force[ind_CDIlmn]*expo_a[1]);
                  bout_a_quart[ind_CDImni] += 1*(coeff_force[ind_RDIlmn]*expo_a[1] + coeff_force[ind_CDIlmn]*expo_a[0]);
                }
            }
      ///////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////
      for (int n=0; n<2*nn+1; n++)
        for (int i=0; i<ni; i++)
          for (int j=0; j<nj; j++)
            {
              int ind_nij = ( (j*ni_o + i) * (2*nn+1) + n );
              // pos[1] = Oy + j*Ly/(nj) - dy/2.; // Localise aux faces
              // pos[1] = Oy + j*Ly/(nj+j_offset) - dy/2.; // Localise aux faces
              pos[1] = pos_y_fy[j]; // Localise aux faces
              for (int m=0; m<2*nm; m++)
                {
                  int ind_mni = ((i*(2*nn+1) + n) * (2*nm+1) + m );
                  double kappa_m = (- kmax + (m)*(2*kmax)/(2*nm));
                  expo_b[0] = cos(kappa_m*(pos[1]-advection_length[1]));
                  expo_b[1] = sin(kappa_m*(pos[1]-advection_length[1]));

                  int ind_RDInij((0*n_dir+dir)*n_ind_nij+ind_nij);
                  int ind_RDImni((0*n_dir+dir)*n_ind_mni+ind_mni);
                  int ind_CDInij((1*n_dir+dir)*n_ind_nij+ind_nij);
                  int ind_CDImni((1*n_dir+dir)*n_ind_mni+ind_mni);

                  bout_b_quart[ind_RDInij] += 1*(bout_a_quart[ind_RDImni]*expo_b[0] - bout_a_quart[ind_CDImni]*expo_b[1]);
                  bout_b_quart[ind_CDInij] += 1*(bout_a_quart[ind_RDImni]*expo_b[1] + bout_a_quart[ind_CDImni]*expo_b[0]);
                }
            }
      ///////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////
      // GAB, question : la boucle en k est exterieure. C'est le plus optimal ? ?
      for (int k=0; k<nk; k++)
        for (int j=0; j<nj; j++)
          for (int i=0; i<ni; i++)
            {
              // pos[2] = Oz + k*Lz/(nk) - dz/2.;   // Localise aux faces
              // pos[2] = Oz + k*Lz/(nk+k_offset) - dz/2.;   // Localise aux faces
              pos[2] = pos_z_fz[k];   // Localise aux faces
              for (int n=0; n<2*nn+1; n++)
                {
                  int ind_nij = ( (j*ni_o + i) * (2*nn+1) + n );
                  double kappa_n = - kmax + (n)*(2*kmax)/(2*nn);
                  expo_c[0] = cos(kappa_n*(pos[2]-advection_length[2]));
                  expo_c[1] = sin(kappa_n*(pos[2]-advection_length[2]));

                  int ind_RDInij((0*n_dir+dir)*n_ind_nij+ind_nij);
                  int ind_CDInij((1*n_dir+dir)*n_ind_nij+ind_nij);
                  // Symetrie Hermitienne :
                  // force[0][dir][ind_ijk] += 1*(bout_b_quart[ind_RDInij]*expo_c[0] - bout_b_quart[ind_CDInij]*expo_c[1]);
                  // force[1][dir][ind_ijk] += 1*(bout_b_quart[ind_RDInij]*expo_c[1] + bout_b_quart[ind_CDInij]*expo_c[0]);

                  force_[dir](i,j,k) += 1*(bout_b_quart[ind_RDInij]*expo_c[0] - bout_b_quart[ind_CDInij]*expo_c[1]);
                }
            }
      ///////////////////////////////////////////////////////////////
    }

}




void Force_ph::write(std::string nom_fichier_sortie, double t)
{
  //////////////////////////////////////////////////////////////////////////////
  // ATTENTION : Cette fonction ecrit force et pas force_
  //             On a supprime force maintenant. Il faut re ecrire cette fct
  //             si on souhaite l'utiliser
  //////////////////////////////////////////////////////////////////////////////
  std::cout << "test in force_ph/write" << std::endl;
  std::ofstream Physical_flux(nom_fichier_sortie.c_str(), std::ios::app);
  if (Physical_flux)
    {
      int i, j, k, ind;

      Physical_flux << std::endl << "time : " << t << std::endl << std::endl;
      for (k=0; k<nk; k++)
        for (j=0; j<nj; j++)
          for (i=0; i<ni; i++)
            {
              ind = (k*nj + j)*ni + i;
              Physical_flux << i<<","<<j<<","<<k<<"\t : ";
              Physical_flux << force[0][0][ind] << " + i" << force[1][0][ind]<<", \t";
              Physical_flux << force[0][1][ind] << " + i" << force[1][1][ind]<<", \t";
              Physical_flux << force[0][2][ind] << " + i" << force[1][2][ind]<<", \t";
              Physical_flux << std::endl;
            }
    }
}


void Force_ph::write_offset_index_position( const IJK_Splitting& my_splitting)
{
  std::ofstream Offset_flux("/volatile/FFTW/creeping_flow/creeping_flow_bis/src/offset.txt", std::ios::app);
  if (Offset_flux)
    {

      Offset_flux << std::endl << "of_i,of_j,of_k,i,j,k,x,y,z"<< std::endl;
      double dz(Lz/nk), dy(Ly/nj), dx(Lx/ni);
      for (int k=0; k<nk; k++)
        {
          double z = Oz + k*Lz/nk + dz/2.;
          for (int j=0; j<nj; j++)
            {
              double y = Oy + j*Ly/nj + dy/2.;
              for (int i=0; i<ni; i++)
                {
                  double x = Ox + i*Lx/ni + dx/2.;
                  Offset_flux << my_splitting.get_offset_local(DIRECTION_I)<<",";
                  Offset_flux << my_splitting.get_offset_local(DIRECTION_J)<<",";
                  Offset_flux << my_splitting.get_offset_local(DIRECTION_K)<<";";
                  Offset_flux << i<<","<<j<<","<<k<<";";
                  Offset_flux << x<<","<<y<<","<<z<<";";
                  Offset_flux << std::endl;
                }
            }
        }
    }
}




void Force_ph::compute_energie()
{
  int i,j,k,dir,ind;
  for (k=0; k<nk; k++)
    for (j=0; j<nj; j++)
      for (i=0; i<ni; i++)
        {
          ind = (k*nj + j)*ni + i;
          for (dir=0; dir<3; dir++)
            energie += force[0][dir][ind]*force[0][dir][ind] + force[1][dir][ind]*force[1][dir][ind];
        }
  energie *= ((1.)/(ni*nj*nk)) ;
}

double Force_ph::get_energie()
{
  return energie;
}

void Force_ph::write_separate(std::string nom_fichier_sortie, double t)
{
  //////////////////////////////////////////////////////////////////////////////
  // ATTENTION : Cette fonction ecrit force et pas force_
  //             On a supprime force maintenant. Il faut re ecrire cette fct
  //             si on souhaite l'utiliser
  //////////////////////////////////////////////////////////////////////////////
  std::ofstream Physical_flux(nom_fichier_sortie.c_str());
  if (Physical_flux)
    {
      int i, j, k, ind;
      double x,y,z;
      Physical_flux << std::endl << "i,j,k, \t x,y,z, \t  rf_x, \t cf_x, \t\t rf_y, \t cf_y, \t\t rf_z, \t cf_z \t";
      Physical_flux << std::endl;

      for (k=0; k<nk; k++)
        for (j=0; j<nj; j++)
          for (i=0; i<ni; i++)
            {
              // ind = (k*nk + j)*nj +i;
              ind = (k*nj + j)*ni + i;

              x = -Lx/2 + i*Lx/ni;
              y = -Ly/2 + j*Ly/nj;
              z = -Lz/2 + k*Lz/nk;

              Physical_flux << i<<","<<j<<","<<k<<"\t,";
              Physical_flux << x << ",\t" << y << ",\t" << z <<", \t\t";
              Physical_flux << force[0][0][ind] << ",\t" << force[1][0][ind]<<", \t\t";
              Physical_flux << force[0][1][ind] << ",\t" << force[1][1][ind]<<", \t\t";
              Physical_flux << force[0][2][ind] << ",\t" << force[1][2][ind]<<"\t";
              Physical_flux << std::endl;
            }
    }
}

FixedVector<IJK_Field_double, 3> Force_ph::get_force_attribute()
{
  return force_;
}

FixedVector<IJK_Field_double, 3>& Force_ph::get_force_attribute2()
{
  return force_;
}


std::vector< std::vector< std:: vector <double >>> set_dimensions(std::vector< std::vector< std:: vector <double >>> the_vector,
                                                                  int dim_one, int dim_two, int dim_three)
{
  int i,j;
  the_vector.resize(dim_one);
  for (i=0; i<dim_one; i++)
    {
      the_vector[i].resize(dim_two);
      for (j=0; j<dim_two; j++)
        {
          the_vector[i][j].resize(dim_three,0.0);
        }
    }
  return the_vector;
}
