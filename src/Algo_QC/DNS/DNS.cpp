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
// File      : DNS_QC_double.cpp
// Directory : $NEW_ALGO_QC_ROOT/src/DNS_QC_double
//
/////////////////////////////////////////////////////////////////////////////
#include <DNS.h>
#include <Param.h>
#include <Interprete_bloc.h>
#include <SFichier.h>
#include <algorithm>
#include <alloca.h>
#include <stat_counters.h>
#include <IJK_Lata_writer.h>
#include <IJK_Navier_Stokes_tools.h>
#include <communications.h>
#include <LecFicDiffuse_JDD.h>
#include <SFichier.h>
#include <EFichier.h>
#include <Boundary_Conditions.h>
#include <DebogIJK.h>
#include <Statistiques.h>
#include <stat_counters.h>
#include <Turbulent_viscosity.h>
#include <Filter_kernel.h>
using std::min;
using std::max;
// Reste a faire:
//  ajouter terme u_div_rho_u dans convection vitesse
//  verifier quick en paroi (=> passer a l'ordre 1)

Implemente_instanciable(DNS_QC_double, "DNS_QC_double", Interprete);

Sortie& DNS_QC_double::printOn(Sortie& s) const
{
  return s;
}
Entree& DNS_QC_double::readOn(Entree& s)
{
  return s;
}


// Programme principal pour le calcul DNS QC sur SuperMUC

// formule pour gnuplot
// lambda(x)=(((((x*(-5.05628e-18)+2.469e-14 ) *x-4.98344e-11)*x)+7.06714e-08)*x+1.0894e-06)*1.93198026315789000e-3 / 1.461e-6
// mu(x)=(((((x*(-5.05628e-18)+2.469e-14 ) *x-4.98344e-11)*x)+7.06714e-08)*x+1.0894e-06)
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
}

static inline void choix_filter_kernel(const int ghost_size,
                                       const Nom& filter_kernel_name,
                                       Filter_kernel_base*& filter_kernel)
{
  // Choose a filter_kernel to compute the filter.
  // Note: THE CALLER _MUdouble_ FREE THE MEMORY USING DELETE.
  // example of use :
  //          Filter_kernel_base* filter_kernel = NULL;
  //          choix_filter_kernel(ghost_size, filter_kernel_name, filter_kernel);
  //          /* do something with filter_kernel */
  //          delete filter_kernel;

  if ( filter_kernel_name == Nom("box") )
    {
      filter_kernel = new Filter_kernel_box(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_13_13_pondere") )
    {
      filter_kernel = new Filter_kernel_weight_13_13_pondere(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_13_13_sansponderation") )
    {
      filter_kernel = new Filter_kernel_weight_13_13_sansponderation(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_13_13_conservatif") )
    {
      filter_kernel = new Filter_kernel_weight_13_13_conservatif(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_12_14_pondere") )
    {
      filter_kernel = new Filter_kernel_weight_12_14_pondere(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_12_14_sansponderation") )
    {
      filter_kernel = new Filter_kernel_weight_12_14_sansponderation(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_12_14_conservatif") )
    {
      filter_kernel = new Filter_kernel_weight_12_14_conservatif(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_23_16_pondere") )
    {
      filter_kernel = new Filter_kernel_weight_23_16_pondere(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_23_16_sansponderation") )
    {
      filter_kernel = new Filter_kernel_weight_23_16_sansponderation(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_23_16_conservatif") )
    {
      filter_kernel = new Filter_kernel_weight_23_16_conservatif(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_14_14_18_pondere") )
    {
      filter_kernel = new Filter_kernel_weight_14_14_18_pondere(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_14_14_18_sansponderation") )
    {
      filter_kernel = new Filter_kernel_weight_14_14_18_sansponderation(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_14_14_18_conservatif") )
    {
      filter_kernel = new Filter_kernel_weight_14_14_18_conservatif(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_15_15_15_pondere") )
    {
      filter_kernel = new Filter_kernel_weight_15_15_15_pondere(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_15_15_15_sansponderation") )
    {
      filter_kernel = new Filter_kernel_weight_15_15_15_sansponderation(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_15_15_15_conservatif") )
    {
      filter_kernel = new Filter_kernel_weight_15_15_15_conservatif(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_16_16_16_112_pondere") )
    {
      filter_kernel = new Filter_kernel_weight_16_16_16_112_pondere(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_16_16_16_112_sansponderation") )
    {
      filter_kernel = new Filter_kernel_weight_16_16_16_112_sansponderation(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_16_16_16_112_conservatif") )
    {
      filter_kernel = new Filter_kernel_weight_16_16_16_112_conservatif(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_14_38_pondere") )
    {
      filter_kernel = new Filter_kernel_weight_14_38_pondere(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_14_38_sansponderation") )
    {
      filter_kernel = new Filter_kernel_weight_14_38_sansponderation(ghost_size);
    }
  else if ( filter_kernel_name == Nom("weight_14_38_conservatif") )
    {
      filter_kernel = new Filter_kernel_weight_14_38_conservatif(ghost_size);
    }
  else if ( filter_kernel_name == Nom("laplacian") )
    {
      filter_kernel = new Filter_kernel_laplacian(ghost_size);
    }
  else
    {
      Cerr << "Error: (Inconsistent parameters) "
           << "The large eddy simulation model requires filtering but the filter kernel name is unknown or unspecified. "
           << "To specify a filter kernel, you may use the keyword FILTER_KERNEL_NAME with the parameters BOX, GAUSSIAN, GAUSSIAN2, LAPLACIAN_SIMPLE, LAPLACIAN or NONE." << finl;
      Process::exit();
    }
}

static inline void choix_modele(const Nom& turbulent_viscosity_model,
                                Turbulent_viscosity_base*& model)
{
  // Choose a model to compute the turbulent viscosity.
  // Note: THE CALLER _MUdouble_ FREE THE MEMORY USING DELETE.
  // example of use :
  //          Turbulent_viscosity_base* model = NULL;
  //          choix_modele(turbulent_viscosity_model, model);
  //          /* do something with model */
  //          delete model;

  if ( turbulent_viscosity_model == Nom("constant") )
    {
      model = new Turbulent_viscosity_constant;
    }
  else if ( turbulent_viscosity_model == Nom("unsrho") )
    {
      model = new Turbulent_viscosity_unsrho;
    }
  else if ( turbulent_viscosity_model == Nom("smagorinsky") )
    {
      model = new Turbulent_viscosity_smagorinsky;
    }
  else if ( turbulent_viscosity_model == Nom("vreman") )
    {
      model = new Turbulent_viscosity_vreman;
    }
  else if ( turbulent_viscosity_model == Nom("sigma") )
    {
      model = new Turbulent_viscosity_sigma;
    }
  else if ( turbulent_viscosity_model == Nom("wale") )
    {
      model = new Turbulent_viscosity_wale;
    }
  else if ( turbulent_viscosity_model == Nom("amd") )
    {
      model = new Turbulent_viscosity_amd;
    }
  else if ( turbulent_viscosity_model == Nom("amd_comp") )
    {
      model = new Turbulent_viscosity_amd_comp;
    }
  else if ( turbulent_viscosity_model == Nom("amdnoclip") )
    {
      model = new Turbulent_viscosity_amdnoclip;
    }
  else if ( turbulent_viscosity_model == Nom("amdscalar") )
    {
      model = new Turbulent_viscosity_amdscalar;
    }
  else if ( turbulent_viscosity_model == Nom("amdscalarnoclip") )
    {
      model = new Turbulent_viscosity_amdscalarnoclip;
    }
  else if ( turbulent_viscosity_model == Nom("rds") )
    {
      model = new Turbulent_viscosity_rds;
    }
  else if ( turbulent_viscosity_model == Nom("vss") )
    {
      model = new Turbulent_viscosity_vss;
    }
  else if ( turbulent_viscosity_model == Nom("kobayashi") )
    {
      model = new Turbulent_viscosity_kobayashi;
    }
  else
    {
      Cerr << "The turbulent diffusion model name is unknown." << finl;
      Process::exit();
    }
}

void calculer_delta_z_pour_delta(const IJK_Splitting& splitting,
                                 const IJK_Grid_Geometry& geom_pour_delta,
                                 const int ghost_size,
                                 ArrOfDouble_with_ghost& delta_z_pour_delta)
{
  const int nktot = splitting.get_grid_geometry().get_nb_elem_tot(DIRECTION_K);
  const ArrOfDouble& coord_z = splitting.get_grid_geometry().get_node_coordinates(DIRECTION_K);
  ArrOfDouble elem_coord(nktot);
  for (int i = 0; i < nktot; i++)
    {
      elem_coord[i] = 0.5 * (coord_z[i] + coord_z[i+1]);
    }

  const int nktot_pour_delta = geom_pour_delta.get_nb_elem_tot(DIRECTION_K);
  const ArrOfDouble& global_delta_pour_delta = geom_pour_delta.get_delta(DIRECTION_K);
  const ArrOfDouble& coord_z_pour_delta = geom_pour_delta.get_node_coordinates(DIRECTION_K);
  ArrOfDouble elem_coord_pour_delta(nktot_pour_delta);
  for (int i = 0; i < nktot_pour_delta; i++)
    {
      elem_coord_pour_delta[i] = 0.5 * (coord_z_pour_delta[i] + coord_z_pour_delta[i+1]);
    }

  const int offset = splitting.get_offset_local(DIRECTION_K);
  const int nk = splitting.get_nb_elem_local(DIRECTION_K);
  delta_z_pour_delta.resize(nk, ghost_size);

  double d0 = 0.;
  double d1 = 0.;
  double coord_0 = 0.;
  double coord_1 = 0.;
  int kg_1 = - ghost_size;
  for (int k = -ghost_size; k < nk + ghost_size; k++)
    {
      const int kg = k + offset;

      double coord;
      if (kg < 0)
        {
          coord = elem_coord[0];
        }
      else if (kg >= nktot)
        {
          coord = elem_coord[nktot - 1];
        }
      else
        {
          coord = elem_coord[kg];
        }

      while (coord_1 < coord)
        {
          d0 = d1;
          coord_0 = coord_1;

          if (kg_1 < 0)
            {
              d1 = 0.;
              coord_1 = coord_z_pour_delta[0];
            }
          else if (kg_1 >= nktot_pour_delta)
            {
              d1 = 0.;
              coord_1 = coord_z_pour_delta[nktot_pour_delta];
            }
          else
            {
              d1 = global_delta_pour_delta[kg_1];
              coord_1 = elem_coord_pour_delta[kg_1];
            }

          kg_1++;
        }

      if (coord_1 == coord_0)
        {
          if (coord_1 != coord || d0 != d1)
            {
              Cerr << "Ce n'est pas normal." << finl;
              Process::exit();
            }
          else
            {
              delta_z_pour_delta[k] = d0;
            }
        }
      else
        {
          delta_z_pour_delta[k] = d0 + (d1 - d0)/(coord_1 - coord_0)*(coord - coord_0);
        }
    }
}

void calculer_delta_z_filtre_identique(const IJK_Splitting& splitting,
                                       const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                       ArrOfDouble_with_ghost& delta_z_filtre)
{
  const int nk = splitting.get_nb_items_local(IJK_Splitting::FACES_K, DIRECTION_K);
  const int nktot = splitting.get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = splitting.get_offset_local(DIRECTION_K);

  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;
      const double dz_pour_delta = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];
      delta_z_filtre[k] = sqrt(dz_pour_delta*dz_pour_delta + dz_pour_delta*dz_pour_delta);
    }
}

void calculer_delta_z_filtre_n_mailles(const int taille_filtre,
                                       const IJK_Splitting& splitting,
                                       const ArrOfDouble_with_ghost& delta_z,
                                       const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                       ArrOfDouble_with_ghost& delta_z_filtre)
{
  Cerr << "Taille du filtre explicite fixee a " << taille_filtre << " mailles." << finl;

  const int nk = splitting.get_nb_items_local(IJK_Splitting::FACES_K, DIRECTION_K);
  const int nktot = splitting.get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = splitting.get_offset_local(DIRECTION_K);

  const bool impair = (taille_filtre%2 != 0);
  const int kp_min = impair ? (taille_filtre-1)/2 : (taille_filtre-2)/2;
  const int kp_max = impair ? (taille_filtre-1)/2 : (taille_filtre-2)/2;

  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      double longueur_elem = 0.;
      for (int kp = -kp_min; kp < kp_max+1; kp++)
        {
          const int kpg = kg + kp;
          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
          longueur_elem += dz;
        }

      if (!impair)
        {
          {
            const int kp = -kp_min - 1;
            const int kpg = kg + kp;
            const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
            longueur_elem += dz;
          }
          {
            const int kp = kp_max + 1;
            const int kpg = kg + kp;
            const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
            longueur_elem += dz;
          }
        }
      const double dz_pour_delta = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];
      delta_z_filtre[k] = sqrt(longueur_elem*longueur_elem + dz_pour_delta*dz_pour_delta);
    }
}

// On suppose que rho a un espace virtuel a jour sur au moins l'epaisseur de joint demandee.
void calculer_temperature_mu_lambda_air(const double P_th,
                                        const double constante_specifique_gaz,
                                        const IJK_Field_double& rho,
                                        IJK_Field_double& temperature,
                                        IJK_Field_double& molecular_mu,
                                        IJK_Field_double& molecular_lambda,
                                        int epaisseur_joint)
{
  const int imax = rho.ni() + epaisseur_joint;
  const int jmax = rho.nj() + epaisseur_joint;
  const int kmax = rho.nk() + epaisseur_joint;

  const double fac_a = -5.05628e-18;
  const double fac_b = 2.469e-14;
  const double fac_c = -4.98344e-11;
  const double fac_d = 7.06714e-08;
  const double fac_e = 1.0894e-06;

  const double facteur = 1.93198026315789000e-3 / 1.461e-6;

  const double facteur_temperature = P_th / constante_specifique_gaz;

  for (int k = -epaisseur_joint; k < kmax; k++)
    {
      for (int j = -epaisseur_joint; j < jmax; j++)
        {
          for (int i = -epaisseur_joint; i < imax; i++)
            {
              const double r = rho(i, j, k);
              const double temp = facteur_temperature / r;
              temperature(i,j,k) = temp;
              double calc = temp * fac_a + fac_b;
              calc = temp * calc  + fac_c;
              calc = temp * calc  + fac_d;
              calc = temp * calc  + fac_e;
              molecular_lambda(i,j,k) = calc * facteur;
              molecular_mu(i,j,k) = calc;
            }
        }
    }
}

void multiplier_champ(const double coefficient,
                      const IJK_Field_double& champ_in,
                      IJK_Field_double& champ_out)
{
  const int ni = champ_out.ni();
  const int nj = champ_out.nj();
  const int nk = champ_out.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              champ_out(i,j,k) = coefficient * champ_in(i,j,k);
            }
        }
    }
}

void multiplier_champ(const IJK_Field_double& champ_rho,
                      IJK_Field_double& champ)
{
  const int ni = champ.ni();
  const int nj = champ.nj();
  const int nk = champ.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              champ(i,j,k) = champ_rho(i,j,k) * champ(i,j,k);
            }
        }
    }
}



// Adds field1 to field2 by reference
void add_fields(const IJK_Field_double& field1, IJK_Field_double& field2)
{
  const int ni = field2.ni();
  const int nj = field2.nj();
  const int nk = field2.nk();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              field2(i,j,k) += field1(i,j,k);
            }
        }
    }
}

// Adds field1 to field2 by reference, on the k^{th} layer of the channel
void add_fields_k(const IJK_Field_double& field1, IJK_Field_double& field2, const int k)
{
  const int ni = field2.ni();
  const int nj = field2.nj();
  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          field2(i,j,k) += field1(i,j,k);
        }
    }
}

void multiplier_champ_rho_arete_ij(const IJK_Field_double& champ_rho,
                                   IJK_Field_double& champ)
{
  const int ni = champ.ni();
  const int nj = champ.nj();
  const int nk = champ.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const double rho = champ_rho(i,j,k);
              const double rho_im1 = champ_rho(i-1,j,k);
              const double rho_jm1 = champ_rho(i,j-1,k);
              const double rho_im1_jm1 = champ_rho(i-1,j-1,k);
              const double rho_ij = 0.25 * (rho + rho_im1 + rho_jm1 + rho_im1_jm1);
              champ(i,j,k) = rho_ij * champ(i,j,k);
            }
        }
    }
}

void multiplier_champ_rho_arete_ik(const IJK_Field_double& champ_rho,
                                   double rho_kmin,
                                   IJK_Field_double& champ)
{
  const int offset = champ_rho.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = champ.ni();
  const int nj = champ.nj();
  const int nk = champ.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;
              const double rho = champ_rho(i,j,k);
              const double rho_im1 = champ_rho(i-1,j,k);
              const double rho_ik = kg==0 ? rho_kmin : 0.25 * (rho + rho_im1 + champ_rho(i,j,k-1) + champ_rho(i-1,j,k-1));
              champ(i,j,k) = rho_ik * champ(i,j,k);
            }
        }
    }
}

void multiplier_champ_rho_arete_jk(const IJK_Field_double& champ_rho,
                                   double rho_kmin,
                                   IJK_Field_double& champ)
{
  const int offset = champ_rho.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = champ.ni();
  const int nj = champ.nj();
  const int nk = champ.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;
              const double rho = champ_rho(i,j,k);
              const double rho_jm1 = champ_rho(i,j-1,k);
              const double rho_jk = kg==0 ? rho_kmin : 0.25 * (rho + rho_jm1 + champ_rho(i,j,k-1) + champ_rho(i,j-1,k-1));
              champ(i,j,k) = rho_jk * champ(i,j,k);
            }
        }
    }
}

void multiplier_champ(const IJK_Field_double& champ_rho,
                      const double& constant,
                      IJK_Field_double& champ)
{
  const int ni = champ.ni();
  const int nj = champ.nj();
  const int nk = champ.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              champ(i,j,k) = champ_rho(i,j,k) * champ(i,j,k) * constant;
            }
        }
    }
}

void multiplier_champ_rho_face_i(const int flag_add,
                                 const IJK_Field_double& champ_rho,
                                 const double& constant,
                                 const IJK_Field_double& champ_in,
                                 IJK_Field_double& champ_out)
{
  const int ni = champ_in.ni();
  const int nj = champ_in.nj();
  const int nk = champ_in.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const double rho = champ_rho(i,j,k);
              const double rho_im1 = champ_rho(i-1,j,k);
              const double rhof_i = 0.5 * (rho + rho_im1);
              if (flag_add)
                {
                  champ_out(i,j,k) += rhof_i * champ_in(i,j,k) * constant;
                }
              else
                {
                  champ_out(i,j,k) = rhof_i * champ_in(i,j,k) * constant;
                }
            }
        }
    }
}

void multiplier_champ_rho_face_j(const int flag_add,
                                 const IJK_Field_double& champ_rho,
                                 const double& constant,
                                 const IJK_Field_double& champ_in,
                                 IJK_Field_double& champ_out)
{
  const int ni = champ_in.ni();
  const int nj = champ_in.nj();
  const int nk = champ_in.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const double rho = champ_rho(i,j,k);
              const double rho_jm1 = champ_rho(i,j-1,k);
              const double rhof_j = 0.5 * (rho + rho_jm1);
              if (flag_add)
                {
                  champ_out(i,j,k) += rhof_j * champ_in(i,j,k) * constant;
                }
              else
                {
                  champ_out(i,j,k) = rhof_j * champ_in(i,j,k) * constant;
                }
            }
        }
    }
}

void multiplier_champ_rho_face_k(const int flag_add,
                                 const IJK_Field_double& champ_rho,
                                 const double& constant,
                                 double rho_kmin,
                                 const IJK_Field_double& champ_in,
                                 IJK_Field_double& champ_out)
{
  const int offset = champ_rho.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = champ_in.ni();
  const int nj = champ_in.nj();
  const int nk = champ_in.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;
              const double rho = champ_rho(i,j,k);
              const double rhof_k = kg==0 ? rho_kmin : 0.5 * (rho + champ_rho(i,j,k-1));
              if (flag_add)
                {
                  champ_out(i,j,k) += rhof_k * champ_in(i,j,k) * constant;
                }
              else
                {
                  champ_out(i,j,k) = rhof_k * champ_in(i,j,k) * constant;
                }
            }
        }
    }
}

inline void multiplier_champ_rho_face_i(const IJK_Field_double& champ_rho,
                                        const double& constant,
                                        IJK_Field_double& champ)
{
  multiplier_champ_rho_face_i(false, champ_rho, constant, champ, champ);
}

inline void multiplier_champ_rho_face_j(const IJK_Field_double& champ_rho,
                                        const double& constant,
                                        IJK_Field_double& champ)
{
  multiplier_champ_rho_face_j(false, champ_rho, constant, champ, champ);
}

inline void multiplier_champ_rho_face_k(const IJK_Field_double& champ_rho,
                                        const double& constant,
                                        double rho_kmin,
                                        IJK_Field_double& champ)
{
  multiplier_champ_rho_face_k(false, champ_rho, constant, rho_kmin, champ, champ);
}

template<class T>
void filtrer_champ_elem(const int flag_add,
                        const IJK_Field_double& champ,
                        const ArrOfDouble_with_ghost& delta_z,
                        const double facteur_delta_x,
                        const double facteur_delta_y,
                        const ArrOfDouble_with_ghost& delta_z_pour_delta,
                        T& kernel,
                        FixedVector<IJK_Field_local_double, 18>& tmp_b,
                        FixedVector<IJK_Field_local_double, 18>& tmp_a,
                        IJK_Field_double& champ_filtre)
{
  const double dx = champ.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = champ.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = champ.ni();
  const int nj = champ.nj();
  const int nk = champ.nk();

  const int nktot = champ.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = champ.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& b = tmp_b[0];
  IJK_Field_local_double& a = tmp_a[0];

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];

      const FixedVector<double, 21> filter_kernel_z = kernel->inhomogeneous(true, k, kg, nktot, dz_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_elem = kernel->size_k_elem(kg, nktot);
      const int shift_k_elem = kernel->shift_k_elem(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_elem = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_elem = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<-1 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double filter_coef_z = filter_kernel_z[kp+10];
                  longueur_elem += filter_coef_z * dz;
                }
              facteur_elem = 1./longueur_elem;
            }
          else
            {
              facteur_elem = dz_glo==0. ? 0. : 1./dz_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              b(i, j, 0) = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_elem(kpg, nktot)))
                    {
                      const double filter_coef_z = filter_kernel_z[kp+10];
                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          b(i, j, 0) += champ(i,j,k+kp) * filter_coef_z * dz * facteur_elem;
                        }
                      else
                        {
                          b(i, j, 0) += champ(i,j,k+kp) * filter_coef_z;
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              a(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  a(i, 0, 0) += b(i, j+jp, 0) * filter_coef_y;
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double r = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  r += a(i+ip, 0, 0) * filter_coef_x;
                }

              if (flag_add)
                {
                  champ_filtre(i,j,k) += r;
                }
              else
                {
                  champ_filtre(i,j,k) = r;
                }
            }
        }
    }
}

template<class T>
void filtrer_champ_face(const int flag_add,
                        const IJK_Field_double& champ,
                        const ArrOfDouble_with_ghost& delta_z,
                        const double facteur_delta_x,
                        const double facteur_delta_y,
                        const ArrOfDouble_with_ghost& delta_z_pour_delta,
                        T& kernel,
                        FixedVector<IJK_Field_local_double, 18>& tmp_b,
                        FixedVector<IJK_Field_local_double, 18>& tmp_a,
                        IJK_Field_double& champ_filtre)
{
  const double dx = champ.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = champ.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = champ.ni();
  const int nj = champ.nj();
  const int nk = champ.nk();

  const int nktot = champ.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = champ.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& b = tmp_b[0];
  IJK_Field_local_double& a = tmp_a[0];

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_m1_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z[k-1];
      const double delta_m_glo = kg==0 ? 0.5*dz_glo : 0.5*(dz_glo + dz_m1_glo);
      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];
      const double dz_m1_pour_delta_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z_pour_delta[k-1];
      const double delta_m_pour_delta_glo = (kg-1<0 || kg>(nktot-1)) ? 0. : 0.5*(dz_pour_delta_glo + dz_m1_pour_delta_glo);

      const FixedVector<double, 21> filter_kernel_z_face = kernel->inhomogeneous(false, k, kg, nktot, delta_m_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_face = kernel->size_k_face(kg, nktot);
      const int shift_k_face = kernel->shift_k_face(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_face = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_face = 0.;
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<0 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                  const double dzf = 0.5*(dz + dz_m1);
                  const double filter_coef_z_face = filter_kernel_z_face[kp+10];
                  longueur_face += filter_coef_z_face * dzf;
                }
              facteur_face = 1./longueur_face;
            }
          else
            {
              facteur_face = delta_m_glo==0. ? 0. : 1./delta_m_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              b(i, j, 0) = 0.;
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_face(kpg, nktot)))
                    {
                      const double filter_coef_z_face = filter_kernel_z_face[kp+10];
                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                          const double dzf = 0.5*(dz + dz_m1);
                          b(i, j, 0) += champ(i,j,k+kp) * filter_coef_z_face * dzf * facteur_face;
                        }
                      else
                        {
                          b(i, j, 0) += champ(i,j,k+kp) * filter_coef_z_face;
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              a(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  a(i, 0, 0) += b(i, j+jp, 0) * filter_coef_y;
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double r = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  r += a(i+ip, 0, 0) * filter_coef_x;
                }

              if (flag_add)
                {
                  champ_filtre(i,j,k) += r;
                }
              else
                {
                  champ_filtre(i,j,k) = r;
                }
            }
        }
    }
}

static inline void calculer_g(int i, int j, int k,
                              const FixedVector<IJK_Field_double, 3>& velocity,
                              const ArrOfDouble_with_ghost& delta_z_maillage,
                              const double dx_delta,
                              const double dy_delta,
                              const double dz_delta,
                              const double delta_m_delta,
                              const double delta_p_delta,
                              FixedVector<FixedVector<double, 3>, 3>& g)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double deltaunsurdx = dx_delta * 1./dx;
  const double deltaunsurdy = dy_delta * 1./dy;

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int kg = k + offset;
  const double dz = delta_z_maillage[k];
  const double delta_m = kg==0 ? 0.5*dz : 0.5*(dz + delta_z_maillage[k-1]);
  const double delta_p = kg==(nktot-1) ? 0.5*dz : 0.5*(dz + delta_z_maillage[k+1]);

  const double deltaunsurdz = dz_delta * 1./dz;
  const double deltaunsurdelta_m = delta_m_delta * 1./delta_m;
  const double deltaunsurdelta_p = delta_p_delta * 1./delta_p;

  const double uf_i = vitesse_i(i,j,k);
  const double uf_ip1 = vitesse_i(i+1,j,k);
  const double uf_i_jm1 = vitesse_i(i,j-1,k);
  const double uf_i_jp1 = vitesse_i(i,j+1,k);
  const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);
  const double uf_i_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i,j,k+1);
  const double uf_ip1_jp1 = vitesse_i(i+1,j+1,k);
  const double uf_ip1_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i+1,j,k+1);
  const double uf_ip1_jm1 = vitesse_i(i+1,j-1,k);
  const double uf_ip1_km1 = kg==0 ? 0. : vitesse_i(i+1,j,k-1);

  const double vf_j = vitesse_j(i,j,k);
  const double vf_j_im1 = vitesse_j(i-1,j,k);
  const double vf_j_ip1 = vitesse_j(i+1,j,k);
  const double vf_jp1 = vitesse_j(i,j+1,k);
  const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);
  const double vf_j_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j,k+1);
  const double vf_jp1_ip1 = vitesse_j(i+1,j+1,k);
  const double vf_jp1_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j+1,k+1);
  const double vf_jp1_im1 = vitesse_j(i-1,j+1,k);
  const double vf_jp1_km1 = kg==0 ? 0. : vitesse_j(i,j+1,k-1);

  const double wf_k = vitesse_k(i,j,k);
  const double wf_k_im1 = vitesse_k(i-1,j,k);
  const double wf_k_ip1 = vitesse_k(i+1,j,k);
  const double wf_k_jm1 = vitesse_k(i,j-1,k);
  const double wf_k_jp1 = vitesse_k(i,j+1,k);
  const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);
  const double wf_kp1_ip1 = kg==(nktot-1) ? 0. : vitesse_k(i+1,j,k+1);
  const double wf_kp1_jp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j+1,k+1);
  const double wf_kp1_im1 = kg==(nktot-1) ? 0. : vitesse_k(i-1,j,k+1);
  const double wf_kp1_jm1 = kg==(nktot-1) ? 0. : vitesse_k(i,j-1,k+1);

  const double duidx = deltaunsurdx * (uf_ip1 - uf_i);
  const double dujdy = deltaunsurdy * (vf_jp1 - vf_j);
  const double dukdz = deltaunsurdz * (wf_kp1 - wf_k);

  const double duidy_ij = deltaunsurdy * (uf_i - uf_i_jm1);
  const double duidy_ip1j = deltaunsurdy * (uf_ip1 - uf_ip1_jm1);
  const double duidy_ijp1 = deltaunsurdy * (uf_i_jp1 - uf_i);
  const double duidy_ip1jp1 = deltaunsurdy * (uf_ip1_jp1 - uf_ip1);

  const double duidz_ik = deltaunsurdelta_m * (uf_i - uf_i_km1);
  const double duidz_ip1k = deltaunsurdelta_m * (uf_ip1 - uf_ip1_km1);
  const double duidz_ikp1 = deltaunsurdelta_p * (uf_i_kp1 - uf_i);
  const double duidz_ip1kp1 = deltaunsurdelta_p * (uf_ip1_kp1 - uf_ip1);

  const double dujdx_ij = deltaunsurdx * (vf_j - vf_j_im1);
  const double dujdx_ip1j = deltaunsurdx * (vf_j_ip1 - vf_j);
  const double dujdx_ijp1 = deltaunsurdx * (vf_jp1 - vf_jp1_im1);
  const double dujdx_ip1jp1 = deltaunsurdx * (vf_jp1_ip1 - vf_jp1);

  const double dujdz_jk = deltaunsurdelta_m * (vf_j - vf_j_km1);
  const double dujdz_jp1k = deltaunsurdelta_m * (vf_jp1 - vf_jp1_km1);
  const double dujdz_jkp1 = deltaunsurdelta_p * (vf_j_kp1 - vf_j);
  const double dujdz_jp1kp1 = deltaunsurdelta_p * (vf_jp1_kp1 - vf_jp1);

  const double dukdx_ik = deltaunsurdx * (wf_k - wf_k_im1);
  const double dukdx_ip1k = deltaunsurdx * (wf_k_ip1 - wf_k);
  const double dukdx_ikp1 = deltaunsurdx * (wf_kp1 - wf_kp1_im1);
  const double dukdx_ip1kp1 = deltaunsurdx * (wf_kp1_ip1 - wf_kp1);

  const double dukdy_jk = deltaunsurdy * (wf_k - wf_k_jm1);
  const double dukdy_jp1k = deltaunsurdy * (wf_k_jp1 - wf_k);
  const double dukdy_jkp1 = deltaunsurdy * (wf_kp1 - wf_kp1_jm1);
  const double dukdy_jp1kp1 = deltaunsurdy * (wf_kp1_jp1 - wf_kp1);

  g[0][0] = duidx;
  g[0][1] = 0.25 * (duidy_ip1jp1 + duidy_ijp1 + duidy_ip1j + duidy_ij);
  g[0][2] = 0.25 * (duidz_ip1kp1 + duidz_ikp1 + duidz_ip1k + duidz_ik);
  g[1][0] = 0.25 * (dujdx_ip1jp1 + dujdx_ip1j + dujdx_ijp1 + dujdx_ij);
  g[1][1] = dujdy;
  g[1][2] = 0.25 * (dujdz_jp1kp1 + dujdz_jkp1 + dujdz_jp1k + dujdz_jk);
  g[2][0] = 0.25 * (dukdx_ip1kp1 + dukdx_ip1k + dukdx_ikp1 + dukdx_ik);
  g[2][1] = 0.25 * (dukdy_jp1kp1 + dukdy_jp1k + dukdy_jkp1 + dukdy_jk);
  g[2][2] = dukdz;
}

static inline void calculer_tau(int i, int j, int k,
                                const FixedVector<IJK_Field_double, 3>& velocity,
                                const IJK_Field_double& turbulent_mu_xx,
                                const IJK_Field_double& turbulent_mu_xy,
                                const IJK_Field_double& turbulent_mu_xz,
                                const IJK_Field_double& turbulent_mu_yy,
                                const IJK_Field_double& turbulent_mu_yz,
                                const IJK_Field_double& turbulent_mu_zz,
                                const ArrOfDouble_with_ghost& delta_z_maillage,
                                const double dx_delta,
                                const double dy_delta,
                                const double dz_delta,
                                const double delta_m_delta,
                                const double delta_p_delta,
                                FixedVector<FixedVector<double, 3>, 3>& f)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double deltaunsurdx = dx_delta * 1./dx;
  const double deltaunsurdy = dy_delta * 1./dy;

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int kg = k + offset;
  const double dz = delta_z_maillage[k];
  const double delta_m = kg==0 ? 0.5*dz : 0.5*(dz + delta_z_maillage[k-1]);
  const double delta_p = kg==(nktot-1) ? 0.5*dz : 0.5*(dz + delta_z_maillage[k+1]);

  const double deltaunsurdz = dz_delta * 1./dz;
  const double deltaunsurdelta_m = delta_m_delta * 1./delta_m;
  const double deltaunsurdelta_p = delta_p_delta * 1./delta_p;

  const double nu_xx = turbulent_mu_xx(i,j,k);

  const double nu_xy = turbulent_mu_xy(i,j,k);
  const double nu_xy_im1 = turbulent_mu_xy(i-1,j,k);
  const double nu_xy_jm1 = turbulent_mu_xy(i,j-1,k);
  const double nu_xy_im1_jm1 = turbulent_mu_xy(i-1,j-1,k);
  const double nu_xy_ip1 = turbulent_mu_xy(i+1,j,k);
  const double nu_xy_jp1 = turbulent_mu_xy(i,j+1,k);
  const double nu_xy_ip1_jp1 = turbulent_mu_xy(i+1,j+1,k);
  const double nu_xy_im1_jp1 = turbulent_mu_xy(i-1,j+1,k);
  const double nu_xy_ip1_jm1 = turbulent_mu_xy(i+1,j-1,k);

  const double nu_xz = turbulent_mu_xz(i,j,k);
  const double nu_xz_im1 = turbulent_mu_xz(i-1,j,k);
  const double nu_xz_km1 = kg==0 ? 0. : turbulent_mu_xz(i,j,k-1);
  const double nu_xz_im1_km1 = kg==0 ? 0. : turbulent_mu_xz(i-1,j,k-1);
  const double nu_xz_ip1 = turbulent_mu_xz(i+1,j,k);
  const double nu_xz_kp1 = kg==(nktot-1) ? 0. : turbulent_mu_xz(i,j,k+1);
  const double nu_xz_ip1_kp1 = kg==(nktot-1) ? 0. : turbulent_mu_xz(i+1,j,k+1);
  const double nu_xz_im1_kp1 = kg==(nktot-1) ? 0. : turbulent_mu_xz(i-1,j,k+1);
  const double nu_xz_ip1_km1 = kg==0 ? 0. : turbulent_mu_xz(i+1,j,k-1);

  const double nu_yy = turbulent_mu_yy(i,j,k);

  const double nu_yz = turbulent_mu_yz(i,j,k);
  const double nu_yz_jm1 = turbulent_mu_yz(i,j-1,k);
  const double nu_yz_km1 = kg==0 ? 0. : turbulent_mu_yz(i,j,k-1);
  const double nu_yz_jm1_km1 = kg==0 ? 0. : turbulent_mu_yz(i,j-1,k-1);
  const double nu_yz_jp1 = turbulent_mu_yz(i,j+1,k);
  const double nu_yz_kp1 = kg==(nktot-1) ? 0. : turbulent_mu_yz(i,j,k+1);
  const double nu_yz_jp1_kp1 = kg==(nktot-1) ? 0. : turbulent_mu_yz(i,j+1,k+1);
  const double nu_yz_jm1_kp1 = kg==(nktot-1) ? 0. : turbulent_mu_yz(i,j-1,k+1);
  const double nu_yz_jp1_km1 = kg==0 ? 0. : turbulent_mu_yz(i,j+1,k-1);

  const double nu_xy_ij = 0.25 * (nu_xy + nu_xy_im1 + nu_xy_jm1 + nu_xy_im1_jm1);
  const double nu_xy_ip1j = 0.25 * (nu_xy_ip1 + nu_xy + nu_xy_ip1_jm1 + nu_xy_jm1);
  const double nu_xy_ijp1 = 0.25 * (nu_xy_jp1 + nu_xy_im1_jp1 + nu_xy + nu_xy_im1);
  const double nu_xy_ip1jp1 = 0.25 * (nu_xy_ip1_jp1 + nu_xy_jp1 + nu_xy_ip1 + nu_xy);
  const double nu_xz_ik = 0.25 * (nu_xz + nu_xz_im1 + nu_xz_km1 + nu_xz_im1_km1);
  const double nu_xz_ip1k = 0.25 * (nu_xz_ip1 + nu_xz + nu_xz_ip1_km1 + nu_xz_km1);
  const double nu_xz_ikp1 = 0.25 * (nu_xz_kp1 + nu_xz_im1_kp1 + nu_xz + nu_xz_im1);
  const double nu_xz_ip1kp1 = 0.25 * (nu_xz_ip1_kp1 + nu_xz_kp1 + nu_xz_ip1 + nu_xz);
  const double nu_yz_jk = 0.25 * (nu_yz + nu_yz_jm1 + nu_yz_km1 + nu_yz_jm1_km1);
  const double nu_yz_jp1k = 0.25 * (nu_yz_jp1 + nu_yz + nu_yz_jp1_km1 + nu_yz_km1);
  const double nu_yz_jkp1 = 0.25 * (nu_yz_kp1 + nu_yz_jm1_kp1 + nu_yz + nu_yz_jm1);
  const double nu_yz_jp1kp1 = 0.25 * (nu_yz_jp1_kp1 + nu_yz_kp1 + nu_yz_jp1 + nu_yz);

  const double nu_zz = turbulent_mu_zz(i,j,k);

  const double uf_i = vitesse_i(i,j,k);
  const double uf_i_jm1 = vitesse_i(i,j-1,k);
  const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);

  const double vf_j = vitesse_j(i,j,k);
  const double vf_j_im1 = vitesse_j(i-1,j,k);
  const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);

  const double wf_k = vitesse_k(i,j,k);
  const double wf_k_im1 = vitesse_k(i-1,j,k);
  const double wf_k_jm1 = vitesse_k(i,j-1,k);

  const double uf_ip1 = vitesse_i(i+1,j,k);
  const double uf_i_jp1 = vitesse_i(i,j+1,k);
  const double uf_i_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i,j,k+1);
  const double uf_ip1_jp1 = vitesse_i(i+1,j+1,k);
  const double uf_ip1_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i+1,j,k+1);
  const double uf_ip1_jm1 = vitesse_i(i+1,j-1,k);
  const double uf_ip1_km1 = kg==0 ? 0. : vitesse_i(i+1,j,k-1);

  const double vf_j_ip1 = vitesse_j(i+1,j,k);
  const double vf_jp1 = vitesse_j(i,j+1,k);
  const double vf_j_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j,k+1);
  const double vf_jp1_ip1 = vitesse_j(i+1,j+1,k);
  const double vf_jp1_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j+1,k+1);
  const double vf_jp1_im1 = vitesse_j(i-1,j+1,k);
  const double vf_jp1_km1 = kg==0 ? 0. : vitesse_j(i,j+1,k-1);

  const double wf_k_ip1 = vitesse_k(i+1,j,k);
  const double wf_k_jp1 = vitesse_k(i,j+1,k);
  const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);
  const double wf_kp1_ip1 = kg==(nktot-1) ? 0. : vitesse_k(i+1,j,k+1);
  const double wf_kp1_jp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j+1,k+1);
  const double wf_kp1_im1 = kg==(nktot-1) ? 0. : vitesse_k(i-1,j,k+1);
  const double wf_kp1_jm1 = kg==(nktot-1) ? 0. : vitesse_k(i,j-1,k+1);

  const double duidx = deltaunsurdx * (uf_ip1 - uf_i);
  const double dujdy = deltaunsurdy * (vf_jp1 - vf_j);
  const double dukdz = deltaunsurdz * (wf_kp1 - wf_k);

  const double duidy_ij = deltaunsurdy * (uf_i - uf_i_jm1);
  const double duidy_ip1j = deltaunsurdy * (uf_ip1 - uf_ip1_jm1);
  const double duidy_ijp1 = deltaunsurdy * (uf_i_jp1 - uf_i);
  const double duidy_ip1jp1 = deltaunsurdy * (uf_ip1_jp1 - uf_ip1);

  const double duidz_ik = deltaunsurdelta_m * (uf_i - uf_i_km1);
  const double duidz_ip1k = deltaunsurdelta_m * (uf_ip1 - uf_ip1_km1);
  const double duidz_ikp1 = deltaunsurdelta_p * (uf_i_kp1 - uf_i);
  const double duidz_ip1kp1 = deltaunsurdelta_p * (uf_ip1_kp1 - uf_ip1);

  const double dujdx_ij = deltaunsurdx * (vf_j - vf_j_im1);
  const double dujdx_ip1j = deltaunsurdx * (vf_j_ip1 - vf_j);
  const double dujdx_ijp1 = deltaunsurdx * (vf_jp1 - vf_jp1_im1);
  const double dujdx_ip1jp1 = deltaunsurdx * (vf_jp1_ip1 - vf_jp1);

  const double dujdz_jk = deltaunsurdelta_m * (vf_j - vf_j_km1);
  const double dujdz_jp1k = deltaunsurdelta_m * (vf_jp1 - vf_jp1_km1);
  const double dujdz_jkp1 = deltaunsurdelta_p * (vf_j_kp1 - vf_j);
  const double dujdz_jp1kp1 = deltaunsurdelta_p * (vf_jp1_kp1 - vf_jp1);

  const double dukdx_ik = deltaunsurdx * (wf_k - wf_k_im1);
  const double dukdx_ip1k = deltaunsurdx * (wf_k_ip1 - wf_k);
  const double dukdx_ikp1 = deltaunsurdx * (wf_kp1 - wf_kp1_im1);
  const double dukdx_ip1kp1 = deltaunsurdx * (wf_kp1_ip1 - wf_kp1);

  const double dukdy_jk = deltaunsurdy * (wf_k - wf_k_jm1);
  const double dukdy_jp1k = deltaunsurdy * (wf_k_jp1 - wf_k);
  const double dukdy_jkp1 = deltaunsurdy * (wf_kp1 - wf_kp1_jm1);
  const double dukdy_jp1kp1 = deltaunsurdy * (wf_kp1_jp1 - wf_kp1);

  const double s_ii = 0.5 * (duidx + duidx);
  const double s_jj = 0.5 * (dujdy + dujdy);
  const double s_kk = 0.5 * (dukdz + dukdz);

  const double s_ij = 0.5 * (duidy_ij + dujdx_ij);
  const double s_ik = 0.5 * (duidz_ik + dukdx_ik);
  const double s_jk = 0.5 * (dujdz_jk + dukdy_jk);

  const double s_ip1j = 0.5 * (duidy_ip1j + dujdx_ip1j);
  const double s_ip1k = 0.5 * (duidz_ip1k + dukdx_ip1k);
  const double s_jp1k = 0.5 * (dujdz_jp1k + dukdy_jp1k);

  const double s_ijp1 = 0.5 * (duidy_ijp1 + dujdx_ijp1);
  const double s_ikp1 = 0.5 * (duidz_ikp1 + dukdx_ikp1);
  const double s_jkp1 = 0.5 * (dujdz_jkp1 + dukdy_jkp1);

  const double s_ip1jp1 = 0.5 * (duidy_ip1jp1 + dujdx_ip1jp1);
  const double s_ip1kp1 = 0.5 * (duidz_ip1kp1 + dukdx_ip1kp1);
  const double s_jp1kp1 = 0.5 * (dujdz_jp1kp1 + dukdy_jp1kp1);

  f[0][0] = - 2. * nu_xx * s_ii;
  f[0][1] = - 2. * 0.25 * (nu_xy_ip1jp1*s_ip1jp1 + nu_xy_ijp1*s_ijp1 + nu_xy_ip1j*s_ip1j + nu_xy_ij*s_ij);
  f[0][2] = - 2. * 0.25 * (nu_xz_ip1kp1*s_ip1kp1 + nu_xz_ikp1*s_ikp1 + nu_xz_ip1k*s_ip1k + nu_xz_ik*s_ik);
  f[1][0] = - 2. * 0.25 * (nu_xy_ip1jp1*s_ip1jp1 + nu_xy_ijp1*s_ijp1 + nu_xy_ip1j*s_ip1j + nu_xy_ij*s_ij);
  f[1][1] = - 2. * nu_yy * s_jj;
  f[1][2] = - 2. * 0.25 * (nu_yz_jp1kp1*s_jp1kp1 + nu_yz_jkp1*s_jkp1 + nu_yz_jp1k*s_jp1k + nu_yz_jk*s_jk);
  f[2][0] = - 2. * 0.25 * (nu_xz_ip1kp1*s_ip1kp1 + nu_xz_ikp1*s_ikp1 + nu_xz_ip1k*s_ip1k + nu_xz_ik*s_ik);
  f[2][1] = - 2. * 0.25 * (nu_yz_jp1kp1*s_jp1kp1 + nu_yz_jkp1*s_jkp1 + nu_yz_jp1k*s_jp1k + nu_yz_jk*s_jk);
  f[2][2] = - 2. * nu_zz * s_kk;
}

static inline void calculer_q(int i, int j, int k,
                              const FixedVector<IJK_Field_double, 3>& velocity,
                              const IJK_Field_double& champ_scalar,
                              double scalar_kmin, double scalar_kmax,
                              const ArrOfDouble_with_ghost& delta_z_maillage,
                              const double dx_delta,
                              const double dy_delta,
                              const double dz_delta,
                              const double delta_m_delta,
                              const double delta_p_delta,
                              FixedVector<double, 3>& q)
{
  const IJK_Field_double& vitesse_k = velocity[2];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double deltaunsurdx = dx_delta * 1./dx;
  const double deltaunsurdy = dy_delta * 1./dy;

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int kg = k + offset;
  const double dz = delta_z_maillage[k];
  const double delta_m = kg==0 ? 0.5*dz : 0.5*(dz + delta_z_maillage[k-1]);
  const double delta_p = kg==(nktot-1) ? 0.5*dz : 0.5*(dz + delta_z_maillage[k+1]);

  const double deltaunsurdelta_m = delta_m_delta * 1./delta_m;
  const double deltaunsurdelta_p = delta_p_delta * 1./delta_p;

  const double scalar = champ_scalar(i,j,k);

  const double scalar_im1 = champ_scalar(i-1,j,k);
  const double scalar_ip1 = champ_scalar(i+1,j,k);

  const double scalar_jm1 = champ_scalar(i,j-1,k);
  const double scalar_jp1 = champ_scalar(i,j+1,k);

  const double scalar_km1 = kg==0 ? scalar_kmin : champ_scalar(i,j,k-1);
  const double scalar_kp1 = kg==(nktot-1) ? scalar_kmax : champ_scalar(i,j,k+1);

  const double dscalardxf_i = deltaunsurdx * (scalar - scalar_im1);
  const double dscalardxf_ip1 = deltaunsurdx * (scalar_ip1 - scalar);

  const double dscalardyf_j = deltaunsurdy * (scalar - scalar_jm1);
  const double dscalardyf_jp1 = deltaunsurdy * (scalar_jp1 - scalar);

  const double dscalardzf_k = deltaunsurdelta_m * (scalar - scalar_km1);
  const double dscalardzf_kp1 = deltaunsurdelta_p * (scalar_kp1 - scalar);

  q[0] = 0.5 * (dscalardxf_i + dscalardxf_ip1);
  q[1] = 0.5 * (dscalardyf_j + dscalardyf_jp1);
  q[2] = 0.5 * (dscalardzf_k + dscalardzf_kp1);
}

static inline void calculer_pi(int i, int j, int k,
                               const FixedVector<IJK_Field_double, 3>& velocity,
                               const IJK_Field_double& champ_scalar,
                               double scalar_kmin, double scalar_kmax,
                               const IJK_Field_double& turbulent_mu_x,
                               const IJK_Field_double& turbulent_mu_y,
                               const IJK_Field_double& turbulent_mu_z,
                               const ArrOfDouble_with_ghost& delta_z_maillage,
                               const double dx_delta,
                               const double dy_delta,
                               const double dz_delta,
                               const double delta_m_delta,
                               const double delta_p_delta,
                               FixedVector<double, 3>& f)
{
  const IJK_Field_double& vitesse_k = velocity[2];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double deltaunsurdx = dx_delta * 1./dx;
  const double deltaunsurdy = dy_delta * 1./dy;

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int kg = k + offset;
  const double dz = delta_z_maillage[k];
  const double delta_m = kg==0 ? 0.5*dz : 0.5*(dz + delta_z_maillage[k-1]);
  const double delta_p = kg==(nktot-1) ? 0.5*dz : 0.5*(dz + delta_z_maillage[k+1]);

  const double deltaunsurdelta_m = delta_m_delta * 1./delta_m;
  const double deltaunsurdelta_p = delta_p_delta * 1./delta_p;

  const double nu_x = turbulent_mu_x(i,j,k);
  const double nu_x_im1 = turbulent_mu_x(i-1,j,k);
  const double nu_x_ip1 = turbulent_mu_x(i+1,j,k);

  const double nu_y = turbulent_mu_y(i,j,k);
  const double nu_y_jm1 = turbulent_mu_y(i,j-1,k);
  const double nu_y_jp1 = turbulent_mu_y(i,j+1,k);

  const double nu_z = turbulent_mu_z(i,j,k);
  const double nu_z_km1 = kg==0 ? 0. : turbulent_mu_z(i,j,k-1);
  const double nu_z_kp1 = kg==(nktot-1) ? 0. : turbulent_mu_z(i,j,k+1);

  const double nu_xf_i = 0.5 * (nu_x + nu_x_im1);
  const double nu_yf_j = 0.5 * (nu_y + nu_y_jm1);
  const double nu_zf_k = 0.5 * (nu_z + nu_z_km1);
  const double nu_xf_ip1 = 0.5 * (nu_x + nu_x_ip1);
  const double nu_yf_jp1 = 0.5 * (nu_y + nu_y_jp1);
  const double nu_zf_kp1 = 0.5 * (nu_z + nu_z_kp1);

  const double scalar = champ_scalar(i,j,k);

  const double scalar_im1 = champ_scalar(i-1,j,k);
  const double scalar_ip1 = champ_scalar(i+1,j,k);

  const double scalar_jm1 = champ_scalar(i,j-1,k);
  const double scalar_jp1 = champ_scalar(i,j+1,k);

  const double scalar_km1 = kg==0 ? scalar_kmin : champ_scalar(i,j,k-1);
  const double scalar_kp1 = kg==(nktot-1) ? scalar_kmax : champ_scalar(i,j,k+1);

  const double dscalardxf_i = deltaunsurdx * (scalar - scalar_im1);
  const double dscalardxf_ip1 = deltaunsurdx * (scalar_ip1 - scalar);

  const double dscalardyf_j = deltaunsurdy * (scalar - scalar_jm1);
  const double dscalardyf_jp1 = deltaunsurdy * (scalar_jp1 - scalar);

  const double dscalardzf_k = deltaunsurdelta_m * (scalar - scalar_km1);
  const double dscalardzf_kp1 = deltaunsurdelta_p * (scalar_kp1 - scalar);

  f[0] = - 0.5 * (nu_xf_i*dscalardxf_i + nu_xf_ip1*dscalardxf_ip1);
  f[1] = - 0.5 * (nu_yf_j*dscalardyf_j + nu_yf_jp1*dscalardyf_jp1);
  f[2] = - 0.5 * (nu_zf_k*dscalardzf_k + nu_zf_kp1*dscalardzf_kp1);
}

template<class T>
void calculer_turbulent_mu(const bool anisotropic,
                           T& model,
                           const double turbulent_viscosity_model_constant,
                           const double variation_cste_modele_fonctionnel_,
                           const double smoothing_center_fr_,
                           const double smoothing_factor_fr_,
                           const double Re_tau_fr_,
                           const double Re_tau_ch_,
                           const double pond_fr_,
                           const double pond_ch_,
                           const double center_constant_,
                           const double Lz_tot_,
                           const FixedVector<IJK_Field_double, 3>& velocity,
                           const IJK_Field_double& champ_rho,
                           const IJK_Field_double& champ_scalar,
                           double scalar_kmin, double scalar_kmax,
                           const ArrOfDouble_with_ghost& delta_z_maillage,
                           const double facteur_x_pour_delta,
                           const double facteur_y_pour_delta,
                           const ArrOfDouble_with_ghost& delta_z_pour_delta,
                           IJK_Field_double& turbulent_mu,
                           const IJK_Splitting& splitting)
{
  const IJK_Field_double& vitesse_k = velocity[2];

  const double dx_pour_delta = facteur_x_pour_delta * vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy_pour_delta = facteur_y_pour_delta * vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  const ArrOfDouble& coord_z = splitting.get_grid_geometry().get_node_coordinates(DIRECTION_K);
  ArrOfDouble elem_coord(nk);
  for (int i = 0; i < nktot; i++)
    {
      elem_coord[i] = 0.5 * (coord_z[i] + coord_z[i+1]);
    }
  elem_coord[nktot]=Lz_tot_;

  FixedVector<FixedVector<double, 3>, 3> g;
  FixedVector<double, 3> q;
// Modif Martin
//  const int nktot_elem = splitting.get_grid_geometry().get_nb_elem_tot(DIRECTION_K);
  // Quick fix using the C++ std::vector
  // double turbulent_viscosity_model_constant_tab[nk];
  std::vector<double> turbulent_viscosity_model_constant_tab(nk);
  if ( !variation_cste_modele_fonctionnel_ )
    {
      if ( ( smoothing_factor_fr_ != 0. ) || ( smoothing_center_fr_ != 0. ) || pond_ch_ != 0. || pond_fr_ != 0. )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A non null smoothing_factor_fr or smoothing_center_fr_ or pond_ch_ or pond_fr_ has been specified while the variation_constante_modele_fonctionnel flag has not been activated."
               << "" << finl;
          Process::exit();
        }
      for (int k = 0; k <= nktot; k++)
        {
          turbulent_viscosity_model_constant_tab[k]=turbulent_viscosity_model_constant;
        }
    }
  else
    {
      if ( smoothing_factor_fr_  == 0. )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A null smoothing factor has been specified while the variation_constante_modele_fonctionnel flag has been activated. Check the 5 parameters used with this flag : smoothing_factor_fr_ smoothing_center_fr_, center_constant_, pond_fr_, pond_ch_."
               << "" << finl;
          Process::exit();
        }
      else if ((pond_fr_ == 0.) || (pond_ch_ == 0.))
        {
          Cerr << "Warning: pond_fr_ and/or pond_ch_ are equal to 0. The functional model constant is then not varying in the wall normal direction. Set pond_fr_ = 1. and pond_ch_ = 1. to obtain a symmetric variation of the model constant." ;
        }
      for (int k = 0; k < nktot/2; k++)
        {
          turbulent_viscosity_model_constant_tab[k]=turbulent_viscosity_model_constant*pond_fr_+(0.5+0.5*tanh((elem_coord[k]-smoothing_center_fr_)/(smoothing_factor_fr_)))*(center_constant_-turbulent_viscosity_model_constant*pond_fr_);
        }
      for (int k = nktot; k > nktot/2-1; --k)
        {
          turbulent_viscosity_model_constant_tab[k]=turbulent_viscosity_model_constant*pond_ch_+(0.5+0.5*tanh((Lz_tot_-elem_coord[k]-smoothing_center_fr_ * Re_tau_fr_/Re_tau_ch_)/(smoothing_factor_fr_ * Re_tau_fr_/Re_tau_ch_)))*(center_constant_-turbulent_viscosity_model_constant*pond_ch_);
        }
    }
// Fin modif Martin

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;
              const double dz_pour_delta = delta_z_pour_delta[k];

              const double rho = champ_rho(i,j,k);

              if (!anisotropic)
                {
                  calculer_g(i, j, k, velocity, delta_z_maillage, 1., 1., 1., 1., 1., g);
                  calculer_q(i, j, k, velocity, champ_scalar, scalar_kmin, scalar_kmax, delta_z_maillage, 1., 1., 1., 1., 1., q);
// Modif Martin
                  turbulent_mu(i,j,k) = model(turbulent_viscosity_model_constant_tab[k], dx_pour_delta, dy_pour_delta, dz_pour_delta, rho, g, q);
// Fin modif Martin
                }
              else
                {
                  const double delta_m_pour_delta = kg==0 ? 0. : 0.5*(dz_pour_delta + delta_z_pour_delta[k-1]);
                  const double delta_p_pour_delta = kg==(nktot-1) ? 0. : 0.5*(dz_pour_delta + delta_z_pour_delta[k+1]);
                  calculer_g(i, j, k, velocity, delta_z_maillage, dx_pour_delta, dy_pour_delta, dz_pour_delta, delta_m_pour_delta, delta_p_pour_delta, g);
                  calculer_q(i, j, k, velocity, champ_scalar, scalar_kmin, scalar_kmax, delta_z_maillage, dx_pour_delta, dy_pour_delta, dz_pour_delta, delta_m_pour_delta, delta_p_pour_delta, q);
// Modif Martin
                  turbulent_mu(i,j,k) = model(turbulent_viscosity_model_constant_tab[k], 1., 1., 1., rho, g, q);
// Fin modif Martin
                }
            }
        }
    }
}

template<class T>
void calculer_ml_dynamic_uu_tensor(const bool anisotropic,
                                   const bool tensorial,
                                   const FixedVector<IJK_Field_double, 3>& velocity,
                                   const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                   const int turbulent_viscosity,
                                   const IJK_Field_double& turbulent_mu_xx,
                                   const IJK_Field_double& turbulent_mu_xy,
                                   const IJK_Field_double& turbulent_mu_xz,
                                   const IJK_Field_double& turbulent_mu_yy,
                                   const IJK_Field_double& turbulent_mu_yz,
                                   const IJK_Field_double& turbulent_mu_zz,
                                   const IJK_Field_double& turbulent_mu_filtre_xx,
                                   const IJK_Field_double& turbulent_mu_filtre_xy,
                                   const IJK_Field_double& turbulent_mu_filtre_xz,
                                   const IJK_Field_double& turbulent_mu_filtre_yy,
                                   const IJK_Field_double& turbulent_mu_filtre_yz,
                                   const IJK_Field_double& turbulent_mu_filtre_zz,
                                   const int structural_uu,
                                   const FixedVector<IJK_Field_double, 6>& structural_uu_tensor,
                                   const FixedVector<IJK_Field_double, 6>& structural_uu_filtre_tensor,
                                   const ArrOfDouble_with_ghost& delta_z,
                                   const double facteur_delta_x,
                                   const double facteur_delta_y,
                                   const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                   const double facteur_delta_filtre_x,
                                   const double facteur_delta_filtre_y,
                                   const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                   T& kernel,
                                   FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                   FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                   FixedVector<FixedVector<ArrOfDouble, 7>, 8>& ml)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const IJK_Field_double& vitesse_i_filtre = velocity_filtre[0];
  const IJK_Field_double& vitesse_j_filtre = velocity_filtre[1];
  const IJK_Field_double& vitesse_k_filtre = velocity_filtre[2];

  const IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  const IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  const IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  const IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  const IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  const IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const IJK_Field_double& structural_uu_filtre_xx = structural_uu_filtre_tensor[0];
  const IJK_Field_double& structural_uu_filtre_xy = structural_uu_filtre_tensor[1];
  const IJK_Field_double& structural_uu_filtre_xz = structural_uu_filtre_tensor[2];
  const IJK_Field_double& structural_uu_filtre_yy = structural_uu_filtre_tensor[3];
  const IJK_Field_double& structural_uu_filtre_yz = structural_uu_filtre_tensor[4];
  const IJK_Field_double& structural_uu_filtre_zz = structural_uu_filtre_tensor[5];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& bf_ii = tmp_b[0];
  IJK_Field_local_double& bf_ij = tmp_b[1];
  IJK_Field_local_double& bf_ik = tmp_b[2];
  IJK_Field_local_double& bf_jj = tmp_b[3];
  IJK_Field_local_double& bf_jk = tmp_b[4];
  IJK_Field_local_double& bf_kk = tmp_b[5];
  IJK_Field_local_double& buu_ii = tmp_b[6];
  IJK_Field_local_double& buu_ij = tmp_b[7];
  IJK_Field_local_double& buu_ik = tmp_b[8];
  IJK_Field_local_double& buu_jj = tmp_b[9];
  IJK_Field_local_double& buu_jk = tmp_b[10];
  IJK_Field_local_double& buu_kk = tmp_b[11];
  IJK_Field_local_double& bg_ii = tmp_b[12];
  IJK_Field_local_double& bg_ij = tmp_b[13];
  IJK_Field_local_double& bg_ik = tmp_b[14];
  IJK_Field_local_double& bg_jj = tmp_b[15];
  IJK_Field_local_double& bg_jk = tmp_b[16];
  IJK_Field_local_double& bg_kk = tmp_b[17];

  IJK_Field_local_double& af_ii = tmp_a[0];
  IJK_Field_local_double& af_ij = tmp_a[1];
  IJK_Field_local_double& af_ik = tmp_a[2];
  IJK_Field_local_double& af_jj = tmp_a[3];
  IJK_Field_local_double& af_jk = tmp_a[4];
  IJK_Field_local_double& af_kk = tmp_a[5];
  IJK_Field_local_double& auu_ii = tmp_a[6];
  IJK_Field_local_double& auu_ij = tmp_a[7];
  IJK_Field_local_double& auu_ik = tmp_a[8];
  IJK_Field_local_double& auu_jj = tmp_a[9];
  IJK_Field_local_double& auu_jk = tmp_a[10];
  IJK_Field_local_double& auu_kk = tmp_a[11];
  IJK_Field_local_double& ag_ii = tmp_a[12];
  IJK_Field_local_double& ag_ij = tmp_a[13];
  IJK_Field_local_double& ag_ik = tmp_a[14];
  IJK_Field_local_double& ag_jj = tmp_a[15];
  IJK_Field_local_double& ag_jk = tmp_a[16];
  IJK_Field_local_double& ag_kk = tmp_a[17];

  ArrOfDouble& moy_lij_xx = ml[0][0];
  ArrOfDouble& moy_lij_xy = ml[0][1];
  ArrOfDouble& moy_lij_xz = ml[0][2];
  ArrOfDouble& moy_lij_yy = ml[0][3];
  ArrOfDouble& moy_lij_yz = ml[0][4];
  ArrOfDouble& moy_lij_zz = ml[0][5];

  ArrOfDouble& moy_mij_xx = ml[1][0];
  ArrOfDouble& moy_mij_xy = ml[1][1];
  ArrOfDouble& moy_mij_xz = ml[1][2];
  ArrOfDouble& moy_mij_yy = ml[1][3];
  ArrOfDouble& moy_mij_yz = ml[1][4];
  ArrOfDouble& moy_mij_zz = ml[1][5];

  ArrOfDouble& moy_hij_xx = ml[2][0];
  ArrOfDouble& moy_hij_xy = ml[2][1];
  ArrOfDouble& moy_hij_xz = ml[2][2];
  ArrOfDouble& moy_hij_yy = ml[2][3];
  ArrOfDouble& moy_hij_yz = ml[2][4];
  ArrOfDouble& moy_hij_zz = ml[2][5];

  ArrOfDouble& moy_mijmij_xx = ml[3][0];
  ArrOfDouble& moy_mijmij_xy = ml[3][1];
  ArrOfDouble& moy_mijmij_xz = ml[3][2];
  ArrOfDouble& moy_mijmij_yy = ml[3][3];
  ArrOfDouble& moy_mijmij_yz = ml[3][4];
  ArrOfDouble& moy_mijmij_zz = ml[3][5];
  ArrOfDouble& moy_mijmij = ml[3][6];

  ArrOfDouble& moy_hijhij_xx = ml[4][0];
  ArrOfDouble& moy_hijhij_xy = ml[4][1];
  ArrOfDouble& moy_hijhij_xz = ml[4][2];
  ArrOfDouble& moy_hijhij_yy = ml[4][3];
  ArrOfDouble& moy_hijhij_yz = ml[4][4];
  ArrOfDouble& moy_hijhij_zz = ml[4][5];
  ArrOfDouble& moy_hijhij = ml[4][6];

  ArrOfDouble& moy_mijlij_xx = ml[5][0];
  ArrOfDouble& moy_mijlij_xy = ml[5][1];
  ArrOfDouble& moy_mijlij_xz = ml[5][2];
  ArrOfDouble& moy_mijlij_yy = ml[5][3];
  ArrOfDouble& moy_mijlij_yz = ml[5][4];
  ArrOfDouble& moy_mijlij_zz = ml[5][5];
  ArrOfDouble& moy_mijlij = ml[5][6];

  ArrOfDouble& moy_hijlij_xx = ml[6][0];
  ArrOfDouble& moy_hijlij_xy = ml[6][1];
  ArrOfDouble& moy_hijlij_xz = ml[6][2];
  ArrOfDouble& moy_hijlij_yy = ml[6][3];
  ArrOfDouble& moy_hijlij_yz = ml[6][4];
  ArrOfDouble& moy_hijlij_zz = ml[6][5];
  ArrOfDouble& moy_hijlij = ml[6][6];

  ArrOfDouble& moy_mijhij_xx = ml[7][0];
  ArrOfDouble& moy_mijhij_xy = ml[7][1];
  ArrOfDouble& moy_mijhij_xz = ml[7][2];
  ArrOfDouble& moy_mijhij_yy = ml[7][3];
  ArrOfDouble& moy_mijhij_yz = ml[7][4];
  ArrOfDouble& moy_mijhij_zz = ml[7][5];
  ArrOfDouble& moy_mijhij = ml[7][6];

  double m_ii = 0.;
  double m_ij = 0.;
  double m_ik = 0.;
  double m_jj = 0.;
  double m_jk = 0.;
  double m_kk = 0.;

  double h_ii = 0.;
  double h_ij = 0.;
  double h_ik = 0.;
  double h_jj = 0.;
  double h_jk = 0.;
  double h_kk = 0.;

  FixedVector<FixedVector<double, 3>, 3> f;
  FixedVector<FixedVector<double, 3>, 3> f_filtre;

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      if (tensorial)
        {
          moy_lij_xx[kg] = 0;
          moy_lij_xy[kg] = 0;
          moy_lij_xz[kg] = 0;
          moy_lij_yy[kg] = 0;
          moy_lij_yz[kg] = 0;
          moy_lij_zz[kg] = 0;

          moy_mij_xx[kg] = 0;
          moy_mij_xy[kg] = 0;
          moy_mij_xz[kg] = 0;
          moy_mij_yy[kg] = 0;
          moy_mij_yz[kg] = 0;
          moy_mij_zz[kg] = 0;

          moy_hij_xx[kg] = 0;
          moy_hij_xy[kg] = 0;
          moy_hij_xz[kg] = 0;
          moy_hij_yy[kg] = 0;
          moy_hij_yz[kg] = 0;
          moy_hij_zz[kg] = 0;

          moy_mijmij_xx[kg] = 0;
          moy_mijmij_xy[kg] = 0;
          moy_mijmij_xz[kg] = 0;
          moy_mijmij_yy[kg] = 0;
          moy_mijmij_yz[kg] = 0;
          moy_mijmij_zz[kg] = 0;

          moy_hijhij_xx[kg] = 0;
          moy_hijhij_xy[kg] = 0;
          moy_hijhij_xz[kg] = 0;
          moy_hijhij_yy[kg] = 0;
          moy_hijhij_yz[kg] = 0;
          moy_hijhij_zz[kg] = 0;

          moy_mijlij_xx[kg] = 0;
          moy_mijlij_xy[kg] = 0;
          moy_mijlij_xz[kg] = 0;
          moy_mijlij_yy[kg] = 0;
          moy_mijlij_yz[kg] = 0;
          moy_mijlij_zz[kg] = 0;

          moy_hijlij_xx[kg] = 0;
          moy_hijlij_xy[kg] = 0;
          moy_hijlij_xz[kg] = 0;
          moy_hijlij_yy[kg] = 0;
          moy_hijlij_yz[kg] = 0;
          moy_hijlij_zz[kg] = 0;

          moy_mijhij_xx[kg] = 0;
          moy_mijhij_xy[kg] = 0;
          moy_mijhij_xz[kg] = 0;
          moy_mijhij_yy[kg] = 0;
          moy_mijhij_yz[kg] = 0;
          moy_mijhij_zz[kg] = 0;
        }

      moy_mijmij[kg] = 0;
      moy_hijhij[kg] = 0;
      moy_mijlij[kg] = 0;
      moy_hijlij[kg] = 0;
      moy_mijhij[kg] = 0;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];

      const FixedVector<double, 21> filter_kernel_z = kernel->inhomogeneous(true, k, kg, nktot, dz_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_elem = kernel->size_k_elem(kg, nktot);
      const int shift_k_elem = kernel->shift_k_elem(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_elem = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_elem = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<-1 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double filter_coef_z = filter_kernel_z[kp+10];
                  longueur_elem += filter_coef_z * dz;
                }
              facteur_elem = 1./longueur_elem;
            }
          else
            {
              facteur_elem = dz_glo==0. ? 0. : 1./dz_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              buu_ii(i, j, 0) = 0.;
              buu_ij(i, j, 0) = 0.;
              buu_ik(i, j, 0) = 0.;
              buu_jj(i, j, 0) = 0.;
              buu_jk(i, j, 0) = 0.;
              buu_kk(i, j, 0) = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_elem(kpg, nktot)))
                    {
                      const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                      const double filter_coef_z = filter_kernel_z[kp+10];
                      const double uf_i = vitesse_i(i,j,k+kp);
                      const double vf_j = vitesse_j(i,j,k+kp);
                      const double wf_k = vitesse_k(i,j,k+kp);

                      const double uf_ip1 = vitesse_i(i+1,j,k+kp);
                      const double vf_jp1 = vitesse_j(i,j+1,k+kp);
                      const double wf_kp1 = kpg==(nktot-1) ? 0. : vitesse_k(i,j,k+1+kp);

                      const double ue = 0.5 * (uf_i + uf_ip1);
                      const double ve = 0.5 * (vf_j + vf_jp1);
                      const double we = 0.5 * (wf_k + wf_kp1);

                      const double uu_e = ue * ue;
                      const double uv_e = ue * ve;
                      const double uw_e = ue * we;
                      const double vv_e = ve * ve;
                      const double vw_e = ve * we;
                      const double ww_e = we * we;

                      if (ponderation_filter_kernel)
                        {
                          buu_ii(i, j, 0) += uu_e * filter_coef_z * dz * facteur_elem;
                          buu_ij(i, j, 0) += uv_e * filter_coef_z * dz * facteur_elem;
                          buu_ik(i, j, 0) += uw_e * filter_coef_z * dz * facteur_elem;
                          buu_jj(i, j, 0) += vv_e * filter_coef_z * dz * facteur_elem;
                          buu_jk(i, j, 0) += vw_e * filter_coef_z * dz * facteur_elem;
                          buu_kk(i, j, 0) += ww_e * filter_coef_z * dz * facteur_elem;
                        }
                      else
                        {
                          buu_ii(i, j, 0) += uu_e * filter_coef_z;
                          buu_ij(i, j, 0) += uv_e * filter_coef_z;
                          buu_ik(i, j, 0) += uw_e * filter_coef_z;
                          buu_jj(i, j, 0) += vv_e * filter_coef_z;
                          buu_jk(i, j, 0) += vw_e * filter_coef_z;
                          buu_kk(i, j, 0) += ww_e * filter_coef_z;
                        }
                    }
                }

              if (turbulent_viscosity)
                {
                  bf_ii(i, j, 0) = 0.;
                  bf_ij(i, j, 0) = 0.;
                  bf_ik(i, j, 0) = 0.;
                  bf_jj(i, j, 0) = 0.;
                  bf_jk(i, j, 0) = 0.;
                  bf_kk(i, j, 0) = 0.;
                  for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                    {
                      const int kpg = kg + kp;
                      if (!(kernel->is_at_wall_elem(kpg, nktot)))
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double filter_coef_z = filter_kernel_z[kp+10];
                          if (!anisotropic)
                            {
                              calculer_tau(i, j, k+kp, velocity, turbulent_mu_xx, turbulent_mu_xy, turbulent_mu_xz, turbulent_mu_yy, turbulent_mu_yz, turbulent_mu_zz, delta_z, 1., 1., 1., 1., 1., f);
                            }
                          else
                            {
                              const double dz_pour_delta = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z_pour_delta[k+kp];
                              const double dz_m1_pour_delta = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z_pour_delta[k-1+kp];
                              const double dz_p1_pour_delta = (kpg+1<0 || kpg+1>(nktot-1)) ? 0. : delta_z_pour_delta[k+1+kp];
                              const double delta_m_pour_delta = (kpg-1<0 || kpg>(nktot-1)) ? 0. : 0.5*(dz_pour_delta + dz_m1_pour_delta);
                              const double delta_p_pour_delta = (kpg<0 || kpg+1>(nktot-1)) ? 0. : 0.5*(dz_pour_delta + dz_p1_pour_delta);
                              calculer_tau(i, j, k+kp, velocity, turbulent_mu_xx, turbulent_mu_xy, turbulent_mu_xz, turbulent_mu_yy, turbulent_mu_yz, turbulent_mu_zz, delta_z, dx_pour_delta, dy_pour_delta, dz_pour_delta, delta_m_pour_delta, delta_p_pour_delta, f);
                            }

                          if (ponderation_filter_kernel)
                            {
                              bf_ii(i, j, 0) += f[0][0] * filter_coef_z * dz * facteur_elem;
                              bf_ij(i, j, 0) += f[0][1] * filter_coef_z * dz * facteur_elem;
                              bf_ik(i, j, 0) += f[0][2] * filter_coef_z * dz * facteur_elem;
                              bf_jj(i, j, 0) += f[1][1] * filter_coef_z * dz * facteur_elem;
                              bf_jk(i, j, 0) += f[1][2] * filter_coef_z * dz * facteur_elem;
                              bf_kk(i, j, 0) += f[2][2] * filter_coef_z * dz * facteur_elem;
                            }
                          else
                            {
                              bf_ii(i, j, 0) += f[0][0] * filter_coef_z;
                              bf_ij(i, j, 0) += f[0][1] * filter_coef_z;
                              bf_ik(i, j, 0) += f[0][2] * filter_coef_z;
                              bf_jj(i, j, 0) += f[1][1] * filter_coef_z;
                              bf_jk(i, j, 0) += f[1][2] * filter_coef_z;
                              bf_kk(i, j, 0) += f[2][2] * filter_coef_z;
                            }
                        }
                    }
                }

              if (structural_uu)
                {
                  bg_ii(i, j, 0) = 0.;
                  bg_ij(i, j, 0) = 0.;
                  bg_ik(i, j, 0) = 0.;
                  bg_jj(i, j, 0) = 0.;
                  bg_jk(i, j, 0) = 0.;
                  bg_kk(i, j, 0) = 0.;
                  for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                    {
                      const int kpg = kg + kp;
                      if (!(kernel->is_at_wall_elem(kpg, nktot)))
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double filter_coef_z = filter_kernel_z[kp+10];
                          const double gxx_e = -structural_uu_xx(i,j,k+kp);
                          const double gyy_e = -structural_uu_yy(i,j,k+kp);
                          const double gzz_e = -structural_uu_zz(i,j,k+kp);

                          const double gxy_ij = -structural_uu_xy(i,j,k+kp);
                          const double gxz_ik = -structural_uu_xz(i,j,k+kp);
                          const double gyz_jk = -structural_uu_yz(i,j,k+kp);

                          const double gxy_ip1j = -structural_uu_xy(i+1,j,k+kp);
                          const double gxz_ip1k = -structural_uu_xz(i+1,j,k+kp);
                          const double gyz_jp1k = -structural_uu_yz(i,j+1,k+kp);

                          const double gxy_ijp1 = -structural_uu_xy(i,j+1,k+kp);
                          const double gxz_ikp1 = kpg==(nktot-1) ? 0. : -structural_uu_xz(i,j,k+1+kp);
                          const double gyz_jkp1 = kpg==(nktot-1) ? 0. : -structural_uu_yz(i,j,k+1+kp);

                          const double gxy_ip1jp1 = -structural_uu_xy(i+1,j+1,k+kp);
                          const double gxz_ip1kp1 = kpg==(nktot-1) ? 0. : -structural_uu_xz(i+1,j,k+1+kp);
                          const double gyz_jp1kp1 = kpg==(nktot-1) ? 0. : -structural_uu_yz(i,j+1,k+1+kp);

                          const double gxy_e = 0.25 * (gxy_ij + gxy_ip1j + gxy_ijp1 + gxy_ip1jp1);
                          const double gxz_e = 0.25 * (gxz_ik + gxz_ip1k + gxz_ikp1 + gxz_ip1kp1);
                          const double gyz_e = 0.25 * (gyz_jk + gyz_jp1k + gyz_jkp1 + gyz_jp1kp1);

                          if (ponderation_filter_kernel)
                            {
                              bg_ii(i, j, 0) += gxx_e * filter_coef_z * dz * facteur_elem;
                              bg_ij(i, j, 0) += gxy_e * filter_coef_z * dz * facteur_elem;
                              bg_ik(i, j, 0) += gxz_e * filter_coef_z * dz * facteur_elem;
                              bg_jj(i, j, 0) += gyy_e * filter_coef_z * dz * facteur_elem;
                              bg_jk(i, j, 0) += gyz_e * filter_coef_z * dz * facteur_elem;
                              bg_kk(i, j, 0) += gzz_e * filter_coef_z * dz * facteur_elem;
                            }
                          else
                            {
                              bg_ii(i, j, 0) += gxx_e * filter_coef_z;
                              bg_ij(i, j, 0) += gxy_e * filter_coef_z;
                              bg_ik(i, j, 0) += gxz_e * filter_coef_z;
                              bg_jj(i, j, 0) += gyy_e * filter_coef_z;
                              bg_jk(i, j, 0) += gyz_e * filter_coef_z;
                              bg_kk(i, j, 0) += gzz_e * filter_coef_z;
                            }
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              auu_ii(i, 0, 0) = 0.;
              auu_ij(i, 0, 0) = 0.;
              auu_ik(i, 0, 0) = 0.;
              auu_jj(i, 0, 0) = 0.;
              auu_jk(i, 0, 0) = 0.;
              auu_kk(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  auu_ii(i, 0, 0) += buu_ii(i, j+jp, 0) * filter_coef_y;
                  auu_ij(i, 0, 0) += buu_ij(i, j+jp, 0) * filter_coef_y;
                  auu_ik(i, 0, 0) += buu_ik(i, j+jp, 0) * filter_coef_y;
                  auu_jj(i, 0, 0) += buu_jj(i, j+jp, 0) * filter_coef_y;
                  auu_jk(i, 0, 0) += buu_jk(i, j+jp, 0) * filter_coef_y;
                  auu_kk(i, 0, 0) += buu_kk(i, j+jp, 0) * filter_coef_y;
                }
              if (turbulent_viscosity)
                {
                  af_ii(i, 0, 0) = 0.;
                  af_ij(i, 0, 0) = 0.;
                  af_ik(i, 0, 0) = 0.;
                  af_jj(i, 0, 0) = 0.;
                  af_jk(i, 0, 0) = 0.;
                  af_kk(i, 0, 0) = 0.;
                  for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                    {
                      const double filter_coef_y = filter_kernel_y[jp+10];
                      af_ii(i, 0, 0) += bf_ii(i, j+jp, 0) * filter_coef_y;
                      af_ij(i, 0, 0) += bf_ij(i, j+jp, 0) * filter_coef_y;
                      af_ik(i, 0, 0) += bf_ik(i, j+jp, 0) * filter_coef_y;
                      af_jj(i, 0, 0) += bf_jj(i, j+jp, 0) * filter_coef_y;
                      af_jk(i, 0, 0) += bf_jk(i, j+jp, 0) * filter_coef_y;
                      af_kk(i, 0, 0) += bf_kk(i, j+jp, 0) * filter_coef_y;
                    }
                }
              if (structural_uu)
                {
                  ag_ii(i, 0, 0) = 0.;
                  ag_ij(i, 0, 0) = 0.;
                  ag_ik(i, 0, 0) = 0.;
                  ag_jj(i, 0, 0) = 0.;
                  ag_jk(i, 0, 0) = 0.;
                  ag_kk(i, 0, 0) = 0.;
                  for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                    {
                      const double filter_coef_y = filter_kernel_y[jp+10];
                      ag_ii(i, 0, 0) += bg_ii(i, j+jp, 0) * filter_coef_y;
                      ag_ij(i, 0, 0) += bg_ij(i, j+jp, 0) * filter_coef_y;
                      ag_ik(i, 0, 0) += bg_ik(i, j+jp, 0) * filter_coef_y;
                      ag_jj(i, 0, 0) += bg_jj(i, j+jp, 0) * filter_coef_y;
                      ag_jk(i, 0, 0) += bg_jk(i, j+jp, 0) * filter_coef_y;
                      ag_kk(i, 0, 0) += bg_kk(i, j+jp, 0) * filter_coef_y;
                    }
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double ruu_ii = 0.;
              double ruu_ij = 0.;
              double ruu_ik = 0.;
              double ruu_jj = 0.;
              double ruu_jk = 0.;
              double ruu_kk = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  ruu_ii += auu_ii(i+ip, 0, 0) * filter_coef_x;
                  ruu_ij += auu_ij(i+ip, 0, 0) * filter_coef_x;
                  ruu_ik += auu_ik(i+ip, 0, 0) * filter_coef_x;
                  ruu_jj += auu_jj(i+ip, 0, 0) * filter_coef_x;
                  ruu_jk += auu_jk(i+ip, 0, 0) * filter_coef_x;
                  ruu_kk += auu_kk(i+ip, 0, 0) * filter_coef_x;
                }

              const double uf_i = vitesse_i_filtre(i,j,k);
              const double vf_j = vitesse_j_filtre(i,j,k);
              const double wf_k = vitesse_k_filtre(i,j,k);

              const double uf_ip1 = vitesse_i_filtre(i+1,j,k);
              const double vf_jp1 = vitesse_j_filtre(i,j+1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k_filtre(i,j,k+1);

              const double ue = 0.5 * (uf_i + uf_ip1);
              const double ve = 0.5 * (vf_j + vf_jp1);
              const double we = 0.5 * (wf_k + wf_kp1);

              const double l_ii = ruu_ii - ue*ue;
              const double l_ij = ruu_ij - ue*ve;
              const double l_ik = ruu_ik - ue*we;
              const double l_jj = ruu_jj - ve*ve;
              const double l_jk = ruu_jk - ve*we;
              const double l_kk = ruu_kk - we*we;

              if (tensorial)
                {
                  moy_lij_xx[kg] += l_ii;
                  moy_lij_xy[kg] += l_ij;
                  moy_lij_xz[kg] += l_ik;
                  moy_lij_yy[kg] += l_jj;
                  moy_lij_yz[kg] += l_jk;
                  moy_lij_zz[kg] += l_kk;
                }

              if (turbulent_viscosity)
                {
                  double rf_ii = 0.;
                  double rf_ij = 0.;
                  double rf_ik = 0.;
                  double rf_jj = 0.;
                  double rf_jk = 0.;
                  double rf_kk = 0.;
                  for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                    {
                      const double filter_coef_x = filter_kernel_x[ip+10];
                      rf_ii += af_ii(i+ip, 0, 0) * filter_coef_x;
                      rf_ij += af_ij(i+ip, 0, 0) * filter_coef_x;
                      rf_ik += af_ik(i+ip, 0, 0) * filter_coef_x;
                      rf_jj += af_jj(i+ip, 0, 0) * filter_coef_x;
                      rf_jk += af_jk(i+ip, 0, 0) * filter_coef_x;
                      rf_kk += af_kk(i+ip, 0, 0) * filter_coef_x;
                    }

                  if (!anisotropic)
                    {
                      calculer_tau(i, j, k, velocity_filtre, turbulent_mu_filtre_xx, turbulent_mu_filtre_xy, turbulent_mu_filtre_xz, turbulent_mu_filtre_yy, turbulent_mu_filtre_yz, turbulent_mu_filtre_zz, delta_z, 1., 1., 1., 1., 1., f_filtre);
                    }
                  else
                    {
                      const double dx_filtre = facteur_delta_filtre_x*dx;
                      const double dy_filtre = facteur_delta_filtre_y*dy;
                      const double dz_filtre = delta_z_pour_delta_filtre[k];
                      const double delta_m_filtre = kg==0 ? 0. : 0.5*(dz_filtre + delta_z_pour_delta_filtre[k-1]);
                      const double delta_p_filtre = kg==(nktot-1) ? 0. : 0.5*(dz_filtre + delta_z_pour_delta_filtre[k+1]);
                      calculer_tau(i, j, k, velocity_filtre, turbulent_mu_filtre_xx, turbulent_mu_filtre_xy, turbulent_mu_filtre_xz, turbulent_mu_filtre_yy, turbulent_mu_filtre_yz, turbulent_mu_filtre_zz, delta_z, dx_filtre, dy_filtre, dz_filtre, delta_m_filtre, delta_p_filtre, f_filtre);
                    }

                  m_ii = f_filtre[0][0] - rf_ii;
                  m_ij = f_filtre[0][1] - rf_ij;
                  m_ik = f_filtre[0][2] - rf_ik;
                  m_jj = f_filtre[1][1] - rf_jj;
                  m_jk = f_filtre[1][2] - rf_jk;
                  m_kk = f_filtre[2][2] - rf_kk;

                  const double  mijmij = m_ii * m_ii
                                         + m_ij * m_ij
                                         + m_ik * m_ik
                                         + m_ij * m_ij
                                         + m_jj * m_jj
                                         + m_jk * m_jk
                                         + m_ik * m_ik
                                         + m_jk * m_jk
                                         + m_kk * m_kk;

                  const double mijlij = m_ii * l_ii
                                        + m_ij * l_ij
                                        + m_ik * l_ik
                                        + m_ij * l_ij
                                        + m_jj * l_jj
                                        + m_jk * l_jk
                                        + m_ik * l_ik
                                        + m_jk * l_jk
                                        + m_kk * l_kk;

                  moy_mijmij[kg] += mijmij;
                  moy_mijlij[kg] += mijlij;

                  if (tensorial)
                    {
                      moy_mij_xx[kg] += m_ii;
                      moy_mij_xy[kg] += m_ij;
                      moy_mij_xz[kg] += m_ik;
                      moy_mij_yy[kg] += m_jj;
                      moy_mij_yz[kg] += m_jk;
                      moy_mij_zz[kg] += m_kk;

                      moy_mijmij_xx[kg] += m_ii * m_ii;
                      moy_mijmij_xy[kg] += m_ij * m_ij;
                      moy_mijmij_xz[kg] += m_ik * m_ik;
                      moy_mijmij_yy[kg] += m_jj * m_jj;
                      moy_mijmij_yz[kg] += m_jk * m_jk;
                      moy_mijmij_zz[kg] += m_kk * m_kk;

                      moy_mijlij_xx[kg] += m_ii * l_ii;
                      moy_mijlij_xy[kg] += m_ij * l_ij;
                      moy_mijlij_xz[kg] += m_ik * l_ik;
                      moy_mijlij_yy[kg] += m_jj * l_jj;
                      moy_mijlij_yz[kg] += m_jk * l_jk;
                      moy_mijlij_zz[kg] += m_kk * l_kk;
                    }
                }

              if (structural_uu)
                {
                  double rg_ii = 0.;
                  double rg_ij = 0.;
                  double rg_ik = 0.;
                  double rg_jj = 0.;
                  double rg_jk = 0.;
                  double rg_kk = 0.;
                  for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                    {
                      const double filter_coef_x = filter_kernel_x[ip+10];
                      rg_ii += ag_ii(i+ip, 0, 0) * filter_coef_x;
                      rg_ij += ag_ij(i+ip, 0, 0) * filter_coef_x;
                      rg_ik += ag_ik(i+ip, 0, 0) * filter_coef_x;
                      rg_jj += ag_jj(i+ip, 0, 0) * filter_coef_x;
                      rg_jk += ag_jk(i+ip, 0, 0) * filter_coef_x;
                      rg_kk += ag_kk(i+ip, 0, 0) * filter_coef_x;
                    }

                  const double gxx_e = -structural_uu_filtre_xx(i,j,k);
                  const double gyy_e = -structural_uu_filtre_yy(i,j,k);
                  const double gzz_e = -structural_uu_filtre_zz(i,j,k);

                  const double gxy_ij = -structural_uu_filtre_xy(i,j,k);
                  const double gxz_ik = -structural_uu_filtre_xz(i,j,k);
                  const double gyz_jk = -structural_uu_filtre_yz(i,j,k);

                  const double gxy_ip1j = -structural_uu_filtre_xy(i+1,j,k);
                  const double gxz_ip1k = -structural_uu_filtre_xz(i+1,j,k);
                  const double gyz_jp1k = -structural_uu_filtre_yz(i,j+1,k);

                  const double gxy_ijp1 = -structural_uu_filtre_xy(i,j+1,k);
                  const double gxz_ikp1 = kg==(nktot-1) ? 0. : -structural_uu_filtre_xz(i,j,k+1);
                  const double gyz_jkp1 = kg==(nktot-1) ? 0. : -structural_uu_filtre_yz(i,j,k+1);

                  const double gxy_ip1jp1 = -structural_uu_filtre_xy(i+1,j+1,k);
                  const double gxz_ip1kp1 = kg==(nktot-1) ? 0. : -structural_uu_filtre_xz(i+1,j,k+1);
                  const double gyz_jp1kp1 = kg==(nktot-1) ? 0. : -structural_uu_filtre_yz(i,j+1,k+1);

                  const double gxy_e = 0.25 * (gxy_ij + gxy_ip1j + gxy_ijp1 + gxy_ip1jp1);
                  const double gxz_e = 0.25 * (gxz_ik + gxz_ip1k + gxz_ikp1 + gxz_ip1kp1);
                  const double gyz_e = 0.25 * (gyz_jk + gyz_jp1k + gyz_jkp1 + gyz_jp1kp1);

                  h_ii = gxx_e - rg_ii;
                  h_ij = gxy_e - rg_ij;
                  h_ik = gxz_e - rg_ik;
                  h_jj = gyy_e - rg_jj;
                  h_jk = gyz_e - rg_jk;
                  h_kk = gzz_e - rg_kk;

                  const double hijhij = h_ii * h_ii
                                        + h_ij * h_ij
                                        + h_ik * h_ik
                                        + h_ij * h_ij
                                        + h_jj * h_jj
                                        + h_jk * h_jk
                                        + h_ik * h_ik
                                        + h_jk * h_jk
                                        + h_kk * h_kk;

                  const double hijlij = h_ii * l_ii
                                        + h_ij * l_ij
                                        + h_ik * l_ik
                                        + h_ij * l_ij
                                        + h_jj * l_jj
                                        + h_jk * l_jk
                                        + h_ik * l_ik
                                        + h_jk * l_jk
                                        + h_kk * l_kk;

                  moy_hijhij[kg] += hijhij;
                  moy_hijlij[kg] += hijlij;

                  if (tensorial)
                    {
                      moy_hij_xx[kg] += h_ii;
                      moy_hij_xy[kg] += h_ij;
                      moy_hij_xz[kg] += h_ik;
                      moy_hij_yy[kg] += h_jj;
                      moy_hij_yz[kg] += h_jk;
                      moy_hij_zz[kg] += h_kk;

                      moy_hijhij_xx[kg] += h_ii * h_ii;
                      moy_hijhij_xy[kg] += h_ij * h_ij;
                      moy_hijhij_xz[kg] += h_ik * h_ik;
                      moy_hijhij_yy[kg] += h_jj * h_jj;
                      moy_hijhij_yz[kg] += h_jk * h_jk;
                      moy_hijhij_zz[kg] += h_kk * h_kk;

                      moy_hijlij_xx[kg] += h_ii * l_ii;
                      moy_hijlij_xy[kg] += h_ij * l_ij;
                      moy_hijlij_xz[kg] += h_ik * l_ik;
                      moy_hijlij_yy[kg] += h_jj * l_jj;
                      moy_hijlij_yz[kg] += h_jk * l_jk;
                      moy_hijlij_zz[kg] += h_kk * l_kk;
                    }

                  if (turbulent_viscosity)
                    {
                      const double mijhij = m_ii * h_ii
                                            + m_ij * h_ij
                                            + m_ik * h_ik
                                            + m_ij * h_ij
                                            + m_jj * h_jj
                                            + m_jk * h_jk
                                            + m_ik * h_ik
                                            + m_jk * h_jk
                                            + m_kk * h_kk;

                      moy_mijhij[kg] += mijhij;

                      if (tensorial)
                        {
                          moy_mijhij_xx[kg] += m_ii * h_ii;
                          moy_mijhij_xy[kg] += m_ij * h_ij;
                          moy_mijhij_xz[kg] += m_ik * h_ik;
                          moy_mijhij_yy[kg] += m_jj * h_jj;
                          moy_mijhij_yz[kg] += m_jk * h_jk;
                          moy_mijhij_zz[kg] += m_kk * h_kk;
                        }
                    }
                }
            }
        }
    }
}

void calculer_constante_direct(const bool global,
                               const bool clipping,
                               const Nom& description,
                               const ArrOfDouble& moy_mij_0,
                               const ArrOfDouble& moy_lij_0,
                               const bool flag_mixte,
                               const ArrOfDouble& moy_hij_0,
                               const FixedVector<IJK_Field_double, 3>& velocity,
                               const ArrOfDouble_with_ghost& delta_z,
                               ArrOfDouble_with_ghost& constante_modele)
{
  const IJK_Field_double& vitesse_k = velocity[2];

  // Copie des tableaux
  ArrOfDouble moy_mij = moy_mij_0;
  ArrOfDouble moy_lij = moy_lij_0;
  ArrOfDouble moy_hij = moy_hij_0;

  const int nk = vitesse_k.nk();
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  Cout <<
       "Dynamic correction of constant [ " << description <<
       " ] type= direct global= " << (int)global <<
       " clipping= " << (int)clipping <<
       " mixte= " << (int)flag_mixte << " :"
       << finl;

  if (global)
    {
      double moy_mij_global = 0;
      double moy_lij_global = 0;
      double moy_hij_global = 0;
      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          moy_mij_global += moy_mij[kg]*delta_z[k];
          moy_lij_global += moy_lij[kg]*delta_z[k];
          if (flag_mixte)
            {
              moy_hij_global += moy_hij[kg]*delta_z[k];
            }
        }

      ArrOfDouble tmp(3);
      tmp[0] = moy_mij_global;
      tmp[1] = moy_lij_global;
      tmp[2] = moy_hij_global;
      mp_sum_for_each_item(tmp);

      double c_global = (tmp[0]==0.) ? 1. : (tmp[1] - tmp[2])/tmp[0];
      if (clipping)
        {
          c_global = max(0., c_global);
        }
      Cout << "(global) constante_modele= " << c_global << finl;

      for (int k = 0; k < nk; k++)
        {
          constante_modele[k] = c_global;
        }
    }
  else
    {
      mp_sum_for_each_item(moy_mij);
      mp_sum_for_each_item(moy_lij);
      if (flag_mixte)
        {
          mp_sum_for_each_item(moy_hij);
        }

      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          if (flag_mixte)
            {
              constante_modele[k] = (moy_mij[kg]==0.) ? 1. : (moy_lij[kg] - moy_hij[kg])/moy_mij[kg];
            }
          else
            {
              constante_modele[k] = (moy_mij[kg]==0.) ? 1. : moy_lij[kg]/moy_mij[kg];
            }
          if (clipping)
            {
              constante_modele[k] = max(0., constante_modele[k]);
            }
          Cout << "k= " << kg << " constante_modele= " << constante_modele[k] << finl;
        }
    }
}

void calculer_constante_lilly(const bool global,
                              const bool clipping,
                              const Nom& description,
                              const ArrOfDouble& moy_mijmij_0,
                              const ArrOfDouble& moy_mijlij_0,
                              const bool flag_mixte,
                              const ArrOfDouble& moy_mijhij_0,
                              const FixedVector<IJK_Field_double, 3>& velocity,
                              const ArrOfDouble_with_ghost& delta_z,
                              ArrOfDouble_with_ghost& constante_modele)
{
  const IJK_Field_double& vitesse_k = velocity[2];

  // Copie des tableaux
  ArrOfDouble moy_mijmij = moy_mijmij_0;
  ArrOfDouble moy_mijlij = moy_mijlij_0;
  ArrOfDouble moy_mijhij = moy_mijhij_0;

  const int nk = vitesse_k.nk();
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  Cout << "Dynamic correction of constant [ " << description << " ] type= lilly global= " << (int)global << " clipping= " << (int)clipping << " mixte= " << (int)flag_mixte << " :" << finl;

  if (global)
    {
      double moy_mijmij_global = 0;
      double moy_mijlij_global = 0;
      double moy_mijhij_global = 0;
      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          moy_mijmij_global += moy_mijmij[kg]*delta_z[k];
          moy_mijlij_global += moy_mijlij[kg]*delta_z[k];
          if (flag_mixte)
            {
              moy_mijhij_global += moy_mijhij[kg]*delta_z[k];
            }
        }

      ArrOfDouble tmp(3);
      tmp[0] = moy_mijmij_global;
      tmp[1] = moy_mijlij_global;
      tmp[2] = moy_mijhij_global;
      mp_sum_for_each_item(tmp);

      double c_global = (tmp[0]==0.) ? 1. : (tmp[1] - tmp[2])/tmp[0];
      if (clipping)
        {
          c_global = max(0., c_global);
        }
      Cout << "(global) constante_modele= " << c_global << finl;

      for (int k = 0; k < nk; k++)
        {
          constante_modele[k] = c_global;
        }
    }
  else
    {
      mp_sum_for_each_item(moy_mijmij);
      mp_sum_for_each_item(moy_mijlij);
      if (flag_mixte)
        {
          mp_sum_for_each_item(moy_mijhij);
        }

      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          if (flag_mixte)
            {
              constante_modele[k] = (moy_mijmij[kg]==0.) ? 1. : (moy_mijlij[kg] - moy_mijhij[kg])/moy_mijmij[kg];
            }
          else
            {
              constante_modele[k] = (moy_mijmij[kg]==0.) ? 1. : moy_mijlij[kg]/moy_mijmij[kg];
            }
          if (clipping)
            {
              constante_modele[k] = max(0., constante_modele[k]);
            }
          Cout << "k= " << kg << " constante_modele= " << constante_modele[k] << finl;
        }
    }
}

void calculer_constante_twoparameters(const bool global,
                                      const bool clipping,
                                      const Nom& description,
                                      const ArrOfDouble& moy_mijmij_0,
                                      const ArrOfDouble& moy_hijhij_0,
                                      const ArrOfDouble& moy_mijlij_0,
                                      const ArrOfDouble& moy_hijlij_0,
                                      const ArrOfDouble& moy_mijhij_0,
                                      const FixedVector<IJK_Field_double, 3>& velocity,
                                      const ArrOfDouble_with_ghost& delta_z,
                                      ArrOfDouble_with_ghost& constante_modele)
{
  const IJK_Field_double& vitesse_k = velocity[2];

  // Copie des tableaux
  ArrOfDouble moy_mijmij = moy_mijmij_0;
  ArrOfDouble moy_hijhij = moy_hijhij_0;
  ArrOfDouble moy_mijlij = moy_mijlij_0;
  ArrOfDouble moy_hijlij = moy_hijlij_0;
  ArrOfDouble moy_mijhij = moy_mijhij_0;

  const int nk = vitesse_k.nk();
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  Cout << "Dynamic correction of constant [ " << description << " ] type= twoparameters global= " << (int)global << " clipping= " << (int)clipping << " :" << finl;

  if (global)
    {
      double moy_mijmij_global = 0;
      double moy_hijhij_global = 0;
      double moy_mijlij_global = 0;
      double moy_hijlij_global = 0;
      double moy_mijhij_global = 0;
      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          moy_mijmij_global += moy_mijmij[kg]*delta_z[k];
          moy_hijhij_global += moy_hijhij[kg]*delta_z[k];
          moy_mijlij_global += moy_mijlij[kg]*delta_z[k];
          moy_hijlij_global += moy_hijlij[kg]*delta_z[k];
          moy_mijhij_global += moy_mijhij[kg]*delta_z[k];
        }

      ArrOfDouble tmp(5);
      tmp[0] = moy_mijmij_global;
      tmp[1] = moy_hijhij_global;
      tmp[2] = moy_mijlij_global;
      tmp[3] = moy_hijlij_global;
      tmp[4] = moy_mijhij_global;
      mp_sum_for_each_item(tmp);

      double c_global = ((tmp[0]*tmp[1] - tmp[4]*tmp[4])==0.) ? 1. :
                        (tmp[2]*tmp[1] - tmp[3]*tmp[4])/(tmp[0]*tmp[1] - tmp[4]*tmp[4]);
      if (clipping)
        {
          c_global = max(0., c_global);
        }
      Cout << "(global) constante_modele= " << c_global << finl;

      for (int k = 0; k < nk; k++)
        {
          constante_modele[k] = c_global;
        }
    }
  else
    {
      mp_sum_for_each_item(moy_mijmij);
      mp_sum_for_each_item(moy_hijhij);
      mp_sum_for_each_item(moy_mijlij);
      mp_sum_for_each_item(moy_hijlij);
      mp_sum_for_each_item(moy_mijhij);

      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          constante_modele[k] = ((moy_mijmij[kg]*moy_hijhij[kg] - moy_mijhij[kg]*moy_mijhij[kg])==0.) ? 1. :
                                (moy_mijlij[kg]*moy_hijhij[kg] - moy_hijlij[kg]*moy_mijhij[kg])/(moy_mijmij[kg]*moy_hijhij[kg] - moy_mijhij[kg]*moy_mijhij[kg]);
          if (clipping)
            {
              constante_modele[k] = max(0., constante_modele[k]);
            }
          Cout << "k= " << kg << " constante_modele= " << constante_modele[k] << finl;
        }
    }
}

void calculer_constante_twonoerror(const bool global,
                                   const bool clipping,
                                   const Nom& description,
                                   const ArrOfDouble& moy_lij_0,
                                   const ArrOfDouble& moy_mij_0,
                                   const ArrOfDouble& moy_hij_0,
                                   const ArrOfDouble& moy_mijmij_0,
                                   const ArrOfDouble& moy_hijhij_0,
                                   const ArrOfDouble& moy_mijlij_0,
                                   const ArrOfDouble& moy_hijlij_0,
                                   const ArrOfDouble& moy_mijhij_0,
                                   const FixedVector<IJK_Field_double, 3>& velocity,
                                   const ArrOfDouble_with_ghost& delta_z,
                                   ArrOfDouble_with_ghost& constante_modele)
{
  const IJK_Field_double& vitesse_k = velocity[2];

  // Copie des tableaux
  ArrOfDouble moy_lij = moy_lij_0;
  ArrOfDouble moy_mij = moy_mij_0;
  ArrOfDouble moy_hij = moy_hij_0;
  ArrOfDouble moy_mijmij = moy_mijmij_0;
  ArrOfDouble moy_hijhij = moy_hijhij_0;
  ArrOfDouble moy_mijlij = moy_mijlij_0;
  ArrOfDouble moy_hijlij = moy_hijlij_0;
  ArrOfDouble moy_mijhij = moy_mijhij_0;

  const int nk = vitesse_k.nk();
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  Cout << "Dynamic correction of constant [ " << description << " ] type= twonoerror global= " << (int)global << " clipping= " << (int)clipping << " :" << finl;

  if (global)
    {
      double moy_lij_global = 0;
      double moy_mij_global = 0;
      double moy_hij_global = 0;
      double moy_mijmij_global = 0;
      double moy_hijhij_global = 0;
      double moy_mijlij_global = 0;
      double moy_hijlij_global = 0;
      double moy_mijhij_global = 0;
      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          moy_lij_global += moy_lij[kg]*delta_z[k];
          moy_mij_global += moy_mij[kg]*delta_z[k];
          moy_hij_global += moy_hij[kg]*delta_z[k];
          moy_mijmij_global += moy_mijmij[kg]*delta_z[k];
          moy_hijhij_global += moy_hijhij[kg]*delta_z[k];
          moy_mijlij_global += moy_mijlij[kg]*delta_z[k];
          moy_hijlij_global += moy_hijlij[kg]*delta_z[k];
          moy_mijhij_global += moy_mijhij[kg]*delta_z[k];
        }

      ArrOfDouble tmp(8);
      tmp[0] = moy_lij_global;
      tmp[1] = moy_mij_global;
      tmp[2] = moy_hij_global;
      tmp[3] = moy_mijmij_global;
      tmp[4] = moy_hijhij_global;
      tmp[5] = moy_mijlij_global;
      tmp[6] = moy_hijlij_global;
      tmp[7] = moy_mijhij_global;
      mp_sum_for_each_item(tmp);

      double c_global = ((tmp[2]*tmp[3] - tmp[1]*tmp[7] + tmp[2]*tmp[7] - tmp[1]*tmp[4])==0.) ? 1. :
                        (tmp[3]*tmp[0] - tmp[1]*tmp[5] + tmp[7]*tmp[0] - tmp[1]*tmp[6])/(tmp[2]*tmp[3] - tmp[1]*tmp[7] + tmp[2]*tmp[7] - tmp[1]*tmp[4]);
      if (clipping)
        {
          c_global = max(0., c_global);
        }
      Cout << "(global) constante_modele= " << c_global << finl;

      for (int k = 0; k < nk; k++)
        {
          constante_modele[k] = c_global;
        }
    }
  else
    {
      mp_sum_for_each_item(moy_lij);
      mp_sum_for_each_item(moy_mij);
      mp_sum_for_each_item(moy_hij);
      mp_sum_for_each_item(moy_mijmij);
      mp_sum_for_each_item(moy_hijhij);
      mp_sum_for_each_item(moy_mijlij);
      mp_sum_for_each_item(moy_hijlij);
      mp_sum_for_each_item(moy_mijhij);

      for (int k = 0; k < nk; k++)
        {
          const int kg = k + offset;

          constante_modele[k] = ((moy_hij[kg]*moy_mijmij[kg] - moy_mij[kg]*moy_mijhij[kg] + moy_hij[kg]*moy_mijhij[kg] - moy_mij[kg]*moy_hijhij[kg])==0.) ? 1. :
                                (moy_mijmij[kg]*moy_lij[kg] - moy_mij[kg]*moy_mijlij[kg] + moy_mijhij[kg]*moy_lij[kg] - moy_mij[kg]*moy_hijlij[kg])/(moy_hij[kg]*moy_mijmij[kg] - moy_mij[kg]*moy_mijhij[kg] + moy_hij[kg]*moy_mijhij[kg] - moy_mij[kg]*moy_hijhij[kg]);
          if (clipping)
            {
              constante_modele[k] = max(0., constante_modele[k]);
            }
          Cout << "k= " << kg << " constante_modele= " << constante_modele[k] << finl;
        }
    }
}

void multiplier_par_constante(const bool face,
                              ArrOfDouble_with_ghost& constante_modele,
                              const FixedVector<IJK_Field_double, 3>& velocity,
                              IJK_Field_double& turbulent_mu)
{
  const IJK_Field_double& vitesse_k = velocity[2];
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  if (face)
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  const int kg = k + offset;
                  const double constante_face = kg==0 ? 0. : 0.5*(constante_modele[k] + constante_modele[k-1]);
                  turbulent_mu(i,j,k) = constante_face * turbulent_mu(i,j,k);
                }
            }
        }
    }
  else
    {
      for (int k = 0; k < nk; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  turbulent_mu(i,j,k) = constante_modele[k] * turbulent_mu(i,j,k);
                }
            }
        }
    }
}

bool calculer_constante_modele(const Nom& turbulent_viscosity_dynamic_type,
                               const Nom& description,
                               const ArrOfDouble& moy_lij,
                               const ArrOfDouble& moy_mij,
                               const ArrOfDouble& moy_hij,
                               const ArrOfDouble& moy_mijmij,
                               const ArrOfDouble& moy_hijhij,
                               const ArrOfDouble& moy_mijlij,
                               const ArrOfDouble& moy_hijlij,
                               const ArrOfDouble& moy_mijhij,
                               const FixedVector<IJK_Field_double, 3>& velocity,
                               const ArrOfDouble_with_ghost& delta_z,
                               ArrOfDouble_with_ghost& constante_modele)
{
  bool mot_cle_reconnu = 0;

  const bool clipping = turbulent_viscosity_dynamic_type.finit_par("clipping");
  Nom dynamic_type = clipping ? turbulent_viscosity_dynamic_type.getPrefix("clipping") : turbulent_viscosity_dynamic_type;

  const bool global = dynamic_type.finit_par("global");
  dynamic_type = global ? dynamic_type.getPrefix("global") : dynamic_type;

  if ( dynamic_type == Nom("direct")
       || dynamic_type == Nom("directmixte") )
    {
      mot_cle_reconnu = 1;
      const bool mixte = ( dynamic_type == Nom("directmixte") );

      calculer_constante_direct(global, clipping,
                                description,
                                moy_mij, moy_lij,
                                mixte, moy_hij,
                                velocity, delta_z,
                                constante_modele);
    }
  else if ( dynamic_type == Nom("lilly")
            || dynamic_type == Nom("lillymixte") )
    {
      mot_cle_reconnu = 1;
      const bool mixte = ( dynamic_type == Nom("lillymixte") );

      calculer_constante_lilly(global, clipping,
                               description,
                               moy_mijmij, moy_mijlij,
                               mixte, moy_mijhij,
                               velocity, delta_z,
                               constante_modele);
    }
  else if ( dynamic_type == Nom("twoparameters") )
    {
      mot_cle_reconnu = 1;
      calculer_constante_twoparameters(global, clipping,
                                       description,
                                       moy_mijmij, moy_hijhij,
                                       moy_mijlij, moy_hijlij,
                                       moy_mijhij,
                                       velocity, delta_z,
                                       constante_modele);
    }
  else if ( dynamic_type == Nom("twonoerror") )
    {
      mot_cle_reconnu = 1;
      calculer_constante_twonoerror(global, clipping,
                                    description,
                                    moy_lij, moy_mij, moy_hij,
                                    moy_mijmij, moy_hijhij,
                                    moy_mijlij, moy_hijlij,
                                    moy_mijhij,
                                    velocity, delta_z,
                                    constante_modele);
    }
  return mot_cle_reconnu;
}

template<class T>
void calculer_ml_dynamic_uscalar_vector(const bool anisotropic,
                                        const bool vectorial,
                                        const FixedVector<IJK_Field_double, 3>& velocity,
                                        const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                        const IJK_Field_double& champ_scalar,
                                        const IJK_Field_double& champ_scalar_filtre,
                                        double scalar_kmin, double scalar_kmax,
                                        const int turbulent_viscosity,
                                        const IJK_Field_double& turbulent_mu_x,
                                        const IJK_Field_double& turbulent_mu_y,
                                        const IJK_Field_double& turbulent_mu_z,
                                        const IJK_Field_double& turbulent_mu_filtre_x,
                                        const IJK_Field_double& turbulent_mu_filtre_y,
                                        const IJK_Field_double& turbulent_mu_filtre_z,
                                        const int structural_uscalar,
                                        const FixedVector<IJK_Field_double, 3>& structural_uscalar_vector,
                                        const FixedVector<IJK_Field_double, 3>& structural_uscalar_filtre_vector,
                                        const ArrOfDouble_with_ghost& delta_z,
                                        const double facteur_delta_x,
                                        const double facteur_delta_y,
                                        const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                        const double facteur_delta_filtre_x,
                                        const double facteur_delta_filtre_y,
                                        const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                        T& kernel,
                                        FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                        FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                        FixedVector<FixedVector<ArrOfDouble, 7>, 8>& ml)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const IJK_Field_double& vitesse_i_filtre = velocity_filtre[0];
  const IJK_Field_double& vitesse_j_filtre = velocity_filtre[1];
  const IJK_Field_double& vitesse_k_filtre = velocity_filtre[2];

  const IJK_Field_double& structural_uscalar_x = structural_uscalar_vector[0];
  const IJK_Field_double& structural_uscalar_y = structural_uscalar_vector[1];
  const IJK_Field_double& structural_uscalar_z = structural_uscalar_vector[2];

  const IJK_Field_double& structural_uscalar_filtre_x = structural_uscalar_filtre_vector[0];
  const IJK_Field_double& structural_uscalar_filtre_y = structural_uscalar_filtre_vector[1];
  const IJK_Field_double& structural_uscalar_filtre_z = structural_uscalar_filtre_vector[2];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& bf_i = tmp_b[0];
  IJK_Field_local_double& bf_j = tmp_b[1];
  IJK_Field_local_double& bf_k = tmp_b[2];
  IJK_Field_local_double& bus_i = tmp_b[3];
  IJK_Field_local_double& bus_j = tmp_b[4];
  IJK_Field_local_double& bus_k = tmp_b[5];
  IJK_Field_local_double& bg_i = tmp_b[6];
  IJK_Field_local_double& bg_j = tmp_b[7];
  IJK_Field_local_double& bg_k = tmp_b[8];

  IJK_Field_local_double& af_i = tmp_a[0];
  IJK_Field_local_double& af_j = tmp_a[1];
  IJK_Field_local_double& af_k = tmp_a[2];
  IJK_Field_local_double& aus_i = tmp_a[3];
  IJK_Field_local_double& aus_j = tmp_a[4];
  IJK_Field_local_double& aus_k = tmp_a[5];
  IJK_Field_local_double& ag_i = tmp_a[6];
  IJK_Field_local_double& ag_j = tmp_a[7];
  IJK_Field_local_double& ag_k = tmp_a[8];

  ArrOfDouble& moy_li_x = ml[0][0];
  ArrOfDouble& moy_li_y = ml[0][1];
  ArrOfDouble& moy_li_z = ml[0][2];

  ArrOfDouble& moy_mi_x = ml[1][0];
  ArrOfDouble& moy_mi_y = ml[1][1];
  ArrOfDouble& moy_mi_z = ml[1][2];

  ArrOfDouble& moy_hi_x = ml[2][0];
  ArrOfDouble& moy_hi_y = ml[2][1];
  ArrOfDouble& moy_hi_z = ml[2][2];

  ArrOfDouble& moy_mimi_x = ml[3][0];
  ArrOfDouble& moy_mimi_y = ml[3][1];
  ArrOfDouble& moy_mimi_z = ml[3][2];
  ArrOfDouble& moy_mimi = ml[3][6];

  ArrOfDouble& moy_hihi_x = ml[4][0];
  ArrOfDouble& moy_hihi_y = ml[4][1];
  ArrOfDouble& moy_hihi_z = ml[4][2];
  ArrOfDouble& moy_hihi = ml[4][6];

  ArrOfDouble& moy_mili_x = ml[5][0];
  ArrOfDouble& moy_mili_y = ml[5][1];
  ArrOfDouble& moy_mili_z = ml[5][2];
  ArrOfDouble& moy_mili = ml[5][6];

  ArrOfDouble& moy_hili_x = ml[6][0];
  ArrOfDouble& moy_hili_y = ml[6][1];
  ArrOfDouble& moy_hili_z = ml[6][2];
  ArrOfDouble& moy_hili = ml[6][6];

  ArrOfDouble& moy_mihi_x = ml[7][0];
  ArrOfDouble& moy_mihi_y = ml[7][1];
  ArrOfDouble& moy_mihi_z = ml[7][2];
  ArrOfDouble& moy_mihi = ml[7][6];

  double m_i = 0.;
  double m_j = 0.;
  double m_k = 0.;

  double h_i = 0.;
  double h_j = 0.;
  double h_k = 0.;

  FixedVector<double, 3> f;
  FixedVector<double, 3> f_filtre;

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      if (vectorial)
        {
          moy_li_x[kg] = 0;
          moy_li_y[kg] = 0;
          moy_li_z[kg] = 0;

          moy_mi_x[kg] = 0;
          moy_mi_y[kg] = 0;
          moy_mi_z[kg] = 0;

          moy_hi_x[kg] = 0;
          moy_hi_y[kg] = 0;
          moy_hi_z[kg] = 0;

          moy_mimi_x[kg] = 0;
          moy_mimi_y[kg] = 0;
          moy_mimi_z[kg] = 0;

          moy_hihi_x[kg] = 0;
          moy_hihi_y[kg] = 0;
          moy_hihi_z[kg] = 0;

          moy_mili_x[kg] = 0;
          moy_mili_y[kg] = 0;
          moy_mili_z[kg] = 0;

          moy_hili_x[kg] = 0;
          moy_hili_y[kg] = 0;
          moy_hili_z[kg] = 0;

          moy_mihi_x[kg] = 0;
          moy_mihi_y[kg] = 0;
          moy_mihi_z[kg] = 0;
        }

      moy_mimi[kg] = 0;
      moy_hihi[kg] = 0;
      moy_mili[kg] = 0;
      moy_hili[kg] = 0;
      moy_mihi[kg] = 0;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];

      const FixedVector<double, 21> filter_kernel_z = kernel->inhomogeneous(true, k, kg, nktot, dz_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_elem = kernel->size_k_elem(kg, nktot);
      const int shift_k_elem = kernel->shift_k_elem(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_elem = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_elem = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<-1 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double filter_coef_z = filter_kernel_z[kp+10];
                  longueur_elem += filter_coef_z * dz;
                }
              facteur_elem = 1./longueur_elem;
            }
          else
            {
              facteur_elem = dz_glo==0. ? 0. : 1./dz_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              bus_i(i, j, 0) = 0.;
              bus_j(i, j, 0) = 0.;
              bus_k(i, j, 0) = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_elem(kpg, nktot)))
                    {
                      const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                      const double filter_coef_z = filter_kernel_z[kp+10];
                      const double scalar = champ_scalar(i,j,k+kp);

                      const double uf_i = vitesse_i(i,j,k+kp);
                      const double vf_j = vitesse_j(i,j,k+kp);
                      const double wf_k = vitesse_k(i,j,k+kp);

                      const double uf_ip1 = vitesse_i(i+1,j,k+kp);
                      const double vf_jp1 = vitesse_j(i,j+1,k+kp);
                      const double wf_kp1 = kpg==(nktot-1) ? 0. : vitesse_k(i,j,k+1+kp);

                      const double ue = 0.5 * (uf_i + uf_ip1);
                      const double ve = 0.5 * (vf_j + vf_jp1);
                      const double we = 0.5 * (wf_k + wf_kp1);

                      if (ponderation_filter_kernel)
                        {
                          bus_i(i, j, 0) += ue*scalar * filter_coef_z * dz * facteur_elem;
                          bus_j(i, j, 0) += ve*scalar * filter_coef_z * dz * facteur_elem;
                          bus_k(i, j, 0) += we*scalar * filter_coef_z * dz * facteur_elem;
                        }
                      else
                        {
                          bus_i(i, j, 0) += ue*scalar * filter_coef_z;
                          bus_j(i, j, 0) += ve*scalar * filter_coef_z;
                          bus_k(i, j, 0) += we*scalar * filter_coef_z;
                        }
                    }
                }

              if (turbulent_viscosity)
                {
                  bf_i(i, j, 0) = 0.;
                  bf_j(i, j, 0) = 0.;
                  bf_k(i, j, 0) = 0.;
                  for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                    {
                      const int kpg = kg + kp;
                      if (!(kernel->is_at_wall_elem(kpg, nktot)))
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double filter_coef_z = filter_kernel_z[kp+10];
                          if (!anisotropic)
                            {
                              calculer_pi(i, j, k+kp, velocity, champ_scalar, scalar_kmin, scalar_kmax, turbulent_mu_x, turbulent_mu_y, turbulent_mu_z, delta_z, 1., 1., 1., 1., 1., f);
                            }
                          else
                            {
                              const double dz_pour_delta = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z_pour_delta[k+kp];
                              const double dz_m1_pour_delta = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z_pour_delta[k-1+kp];
                              const double dz_p1_pour_delta = (kpg+1<0 || kpg+1>(nktot-1)) ? 0. : delta_z_pour_delta[k+1+kp];
                              const double delta_m_pour_delta = (kpg-1<0 || kpg>(nktot-1)) ? 0. : 0.5*(dz_pour_delta + dz_m1_pour_delta);
                              const double delta_p_pour_delta = (kpg<0 || kpg+1>(nktot-1)) ? 0. : 0.5*(dz_pour_delta + dz_p1_pour_delta);
                              calculer_pi(i, j, k+kp, velocity, champ_scalar, scalar_kmin, scalar_kmax, turbulent_mu_x, turbulent_mu_y, turbulent_mu_z, delta_z, dx_pour_delta, dy_pour_delta, dz_pour_delta, delta_m_pour_delta, delta_p_pour_delta, f);
                            }

                          if (ponderation_filter_kernel)
                            {
                              bf_i(i, j, 0) += f[0] * filter_coef_z * dz * facteur_elem;
                              bf_j(i, j, 0) += f[1] * filter_coef_z * dz * facteur_elem;
                              bf_k(i, j, 0) += f[2] * filter_coef_z * dz * facteur_elem;
                            }
                          else
                            {
                              bf_i(i, j, 0) += f[0] * filter_coef_z;
                              bf_j(i, j, 0) += f[1] * filter_coef_z;
                              bf_k(i, j, 0) += f[2] * filter_coef_z;
                            }
                        }
                    }
                }

              if (structural_uscalar)
                {
                  bg_i(i, j, 0) = 0.;
                  bg_j(i, j, 0) = 0.;
                  bg_k(i, j, 0) = 0.;
                  for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                    {
                      const int kpg = kg + kp;
                      if (!(kernel->is_at_wall_elem(kpg, nktot)))
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double filter_coef_z = filter_kernel_z[kp+10];
                          const double gxf_i = -structural_uscalar_x(i,j,k+kp);
                          const double gyf_j = -structural_uscalar_y(i,j,k+kp);
                          const double gzf_k = -structural_uscalar_z(i,j,k+kp);

                          const double gxf_ip1 = -structural_uscalar_x(i+1,j,k+kp);
                          const double gyf_jp1 = -structural_uscalar_y(i,j+1,k+kp);
                          const double gzf_kp1 = kpg==(nktot-1) ? 0. : -structural_uscalar_z(i,j,k+1+kp);

                          const double gx_e = 0.5 * (gxf_i + gxf_ip1);
                          const double gy_e = 0.5 * (gyf_j + gyf_jp1);
                          const double gz_e = 0.5 * (gzf_k + gzf_kp1);

                          if (ponderation_filter_kernel)
                            {
                              bg_i(i, j, 0) += gx_e * filter_coef_z * dz * facteur_elem;
                              bg_j(i, j, 0) += gy_e * filter_coef_z * dz * facteur_elem;
                              bg_k(i, j, 0) += gz_e * filter_coef_z * dz * facteur_elem;
                            }
                          else
                            {
                              bg_i(i, j, 0) += gx_e * filter_coef_z;
                              bg_j(i, j, 0) += gy_e * filter_coef_z;
                              bg_k(i, j, 0) += gz_e * filter_coef_z;
                            }
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              aus_i(i, 0, 0) = 0.;
              aus_j(i, 0, 0) = 0.;
              aus_k(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  aus_i(i, 0, 0) += bus_i(i, j+jp, 0) * filter_coef_y;
                  aus_j(i, 0, 0) += bus_j(i, j+jp, 0) * filter_coef_y;
                  aus_k(i, 0, 0) += bus_k(i, j+jp, 0) * filter_coef_y;
                }
              if (turbulent_viscosity)
                {
                  af_i(i, 0, 0) = 0.;
                  af_j(i, 0, 0) = 0.;
                  af_k(i, 0, 0) = 0.;
                  for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                    {
                      const double filter_coef_y = filter_kernel_y[jp+10];
                      af_i(i, 0, 0) += bf_i(i, j+jp, 0) * filter_coef_y;
                      af_j(i, 0, 0) += bf_j(i, j+jp, 0) * filter_coef_y;
                      af_k(i, 0, 0) += bf_k(i, j+jp, 0) * filter_coef_y;
                    }
                }
              if (structural_uscalar)
                {
                  ag_i(i, 0, 0) = 0.;
                  ag_j(i, 0, 0) = 0.;
                  ag_k(i, 0, 0) = 0.;
                  for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                    {
                      const double filter_coef_y = filter_kernel_y[jp+10];
                      ag_i(i, 0, 0) += bg_i(i, j+jp, 0) * filter_coef_y;
                      ag_j(i, 0, 0) += bg_j(i, j+jp, 0) * filter_coef_y;
                      ag_k(i, 0, 0) += bg_k(i, j+jp, 0) * filter_coef_y;
                    }
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double rus_i = 0.;
              double rus_j = 0.;
              double rus_k = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  rus_i += aus_i(i+ip, 0, 0) * filter_coef_x;
                  rus_j += aus_j(i+ip, 0, 0) * filter_coef_x;
                  rus_k += aus_k(i+ip, 0, 0) * filter_coef_x;
                }

              const double scalar = champ_scalar_filtre(i,j,k);

              const double uf_i = vitesse_i_filtre(i,j,k);
              const double vf_j = vitesse_j_filtre(i,j,k);
              const double wf_k = vitesse_k_filtre(i,j,k);

              const double uf_ip1 = vitesse_i_filtre(i+1,j,k);
              const double vf_jp1 = vitesse_j_filtre(i,j+1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k_filtre(i,j,k+1);

              const double ue = 0.5 * (uf_i + uf_ip1);
              const double ve = 0.5 * (vf_j + vf_jp1);
              const double we = 0.5 * (wf_k + wf_kp1);

              const double l_i = rus_i - ue*scalar;
              const double l_j = rus_j - ve*scalar;
              const double l_k = rus_k - we*scalar;

              if (vectorial)
                {
                  moy_li_x[kg] += l_i;
                  moy_li_y[kg] += l_j;
                  moy_li_z[kg] += l_k;
                }

              if (turbulent_viscosity)
                {
                  double rf_i = 0.;
                  double rf_j = 0.;
                  double rf_k = 0.;
                  for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                    {
                      const double filter_coef_x = filter_kernel_x[ip+10];
                      rf_i += af_i(i+ip, 0, 0) * filter_coef_x;
                      rf_j += af_j(i+ip, 0, 0) * filter_coef_x;
                      rf_k += af_k(i+ip, 0, 0) * filter_coef_x;
                    }

                  if (!anisotropic)
                    {
                      calculer_pi(i, j, k, velocity_filtre, champ_scalar_filtre, scalar_kmin, scalar_kmax, turbulent_mu_filtre_x, turbulent_mu_filtre_y, turbulent_mu_filtre_z, delta_z, 1., 1., 1., 1., 1., f_filtre);
                    }
                  else
                    {
                      const double dx_filtre = facteur_delta_filtre_x*dx;
                      const double dy_filtre = facteur_delta_filtre_y*dy;
                      const double dz_filtre = delta_z_pour_delta_filtre[k];
                      const double delta_m_filtre = kg==0 ? 0. : 0.5*(dz_filtre + delta_z_pour_delta_filtre[k-1]);
                      const double delta_p_filtre = kg==(nktot-1) ? 0. : 0.5*(dz_filtre + delta_z_pour_delta_filtre[k+1]);
                      calculer_pi(i, j, k, velocity_filtre, champ_scalar_filtre, scalar_kmin, scalar_kmax, turbulent_mu_filtre_x, turbulent_mu_filtre_y, turbulent_mu_filtre_z, delta_z, dx_filtre, dy_filtre, dz_filtre, delta_m_filtre, delta_p_filtre, f_filtre);
                    }

                  m_i = f_filtre[0] - rf_i;
                  m_j = f_filtre[1] - rf_j;
                  m_k = f_filtre[2] - rf_k;

                  const double mimi = m_i * m_i
                                      + m_j * m_j
                                      + m_k * m_k;

                  const double mili = m_i * l_i
                                      + m_j * l_j
                                      + m_k * l_k;

                  moy_mimi[kg] += mimi;
                  moy_mili[kg] += mili;

                  if (vectorial)
                    {
                      moy_mi_x[kg] += m_i;
                      moy_mi_y[kg] += m_j;
                      moy_mi_z[kg] += m_k;

                      moy_mimi_x[kg] += m_i * m_i;
                      moy_mimi_y[kg] += m_j * m_j;
                      moy_mimi_z[kg] += m_k * m_k;

                      moy_mili_x[kg] += m_i * l_i;
                      moy_mili_y[kg] += m_j * l_j;
                      moy_mili_z[kg] += m_k * l_k;
                    }
                }

              if (structural_uscalar)
                {
                  double rg_i = 0.;
                  double rg_j = 0.;
                  double rg_k = 0.;
                  for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                    {
                      const double filter_coef_x = filter_kernel_x[ip+10];
                      rg_i += ag_i(i+ip, 0, 0) * filter_coef_x;
                      rg_j += ag_j(i+ip, 0, 0) * filter_coef_x;
                      rg_k += ag_k(i+ip, 0, 0) * filter_coef_x;
                    }

                  const double gxf_i = -structural_uscalar_filtre_x(i,j,k);
                  const double gyf_j = -structural_uscalar_filtre_y(i,j,k);
                  const double gzf_k = -structural_uscalar_filtre_z(i,j,k);

                  const double gxf_ip1 = -structural_uscalar_filtre_x(i+1,j,k);
                  const double gyf_jp1 = -structural_uscalar_filtre_y(i,j+1,k);
                  const double gzf_kp1 = kg==(nktot-1) ? 0. : -structural_uscalar_filtre_z(i,j,k+1);

                  const double gx_e = 0.5 * (gxf_i + gxf_ip1);
                  const double gy_e = 0.5 * (gyf_j + gyf_jp1);
                  const double gz_e = 0.5 * (gzf_k + gzf_kp1);

                  h_i = gx_e - rg_i;
                  h_j = gy_e - rg_j;
                  h_k = gz_e - rg_k;

                  const double hihi = h_i * h_i
                                      + h_j * h_j
                                      + h_k * h_k;

                  const double hili = h_i * l_i
                                      + h_j * l_j
                                      + h_k * l_k;

                  moy_hihi[kg] += hihi;
                  moy_hili[kg] += hili;

                  if (vectorial)
                    {
                      moy_hi_x[kg] += h_i;
                      moy_hi_y[kg] += h_j;
                      moy_hi_z[kg] += h_k;

                      moy_hihi_x[kg] += h_i * h_i;
                      moy_hihi_y[kg] += h_j * h_j;
                      moy_hihi_z[kg] += h_k * h_k;

                      moy_hili_x[kg] += h_i * l_i;
                      moy_hili_y[kg] += h_j * l_j;
                      moy_hili_z[kg] += h_k * l_k;
                    }

                  if (turbulent_viscosity)
                    {
                      const double mihi = m_i * h_i
                                          + m_j * h_j
                                          + m_k * h_k;

                      moy_mihi[kg] += mihi;

                      if (vectorial)
                        {
                          moy_mihi_x[kg] += m_i * h_i;
                          moy_mihi_y[kg] += m_j * h_j;
                          moy_mihi_z[kg] += m_k * h_k;
                        }
                    }
                }
            }
        }
    }
}

template<class T>
void calculer_turbulent_mu_scalar(const bool anisotropic,
                                  const Nom& turbulent_viscosity_model,
                                  const double turbulent_viscosity_model_constant,
                                  const double variation_cste_modele_fonctionnel_,
                                  const double smoothing_center_fr_,
                                  const double smoothing_factor_fr_,
                                  const double Re_tau_fr_,
                                  const double Re_tau_ch_,
                                  const double pond_fr_,
                                  const double pond_ch_,
                                  const double center_constant_,
                                  const double Lz_tot_,
                                  FixedVector<IJK_Field_double, 3>& velocity,
                                  FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                  IJK_Field_double& rho,
                                  IJK_Field_double& rho_filtre,
                                  double rho_kmin, double rho_kmax,
                                  IJK_Field_double& scalar,
                                  IJK_Field_double& scalar_filtre,
                                  double scalar_kmin, double scalar_kmax,
                                  const ArrOfDouble_with_ghost& delta_z,
                                  const double facteur_delta_x,
                                  const double facteur_delta_y,
                                  const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                  const double facteur_delta_filtre_x,
                                  const double facteur_delta_filtre_y,
                                  const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                  T& kernel,
                                  FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                  FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                  const bool flag_turbulent_mu_filtre,
                                  IJK_Field_double& turbulent_mu_filtre,
                                  IJK_Field_double& turbulent_mu,
                                  const IJK_Splitting& splitting)
{
  Turbulent_viscosity_base* model = NULL;

  choix_modele(turbulent_viscosity_model, model);

  calculer_turbulent_mu(anisotropic,
                        *model, turbulent_viscosity_model_constant, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                        velocity, rho, scalar, scalar_kmin, scalar_kmax,
                        delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                        turbulent_mu, splitting);
  int ghost_size_filter;
  int ghost_size_velocity;
  int ghost_size_scalar;
  if (flag_turbulent_mu_filtre)
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      ghost_size_scalar = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);
      rho.echange_espace_virtuel(ghost_size_scalar);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);
      filtrer_champ_elem(flag_add, rho, delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, rho_filtre);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);
      rho_filtre.echange_espace_virtuel(2);

      // si scalar et rho ne sont pas le meme objet
      if (&scalar != &rho)
        {
          scalar.echange_espace_virtuel(ghost_size_scalar);
          filtrer_champ_elem(flag_add, scalar, delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, scalar_filtre);
          scalar_filtre.echange_espace_virtuel(2);
        }

      calculer_turbulent_mu(anisotropic,
                            *model, turbulent_viscosity_model_constant, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                            velocity_filtre, rho_filtre, scalar_filtre, scalar_kmin, scalar_kmax,
                            delta_z, facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                            turbulent_mu_filtre, splitting);
    }

  delete model;
}

template<class T>
void calculer_turbulent_mu_tensor(const bool anisotropic,
                                  const Nom& turbulent_viscosity_model,
                                  const double turbulent_viscosity_model_constant,
                                  const ArrOfDouble& turbulent_viscosity_tensor_coefficients,
                                  const double variation_cste_modele_fonctionnel_,
                                  const double smoothing_center_fr_,
                                  const double smoothing_factor_fr_,
                                  const double Re_tau_fr_,
                                  const double Re_tau_ch_,
                                  const double pond_fr_,
                                  const double pond_ch_,
                                  const double center_constant_,
                                  const double Lz_tot_,
                                  FixedVector<IJK_Field_double, 3>& velocity,
                                  FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                  IJK_Field_double& rho,
                                  IJK_Field_double& rho_filtre,
                                  double rho_kmin, double rho_kmax,
                                  IJK_Field_double& scalar,
                                  IJK_Field_double& scalar_filtre,
                                  double scalar_kmin, double scalar_kmax,
                                  const ArrOfDouble_with_ghost& delta_z,
                                  const double facteur_delta_x,
                                  const double facteur_delta_y,
                                  const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                  const double facteur_delta_filtre_x,
                                  const double facteur_delta_filtre_y,
                                  const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                  T& kernel,
                                  FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                  FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                  const bool flag_turbulent_mu_filtre,
                                  FixedVector<IJK_Field_double, 6>& turbulent_mu_filtre_tensor,
                                  FixedVector<IJK_Field_double, 6>& turbulent_mu_tensor,
                                  const IJK_Splitting& splitting)
{
  Turbulent_viscosity_base* model = NULL;

  const double& coefficient_xx = turbulent_viscosity_tensor_coefficients[0];
  const double& coefficient_xy = turbulent_viscosity_tensor_coefficients[1];
  const double& coefficient_xz = turbulent_viscosity_tensor_coefficients[2];
  const double& coefficient_yy = turbulent_viscosity_tensor_coefficients[3];
  const double& coefficient_yz = turbulent_viscosity_tensor_coefficients[4];
  const double& coefficient_zz = turbulent_viscosity_tensor_coefficients[5];

  IJK_Field_double& turbulent_mu_xx = turbulent_mu_tensor[0];
  IJK_Field_double& turbulent_mu_xy = turbulent_mu_tensor[1];
  IJK_Field_double& turbulent_mu_xz = turbulent_mu_tensor[2];
  IJK_Field_double& turbulent_mu_yy = turbulent_mu_tensor[3];
  IJK_Field_double& turbulent_mu_yz = turbulent_mu_tensor[4];
  IJK_Field_double& turbulent_mu_zz = turbulent_mu_tensor[5];

  IJK_Field_double& turbulent_mu_filtre_xx = turbulent_mu_filtre_tensor[0];
  IJK_Field_double& turbulent_mu_filtre_xy = turbulent_mu_filtre_tensor[1];
  IJK_Field_double& turbulent_mu_filtre_xz = turbulent_mu_filtre_tensor[2];
  IJK_Field_double& turbulent_mu_filtre_yy = turbulent_mu_filtre_tensor[3];
  IJK_Field_double& turbulent_mu_filtre_yz = turbulent_mu_filtre_tensor[4];
  IJK_Field_double& turbulent_mu_filtre_zz = turbulent_mu_filtre_tensor[5];

  choix_modele(turbulent_viscosity_model, model);

  if (0)
    {
      // do nothing
    }
  else
    {
      calculer_turbulent_mu_scalar(anisotropic,
                                   turbulent_viscosity_model, turbulent_viscosity_model_constant, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                   velocity, velocity_filtre,
                                   rho, rho_filtre, rho_kmin, rho_kmax,
                                   scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                   delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                   facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                   kernel, tmp_b, tmp_a,
                                   flag_turbulent_mu_filtre, turbulent_mu_filtre_xx,
                                   turbulent_mu_xx, splitting);

      multiplier_champ(coefficient_xy, turbulent_mu_xx, turbulent_mu_xy);
      multiplier_champ(coefficient_xz, turbulent_mu_xx, turbulent_mu_xz);
      multiplier_champ(coefficient_yy, turbulent_mu_xx, turbulent_mu_yy);
      multiplier_champ(coefficient_yz, turbulent_mu_xx, turbulent_mu_yz);
      multiplier_champ(coefficient_zz, turbulent_mu_xx, turbulent_mu_zz);
      multiplier_champ(coefficient_xx, turbulent_mu_xx, turbulent_mu_xx);

      if (flag_turbulent_mu_filtre)
        {
          multiplier_champ(coefficient_xy, turbulent_mu_filtre_xx, turbulent_mu_filtre_xy);
          multiplier_champ(coefficient_xz, turbulent_mu_filtre_xx, turbulent_mu_filtre_xz);
          multiplier_champ(coefficient_yy, turbulent_mu_filtre_xx, turbulent_mu_filtre_yy);
          multiplier_champ(coefficient_yz, turbulent_mu_filtre_xx, turbulent_mu_filtre_yz);
          multiplier_champ(coefficient_zz, turbulent_mu_filtre_xx, turbulent_mu_filtre_zz);
          multiplier_champ(coefficient_xx, turbulent_mu_filtre_xx, turbulent_mu_filtre_xx);
        }
    }

  delete model;
}

template<class T>
void calculer_turbulent_mu_vector(const bool anisotropic,
                                  const Nom& turbulent_viscosity_model,
                                  const double turbulent_viscosity_model_constant,
                                  const ArrOfDouble& turbulent_diffusivity_vector_coefficients,
                                  const double variation_cste_modele_fonctionnel_,
                                  const double smoothing_center_fr_,
                                  const double smoothing_factor_fr_,
                                  const double Re_tau_fr_,
                                  const double Re_tau_ch_,
                                  const double pond_fr_,
                                  const double pond_ch_,
                                  const double center_constant_,
                                  const double Lz_tot_,
                                  FixedVector<IJK_Field_double, 3>& velocity,
                                  FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                  IJK_Field_double& rho,
                                  IJK_Field_double& rho_filtre,
                                  double rho_kmin, double rho_kmax,
                                  IJK_Field_double& scalar,
                                  IJK_Field_double& scalar_filtre,
                                  double scalar_kmin, double scalar_kmax,
                                  const ArrOfDouble_with_ghost& delta_z,
                                  const double facteur_delta_x,
                                  const double facteur_delta_y,
                                  const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                  const double facteur_delta_filtre_x,
                                  const double facteur_delta_filtre_y,
                                  const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                  T& kernel,
                                  FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                  FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                  const bool flag_turbulent_mu_filtre,
                                  FixedVector<IJK_Field_double, 3>& turbulent_mu_filtre_vector,
                                  FixedVector<IJK_Field_double, 3>& turbulent_mu_vector,
                                  const IJK_Splitting& splitting)
{
  Turbulent_viscosity_base* model = NULL;

  const double& coefficient_x = turbulent_diffusivity_vector_coefficients[0];
  const double& coefficient_y = turbulent_diffusivity_vector_coefficients[1];
  const double& coefficient_z = turbulent_diffusivity_vector_coefficients[2];

  IJK_Field_double& turbulent_mu_x = turbulent_mu_vector[0];
  IJK_Field_double& turbulent_mu_y = turbulent_mu_vector[1];
  IJK_Field_double& turbulent_mu_z = turbulent_mu_vector[2];

  IJK_Field_double& turbulent_mu_filtre_x = turbulent_mu_filtre_vector[0];
  IJK_Field_double& turbulent_mu_filtre_y = turbulent_mu_filtre_vector[1];
  IJK_Field_double& turbulent_mu_filtre_z = turbulent_mu_filtre_vector[2];

  choix_modele(turbulent_viscosity_model, model);

  if (0)
    {
      // do nothing
    }
  else
    {
      calculer_turbulent_mu_scalar(anisotropic,
                                   turbulent_viscosity_model, turbulent_viscosity_model_constant, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                   velocity, velocity_filtre,
                                   rho, rho_filtre, rho_kmin, rho_kmax,
                                   scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                   delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                   facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                   kernel, tmp_b, tmp_a,
                                   flag_turbulent_mu_filtre, turbulent_mu_filtre_x,
                                   turbulent_mu_x, splitting);

      multiplier_champ(coefficient_y, turbulent_mu_x, turbulent_mu_y);
      multiplier_champ(coefficient_z, turbulent_mu_x, turbulent_mu_z);
      multiplier_champ(coefficient_x, turbulent_mu_x, turbulent_mu_x);

      if (flag_turbulent_mu_filtre)
        {
          multiplier_champ(coefficient_y, turbulent_mu_filtre_x, turbulent_mu_filtre_y);
          multiplier_champ(coefficient_z, turbulent_mu_filtre_x, turbulent_mu_filtre_z);
          multiplier_champ(coefficient_x, turbulent_mu_filtre_x, turbulent_mu_filtre_x);
        }
    }

  delete model;
}

void calculer_structural_uu_gradient(const double structural_uu_model_constant,
                                     const ArrOfDouble& structural_uu_tensor_coefficients,
                                     const FixedVector<IJK_Field_double, 3>& velocity,
                                     const ArrOfDouble_with_ghost& delta_z_maillage,
                                     const double facteur_x_pour_delta,
                                     const double facteur_y_pour_delta,
                                     const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                     FixedVector<IJK_Field_double, 6>& structural_uu_tensor)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const double& coefficient_xx = structural_uu_tensor_coefficients[0];
  const double& coefficient_xy = structural_uu_tensor_coefficients[1];
  const double& coefficient_xz = structural_uu_tensor_coefficients[2];
  const double& coefficient_yy = structural_uu_tensor_coefficients[3];
  const double& coefficient_yz = structural_uu_tensor_coefficients[4];
  const double& coefficient_zz = structural_uu_tensor_coefficients[5];

  const double deltaunsurdx = facteur_x_pour_delta;
  const double deltaunsurdy = facteur_y_pour_delta;

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;

              const double dz = delta_z_maillage[k];
              const double dz_m1 = kg==0 ? 0. : delta_z_maillage[k-1];
              const double delta_m = kg==0 ? 0.5*dz : 0.5*(dz + delta_z_maillage[k-1]);
              const double delta_p = kg==(nktot-1) ? 0.5*dz : 0.5*(dz + delta_z_maillage[k+1]);

              const double deltaunsurdz = delta_z_pour_delta[k] * 1./dz;
              const double deltaunsurdz_m1 = kg==0 ? 0. : delta_z_pour_delta[k-1] * 1./dz_m1;
              const double deltaunsurdelta_m = kg==0 ? 0. : 0.5*(delta_z_pour_delta[k] + delta_z_pour_delta[k-1]) * 1./delta_m;
              const double deltaunsurdelta_p = kg==(nktot-1) ? 0. : 0.5*(delta_z_pour_delta[k] + delta_z_pour_delta[k+1]) * 1./delta_p;

              const double uf_i = vitesse_i(i,j,k);
              const double uf_ip1 = vitesse_i(i+1,j,k);
              const double uf_im1 = vitesse_i(i-1,j,k);
              const double uf_i_jm1 = vitesse_i(i,j-1,k);
              const double uf_i_jp1 = vitesse_i(i,j+1,k);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);
              const double uf_i_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i,j,k+1);
              const double uf_i_jm1_km1 = kg==0 ? 0. : vitesse_i(i,j-1,k-1);
              const double uf_i_jp1_km1 = kg==0 ? 0. : vitesse_i(i,j+1,k-1);
              const double uf_i_jm1_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i,j-1,k+1);
              const double uf_ip1_jp1 = vitesse_i(i+1,j+1,k);
              const double uf_ip1_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i+1,j,k+1);
              const double uf_ip1_jm1 = vitesse_i(i+1,j-1,k);
              const double uf_ip1_km1 = kg==0 ? 0. : vitesse_i(i+1,j,k-1);
              const double uf_im1_jm1 = vitesse_i(i-1,j-1,k);
              const double uf_im1_km1 = kg==0 ? 0. : vitesse_i(i-1,j,k-1);

              const double vf_j = vitesse_j(i,j,k);
              const double vf_j_im1 = vitesse_j(i-1,j,k);
              const double vf_j_ip1 = vitesse_j(i+1,j,k);
              const double vf_jp1 = vitesse_j(i,j+1,k);
              const double vf_jm1 = vitesse_j(i,j-1,k);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);
              const double vf_j_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j,k+1);
              const double vf_j_im1_km1 = kg==0 ? 0. : vitesse_j(i-1,j,k-1);
              const double vf_j_ip1_km1 = kg==0 ? 0. : vitesse_j(i+1,j,k-1);
              const double vf_j_im1_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i-1,j,k+1);
              const double vf_jp1_ip1 = vitesse_j(i+1,j+1,k);
              const double vf_jp1_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j+1,k+1);
              const double vf_jp1_im1 = vitesse_j(i-1,j+1,k);
              const double vf_jp1_km1 = kg==0 ? 0. : vitesse_j(i,j+1,k-1);
              const double vf_jm1_im1 = vitesse_j(i-1,j-1,k);
              const double vf_jm1_km1 = kg==0 ? 0. : vitesse_j(i,j-1,k-1);

              const double wf_k = vitesse_k(i,j,k);
              const double wf_k_im1 = vitesse_k(i-1,j,k);
              const double wf_k_ip1 = vitesse_k(i+1,j,k);
              const double wf_k_jm1 = vitesse_k(i,j-1,k);
              const double wf_k_jp1 = vitesse_k(i,j+1,k);
              const double wf_k_im1_jm1 = vitesse_k(i-1,j-1,k);
              const double wf_k_ip1_jm1 = vitesse_k(i+1,j-1,k);
              const double wf_k_im1_jp1 = vitesse_k(i-1,j+1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);
              const double wf_km1 = kg==0 ? 0. : vitesse_k(i,j,k-1);
              const double wf_kp1_ip1 = kg==(nktot-1) ? 0. : vitesse_k(i+1,j,k+1);
              const double wf_kp1_jp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j+1,k+1);
              const double wf_kp1_im1 = kg==(nktot-1) ? 0. : vitesse_k(i-1,j,k+1);
              const double wf_kp1_jm1 = kg==(nktot-1) ? 0. : vitesse_k(i,j-1,k+1);
              const double wf_km1_im1 = kg==0 ? 0. : vitesse_k(i-1,j,k-1);
              const double wf_km1_jm1 = kg==0 ? 0. : vitesse_k(i,j-1,k-1);

              const double duidx = deltaunsurdx * (uf_ip1 - uf_i);
              const double dujdy = deltaunsurdy * (vf_jp1 - vf_j);
              const double dukdz = deltaunsurdz * (wf_kp1 - wf_k);

              const double duidx_im1 = deltaunsurdx * (uf_i - uf_im1);
              const double dujdy_im1 = deltaunsurdy * (vf_jp1_im1 - vf_j_im1);
              const double dukdz_im1 = deltaunsurdz * (wf_kp1_im1 - wf_k_im1);

              const double duidx_jm1 = deltaunsurdx * (uf_ip1_jm1 - uf_i_jm1);
              const double dujdy_jm1 = deltaunsurdy * (vf_j - vf_jm1);
              const double dukdz_jm1 = deltaunsurdz * (wf_kp1_jm1 - wf_k_jm1);

              const double duidx_km1 = deltaunsurdx * (uf_ip1_km1 - uf_i_km1);
              const double dujdy_km1 = deltaunsurdy * (vf_jp1_km1 - vf_j_km1);
              const double dukdz_km1 = deltaunsurdz_m1 * (wf_k - wf_km1);

              const double duidx_im1_jm1 = deltaunsurdx * (uf_i_jm1 - uf_im1_jm1);
              const double dujdy_im1_jm1 = deltaunsurdy * (vf_j_im1 - vf_jm1_im1);

              const double dujdy_jm1_km1 = deltaunsurdy * (vf_j_km1 - vf_jm1_km1);
              const double dukdz_jm1_km1 = deltaunsurdz * (wf_k_jm1 - wf_km1_jm1);

              const double duidx_im1_km1 = deltaunsurdx * (uf_i_km1 - uf_im1_km1);
              const double dukdz_im1_km1 = deltaunsurdz * (wf_k_im1 - wf_km1_im1);

              const double duidy_ij = deltaunsurdy * (uf_i - uf_i_jm1);
              const double duidy_ip1j = deltaunsurdy * (uf_ip1 - uf_ip1_jm1);
              const double duidy_ijp1 = deltaunsurdy * (uf_i_jp1 - uf_i);
              const double duidy_ip1jp1 = deltaunsurdy * (uf_ip1_jp1 - uf_ip1);
              const double duidy_ij_km1 = deltaunsurdy * (uf_i_km1 - uf_i_jm1_km1);
              const double duidy_ijp1_km1 = deltaunsurdy * (uf_i_jp1_km1 - uf_i_km1);

              const double duidz_ik = deltaunsurdelta_m * (uf_i - uf_i_km1);
              const double duidz_ip1k = deltaunsurdelta_m * (uf_ip1 - uf_ip1_km1);
              const double duidz_ikp1 = deltaunsurdelta_p * (uf_i_kp1 - uf_i);
              const double duidz_ip1kp1 = deltaunsurdelta_p * (uf_ip1_kp1 - uf_ip1);
              const double duidz_ik_jm1 = deltaunsurdelta_m * (uf_i_jm1 - uf_i_jm1_km1);
              const double duidz_ikp1_jm1 = deltaunsurdelta_p * (uf_i_jm1_kp1 - uf_i_jm1);

              const double dujdx_ij = deltaunsurdx * (vf_j - vf_j_im1);
              const double dujdx_ip1j = deltaunsurdx * (vf_j_ip1 - vf_j);
              const double dujdx_ijp1 = deltaunsurdx * (vf_jp1 - vf_jp1_im1);
              const double dujdx_ip1jp1 = deltaunsurdx * (vf_jp1_ip1 - vf_jp1);
              const double dujdx_ij_km1 = deltaunsurdx * (vf_j_km1 - vf_j_im1_km1);
              const double dujdx_ip1j_km1 = deltaunsurdx * (vf_j_ip1_km1 - vf_j_km1);

              const double dujdz_jk = deltaunsurdelta_m * (vf_j - vf_j_km1);
              const double dujdz_jp1k = deltaunsurdelta_m * (vf_jp1 - vf_jp1_km1);
              const double dujdz_jkp1 = deltaunsurdelta_p * (vf_j_kp1 - vf_j);
              const double dujdz_jp1kp1 = deltaunsurdelta_p * (vf_jp1_kp1 - vf_jp1);
              const double dujdz_jk_im1 = deltaunsurdelta_m * (vf_j_im1 - vf_j_im1_km1);
              const double dujdz_jkp1_im1 = deltaunsurdelta_p * (vf_j_im1_kp1 - vf_j_im1);

              const double dukdx_ik = deltaunsurdx * (wf_k - wf_k_im1);
              const double dukdx_ip1k = deltaunsurdx * (wf_k_ip1 - wf_k);
              const double dukdx_ikp1 = deltaunsurdx * (wf_kp1 - wf_kp1_im1);
              const double dukdx_ip1kp1 = deltaunsurdx * (wf_kp1_ip1 - wf_kp1);
              const double dukdx_ik_jm1 = deltaunsurdx * (wf_k_jm1 - wf_k_im1_jm1);
              const double dukdx_ip1k_jm1 = deltaunsurdx * (wf_k_ip1_jm1 - wf_k_jm1);

              const double dukdy_jk = deltaunsurdy * (wf_k - wf_k_jm1);
              const double dukdy_jp1k = deltaunsurdy * (wf_k_jp1 - wf_k);
              const double dukdy_jkp1 = deltaunsurdy * (wf_kp1 - wf_kp1_jm1);
              const double dukdy_jp1kp1 = deltaunsurdy * (wf_kp1_jp1 - wf_kp1);
              const double dukdy_jk_im1 = deltaunsurdy * (wf_k_im1 - wf_k_im1_jm1);
              const double dukdy_jp1k_im1 = deltaunsurdy * (wf_k_im1_jp1 - wf_k_im1);

              const double g_e_ii = duidx;
              const double g_e_ij = 0.25 * (duidy_ip1jp1 + duidy_ijp1 + duidy_ip1j + duidy_ij);
              const double g_e_ik = 0.25 * (duidz_ip1kp1 + duidz_ikp1 + duidz_ip1k + duidz_ik);
              const double g_e_ji = 0.25 * (dujdx_ip1jp1 + dujdx_ip1j + dujdx_ijp1 + dujdx_ij);
              const double g_e_jj = dujdy;
              const double g_e_jk = 0.25 * (dujdz_jp1kp1 + dujdz_jkp1 + dujdz_jp1k + dujdz_jk);
              const double g_e_ki = 0.25 * (dukdx_ip1kp1 + dukdx_ip1k + dukdx_ikp1 + dukdx_ik);
              const double g_e_kj = 0.25 * (dukdy_jp1kp1 + dukdy_jp1k + dukdy_jkp1 + dukdy_jk);
              const double g_e_kk = dukdz;

              const double g_aij_ii = 0.25 * (duidx_im1_jm1 + duidx_jm1 + duidx_im1 + duidx);
              const double g_aij_ij = duidy_ij;
              const double g_aij_ik = 0.25 * (duidz_ikp1_jm1 + duidz_ikp1 + duidz_ik_jm1 + duidz_ik);
              const double g_aij_ji = dujdx_ij;
              const double g_aij_jj = 0.25 * (dujdy_im1_jm1 + dujdy_jm1 + dujdy_im1 + dujdy);
              const double g_aij_jk = 0.25 * (dujdz_jkp1_im1 + dujdz_jkp1 + dujdz_jk_im1 + dujdz_jk);

              const double g_aik_ii = 0.25 * (duidx_im1_km1 + duidx_km1 + duidx_im1 + duidx);
              const double g_aik_ij = 0.25 * (duidy_ijp1_km1 + duidy_ijp1 + duidy_ij_km1 + duidy_ij);
              const double g_aik_ik = duidz_ik;
              const double g_aik_ki = dukdx_ik;
              const double g_aik_kj = 0.25 * (dukdy_jp1k_im1 + dukdy_jp1k + dukdy_jk_im1 + dukdy_jk);
              const double g_aik_kk = 0.25 * (dukdz_im1_km1 + dukdz_km1 + dukdz_im1 + dukdz);

              const double g_ajk_ji = 0.25 * (dujdx_ip1j_km1 + dujdx_ip1j + dujdx_ij_km1 + dujdx_ij);
              const double g_ajk_jj = 0.25 * (dujdy_jm1_km1 + dujdy_km1 + dujdy_jm1 + dujdy);
              const double g_ajk_jk = dujdz_jk;
              const double g_ajk_ki = 0.25 * (dukdx_ip1k_jm1 + dukdx_ip1k + dukdx_ik_jm1 + dukdx_ik);
              const double g_ajk_kj = dukdy_jk;
              const double g_ajk_kk = 0.25 * (dukdz_jm1_km1 + dukdz_km1 + dukdz_jm1 + dukdz);

              const double c_ii = g_e_ii*g_e_ii + g_e_ij*g_e_ij + g_e_ik*g_e_ik;
              const double c_ij = g_aij_ii*g_aij_ji + g_aij_ij*g_aij_jj + g_aij_ik*g_aij_jk;
              const double c_ik = g_aik_ii*g_aik_ki + g_aik_ij*g_aik_kj + g_aik_ik*g_aik_kk;
              const double c_jj = g_e_ji*g_e_ji + g_e_jj*g_e_jj + g_e_jk*g_e_jk;
              const double c_jk = g_ajk_ji*g_ajk_ki + g_ajk_jj*g_ajk_kj + g_ajk_jk*g_ajk_kk;
              const double c_kk = g_e_ki*g_e_ki + g_e_kj*g_e_kj + g_e_kk*g_e_kk;

              structural_uu_xx(i,j,k) = - coefficient_xx * structural_uu_model_constant * c_ii/12.;
              structural_uu_xy(i,j,k) = - coefficient_xy * structural_uu_model_constant * c_ij/12.;
              structural_uu_xz(i,j,k) = - coefficient_xz * structural_uu_model_constant * c_ik/12.;
              structural_uu_yy(i,j,k) = - coefficient_yy * structural_uu_model_constant * c_jj/12.;
              structural_uu_yz(i,j,k) = - coefficient_yz * structural_uu_model_constant * c_jk/12.;
              structural_uu_zz(i,j,k) = - coefficient_zz * structural_uu_model_constant * c_kk/12.;
            }
        }
    }
}

void calculer_laplacien_u(const FixedVector<IJK_Field_double, 3>& velocity,
                          const ArrOfDouble_with_ghost& delta_z_maillage,
                          const double facteur_x_pour_delta,
                          const double facteur_y_pour_delta,
                          const ArrOfDouble_with_ghost& delta_z_pour_delta,
                          FixedVector<IJK_Field_double, 3>& laplacien_velocity)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  IJK_Field_double& laplacien_i = laplacien_velocity[0];
  IJK_Field_double& laplacien_j = laplacien_velocity[1];
  IJK_Field_double& laplacien_k = laplacien_velocity[2];

  const double deltaunsurdx = facteur_x_pour_delta;
  const double deltaunsurdy = facteur_y_pour_delta;

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;

              const double dz = delta_z_maillage[k];
              const double dz_m1 = kg==0 ? 0. : delta_z_maillage[k-1];
              const double delta_m = kg==0 ? 0.5*dz : 0.5*(dz + delta_z_maillage[k-1]);
              const double delta_p = kg==(nktot-1) ? 0.5*dz : 0.5*(dz + delta_z_maillage[k+1]);

              const double deltaunsurdz = delta_z_pour_delta[k] * 1./dz;
              const double deltaunsurdz_m1 = kg==0 ? 0. : delta_z_pour_delta[k-1] * 1./dz_m1;
              const double deltaunsurdelta_m = kg==0 ? 0. : 0.5*(delta_z_pour_delta[k] + delta_z_pour_delta[k-1]) * 1./delta_m;
              const double deltaunsurdelta_p = kg==(nktot-1) ? 0. : 0.5*(delta_z_pour_delta[k] + delta_z_pour_delta[k+1]) * 1./delta_p;

              const double uf_i = vitesse_i(i,j,k);
              const double uf_ip1 = vitesse_i(i+1,j,k);
              const double uf_im1 = vitesse_i(i-1,j,k);
              const double uf_i_jm1 = vitesse_i(i,j-1,k);
              const double uf_i_jp1 = vitesse_i(i,j+1,k);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);
              const double uf_i_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i,j,k+1);

              const double vf_j = vitesse_j(i,j,k);
              const double vf_j_im1 = vitesse_j(i-1,j,k);
              const double vf_j_ip1 = vitesse_j(i+1,j,k);
              const double vf_jp1 = vitesse_j(i,j+1,k);
              const double vf_jm1 = vitesse_j(i,j-1,k);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);
              const double vf_j_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j,k+1);

              const double wf_k = vitesse_k(i,j,k);
              const double wf_k_im1 = vitesse_k(i-1,j,k);
              const double wf_k_ip1 = vitesse_k(i+1,j,k);
              const double wf_k_jm1 = vitesse_k(i,j-1,k);
              const double wf_k_jp1 = vitesse_k(i,j+1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);
              const double wf_km1 = kg==0 ? 0. : vitesse_k(i,j,k-1);

              const double duidx = deltaunsurdx * (uf_ip1 - uf_i);
              const double dujdy = deltaunsurdy * (vf_jp1 - vf_j);
              const double dukdz = deltaunsurdz * (wf_kp1 - wf_k);

              const double duidx_im1 = deltaunsurdx * (uf_i - uf_im1);
              const double dujdy_jm1 = deltaunsurdy * (vf_j - vf_jm1);
              const double dukdz_km1 = deltaunsurdz_m1 * (wf_k - wf_km1);

              const double duidy_ij = deltaunsurdy * (uf_i - uf_i_jm1);
              const double duidy_ijp1 = deltaunsurdy * (uf_i_jp1 - uf_i);

              const double duidz_ik = deltaunsurdelta_m * (uf_i - uf_i_km1);
              const double duidz_ikp1 = deltaunsurdelta_p * (uf_i_kp1 - uf_i);

              const double dujdx_ij = deltaunsurdx * (vf_j - vf_j_im1);
              const double dujdx_ip1j = deltaunsurdx * (vf_j_ip1 - vf_j);

              const double dujdz_jk = deltaunsurdelta_m * (vf_j - vf_j_km1);
              const double dujdz_jkp1 = deltaunsurdelta_p * (vf_j_kp1 - vf_j);

              const double dukdx_ik = deltaunsurdx * (wf_k - wf_k_im1);
              const double dukdx_ip1k = deltaunsurdx * (wf_k_ip1 - wf_k);

              const double dukdy_jk = deltaunsurdy * (wf_k - wf_k_jm1);
              const double dukdy_jp1k = deltaunsurdy * (wf_k_jp1 - wf_k);

              const double d2uidx2f_i = deltaunsurdx * (duidx - duidx_im1);
              const double d2ujdy2f_j = deltaunsurdy * (dujdy - dujdy_jm1);
              const double d2ukdz2f_k = deltaunsurdelta_m * (dukdz - dukdz_km1);

              const double d2uidy2f_i = deltaunsurdy * (duidy_ijp1 - duidy_ij);
              const double d2uidz2f_i = deltaunsurdz * (duidz_ikp1 - duidz_ik);

              const double d2ujdx2f_j = deltaunsurdx * (dujdx_ip1j - dujdx_ij);
              const double d2ujdz2f_j = deltaunsurdz * (dujdz_jkp1 - dujdz_jk);

              const double d2ukdx2f_k = deltaunsurdx * (dukdx_ip1k - dukdx_ik);
              const double d2ukdy2f_k = deltaunsurdy * (dukdy_jp1k - dukdy_jk);

              const double laplacien_uf_i = d2uidx2f_i + d2uidy2f_i + d2uidz2f_i;
              const double laplacien_vf_j = d2ujdx2f_j + d2ujdy2f_j + d2ujdz2f_j;
              const double laplacien_wf_k = d2ukdx2f_k + d2ukdy2f_k + d2ukdz2f_k;

              laplacien_i(i,j,k) = laplacien_uf_i;
              laplacien_j(i,j,k) = laplacien_vf_j;
              laplacien_k(i,j,k) = laplacien_wf_k;
            }
        }
    }
}

void calculer_structural_uu_su_laplacien_u(const double structural_uu_model_constant,
                                           const ArrOfDouble& structural_uu_tensor_coefficients,
                                           const FixedVector<IJK_Field_double, 3>& velocity,
                                           const FixedVector<IJK_Field_double, 3>& laplacien_velocity,
                                           FixedVector<IJK_Field_double, 6>& structural_uu_tensor)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const IJK_Field_double& laplacien_i = laplacien_velocity[0];
  const IJK_Field_double& laplacien_j = laplacien_velocity[1];
  const IJK_Field_double& laplacien_k = laplacien_velocity[2];

  IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const double& coefficient_xx = structural_uu_tensor_coefficients[0];
  const double& coefficient_xy = structural_uu_tensor_coefficients[1];
  const double& coefficient_xz = structural_uu_tensor_coefficients[2];
  const double& coefficient_yy = structural_uu_tensor_coefficients[3];
  const double& coefficient_yz = structural_uu_tensor_coefficients[4];
  const double& coefficient_zz = structural_uu_tensor_coefficients[5];

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;

              const double uf_i = vitesse_i(i,j,k);
              const double uf_ip1 = vitesse_i(i+1,j,k);
              const double uf_i_jm1 = vitesse_i(i,j-1,k);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);

              const double vf_j = vitesse_j(i,j,k);
              const double vf_j_im1 = vitesse_j(i-1,j,k);
              const double vf_jp1 = vitesse_j(i,j+1,k);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);

              const double wf_k = vitesse_k(i,j,k);
              const double wf_k_im1 = vitesse_k(i-1,j,k);
              const double wf_k_jm1 = vitesse_k(i,j-1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);

              const double u_ij = 0.5*(uf_i + uf_i_jm1);
              const double u_ik = 0.5*(uf_i + uf_i_km1);
              const double v_ij = 0.5*(vf_j + vf_j_im1);
              const double v_jk = 0.5*(vf_j + vf_j_km1);
              const double w_ik = 0.5*(wf_k + wf_k_im1);
              const double w_jk = 0.5*(wf_k + wf_k_jm1);

              const double ue = 0.5 * (uf_i + uf_ip1);
              const double ve = 0.5 * (vf_j + vf_jp1);
              const double we = 0.5 * (wf_k + wf_kp1);

              const double laplacien_uf_i = laplacien_i(i,j,k);
              const double laplacien_vf_j = laplacien_j(i,j,k);
              const double laplacien_wf_k = laplacien_k(i,j,k);

              const double laplacien_uf_ip1 = laplacien_i(i+1,j,k);
              const double laplacien_vf_jp1 = laplacien_j(i,j+1,k);
              const double laplacien_wf_kp1 = kg==(nktot-1) ? 0. : laplacien_k(i,j,k+1);

              const double laplacien_uf_i_jm1 = laplacien_i(i,j-1,k);
              const double laplacien_uf_i_km1 = kg==0 ? 0. : laplacien_i(i,j,k-1);
              const double laplacien_vf_j_im1 = laplacien_j(i-1,j,k);
              const double laplacien_vf_j_km1 = kg==0 ? 0. : laplacien_j(i,j,k-1);
              const double laplacien_wf_k_im1 = laplacien_k(i-1,j,k);
              const double laplacien_wf_k_jm1 = laplacien_k(i,j-1,k);

              const double laplacien_u_aij = 0.5*(laplacien_uf_i + laplacien_uf_i_jm1);
              const double laplacien_u_aik = 0.5*(laplacien_uf_i + laplacien_uf_i_km1);
              const double laplacien_v_aij = 0.5*(laplacien_vf_j + laplacien_vf_j_im1);
              const double laplacien_v_ajk = 0.5*(laplacien_vf_j + laplacien_vf_j_km1);
              const double laplacien_w_aik = 0.5*(laplacien_wf_k + laplacien_wf_k_im1);
              const double laplacien_w_ajk = 0.5*(laplacien_wf_k + laplacien_wf_k_jm1);

              const double laplacien_u_e = 0.5 * (laplacien_uf_i + laplacien_uf_ip1);
              const double laplacien_v_e = 0.5 * (laplacien_vf_j + laplacien_vf_jp1);
              const double laplacien_w_e = 0.5 * (laplacien_wf_k + laplacien_wf_kp1);

              const double udu_ii = ue*laplacien_u_e/24.;
              const double udu_ij = u_ij*laplacien_v_aij/24.;
              const double udu_ik = u_ik*laplacien_w_aik/24.;
              const double udu_ji = v_ij*laplacien_u_aij/24.;
              const double udu_jj = ve*laplacien_v_e/24.;
              const double udu_jk = v_jk*laplacien_w_ajk/24.;
              const double udu_ki = w_ik*laplacien_u_aik/24.;
              const double udu_kj = w_jk*laplacien_v_ajk/24.;
              const double udu_kk = we*laplacien_w_e/24.;

              const double su_laplacien_u_ii = - udu_ii - udu_ii;
              const double su_laplacien_u_ij = - udu_ij - udu_ji;
              const double su_laplacien_u_ik = - udu_ik - udu_ki;
              const double su_laplacien_u_jj = - udu_jj - udu_jj;
              const double su_laplacien_u_jk = - udu_jk - udu_kj;
              const double su_laplacien_u_kk = - udu_kk - udu_kk;

              structural_uu_xx(i,j,k) = - coefficient_xx * structural_uu_model_constant * su_laplacien_u_ii;
              structural_uu_xy(i,j,k) = - coefficient_xy * structural_uu_model_constant * su_laplacien_u_ij;
              structural_uu_xz(i,j,k) = - coefficient_xz * structural_uu_model_constant * su_laplacien_u_ik;
              structural_uu_yy(i,j,k) = - coefficient_yy * structural_uu_model_constant * su_laplacien_u_jj;
              structural_uu_yz(i,j,k) = - coefficient_yz * structural_uu_model_constant * su_laplacien_u_jk;
              structural_uu_zz(i,j,k) = - coefficient_zz * structural_uu_model_constant * su_laplacien_u_kk;
            }
        }
    }
}

void calculer_structural_uu_convection(const double structural_uu_model_constant,
                                       const ArrOfDouble& structural_uu_tensor_coefficients,
                                       const FixedVector<IJK_Field_double, 3>& velocity,
                                       FixedVector<IJK_Field_double, 6>& structural_uu_tensor)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const double& coefficient_xx = structural_uu_tensor_coefficients[0];
  const double& coefficient_xy = structural_uu_tensor_coefficients[1];
  const double& coefficient_xz = structural_uu_tensor_coefficients[2];
  const double& coefficient_yy = structural_uu_tensor_coefficients[3];
  const double& coefficient_yz = structural_uu_tensor_coefficients[4];
  const double& coefficient_zz = structural_uu_tensor_coefficients[5];

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;

              const double uf_i = vitesse_i(i,j,k);
              const double uf_ip1 = vitesse_i(i+1,j,k);
              const double uf_i_jm1 = vitesse_i(i,j-1,k);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);

              const double vf_j = vitesse_j(i,j,k);
              const double vf_j_im1 = vitesse_j(i-1,j,k);
              const double vf_jp1 = vitesse_j(i,j+1,k);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);

              const double wf_k = vitesse_k(i,j,k);
              const double wf_k_im1 = vitesse_k(i-1,j,k);
              const double wf_k_jm1 = vitesse_k(i,j-1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);

              const double u_ij = 0.5*(uf_i + uf_i_jm1);
              const double u_ik = 0.5*(uf_i + uf_i_km1);
              const double v_ij = 0.5*(vf_j + vf_j_im1);
              const double v_jk = 0.5*(vf_j + vf_j_km1);
              const double w_ik = 0.5*(wf_k + wf_k_im1);
              const double w_jk = 0.5*(wf_k + wf_k_jm1);

              const double ue = 0.5 * (uf_i + uf_ip1);
              const double ve = 0.5 * (vf_j + vf_jp1);
              const double we = 0.5 * (wf_k + wf_kp1);

              const double uu_ii = ue*ue;
              const double uu_ij = u_ij*v_ij;
              const double uu_ik = u_ik*w_ik;
              const double uu_jj = ve*ve;
              const double uu_jk = v_jk*w_jk;
              const double uu_kk = we*we;

              structural_uu_xx(i,j,k) = - coefficient_xx * structural_uu_model_constant * uu_ii;
              structural_uu_xy(i,j,k) = - coefficient_xy * structural_uu_model_constant * uu_ij;
              structural_uu_xz(i,j,k) = - coefficient_xz * structural_uu_model_constant * uu_ik;
              structural_uu_yy(i,j,k) = - coefficient_yy * structural_uu_model_constant * uu_jj;
              structural_uu_yz(i,j,k) = - coefficient_yz * structural_uu_model_constant * uu_jk;
              structural_uu_zz(i,j,k) = - coefficient_zz * structural_uu_model_constant * uu_kk;
            }
        }
    }
}
template<class T>
void calculer_structural_uu_similarity_comp(const double structural_uu_model_constant,
                                            const ArrOfDouble& structural_uu_tensor_coefficients,
                                            const IJK_Field_double& rho,
                                            const FixedVector<IJK_Field_double, 3>& velocity,
                                            const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                            const ArrOfDouble_with_ghost& delta_z,
                                            const double facteur_delta_x,
                                            const double facteur_delta_y,
                                            const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                            T& kernel,
                                            FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                            FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                            FixedVector<IJK_Field_double, 6>& structural_uu_tensor)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];
  const IJK_Field_double& masse_vol_ijk = rho;
  IJK_Field_local_double& masse_vol_ijk_filtre = tmp_b[15];
  IJK_Field_local_double&   aa_rho_ijk_filtre = tmp_a[15];

  // const IJK_Field_double& vitesse_i_filtre = velocity_filtre[0];
  // const IJK_Field_double& vitesse_j_filtre = velocity_filtre[1];
  // const IJK_Field_double& vitesse_k_filtre = velocity_filtre[2];

  IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const double& coefficient_xx = structural_uu_tensor_coefficients[0];
  const double& coefficient_xy = structural_uu_tensor_coefficients[1];
  const double& coefficient_xz = structural_uu_tensor_coefficients[2];
  const double& coefficient_yy = structural_uu_tensor_coefficients[3];
  const double& coefficient_yz = structural_uu_tensor_coefficients[4];
  const double& coefficient_zz = structural_uu_tensor_coefficients[5];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();
  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& b_ii = tmp_b[0];
  IJK_Field_local_double& b_ij = tmp_b[1];
  IJK_Field_local_double& b_ik = tmp_b[2];
  IJK_Field_local_double& b_jj = tmp_b[3];
  IJK_Field_local_double& b_jk = tmp_b[4];
  IJK_Field_local_double& b_kk = tmp_b[5];
  IJK_Field_local_double& bb_uii = tmp_b[6];
  IJK_Field_local_double& bb_uij = tmp_b[7];
  IJK_Field_local_double& bb_vjj = tmp_b[8];
  IJK_Field_local_double& bb_vij = tmp_b[9];
  IJK_Field_local_double& bb_wkk = tmp_b[10];
  IJK_Field_local_double& bb_uik = tmp_b[11];
  IJK_Field_local_double& bb_vjk = tmp_b[12];
  IJK_Field_local_double& bb_wik = tmp_b[13];
  IJK_Field_local_double& bb_wjk = tmp_b[14];

  IJK_Field_local_double& a_ii = tmp_a[0];
  IJK_Field_local_double& a_ij = tmp_a[1];
  IJK_Field_local_double& a_ik = tmp_a[2];
  IJK_Field_local_double& a_jj = tmp_a[3];
  IJK_Field_local_double& a_jk = tmp_a[4];
  IJK_Field_local_double& a_kk = tmp_a[5];
  IJK_Field_local_double& aa_uii = tmp_a[6];
  IJK_Field_local_double& aa_uij = tmp_a[7];
  IJK_Field_local_double& aa_vjj = tmp_a[8];
  IJK_Field_local_double& aa_vij = tmp_a[9];
  IJK_Field_local_double& aa_wkk = tmp_a[10];
  IJK_Field_local_double& aa_uik = tmp_a[11];
  IJK_Field_local_double& aa_vjk = tmp_a[12];
  IJK_Field_local_double& aa_wik = tmp_a[13];
  IJK_Field_local_double& aa_wjk = tmp_a[14];

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_m1_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z[k-1];
      const double delta_m_glo = kg==0 ? 0.5*dz_glo : 0.5*(dz_glo + dz_m1_glo);
      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];
      const double dz_m1_pour_delta_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z_pour_delta[k-1];
      const double delta_m_pour_delta_glo = (kg-1<0 || kg>(nktot-1)) ? 0. : 0.5*(dz_pour_delta_glo + dz_m1_pour_delta_glo);

      const FixedVector<double, 21> filter_kernel_z = kernel->inhomogeneous(true, k, kg, nktot, dz_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_z_face = kernel->inhomogeneous(false, k, kg, nktot, delta_m_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_elem = kernel->size_k_elem(kg, nktot);
      const int size_k_face = kernel->size_k_face(kg, nktot);
      const int shift_k_elem = kernel->shift_k_elem(kg);
      const int shift_k_face = kernel->shift_k_face(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_elem = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_elem = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<-1 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double filter_coef_z = filter_kernel_z[kp+10];
                  longueur_elem += filter_coef_z * dz;
                }
              facteur_elem = 1./longueur_elem;
            }
          else
            {
              facteur_elem = dz_glo==0. ? 0. : 1./dz_glo;
            }
        }
      double facteur_face = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_face = 0.;
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<0 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                  const double dzf = 0.5*(dz + dz_m1);
                  const double filter_coef_z_face = filter_kernel_z_face[kp+10];
                  longueur_face += filter_coef_z_face * dzf;
                }
              facteur_face = 1./longueur_face;
            }
          else
            {
              facteur_face = delta_m_glo==0. ? 0. : 1./delta_m_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              b_ii(i, j, 0) = 0.;
              bb_uii(i, j, 0) = 0.;
              bb_uij(i, j, 0) = 0.;
              b_ij(i, j, 0) = 0.;
              b_ik(i, j, 0) = 0.;
              b_jj(i, j, 0) = 0.;
              bb_vjj(i, j, 0) = 0.;
              bb_vij(i, j, 0) = 0.;
              b_jk(i, j, 0) = 0.;
              b_kk(i, j, 0) = 0.;
              bb_wkk(i, j, 0) = 0.;
              bb_uik(i, j, 0) = 0.;
              bb_vjk(i, j, 0) = 0.;
              bb_wik(i, j, 0) = 0.;
              bb_wjk(i, j, 0) = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_elem(kpg, nktot)))
                    {
                      const double filter_coef_z = filter_kernel_z[kp+10];

                      const double rho_ijk = masse_vol_ijk(i,j,k+kp);
                      const double uf_i = vitesse_i(i,j,k+kp);
                      const double vf_j = vitesse_j(i,j,k+kp);
                      const double wf_k = vitesse_k(i,j,k+kp);

                      const double uf_ip1 = vitesse_i(i+1,j,k+kp);
                      const double vf_jp1 = vitesse_j(i,j+1,k+kp);
                      const double wf_kp1 = kpg==(nktot-1) ? 0. : vitesse_k(i,j,k+1+kp);

                      const double uf_i_jm1 = vitesse_i(i,j-1,k+kp);
                      const double vf_j_im1 = vitesse_j(i-1,j,k+kp);

                      const double u_ij = 0.5*(uf_i + uf_i_jm1);
                      const double v_ij = 0.5*(vf_j + vf_j_im1);

                      const double ue = 0.5 * (uf_i + uf_ip1);
                      const double ve = 0.5 * (vf_j + vf_jp1);
                      const double we = 0.5 * (wf_k + wf_kp1);

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          b_ii(i, j, 0) += rho_ijk*ue*ue     * filter_coef_z      * dz  * facteur_elem;
                          b_ij(i, j, 0) += rho_ijk*u_ij*v_ij * filter_coef_z      * dz  * facteur_elem;
                          b_jj(i, j, 0) += rho_ijk*ve*ve     * filter_coef_z      * dz  * facteur_elem;
                          b_kk(i, j, 0) += rho_ijk*we*we     * filter_coef_z      * dz  * facteur_elem;
                          bb_uii(i, j, 0) += rho_ijk*ue     * filter_coef_z      * dz  * facteur_elem;
                          bb_uij(i, j, 0) += rho_ijk*u_ij     * filter_coef_z      * dz  * facteur_elem;
                          bb_vij(i, j, 0) += rho_ijk*v_ij     * filter_coef_z      * dz  * facteur_elem;
                          bb_vjj(i, j, 0) += rho_ijk*ve     * filter_coef_z      * dz  * facteur_elem;
                          bb_wkk(i, j, 0) += rho_ijk*we     * filter_coef_z      * dz  * facteur_elem;
                          masse_vol_ijk_filtre(i, j, 0) += rho_ijk     * filter_coef_z      * dz  * facteur_elem;
                        }
                      else
                        {
                          b_ii(i, j, 0) += rho_ijk*ue*ue     * filter_coef_z;
                          b_ij(i, j, 0) += rho_ijk*u_ij*v_ij * filter_coef_z;
                          b_jj(i, j, 0) += rho_ijk*ve*ve     * filter_coef_z;
                          b_kk(i, j, 0) += rho_ijk*we*we     * filter_coef_z;
                          bb_uii(i, j, 0) += rho_ijk*ue     * filter_coef_z;
                          bb_uij(i, j, 0) += rho_ijk*u_ij     * filter_coef_z;
                          bb_vij(i, j, 0) += rho_ijk*v_ij     * filter_coef_z;
                          bb_vjj(i, j, 0) += rho_ijk*ve     * filter_coef_z;
                          bb_wkk(i, j, 0) += rho_ijk*we     * filter_coef_z;
                          masse_vol_ijk_filtre(i, j, 0) += rho_ijk     * filter_coef_z;
                        }
                    }
                }
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_face(kpg, nktot)))
                    {
                      const double filter_coef_z_face = filter_kernel_z_face[kp+10];

                      const double rho_ijk = masse_vol_ijk(i,j,k+kp);
                      const double uf_i = vitesse_i(i,j,k+kp);
                      const double vf_j = vitesse_j(i,j,k+kp);
                      const double wf_k = vitesse_k(i,j,k+kp);

                      const double uf_i_km1 = kpg==0 ? 0. : vitesse_i(i,j,k-1+kp);
                      const double vf_j_km1 = kpg==0 ? 0. : vitesse_j(i,j,k-1+kp);
                      const double wf_k_im1 = vitesse_k(i-1,j,k+kp);
                      const double wf_k_jm1 = vitesse_k(i,j-1,k+kp);

                      const double u_ik = 0.5*(uf_i + uf_i_km1);
                      const double v_jk = 0.5*(vf_j + vf_j_km1);
                      const double w_ik = 0.5*(wf_k + wf_k_im1);
                      const double w_jk = 0.5*(wf_k + wf_k_jm1);

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                          const double dzf = 0.5*(dz + dz_m1);
                          b_ik(i, j, 0) += rho_ijk*u_ik*w_ik * filter_coef_z_face * dzf * facteur_face;
                          b_jk(i, j, 0) += rho_ijk*v_jk*w_jk * filter_coef_z_face * dzf * facteur_face;
                          bb_uik(i, j, 0) += rho_ijk*u_ik * filter_coef_z_face * dzf * facteur_face;
                          bb_vjk(i, j, 0) += rho_ijk*v_jk * filter_coef_z_face * dzf * facteur_face;
                          bb_wik(i, j, 0) += rho_ijk*w_ik * filter_coef_z_face * dzf * facteur_face;
                          bb_wjk(i, j, 0) += rho_ijk*w_jk * filter_coef_z_face * dzf * facteur_face;
                        }
                      else
                        {
                          b_ik(i, j, 0) += rho_ijk*u_ik*w_ik * filter_coef_z_face;
                          b_jk(i, j, 0) += rho_ijk*v_jk*w_jk * filter_coef_z_face;
                          bb_uik(i, j, 0) += rho_ijk*u_ik * filter_coef_z_face;
                          bb_vjk(i, j, 0) += rho_ijk*v_jk * filter_coef_z_face;
                          bb_wik(i, j, 0) += rho_ijk*w_ik * filter_coef_z_face;
                          bb_wjk(i, j, 0) += rho_ijk*w_jk * filter_coef_z_face;
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              a_ii(i, 0, 0) = 0.;
              aa_uii(i, 0, 0) = 0.;
              aa_uij(i, 0, 0) = 0.;
              a_ij(i, 0, 0) = 0.;
              a_ik(i, 0, 0) = 0.;
              a_jj(i, 0, 0) = 0.;
              aa_vjj(i, 0, 0) = 0.;
              aa_vij(i, 0, 0) = 0.;
              a_jk(i, 0, 0) = 0.;
              a_kk(i, 0, 0) = 0.;
              aa_wkk(i, 0, 0) = 0.;
              aa_uik(i, 0, 0) = 0.;
              aa_vjk(i, 0, 0) = 0.;
              aa_wik(i, 0, 0) = 0.;
              aa_wjk(i, 0, 0) = 0.;
              aa_rho_ijk_filtre(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  a_ii(i, 0, 0) += b_ii(i, j+jp, 0) * filter_coef_y;
                  aa_uii(i, 0, 0) += bb_uii(i, j+jp, 0) * filter_coef_y;
                  aa_uij(i, 0, 0) += bb_uij(i, j+jp, 0) * filter_coef_y;
                  a_ij(i, 0, 0) += b_ij(i, j+jp, 0) * filter_coef_y;
                  a_ik(i, 0, 0) += b_ik(i, j+jp, 0) * filter_coef_y;
                  a_jj(i, 0, 0) += b_jj(i, j+jp, 0) * filter_coef_y;
                  aa_vjj(i, 0, 0) += bb_vjj(i, j+jp, 0) * filter_coef_y;
                  aa_vij(i, 0, 0) += bb_vij(i, j+jp, 0) * filter_coef_y;
                  a_jk(i, 0, 0) += b_jk(i, j+jp, 0) * filter_coef_y;
                  a_kk(i, 0, 0) += b_kk(i, j+jp, 0) * filter_coef_y;
                  aa_wkk(i, 0, 0) += bb_wkk(i, j+jp, 0) * filter_coef_y;
                  aa_uik(i, 0, 0) += bb_uik(i, j+jp, 0) * filter_coef_y;
                  aa_vjk(i, 0, 0) += bb_vjk(i, j+jp, 0) * filter_coef_y;
                  aa_wik(i, 0, 0) += bb_wik(i, j+jp, 0) * filter_coef_y;
                  aa_wjk(i, 0, 0) += bb_wjk(i, j+jp, 0) * filter_coef_y;
                  aa_rho_ijk_filtre(i, 0, 0) += masse_vol_ijk_filtre(i, j+jp, 0)* filter_coef_y;
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double r_ii = 0.;
              double r_ij = 0.;
              double r_ik = 0.;
              double r_jj = 0.;
              double r_jk = 0.;
              double r_kk = 0.;
              double rho_ue = 0.;
              double rho_ve = 0.;
              double rho_we = 0.;
              double rho_uij = 0.;
              double rho_uik = 0.;
              double rho_vij = 0.;
              double rho_vjk = 0.;
              double rho_wik = 0.;
              double rho_wjk = 0.;
              double rho_filtre = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  r_ii += a_ii(i+ip, 0, 0) * filter_coef_x;
                  r_ij += a_ij(i+ip, 0, 0) * filter_coef_x;
                  r_ik += a_ik(i+ip, 0, 0) * filter_coef_x;
                  r_jj += a_jj(i+ip, 0, 0) * filter_coef_x;
                  r_jk += a_jk(i+ip, 0, 0) * filter_coef_x;
                  r_kk += a_kk(i+ip, 0, 0) * filter_coef_x;
                  rho_ue += aa_uii(i+ip, 0, 0) * filter_coef_x;
                  rho_ve += aa_vjj(i+ip, 0, 0) * filter_coef_x;
                  rho_we += aa_wkk(i+ip, 0, 0) * filter_coef_x;
                  rho_uij += aa_uij(i+ip, 0, 0) * filter_coef_x;
                  rho_uik += aa_uik(i+ip, 0, 0) * filter_coef_x;
                  rho_vij += aa_vij(i+ip, 0, 0) * filter_coef_x;
                  rho_vjk += aa_vjk(i+ip, 0, 0) * filter_coef_x;
                  rho_wik += aa_wik(i+ip, 0, 0) * filter_coef_x;
                  rho_wjk += aa_wjk(i+ip, 0, 0) * filter_coef_x;
                  rho_filtre += aa_rho_ijk_filtre(i+ip, 0, 0)* filter_coef_x;
                }

              const double c_ii = (r_ii)/rho_filtre - (rho_ue*rho_ue)/(rho_filtre*rho_filtre);
              const double c_ij = (r_ij)/rho_filtre - (rho_uij*rho_vij)/(rho_filtre*rho_filtre);
              const double c_ik = (r_ik)/rho_filtre - (rho_uik*rho_wik)/(rho_filtre*rho_filtre);
              const double c_jj = (r_jj)/rho_filtre - (rho_ve*rho_ve)/(rho_filtre*rho_filtre);
              const double c_jk = (r_jk)/rho_filtre - (rho_vjk*rho_wjk)/(rho_filtre*rho_filtre);
              const double c_kk = (r_kk)/rho_filtre - (rho_we*rho_we)/(rho_filtre*rho_filtre);

              structural_uu_xx(i,j,k) = - coefficient_xx * structural_uu_model_constant * c_ii;
              structural_uu_xy(i,j,k) = - coefficient_xy * structural_uu_model_constant * c_ij;
              structural_uu_xz(i,j,k) = - coefficient_xz * structural_uu_model_constant * c_ik;
              structural_uu_yy(i,j,k) = - coefficient_yy * structural_uu_model_constant * c_jj;
              structural_uu_yz(i,j,k) = - coefficient_yz * structural_uu_model_constant * c_jk;
              structural_uu_zz(i,j,k) = - coefficient_zz * structural_uu_model_constant * c_kk;
            }
        }
    }
}

template<class T>
void calculer_structural_uu_similarity(const double structural_uu_model_constant,
                                       const ArrOfDouble& structural_uu_tensor_coefficients,
                                       const FixedVector<IJK_Field_double, 3>& velocity,
                                       const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                       const ArrOfDouble_with_ghost& delta_z,
                                       const double facteur_delta_x,
                                       const double facteur_delta_y,
                                       const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                       T& kernel,
                                       FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                       FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                       FixedVector<IJK_Field_double, 6>& structural_uu_tensor)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const IJK_Field_double& vitesse_i_filtre = velocity_filtre[0];
  const IJK_Field_double& vitesse_j_filtre = velocity_filtre[1];
  const IJK_Field_double& vitesse_k_filtre = velocity_filtre[2];

  IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const double& coefficient_xx = structural_uu_tensor_coefficients[0];
  const double& coefficient_xy = structural_uu_tensor_coefficients[1];
  const double& coefficient_xz = structural_uu_tensor_coefficients[2];
  const double& coefficient_yy = structural_uu_tensor_coefficients[3];
  const double& coefficient_yz = structural_uu_tensor_coefficients[4];
  const double& coefficient_zz = structural_uu_tensor_coefficients[5];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& b_ii = tmp_b[0];
  IJK_Field_local_double& b_ij = tmp_b[1];
  IJK_Field_local_double& b_ik = tmp_b[2];
  IJK_Field_local_double& b_jj = tmp_b[3];
  IJK_Field_local_double& b_jk = tmp_b[4];
  IJK_Field_local_double& b_kk = tmp_b[5];

  IJK_Field_local_double& a_ii = tmp_a[0];
  IJK_Field_local_double& a_ij = tmp_a[1];
  IJK_Field_local_double& a_ik = tmp_a[2];
  IJK_Field_local_double& a_jj = tmp_a[3];
  IJK_Field_local_double& a_jk = tmp_a[4];
  IJK_Field_local_double& a_kk = tmp_a[5];

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_m1_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z[k-1];
      const double delta_m_glo = kg==0 ? 0.5*dz_glo : 0.5*(dz_glo + dz_m1_glo);
      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];
      const double dz_m1_pour_delta_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z_pour_delta[k-1];
      const double delta_m_pour_delta_glo = (kg-1<0 || kg>(nktot-1)) ? 0. : 0.5*(dz_pour_delta_glo + dz_m1_pour_delta_glo);

      const FixedVector<double, 21> filter_kernel_z = kernel->inhomogeneous(true, k, kg, nktot, dz_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_z_face = kernel->inhomogeneous(false, k, kg, nktot, delta_m_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_elem = kernel->size_k_elem(kg, nktot);
      const int size_k_face = kernel->size_k_face(kg, nktot);
      const int shift_k_elem = kernel->shift_k_elem(kg);
      const int shift_k_face = kernel->shift_k_face(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_elem = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_elem = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<-1 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double filter_coef_z = filter_kernel_z[kp+10];
                  longueur_elem += filter_coef_z * dz;
                }
              facteur_elem = 1./longueur_elem;
            }
          else
            {
              facteur_elem = dz_glo==0. ? 0. : 1./dz_glo;
            }
        }
      double facteur_face = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_face = 0.;
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<0 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                  const double dzf = 0.5*(dz + dz_m1);
                  const double filter_coef_z_face = filter_kernel_z_face[kp+10];
                  longueur_face += filter_coef_z_face * dzf;
                }
              facteur_face = 1./longueur_face;
            }
          else
            {
              facteur_face = delta_m_glo==0. ? 0. : 1./delta_m_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              b_ii(i, j, 0) = 0.;
              b_ij(i, j, 0) = 0.;
              b_ik(i, j, 0) = 0.;
              b_jj(i, j, 0) = 0.;
              b_jk(i, j, 0) = 0.;
              b_kk(i, j, 0) = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_elem(kpg, nktot)))
                    {
                      const double filter_coef_z = filter_kernel_z[kp+10];

                      const double uf_i = vitesse_i(i,j,k+kp);
                      const double vf_j = vitesse_j(i,j,k+kp);
                      const double wf_k = vitesse_k(i,j,k+kp);

                      const double uf_ip1 = vitesse_i(i+1,j,k+kp);
                      const double vf_jp1 = vitesse_j(i,j+1,k+kp);
                      const double wf_kp1 = kpg==(nktot-1) ? 0. : vitesse_k(i,j,k+1+kp);

                      const double uf_i_jm1 = vitesse_i(i,j-1,k+kp);
                      const double vf_j_im1 = vitesse_j(i-1,j,k+kp);

                      const double u_ij = 0.5*(uf_i + uf_i_jm1);
                      const double v_ij = 0.5*(vf_j + vf_j_im1);

                      const double ue = 0.5 * (uf_i + uf_ip1);
                      const double ve = 0.5 * (vf_j + vf_jp1);
                      const double we = 0.5 * (wf_k + wf_kp1);

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          b_ii(i, j, 0) += ue*ue     * filter_coef_z      * dz  * facteur_elem;
                          b_ij(i, j, 0) += u_ij*v_ij * filter_coef_z      * dz  * facteur_elem;
                          b_jj(i, j, 0) += ve*ve     * filter_coef_z      * dz  * facteur_elem;
                          b_kk(i, j, 0) += we*we     * filter_coef_z      * dz  * facteur_elem;
                        }
                      else
                        {
                          b_ii(i, j, 0) += ue*ue     * filter_coef_z;
                          b_ij(i, j, 0) += u_ij*v_ij * filter_coef_z;
                          b_jj(i, j, 0) += ve*ve     * filter_coef_z;
                          b_kk(i, j, 0) += we*we     * filter_coef_z;
                        }
                    }
                }
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_face(kpg, nktot)))
                    {
                      const double filter_coef_z_face = filter_kernel_z_face[kp+10];

                      const double uf_i = vitesse_i(i,j,k+kp);
                      const double vf_j = vitesse_j(i,j,k+kp);
                      const double wf_k = vitesse_k(i,j,k+kp);

                      const double uf_i_km1 = kpg==0 ? 0. : vitesse_i(i,j,k-1+kp);
                      const double vf_j_km1 = kpg==0 ? 0. : vitesse_j(i,j,k-1+kp);
                      const double wf_k_im1 = vitesse_k(i-1,j,k+kp);
                      const double wf_k_jm1 = vitesse_k(i,j-1,k+kp);

                      const double u_ik = 0.5*(uf_i + uf_i_km1);
                      const double v_jk = 0.5*(vf_j + vf_j_km1);
                      const double w_ik = 0.5*(wf_k + wf_k_im1);
                      const double w_jk = 0.5*(wf_k + wf_k_jm1);

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                          const double dzf = 0.5*(dz + dz_m1);
                          b_ik(i, j, 0) += u_ik*w_ik * filter_coef_z_face * dzf * facteur_face;
                          b_jk(i, j, 0) += v_jk*w_jk * filter_coef_z_face * dzf * facteur_face;
                        }
                      else
                        {
                          b_ik(i, j, 0) += u_ik*w_ik * filter_coef_z_face;
                          b_jk(i, j, 0) += v_jk*w_jk * filter_coef_z_face;
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              a_ii(i, 0, 0) = 0.;
              a_ij(i, 0, 0) = 0.;
              a_ik(i, 0, 0) = 0.;
              a_jj(i, 0, 0) = 0.;
              a_jk(i, 0, 0) = 0.;
              a_kk(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  a_ii(i, 0, 0) += b_ii(i, j+jp, 0) * filter_coef_y;
                  a_ij(i, 0, 0) += b_ij(i, j+jp, 0) * filter_coef_y;
                  a_ik(i, 0, 0) += b_ik(i, j+jp, 0) * filter_coef_y;
                  a_jj(i, 0, 0) += b_jj(i, j+jp, 0) * filter_coef_y;
                  a_jk(i, 0, 0) += b_jk(i, j+jp, 0) * filter_coef_y;
                  a_kk(i, 0, 0) += b_kk(i, j+jp, 0) * filter_coef_y;
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double r_ii = 0.;
              double r_ij = 0.;
              double r_ik = 0.;
              double r_jj = 0.;
              double r_jk = 0.;
              double r_kk = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  r_ii += a_ii(i+ip, 0, 0) * filter_coef_x;
                  r_ij += a_ij(i+ip, 0, 0) * filter_coef_x;
                  r_ik += a_ik(i+ip, 0, 0) * filter_coef_x;
                  r_jj += a_jj(i+ip, 0, 0) * filter_coef_x;
                  r_jk += a_jk(i+ip, 0, 0) * filter_coef_x;
                  r_kk += a_kk(i+ip, 0, 0) * filter_coef_x;
                }

              const double uf_i = vitesse_i_filtre(i,j,k);
              const double vf_j = vitesse_j_filtre(i,j,k);
              const double wf_k = vitesse_k_filtre(i,j,k);

              const double uf_ip1 = vitesse_i_filtre(i+1,j,k);
              const double vf_jp1 = vitesse_j_filtre(i,j+1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k_filtre(i,j,k+1);

              const double uf_i_jm1 = vitesse_i_filtre(i,j-1,k);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i_filtre(i,j,k-1);
              const double vf_j_im1 = vitesse_j_filtre(i-1,j,k);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j_filtre(i,j,k-1);
              const double wf_k_im1 = vitesse_k_filtre(i-1,j,k);
              const double wf_k_jm1 = vitesse_k_filtre(i,j-1,k);

              const double u_ij = 0.5*(uf_i + uf_i_jm1);
              const double u_ik = 0.5*(uf_i + uf_i_km1);
              const double v_ij = 0.5*(vf_j + vf_j_im1);
              const double v_jk = 0.5*(vf_j + vf_j_km1);
              const double w_ik = 0.5*(wf_k + wf_k_im1);
              const double w_jk = 0.5*(wf_k + wf_k_jm1);

              const double ue = 0.5 * (uf_i + uf_ip1);
              const double ve = 0.5 * (vf_j + vf_jp1);
              const double we = 0.5 * (wf_k + wf_kp1);

              const double c_ii = r_ii - ue*ue;
              const double c_ij = r_ij - u_ij*v_ij;
              const double c_ik = r_ik - u_ik*w_ik;
              const double c_jj = r_jj - ve*ve;
              const double c_jk = r_jk - v_jk*w_jk;
              const double c_kk = r_kk - we*we;

              structural_uu_xx(i,j,k) = - coefficient_xx * structural_uu_model_constant * c_ii;
              structural_uu_xy(i,j,k) = - coefficient_xy * structural_uu_model_constant * c_ij;
              structural_uu_xz(i,j,k) = - coefficient_xz * structural_uu_model_constant * c_ik;
              structural_uu_yy(i,j,k) = - coefficient_yy * structural_uu_model_constant * c_jj;
              structural_uu_yz(i,j,k) = - coefficient_yz * structural_uu_model_constant * c_jk;
              structural_uu_zz(i,j,k) = - coefficient_zz * structural_uu_model_constant * c_kk;
            }
        }
    }
}


//#####################
//#####Martin 092021###
//#####################

template<class T>
void calculer_structural_uu_similarity_Streher(const double structural_uu_model_constant,
                                               const FixedVector<IJK_Field_double, 3>& velocity,
                                               T& kernel,
                                               FixedVector<IJK_Field_double, 6>& structural_uu_tensor)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const double uf_im1 = vitesse_i(i-1,j,k);
              const double uf_i = vitesse_i(i,j,k);
              const double uf_i_jm2 = vitesse_i(i,j-2,k);
              const double uf_i_jm1 = vitesse_i(i,j-1,k);
              const double uf_i_jp1 = vitesse_i(i,j+1,k);
              const double uf_i_km2 = kg<=1 ? 0. : vitesse_i(i,j,k-2);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);
              const double uf_i_kp1 = kg>=(nktot-1) ? 0. : vitesse_i(i,j,k+1);
              const double uf_ip1 = vitesse_i(i+1,j,k);
              const double uf_ip2 = vitesse_i(i+2,j,k);

              const double vf_jm1 = vitesse_j(i,j-1,k);
              const double vf_j = vitesse_j(i,j,k);
              const double vf_j_im2 = vitesse_j(i-2,j,k);
              const double vf_j_im1 = vitesse_j(i-1,j,k);
              const double vf_j_ip1 = vitesse_j(i+1,j,k);
              const double vf_j_km2 = kg<=1 ? 0. : vitesse_j(i,j,k-2);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);
              const double vf_j_kp1 = kg>=(nktot-1) ? 0. : vitesse_j(i,j,k+1);
              const double vf_jp1 = vitesse_j(i,j+1,k);
              const double vf_jp2 = vitesse_j(i,j+2,k);

              const double wf_km1 = kg<=1 ? 0. : vitesse_k(i,j,k-1);
              const double wf_k = vitesse_k(i,j,k);
              const double wf_k_im2 = vitesse_k(i-2,j,k);
              const double wf_k_im1 = vitesse_k(i-1,j,k);
              const double wf_k_ip1 = vitesse_k(i+1,j,k);
              const double wf_k_jm2 = vitesse_k(i,j-2,k);
              const double wf_k_jm1 = vitesse_k(i,j-1,k);
              const double wf_k_jp1 = vitesse_k(i,j+1,k);
              const double wf_kp1 = kg>=(nktot-1) ? 0. : vitesse_k(i,j,k+1);
              const double wf_kp2 = kg>=(nktot-2) ? 0. : vitesse_k(i,j,k+2);

              const double c_ii = (uf_i+uf_ip1)*(uf_i+uf_ip1)/4.-(uf_im1+uf_i+uf_ip1+uf_ip2)*(uf_im1+uf_i+uf_ip1+uf_ip2)/16.;
              const double c_jj = (vf_j+vf_jp1)*(vf_j+vf_jp1)/4.-(vf_jm1+vf_j+vf_jp1+vf_jp2)*(vf_jm1+vf_j+vf_jp1+vf_jp2)/16.;
              const double c_kk = (wf_k+wf_kp1)*(wf_k+wf_kp1)/4.-(wf_km1+wf_k+wf_kp1+wf_kp2)*(wf_km1+wf_k+wf_kp1+wf_kp2)/16.;
              const double c_ij = (uf_i_jm1+uf_i)*(vf_j_im1+vf_j)/4.-(uf_i_jp1+uf_i+uf_i_jm1+uf_i_jm2)*(vf_j_ip1+vf_j+vf_j_im1+vf_j_im2)/16.;
              const double c_ik = (uf_i_km1+uf_i)*(wf_k_im1+wf_k)/4.-(uf_i_kp1+uf_i+uf_i_km1+uf_i_km2)*(wf_k_ip1+wf_k+wf_k_im1+wf_k_im2)/16.;
              const double c_jk = (vf_j_km1+vf_j)*(wf_k_jm1+wf_k)/4.-(vf_j_kp1+vf_j+vf_j_km1+vf_j_km2)*(wf_k_jp1+wf_k+wf_k_jm1+wf_k_jm2)/16.;

              structural_uu_xx(i,j,k) = - c_ii * structural_uu_model_constant;
              structural_uu_xy(i,j,k) = - c_ij * structural_uu_model_constant;
              structural_uu_xz(i,j,k) = - c_ik * structural_uu_model_constant;
              structural_uu_yy(i,j,k) = - c_jj * structural_uu_model_constant;
              structural_uu_yz(i,j,k) = - c_jk * structural_uu_model_constant;
              structural_uu_zz(i,j,k) = - c_kk * structural_uu_model_constant;
            }
        }
    }
}

template<class T>
void calculer_structural_uu(const Nom& structural_uu_model,
                            const double structural_uu_model_constant,
                            const ArrOfDouble& structural_uu_tensor_coefficients,
                            IJK_Field_double& rho,
                            FixedVector<IJK_Field_double, 3>& velocity,
                            FixedVector<IJK_Field_double, 3>& velocity_filtre,
                            const ArrOfDouble_with_ghost& delta_z,
                            const double facteur_delta_x,
                            const double facteur_delta_y,
                            const ArrOfDouble_with_ghost& delta_z_pour_delta,
                            const double facteur_delta_filtre_x,
                            const double facteur_delta_filtre_y,
                            const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                            T& kernel,
                            FixedVector<IJK_Field_local_double, 18>& tmp_b,
                            FixedVector<IJK_Field_local_double, 18>& tmp_a,
                            FixedVector<IJK_Field_double, 6>& structural_uu_tmp_tensor,
                            const bool flag_structural_uu_filtre,
                            FixedVector<IJK_Field_double, 6>& structural_uu_filtre_tensor,
                            FixedVector<IJK_Field_double, 6>& structural_uu_tensor)
{
  int ghost_size_filter;
  int ghost_size_velocity;
  int ghost_size_structural_uu_tmp;

  if ( structural_uu_model == Nom("gradient") )
    {
      calculer_structural_uu_gradient(structural_uu_model_constant,
                                      structural_uu_tensor_coefficients,
                                      velocity,
                                      delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                      structural_uu_tensor);
    }
  else if ( structural_uu_model == Nom("gradient_filtre") )
    {
      calculer_structural_uu_gradient(structural_uu_model_constant,
                                      structural_uu_tensor_coefficients,
                                      velocity,
                                      delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                      structural_uu_tmp_tensor);

      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_structural_uu_tmp = max((int)2, ghost_size_filter);
      structural_uu_tmp_tensor[0].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[1].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[2].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[3].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[4].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[5].echange_espace_virtuel(ghost_size_structural_uu_tmp);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[0]);
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[1]);
      filtrer_champ_face(flag_add, structural_uu_tmp_tensor[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[2]);
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[3], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[3]);
      filtrer_champ_face(flag_add, structural_uu_tmp_tensor[4], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[4]);
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[5], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[5]);
    }
  else if ( structural_uu_model == Nom("su_laplacien_u") )
    {
      velocity[0].echange_espace_virtuel(2);
      velocity[1].echange_espace_virtuel(2);
      velocity[2].echange_espace_virtuel(2);

      calculer_laplacien_u(velocity,
                           delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                           velocity_filtre);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);

      calculer_structural_uu_su_laplacien_u(structural_uu_model_constant,
                                            structural_uu_tensor_coefficients,
                                            velocity,
                                            velocity_filtre,
                                            structural_uu_tensor);
    }
  else if ( structural_uu_model == Nom("convection") )
    {
      calculer_structural_uu_convection(structural_uu_model_constant,
                                        structural_uu_tensor_coefficients,
                                        velocity,
                                        structural_uu_tensor);
    }
  else if ( structural_uu_model == Nom("convection_filtre") )
    {
      calculer_structural_uu_convection(structural_uu_model_constant,
                                        structural_uu_tensor_coefficients,
                                        velocity,
                                        structural_uu_tmp_tensor);

      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_structural_uu_tmp = max((int)2, ghost_size_filter);
      structural_uu_tmp_tensor[0].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[1].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[2].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[3].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[4].echange_espace_virtuel(ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor[5].echange_espace_virtuel(ghost_size_structural_uu_tmp);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[0]);
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[1]);
      filtrer_champ_face(flag_add, structural_uu_tmp_tensor[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[2]);
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[3], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[3]);
      filtrer_champ_face(flag_add, structural_uu_tmp_tensor[4], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[4]);
      filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[5], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_tensor[5]);
    }
  else if ( structural_uu_model == Nom("similarity") )
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);

      calculer_structural_uu_similarity(structural_uu_model_constant,
                                        structural_uu_tensor_coefficients,
                                        velocity, velocity_filtre,
                                        delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                        kernel, tmp_b, tmp_a,
                                        structural_uu_tensor);
    }
  else if ( structural_uu_model == Nom("similarity_comp") )
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);

      calculer_structural_uu_similarity_comp(structural_uu_model_constant,
                                             structural_uu_tensor_coefficients,
                                             rho,
                                             velocity, velocity_filtre,
                                             delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                             kernel, tmp_b, tmp_a,
                                             structural_uu_tensor);
    }
  else if ( structural_uu_model == Nom("similarity_streher") )
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);

      calculer_structural_uu_similarity_Streher(structural_uu_model_constant,
                                                velocity,
                                                kernel,
                                                structural_uu_tensor);
    }
  else
    {
      Cerr << "The name of the structural model for uu is unknown." << finl;
      Process::exit();
    }

  if (flag_structural_uu_filtre)
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);


      if ( structural_uu_model == Nom("gradient") )
        {
          calculer_structural_uu_gradient(structural_uu_model_constant,
                                          structural_uu_tensor_coefficients,
                                          velocity_filtre,
                                          delta_z, facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                          structural_uu_filtre_tensor);
        }
      else if ( structural_uu_model == Nom("gradient_filtre") )
        {
          calculer_structural_uu_gradient(structural_uu_model_constant,
                                          structural_uu_tensor_coefficients,
                                          velocity_filtre,
                                          delta_z, facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                          structural_uu_tmp_tensor);

          ghost_size_filter = 1 + kernel->ghost_size();
          ghost_size_structural_uu_tmp = max((int)2, ghost_size_filter);
          structural_uu_tmp_tensor[0].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[1].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[2].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[3].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[4].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[5].echange_espace_virtuel(ghost_size_structural_uu_tmp);

          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[0]);
          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[1]);
          filtrer_champ_face(flag_add, structural_uu_tmp_tensor[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[2]);
          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[3], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[3]);
          filtrer_champ_face(flag_add, structural_uu_tmp_tensor[4], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[4]);
          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[5], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[5]);
        }
      else if ( structural_uu_model == Nom("su_laplacien_u") )
        {
          calculer_laplacien_u(velocity_filtre,
                               delta_z, facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                               velocity_filtre);

          velocity_filtre[0].echange_espace_virtuel(2);
          velocity_filtre[1].echange_espace_virtuel(2);
          velocity_filtre[2].echange_espace_virtuel(2);

          calculer_structural_uu_su_laplacien_u(structural_uu_model_constant,
                                                structural_uu_tensor_coefficients,
                                                velocity_filtre,
                                                velocity_filtre,
                                                structural_uu_filtre_tensor);
        }
      else if ( structural_uu_model == Nom("convection") )
        {
          calculer_structural_uu_convection(structural_uu_model_constant,
                                            structural_uu_tensor_coefficients,
                                            velocity_filtre,
                                            structural_uu_filtre_tensor);
        }
      else if ( structural_uu_model == Nom("convection_filtre") )
        {
          calculer_structural_uu_convection(structural_uu_model_constant,
                                            structural_uu_tensor_coefficients,
                                            velocity_filtre,
                                            structural_uu_tmp_tensor);

          ghost_size_filter = 1 + kernel->ghost_size();
          ghost_size_structural_uu_tmp = max((int)2, ghost_size_filter);
          structural_uu_tmp_tensor[0].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[1].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[2].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[3].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[4].echange_espace_virtuel(ghost_size_structural_uu_tmp);
          structural_uu_tmp_tensor[5].echange_espace_virtuel(ghost_size_structural_uu_tmp);

          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[0]);
          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[1]);
          filtrer_champ_face(flag_add, structural_uu_tmp_tensor[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[2]);
          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[3], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[3]);
          filtrer_champ_face(flag_add, structural_uu_tmp_tensor[4], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[4]);
          filtrer_champ_elem(flag_add, structural_uu_tmp_tensor[5], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uu_filtre_tensor[5]);
        }
      else if ( structural_uu_model == Nom("similarity") )
        {
          Cerr << "The use of a dynamic constant with the SIMILARITY structural model is not allowed." << finl;
          Process::exit();
        }
    }
}

void calculer_structural_uscalar_gradient(const double structural_uscalar_model_constant,
                                          const ArrOfDouble& structural_uscalar_vector_coefficients,
                                          const FixedVector<IJK_Field_double, 3>& velocity,
                                          const IJK_Field_double& champ_scalar,
                                          double scalar_kmin, double scalar_kmax,
                                          const ArrOfDouble_with_ghost& delta_z_maillage,
                                          const double facteur_x_pour_delta,
                                          const double facteur_y_pour_delta,
                                          const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                          FixedVector<IJK_Field_double, 3>& structural_uscalar_vector)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  IJK_Field_double& structural_uscalar_x = structural_uscalar_vector[0];
  IJK_Field_double& structural_uscalar_y = structural_uscalar_vector[1];
  IJK_Field_double& structural_uscalar_z = structural_uscalar_vector[2];

  const double& coefficient_x = structural_uscalar_vector_coefficients[0];
  const double& coefficient_y = structural_uscalar_vector_coefficients[1];
  const double& coefficient_z = structural_uscalar_vector_coefficients[2];

  const double deltaunsurdx = facteur_x_pour_delta;
  const double deltaunsurdy = facteur_y_pour_delta;

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              const int kg = k + offset;

              const double dz = delta_z_maillage[k];
              const double dz_m1 = kg==0 ? 0. : delta_z_maillage[k-1];
              const double delta_m = kg==0 ? 0.5*dz : 0.5*(dz + delta_z_maillage[k-1]);
              const double delta_p = kg==(nktot-1) ? 0.5*dz : 0.5*(dz + delta_z_maillage[k+1]);

              const double deltaunsurdz = delta_z_pour_delta[k] * 1./dz;
              const double deltaunsurdz_m1 = kg==0 ? 0. : delta_z_pour_delta[k-1] * 1./dz_m1;
              const double deltaunsurdelta_m = kg==0 ? 0. : 0.5*(delta_z_pour_delta[k] + delta_z_pour_delta[k-1]) * 1./delta_m;
              const double deltaunsurdelta_p = kg==(nktot-1) ? 0. : 0.5*(delta_z_pour_delta[k] + delta_z_pour_delta[k+1]) * 1./delta_p;

              const double uf_i = vitesse_i(i,j,k);
              const double uf_ip1 = vitesse_i(i+1,j,k);
              const double uf_im1 = vitesse_i(i-1,j,k);
              const double uf_i_jm1 = vitesse_i(i,j-1,k);
              const double uf_i_jp1 = vitesse_i(i,j+1,k);
              const double uf_i_km1 = kg==0 ? 0. : vitesse_i(i,j,k-1);
              const double uf_i_kp1 = kg==(nktot-1) ? 0. : vitesse_i(i,j,k+1);

              const double vf_j = vitesse_j(i,j,k);
              const double vf_j_im1 = vitesse_j(i-1,j,k);
              const double vf_j_ip1 = vitesse_j(i+1,j,k);
              const double vf_jp1 = vitesse_j(i,j+1,k);
              const double vf_jm1 = vitesse_j(i,j-1,k);
              const double vf_j_km1 = kg==0 ? 0. : vitesse_j(i,j,k-1);
              const double vf_j_kp1 = kg==(nktot-1) ? 0. : vitesse_j(i,j,k+1);

              const double wf_k = vitesse_k(i,j,k);
              const double wf_k_im1 = vitesse_k(i-1,j,k);
              const double wf_k_ip1 = vitesse_k(i+1,j,k);
              const double wf_k_jm1 = vitesse_k(i,j-1,k);
              const double wf_k_jp1 = vitesse_k(i,j+1,k);
              const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);
              const double wf_km1 = kg==0 ? 0. : vitesse_k(i,j,k-1);

              const double duidx = deltaunsurdx * (uf_ip1 - uf_i);
              const double dujdy = deltaunsurdy * (vf_jp1 - vf_j);
              const double dukdz = deltaunsurdz * (wf_kp1 - wf_k);

              const double duidx_im1 = deltaunsurdx * (uf_i - uf_im1);
              const double dujdy_jm1 = deltaunsurdy * (vf_j - vf_jm1);
              const double dukdz_km1 = deltaunsurdz_m1 * (wf_k - wf_km1);

              const double duidy_ij = deltaunsurdy * (uf_i - uf_i_jm1);
              const double duidy_ijp1 = deltaunsurdy * (uf_i_jp1 - uf_i);

              const double duidz_ik = deltaunsurdelta_m * (uf_i - uf_i_km1);
              const double duidz_ikp1 = deltaunsurdelta_p * (uf_i_kp1 - uf_i);

              const double dujdx_ij = deltaunsurdx * (vf_j - vf_j_im1);
              const double dujdx_ip1j = deltaunsurdx * (vf_j_ip1 - vf_j);

              const double dujdz_jk = deltaunsurdelta_m * (vf_j - vf_j_km1);
              const double dujdz_jkp1 = deltaunsurdelta_p * (vf_j_kp1 - vf_j);

              const double dukdx_ik = deltaunsurdx * (wf_k - wf_k_im1);
              const double dukdx_ip1k = deltaunsurdx * (wf_k_ip1 - wf_k);

              const double dukdy_jk = deltaunsurdy * (wf_k - wf_k_jm1);
              const double dukdy_jp1k = deltaunsurdy * (wf_k_jp1 - wf_k);

              const double g_fi_ii = 0.5 * (duidx_im1 + duidx);
              const double g_fi_ij = 0.5 * (duidy_ijp1 + duidy_ij);
              const double g_fi_ik = 0.5 * (duidz_ikp1 + duidz_ik);
              const double g_fj_ji = 0.5 * (dujdx_ip1j + dujdx_ij);
              const double g_fj_jj = 0.5 * (dujdy_jm1 + dujdy);
              const double g_fj_jk = 0.5 * (dujdz_jkp1 + dujdz_jk);
              const double g_fk_ki = 0.5 * (dukdx_ip1k + dukdx_ik);
              const double g_fk_kj = 0.5 * (dukdy_jp1k + dukdy_jk);
              const double g_fk_kk = 0.5 * (dukdz_km1 + dukdz);

              const double scalar = champ_scalar(i,j,k);

              const double scalar_im1 = champ_scalar(i-1,j,k);
              const double scalar_ip1 = champ_scalar(i+1,j,k);
              const double scalar_im1_jm1 = champ_scalar(i-1,j-1,k);
              const double scalar_ip1_jm1 = champ_scalar(i+1,j-1,k);
              const double scalar_im1_km1 = kg==0 ? scalar_kmin : champ_scalar(i-1,j,k-1);
              const double scalar_ip1_km1 = kg==0 ? scalar_kmin : champ_scalar(i+1,j,k-1);

              const double scalar_jm1 = champ_scalar(i,j-1,k);
              const double scalar_jp1 = champ_scalar(i,j+1,k);
              const double scalar_jm1_km1 = kg==0 ? scalar_kmin : champ_scalar(i,j-1,k-1);
              const double scalar_jp1_km1 = kg==0 ? scalar_kmin : champ_scalar(i,j+1,k-1);

              const double scalar_km1 = kg==0 ? scalar_kmin : champ_scalar(i,j,k-1);
              const double scalar_kp1 = kg==(nktot-1) ? scalar_kmax : champ_scalar(i,j,k+1);

              const double scalar_im1_jp1 = champ_scalar(i-1,j+1,k);
              const double scalar_im1_kp1 = kg==(nktot-1) ? scalar_kmax : champ_scalar(i-1,j,k+1);
              const double scalar_jm1_kp1 = kg==(nktot-1) ? scalar_kmax : champ_scalar(i,j-1,k+1);

              const double dscalardxf_i = deltaunsurdx * (scalar - scalar_im1);
              const double dscalardxf_ip1 = deltaunsurdx * (scalar_ip1 - scalar);
              const double dscalardxf_i_jm1 = deltaunsurdx * (scalar_jm1 - scalar_im1_jm1);
              const double dscalardxf_ip1_jm1 = deltaunsurdx * (scalar_ip1_jm1 - scalar_jm1);
              const double dscalardxf_i_km1 = deltaunsurdx * (scalar_km1 - scalar_im1_km1);
              const double dscalardxf_ip1_km1 = deltaunsurdx * (scalar_ip1_km1 - scalar_km1);

              const double dscalardyf_j = deltaunsurdy * (scalar - scalar_jm1);
              const double dscalardyf_jp1 = deltaunsurdy * (scalar_jp1 - scalar);
              const double dscalardyf_j_im1 = deltaunsurdy * (scalar_im1 - scalar_im1_jm1);
              const double dscalardyf_jp1_im1 = deltaunsurdy * (scalar_im1_jp1 - scalar_im1);
              const double dscalardyf_j_km1 = deltaunsurdy * (scalar_km1 - scalar_jm1_km1);
              const double dscalardyf_jp1_km1 = deltaunsurdy * (scalar_jp1_km1 - scalar_km1);

              const double dscalardzf_k = deltaunsurdelta_m * (scalar - scalar_km1);
              const double dscalardzf_kp1 = deltaunsurdelta_p * (scalar_kp1 - scalar);
              const double dscalardzf_k_im1 = deltaunsurdelta_m * (scalar_im1 - scalar_im1_km1);
              const double dscalardzf_kp1_im1 = deltaunsurdelta_p * (scalar_im1_kp1 - scalar_im1);
              const double dscalardzf_k_jm1 = deltaunsurdelta_m * (scalar_jm1 - scalar_jm1_km1);
              const double dscalardzf_kp1_jm1 = deltaunsurdelta_p * (scalar_jm1_kp1 - scalar_jm1);

              const double q_fi_i = dscalardxf_i;
              const double q_fi_j = 0.25 * (dscalardyf_jp1_im1 + dscalardyf_jp1 + dscalardyf_j_im1 + dscalardyf_j);
              const double q_fi_k = 0.25 * (dscalardzf_kp1_im1 + dscalardzf_kp1 + dscalardzf_k_im1 + dscalardzf_k);

              const double q_fj_i = 0.25 * (dscalardxf_ip1_jm1 + dscalardxf_ip1 + dscalardxf_i_jm1 + dscalardxf_i);
              const double q_fj_j = dscalardyf_j;
              const double q_fj_k = 0.25 * (dscalardzf_kp1_jm1 + dscalardzf_kp1 + dscalardzf_k_jm1 + dscalardzf_k);

              const double q_fk_i = 0.25 * (dscalardxf_ip1_km1 + dscalardxf_ip1 + dscalardxf_i_km1 + dscalardxf_i);
              const double q_fk_j = 0.25 * (dscalardyf_jp1_km1 + dscalardyf_jp1 + dscalardyf_j_km1 + dscalardyf_j);
              const double q_fk_k = dscalardzf_k;

              const double c_i = g_fi_ii*q_fi_i + g_fi_ij*q_fi_j + g_fi_ik*q_fi_k;
              const double c_j = g_fj_ji*q_fj_i + g_fj_jj*q_fj_j + g_fj_jk*q_fj_k;
              const double c_k = g_fk_ki*q_fk_i + g_fk_kj*q_fk_j + g_fk_kk*q_fk_k;

              structural_uscalar_x(i,j,k) = - coefficient_x * structural_uscalar_model_constant * c_i/12.;
              structural_uscalar_y(i,j,k) = - coefficient_y * structural_uscalar_model_constant * c_j/12.;
              structural_uscalar_z(i,j,k) = - coefficient_z * structural_uscalar_model_constant * c_k/12.;
            }
        }
    }
}


template<class T>
void calculer_structural_uscalar_similarity_comp(const double structural_uscalar_model_constant,
                                                 const ArrOfDouble& structural_uscalar_vector_coefficients,
                                                 const IJK_Field_double& rho,
                                                 const FixedVector<IJK_Field_double, 3>& velocity,
                                                 const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                                 const IJK_Field_double& champ_scalar,
                                                 const IJK_Field_double& champ_scalar_filtre,
                                                 double scalar_kmin, double scalar_kmax,
                                                 const ArrOfDouble_with_ghost& delta_z,
                                                 const double facteur_delta_x,
                                                 const double facteur_delta_y,
                                                 const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                                 T& kernel,
                                                 FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                                 FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                                 FixedVector<IJK_Field_double, 3>& structural_uscalar_vector)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];
  const IJK_Field_double& masse_vol_ijk = rho;
  IJK_Field_local_double& masse_vol_ijk_filtre = tmp_b[9];
  IJK_Field_local_double&   aa_rho_ijk_filtre = tmp_a[9];

  // const IJK_Field_double& vitesse_i_filtre = velocity_filtre[0];
  // const IJK_Field_double& vitesse_j_filtre = velocity_filtre[1];
  // const IJK_Field_double& vitesse_k_filtre = velocity_filtre[2];

  IJK_Field_double& structural_uscalar_x = structural_uscalar_vector[0];
  IJK_Field_double& structural_uscalar_y = structural_uscalar_vector[1];
  IJK_Field_double& structural_uscalar_z = structural_uscalar_vector[2];

  const double& coefficient_x = structural_uscalar_vector_coefficients[0];
  const double& coefficient_y = structural_uscalar_vector_coefficients[1];
  const double& coefficient_z = structural_uscalar_vector_coefficients[2];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& b_is = tmp_b[0];
  IJK_Field_local_double& b_js = tmp_b[1];
  IJK_Field_local_double& b_ks = tmp_b[2];
  IJK_Field_local_double& bb_is = tmp_b[3];
  IJK_Field_local_double& bb_js = tmp_b[4];
  IJK_Field_local_double& bb_ks = tmp_b[5];
  IJK_Field_local_double& bbT_i = tmp_b[6];
  IJK_Field_local_double& bbT_j = tmp_b[7];
  IJK_Field_local_double& bbT_k = tmp_b[8];

  IJK_Field_local_double& a_is = tmp_a[0];
  IJK_Field_local_double& a_js = tmp_a[1];
  IJK_Field_local_double& a_ks = tmp_a[2];
  IJK_Field_local_double& aa_is = tmp_a[3];
  IJK_Field_local_double& aa_js = tmp_a[4];
  IJK_Field_local_double& aa_ks = tmp_a[5];
  IJK_Field_local_double& aaT_i = tmp_a[6];
  IJK_Field_local_double& aaT_j = tmp_a[7];
  IJK_Field_local_double& aaT_k = tmp_a[8];

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_m1_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z[k-1];
      const double delta_m_glo = kg==0 ? 0.5*dz_glo : 0.5*(dz_glo + dz_m1_glo);

      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];
      const double dz_m1_pour_delta_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z_pour_delta[k-1];
      const double delta_m_pour_delta_glo = (kg-1<0 || kg>(nktot-1)) ? 0. : 0.5*(dz_pour_delta_glo + dz_m1_pour_delta_glo);

      const FixedVector<double, 21> filter_kernel_z = kernel->inhomogeneous(true, k, kg, nktot, dz_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_z_face = kernel->inhomogeneous(false, k, kg, nktot, delta_m_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_elem = kernel->size_k_elem(kg, nktot);
      const int size_k_face = kernel->size_k_face(kg, nktot);
      const int shift_k_elem = kernel->shift_k_elem(kg);
      const int shift_k_face = kernel->shift_k_face(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_elem = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_elem = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<-1 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double filter_coef_z = filter_kernel_z[kp+10];
                  longueur_elem += filter_coef_z * dz;
                }
              facteur_elem = 1./longueur_elem;
            }
          else
            {
              facteur_elem = dz_glo==0. ? 0. : 1./dz_glo;
            }
        }
      double facteur_face = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_face = 0.;
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<0 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                  const double dzf = 0.5*(dz + dz_m1);
                  const double filter_coef_z_face = filter_kernel_z_face[kp+10];
                  longueur_face += filter_coef_z_face * dzf;
                }
              facteur_face = 1./longueur_face;
            }
          else
            {
              facteur_face = delta_m_glo==0. ? 0. : 1./delta_m_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              b_is(i, j, 0) = 0.;
              b_js(i, j, 0) = 0.;
              b_ks(i, j, 0) = 0.;
              bb_is(i, j, 0) = 0.;
              bb_js(i, j, 0) = 0.;
              bb_ks(i, j, 0) = 0.;
              bbT_i(i, j, 0) = 0.;
              bbT_j(i, j, 0) = 0.;
              bbT_k(i, j, 0) = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_elem(kpg, nktot)))
                    {
                      const double filter_coef_z = filter_kernel_z[kp+10];

                      const double rho_ijk = masse_vol_ijk(i,j,k+kp);
                      const double scalar = champ_scalar(i,j,k+kp);
                      // const double scalar_im1 = champ_scalar(i-1,j,k+kp);
                      // const double scalar_jm1 = champ_scalar(i,j-1,k+kp);
                      const double scalarf_i = (scalar);
                      const double scalarf_j = (scalar);

                      const double uf_i = 0.5*(vitesse_i(i,j,k+kp)+vitesse_i(i+1,j,k+kp));
                      const double vf_j = 0.5*(vitesse_j(i,j,k+kp)+vitesse_j(i,j+1,k+kp));

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          b_is(i, j, 0) += rho_ijk*uf_i*scalarf_i * filter_coef_z      * dz  * facteur_elem;
                          b_js(i, j, 0) += rho_ijk*vf_j*scalarf_j * filter_coef_z      * dz  * facteur_elem;
                          bb_is(i, j, 0) += rho_ijk*uf_i * filter_coef_z      * dz  * facteur_elem;
                          bb_js(i, j, 0) += rho_ijk*vf_j * filter_coef_z      * dz  * facteur_elem;
                          bbT_i(i, j, 0) += rho_ijk*scalarf_i * filter_coef_z      * dz  * facteur_elem;
                          bbT_j(i, j, 0) += rho_ijk*scalarf_j * filter_coef_z      * dz  * facteur_elem;
                          masse_vol_ijk_filtre(i, j, 0) += rho_ijk     * filter_coef_z      * dz  * facteur_elem;
                        }
                      else
                        {
                          b_is(i, j, 0) += rho_ijk*uf_i*scalarf_i * filter_coef_z;
                          b_js(i, j, 0) += rho_ijk*vf_j*scalarf_j * filter_coef_z;
                          bb_is(i, j, 0) += rho_ijk*uf_i * filter_coef_z;
                          bb_js(i, j, 0) += rho_ijk*vf_j * filter_coef_z;
                          bbT_i(i, j, 0) += rho_ijk*scalarf_i * filter_coef_z;
                          bbT_j(i, j, 0) += rho_ijk*scalarf_j * filter_coef_z;
                          masse_vol_ijk_filtre(i, j, 0) += rho_ijk     * filter_coef_z;
                        }
                    }
                }
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_face(kpg, nktot)))
                    {
                      const double filter_coef_z_face = filter_kernel_z_face[kp+10];

                      const double rho_ijk = masse_vol_ijk(i,j,k+kp);
                      const double scalar = champ_scalar(i,j,k+kp);
                      // const double scalar_km1 = kpg==0 ? scalar_kmin : champ_scalar(i,j,k-1+kp);
                      const double scalarf_k = (scalar);

                      const double wf_k = kpg==(nktot-1) ? 0. : 0.5*(vitesse_k(i,j,k+kp)+vitesse_k(i,j,k+kp+1));

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                          const double dzf = 0.5*(dz + dz_m1);
                          b_ks(i, j, 0) += rho_ijk*wf_k*scalarf_k * filter_coef_z_face * dzf * facteur_face;
                          bb_ks(i, j, 0) += rho_ijk*wf_k*scalarf_k * filter_coef_z_face * dzf * facteur_face;
                          bbT_k(i, j, 0) += rho_ijk*scalarf_k * filter_coef_z_face * dzf * facteur_face;
                        }
                      else
                        {
                          b_ks(i, j, 0) += rho_ijk*wf_k*scalarf_k * filter_coef_z_face;
                          bb_ks(i, j, 0) += rho_ijk*wf_k*scalarf_k * filter_coef_z_face;
                          bbT_k(i, j, 0) += rho_ijk*scalarf_k * filter_coef_z_face;
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              a_is(i, 0, 0) = 0.;
              a_js(i, 0, 0) = 0.;
              a_ks(i, 0, 0) = 0.;
              aa_is(i, 0, 0) = 0.;
              aa_js(i, 0, 0) = 0.;
              aa_ks(i, 0, 0) = 0.;
              aaT_i(i, 0, 0) = 0.;
              aaT_j(i, 0, 0) = 0.;
              aaT_k(i, 0, 0) = 0.;
              aa_rho_ijk_filtre(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  a_is(i, 0, 0) += b_is(i, j+jp, 0) * filter_coef_y;
                  a_js(i, 0, 0) += b_js(i, j+jp, 0) * filter_coef_y;
                  a_ks(i, 0, 0) += b_ks(i, j+jp, 0) * filter_coef_y;
                  aa_is(i, 0, 0) += bb_is(i, j+jp, 0) * filter_coef_y;
                  aa_js(i, 0, 0) += bb_js(i, j+jp, 0) * filter_coef_y;
                  aa_ks(i, 0, 0) += bb_ks(i, j+jp, 0) * filter_coef_y;
                  aaT_i(i, 0, 0) += bbT_i(i, j+jp, 0) * filter_coef_y;
                  aaT_j(i, 0, 0) += bbT_j(i, j+jp, 0) * filter_coef_y;
                  aaT_k(i, 0, 0) += bbT_k(i, j+jp, 0) * filter_coef_y;
                  aa_rho_ijk_filtre(i, 0, 0) += masse_vol_ijk_filtre(i, j+jp, 0)* filter_coef_y;
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double r_is = 0.;
              double r_js = 0.;
              double r_ks = 0.;
              double rho_uf_i = 0.;
              double rho_vf_j = 0.;
              double rho_wf_k = 0.;
              double rho_scalarf_i = 0.;
              double rho_scalarf_j = 0.;
              double rho_scalarf_k = 0.;
              double rho_filtre = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  r_is += a_is(i+ip, 0, 0) * filter_coef_x;
                  r_js += a_js(i+ip, 0, 0) * filter_coef_x;
                  r_ks += a_ks(i+ip, 0, 0) * filter_coef_x;
                  rho_uf_i += aa_is(i+ip, 0, 0) * filter_coef_x;
                  rho_vf_j += aa_js(i+ip, 0, 0) * filter_coef_x;
                  rho_wf_k += aa_ks(i+ip, 0, 0) * filter_coef_x;
                  rho_scalarf_i += aaT_i(i+ip, 0, 0) * filter_coef_x;
                  rho_scalarf_j += aaT_j(i+ip, 0, 0) * filter_coef_x;
                  rho_scalarf_k += aaT_k(i+ip, 0, 0) * filter_coef_x;
                  rho_filtre += aa_rho_ijk_filtre(i+ip, 0, 0)* filter_coef_x;
                }

              Cerr << "i=" << i << finl;
              Cerr << "j=" << j << finl;
              Cerr << "k=" << k << finl;
              Cerr << "rho_filtre=" << rho_filtre << finl;
              Cerr << "r_is=" << r_is << finl;
              Cerr << "rho_uf_i=" << rho_uf_i << finl;
              Cerr << "rho_scalarf_i=" << rho_scalarf_i << finl;
              Cerr << "-----------------------" << finl;

              const double c_is = (r_is)/rho_filtre - (rho_uf_i*rho_scalarf_i)/(rho_filtre*rho_filtre);
              const double c_js = (r_js)/rho_filtre - (rho_vf_j*rho_scalarf_j)/(rho_filtre*rho_filtre);
              const double c_ks = (r_ks)/rho_filtre - (rho_wf_k*rho_scalarf_k)/(rho_filtre*rho_filtre);

              structural_uscalar_x(i,j,k) = - coefficient_x * structural_uscalar_model_constant * c_is;
              structural_uscalar_y(i,j,k) = - coefficient_y * structural_uscalar_model_constant * c_js;
              structural_uscalar_z(i,j,k) = - coefficient_z * structural_uscalar_model_constant * c_ks;
            }
        }
    }
}

template<class T>
void calculer_structural_uscalar_similarity(const double structural_uscalar_model_constant,
                                            const ArrOfDouble& structural_uscalar_vector_coefficients,
                                            const FixedVector<IJK_Field_double, 3>& velocity,
                                            const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                            const IJK_Field_double& champ_scalar,
                                            const IJK_Field_double& champ_scalar_filtre,
                                            double scalar_kmin, double scalar_kmax,
                                            const ArrOfDouble_with_ghost& delta_z,
                                            const double facteur_delta_x,
                                            const double facteur_delta_y,
                                            const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                            T& kernel,
                                            FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                            FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                            FixedVector<IJK_Field_double, 3>& structural_uscalar_vector)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const IJK_Field_double& vitesse_i_filtre = velocity_filtre[0];
  const IJK_Field_double& vitesse_j_filtre = velocity_filtre[1];
  const IJK_Field_double& vitesse_k_filtre = velocity_filtre[2];

  IJK_Field_double& structural_uscalar_x = structural_uscalar_vector[0];
  IJK_Field_double& structural_uscalar_y = structural_uscalar_vector[1];
  IJK_Field_double& structural_uscalar_z = structural_uscalar_vector[2];

  const double& coefficient_x = structural_uscalar_vector_coefficients[0];
  const double& coefficient_y = structural_uscalar_vector_coefficients[1];
  const double& coefficient_z = structural_uscalar_vector_coefficients[2];

  const double dx = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vitesse_k.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
  const double dx_pour_delta = facteur_delta_x*dx;
  const double dy_pour_delta = facteur_delta_y*dy;

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  IJK_Field_local_double& b_is = tmp_b[0];
  IJK_Field_local_double& b_js = tmp_b[1];
  IJK_Field_local_double& b_ks = tmp_b[2];

  IJK_Field_local_double& a_is = tmp_a[0];
  IJK_Field_local_double& a_js = tmp_a[1];
  IJK_Field_local_double& a_ks = tmp_a[2];

  const int ghost_size_filter = kernel->ghost_size();
  const int size_uniform = kernel->size_uniform();
  const int shift_uniform = kernel->shift_uniform();
  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;

      const double dz_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z[k];
      const double dz_m1_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z[k-1];
      const double delta_m_glo = kg==0 ? 0.5*dz_glo : 0.5*(dz_glo + dz_m1_glo);

      const double dz_pour_delta_glo = (kg<0 || kg>(nktot-1)) ? 0. : delta_z_pour_delta[k];
      const double dz_m1_pour_delta_glo = (kg-1<0 || kg-1>(nktot-1)) ? 0. : delta_z_pour_delta[k-1];
      const double delta_m_pour_delta_glo = (kg-1<0 || kg>(nktot-1)) ? 0. : 0.5*(dz_pour_delta_glo + dz_m1_pour_delta_glo);

      const FixedVector<double, 21> filter_kernel_z = kernel->inhomogeneous(true, k, kg, nktot, dz_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_z_face = kernel->inhomogeneous(false, k, kg, nktot, delta_m_pour_delta_glo, delta_z);
      const FixedVector<double, 21> filter_kernel_x = kernel->uniform(dx_pour_delta, dx);
      const FixedVector<double, 21> filter_kernel_y = kernel->uniform(dy_pour_delta, dy);
      const int size_k_elem = kernel->size_k_elem(kg, nktot);
      const int size_k_face = kernel->size_k_face(kg, nktot);
      const int shift_k_elem = kernel->shift_k_elem(kg);
      const int shift_k_face = kernel->shift_k_face(kg);
      const bool ponderation_filter_kernel = kernel->ponderation();
      const bool normalisation_filter_kernel = kernel->normalisation();

      double facteur_elem = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_elem = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<-1 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double filter_coef_z = filter_kernel_z[kp+10];
                  longueur_elem += filter_coef_z * dz;
                }
              facteur_elem = 1./longueur_elem;
            }
          else
            {
              facteur_elem = dz_glo==0. ? 0. : 1./dz_glo;
            }
        }
      double facteur_face = 0.;
      if (ponderation_filter_kernel)
        {
          if (normalisation_filter_kernel)
            {
              double longueur_face = 0.;
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (kpg<0 || kpg>nktot)
                    {
                      Cerr << "This should not happen." << finl;
                      Process::exit();
                    }
                  const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                  const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                  const double dzf = 0.5*(dz + dz_m1);
                  const double filter_coef_z_face = filter_kernel_z_face[kp+10];
                  longueur_face += filter_coef_z_face * dzf;
                }
              facteur_face = 1./longueur_face;
            }
          else
            {
              facteur_face = delta_m_glo==0. ? 0. : 1./delta_m_glo;
            }
        }

      for (int j = -ghost_size_filter; j < nj+ghost_size_filter; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              b_is(i, j, 0) = 0.;
              b_js(i, j, 0) = 0.;
              b_ks(i, j, 0) = 0.;
              for (int kp = -shift_k_elem; kp < size_k_elem-shift_k_elem; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_elem(kpg, nktot)))
                    {
                      const double filter_coef_z = filter_kernel_z[kp+10];

                      const double scalar = champ_scalar(i,j,k+kp);
                      const double scalar_im1 = champ_scalar(i-1,j,k+kp);
                      const double scalar_jm1 = champ_scalar(i,j-1,k+kp);
                      const double scalarf_i = 0.5*(scalar + scalar_im1);
                      const double scalarf_j = 0.5*(scalar + scalar_jm1);

                      const double uf_i = vitesse_i(i,j,k+kp);
                      const double vf_j = vitesse_j(i,j,k+kp);

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          b_is(i, j, 0) += uf_i*scalarf_i * filter_coef_z      * dz  * facteur_elem;
                          b_js(i, j, 0) += vf_j*scalarf_j * filter_coef_z      * dz  * facteur_elem;
                        }
                      else
                        {
                          b_is(i, j, 0) += uf_i*scalarf_i * filter_coef_z;
                          b_js(i, j, 0) += vf_j*scalarf_j * filter_coef_z;
                        }
                    }
                }
              for (int kp = -shift_k_face; kp < size_k_face-shift_k_face; kp++)
                {
                  const int kpg = kg + kp;
                  if (!(kernel->is_at_wall_face(kpg, nktot)))
                    {
                      const double filter_coef_z_face = filter_kernel_z_face[kp+10];

                      const double scalar = champ_scalar(i,j,k+kp);
                      const double scalar_km1 = kpg==0 ? scalar_kmin : champ_scalar(i,j,k-1+kp);
                      const double scalarf_k = 0.5*(scalar + scalar_km1);

                      const double wf_k = vitesse_k(i,j,k+kp);

                      if (ponderation_filter_kernel)
                        {
                          const double dz = (kpg<0 || kpg>(nktot-1)) ? 0. : delta_z[k+kp];
                          const double dz_m1 = (kpg-1<0 || kpg-1>(nktot-1)) ? 0. : delta_z[k-1+kp];
                          const double dzf = 0.5*(dz + dz_m1);
                          b_ks(i, j, 0) += wf_k*scalarf_k * filter_coef_z_face * dzf * facteur_face;
                        }
                      else
                        {
                          b_ks(i, j, 0) += wf_k*scalarf_k * filter_coef_z_face;
                        }
                    }
                }
            }
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = -ghost_size_filter; i < ni+ghost_size_filter; i++)
            {
              a_is(i, 0, 0) = 0.;
              a_js(i, 0, 0) = 0.;
              a_ks(i, 0, 0) = 0.;
              for (int jp = -shift_uniform; jp < size_uniform-shift_uniform; jp++)
                {
                  const double filter_coef_y = filter_kernel_y[jp+10];
                  a_is(i, 0, 0) += b_is(i, j+jp, 0) * filter_coef_y;
                  a_js(i, 0, 0) += b_js(i, j+jp, 0) * filter_coef_y;
                  a_ks(i, 0, 0) += b_ks(i, j+jp, 0) * filter_coef_y;
                }
            }

          for (int i = 0; i < ni; i++)
            {
              double r_is = 0.;
              double r_js = 0.;
              double r_ks = 0.;
              for (int ip = -shift_uniform; ip < size_uniform-shift_uniform; ip++)
                {
                  const double filter_coef_x = filter_kernel_x[ip+10];
                  r_is += a_is(i+ip, 0, 0) * filter_coef_x;
                  r_js += a_js(i+ip, 0, 0) * filter_coef_x;
                  r_ks += a_ks(i+ip, 0, 0) * filter_coef_x;
                }

              const double scalar = champ_scalar_filtre(i,j,k);
              const double scalar_im1 = champ_scalar_filtre(i-1,j,k);
              const double scalar_jm1 = champ_scalar_filtre(i,j-1,k);
              const double scalar_km1 = kg==0 ? scalar_kmin : champ_scalar_filtre(i,j,k-1);
              const double scalarf_i = 0.5*(scalar + scalar_im1);
              const double scalarf_j = 0.5*(scalar + scalar_jm1);
              const double scalarf_k = 0.5*(scalar + scalar_km1);

              const double uf_i = vitesse_i_filtre(i,j,k);
              const double vf_j = vitesse_j_filtre(i,j,k);
              const double wf_k = vitesse_k_filtre(i,j,k);

              const double c_is = r_is - uf_i*scalarf_i;
              const double c_js = r_js - vf_j*scalarf_j;
              const double c_ks = r_ks - wf_k*scalarf_k;

              structural_uscalar_x(i,j,k) = - coefficient_x * structural_uscalar_model_constant * c_is;
              structural_uscalar_y(i,j,k) = - coefficient_y * structural_uscalar_model_constant * c_js;
              structural_uscalar_z(i,j,k) = - coefficient_z * structural_uscalar_model_constant * c_ks;
            }
        }
    }
}

void calculer_structural_uscalar_similarity_Streher(const double structural_uscalar_model_constant,
                                                    const ArrOfDouble& structural_uscalar_vector_coefficients,
                                                    const FixedVector<IJK_Field_double, 3>& velocity,
                                                    const IJK_Field_double& champ_scalar,
                                                    double scalar_kmin, double scalar_kmax,
                                                    FixedVector<IJK_Field_double, 3>& structural_uscalar_vector)
{
  const IJK_Field_double& vitesse_i = velocity[0];
  const IJK_Field_double& vitesse_j = velocity[1];
  const IJK_Field_double& vitesse_k = velocity[2];

  const int ni = vitesse_k.ni();
  const int nj = vitesse_k.nj();
  const int nk = vitesse_k.nk();

  IJK_Field_double& structural_uscalar_x = structural_uscalar_vector[0];
  IJK_Field_double& structural_uscalar_y = structural_uscalar_vector[1];
  IJK_Field_double& structural_uscalar_z = structural_uscalar_vector[2];

  const int nktot = vitesse_k.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
  const int offset = vitesse_k.get_splitting().get_offset_local(DIRECTION_K);

  for (int k = 0; k < nk; k++)
    {
      const int kg = k + offset;
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {

              // const double uf_im1 = vitesse_i(i-1,j,k);
              const double uf_i = vitesse_i(i,j,k);
              // const double uf_ip1 = vitesse_i(i+1,j,k);
              // const double uf_ip2 = vitesse_i(i+2,j,k);

              // const double vf_jm1 = vitesse_j(i,j-1,k);
              const double vf_j = vitesse_j(i,j,k);
              // const double vf_jp1 = vitesse_j(i,j+1,k);
              // const double vf_jp2 = vitesse_j(i,j+2,k);

              // const double wf_km1 = kg<=1 ? 0. : vitesse_k(i,j,k-1);
              const double wf_k = vitesse_k(i,j,k);
              // const double wf_kp1 = kg==(nktot-1) ? 0. : vitesse_k(i,j,k+1);
              // const double wf_kp2 = kg>=(nktot-2) ? 0. : vitesse_k(i,j,k+2);

              const double sf_im2 = champ_scalar(i-2,j,k);
              const double sf_im1 = champ_scalar(i-1,j,k);
              const double sf_i = champ_scalar(i,j,k);
              const double sf_ip1 = champ_scalar(i+1,j,k);
              // const double sf_ip2 = champ_scalar(i+2,j,k);

              const double sf_jm2 = champ_scalar(i,j-2,k);
              const double sf_jm1 = champ_scalar(i,j-1,k);
              const double sf_j = champ_scalar(i,j,k);
              const double sf_jp1 = champ_scalar(i,j+1,k);
              // const double sf_jp2 = champ_scalar(i,j+2,k);

              const double sf_km2 = kg<=1 ? scalar_kmin : champ_scalar(i,j,k-2);
              const double sf_km1 = kg==0 ? scalar_kmin : champ_scalar(i,j,k-1);
              const double sf_k = kg >=(nktot-1) ? scalar_kmax : champ_scalar(i,j,k);
              const double sf_kp1 = kg >=(nktot-2) ? scalar_kmax : champ_scalar(i,j,k+1);
              // const double sf_kp2 = kg >=(nktot-3) ? scalar_kmax : champ_scalar(i,j,k+2);


              const double c_is = uf_i*(sf_i+sf_im1)/2. - (uf_i)*(sf_im2+sf_im1+sf_i+sf_ip1)/4.;
              const double c_js = vf_j*(sf_j+sf_jm1)/2. - (vf_j)*(sf_jm2+sf_jm1+sf_j+sf_jp1)/4.;
              const double c_ks = wf_k*(sf_k+sf_km1)/2. - (wf_k)*(sf_km2+sf_km1+sf_k+sf_kp1)/4.;

              structural_uscalar_x(i,j,k) = - structural_uscalar_model_constant * c_is;
              structural_uscalar_y(i,j,k) = - structural_uscalar_model_constant * c_js;
              structural_uscalar_z(i,j,k) = - structural_uscalar_model_constant * c_ks;

            }
        }
    }
}

template<class T>
void calculer_structural_uscalar(const Nom& structural_uscalar_model,
                                 const double structural_uscalar_model_constant,
                                 const ArrOfDouble& structural_uscalar_vector_coefficients,
                                 IJK_Field_double& rho,
                                 FixedVector<IJK_Field_double, 3>& velocity,
                                 FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                 IJK_Field_double& scalar,
                                 IJK_Field_double& scalar_filtre,
                                 double scalar_kmin, double scalar_kmax,
                                 const ArrOfDouble_with_ghost& delta_z,
                                 const double facteur_delta_x,
                                 const double facteur_delta_y,
                                 const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                 const double facteur_delta_filtre_x,
                                 const double facteur_delta_filtre_y,
                                 const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                 T& kernel,
                                 FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                 FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                 FixedVector<IJK_Field_double, 3>& structural_uscalar_tmp_vector,
                                 const bool flag_structural_uscalar_filtre,
                                 FixedVector<IJK_Field_double, 3>& structural_uscalar_filtre_vector,
                                 FixedVector<IJK_Field_double, 3>& structural_uscalar_vector)
{
  int ghost_size_filter;
  int ghost_size_velocity;
  int ghost_size_scalar;
  int ghost_size_structural_uscalar_tmp;


  if ( structural_uscalar_model == Nom("gradient") )
    {
      calculer_structural_uscalar_gradient(structural_uscalar_model_constant,
                                           structural_uscalar_vector_coefficients,
                                           velocity, scalar, scalar_kmin, scalar_kmax,
                                           delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                           structural_uscalar_vector);
    }
  else if ( structural_uscalar_model == Nom("gradient_filtre") )
    {
      calculer_structural_uscalar_gradient(structural_uscalar_model_constant,
                                           structural_uscalar_vector_coefficients,
                                           velocity, scalar, scalar_kmin, scalar_kmax,
                                           delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                           structural_uscalar_tmp_vector);

      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_structural_uscalar_tmp = max((int)2, ghost_size_filter);
      structural_uscalar_tmp_vector[0].echange_espace_virtuel(ghost_size_structural_uscalar_tmp);
      structural_uscalar_tmp_vector[1].echange_espace_virtuel(ghost_size_structural_uscalar_tmp);
      structural_uscalar_tmp_vector[2].echange_espace_virtuel(ghost_size_structural_uscalar_tmp);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, structural_uscalar_tmp_vector[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uscalar_vector[0]);
      filtrer_champ_elem(flag_add, structural_uscalar_tmp_vector[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uscalar_vector[1]);
      filtrer_champ_face(flag_add, structural_uscalar_tmp_vector[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uscalar_vector[2]);
    }
  else if ( structural_uscalar_model == Nom("similarity") )
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      ghost_size_scalar = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);
      scalar.echange_espace_virtuel(ghost_size_scalar);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);
      filtrer_champ_elem(flag_add, scalar, delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, scalar_filtre);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);
      scalar_filtre.echange_espace_virtuel(2);

      calculer_structural_uscalar_similarity(structural_uscalar_model_constant,
                                             structural_uscalar_vector_coefficients,
                                             velocity, velocity_filtre,
                                             scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                             delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                             kernel, tmp_b, tmp_a,
                                             structural_uscalar_vector);
    }
  else if ( structural_uscalar_model == Nom("similarity_comp") )
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      ghost_size_scalar = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);
      scalar.echange_espace_virtuel(ghost_size_scalar);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);
      filtrer_champ_elem(flag_add, scalar, delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, scalar_filtre);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);
      scalar_filtre.echange_espace_virtuel(2);

      calculer_structural_uscalar_similarity_comp(structural_uscalar_model_constant,
                                                  structural_uscalar_vector_coefficients,
                                                  rho,
                                                  velocity, velocity_filtre,
                                                  scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                                  delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                                  kernel, tmp_b, tmp_a,
                                                  structural_uscalar_vector);
    }
  else if ( structural_uscalar_model == Nom("similarity_streher") )
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      ghost_size_scalar = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);
      scalar.echange_espace_virtuel(ghost_size_scalar);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);
      filtrer_champ_elem(flag_add, scalar, delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, scalar_filtre);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);
      scalar_filtre.echange_espace_virtuel(2);

      calculer_structural_uscalar_similarity_Streher(structural_uscalar_model_constant,
                                                     structural_uscalar_vector_coefficients,
                                                     velocity,
                                                     scalar, scalar_kmin, scalar_kmax,
                                                     structural_uscalar_vector);
    }
  else
    {
      Cerr << "The name of the structural model for uscalar is unknown." << finl;
      Process::exit();
    }

  if (flag_structural_uscalar_filtre)
    {
      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_velocity = max((int)2, ghost_size_filter);
      ghost_size_scalar = max((int)2, ghost_size_filter);
      velocity[0].echange_espace_virtuel(ghost_size_velocity);
      velocity[1].echange_espace_virtuel(ghost_size_velocity);
      velocity[2].echange_espace_virtuel(ghost_size_velocity);
      scalar.echange_espace_virtuel(ghost_size_scalar);

      const int flag_add = 0;
      filtrer_champ_elem(flag_add, velocity[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[0]);
      filtrer_champ_elem(flag_add, velocity[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[1]);
      filtrer_champ_face(flag_add, velocity[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, velocity_filtre[2]);
      filtrer_champ_elem(flag_add, scalar, delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, scalar_filtre);

      velocity_filtre[0].echange_espace_virtuel(2);
      velocity_filtre[1].echange_espace_virtuel(2);
      velocity_filtre[2].echange_espace_virtuel(2);
      scalar_filtre.echange_espace_virtuel(2);

      if ( structural_uscalar_model == Nom("gradient") )
        {
          calculer_structural_uscalar_gradient(structural_uscalar_model_constant,
                                               structural_uscalar_vector_coefficients,
                                               velocity_filtre, scalar_filtre, scalar_kmin, scalar_kmax,
                                               delta_z, facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                               structural_uscalar_filtre_vector);
        }
      else if ( structural_uscalar_model == Nom("gradient_filtre") )
        {
          calculer_structural_uscalar_gradient(structural_uscalar_model_constant,
                                               structural_uscalar_vector_coefficients,
                                               velocity_filtre, scalar_filtre, scalar_kmin, scalar_kmax,
                                               delta_z, facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                               structural_uscalar_tmp_vector);

          ghost_size_filter = 1 + kernel->ghost_size();
          ghost_size_structural_uscalar_tmp = max((int)2, ghost_size_filter);
          structural_uscalar_tmp_vector[0].echange_espace_virtuel(ghost_size_structural_uscalar_tmp);
          structural_uscalar_tmp_vector[1].echange_espace_virtuel(ghost_size_structural_uscalar_tmp);
          structural_uscalar_tmp_vector[2].echange_espace_virtuel(ghost_size_structural_uscalar_tmp);

          filtrer_champ_elem(flag_add, structural_uscalar_tmp_vector[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uscalar_filtre_vector[0]);
          filtrer_champ_elem(flag_add, structural_uscalar_tmp_vector[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uscalar_filtre_vector[1]);
          filtrer_champ_face(flag_add, structural_uscalar_tmp_vector[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, structural_uscalar_filtre_vector[2]);
        }
      else if ( structural_uscalar_model == Nom("similarity") )
        {
          Cerr << "The use of a dynamic constant with the SIMILARITY structural model is not allowed." << finl;
          Process::exit();
        }
    }
}

template<class T>
void modification_modele_dynamic_uu_scalar(const bool anisotropic,
                                           const Nom& turbulent_viscosity_dynamic_type,
                                           const Nom& structural_uu_dynamic_type,
                                           const FixedVector<IJK_Field_double, 3>& velocity,
                                           const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                           const IJK_Field_double& scalar,
                                           const IJK_Field_double& scalar_filtre,
                                           double scalar_kmin, double scalar_kmax,
                                           const ArrOfDouble_with_ghost& delta_z,
                                           const double facteur_delta_x,
                                           const double facteur_delta_y,
                                           const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                           const double facteur_delta_filtre_x,
                                           const double facteur_delta_filtre_y,
                                           const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                           T& kernel,
                                           FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                           FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                           ArrOfDouble_with_ghost& constante_modele,
                                           FixedVector<FixedVector<ArrOfDouble, 7>, 8>& ml,
                                           const int turbulent_viscosity,
                                           IJK_Field_double& turbulent_mu,
                                           IJK_Field_double& turbulent_mu_filtre,
                                           const int structural_uu,
                                           FixedVector<IJK_Field_double, 6>& structural_uu_tensor,
                                           FixedVector<IJK_Field_double, 6>& structural_uu_filtre_tensor)
{
  int ghost_size_filter;
  int ghost_size_turbulent_mu;
  int ghost_size_structural_uu;
  if ( turbulent_viscosity_dynamic_type == Nom("not_dynamic")
       && structural_uu_dynamic_type == Nom("not_dynamic") )
    {
      // do nothing
    }
  else
    {
      if (turbulent_viscosity)
        {
          ghost_size_filter = 1 + kernel->ghost_size();
          ghost_size_turbulent_mu = max((int)2, ghost_size_filter);
          turbulent_mu.echange_espace_virtuel(ghost_size_turbulent_mu);
          turbulent_mu_filtre.echange_espace_virtuel(2);
        }

      if (structural_uu)
        {
          for (int j=0 ; j<6 ; j++)
            {
              ghost_size_filter = 1 + kernel->ghost_size();
              ghost_size_structural_uu = max((int)2, ghost_size_filter);
              structural_uu_tensor[j].echange_espace_virtuel(ghost_size_structural_uu);
              structural_uu_filtre_tensor[j].echange_espace_virtuel(2);
            }
        }

      const bool tensorial_struct = structural_uu_dynamic_type.finit_par("tensorial");

      calculer_ml_dynamic_uu_tensor(anisotropic, tensorial_struct,
                                    velocity, velocity_filtre,
                                    turbulent_viscosity,
                                    turbulent_mu,
                                    turbulent_mu,
                                    turbulent_mu,
                                    turbulent_mu,
                                    turbulent_mu,
                                    turbulent_mu,
                                    turbulent_mu_filtre,
                                    turbulent_mu_filtre,
                                    turbulent_mu_filtre,
                                    turbulent_mu_filtre,
                                    turbulent_mu_filtre,
                                    turbulent_mu_filtre,
                                    structural_uu,
                                    structural_uu_tensor,
                                    structural_uu_filtre_tensor,
                                    delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                    facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                    kernel, tmp_b, tmp_a,
                                    ml);

      FixedVector<ArrOfDouble, 7>& lij = ml[0];
      FixedVector<ArrOfDouble, 7>& mij = ml[1];
      FixedVector<ArrOfDouble, 7>& hij = ml[2];
      FixedVector<ArrOfDouble, 7>& mijmij = ml[3];
      FixedVector<ArrOfDouble, 7>& hijhij = ml[4];
      FixedVector<ArrOfDouble, 7>& mijlij = ml[5];
      FixedVector<ArrOfDouble, 7>& hijlij = ml[6];
      FixedVector<ArrOfDouble, 7>& mijhij = ml[7];

      ArrOfDouble& moy_lij = ml[0][6];
      ArrOfDouble& moy_mij = ml[1][6];
      ArrOfDouble& moy_hij = ml[2][6];
      ArrOfDouble& moy_mijmij = ml[3][6];
      ArrOfDouble& moy_hijhij = ml[4][6];
      ArrOfDouble& moy_mijlij = ml[5][6];
      ArrOfDouble& moy_hijlij = ml[6][6];
      ArrOfDouble& moy_mijhij = ml[7][6];

      const bool visc_reconnu = calculer_constante_modele(turbulent_viscosity_dynamic_type,
                                                          Nom("uu, viscosity"),
                                                          moy_lij, moy_mij, moy_hij,
                                                          moy_mijmij, moy_hijhij,
                                                          moy_mijlij, moy_hijlij,
                                                          moy_mijhij,
                                                          velocity, delta_z,
                                                          constante_modele);
      if (visc_reconnu)
        {
          const bool flag_face = false;
          multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu);
          multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_filtre);
        }

      const bool struct_reconnu = calculer_constante_modele(structural_uu_dynamic_type,
                                                            Nom("uu, structural"),
                                                            moy_lij, moy_hij, moy_mij,
                                                            moy_hijhij, moy_mijmij,
                                                            moy_hijlij, moy_mijlij,
                                                            moy_mijhij,
                                                            velocity, delta_z,
                                                            constante_modele);
      if (struct_reconnu)
        {
          for (int j=0 ; j<6 ; j++)
            {
              const bool flag_face = (j==2||j==4);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_tensor[j]);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_filtre_tensor[j]);
            }
        }

      if (tensorial_struct)
        {
          Nom struct_dynamic_type = structural_uu_dynamic_type.getPrefix("tensorial");
          for (int j=0 ; j<6 ; j++)
            {
              const bool reconnu = calculer_constante_modele(struct_dynamic_type,
                                                             Nom("uu, structural, tensorial, j= ") + Nom(j),
                                                             lij[j], hij[j], mij[j],
                                                             hijhij[j], mijmij[j],
                                                             hijlij[j], mijlij[j],
                                                             mijhij[j],
                                                             velocity, delta_z,
                                                             constante_modele);
              if (reconnu)
                {
                  const bool flag_face = (j==2||j==4);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_tensor[j]);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_filtre_tensor[j]);
                }
              else
                {
                  Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
                  Process::exit();
                }
            }
        }

      const bool twopass_visc = turbulent_viscosity_dynamic_type.finit_par("_twopass");
      const bool twopass_struct = structural_uu_dynamic_type.finit_par("_twopass");
      if (twopass_visc || twopass_struct)
        {
          if (!turbulent_viscosity)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'twopass' en ayant pas de viscosite turbulente" << finl;
              Process::exit();
            }
          if (!structural_uu)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'twopass' en ayant pas de modele structurel" << finl;
              Process::exit();
            }

          Nom visc_dynamic_type;
          Nom struct_dynamic_type;
          if (twopass_visc)
            {
              Cout << "Second pass, dynamic correction of constant: uu, viscosity" << finl;
              visc_dynamic_type = turbulent_viscosity_dynamic_type.getPrefix("_twopass");
              struct_dynamic_type = Nom("not_dynamic");
            }
          else
            {
              Cout << "Second pass, dynamic correction of constant: uu, structural" << finl;
              visc_dynamic_type = Nom("not_dynamic");
              struct_dynamic_type = structural_uu_dynamic_type.getPrefix("_twopass");
            }

          modification_modele_dynamic_uu_scalar(anisotropic,
                                                visc_dynamic_type,
                                                struct_dynamic_type,
                                                velocity, velocity_filtre,
                                                scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                                delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                                facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                                kernel, tmp_b, tmp_a,
                                                constante_modele, ml,
                                                turbulent_viscosity, turbulent_mu, turbulent_mu_filtre,
                                                structural_uu, structural_uu_tensor, structural_uu_filtre_tensor);
        }

      if ((!visc_reconnu) && (!struct_reconnu) && (!twopass_visc) && (!tensorial_struct) && (!twopass_struct))
        {
          Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
          Process::exit();
        }
    }
}

template<class T>
void modification_modele_dynamic_uu_tensor(const bool anisotropic,
                                           const Nom& turbulent_viscosity_dynamic_type,
                                           const Nom& structural_uu_dynamic_type,
                                           const FixedVector<IJK_Field_double, 3>& velocity,
                                           const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                           const IJK_Field_double& scalar,
                                           const IJK_Field_double& scalar_filtre,
                                           double scalar_kmin, double scalar_kmax,
                                           const ArrOfDouble_with_ghost& delta_z,
                                           const double facteur_delta_x,
                                           const double facteur_delta_y,
                                           const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                           const double facteur_delta_filtre_x,
                                           const double facteur_delta_filtre_y,
                                           const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                           T& kernel,
                                           FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                           FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                           ArrOfDouble_with_ghost& constante_modele,
                                           FixedVector<FixedVector<ArrOfDouble, 7>, 8>& ml,
                                           const int turbulent_viscosity,
                                           FixedVector<IJK_Field_double, 6>& turbulent_mu_tensor,
                                           FixedVector<IJK_Field_double, 6>& turbulent_mu_filtre_tensor,
                                           const int structural_uu,
                                           FixedVector<IJK_Field_double, 6>& structural_uu_tensor,
                                           FixedVector<IJK_Field_double, 6>& structural_uu_filtre_tensor)
{
  IJK_Field_double& turbulent_mu_xx = turbulent_mu_tensor[0];
  IJK_Field_double& turbulent_mu_xy = turbulent_mu_tensor[1];
  IJK_Field_double& turbulent_mu_xz = turbulent_mu_tensor[2];
  IJK_Field_double& turbulent_mu_yy = turbulent_mu_tensor[3];
  IJK_Field_double& turbulent_mu_yz = turbulent_mu_tensor[4];
  IJK_Field_double& turbulent_mu_zz = turbulent_mu_tensor[5];

  IJK_Field_double& turbulent_mu_filtre_xx = turbulent_mu_filtre_tensor[0];
  IJK_Field_double& turbulent_mu_filtre_xy = turbulent_mu_filtre_tensor[1];
  IJK_Field_double& turbulent_mu_filtre_xz = turbulent_mu_filtre_tensor[2];
  IJK_Field_double& turbulent_mu_filtre_yy = turbulent_mu_filtre_tensor[3];
  IJK_Field_double& turbulent_mu_filtre_yz = turbulent_mu_filtre_tensor[4];
  IJK_Field_double& turbulent_mu_filtre_zz = turbulent_mu_filtre_tensor[5];

  int ghost_size_filter, ghost_size_turbulent_mu, ghost_size_structural_uu;

  if ( turbulent_viscosity_dynamic_type == Nom("not_dynamic")
       && structural_uu_dynamic_type == Nom("not_dynamic") )
    {
      // do nothing
    }
  else
    {
      if (turbulent_viscosity)
        {
          for (int j=0 ; j<6 ; j++)
            {
              ghost_size_filter = 1 + kernel->ghost_size();
              ghost_size_turbulent_mu = max((int)2, ghost_size_filter);
              turbulent_mu_tensor[j].echange_espace_virtuel(ghost_size_turbulent_mu);
              turbulent_mu_filtre_tensor[j].echange_espace_virtuel(2);
            }
        }

      if (structural_uu)
        {
          for (int j=0 ; j<6 ; j++)
            {
              ghost_size_filter = 1 + kernel->ghost_size();
              ghost_size_structural_uu = max((int)2, ghost_size_filter);
              structural_uu_tensor[j].echange_espace_virtuel(ghost_size_structural_uu);
              structural_uu_filtre_tensor[j].echange_espace_virtuel(2);
            }
        }

      const bool tensorial_visc = turbulent_viscosity_dynamic_type.finit_par("tensorial");
      const bool tensorial_struct = structural_uu_dynamic_type.finit_par("tensorial");
      const bool tensorial = tensorial_visc || tensorial_struct;

      calculer_ml_dynamic_uu_tensor(anisotropic, tensorial,
                                    velocity, velocity_filtre,
                                    turbulent_viscosity,
                                    turbulent_mu_xx,
                                    turbulent_mu_xy,
                                    turbulent_mu_xz,
                                    turbulent_mu_yy,
                                    turbulent_mu_yz,
                                    turbulent_mu_zz,
                                    turbulent_mu_filtre_xx,
                                    turbulent_mu_filtre_xy,
                                    turbulent_mu_filtre_xz,
                                    turbulent_mu_filtre_yy,
                                    turbulent_mu_filtre_yz,
                                    turbulent_mu_filtre_zz,
                                    structural_uu,
                                    structural_uu_tensor,
                                    structural_uu_filtre_tensor,
                                    delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                    facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                    kernel, tmp_b, tmp_a,
                                    ml);

      FixedVector<ArrOfDouble, 7>& lij = ml[0];
      FixedVector<ArrOfDouble, 7>& mij = ml[1];
      FixedVector<ArrOfDouble, 7>& hij = ml[2];
      FixedVector<ArrOfDouble, 7>& mijmij = ml[3];
      FixedVector<ArrOfDouble, 7>& hijhij = ml[4];
      FixedVector<ArrOfDouble, 7>& mijlij = ml[5];
      FixedVector<ArrOfDouble, 7>& hijlij = ml[6];
      FixedVector<ArrOfDouble, 7>& mijhij = ml[7];

      ArrOfDouble& moy_lij = ml[0][6];
      ArrOfDouble& moy_mij = ml[1][6];
      ArrOfDouble& moy_hij = ml[2][6];
      ArrOfDouble& moy_mijmij = ml[3][6];
      ArrOfDouble& moy_hijhij = ml[4][6];
      ArrOfDouble& moy_mijlij = ml[5][6];
      ArrOfDouble& moy_hijlij = ml[6][6];
      ArrOfDouble& moy_mijhij = ml[7][6];

      const bool visc_reconnu = calculer_constante_modele(turbulent_viscosity_dynamic_type,
                                                          Nom("uu, viscosity"),
                                                          moy_lij, moy_mij, moy_hij,
                                                          moy_mijmij, moy_hijhij,
                                                          moy_mijlij, moy_hijlij,
                                                          moy_mijhij,
                                                          velocity, delta_z,
                                                          constante_modele);
      if (visc_reconnu)
        {
          for (int j=0 ; j<6 ; j++)
            {
              const bool flag_face = false;
              multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_tensor[j]);
              multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_filtre_tensor[j]);
            }
        }

      const bool struct_reconnu = calculer_constante_modele(structural_uu_dynamic_type,
                                                            Nom("uu, structural"),
                                                            moy_lij, moy_hij, moy_mij,
                                                            moy_hijhij, moy_mijmij,
                                                            moy_hijlij, moy_mijlij,
                                                            moy_mijhij,
                                                            velocity, delta_z,
                                                            constante_modele);
      if (struct_reconnu)
        {
          for (int j=0 ; j<6 ; j++)
            {
              const bool flag_face = (j==2||j==4);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_tensor[j]);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_filtre_tensor[j]);
            }
        }

      if (tensorial_visc)
        {
          Nom visc_dynamic_type = turbulent_viscosity_dynamic_type.getPrefix("tensorial");
          for (int j=0 ; j<6 ; j++)
            {
              const bool reconnu = calculer_constante_modele(visc_dynamic_type,
                                                             Nom("uu, viscosity, tensorial, j= ") + Nom(j),
                                                             lij[j], mij[j], hij[j],
                                                             mijmij[j], hijhij[j],
                                                             mijlij[j], hijlij[j],
                                                             mijhij[j],
                                                             velocity, delta_z,
                                                             constante_modele);
              if (reconnu)
                {
                  const bool flag_face = false;
                  multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_tensor[j]);
                  multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_filtre_tensor[j]);
                }
              else
                {
                  Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
                  Process::exit();
                }
            }
        }

      if (tensorial_struct)
        {
          Nom struct_dynamic_type = structural_uu_dynamic_type.getPrefix("tensorial");
          for (int j=0 ; j<6 ; j++)
            {
              const bool reconnu = calculer_constante_modele(struct_dynamic_type,
                                                             Nom("uu, structural, tensorial, j= ") + Nom(j),
                                                             lij[j], hij[j], mij[j],
                                                             hijhij[j], mijmij[j],
                                                             hijlij[j], mijlij[j],
                                                             mijhij[j],
                                                             velocity, delta_z,
                                                             constante_modele);
              if (reconnu)
                {
                  const bool flag_face = (j==2||j==4);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_tensor[j]);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uu_filtre_tensor[j]);
                }
              else
                {
                  Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
                  Process::exit();
                }
            }
        }

      const bool twopass_visc = turbulent_viscosity_dynamic_type.finit_par("_twopass");
      const bool twopass_struct = structural_uu_dynamic_type.finit_par("_twopass");
      if (twopass_visc || twopass_struct)
        {
          if (!turbulent_viscosity)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'twopass' en ayant pas de viscosite turbulente" << finl;
              Process::exit();
            }
          if (!structural_uu)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'twopass' en ayant pas de modele structurel" << finl;
              Process::exit();
            }

          Nom visc_dynamic_type;
          Nom struct_dynamic_type;
          if (twopass_visc)
            {
              Cout << "Second pass, dynamic correction of constant: uu, viscosity" << finl;
              visc_dynamic_type = turbulent_viscosity_dynamic_type.getPrefix("_twopass");
              struct_dynamic_type = Nom("not_dynamic");
            }
          else
            {
              Cout << "Second pass, dynamic correction of constant: uu, structural" << finl;
              visc_dynamic_type = Nom("not_dynamic");
              struct_dynamic_type = structural_uu_dynamic_type.getPrefix("_twopass");
            }

          modification_modele_dynamic_uu_tensor(anisotropic,
                                                visc_dynamic_type,
                                                struct_dynamic_type,
                                                velocity, velocity_filtre,
                                                scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                                delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                                facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                                kernel, tmp_b, tmp_a,
                                                constante_modele, ml,
                                                turbulent_viscosity, turbulent_mu_tensor, turbulent_mu_filtre_tensor,
                                                structural_uu, structural_uu_tensor, structural_uu_filtre_tensor);
        }

      if ((!visc_reconnu) && (!struct_reconnu) && (!tensorial_visc) && (!tensorial_struct) && (!twopass_visc) && (!twopass_struct))
        {
          Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
          Process::exit();
        }
    }
}

template<class T>
void modification_modele_dynamic_uscalar_scalar(const bool anisotropic,
                                                const Nom& turbulent_viscosity_dynamic_type,
                                                const Nom& structural_uscalar_dynamic_type,
                                                const FixedVector<IJK_Field_double, 3>& velocity,
                                                const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                                const IJK_Field_double& scalar,
                                                const IJK_Field_double& scalar_filtre,
                                                double scalar_kmin, double scalar_kmax,
                                                const ArrOfDouble_with_ghost& delta_z,
                                                const double facteur_delta_x,
                                                const double facteur_delta_y,
                                                const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                                const double facteur_delta_filtre_x,
                                                const double facteur_delta_filtre_y,
                                                const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                                T& kernel,
                                                FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                                FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                                ArrOfDouble_with_ghost& constante_modele,
                                                FixedVector<FixedVector<ArrOfDouble, 7>, 8>& ml,
                                                const int turbulent_viscosity,
                                                IJK_Field_double& turbulent_mu,
                                                IJK_Field_double& turbulent_mu_filtre,
                                                const int structural_uscalar,
                                                FixedVector<IJK_Field_double, 3>& structural_uscalar_vector,
                                                FixedVector<IJK_Field_double, 3>& structural_uscalar_filtre_vector)
{
  int ghost_size_filter, ghost_size_turbulent_mu, ghost_size_structural_uu;
  if ( turbulent_viscosity_dynamic_type == Nom("not_dynamic")
       && structural_uscalar_dynamic_type == Nom("not_dynamic") )
    {
      // do nothing
    }
  else
    {
      if (turbulent_viscosity)
        {
          ghost_size_filter = 1 + kernel->ghost_size();
          ghost_size_turbulent_mu = max((int)2, ghost_size_filter);
          turbulent_mu.echange_espace_virtuel(ghost_size_turbulent_mu);
          turbulent_mu_filtre.echange_espace_virtuel(2);
        }

      if (structural_uscalar)
        {
          for (int j=0 ; j<3 ; j++)
            {
              ghost_size_filter = 1 + kernel->ghost_size();
              ghost_size_structural_uu = max((int)2, ghost_size_filter);
              structural_uscalar_vector[j].echange_espace_virtuel(ghost_size_structural_uu);
              structural_uscalar_filtre_vector[j].echange_espace_virtuel(2);
            }
        }

      const bool vectorial_struct = structural_uscalar_dynamic_type.finit_par("vectorial");

      calculer_ml_dynamic_uscalar_vector(anisotropic, vectorial_struct,
                                         velocity, velocity_filtre,
                                         scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                         turbulent_viscosity,
                                         turbulent_mu,
                                         turbulent_mu,
                                         turbulent_mu,
                                         turbulent_mu_filtre,
                                         turbulent_mu_filtre,
                                         turbulent_mu_filtre,
                                         structural_uscalar,
                                         structural_uscalar_vector,
                                         structural_uscalar_filtre_vector,
                                         delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                         facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                         kernel, tmp_b, tmp_a,
                                         ml);

      FixedVector<ArrOfDouble, 7>& li = ml[0];
      FixedVector<ArrOfDouble, 7>& mi = ml[1];
      FixedVector<ArrOfDouble, 7>& hi = ml[2];
      FixedVector<ArrOfDouble, 7>& mimi = ml[3];
      FixedVector<ArrOfDouble, 7>& hihi = ml[4];
      FixedVector<ArrOfDouble, 7>& mili = ml[5];
      FixedVector<ArrOfDouble, 7>& hili = ml[6];
      FixedVector<ArrOfDouble, 7>& mihi = ml[7];

      ArrOfDouble& moy_li = ml[0][6];
      ArrOfDouble& moy_mi = ml[1][6];
      ArrOfDouble& moy_hi = ml[2][6];
      ArrOfDouble& moy_mimi = ml[3][6];
      ArrOfDouble& moy_hihi = ml[4][6];
      ArrOfDouble& moy_mili = ml[5][6];
      ArrOfDouble& moy_hili = ml[6][6];
      ArrOfDouble& moy_mihi = ml[7][6];

      const bool visc_reconnu = calculer_constante_modele(turbulent_viscosity_dynamic_type,
                                                          Nom("uscalar, viscosity"),
                                                          moy_li, moy_mi, moy_hi,
                                                          moy_mimi, moy_hihi,
                                                          moy_mili, moy_hili,
                                                          moy_mihi,
                                                          velocity, delta_z,
                                                          constante_modele);
      if (visc_reconnu)
        {
          const bool flag_face = false;
          multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu);
          multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_filtre);
        }

      const bool struct_reconnu = calculer_constante_modele(structural_uscalar_dynamic_type,
                                                            Nom("uscalar, structural"),
                                                            moy_li, moy_hi, moy_mi,
                                                            moy_hihi, moy_mimi,
                                                            moy_hili, moy_mili,
                                                            moy_mihi,
                                                            velocity, delta_z,
                                                            constante_modele);
      if (struct_reconnu)
        {
          for (int j=0 ; j<3 ; j++)
            {
              const bool flag_face = (j==2);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_vector[j]);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_filtre_vector[j]);
            }
        }

      if (vectorial_struct)
        {
          Nom struct_dynamic_type = structural_uscalar_dynamic_type.getPrefix("vectorial");
          for (int j=0 ; j<3 ; j++)
            {
              const bool reconnu = calculer_constante_modele(struct_dynamic_type,
                                                             Nom("uscalar, structural, vectorial, j= ") + Nom(j),
                                                             li[j], hi[j], mi[j],
                                                             hihi[j], mimi[j],
                                                             hili[j], mili[j],
                                                             mihi[j],
                                                             velocity, delta_z,
                                                             constante_modele);
              if (reconnu)
                {
                  const bool flag_face = (j==2);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_vector[j]);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_filtre_vector[j]);
                }
              else
                {
                  Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
                  Process::exit();
                }
            }
        }

      const bool twopass_visc = turbulent_viscosity_dynamic_type.finit_par("_twopass");
      const bool twopass_struct = structural_uscalar_dynamic_type.finit_par("_twopass");
      if (twopass_visc || twopass_struct)
        {
          if (!turbulent_viscosity)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'lillymixtedynamic' en ayant pas de viscosite turbulente" << finl;
              Process::exit();
            }
          if (!structural_uscalar)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'lillymixtedynamic' en ayant pas de modele structurel" << finl;
              Process::exit();
            }

          Nom visc_dynamic_type;
          Nom struct_dynamic_type;
          if (twopass_visc)
            {
              Cout << "Second pass, dynamic correction of constant: uscalar, viscosity" << finl;
              visc_dynamic_type = turbulent_viscosity_dynamic_type.getPrefix("_twopass");
              struct_dynamic_type = Nom("not_dynamic");
            }
          else
            {
              Cout << "Second pass, dynamic correction of constant: uscalar, structural" << finl;
              visc_dynamic_type = Nom("not_dynamic");
              struct_dynamic_type = structural_uscalar_dynamic_type.getPrefix("_twopass");
            }

          modification_modele_dynamic_uscalar_scalar(anisotropic,
                                                     visc_dynamic_type,
                                                     struct_dynamic_type,
                                                     velocity, velocity_filtre,
                                                     scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                                     delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                                     facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                                     kernel, tmp_b, tmp_a,
                                                     constante_modele, ml,
                                                     turbulent_viscosity, turbulent_mu, turbulent_mu_filtre,
                                                     structural_uscalar, structural_uscalar_vector, structural_uscalar_filtre_vector);
        }

      if ((!visc_reconnu) && (!struct_reconnu) && (!vectorial_struct) && (!twopass_visc) && (!twopass_struct))
        {
          Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
          Process::exit();
        }
    }
}

template<class T>
void modification_modele_dynamic_uscalar_vector(const bool anisotropic,
                                                const Nom& turbulent_viscosity_dynamic_type,
                                                const Nom& structural_uscalar_dynamic_type,
                                                const FixedVector<IJK_Field_double, 3>& velocity,
                                                const FixedVector<IJK_Field_double, 3>& velocity_filtre,
                                                const IJK_Field_double& scalar,
                                                const IJK_Field_double& scalar_filtre,
                                                double scalar_kmin, double scalar_kmax,
                                                const ArrOfDouble_with_ghost& delta_z,
                                                const double facteur_delta_x,
                                                const double facteur_delta_y,
                                                const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                                const double facteur_delta_filtre_x,
                                                const double facteur_delta_filtre_y,
                                                const ArrOfDouble_with_ghost& delta_z_pour_delta_filtre,
                                                T& kernel,
                                                FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                                FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                                ArrOfDouble_with_ghost& constante_modele,
                                                FixedVector<FixedVector<ArrOfDouble, 7>, 8>& ml,
                                                const int turbulent_viscosity,
                                                FixedVector<IJK_Field_double, 3>& turbulent_mu_vector,
                                                FixedVector<IJK_Field_double, 3>& turbulent_mu_filtre_vector,
                                                const int structural_uscalar,
                                                FixedVector<IJK_Field_double, 3>& structural_uscalar_vector,
                                                FixedVector<IJK_Field_double, 3>& structural_uscalar_filtre_vector)
{
  IJK_Field_double& turbulent_mu_x = turbulent_mu_vector[0];
  IJK_Field_double& turbulent_mu_y = turbulent_mu_vector[1];
  IJK_Field_double& turbulent_mu_z = turbulent_mu_vector[2];

  IJK_Field_double& turbulent_mu_filtre_x = turbulent_mu_filtre_vector[0];
  IJK_Field_double& turbulent_mu_filtre_y = turbulent_mu_filtre_vector[1];
  IJK_Field_double& turbulent_mu_filtre_z = turbulent_mu_filtre_vector[2];

  if ( turbulent_viscosity_dynamic_type == Nom("not_dynamic")
       && structural_uscalar_dynamic_type == Nom("not_dynamic") )
    {
      // do nothing
    }
  else
    {
      int ghost_size_filter, ghost_size_turbulent_mu, ghost_size_structural_uu;
      if (turbulent_viscosity)
        {
          for (int j=0 ; j<3 ; j++)
            {
              ghost_size_filter = 1 + kernel->ghost_size();
              ghost_size_turbulent_mu = max((int)2, ghost_size_filter);
              turbulent_mu_vector[j].echange_espace_virtuel(ghost_size_turbulent_mu);
              turbulent_mu_filtre_vector[j].echange_espace_virtuel(2);
            }
        }

      if (structural_uscalar)
        {
          for (int j=0 ; j<3 ; j++)
            {
              ghost_size_filter = 1 + kernel->ghost_size();
              ghost_size_structural_uu = max((int)2, ghost_size_filter);
              structural_uscalar_vector[j].echange_espace_virtuel(ghost_size_structural_uu);
              structural_uscalar_filtre_vector[j].echange_espace_virtuel(2);
            }
        }

      const bool vectorial_visc = turbulent_viscosity_dynamic_type.finit_par("vectorial");
      const bool vectorial_struct = structural_uscalar_dynamic_type.finit_par("vectorial");
      const bool vectorial = vectorial_visc || vectorial_struct;

      calculer_ml_dynamic_uscalar_vector(anisotropic, vectorial,
                                         velocity, velocity_filtre,
                                         scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                         turbulent_viscosity,
                                         turbulent_mu_x,
                                         turbulent_mu_y,
                                         turbulent_mu_z,
                                         turbulent_mu_filtre_x,
                                         turbulent_mu_filtre_y,
                                         turbulent_mu_filtre_z,
                                         structural_uscalar,
                                         structural_uscalar_vector,
                                         structural_uscalar_filtre_vector,
                                         delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                         facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                         kernel, tmp_b, tmp_a,
                                         ml);

      FixedVector<ArrOfDouble, 7>& li = ml[0];
      FixedVector<ArrOfDouble, 7>& mi = ml[1];
      FixedVector<ArrOfDouble, 7>& hi = ml[2];
      FixedVector<ArrOfDouble, 7>& mimi = ml[3];
      FixedVector<ArrOfDouble, 7>& hihi = ml[4];
      FixedVector<ArrOfDouble, 7>& mili = ml[5];
      FixedVector<ArrOfDouble, 7>& hili = ml[6];
      FixedVector<ArrOfDouble, 7>& mihi = ml[7];

      ArrOfDouble& moy_li = ml[0][6];
      ArrOfDouble& moy_mi = ml[1][6];
      ArrOfDouble& moy_hi = ml[2][6];
      ArrOfDouble& moy_mimi = ml[3][6];
      ArrOfDouble& moy_hihi = ml[4][6];
      ArrOfDouble& moy_mili = ml[5][6];
      ArrOfDouble& moy_hili = ml[6][6];
      ArrOfDouble& moy_mihi = ml[7][6];

      const bool visc_reconnu = calculer_constante_modele(turbulent_viscosity_dynamic_type,
                                                          Nom("uscalar, viscosity"),
                                                          moy_li, moy_mi, moy_hi,
                                                          moy_mimi, moy_hihi,
                                                          moy_mili, moy_hili,
                                                          moy_mihi,
                                                          velocity, delta_z,
                                                          constante_modele);
      if (visc_reconnu)
        {
          for (int j=0 ; j<3 ; j++)
            {
              const bool flag_face = false;
              multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_vector[j]);
              multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_filtre_vector[j]);
            }
        }

      const bool struct_reconnu = calculer_constante_modele(structural_uscalar_dynamic_type,
                                                            Nom("uscalar, structural"),
                                                            moy_li, moy_hi, moy_mi,
                                                            moy_hihi, moy_mimi,
                                                            moy_hili, moy_mili,
                                                            moy_mihi,
                                                            velocity, delta_z,
                                                            constante_modele);
      if (struct_reconnu)
        {
          for (int j=0 ; j<3 ; j++)
            {
              const bool flag_face = (j==2);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_vector[j]);
              multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_filtre_vector[j]);
            }
        }

      if (vectorial_visc)
        {
          Nom visc_dynamic_type = turbulent_viscosity_dynamic_type.getPrefix("vectorial");
          for (int j=0 ; j<3 ; j++)
            {
              const bool reconnu = calculer_constante_modele(visc_dynamic_type,
                                                             Nom("uscalar, viscosity, vectorial, j= ") + Nom(j),
                                                             li[j], mi[j], hi[j],
                                                             mimi[j], hihi[j],
                                                             mili[j], hili[j],
                                                             mihi[j],
                                                             velocity, delta_z,
                                                             constante_modele);
              if (reconnu)
                {
                  const bool flag_face = false;
                  multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_vector[j]);
                  multiplier_par_constante(flag_face, constante_modele, velocity, turbulent_mu_filtre_vector[j]);
                }
              else
                {
                  Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
                  Process::exit();
                }
            }
        }

      if (vectorial_struct)
        {
          Nom struct_dynamic_type = structural_uscalar_dynamic_type.getPrefix("vectorial");
          for (int j=0 ; j<3 ; j++)
            {
              const bool reconnu = calculer_constante_modele(struct_dynamic_type,
                                                             Nom("uscalar, structural, vectorial, j= ") + Nom(j),
                                                             li[j], hi[j], mi[j],
                                                             hihi[j], mimi[j],
                                                             hili[j], mili[j],
                                                             mihi[j],
                                                             velocity, delta_z,
                                                             constante_modele);
              if (reconnu)
                {
                  const bool flag_face = (j==2);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_vector[j]);
                  multiplier_par_constante(flag_face, constante_modele, velocity, structural_uscalar_filtre_vector[j]);
                }
              else
                {
                  Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
                  Process::exit();
                }
            }
        }

      const bool twopass_visc = turbulent_viscosity_dynamic_type.finit_par("_twopass");
      const bool twopass_struct = structural_uscalar_dynamic_type.finit_par("_twopass");
      if (twopass_visc || twopass_struct)
        {
          if (!turbulent_viscosity)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'lillymixtedynamic' en ayant pas de viscosite turbulente" << finl;
              Process::exit();
            }
          if (!structural_uscalar)
            {
              Cerr << "Erreur : On ne devrait pas pouvoir rentrer dans une evaluation dynamique de la constante de type 'lillymixtedynamic' en ayant pas de modele structurel" << finl;
              Process::exit();
            }

          Nom visc_dynamic_type;
          Nom struct_dynamic_type;
          if (twopass_visc)
            {
              Cout << "Second pass, dynamic correction of constant: uscalar, viscosity" << finl;
              visc_dynamic_type = turbulent_viscosity_dynamic_type.getPrefix("_twopass");
              struct_dynamic_type = Nom("not_dynamic");
            }
          else
            {
              Cout << "Second pass, dynamic correction of constant: uscalar, structural" << finl;
              visc_dynamic_type = Nom("not_dynamic");
              struct_dynamic_type = structural_uscalar_dynamic_type.getPrefix("_twopass");
            }

          modification_modele_dynamic_uscalar_vector(anisotropic,
                                                     visc_dynamic_type,
                                                     struct_dynamic_type,
                                                     velocity, velocity_filtre,
                                                     scalar, scalar_filtre, scalar_kmin, scalar_kmax,
                                                     delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta,
                                                     facteur_delta_filtre_x, facteur_delta_filtre_y, delta_z_pour_delta_filtre,
                                                     kernel, tmp_b, tmp_a,
                                                     constante_modele, ml,
                                                     turbulent_viscosity, turbulent_mu_vector, turbulent_mu_filtre_vector,
                                                     structural_uscalar, structural_uscalar_vector, structural_uscalar_filtre_vector);
        }

      if ((!visc_reconnu) && (!struct_reconnu) && (!vectorial_visc) && (!vectorial_struct) && (!twopass_visc) && (!twopass_struct))
        {
          Cerr << "Error: The type of the dynamic correction of the constant is invalid." << finl;
          Process::exit();
        }
    }
}


Entree& DNS_QC_double::interpreter(Entree& is)
{
  //M.D. 12-02-2019 ajout possibilite puit
  puit_ = 0.;
  fichier_reprise_rho_ = "??";
  fichier_reprise_vitesse_ = "??";
  timestep_reprise_rho_ = -2;
  timestep_reprise_vitesse_ = -2;
  expression_vitesse_initiale_.dimensionner(3);
  t_debut_statistiques_ = 0.;
  check_divergence_ = 0;
  Param param(que_suis_je());
  Nom ijk_splitting_name;
  reprise_fourier_ = 0;
  dt_post_ = 100;
  timestep_facsec_ = 1.;
  timestep_max_ = 1.;
  timestep_ = 1.e28;
  old_timestep_= 1.e28; // Evite la division par zero
  current_time_ = 0.;
  dt_sauvegarde_ = 2000000000; // jamais
  nom_sauvegarde_ = nom_du_cas() + ".sauv";
  nom_reprise_ = "??"; // pas de reprise
  compteur_post_instantanes_ = 0;
  postraiter_sous_pas_de_temps_ = 0;
  calcul_2d_ = 0;
  projection_initiale_demandee_ = 0;
  // F.A 3/03/14 modification source
  dump_factor_ = 10./3.;
  // F.A 17/03/14 ajout diffs et convections negligeables
  diff_qdm_negligeable_ = 0;
  diff_temp_negligeable_= 0;
  conv_qdm_negligeable_ = 0;
  conv_rho_negligeable_ = 0;
  // DD,2016-10-14: ajout disable_solveur_poisson_
  disable_solveur_poisson_ = 0;
  // F.A dt_start ajout possibilite de dt_start
  dt_start_ = 1.e28; // tres grand

  //oscillating boundary condition
  amplitude_oscillating_boundary_=0;
  frequency_oscillating_boundary_=0;

  /* amplitude_oscillating_boundary_2_=0; */
  /* frequency_oscillating_boundary_2_=0; */
  /* flag_oscillating_boundary_2_=false; */
  /* statistiques spectrales */
  dt_post_spectral_ = -1;
  Nom post_splitting_name = "??";

  /* sauvegarde des lata par plan */
  sauvegarde_splitting_name_ = "??";

  // DD,2016-28-01: statistiques a partir de fichiers lata uniquement
  sauvegarde_post_instantanes_ = 0;
  lecture_post_instantanes_ = 0;
  lecture_post_instantanes_filtrer_u_ = 0;
  lecture_post_instantanes_filtrer_rho_ = 0;
  lecture_post_instantanes_filtrer_p_ = 0;
  lecture_post_instantanes_filtrer_tous_ = 0;

  // DD,2016-06-15: changement des pas de temps de stabilite,
  old_dtstab_ = 0;

  debit_cible_=0;
  aim_t_bulk_=-1;
  dump_factor_2_=3./100.;

  convection_rho_amont_ = 0;
  convection_rho_centre2_ = 0;
  convection_rho_centre4_ = 0;

  convection_velocity_amont_ = 0;
  convection_velocity_centre2_ = 0;
  convection_velocity_quicksharp_ = 0;

  terme_source_acceleration_=0;
  terme_source_acceleration_constant_=0;
  mode_terme_source_impose_=0;

  variation_cste_modele_fonctionnel_ = 0;
  smoothing_center_fr_=0.;
  smoothing_factor_fr_=0.;
  Re_tau_fr_=1.;
  Re_tau_ch_=1.;
  pond_fr_=0.;
  pond_ch_=0.;
  center_constant_=0.;

  // DD,2017-04-27: diffusion modifie en vue de l'ajout de modeles
  type_velocity_diffusion_ = Nom("full");

  type_velocity_turbulent_diffusion_ = Nom("none");
  turbulent_viscosity_model_ = Nom("none");
  turbulent_viscosity_dynamic_type_ = Nom("not_dynamic");
  turbulent_viscosity_model_constant_ = -1;
  turbulent_viscosity_tensor_coefficients_.resize_array(6);
  turbulent_viscosity_tensor_coefficients_ = 1.;
  turbulent_viscosity_ = 0;

  type_scalar_turbulent_diffusion_ = Nom("none");
  turbulent_diffusivity_model_ = Nom("none");
  turbulent_diffusivity_dynamic_type_ = Nom("not_dynamic");
  turbulent_diffusivity_model_constant_ = -1;
  turbulent_diffusivity_vector_coefficients_.resize_array(3);
  turbulent_diffusivity_vector_coefficients_ = 1.;
  turbulent_diffusivity_ = 0;

  structural_uu_model_ = Nom("none");
  structural_uu_dynamic_type_ = Nom("not_dynamic");
  structural_uu_model_constant_ = -1;
  structural_uu_tensor_coefficients_.resize_array(6);
  structural_uu_tensor_coefficients_ = 1.;
  structural_uu_ = 0;

  structural_uscalar_model_ = Nom("none");
  structural_uscalar_dynamic_type_ = Nom("not_dynamic");
  structural_uscalar_model_constant_ = -1;
  structural_uscalar_vector_coefficients_.resize_array(3);
  structural_uscalar_vector_coefficients_ = 1.;
  structural_uscalar_ = 0;

  filter_kernel_name_ = Nom("none");
  flag_filtrage_convection_qdm_ = 0;
  flag_filtrage_turbulent_diffusion_qdm_ = 0;
  flag_filtrage_structural_diffusion_qdm_ = 0;
  flag_convection_qdm_sans_rho_ = 0;
  flag_convection_qdm_sans_divergence_ = 0;

  large_eddy_simulation_formulation_ = Nom("none");
  Nom geom_name_pour_delta = Nom("__identique_au_maillage__");
  formulation_velocity_ = 0;
  formulation_favre_ = 0;

  param.ajouter("ijk_splitting", &ijk_splitting_name, Param::REQUIRED);
  param.ajouter("tinit", &current_time_);
  param.ajouter("timestep", &timestep_max_, Param::REQUIRED);
  param.ajouter("timestep_facsec", &timestep_facsec_);
  param.ajouter("dt_start", &dt_start_);
  param.ajouter("nb_pas_dt_max", &nb_timesteps_, Param::REQUIRED);
  param.ajouter("multigrid_solver", &poisson_solver_, Param::REQUIRED);
  param.ajouter_flag("check_divergence", &check_divergence_);

  param.ajouter("t_paroi_impose_kmin", &T_paroi_impose_kmin_, Param::REQUIRED);
  param.ajouter("t_paroi_impose_kmax", &T_paroi_impose_kmax_, Param::REQUIRED);

  param.ajouter("p_thermo_init", &P_thermodynamique_, Param::REQUIRED);
  param.ajouter("puit", &puit_);
  param.ajouter("smoothing_center_fr", &smoothing_center_fr_);
  param.ajouter("smoothing_factor_fr", &smoothing_factor_fr_);
  param.ajouter("Re_tau_fr", &Re_tau_fr_);
  param.ajouter("Re_tau_ch", &Re_tau_ch_);
  param.ajouter("ponderation_fr", &pond_fr_);
  param.ajouter("ponderation_ch", &pond_ch_);
  param.ajouter("center_constant", &center_constant_);
  param.ajouter("cp", &Cp_gaz_, Param::REQUIRED);
  param.ajouter("gamma", &gamma_, Param::REQUIRED);


  param.ajouter("expression_vx_init", &expression_vitesse_initiale_[0]);
  param.ajouter("expression_vy_init", &expression_vitesse_initiale_[1]);
  param.ajouter("expression_vz_init", &expression_vitesse_initiale_[2]);
  param.ajouter("expression_t_init", &expression_temperature_initiale_);
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
  param.ajouter("fichier_reprise_rho", &fichier_reprise_rho_);
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
  param.ajouter("timestep_reprise_rho", &timestep_reprise_rho_);

  // F.A on change la source utilisee
  param.ajouter("debit_massique", &debit_cible_);
  // Y.Z. forcing T^b by changing H_s (heat sink)
  param.ajouter("t_bulk", &aim_t_bulk_);
  param.ajouter("dump_factor_2", &dump_factor_2_);

  param.ajouter_flag("convection_rho_amont", &convection_rho_amont_);
  param.ajouter_flag("convection_rho_centre2", &convection_rho_centre2_);
  param.ajouter_flag("convection_rho_centre4", &convection_rho_centre4_);

  param.ajouter_flag("convection_velocity_amont", &convection_velocity_amont_);
  param.ajouter_flag("convection_velocity_quicksharp", &convection_velocity_quicksharp_);
  param.ajouter_flag("convection_velocity_centre2", &convection_velocity_centre2_);

  param.ajouter_flag("variation_cste_modele_fonctionnel", &variation_cste_modele_fonctionnel_);

  param.ajouter("dumping_factor", &dump_factor_);

  // param.ajouter("terme_force_init", &terme_source_acceleration_, Param::REQUIRED);
  // param.ajouter("expression_derivee_force", &expression_derivee_acceleration_, Param::REQUIRED);
  param.ajouter("terme_source_acceleration_constant", &terme_source_acceleration_constant_);

  param.ajouter("dt_post", &dt_post_);
  param.ajouter("champs_a_postraiter", &liste_post_instantanes_);

  param.ajouter("dt_sauvegarde", &dt_sauvegarde_);
  param.ajouter_flag("postraiter_sous_pas_de_temps", &postraiter_sous_pas_de_temps_);
  param.ajouter("nom_sauvegarde", &nom_sauvegarde_);
  param.ajouter("nom_reprise", &nom_reprise_);

  param.ajouter("t_debut_statistiques", &t_debut_statistiques_);

  /* stats spectrales */
  param.ajouter("dt_post_spectral", &dt_post_spectral_);
  param.ajouter("spectral_splitting", &post_splitting_name);
  param.ajouter_flag("reprise_spectrale", &reprise_fourier_);
  /* -------------------- */

  /* sauvegarde des lata par plan */
  param.ajouter("sauvegarde_splitting", &sauvegarde_splitting_name_);

  param.ajouter_flag("calcul_2d", &calcul_2d_);
  param.ajouter("check_stop_file", &check_stop_file_);
  param.ajouter_flag("projection_initiale", &projection_initiale_demandee_);

  // Ajout possibilite diffusion ou convection negligeable pour chaque equa  (sauf celle de la chaleur !! ) !!
  param.ajouter_flag("convection_rho_negligeable", &conv_rho_negligeable_);
  param.ajouter_flag("convection_qdm_negligeable", &conv_qdm_negligeable_);
  param.ajouter_flag("diffusion_qdm_negligeable", &diff_qdm_negligeable_);
  param.ajouter_flag("diffusion_temp_negligeable", &diff_temp_negligeable_);
  // DD,2016-10-14: ajout disable_solveur_poisson_
  param.ajouter_flag("disable_solveur_poisson", &disable_solveur_poisson_);

  // DD,2017-04-27: diffusion modifie en vue de l'ajout de modeles
  param.ajouter("type_velocity_diffusion", &type_velocity_diffusion_);

  param.ajouter("large_eddy_simulation_formulation", &large_eddy_simulation_formulation_);
  param.ajouter("ijk_grid_geometry_pour_delta", &geom_name_pour_delta);

  param.ajouter("type_velocity_turbulent_diffusion", &type_velocity_turbulent_diffusion_);
  param.ajouter("turbulent_viscosity_model", &turbulent_viscosity_model_);
  param.ajouter("turbulent_viscosity_dynamic_type", &turbulent_viscosity_dynamic_type_);
  param.ajouter("turbulent_viscosity_model_constant", &turbulent_viscosity_model_constant_);
  param.ajouter("turbulent_viscosity_tensor_coefficients", &turbulent_viscosity_tensor_coefficients_);
  param.ajouter_flag("turbulent_viscosity", &turbulent_viscosity_);

  param.ajouter("type_scalar_turbulent_diffusion", &type_scalar_turbulent_diffusion_);
  param.ajouter("turbulent_diffusivity_model", &turbulent_diffusivity_model_);
  param.ajouter("turbulent_diffusivity_dynamic_type", &turbulent_diffusivity_dynamic_type_);
  param.ajouter("turbulent_diffusivity_model_constant", &turbulent_diffusivity_model_constant_);
  param.ajouter("turbulent_diffusivity_vector_coefficients", &turbulent_diffusivity_vector_coefficients_);
  param.ajouter_flag("turbulent_diffusivity", &turbulent_diffusivity_);

  param.ajouter("structural_uu_model", &structural_uu_model_);
  param.ajouter("structural_uu_dynamic_type", &structural_uu_dynamic_type_);
  param.ajouter("structural_uu_model_constant", &structural_uu_model_constant_);
  param.ajouter("structural_uu_tensor_coefficients", &structural_uu_tensor_coefficients_);
  param.ajouter_flag("structural_uu", &structural_uu_);

  param.ajouter("structural_uscalar_model", &structural_uscalar_model_);
  param.ajouter("structural_uscalar_dynamic_type", &structural_uscalar_dynamic_type_);
  param.ajouter("structural_uscalar_model_constant", &structural_uscalar_model_constant_);
  param.ajouter("structural_uscalar_vector_coefficients", &structural_uscalar_vector_coefficients_);
  param.ajouter_flag("structural_uscalar", &structural_uscalar_);

  param.ajouter("filter_kernel_name", &filter_kernel_name_);
  param.ajouter_flag("filtrage_convection_qdm", &flag_filtrage_convection_qdm_);
  param.ajouter_flag("filtrage_turbulent_diffusion_qdm", &flag_filtrage_turbulent_diffusion_qdm_);
  param.ajouter_flag("filtrage_structural_diffusion_qdm", &flag_filtrage_structural_diffusion_qdm_);
  param.ajouter_flag("convection_qdm_sans_rho", &flag_convection_qdm_sans_rho_);
  param.ajouter_flag("convection_qdm_sans_divergence", &flag_convection_qdm_sans_divergence_);

  // DD,2016-28-01: statistiques a partir de fichiers lata uniquement
  param.ajouter("statlata_namelist", &statlata_namelist_);
  param.ajouter_flag("sauvegarde_post_instantanes", &sauvegarde_post_instantanes_);
  param.ajouter_flag("lecture_post_instantanes", &lecture_post_instantanes_);
  param.ajouter_flag("lecture_post_instantanes_filtrer_u", &lecture_post_instantanes_filtrer_u_);
  param.ajouter_flag("lecture_post_instantanes_filtrer_rho", &lecture_post_instantanes_filtrer_rho_);
  param.ajouter_flag("lecture_post_instantanes_filtrer_p", &lecture_post_instantanes_filtrer_p_);
  param.ajouter_flag("lecture_post_instantanes_filtrer_tous", &lecture_post_instantanes_filtrer_tous_);
  param.ajouter("compteur_post_instantanes", &compteur_post_instantanes_);

  // DD,2016-06-15: changement des pas de temps de stabilite,
  param.ajouter_flag("old_dtstab", &old_dtstab_);

  // Oscillating boundary condition expression
  // param.ajouter("expression_oscillating_boundary", &expression_oscillating_boundary);
  param.ajouter("amplitude_oscillating_boundary", &amplitude_oscillating_boundary_);
  param.ajouter("frequency_oscillating_boundary", &frequency_oscillating_boundary_);

  /* param.ajouter("amplitude_oscillating_boundary_2", &amplitude_oscillating_boundary_2_); */
  /* param.ajouter("frequency_oscillating_boundary_2", &frequency_oscillating_boundary_2_); */


  param.lire_avec_accolades(is);

  // Recuperation des donnees de maillage
  splitting_ = ref_cast(IJK_Splitting, Interprete_bloc::objet_global(ijk_splitting_name));
  // Initialisation de delta_z_local_
  splitting_.get_local_mesh_delta(DIRECTION_K, 2 /* ghost cells */, delta_z_local_);

  if (geom_name_pour_delta == Nom("__identique_au_maillage__"))
    {
      facteur_delta_x_ = 1.;
      facteur_delta_y_ = 1.;
      delta_z_local_pour_delta_ = delta_z_local_;
    }
  else
    {
      Cerr << "Lecture de la taille du filtre depuis le IJK_Grid_Geometry " << geom_name_pour_delta << "." << finl;
      IJK_Grid_Geometry geom_pour_delta = ref_cast(IJK_Grid_Geometry, Interprete_bloc::objet_global(geom_name_pour_delta));
      const double dx_maillage = splitting_.get_grid_geometry().get_constant_delta(DIRECTION_I);
      const double dy_maillage = splitting_.get_grid_geometry().get_constant_delta(DIRECTION_J);
      const double dx_pour_delta = geom_pour_delta.get_constant_delta(DIRECTION_I);
      const double dy_pour_delta = geom_pour_delta.get_constant_delta(DIRECTION_J);
      facteur_delta_x_ = dx_pour_delta/dx_maillage;
      facteur_delta_y_ = dy_pour_delta/dy_maillage;

      calculer_delta_z_pour_delta(splitting_, geom_pour_delta, 2, delta_z_local_pour_delta_);
      for (int i=-2 ; i<splitting_.get_nb_elem_local(DIRECTION_K)+2 ; i++)
        {
          Cerr << i << "      " << delta_z_local_[i] << "        " << delta_z_local_pour_delta_[i] << finl;
        }
    }

  mode_terme_source_impose_ = (terme_source_acceleration_constant_ != 0);
  if ( (debit_cible_ != 0) && mode_terme_source_impose_ )
    {
      Cerr << " Error: (Inconsistent parameters) The keyword DEBIT_MASSIQUE cannot be used with TERME_SOURCE_ACCELERATION_CONdoubleANT." << finl;
      Process::exit();
    }
  else if (mode_terme_source_impose_)
    {
      Cerr << "Simulation en mode : TERME SOURCE IMPOSE." << finl;
    }
  else
    {
      Cerr << "Simulation en mode : DEBIT IMPOSE." << finl;
    }

  if ( (post_splitting_name != "??") || ( dt_post_spectral_ != -1 ) )
    {
      if ( dt_post_spectral_ <= 0 )
        {
          Cerr << " You need to specify an dt_post_spectral X," <<
               " were X is the number of step between 2 post" << finl;
          Process::exit();
        }
      else if ( (post_splitting_name == "??") )
        {
          Cerr << " You need to specify an other splitting like : spectral_splitting post_splitting " << finl;
          Cerr << "in order to compute FFT on planes you must only split in z direction" << finl;
          Process::exit();
        }
      else
        {
#if defined(WITH_FFTW)
          post_splitting_ =ref_cast(IJK_Splitting, Interprete_bloc::objet_global(post_splitting_name));

          const int decoupage_selon_x_ou_y = post_splitting_.get_nprocessor_per_direction(0)*post_splitting_.get_nprocessor_per_direction(1);
          if (decoupage_selon_x_ou_y > 1 )
            {
              Cerr << "In order to compute FFT on planes you must only split in z direction the spectral_splitting" << finl;
              exit();
            }
#else
          Cerr << " Module FFTW desactive impossible de faire la FFT " << finl;
          exit();
#endif
        }
    }

  if (sauvegarde_splitting_name_ != "??")
    {
      sauvegarde_splitting_ = ref_cast(IJK_Splitting, Interprete_bloc::objet_global(sauvegarde_splitting_name_));
      const int decoupage_selon_x_ou_y = sauvegarde_splitting_.get_nprocessor_per_direction(0)*sauvegarde_splitting_.get_nprocessor_per_direction(1);
      // Modif Martin
      if (decoupage_selon_x_ou_y > 1 )
        {
          Cerr << "In order to save lata on planes you must only split in z direction the sauvegarde_splitting" << finl;
          exit();
        }
      // Fin modif Martin

      Cerr << "Initialisation de la sauvegarde des lata par plan" << finl;
      redistribute_to_sauvegarde_splitting_elem_.initialize(splitting_,sauvegarde_splitting_,IJK_Splitting::ELEM);
      redistribute_to_sauvegarde_splitting_faces_[0].initialize(splitting_,sauvegarde_splitting_,IJK_Splitting::FACES_I);
      redistribute_to_sauvegarde_splitting_faces_[1].initialize(splitting_,sauvegarde_splitting_,IJK_Splitting::FACES_J);
      redistribute_to_sauvegarde_splitting_faces_[2].initialize(splitting_,sauvegarde_splitting_,IJK_Splitting::FACES_K);
    }

  if ( (!lecture_post_instantanes_) && (lecture_post_instantanes_filtrer_u_
                                        || lecture_post_instantanes_filtrer_rho_
                                        || lecture_post_instantanes_filtrer_p_
                                        || lecture_post_instantanes_filtrer_tous_) )
    {
      Cerr << " Error: (Inconsistent parameters) The keyword LECTURE_POdouble_INdoubleANTANES_FILTRER_U, LECTURE_POdouble_INdoubleANTANES_FILTRER_RHO, LECTURE_POdouble_INdoubleANTANES_FILTRER_P and LECTURE_POdouble_INdoubleANTANES_FILTRER_TOUS cannot be used without the keyword LECTURE_POdouble_INdoubleANTANES." << finl;
      Process::exit();
    }

  if (turbulent_viscosity_tensor_coefficients_.size_array() != 6)
    {
      Cerr << "Erreur: TURBULENT_VISCOSITY_TENSOR_COEFFICIENTS doit etre un vecteur de 6 composantes (xx, xy, xz, yy, yz, zz)." << finl;
      Process::exit();
    }
  if (turbulent_diffusivity_vector_coefficients_.size_array() != 3)
    {
      Cerr << "Erreur: TURBULENT_VISCOSITY_TENSOR_COEFFICIENTS doit etre un vecteur de 3 composantes (x, y, z)." << finl;
      Process::exit();
    }
  if (structural_uu_tensor_coefficients_.size_array() != 6)
    {
      Cerr << "Erreur: sTRUCTURAL_UU_TENSOR_COEFFICIENTS doit etre un vecteur de 6 composantes (xx, xy, xz, yy, yz, zz)." << finl;
      Process::exit();
    }
  if (structural_uscalar_vector_coefficients_.size_array() != 3)
    {
      Cerr << "Erreur: sTRUCTURAL_USCALAR_VECTOR_COEFFICIENTS doit etre un vecteur de 3 composantes (x, y, z)." << finl;
      Process::exit();
    }

  if ( large_eddy_simulation_formulation_ == Nom("favre") )
    {
      formulation_favre_ = 1;
    }
  else if ( large_eddy_simulation_formulation_ == Nom("velocity") )
    {
      formulation_velocity_ = 1;
    }

  // DD,2017-04-27: diffusion modifie en vue de l'ajout de modeles
  if ( large_eddy_simulation_formulation_ == Nom("none") )
    {
      if (turbulent_viscosity_)
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_VISCOSITY flag was activated but no large eddy simulation formulation was specified. "
               << "To specify a large eddy simulation formulation, you may use the keyword LARGE_EDDY_SIMULATION_FORMULATION with the parameters FAVRE, VELOCITY or NONE." << finl;
          Process::exit();
        }
      else if (turbulent_diffusivity_)
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_DIFFUSIVITY flag was activated but no large eddy simulation formulation was specified. "
               << "To specify a large eddy simulation formulation, you may use the keyword LARGE_EDDY_SIMULATION_FORMULATION with the parameters FAVRE, VELOCITY or NONE." << finl;
          Process::exit();
        }
      else if (structural_uu_)
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The sTRUCTURAL_UU flag was activated but no large eddy simulation formulation was specified. "
               << "To specify a large eddy simulation formulation, you may use the keyword LARGE_EDDY_SIMULATION_FORMULATION with the parameters FAVRE, VELOCITY or NONE." << finl;
          Process::exit();
        }
      else if (structural_uscalar_)
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The sTRUCTURAL_USCALAR flag was activated but no large eddy simulation formulation was specified. "
               << "To specify a large eddy simulation formulation, you may use the keyword LARGE_EDDY_SIMULATION_FORMULATION with the parameters FAVRE, VELOCITY or NONE." << finl;
          Process::exit();
        }
    }
  else
    {
      if (!turbulent_viscosity_ && !turbulent_diffusivity_ && !structural_uu_ && !structural_uscalar_)
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A LARGE_EDDY_SIMULATION_FORMULATION was specified but neither the TURBULENT_VISCOSITY flag, the TURBULENT_DIFFUSIVITY flag, the sTRUCTURAL_UU flag nor the sTRUCTURAL_USCALAR flag have been activated. "
               << "" << finl;
          Process::exit();
        }
    }

  if (!turbulent_viscosity_)
    {
      if ( type_velocity_turbulent_diffusion_ != Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A velocity turbulent diffusion type was specified but the TURBULENT_VISCOSITY flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
      else if ( turbulent_viscosity_model_ != Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A subgrid viscosity model was specified but the TURBULENT_VISCOSITY flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
      else if ( turbulent_viscosity_model_constant_ != -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A subgrid viscosity model constant was specified but the TURBULENT_VISCOSITY flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
    }
  else
    {
      if ( type_velocity_turbulent_diffusion_ == Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_VISCOSITY flag expect a velocity turbulent diffusion type but no velocity turbulent diffusion type was specified. "
               << "To specify a velocity turbulent diffusion type, you may use the keyword TYPE_VELOCITY_TURBULENT_DIFFUSION with the parameters SIMPLE, SIMPLE_WITH_TRANSPOSE, FULL or NONE" << finl;
          Process::exit();
        }
      else if ( turbulent_viscosity_model_ == Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_VISCOSITY flag expect a subgrid viscosity model but no model was specified. "
               << "To specify a model, you may use the keyword TURBULENT_VISCOSITY_MODEL with the parameters CONdoubleANT, UNSRHO, SMAGORINSKY, VREMAN, SIGMA, WALE, AMD, AMD_COMP, AMDNOCLIP, AMDSCALAR, AMDSCALARNOCLIP, RDS, VSS, KOBAYASHI or NONE." << finl;
          Process::exit();
        }
      else if ( turbulent_viscosity_model_constant_ == -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_VISCOSITY flag expect a subgrid viscosity model constant but no constant was specified. "
               << "To specify a constant, you may use the keyword TURBULENT_VISCOSITY_MODEL_CONdoubleANT." << finl;
          Process::exit();
        }
    }

  if (!turbulent_diffusivity_)
    {
      if ( type_scalar_turbulent_diffusion_ != Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A scalar turbulent diffusion type was specified but the TURBULENT_DIFFUSIVITY flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
      if ( turbulent_diffusivity_model_ != Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A subgrid diffusivity model was specified but the TURBULENT_DIFFUSIVITY flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
      if ( turbulent_diffusivity_model_constant_ != -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A subgrid diffusivity model constant was specified but the TURBULENT_DIFFUSIVITY flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
    }
  else
    {
      if ( type_scalar_turbulent_diffusion_ == Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_DIFFUSIVITY flag expect a scalar turbulent diffusion type but no scalar turbulent diffusion type was specified. "
               << "To specify a scalar turbulent diffusion type, you may use the keyword TYPE_SCALAR_TURBULENT_DIFFUSION with the parameters NORMAL, ANISOTROPIC or NONE" << finl;
          Process::exit();
        }
      if ( turbulent_diffusivity_model_ == Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_DIFFUSIVITY flag expect a subgrid diffusivity model but no model was specified. "
               << "To specify a model, you may use the keyword TURBULENT_DIFFUSIVITY_MODEL with the parameters CONdoubleANT, UNSRHO, SMAGORINSKY, VREMAN, SIGMA, WALE, AMD, AMD_COMP, AMDNOCLIP, AMDSCALAR, AMDSCALARNOCLIP, RDS, VSS, KOBAYASHI or NONE." << finl;
          Process::exit();
        }
      else if ( turbulent_diffusivity_model_constant_ == -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The TURBULENT_DIFFUSIVITY flag expect a subgrid diffusivity model constant but no constant was specified. "
               << "To specify a constant, you may use the keyword TURBULENT_DIFFUSIVITY_MODEL_CONdoubleANT." << finl;
          Process::exit();
        }
    }

  if (!structural_uu_)
    {
      if ( structural_uu_model_ != Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A structural model for uu was specified but the sTRUCTURAL_UU flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
      else if ( structural_uu_model_constant_ != -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A structural model for uu constant was specified but the sTRUCTURAL_UU flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
    }
  else
    {
      if ( structural_uu_model_ == Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The sTRUCTURAL_UU flag expect a structural model but no model was specified. "
               << "To specify a model, you may use the keyword sTRUCTURAL_UU_MODEL with the parameters GRADIENT, GRADIENT_FILTRE, SU_LAPLACIEN_U, SIMILARITY or NONE." << finl;
          Process::exit();
        }
      else if ( structural_uu_model_constant_ == -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The sTRUCTURAL_UU flag expect a structural model constant but no constant was specified. "
               << "To specify a constant, you may use the keyword sTRUCTURAL_UU_MODEL_CONdoubleANT." << finl;
          Process::exit();
        }
    }

  if (!structural_uscalar_)
    {
      if ( structural_uscalar_model_ != Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A structural model for uscalar was specified but the sTRUCTURAL_USCALAR flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
      if ( structural_uscalar_model_constant_ != -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "A structural model for uscalar constant was specified but the sTRUCTURAL_USCALAR flag has not been activated. "
               << "" << finl;
          Process::exit();
        }
    }
  else
    {
      if ( structural_uscalar_model_ == Nom("none") )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The sTRUCTURAL_USCALAR flag expect a subgrid diffusivity model but no model was specified. "
               << "To specify a model, you may use the keyword sTRUCTURAL_USCALAR_MODEL with the parameters GRADIENT, GRADIENT_FILTRE, SIMILARITY or NONE." << finl;
          Process::exit();
        }
      else if ( structural_uscalar_model_constant_ == -1 )
        {
          Cerr << "Error: (Inconsistent parameters) "
               << "The sTRUCTURAL_USCALAR flag expect a subgrid diffusivity model constant but no constant was specified. "
               << "To specify a constant, you may use the keyword sTRUCTURAL_USCALAR_MODEL_CONdoubleANT." << finl;
          Process::exit();
        }
    }

  // Note :
  // 'dynamic_uu', 'dynamicglobal_uu', 'dynamic_urho', 'dynamicglobal_urho',
  // 'dynamic_ut' et 'dynamicglobal_ut' sont des options obsoletes
  // equivalentes a 'lilly' et 'lillyglobal'.

  if (turbulent_viscosity_dynamic_type_ == Nom("dynamic_uu"))
    {
      turbulent_viscosity_dynamic_type_ = Nom("lilly");
    }
  if (turbulent_viscosity_dynamic_type_ == Nom("dynamicglobal_uu"))
    {
      turbulent_viscosity_dynamic_type_ = Nom("lillyglobal");
    }
  if ( turbulent_diffusivity_dynamic_type_ == Nom("dynamic_urho")
       || turbulent_diffusivity_dynamic_type_ == Nom("dynamic_ut") )
    {
      turbulent_diffusivity_dynamic_type_ = Nom("lilly");
    }
  if ( turbulent_diffusivity_dynamic_type_ == Nom("dynamicglobal_urho")
       || turbulent_diffusivity_dynamic_type_ == Nom("dynamicglobal_ut") )
    {
      turbulent_diffusivity_dynamic_type_ = Nom("lillyglobal");
    }

  flag_u_filtre_ = lecture_post_instantanes_filtrer_u_ || lecture_post_instantanes_filtrer_tous_
                   || (structural_uu_ && (structural_uu_model_ != Nom("gradient")) )
                   || (structural_uscalar_ && (structural_uscalar_model_ != Nom("gradient")) )
                   || (turbulent_viscosity_dynamic_type_ != Nom("not_dynamic"))
                   || (turbulent_diffusivity_dynamic_type_ != Nom("not_dynamic"))
                   || (structural_uu_dynamic_type_ != Nom("not_dynamic"))
                   || (structural_uscalar_dynamic_type_ != Nom("not_dynamic"));
  flag_rho_filtre_ = lecture_post_instantanes_filtrer_rho_ || lecture_post_instantanes_filtrer_tous_
                     || (formulation_velocity_ && structural_uscalar_ && (structural_uscalar_model_ != Nom("gradient")) )
                     || (turbulent_diffusivity_dynamic_type_ != Nom("not_dynamic"))
                     || (turbulent_viscosity_dynamic_type_ != Nom("not_dynamic"))
                     || (structural_uu_dynamic_type_ != Nom("not_dynamic"))
                     || (structural_uscalar_dynamic_type_ != Nom("not_dynamic"));
  flag_temperature_filtre_ = formulation_favre_ && (
                               (structural_uscalar_ && (structural_uscalar_model_ != Nom("gradient")) )
                               || (turbulent_viscosity_dynamic_type_ != Nom("not_dynamic"))
                               || (turbulent_diffusivity_dynamic_type_ != Nom("not_dynamic"))
                               || (structural_uu_dynamic_type_ != Nom("not_dynamic"))
                               || (structural_uscalar_dynamic_type_ != Nom("not_dynamic"))
                             );
  flag_turbulent_mu_filtre_ = turbulent_viscosity_ && ((turbulent_viscosity_dynamic_type_ != Nom("not_dynamic")) || (structural_uu_dynamic_type_ != Nom("not_dynamic")));
  flag_turbulent_kappa_filtre_ = turbulent_diffusivity_ && ((turbulent_diffusivity_dynamic_type_ != Nom("not_dynamic")) || (structural_uscalar_dynamic_type_ != Nom("not_dynamic")));
  flag_structural_uu_filtre_ = structural_uu_ && ((turbulent_viscosity_dynamic_type_ != Nom("not_dynamic")) || (structural_uu_dynamic_type_ != Nom("not_dynamic")));
  flag_structural_uscalar_filtre_ = structural_uscalar_ && ((turbulent_diffusivity_dynamic_type_ != Nom("not_dynamic")) || (structural_uscalar_dynamic_type_ != Nom("not_dynamic")));

  flag_structural_uu_tmp_ = structural_uu_ && (structural_uu_model_ == Nom("gradient_filtre") || structural_uu_model_ == Nom("convection_filtre"));
  flag_structural_uscalar_tmp_ = structural_uscalar_ && (structural_uscalar_model_ == Nom("gradient_filtre"));
  flag_d_velocity_tmp_ = flag_filtrage_convection_qdm_ || flag_filtrage_turbulent_diffusion_qdm_ || flag_filtrage_structural_diffusion_qdm_;

  flag_nu_anisotropic_ = ( (type_velocity_turbulent_diffusion_ == Nom("simple_anisotropic"))
                           || (type_velocity_turbulent_diffusion_ == Nom("simple_with_transpose_anisotropic"))
                           || (type_velocity_turbulent_diffusion_ == Nom("full_anisotropic")) );
  flag_nu_tensorial_ = ( turbulent_viscosity_dynamic_type_.finit_par("tensorial")
                         || turbulent_viscosity_dynamic_type_.finit_par("tensorial_twopass")
                         || (turbulent_viscosity_tensor_coefficients_[0] != 1.)
                         || (turbulent_viscosity_tensor_coefficients_[1] != 1.)
                         || (turbulent_viscosity_tensor_coefficients_[2] != 1.)
                         || (turbulent_viscosity_tensor_coefficients_[3] != 1.)
                         || (turbulent_viscosity_tensor_coefficients_[4] != 1.)
                         || (turbulent_viscosity_tensor_coefficients_[5] != 1.) );
  flag_kappa_anisotropic_ = (type_scalar_turbulent_diffusion_ == Nom("anisotropic"));
  flag_kappa_vectorial_ = ( turbulent_diffusivity_dynamic_type_.finit_par("vectorial")
                            || turbulent_diffusivity_dynamic_type_.finit_par("vectorial_twopass")
                            || (turbulent_diffusivity_vector_coefficients_[0] != 1.)
                            || (turbulent_diffusivity_vector_coefficients_[1] != 1.)
                            || (turbulent_diffusivity_vector_coefficients_[2] != 1.) );

  // //3/03/14 changement de la source
  // // Preparation de l'expression derivee de l'acceleration
  // String2 tmpstring(expression_derivee_acceleration_);
  // parser_derivee_acceleration_.setString(tmpstring);
  // parser_derivee_acceleration_.setNbVar(3);
  // parser_derivee_acceleration_.addVar("force");
  // parser_derivee_acceleration_.addVar("v_moyen");
  // parser_derivee_acceleration_.addVar("rho_v_moyen");
  // parser_derivee_acceleration_.parseString();
  //

  flag_oscillating_boundary=false;
  if(amplitude_oscillating_boundary_ && frequency_oscillating_boundary_)
    {
      Cerr <<
           "The boundary is oscillating at amplitude "
           << amplitude_oscillating_boundary_
           <<  " and frequency "
           << frequency_oscillating_boundary_
           << finl;
      flag_oscillating_boundary = true;
    }


  flag_t_bulk_forced_=false;
  if(aim_t_bulk_ != -1)
    {
      Cerr <<
           "Simulation in forced T bulk mode with value "
           << aim_t_bulk_
           << " K"
           << finl;
      flag_t_bulk_forced_ = true;
    }



  if (nom_reprise_ != "??")
    reprendre_qc(nom_reprise_);

  run();
  return is;
}

// Attention, il faut 1 epaisseur de joint valide sur rho
static void calculer_debit(const IJK_Field_double& vx, const IJK_Field_double& rho, const ArrOfDouble_with_ghost& delta_z, const double longeur_Lx_tot, double& debit)
{
  const int ni = vx.ni();
  const int nj = vx.nj();
  const int nk = vx.nk();
  const double dx = vx.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = vx.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);

  debit = 0.;
  // Ponderation par le volume des mailles pour traiter le cas ou le maillage est irregulier en k
  for (int k = 0; k < nk; k++)
    {
      double somme_rhov = 0;
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              double rhov = vx(i,j,k) * (rho(i-1,j,k) + rho(i,j,k)) * 0.5; // rho moyen sur la face
              somme_rhov += rhov;
              //Cerr << "somme " << somme_rhov << " i " << i << " j " << j << " k " << k << " vx " << vx(i,j,k) << " rho -1 " << rho(i-1,j,k) << " rho " << rho(i,j,k) << finl;
            }
        }
      // volume d'une maille de ce plan
      double volume = dx * dy * delta_z[k];
      debit += volume * somme_rhov;
    }
  // somme sur tous les processeurs. On regroupe toutes les communications en un seul appel:
  ArrOfDouble tmp(1);
  tmp[0] = debit;
  mp_sum_for_each_item(tmp);

  // ici on calcule la valeur moyenne sur le volume V = Lx * Ly * Lz et on multiplie par la surface debitante S = Ly * Lz
// on doit donc juste diviser par Lx
  debit = tmp[0] / longeur_Lx_tot;
}

// Methode appelee dans run().
// Intialise toutes les variables et champs pour le demarrage du calcul
// (soit avec des valeurs du jdd, soit avec des valeurs lues dans le fichier de reprise)
void DNS_QC_double::initialise()
{
  Cout << "DNS_QC_double::initialise()" << finl;

  // precision de sortie du Cerr !!!
  Cerr.setf(ios::scientific);
  Cerr.precision(20);
  Cout.setf(ios::scientific);
  Cout.precision(20);

  Cout << "precision de sortie du out et du err modifier" << finl;

  // calcul de la constante specifique du gaz
  constante_specifique_gaz_ = Cp_gaz_ *(1.-1./gamma_);
  Cout << " Cp_gaz = " << Cp_gaz_
       << "\n constante_specifique_gaz = " << constante_specifique_gaz_ << finl;

  // calcule les valeurs aux parois
  const double Pth_sur_R = P_thermodynamique_ / constante_specifique_gaz_;
  rho_paroi_impose_kmin_ = Pth_sur_R / T_paroi_impose_kmin_;
  rho_paroi_impose_kmax_ = Pth_sur_R / T_paroi_impose_kmax_;
  lambda_de_t_paroi_kmin_ = calculer_lambda_air(T_paroi_impose_kmin_);
  lambda_de_t_paroi_kmax_ = calculer_lambda_air(T_paroi_impose_kmax_);
  Cout << " T_paroi_impose_kmin_ = " << T_paroi_impose_kmin_
       << "\n T_paroi_impose_kmax_ = " << T_paroi_impose_kmax_
       << "\n lambda_de_t_paroi_kmin = " << lambda_de_t_paroi_kmin_
       << "\n lambda_de_t_paroi_kmax = " << lambda_de_t_paroi_kmax_ << finl;

  // calcul du volume total du domaine
  const IJK_Grid_Geometry& geometry = splitting_.get_grid_geometry();
  const Nom& geom_name = splitting_.get_grid_geometry().le_nom();

  int Nx_tot = geometry.get_nb_elem_tot(0) + 1;
  int Ny_tot = geometry.get_nb_elem_tot(1) + 1;
  int Nz_tot = geometry.get_nb_elem_tot(2) + 1;
  Lx_tot_  = geometry.get_node_coordinates(0)[Nx_tot - 1] - geometry.get_origin(0);
  Ly_tot_  = geometry.get_node_coordinates(1)[Ny_tot - 1] - geometry.get_origin(1);
  Lz_tot_  = geometry.get_node_coordinates(2)[Nz_tot - 1] - geometry.get_origin(2);
  volume_total_domaine_ = Lx_tot_ * Ly_tot_ * Lz_tot_;
  Cout << " Volume total domaine = " << volume_total_domaine_ << finl;

  if ( fichier_reprise_rho_ == "??" )   // si on ne fait pas une reprise on initialise T et rho suivant une methode
    {
      Cout << "Initialisation temperature, expression = " << expression_temperature_initiale_ << finl;
      set_field_data(temperature_, expression_temperature_initiale_);
      // Calcul de rho_ fonction de la temperature:
      const int nk = temperature_.nk();
      const int nj = temperature_.nj();
      const int ni = temperature_.ni();
      for (int k = 0; k < nk; k++ )
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  double t = temperature_(i,j,k);
                  rho_(i,j,k) = Pth_sur_R / t;
                }
            }
        }
    }
  else
    {
      Cout << "Lecture rho initial dans fichier " << fichier_reprise_rho_ << " timestep= " << timestep_reprise_rho_ << finl;
      lire_dans_lata(fichier_reprise_rho_, timestep_reprise_rho_, geom_name, "RHO", rho_);
    }

  if (fichier_reprise_vitesse_ == "??")   // si on ne fait pas une reprise on initialise V
    {
      if (expression_vitesse_initiale_.size() != 3)
        {
          Cerr << "Erreur dans l'initialisation: la vitesse initiale doit etre fournie avec trois expressions" << finl;
          Process::exit();
        }
      Cout << "Initialisation vitesse \nvx = " << expression_vitesse_initiale_[0]
           << "\nvy = " <<  expression_vitesse_initiale_[1]
           << "\nvz = " <<  expression_vitesse_initiale_[2]  << finl;
      for (int i = 0; i < 3; i++)
        set_field_data(velocity_[i], expression_vitesse_initiale_[i]);
    }
  else
    {
      Cout << "Lecture vitesse initiale dans fichier " << fichier_reprise_vitesse_ << " timestep= " << timestep_reprise_vitesse_ << finl;
      lire_dans_lata(fichier_reprise_vitesse_, timestep_reprise_vitesse_, geom_name, "VELOCITY",
                     velocity_[0], velocity_[1], velocity_[2]); // fonction qui lit un champ a partir d'un lata .
    }

  // initialise temperature, lambda et mu
  // (on a besoin de mu(rho) initial pour l'operateur diffusion, voir rk3_sub_step()
  rho_.echange_espace_virtuel(1); // Nous avons besoin d'avoir rho a jour aux bord pour calculer le joints de mu et lambda !!
  calculer_temperature_mu_lambda_air(P_thermodynamique_, constante_specifique_gaz_, rho_, temperature_, molecular_mu_, molecular_lambda_, 1);

  // statistiques...
  Cout << "Initialisation des statistiques. T_debut_statistiques=" << t_debut_statistiques_ << finl;
  statistiques_.initialize(splitting_.get_grid_geometry(), T_paroi_impose_kmax_,T_paroi_impose_kmin_, constante_specifique_gaz_ );

  if ( dt_post_spectral_ > 0 ) // si on a un post traitement spectral
    {
      if (!statistiques_.is_converge())
        {
          Cerr << " les stats ne sont pas converger pas de productions spectrale allowed." << finl;
          exit();
        }

      partie_fourier_.initialize(splitting_ /* maillage de simu */
                                 ,post_splitting_ /* maillage de post traitement */
                                 ,statistiques_.vitesse_moyenne() /* vitesse moyenne utilisee dans les stats standart */
                                 ,statistiques_.masse_volumique_moyenne() /* masse volumique moyenne utilisee dans les stats standart */
                                 ,statistiques_.viscosite_cinematique_moyenne() /* viscosite cinematique moyenne utilisee dans les stats standart */
                                 ,T_paroi_impose_kmax_
                                 ,T_paroi_impose_kmin_
                                 ,reprise_fourier_);
    }

// pour source :
  debit_actuel_ = debit_cible_;
  if (debit_actuel_<0)
    {
      calculer_debit(velocity_[0], rho_, delta_z_local_, Lx_tot_, debit_actuel_);
      debit_cible_=debit_actuel_;
      Cerr<< "debit "<<debit_actuel_<<finl;
    }

  if(flag_t_bulk_forced_)
    actual_t_bulk_ = aim_t_bulk_;
}

static void force_zero_normal_velocity_on_walls(IJK_Field_double& vz)
{
  const int nj = vz.nj();
  const int ni = vz.ni();
  const int kmin = vz.get_splitting().get_offset_local(DIRECTION_K);
  const int nktot = vz.get_splitting().get_nb_items_global(IJK_Splitting::FACES_K, DIRECTION_K);
  if (kmin == 0)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              vz(i,j,0) = 0.;
            }
        }
    }
  if (kmin + vz.nk() == nktot)
    {
      const int k = vz.nk()-1;
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              vz(i,j,k) = 0.;
            }
        }
    }
}

// // Attention, il faut 1 epaisseur de joint valide sur rho
// static void calculer_v_et_rhov_moyen(const IJK_Field_double & vx, const IJK_Field_double & rho, const ArrOfDouble_with_ghost & delta_z, const double volume_total, double & v_moy, double & rho_v_moy)
// {
//   const int ni = vx.ni();
//   const int nj = vx.nj();
//   const int nk = vx.nk();
//   const double dx = vx.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
//   const double dy = vx.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);
//
//   v_moy = 0.;
//   rho_v_moy = 0.;
//
//   // Ponderation par le volume des mailles pour traiter le cas ou le maillage est irregulier en k
//   for (int k = 0; k < nk; k++) {
//     double somme_v = 0.;
//     double somme_rhov = 0;
//     for (int j = 0; j < nj; j++) {
//       for (int i = 0; i < ni; i++) {
//      double v = vx(i,j,k);
//      double rhov = v * (rho(i-1,j,k) + rho(i,j,k)) * 0.5; // rho moyen sur la face
//      somme_v += v;
//      somme_rhov += rhov;
//       }
//     }
//     // volume d'une maille de ce plan
//     double volume = dx * dy * delta_z[k];
//     v_moy += volume * somme_v;
//     rho_v_moy += volume * somme_rhov;
//   }
//   // somme sur tous les processeurs. On regroupe toutes les communications en un seul appel:
//   ArrOfDouble tmp(2);
//   tmp[0] = v_moy;
//   tmp[1] = rho_v_moy;
//   mp_sum_for_each_item(tmp);
//
//   v_moy = tmp[0] / volume_total;
//   rho_v_moy = tmp[1] / volume_total;
// }


static double calculer_dtstab_diffusion_temperature_local(const IJK_Field_double& lambda,
                                                          const IJK_Field_double& rho,
                                                          const double cp_gaz)
{
  const IJK_Splitting& split = lambda.get_splitting();
  const int ni = split.get_nb_elem_local(DIRECTION_I);
  const int nj = split.get_nb_elem_local(DIRECTION_J);
  const int nk = split.get_nb_elem_local(DIRECTION_K);
  const IJK_Grid_Geometry& geom = split.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const ArrOfDouble& delta_z = geom.get_delta(DIRECTION_K);
  const int k_offset = split.get_offset_local(DIRECTION_K);

  double inv_dtstab = 1e-20;
  for (int k = 0; k < nk; k++)
    {
      const double dz = delta_z[k + k_offset];
      const double coeff = 2. * (1./(dx*dx) + 1./(dy*dy) + 1./(dz*dz)) / cp_gaz;
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              double dt = lambda(i,j,k) / rho(i,j,k) * coeff;
              inv_dtstab = max(dt, inv_dtstab);
            }
        }
    }
  return 1. / inv_dtstab;
}

static double calculer_dtstab_diffusion_temperature_local_sans_rho(const bool anisotropic,
                                                                   const IJK_Field_double& lambda_x,
                                                                   const IJK_Field_double& lambda_y,
                                                                   const IJK_Field_double& lambda_z,
                                                                   const double cp_gaz)
{
  const IJK_Splitting& split = lambda_x.get_splitting();
  const int ni = split.get_nb_elem_local(DIRECTION_I);
  const int nj = split.get_nb_elem_local(DIRECTION_J);
  const int nk = split.get_nb_elem_local(DIRECTION_K);
  const IJK_Grid_Geometry& geom = split.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const ArrOfDouble& delta_z = geom.get_delta(DIRECTION_K);
  const int k_offset = split.get_offset_local(DIRECTION_K);

  double inv_dtstab = 1e-20;
  for (int k = 0; k < nk; k++)
    {
      const double dz = delta_z[k + k_offset];
      double coeff_x;
      double coeff_y;
      double coeff_z;
      if (anisotropic)
        {
          coeff_x = 2. * (1./(dx)) / cp_gaz;
          coeff_y = 2. * (1./(dy)) / cp_gaz;
          coeff_z = 2. * (1./(dz)) / cp_gaz;
        }
      else
        {
          coeff_x = 2. * (1./(dx*dx)) / cp_gaz;
          coeff_y = 2. * (1./(dy*dy)) / cp_gaz;
          coeff_z = 2. * (1./(dz*dz)) / cp_gaz;
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              double dt_x = lambda_x(i,j,k) * coeff_x;
              double dt_y = lambda_y(i,j,k) * coeff_y;
              double dt_z = lambda_z(i,j,k) * coeff_z;
              double dt = dt_x + dt_y + dt_z;
              inv_dtstab = max(dt, inv_dtstab);
            }
        }
    }
  return 1. / inv_dtstab;
}

static double calculer_dtstab_diffusion_temperature_local_avec_turbulent_favre(const bool anisotropic,
                                                                               const IJK_Field_double& lambda,
                                                                               const IJK_Field_double& lambda_turbulent_xx,
                                                                               const IJK_Field_double& lambda_turbulent_xy,
                                                                               const IJK_Field_double& lambda_turbulent_xz,
                                                                               const IJK_Field_double& lambda_turbulent_yy,
                                                                               const IJK_Field_double& lambda_turbulent_yz,
                                                                               const IJK_Field_double& lambda_turbulent_zz,
                                                                               const IJK_Field_double& rho,
                                                                               const double cp_gaz)
{
  const IJK_Field_double& lambda_turbulent_yx = lambda_turbulent_xy;
  const IJK_Field_double& lambda_turbulent_zx = lambda_turbulent_xz;
  const IJK_Field_double& lambda_turbulent_zy = lambda_turbulent_yz;

  const IJK_Splitting& split = lambda.get_splitting();
  const int ni = split.get_nb_elem_local(DIRECTION_I);
  const int nj = split.get_nb_elem_local(DIRECTION_J);
  const int nk = split.get_nb_elem_local(DIRECTION_K);
  const IJK_Grid_Geometry& geom = split.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const ArrOfDouble& delta_z = geom.get_delta(DIRECTION_K);
  const int k_offset = split.get_offset_local(DIRECTION_K);

  double inv_dtstab = 1e-20;
  for (int k = 0; k < nk; k++)
    {
      const double dz = delta_z[k + k_offset];
      const double coeff_x = 2. * (1./(dx*dx)) / cp_gaz;
      const double coeff_y = 2. * (1./(dy*dy)) / cp_gaz;
      const double coeff_z = 2. * (1./(dz*dz)) / cp_gaz;
      double fac_x;
      double fac_y;
      double fac_z;
      if (anisotropic)
        {
          fac_x = dx;
          fac_y = dy;
          fac_z = dz;
        }
      else
        {
          fac_x = 1.;
          fac_y = 1.;
          fac_z = 1.;
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              double dt_xx = ( lambda(i,j,k) + lambda_turbulent_xx(i,j,k)*fac_x ) / rho(i,j,k) * coeff_x;
              double dt_xy = ( lambda(i,j,k) + lambda_turbulent_xy(i,j,k)*fac_y ) / rho(i,j,k) * coeff_y;
              double dt_xz = ( lambda(i,j,k) + lambda_turbulent_xz(i,j,k)*fac_z ) / rho(i,j,k) * coeff_z;
              double dt_yx = ( lambda(i,j,k) + lambda_turbulent_yx(i,j,k)*fac_x ) / rho(i,j,k) * coeff_x;
              double dt_yy = ( lambda(i,j,k) + lambda_turbulent_yy(i,j,k)*fac_y ) / rho(i,j,k) * coeff_y;
              double dt_yz = ( lambda(i,j,k) + lambda_turbulent_yz(i,j,k)*fac_z ) / rho(i,j,k) * coeff_z;
              double dt_zx = ( lambda(i,j,k) + lambda_turbulent_zx(i,j,k)*fac_x ) / rho(i,j,k) * coeff_x;
              double dt_zy = ( lambda(i,j,k) + lambda_turbulent_zy(i,j,k)*fac_y ) / rho(i,j,k) * coeff_y;
              double dt_zz = ( lambda(i,j,k) + lambda_turbulent_zz(i,j,k)*fac_z ) / rho(i,j,k) * coeff_z;
              double dt_x = dt_xx + dt_xy + dt_xz;
              double dt_y = dt_yx + dt_yy + dt_yz;
              double dt_z = dt_zx + dt_zy + dt_zz;
              double dt = max(dt_x , max(dt_y, dt_z));
              inv_dtstab = max(dt, inv_dtstab);
            }
        }
    }
  return 1. / inv_dtstab;
}

static double calculer_dtstab_diffusion_temperature_local_avec_turbulent_velocity(const bool anisotropic,
                                                                                  const IJK_Field_double& lambda,
                                                                                  const IJK_Field_double& lambda_turbulent_xx,
                                                                                  const IJK_Field_double& lambda_turbulent_xy,
                                                                                  const IJK_Field_double& lambda_turbulent_xz,
                                                                                  const IJK_Field_double& lambda_turbulent_yy,
                                                                                  const IJK_Field_double& lambda_turbulent_yz,
                                                                                  const IJK_Field_double& lambda_turbulent_zz,
                                                                                  const IJK_Field_double& rho,
                                                                                  const double cp_gaz)
{
  const IJK_Field_double& lambda_turbulent_yx = lambda_turbulent_xy;
  const IJK_Field_double& lambda_turbulent_zx = lambda_turbulent_xz;
  const IJK_Field_double& lambda_turbulent_zy = lambda_turbulent_yz;

  const IJK_Splitting& split = lambda.get_splitting();
  const int ni = split.get_nb_elem_local(DIRECTION_I);
  const int nj = split.get_nb_elem_local(DIRECTION_J);
  const int nk = split.get_nb_elem_local(DIRECTION_K);
  const IJK_Grid_Geometry& geom = split.get_grid_geometry();
  const double dx = geom.get_constant_delta(DIRECTION_I);
  const double dy = geom.get_constant_delta(DIRECTION_J);
  const ArrOfDouble& delta_z = geom.get_delta(DIRECTION_K);
  const int k_offset = split.get_offset_local(DIRECTION_K);

  double inv_dtstab = 1e-20;
  for (int k = 0; k < nk; k++)
    {
      const double dz = delta_z[k + k_offset];
      const double coeff_x = 2. * (1./(dx*dx)) / cp_gaz;
      const double coeff_y = 2. * (1./(dy*dy)) / cp_gaz;
      const double coeff_z = 2. * (1./(dz*dz)) / cp_gaz;
      double fac_x;
      double fac_y;
      double fac_z;
      if (anisotropic)
        {
          fac_x = dx;
          fac_y = dy;
          fac_z = dz;
        }
      else
        {
          fac_x = 1.;
          fac_y = 1.;
          fac_z = 1.;
        }
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              double dt_xx = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_xx(i,j,k)*fac_x ) * coeff_x;
              double dt_xy = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_xy(i,j,k)*fac_y ) * coeff_y;
              double dt_xz = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_xz(i,j,k)*fac_z ) * coeff_z;
              double dt_yx = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_yx(i,j,k)*fac_x ) * coeff_x;
              double dt_yy = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_yy(i,j,k)*fac_y ) * coeff_y;
              double dt_yz = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_yz(i,j,k)*fac_z ) * coeff_z;
              double dt_zx = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_zx(i,j,k)*fac_x ) * coeff_x;
              double dt_zy = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_zy(i,j,k)*fac_y ) * coeff_y;
              double dt_zz = ( lambda(i,j,k)/rho(i,j,k) + lambda_turbulent_zz(i,j,k)*fac_z ) * coeff_z;
              double dt_x = dt_xx + dt_xy + dt_xz;
              double dt_y = dt_yx + dt_yy + dt_yz;
              double dt_z = dt_zx + dt_zy + dt_zz;
              double dt = max(dt_x , max(dt_y, dt_z));
              inv_dtstab = max(dt, inv_dtstab);
            }
        }
    }
  return 1. / inv_dtstab;
}

// Impression dans cerr du flux vertical (direction K) conductif et
void DNS_QC_double::calculer_moyennes_flux()
{
  const int ni = velocity_[2].ni();
  const int nj = velocity_[2].nj();
  const int nk = velocity_[2].nk();

  const double dx = velocity_[2].get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = velocity_[2].get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);

  const double facteur = Cp_gaz_ * P_thermodynamique_ / constante_specifique_gaz_;
  const double surface = 1. / (ni*nj); // calcule directement la surface et la moyenne.
  const int offset = temperature_.get_splitting().get_offset_local(DIRECTION_K);
  const int nktot = velocity_[2].get_splitting().get_nb_items_global(IJK_Splitting::FACES_K, DIRECTION_K);

  // tableau globaux plus facile pour faire un mp_sum.
  ArrOfDouble flux_cv(nktot);
  ArrOfDouble flux_cd(nktot);

  flux_cv =0;
  flux_cd =0;

  for (int k = 0; k < nk; k++)
    {
      double partiel_flux_cv = 0.;
      double partiel_flux_cd = 0.;

      if ( (k + offset) == (nktot-1) )
        {
          partiel_flux_cv =0; // vitesse nulle a la paroi
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  partiel_flux_cd += boundary_flux_kmax_(i,j,0)/ (dx* dy); // somme des flux de bord
                }
            }
        }
      else if ( (k + offset) == 0 )
        {
          partiel_flux_cv =0; // vitesse nulle a la paroi
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  partiel_flux_cd += boundary_flux_kmin_(i,j,0) / (dx* dy); // somme des flux de bord
                }
            }
        }
      else
        {
          const double d0 = delta_z_local_[k-1] * 0.5;
          const double d1 = delta_z_local_[k]   * 0.5;
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  const double v = velocity_[2](i,j,k);
                  const double L0 = molecular_lambda_(i,j,k-1);
                  const double L1 = molecular_lambda_(i,j,k);

                  const double lambda_face = (L0 * L1) / ( d0 * L0 + d1 *L1);
                  // mauvaise methode ici il faut utiliser la CL de flux bord.
                  const double T_inf = temperature_(i,j,k-1);
                  const double T_sup = temperature_(i,j,k);

                  double flux = lambda_face * (  T_inf - T_sup); // on calcule -L dT/dy

                  partiel_flux_cv += v * facteur;
                  partiel_flux_cd += flux;
                }
            }
        }
      flux_cv[k+offset] = partiel_flux_cv;
      flux_cd[k+offset] = partiel_flux_cd;
    }
  flux_cv *= surface;
  flux_cd *= surface;

  // communication
  // je veux ton recevoir sur le 0 pour ecrire.
  // pas besoin d'un envoi general.
  // bon pour le moment je m'en contente.
  mp_sum_for_each_item(flux_cv);
  mp_sum_for_each_item(flux_cd);

  // ici on va imprimmer ce dont on a besoin.
#if 0
  const IJK_Grid_Geometry& geometry = splitting_.get_grid_geometry();
  if(Process::je_suis_maitre())
    {
      for (int k = 0; k < nktot; k++)
        {
          double z = geometry.get_node_coordinates(2)[k];
          Cerr << "iiii " << z << " " << flux_cv[k] << " " << flux_cd[k] << finl;
        }
      Cerr << "iiii " << finl;
    }
#endif
}


void DNS_QC_double::ecrire_fichier_sauv(const char *fichier_sauvegarde, const char *lata_name)
{
  if (Process::je_suis_maitre())
    {
      Cerr << "T= " << current_time_ << " Checkpointing dans le fichier " << fichier_sauvegarde << finl;
      SFichier fichier(fichier_sauvegarde);
      fichier.precision(17);
      fichier.setf(std::ios_base::scientific);
      fichier << "{\n"
              << " tinit " << current_time_ << "\n"
              << " terme_acceleration_init " << terme_source_acceleration_ << "\n"
              << " p_thermo_init " << P_thermodynamique_ << "\n"
              << " fichier_reprise_vitesse " << lata_name << "\n";
      fichier << " fichier_reprise_rho " << lata_name << "\n"
              << " timestep_reprise_vitesse " << (int)1 << "\n"
              << " timestep_reprise_rho " << (int)1 << "\n";
      if (statistiques_.t_integration() > 0.)
        fichier << " statistiques " << statistiques_;
      fichier << "}\n";
    }
}

void DNS_QC_double::sauvegarder_qc(const char *fichier_sauvegarde)
{
  Nom lata_name(fichier_sauvegarde);
  lata_name += ".lata";

  ecrire_fichier_sauv(fichier_sauvegarde, lata_name);
  if (sauvegarde_splitting_name_ != "??")
    {
      for (int i = 0; i < 3; i++)
        {
          redistribute_to_sauvegarde_splitting_faces_[i].redistribute(velocity_[i], velocity_sauvegarde_[i]);
          velocity_sauvegarde_[i].echange_espace_virtuel(velocity_sauvegarde_[i].ghost());
        }
      redistribute_to_sauvegarde_splitting_elem_.redistribute(rho_, rho_sauvegarde_);
      rho_sauvegarde_.echange_espace_virtuel(rho_sauvegarde_.ghost());

      int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
      int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
      int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

      if (ni > 0 || nj > 0 || nk > 0)
        {
          dumplata_header(lata_name, rho_sauvegarde_ /* on passe un champ pour ecrire la geometrie */);
          dumplata_newtime(lata_name, current_time_);
          dumplata_vector_parallele_plan(lata_name, "VELOCITY", velocity_sauvegarde_[0], velocity_sauvegarde_[1], velocity_sauvegarde_[2], 0);
          dumplata_scalar_parallele_plan(lata_name, "RHO", rho_sauvegarde_, 0);
        }
    }
  else
    {
      dumplata_header(lata_name, rho_ /* on passe un champ pour ecrire la geometrie */);
      dumplata_newtime(lata_name,current_time_);
      dumplata_vector(lata_name,"VELOCITY", velocity_[0], velocity_[1], velocity_[2], 0);
      dumplata_scalar(lata_name,"RHO", rho_, 0);
    }

  if (dt_post_spectral_ > 0 )
    partie_fourier_.sauvegarde();

}

void DNS_QC_double::reprendre_qc(const char *fichier_reprise)
{
  Cerr << "Reprise du calcul dans le fichier " << fichier_reprise << finl;
  // Lecture par tous les processeurs, on retire les commentaires etc...
  LecFicDiffuse_JDD fichier(fichier_reprise);
  Param param(que_suis_je());
  param.ajouter("tinit", &current_time_);
  param.ajouter("terme_acceleration_init", &terme_source_acceleration_, Param::REQUIRED);
  param.ajouter("p_thermo_init", &P_thermodynamique_, Param::REQUIRED);
  param.ajouter("fichier_reprise_vitesse", &fichier_reprise_vitesse_);
  param.ajouter("fichier_reprise_rho", &fichier_reprise_rho_);
  param.ajouter("timestep_reprise_vitesse", &timestep_reprise_vitesse_);
  param.ajouter("timestep_reprise_rho", &timestep_reprise_rho_);
  Cerr << "Avant ajout statistiques" << finl;
  param.ajouter("statistiques", &statistiques_);
  Cerr << "Apres ajout statistiques" << finl;
  param.lire_avec_accolades(fichier);
  // Appeler ensuite initialize() pour lire les fichiers lata etc...
  Cerr << "Reprise des donnees a t=" << current_time_ << "\n P_thermo=" << P_thermodynamique_ << finl;
}

// MD 31/07/2019
// Calcul de la vitesse aux elements
void DNS_QC_double::calculer_velocity_elem()
{
  const int ni = velocity_elem_X_.ni();
  const int nj = velocity_elem_Y_.nj();
  const int nk = velocity_elem_Z_.nk();

  for (int k = 0 ; k < nk ; k++)
    {
      for (int j = 0 ; j < nj ; j++)
        {
          for (int i = 0 ; i < ni ; i++)
            {
              velocity_elem_X_(i,j,k) = 0.5*(velocity_[0](i,j,k)+velocity_[0](i+1,j,k));
              velocity_elem_Y_(i,j,k) = 0.5*(velocity_[1](i,j,k)+velocity_[1](i,j+1,k));
              //if (kg == (nktot-1)) {
              //  velocity_elem_Z_(i,j,k) = 0.5*(velocity_[2](i,j,k) + 0.);
              //} else {
              velocity_elem_Z_(i,j,k) = 0.5*(velocity_[2](i,j,k) + velocity_[2](i,j,k+1));
              //}
            }
        }
    }
  velocity_elem_X_.echange_espace_virtuel(1);
  velocity_elem_Y_.echange_espace_virtuel(1);
  velocity_elem_Z_.echange_espace_virtuel(1);
}

void DNS_QC_double::posttraiter_champs_instantanes(const char *lata_name, double current_time)
{
  // DD,2016-28-01: statistiques a partir de fichiers lata uniquement
  if (sauvegarde_post_instantanes_ )
    {
      Nom nom_fichier_sauvegarde(lata_name);
      nom_fichier_sauvegarde += Nom(compteur_post_instantanes_);
      nom_fichier_sauvegarde += Nom(".sauv");

      Nom lata_name_lata(nom_fichier_sauvegarde);
      lata_name_lata += ".lata";

      ecrire_fichier_sauv(nom_fichier_sauvegarde, lata_name_lata);

      dumplata_header(lata_name_lata, rho_ /* on passe un champ pour ecrire la geometrie */);
      dumplata_newtime(lata_name_lata, current_time_);

      if (liste_post_instantanes_.contient_("TOUS"))
        {
          liste_post_instantanes_.dimensionner_force(0);
          liste_post_instantanes_.add("VELOCITY");
          liste_post_instantanes_.add("VELOCITY_ELEM_X");
          liste_post_instantanes_.add("VELOCITY_ELEM_Y");
          liste_post_instantanes_.add("VELOCITY_ELEM_Z");
          liste_post_instantanes_.add("D_VELOCITY");
          liste_post_instantanes_.add("PRESSURE");
          liste_post_instantanes_.add("TEMPERATURE");
          liste_post_instantanes_.add("RHO");
          liste_post_instantanes_.add("LAMBDA");
          liste_post_instantanes_.add("MU");
          liste_post_instantanes_.add("PRESSURE_RHS");
          liste_post_instantanes_.add("DIV_LAMBDA_GRAD_T_VOLUME");
          liste_post_instantanes_.add("U_DIV_RHO_U");
          liste_post_instantanes_.add("DRHO_DT");
          if (turbulent_viscosity_ && (!flag_nu_tensorial_)) liste_post_instantanes_.add("TURBULENT_MU");
          if (turbulent_viscosity_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_XX");
          if (turbulent_viscosity_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_XY");
          if (turbulent_viscosity_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_XZ");
          if (turbulent_viscosity_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_YY");
          if (turbulent_viscosity_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_YZ");
          if (turbulent_viscosity_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_ZZ");
          if (turbulent_diffusivity_ && (!flag_kappa_vectorial_)) liste_post_instantanes_.add("TURBULENT_KAPPA");
          if (turbulent_diffusivity_ && flag_kappa_vectorial_) liste_post_instantanes_.add("TURBULENT_KAPPA_X");
          if (turbulent_diffusivity_ && flag_kappa_vectorial_) liste_post_instantanes_.add("TURBULENT_KAPPA_Y");
          if (turbulent_diffusivity_ && flag_kappa_vectorial_) liste_post_instantanes_.add("TURBULENT_KAPPA_Z");
          if (structural_uu_) liste_post_instantanes_.add("sTRUCTURAL_UU_XX");
          if (structural_uu_) liste_post_instantanes_.add("sTRUCTURAL_UU_XY");
          if (structural_uu_) liste_post_instantanes_.add("sTRUCTURAL_UU_XZ");
          if (structural_uu_) liste_post_instantanes_.add("sTRUCTURAL_UU_YY");
          if (structural_uu_) liste_post_instantanes_.add("sTRUCTURAL_UU_YZ");
          if (structural_uu_) liste_post_instantanes_.add("sTRUCTURAL_UU_ZZ");
          if (structural_uscalar_) liste_post_instantanes_.add("sTRUCTURAL_USCALAR_X");
          if (structural_uscalar_) liste_post_instantanes_.add("sTRUCTURAL_USCALAR_Y");
          if (structural_uscalar_) liste_post_instantanes_.add("sTRUCTURAL_USCALAR_Z");
          if (flag_u_filtre_) liste_post_instantanes_.add("VELOCITY_FILTRE");
          if (flag_rho_filtre_) liste_post_instantanes_.add("RHO_FILTRE");
          if (flag_temperature_filtre_) liste_post_instantanes_.add("TEMPERATURE_FILTRE");
          if (flag_turbulent_mu_filtre_ && (!flag_nu_tensorial_)) liste_post_instantanes_.add("TURBULENT_MU_FILTRE");
          if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_XX");
          if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_FILTRE_XY");
          if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_FILTRE_XZ");
          if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_FILTRE_YY");
          if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_FILTRE_YZ");
          if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_) liste_post_instantanes_.add("TURBULENT_MU_FILTRE_ZZ");
          if (flag_turbulent_kappa_filtre_ && (!flag_kappa_vectorial_)) liste_post_instantanes_.add("TURBULENT_KAPPA_FILTRE");
          if (flag_turbulent_kappa_filtre_ && flag_kappa_vectorial_) liste_post_instantanes_.add("TURBULENT_KAPPA_FILTRE_X");
          if (flag_turbulent_kappa_filtre_ && flag_kappa_vectorial_) liste_post_instantanes_.add("TURBULENT_KAPPA_FILTRE_Y");
          if (flag_turbulent_kappa_filtre_ && flag_kappa_vectorial_) liste_post_instantanes_.add("TURBULENT_KAPPA_FILTRE_Z");
          if (flag_structural_uu_filtre_) liste_post_instantanes_.add("sTRUCTURAL_UU_FILTRE_XX");
          if (flag_structural_uu_filtre_) liste_post_instantanes_.add("sTRUCTURAL_UU_FILTRE_XY");
          if (flag_structural_uu_filtre_) liste_post_instantanes_.add("sTRUCTURAL_UU_FILTRE_XZ");
          if (flag_structural_uu_filtre_) liste_post_instantanes_.add("sTRUCTURAL_UU_FILTRE_YY");
          if (flag_structural_uu_filtre_) liste_post_instantanes_.add("sTRUCTURAL_UU_FILTRE_YZ");
          if (flag_structural_uu_filtre_) liste_post_instantanes_.add("sTRUCTURAL_UU_FILTRE_ZZ");
          if (flag_structural_uscalar_filtre_) liste_post_instantanes_.add("sTRUCTURAL_USCALAR_FILTRE_X");
          if (flag_structural_uscalar_filtre_) liste_post_instantanes_.add("sTRUCTURAL_USCALAR_FILTRE_Y");
          if (flag_structural_uscalar_filtre_) liste_post_instantanes_.add("sTRUCTURAL_USCALAR_FILTRE_Z");
        }
      int n = liste_post_instantanes_.size();
      Cerr << " Martinn n= " << n << finl;

      if (liste_post_instantanes_.contient_("VELOCITY"))
        {
          n=n-1;
          dumplata_vector(lata_name_lata,"VELOCITY", velocity_[0], velocity_[1], velocity_[2], 0);
        }
      if (liste_post_instantanes_.contient_("VELOCITY_ELEM_X"))
        {
          n=n-1;
          calculer_velocity_elem();
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(velocity_elem_X_, velocity_elem_X_sauvegarde_);
              velocity_elem_X_sauvegarde_.echange_espace_virtuel(velocity_elem_X_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "VELOCITY_ELEM_X", velocity_elem_X_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"VELOCITY_ELEM_X", velocity_elem_X_, 0);
            }
        }

      if (liste_post_instantanes_.contient_("VELOCITY_ELEM_Y"))
        {
          n=n-1;
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(velocity_elem_Y_, velocity_elem_Y_sauvegarde_);
              velocity_elem_Y_sauvegarde_.echange_espace_virtuel(velocity_elem_Y_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "VELOCITY_ELEM_Y", velocity_elem_Y_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"VELOCITY_ELEM_Y", velocity_elem_Y_, 0);
            }
        }

      if (liste_post_instantanes_.contient_("VELOCITY_ELEM_Z"))
        {
          n=n-1;
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(velocity_elem_Z_, velocity_elem_Z_sauvegarde_);
              velocity_elem_Z_sauvegarde_.echange_espace_virtuel(velocity_elem_Z_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "VELOCITY_ELEM_Z", velocity_elem_Z_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"VELOCITY_ELEM_Z", velocity_elem_Z_, 0);
            }
        }
      if (liste_post_instantanes_.contient_("D_VELOCITY"))
        n--, dumplata_vector(lata_name_lata,"D_VELOCITY", d_velocity_[0], d_velocity_[1], d_velocity_[2], 0);
      if (liste_post_instantanes_.contient_("PRESSURE"))
        {
          n=n-1;
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(pressure_, pressure_sauvegarde_);
              pressure_sauvegarde_.echange_espace_virtuel(pressure_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "PRESSURE", pressure_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"PRESSURE", pressure_, 0);
            }
        }
      if (liste_post_instantanes_.contient_("TEMPERATURE"))
        {
          n=n-1;
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(temperature_, temperature_sauvegarde_);
              temperature_sauvegarde_.echange_espace_virtuel(temperature_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "TEMPERATURE", temperature_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"TEMPERATURE", temperature_, 0);
            }
        }
      if (liste_post_instantanes_.contient_("RHO"))
        {
          n=n-1;
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(rho_, rho_sauvegarde_);
              rho_sauvegarde_.echange_espace_virtuel(rho_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "RHO", rho_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"RHO", rho_, 0);
            }
        }
      if (liste_post_instantanes_.contient_("LAMBDA"))
        {
          n=n-1;
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(molecular_lambda_, molecular_lambda_sauvegarde_);
              molecular_lambda_sauvegarde_.echange_espace_virtuel(molecular_lambda_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "LAMBDA", molecular_lambda_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"LAMBDA", molecular_lambda_, 0);
            }
        }
      if (liste_post_instantanes_.contient_("MU"))
        {
          n=n-1;
          if (sauvegarde_splitting_name_ != "??")
            {
              redistribute_to_sauvegarde_splitting_elem_.redistribute(molecular_mu_, molecular_mu_sauvegarde_);
              molecular_mu_sauvegarde_.echange_espace_virtuel(molecular_mu_sauvegarde_.ghost());

              int ni = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_I);
              int nj = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_J);
              int nk = sauvegarde_splitting_.get_nb_elem_local(DIRECTION_K);

              if (ni > 0 || nj > 0 || nk > 0)
                {
                  dumplata_scalar_parallele_plan(lata_name_lata, "MU", molecular_mu_sauvegarde_, 0);
                }
            }
          else
            {
              dumplata_scalar(lata_name_lata,"MU", molecular_mu_, 0);
            }
        }
      if (liste_post_instantanes_.contient_("PRESSURE_RHS"))
        n--,dumplata_scalar(lata_name_lata,"PRESSURE_RHS", pressure_rhs_, 0);
      if (liste_post_instantanes_.contient_("DIV_LAMBDA_GRAD_T_VOLUME"))
        n--,dumplata_scalar(lata_name_lata,"DIV_LAMBDA_GRAD_T_VOLUME", div_lambda_grad_T_volume_, 0);
      if (liste_post_instantanes_.contient_("U_DIV_RHO_U"))
        n--,dumplata_scalar(lata_name_lata,"U_DIV_RHO_U", u_div_rho_u_, 0);
      if (liste_post_instantanes_.contient_("DRHO_DT"))
        n--,dumplata_scalar(lata_name_lata,"DRHO_DT", d_rho_, 0);
      if (turbulent_viscosity_ && (!flag_nu_tensorial_) && liste_post_instantanes_.contient_("TURBULENT_MU"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU", turbulent_mu_, 0);
      if (turbulent_viscosity_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_XX"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_XX", turbulent_mu_tensor_[0], 0);
      if (turbulent_viscosity_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_XY"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_XY", turbulent_mu_tensor_[1], 0);
      if (turbulent_viscosity_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_XZ"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_XZ", turbulent_mu_tensor_[2], 0);
      if (turbulent_viscosity_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_YY"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_YY", turbulent_mu_tensor_[3], 0);
      if (turbulent_viscosity_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_YZ"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_YZ", turbulent_mu_tensor_[4], 0);
      if (turbulent_viscosity_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_ZZ"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_ZZ", turbulent_mu_tensor_[5], 0);
      if (turbulent_diffusivity_ && (!flag_kappa_vectorial_) && liste_post_instantanes_.contient_("TURBULENT_KAPPA"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA", turbulent_kappa_, 0);
      if (turbulent_diffusivity_ && flag_kappa_vectorial_ && liste_post_instantanes_.contient_("TURBULENT_KAPPA_X"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA_X", turbulent_kappa_vector_[0], 0);
      if (turbulent_diffusivity_ && flag_kappa_vectorial_ && liste_post_instantanes_.contient_("TURBULENT_KAPPA_Y"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA_Y", turbulent_kappa_vector_[1], 0);
      if (turbulent_diffusivity_ && flag_kappa_vectorial_ && liste_post_instantanes_.contient_("TURBULENT_KAPPA_Z"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA_Z", turbulent_kappa_vector_[2], 0);
      if (structural_uu_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_XX"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_XX", structural_uu_tensor_[0], 0);
      if (structural_uu_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_XY"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_XY", structural_uu_tensor_[1], 0);
      if (structural_uu_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_XZ"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_XZ", structural_uu_tensor_[2], 0);
      if (structural_uu_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_YY"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_YY", structural_uu_tensor_[3], 0);
      if (structural_uu_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_YZ"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_YZ", structural_uu_tensor_[4], 0);
      if (structural_uu_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_ZZ"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_ZZ", structural_uu_tensor_[5], 0);
      if (structural_uscalar_ && liste_post_instantanes_.contient_("sTRUCTURAL_USCALAR_X"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_USCALAR_X", structural_uscalar_vector_[0], 0);
      if (structural_uscalar_ && liste_post_instantanes_.contient_("sTRUCTURAL_USCALAR_Y"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_USCALAR_Y", structural_uscalar_vector_[1], 0);
      if (structural_uscalar_ && liste_post_instantanes_.contient_("sTRUCTURAL_USCALAR_Z"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_USCALAR_Z", structural_uscalar_vector_[2], 0);
      if (flag_u_filtre_ && liste_post_instantanes_.contient_("VELOCITY_FILTRE"))
        n--,dumplata_vector(lata_name_lata,"VELOCITY_FILTRE", velocity_filtre_[0], velocity_filtre_[1], velocity_filtre_[2], 0);
      if (flag_rho_filtre_ && liste_post_instantanes_.contient_("RHO_FILTRE"))
        n--,dumplata_scalar(lata_name_lata,"RHO_FILTRE", rho_filtre_, 0);
      if (flag_temperature_filtre_ && liste_post_instantanes_.contient_("TEMPERATURE_FILTRE"))
        n--,dumplata_scalar(lata_name_lata,"TEMPERATURE_FILTRE", temperature_filtre_, 0);
      if (flag_turbulent_mu_filtre_ && (!flag_nu_tensorial_) && liste_post_instantanes_.contient_("TURBULENT_MU_FILTRE"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_FILTRE", turbulent_mu_filtre_, 0);
      if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_FILTRE_XX"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_FILTRE_XX", turbulent_mu_filtre_tensor_[0], 0);
      if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_FILTRE_XY"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_FILTRE_XY", turbulent_mu_filtre_tensor_[1], 0);
      if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_FILTRE_XZ"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_FILTRE_XZ", turbulent_mu_filtre_tensor_[2], 0);
      if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_FILTRE_YY"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_FILTRE_YY", turbulent_mu_filtre_tensor_[3], 0);
      if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_FILTRE_YZ"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_FILTRE_YZ", turbulent_mu_filtre_tensor_[4], 0);
      if (flag_turbulent_mu_filtre_ && flag_nu_tensorial_ && liste_post_instantanes_.contient_("TURBULENT_MU_FILTRE_ZZ"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_MU_FILTRE_ZZ", turbulent_mu_filtre_tensor_[5], 0);
      if (flag_turbulent_kappa_filtre_ && (!flag_kappa_vectorial_) && liste_post_instantanes_.contient_("TURBULENT_KAPPA_FILTRE"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA_FILTRE", turbulent_kappa_filtre_, 0);
      if (flag_turbulent_kappa_filtre_ && flag_kappa_vectorial_ && liste_post_instantanes_.contient_("TURBULENT_KAPPA_FILTRE_X"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA_FILTRE_X", turbulent_kappa_filtre_vector_[0], 0);
      if (flag_turbulent_kappa_filtre_ && flag_kappa_vectorial_ && liste_post_instantanes_.contient_("TURBULENT_KAPPA_FILTRE_Y"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA_FILTRE_Y", turbulent_kappa_filtre_vector_[1], 0);
      if (flag_turbulent_kappa_filtre_ && flag_kappa_vectorial_ && liste_post_instantanes_.contient_("TURBULENT_KAPPA_FILTRE_Z"))
        n--,dumplata_scalar(lata_name_lata,"TURBULENT_KAPPA_FILTRE_Z", turbulent_kappa_filtre_vector_[2], 0);
      if (flag_structural_uu_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_FILTRE_XX"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_FILTRE_XX", structural_uu_filtre_tensor_[0], 0);
      if (flag_structural_uu_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_FILTRE_XY"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_FILTRE_XY", structural_uu_filtre_tensor_[1], 0);
      if (flag_structural_uu_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_FILTRE_XZ"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_FILTRE_XZ", structural_uu_filtre_tensor_[2], 0);
      if (flag_structural_uu_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_FILTRE_YY"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_FILTRE_YY", structural_uu_filtre_tensor_[3], 0);
      if (flag_structural_uu_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_FILTRE_YZ"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_FILTRE_YZ", structural_uu_filtre_tensor_[4], 0);
      if (flag_structural_uu_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_UU_FILTRE_ZZ"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_UU_FILTRE_ZZ", structural_uu_filtre_tensor_[5], 0);
      if (flag_structural_uscalar_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_USCALAR_FILTRE_X"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_USCALAR_FILTRE_X", structural_uscalar_filtre_vector_[0], 0);
      if (flag_structural_uscalar_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_USCALAR_FILTRE_Y"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_USCALAR_FILTRE_Y", structural_uscalar_filtre_vector_[1], 0);
      if (flag_structural_uscalar_filtre_ && liste_post_instantanes_.contient_("sTRUCTURAL_USCALAR_FILTRE_Z"))
        n--,dumplata_scalar(lata_name_lata,"sTRUCTURAL_USCALAR_FILTRE_Z", structural_uscalar_filtre_vector_[2], 0);
      if (n>0)
        {
          Cerr << " Martinn n2= " << n << finl;
          Cerr << "Il y a des noms de champs a postraiter inconnus ou dupliques dans la liste de champs a postraiter"
               << finl << liste_post_instantanes_ << finl;
          Process::exit();
        }
      compteur_post_instantanes_++;
    }

  if (Process::je_suis_maitre())
    {
      Nom nom_fichier("moyenne_spatiale_");
      nom_fichier += Nom(current_time);
      nom_fichier += Nom(".txt");
      SFichier f(nom_fichier);
      // F.A modification de la precision pour allez chercher les 4 ordres
      f.setf(ios::scientific);
      f.precision(15);
      statistiques_.postraiter(f, 1 /* flag pour ecrire la moyenne instantanee */);
      // modif AT 20/06/2013
      if (statistiques_.check_converge())
        {
          Nom nom_fichier_ec("spatiale_ec_");
          nom_fichier_ec += Nom(current_time);
          nom_fichier_ec += Nom(".txt");
          SFichier fk(nom_fichier_ec);
          // F.A modification de la precision pour allez chercher les 4 ordres
          fk.setf(ios::scientific);
          fk.precision(15);
          statistiques_.postraiter_k(fk,1 /* valeur spatiale */);
        }
    }

  if (Process::je_suis_maitre() && statistiques_.t_integration() > 0.)
    {
      Nom nom_fichier("statistiques_");
      nom_fichier += Nom(current_time);
      nom_fichier += Nom(".txt");
      SFichier fs(nom_fichier);
      // F.A modification de la precision pour allez chercher les 4 ordres
      fs.setf(ios::scientific);
      fs.precision(15);
      statistiques_.postraiter(fs,0 /* moyenne temporelle */);
    }

  if (Process::je_suis_maitre() && statistiques_.t_integration_k() > 0. )
    {
      Nom nom_fichier("stat_ec_");
      nom_fichier += Nom(current_time);
      nom_fichier += Nom(".txt");
      SFichier fk(nom_fichier);
      // F.A modification de la precision pour allez chercher les 4 ordres
      fk.setf(ios::scientific);
      fk.precision(15);
      statistiques_.postraiter_k(fk,0 /* moyenne temporelle */);
    }
}

void DNS_QC_double::run()
{
  Cerr << "IJK_problem_double::run()" << finl;

  if (dt_start_ == 1.e20)
    Cerr << " pas de dt_start choisi, si initialisation d'un calcul il faut faire un dt_start petit ou une projection initiale si isotherme " << finl;

  if(conv_qdm_negligeable_)
    Cerr << " Attention la convection de la vitesse est negligee " << finl;
  if(conv_rho_negligeable_)
    Cerr << " Attention la convection de la masse volumique est negligee " << finl;
  if(diff_qdm_negligeable_)
    Cerr << " Attention la diffusion visqueuse est negligee " << finl;
  if(diff_temp_negligeable_)
    Cerr << " Attention la diffusion thermique est negligee " << finl;

  kernel_ = NULL;
  int ghost_size_filter = 0, ghost_size_d_velocity_tmp;
  if (flag_turbulent_mu_filtre_ || flag_turbulent_kappa_filtre_
      || flag_structural_uu_filtre_ || flag_structural_uscalar_filtre_
      || flag_u_filtre_ || flag_rho_filtre_ || flag_temperature_filtre_
      || flag_structural_uu_tmp_ || flag_structural_uscalar_tmp_
      || flag_d_velocity_tmp_)
    {
      int ghost_size = 0;

      choix_filter_kernel(ghost_size, filter_kernel_name_, kernel_);
      Cerr << "Filter: The ghost size is " << kernel_->ghost_size() << finl;

      ghost_size_filter = 1 + kernel_->ghost_size();
    }

  int ghost_size_velocity = flag_u_filtre_ ? max((int)2, ghost_size_filter) : 2;
  int ghost_size_rho = flag_rho_filtre_ ? max((int)2, ghost_size_filter) : 2;
  int ghost_size_temperature = flag_temperature_filtre_ ? max((int)2, ghost_size_filter) : 2;
  int ghost_size_pressure = (lecture_post_instantanes_filtrer_p_ || lecture_post_instantanes_filtrer_tous_) ? max((int)1, ghost_size_filter) : 1;

  /* allocation des tableaux */
  allocate_velocity(velocity_, splitting_, ghost_size_velocity);
  velocity_elem_X_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  velocity_elem_Y_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  velocity_elem_Z_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  rho_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_rho);
  allocate_velocity(rho_v_, splitting_, 2);
  allocate_velocity(d_velocity_, splitting_, 1);
  allocate_velocity(RK3_F_velocity_, splitting_, 0);
  pressure_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_pressure);
  molecular_mu_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  pressure_rhs_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  d_rho_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  temperature_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_temperature);
  RK3_F_rho_.allocate(splitting_, IJK_Splitting::ELEM, 0);
  molecular_lambda_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  div_lambda_grad_T_volume_.allocate(splitting_, IJK_Splitting::ELEM, 0);
  u_div_rho_u_.allocate(splitting_, IJK_Splitting::ELEM, 1);
  divergence_.allocate(splitting_, IJK_Splitting::ELEM, 1);

  //Modif Martin
  if (sauvegarde_splitting_name_ != "??")
    {
      if (( sauvegarde_splitting_.get_local_slice_index(0) == sauvegarde_splitting_.get_nprocessor_per_direction(0) - 1) || ( sauvegarde_splitting_.get_local_slice_index(1) == sauvegarde_splitting_.get_nprocessor_per_direction(1) - 1) || ( sauvegarde_splitting_.get_local_slice_index(2) == sauvegarde_splitting_.get_nprocessor_per_direction(2) - 1))
        {
          allocate_velocity(velocity_sauvegarde_, sauvegarde_splitting_, 1);
          temperature_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
          molecular_lambda_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
          molecular_mu_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
          rho_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
          pressure_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
          velocity_elem_X_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
          velocity_elem_Y_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
          velocity_elem_Z_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 1);
        }
      else
        {
          allocate_velocity(velocity_sauvegarde_, sauvegarde_splitting_, 0);
          temperature_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
          molecular_lambda_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
          molecular_mu_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
          rho_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
          pressure_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
          velocity_elem_X_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
          velocity_elem_Y_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
          velocity_elem_Z_sauvegarde_.allocate(sauvegarde_splitting_, IJK_Splitting::ELEM, 0);
        }
    }

  /* initialisation des tableaux */
  // Fill with valid floating point data in walls and ghost cells:
  velocity_[0].data() = 0.;
  velocity_[1].data() = 0.;
  velocity_[2].data() = 0.;
  velocity_elem_X_.data()=0.;
  velocity_elem_Y_.data()=0.;
  velocity_elem_Z_.data()=0.;
  rho_.data() = 1.; // not zero=> might want 1./rho and must not crash
  rho_v_[0].data() = 0.;
  rho_v_[1].data() = 0.;
  rho_v_[2].data() = 0.;
  d_velocity_[0].data() = 0.;
  d_velocity_[1].data() = 0.;
  d_velocity_[2].data() = 0.;
  RK3_F_velocity_[0].data() = 0.;
  RK3_F_velocity_[1].data() = 0.;
  RK3_F_velocity_[2].data() = 0.;
  pressure_.data() = 0.;
  molecular_mu_.data() = 0.;
  pressure_rhs_.data() = 0.;
  d_rho_.data() = 0.;
  temperature_.data() = 293.; // like rho: not zero
  RK3_F_rho_.data() = 0.;
  molecular_lambda_.data() = 0.;
  div_lambda_grad_T_volume_.data() = 0.;
  u_div_rho_u_.data() = 0.;
  divergence_.data() = 0.;


  if (flag_turbulent_mu_filtre_ || flag_turbulent_kappa_filtre_
      || flag_structural_uu_filtre_ || flag_structural_uscalar_filtre_
      || flag_u_filtre_ || flag_rho_filtre_ || flag_temperature_filtre_
      || flag_structural_uu_tmp_ || flag_structural_uscalar_tmp_
      || flag_d_velocity_tmp_)
    {
      if (filter_kernel_name_ == Nom("laplacian"))
        {
          Cerr << "Taille du filtre explicite identique a celle du filtre." << finl;
          facteur_delta_filtre_x_ = sqrt(facteur_delta_x_*facteur_delta_x_ + facteur_delta_x_*facteur_delta_x_);
          facteur_delta_filtre_y_ = sqrt(facteur_delta_y_*facteur_delta_y_ + facteur_delta_y_*facteur_delta_y_);
          delta_z_local_pour_delta_filtre_ = delta_z_local_;
          calculer_delta_z_filtre_identique(splitting_, delta_z_local_pour_delta_, delta_z_local_pour_delta_filtre_);
        }
      else
        {
          const double n_mailles = kernel_->n_mailles();
          Cerr << "Taille du filtre explicite fixee a " << n_mailles << " mailles." << finl;

          facteur_delta_filtre_x_ = sqrt(n_mailles*n_mailles + facteur_delta_x_*facteur_delta_x_);
          facteur_delta_filtre_y_ = sqrt(n_mailles*n_mailles + facteur_delta_y_*facteur_delta_y_);
          delta_z_local_pour_delta_filtre_ = delta_z_local_;
          calculer_delta_z_filtre_n_mailles((int)n_mailles, splitting_, delta_z_local_, delta_z_local_pour_delta_, delta_z_local_pour_delta_filtre_);
        }
      constante_modele_ = delta_z_local_;
      for (int i=0 ; i<18 ; i++)
        {
          const int ni = velocity_[2].ni();
          const int nj = velocity_[2].nj();

          // On voudrait
          // tmp_b_[i].allocate(ni, nj, 1, ghost_size_filter);
          // tmp_a_[i].allocate(ni, 1, 1, ghost_size_filter);
          // mais 'allocate' interdit ghost_size < ni, nj ou nk.
          tmp_b_[i].allocate(ni, nj, ghost_size_filter, ghost_size_filter);
          tmp_a_[i].allocate(ni, ghost_size_filter, ghost_size_filter, ghost_size_filter);
        }
      const int nktot = velocity_[2].get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
      for (int i=0; i<8; i++)
        {
          for (int j=0; j<7; j++)
            {
              ml_[i][j].resize_array(nktot+1);
            }
        }
    }

  // DD,2017-04-27: diffusion modifie en vue de l'ajout de modeles
  if (turbulent_viscosity_)
    {
      int ghost_size_turbulent_mu = flag_turbulent_mu_filtre_ ? max((int)2, ghost_size_filter) : 2;

      if (flag_nu_tensorial_)
        {
          turbulent_mu_tensor_[0].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_mu);
          turbulent_mu_tensor_[1].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_mu);
          turbulent_mu_tensor_[2].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_mu);
          turbulent_mu_tensor_[3].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_mu);
          turbulent_mu_tensor_[4].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_mu);
          turbulent_mu_tensor_[5].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_mu);
          turbulent_mu_tensor_[0].data() = 0.;
          turbulent_mu_tensor_[1].data() = 0.;
          turbulent_mu_tensor_[2].data() = 0.;
          turbulent_mu_tensor_[3].data() = 0.;
          turbulent_mu_tensor_[4].data() = 0.;
          turbulent_mu_tensor_[5].data() = 0.;
        }
      else
        {
          turbulent_mu_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_mu);
          turbulent_mu_.data() = 0.;
        }
    }
  if (turbulent_diffusivity_)
    {
      int ghost_size_turbulent_kappa = flag_turbulent_kappa_filtre_ ? max((int)2, ghost_size_filter) : 2;

      if (flag_kappa_vectorial_)
        {
          turbulent_kappa_vector_[0].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_kappa);
          turbulent_kappa_vector_[1].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_kappa);
          turbulent_kappa_vector_[2].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_kappa);
          turbulent_kappa_vector_[0].data() = 0.;
          turbulent_kappa_vector_[1].data() = 0.;
          turbulent_kappa_vector_[2].data() = 0.;
        }
      else
        {
          turbulent_kappa_.allocate(splitting_, IJK_Splitting::ELEM, ghost_size_turbulent_kappa);
          turbulent_kappa_.data() = 0.;
        }
    }

  if (flag_u_filtre_)
    {
      allocate_velocity(velocity_filtre_, splitting_, 2);
      velocity_filtre_[0].data() = 0.;
      velocity_filtre_[1].data() = 0.;
      velocity_filtre_[2].data() = 0.;
    }
  if (structural_uu_)
    {
      int ghost_size_structural_uu = flag_structural_uu_filtre_ ? max((int)2, ghost_size_filter) : 2;

      structural_uu_tensor_[0].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu);
      structural_uu_tensor_[1].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu);
      structural_uu_tensor_[2].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu);
      structural_uu_tensor_[3].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu);
      structural_uu_tensor_[4].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu);
      structural_uu_tensor_[5].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu);
      structural_uu_tensor_[0].data() = 0.;
      structural_uu_tensor_[1].data() = 0.;
      structural_uu_tensor_[2].data() = 0.;
      structural_uu_tensor_[3].data() = 0.;
      structural_uu_tensor_[4].data() = 0.;
      structural_uu_tensor_[5].data() = 0.;
    }

  if (flag_rho_filtre_)
    {
      rho_filtre_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      rho_filtre_.data() = 1.; // not zero=> might want 1./rho and must not crash
    }
  if (flag_temperature_filtre_)
    {
      temperature_filtre_.allocate(splitting_, IJK_Splitting::ELEM, 2);
      temperature_filtre_.data() = 293.; // like rho: not zero
    }
  if (structural_uscalar_)
    {
      int ghost_size_structural_uscalar = flag_structural_uscalar_filtre_ ? max((int)2, ghost_size_filter) : 2;

      structural_uscalar_vector_[0].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uscalar);
      structural_uscalar_vector_[1].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uscalar);
      structural_uscalar_vector_[2].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uscalar);
      structural_uscalar_vector_[0].data() = 0.;
      structural_uscalar_vector_[1].data() = 0.;
      structural_uscalar_vector_[2].data() = 0.;
    }

  if (flag_turbulent_mu_filtre_)
    {
      if (flag_nu_tensorial_)
        {
          turbulent_mu_filtre_tensor_[0].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_mu_filtre_tensor_[1].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_mu_filtre_tensor_[2].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_mu_filtre_tensor_[3].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_mu_filtre_tensor_[4].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_mu_filtre_tensor_[5].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_mu_filtre_tensor_[0].data() = 0.;
          turbulent_mu_filtre_tensor_[1].data() = 0.;
          turbulent_mu_filtre_tensor_[2].data() = 0.;
          turbulent_mu_filtre_tensor_[3].data() = 0.;
          turbulent_mu_filtre_tensor_[4].data() = 0.;
          turbulent_mu_filtre_tensor_[5].data() = 0.;
        }
      else
        {
          turbulent_mu_filtre_.allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_mu_filtre_.data() = 0.;
        }
    }
  if (flag_turbulent_kappa_filtre_)
    {
      if (flag_kappa_vectorial_)
        {
          turbulent_kappa_filtre_vector_[0].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_kappa_filtre_vector_[1].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_kappa_filtre_vector_[2].allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_kappa_filtre_vector_[0].data() = 0.;
          turbulent_kappa_filtre_vector_[1].data() = 0.;
          turbulent_kappa_filtre_vector_[2].data() = 0.;
        }
      else
        {
          turbulent_kappa_filtre_.allocate(splitting_, IJK_Splitting::ELEM, 2);
          turbulent_kappa_filtre_.data() = 0.;
        }
    }

  if (flag_structural_uu_filtre_)
    {
      structural_uu_filtre_tensor_[0].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uu_filtre_tensor_[1].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uu_filtre_tensor_[2].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uu_filtre_tensor_[3].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uu_filtre_tensor_[4].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uu_filtre_tensor_[5].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uu_filtre_tensor_[0].data() = 0.;
      structural_uu_filtre_tensor_[1].data() = 0.;
      structural_uu_filtre_tensor_[2].data() = 0.;
      structural_uu_filtre_tensor_[3].data() = 0.;
      structural_uu_filtre_tensor_[4].data() = 0.;
      structural_uu_filtre_tensor_[5].data() = 0.;
    }
  if (flag_structural_uscalar_filtre_)
    {
      structural_uscalar_filtre_vector_[0].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uscalar_filtre_vector_[1].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uscalar_filtre_vector_[2].allocate(splitting_, IJK_Splitting::ELEM, 2);
      structural_uscalar_filtre_vector_[0].data() = 0.;
      structural_uscalar_filtre_vector_[1].data() = 0.;
      structural_uscalar_filtre_vector_[2].data() = 0.;
    }

  if (flag_structural_uu_tmp_)
    {
      int ghost_size_structural_uu_tmp = flag_structural_uu_tmp_ ? max((int)2, ghost_size_filter) : 2;

      structural_uu_tmp_tensor_[0].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor_[1].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor_[2].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor_[3].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor_[4].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor_[5].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uu_tmp);
      structural_uu_tmp_tensor_[0].data() = 0.;
      structural_uu_tmp_tensor_[1].data() = 0.;
      structural_uu_tmp_tensor_[2].data() = 0.;
      structural_uu_tmp_tensor_[3].data() = 0.;
      structural_uu_tmp_tensor_[4].data() = 0.;
      structural_uu_tmp_tensor_[5].data() = 0.;
    }
  if (flag_structural_uscalar_tmp_)
    {
      int ghost_size_structural_uscalar_tmp = flag_structural_uscalar_tmp_ ? max((int)2, ghost_size_filter) : 2;

      structural_uscalar_tmp_vector_[0].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uscalar_tmp);
      structural_uscalar_tmp_vector_[1].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uscalar_tmp);
      structural_uscalar_tmp_vector_[2].allocate(splitting_, IJK_Splitting::ELEM, ghost_size_structural_uscalar_tmp);
      structural_uscalar_tmp_vector_[0].data() = 0.;
      structural_uscalar_tmp_vector_[1].data() = 0.;
      structural_uscalar_tmp_vector_[2].data() = 0.;
    }
  if (flag_d_velocity_tmp_)
    {
      ghost_size_d_velocity_tmp = flag_d_velocity_tmp_ ? max((int)2, ghost_size_filter) : 2;

      allocate_velocity(d_velocity_tmp_, splitting_, ghost_size_d_velocity_tmp);
      d_velocity_tmp_[0].data() = 0.;
      d_velocity_tmp_[1].data() = 0.;
      d_velocity_tmp_[2].data() = 0.;
    }

  const int nb_allocated_arrays = 21;
  Cerr << " Allocating " << nb_allocated_arrays << " arrays, approx total size= "
       << (double)(unsigned long)molecular_mu_.data().size_array() * sizeof(double) * nb_allocated_arrays * 9.537E-07 << " MB per core" << finl;


  // DD,2017-04-27: diffusion modifie en vue de l'ajout de modeles
  if (type_velocity_diffusion_ == Nom("simple"))
    {
      Cerr << "The diffusion is 'simple': the flux is 'molecular_mu * grad u'" << finl;
      velocity_diffusion_op_simple_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_diffusion_ == Nom("simple_with_transpose"))
    {
      Cerr << "The diffusion is 'simple_with_transpose': the flux is 'molecular_mu * (grad u + grad^T u)'" << finl;
      velocity_diffusion_op_simple_with_transpose_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_diffusion_ == Nom("full"))
    {
      Cerr << "The diffusion is 'full': the flux is 'molecular_mu * (grad u + grad^T u - 2/3 * div u * Id)'" << finl;
      velocity_diffusion_op_full_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_diffusion_ == Nom("none"))
    {
      Cerr << "The diffusion is 'none': no molecular diffusion" << finl;
    }
  else
    {
      Cerr << "Unknown velocity diffusion operator! " << finl;
      Process::exit();
    }

  if (type_velocity_turbulent_diffusion_ == Nom("simple"))
    {
      Cerr << "The velocity turbulent diffusion is 'simple': the flux is 'turbulent_mu * grad u'" << finl;
      velocity_turbulent_diffusion_op_simple_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_turbulent_diffusion_ == Nom("simple_with_transpose"))
    {
      Cerr << "The velocity turbulent diffusion is 'simple_with_transpose': the flux is 'turbulent_mu * (grad u + grad^T u)'" << finl;
      velocity_turbulent_diffusion_op_simple_with_transpose_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_turbulent_diffusion_ == Nom("full"))
    {
      Cerr << "The velocity turbulent diffusion is 'full': the flux is 'turbulent_mu * (grad u + grad^T u - 2/3 * div u * Id)'" << finl;
      velocity_turbulent_diffusion_op_full_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_turbulent_diffusion_ == Nom("simple_anisotropic"))
    {
      Cerr << "The velocity turbulent diffusion is 'simple_anisotropic': the flux is 'turbulent_mu^a * grad^a u' where (grad^a)_i = Delta_i (grad)_i" << finl;
      velocity_turbulent_diffusion_op_simple_anisotropic_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_turbulent_diffusion_ == Nom("simple_with_transpose_anisotropic"))
    {
      Cerr << "The velocity turbulent diffusion is 'simple_with_transpose_anisotropic': the flux is 'turbulent_mu^a * (grad^a u + grad^a^T u)' where (grad^a)_i = Delta_i (grad)_i" << finl;
      velocity_turbulent_diffusion_op_simple_with_transpose_anisotropic_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_turbulent_diffusion_ == Nom("full_anisotropic"))
    {
      Cerr << "The velocity turbulent diffusion is 'full_anisotropic': the flux is 'turbulent_mu^a * (grad^a u + grad^a^T u - 2/3 * div^a u * Id)' where (grad^a)_i = Delta_i (grad)_i" << finl;
      velocity_turbulent_diffusion_op_full_anisotropic_.initialize(splitting_, boundary_conditions_);
    }
  else if (type_velocity_turbulent_diffusion_ == Nom("none"))
    {
      Cerr << "The velocity turbulent diffusion is 'none': no turbulent diffusion" << finl;
    }
  else
    {
      Cerr << "Unknown velocity turbulent diffusion operator! " << finl;
      Process::exit();
    }

  if (structural_uu_)
    {
      Cerr << "A structural model will be added to the velocity turbulent diffusion. The structural uu tensor coefficients are:";
      Cerr << " xx: " << structural_uu_tensor_coefficients_[0];
      Cerr << " xy: " << structural_uu_tensor_coefficients_[1];
      Cerr << " xz: " << structural_uu_tensor_coefficients_[2];
      Cerr << " yy: " << structural_uu_tensor_coefficients_[3];
      Cerr << " yz: " << structural_uu_tensor_coefficients_[4];
      Cerr << " zz: " << structural_uu_tensor_coefficients_[5];
      Cerr << finl;
      velocity_turbulent_diffusion_op_structural_.initialize(splitting_, boundary_conditions_);
    }

  if (structural_uscalar_)
    {
      Cerr << "A structural model will be added to the scalar turbulent diffusion. The structural uscalar vector coefficients are:"
           << " x: " << structural_uscalar_vector_coefficients_[0]
           << " y: " << structural_uscalar_vector_coefficients_[1]
           << " z: " << structural_uscalar_vector_coefficients_[2]
           << finl;
      operateur_diffusion_temperature_structural_.initialize(splitting_);
    }

  if (type_scalar_turbulent_diffusion_ == Nom("normal"))
    {
      Cerr << "The scalar turbulent diffusion is 'normal': the flux is 'lambda * grad s'" << finl;
      operateur_diffusion_turbulent_scalar_.initialize(splitting_);
    }
  else if (type_scalar_turbulent_diffusion_ == Nom("anisotropic"))
    {
      Cerr << "The scalar turbulent diffusion is 'anisotropic': the flux is 'lambda^a * grad^a s' where (grad^a)_i = Delta_i (grad)_i" << finl;
      operateur_diffusion_turbulent_scalar_anisotropic_.initialize(splitting_);
    }
  else if (type_scalar_turbulent_diffusion_ == Nom("none"))
    {
      Cerr << "The scalar turbulent diffusion is 'none': no scalar turbulent diffusion" << finl;
    }
  else
    {
      Cerr << "Unknown scalar turbulent diffusion operator! " << finl;
      Process::exit();
    }

  // On cree des variables ref_turbulent_mu_ij qui sont des references a
  // turbulent_mu_ij si la viscosite est tensorielle et
  // turbulent_mu_ si la viscosite n'est pas tensorielle.
  IJK_Field_double* ptr_turbulent_mu_xx;
  IJK_Field_double* ptr_turbulent_mu_xy;
  IJK_Field_double* ptr_turbulent_mu_xz;
  IJK_Field_double* ptr_turbulent_mu_yy;
  IJK_Field_double* ptr_turbulent_mu_yz;
  IJK_Field_double* ptr_turbulent_mu_zz;

  if (flag_nu_tensorial_)
    {
      ptr_turbulent_mu_xx = &turbulent_mu_tensor_[0];
      ptr_turbulent_mu_xy = &turbulent_mu_tensor_[1];
      ptr_turbulent_mu_xz = &turbulent_mu_tensor_[2];
      ptr_turbulent_mu_yy = &turbulent_mu_tensor_[3];
      ptr_turbulent_mu_yz = &turbulent_mu_tensor_[4];
      ptr_turbulent_mu_zz = &turbulent_mu_tensor_[5];
      Cerr << "The turbulent viscosity is tensorial. The turbulent viscosity tensor coefficients are:"
           << " xx: " << turbulent_viscosity_tensor_coefficients_[0]
           << " xy: " << turbulent_viscosity_tensor_coefficients_[1];
      Cerr << " xz: " << turbulent_viscosity_tensor_coefficients_[2]
           << " yy: " << turbulent_viscosity_tensor_coefficients_[3]
           << " yz: " << turbulent_viscosity_tensor_coefficients_[4];
      Cerr << " zz: " << turbulent_viscosity_tensor_coefficients_[5]
           << finl;
    }
  else
    {
      ptr_turbulent_mu_xx = &turbulent_mu_;
      ptr_turbulent_mu_xy = &turbulent_mu_;
      ptr_turbulent_mu_xz = &turbulent_mu_;
      ptr_turbulent_mu_yy = &turbulent_mu_;
      ptr_turbulent_mu_yz = &turbulent_mu_;
      ptr_turbulent_mu_zz = &turbulent_mu_;
      Cerr << "The turbulent viscosity is not tensorial." << finl;
    }

  IJK_Field_double& ref_turbulent_mu_xx = *ptr_turbulent_mu_xx;
  IJK_Field_double& ref_turbulent_mu_xy = *ptr_turbulent_mu_xy;
  IJK_Field_double& ref_turbulent_mu_xz = *ptr_turbulent_mu_xz;
  IJK_Field_double& ref_turbulent_mu_yy = *ptr_turbulent_mu_yy;
  IJK_Field_double& ref_turbulent_mu_yz = *ptr_turbulent_mu_yz;
  IJK_Field_double& ref_turbulent_mu_zz = *ptr_turbulent_mu_zz;

  // On cree des variables ref_turbulent_kappa_i qui sont des references a
  // turbulent_kappa_i si la diffusivite est vectorielle et
  // turbulent_kappa_ si la diffusivite n'est pas vectorielle.
  IJK_Field_double* ptr_turbulent_kappa_x;
  IJK_Field_double* ptr_turbulent_kappa_y;
  IJK_Field_double* ptr_turbulent_kappa_z;

  if (flag_kappa_vectorial_)
    {
      ptr_turbulent_kappa_x = &turbulent_kappa_vector_[0];
      ptr_turbulent_kappa_y = &turbulent_kappa_vector_[1];
      ptr_turbulent_kappa_z = &turbulent_kappa_vector_[2];
      Cerr << "The turbulent diffusivity is vectorial. The turbulent diffusivity vector coefficients are:"
           << " x: " << turbulent_diffusivity_vector_coefficients_[0]
           << " y: " << turbulent_diffusivity_vector_coefficients_[1]
           << " z: " << turbulent_diffusivity_vector_coefficients_[2]
           << finl;
    }
  else
    {
      ptr_turbulent_kappa_x = &turbulent_kappa_;
      ptr_turbulent_kappa_y = &turbulent_kappa_;
      ptr_turbulent_kappa_z = &turbulent_kappa_;
      Cerr << "The turbulent diffusivity is not tensorial." << finl;
    }

  IJK_Field_double& ref_turbulent_kappa_x = *ptr_turbulent_kappa_x;
  IJK_Field_double& ref_turbulent_kappa_y = *ptr_turbulent_kappa_y;
  IJK_Field_double& ref_turbulent_kappa_z = *ptr_turbulent_kappa_z;

//  velocity_convection_op_.initialize(splitting_);

  if (convection_velocity_amont_)
    {
      velocity_convection_op_amont_.initialize(splitting_);
    }
  else if (convection_velocity_quicksharp_)
    {
      velocity_convection_op_quicksharp_.initialize(splitting_);
    }
  else if (convection_velocity_centre2_)
    {
      velocity_convection_op_centre_2_.initialize(splitting_, boundary_conditions_);
    }
  else
    {
      // Schema centre4 (utilise par defaut)
      velocity_convection_op_.initialize(splitting_, boundary_conditions_);
    }

  if (convection_rho_centre2_)
    {
      rho_convection_op_centre2_.initialize(splitting_);
    }
  else if (convection_rho_amont_)
    {
      rho_convection_op_amont_.initialize(splitting_);
    }
  else if (convection_rho_centre4_)
    {
      // ATTENTION le schema centre 4 OpCentre4IJK a ete modifie pour devenir un centre 2
      rho_convection_op_centre4_.initialize(splitting_, boundary_conditions_);
    }
  else
    {
      // Schema quick (utilise par defaut)
      rho_convection_op_.initialize(splitting_);
    }
  operateur_diffusion_temperature_.initialize(splitting_);
  if (!disable_solveur_poisson_)
    {
      poisson_solver_.initialize(splitting_);
    }

  initialise();
  force_zero_normal_velocity_on_walls(velocity_[2]);

  static Stat_Counter_Id cnt_dtstab = statistiques().new_counter(1, "calcul dtstab QC");
  static Stat_Counter_Id cnt_updtstat = statistiques().new_counter(1, "update statistiques");
  static Stat_Counter_Id cnt_calctermeacc = statistiques().new_counter(1, "calcul terme acceleration");
  static Stat_Counter_Id cnt_sauvegarde = statistiques().new_counter(1, "checkpointing");
  static Stat_Counter_Id fourier_upstat = statistiques().new_counter(1, "TF update");

  // modif FA AT 16/07/2013 necessaire pour le calcul des derivees
  // dans  statistiques_.update_stat_k(...)
  temperature_.echange_espace_virtuel(1);
  molecular_lambda_.echange_espace_virtuel(1);
  molecular_mu_.echange_espace_virtuel(1);
  pressure_.echange_espace_virtuel(1);
  rho_.echange_espace_virtuel(2); // rho est echange sur deux mailles en prevision de rk_step
  velocity_[0].echange_espace_virtuel(2); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
  velocity_[1].echange_espace_virtuel(2); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
  velocity_[2].echange_espace_virtuel(2); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/

  // Projection initiale sur div(u)=0, si demande: (attention, ne pas le faire en reprise)
  if (!disable_solveur_poisson_)
    {
      if (projection_initiale_demandee_)
        {
          Cerr << "*****************************************************************************\n"
               << "  Attention : projection du champ de vitesse initial sur div(u)=0\n"
               << "*****************************************************************************" << finl;

          // pressure_projection_with_rho(rho_,velocity_[0], velocity_[1], velocity_[2],
          //                              pressure_, 1.0 /* dt */, pressure_rhs_,
          //                              1 /* check divergence */,
          //                              poisson_solver_);
          pressure_.data() = 0.;
          pressure_rhs_.data() = 0.;
          // // Echange espace virtuel inutiles car deja fait dans pressure_projection_with_rho
          // velocity_[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          // velocity_[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          // velocity_[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          ///
        }
    }

  const Nom lata_name = nom_du_cas() + Nom("_lata_");

  ArrOfDouble tmp_size3(3);
  // Postraiter la condition initiale:
  // Calcul des moyennes spatiales sur la condition initiale:
  statistiques_.update_stat(velocity_, pressure_, temperature_, rho_, molecular_mu_, molecular_lambda_,
                            delta_z_local_pour_delta_,
                            flag_nu_anisotropic_, turbulent_viscosity_, ref_turbulent_mu_xx, ref_turbulent_mu_xy, ref_turbulent_mu_xz, ref_turbulent_mu_yy, ref_turbulent_mu_yz, ref_turbulent_mu_zz,
                            flag_kappa_anisotropic_, turbulent_diffusivity_, ref_turbulent_kappa_x, ref_turbulent_kappa_y, ref_turbulent_kappa_z,
                            structural_uu_, structural_uu_tensor_,
                            structural_uscalar_, structural_uscalar_vector_,
                            formulation_favre_, formulation_velocity_,
                            Cp_gaz_, P_thermodynamique_, 0.);
  if (statistiques_.check_converge())
    {
      statistiques_.update_stat_k(velocity_, pressure_,rho_, molecular_mu_, P_thermodynamique_, terme_source_acceleration_, 0.); // modif AT 20/06/2013
    }
  posttraiter_champs_instantanes(lata_name, current_time_);

  // ATTENTION PARTIE pour le test et debugage
  if ( 0 )
    {
      partie_fourier_.update(velocity_, pressure_, rho_, molecular_mu_);
      Cerr << " Post-traitement Spectral " << finl;
      Nom Nom_post_test = "Spectrale";
      partie_fourier_.postraiter(Nom_post_test);
      exit();
    }

  int stop = 0;
  double previous_time = 0;

  // DD,2016-28-01: statistiques a partir de fichiers lata uniquement
  if (lecture_post_instantanes_)
    {
      nb_timesteps_ = statlata_namelist_.size();

      Cerr << "*****************************************************************************\n"
           << "  Attention : On ne fait que postraiter les statistiques depuis les lata\n"
           << "*****************************************************************************" << finl;
    }

  for (int tstep = 0; tstep < nb_timesteps_ && stop == 0; tstep++)
    {
      statistiques().begin_count(timestep_counter_);
      DebogIJK::verifier("Vitesse_X init", velocity_[0]);
      DebogIJK::verifier("Vitesse_Y init", velocity_[1]);
      DebogIJK::verifier("Vitesse_Z init", velocity_[2]);
      DebogIJK::verifier("rho init", rho_);

      // DD,2016-28-01: statistiques a partir de fichiers lata uniquement
      if (lecture_post_instantanes_)
        {
          Nom statlata_fichier = Nom(statlata_namelist_[tstep]);
          statlata_fichier += Nom(".sauv");

          previous_time = current_time_;
          reprendre_qc(statlata_fichier); // current_time_ mis a jour
          timestep_ = current_time_ - previous_time;

          Cout << "T= " << current_time_
               << " timestep= " << timestep_ << finl;

          const Nom& geom_name = splitting_.get_grid_geometry().le_nom();

          Cout << "Lecture rho dans fichier " << fichier_reprise_rho_ << " timestep= " << timestep_reprise_rho_ << finl;
          lire_dans_lata(fichier_reprise_rho_, timestep_reprise_rho_, geom_name, "RHO", rho_);

          Cout << "Lecture pression dans fichier " << fichier_reprise_rho_ << " timestep= " << timestep_reprise_rho_ << finl;
          lire_dans_lata(fichier_reprise_rho_, timestep_reprise_rho_, geom_name, "PRESSURE", pressure_);

          Cout << "Lecture vitesse dans fichier " << fichier_reprise_vitesse_ << " timestep= " << timestep_reprise_vitesse_ << finl;
          lire_dans_lata(fichier_reprise_vitesse_, timestep_reprise_vitesse_, geom_name, "VELOCITY",
                         velocity_[0], velocity_[1], velocity_[2]);

          if (lecture_post_instantanes_filtrer_u_ || lecture_post_instantanes_filtrer_tous_)
            {
              const int flag_add = 0;
              velocity_[0].echange_espace_virtuel(ghost_size_velocity);
              velocity_[1].echange_espace_virtuel(ghost_size_velocity);
              velocity_[2].echange_espace_virtuel(ghost_size_velocity);
              filtrer_champ_elem(flag_add, velocity_[0], delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_, kernel_, tmp_b_, tmp_a_, velocity_[0]);
              filtrer_champ_elem(flag_add, velocity_[1], delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_, kernel_, tmp_b_, tmp_a_, velocity_[1]);
              filtrer_champ_face(flag_add, velocity_[2], delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_, kernel_, tmp_b_, tmp_a_, velocity_[2]);
            }

          if (lecture_post_instantanes_filtrer_rho_ || lecture_post_instantanes_filtrer_tous_)
            {
              const int flag_add = 0;
              rho_.echange_espace_virtuel(ghost_size_rho);
              filtrer_champ_elem(flag_add, rho_, delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_, kernel_, tmp_b_, tmp_a_, rho_);
            }

          if (lecture_post_instantanes_filtrer_p_ || lecture_post_instantanes_filtrer_tous_)
            {
              const int flag_add = 0;
              pressure_.echange_espace_virtuel(ghost_size_pressure);
              filtrer_champ_elem(flag_add, pressure_, delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_, kernel_, tmp_b_, tmp_a_, pressure_);
            }

          rho_.echange_espace_virtuel(1);
          calculer_temperature_mu_lambda_air(P_thermodynamique_, constante_specifique_gaz_, rho_, temperature_, molecular_mu_, molecular_lambda_, 1);

          temperature_.echange_espace_virtuel(1);
          molecular_lambda_.echange_espace_virtuel(1);
          molecular_mu_.echange_espace_virtuel(1);
          pressure_.echange_espace_virtuel(1);
          rho_.echange_espace_virtuel(2);
          velocity_[0].echange_espace_virtuel(2);
          velocity_[1].echange_espace_virtuel(2);
          velocity_[2].echange_espace_virtuel(2);
        }
      else
        {
          // Calcul du pas de temps:

          velocity_[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          velocity_[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          velocity_[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/

          statistiques().begin_count(cnt_dtstab);
          old_timestep_ = timestep_;

          double dt_conv    =  conv_qdm_negligeable_  ? 1.e20 : velocity_convection_op_.compute_dtstab_convection_local(velocity_[0], velocity_[1], velocity_[2]);
          double dt_diff    = 1.e20;
          double dt_diff_mu = 1.e20;
          if (old_dtstab_)
            {
              dt_diff_mu =  diff_qdm_negligeable_  ? 1.e20 : calculer_dtstab_diffusion_temperature_local(molecular_mu_, rho_, 1.);
              dt_diff    =  diff_temp_negligeable_ ? 1.e20 : calculer_dtstab_diffusion_temperature_local(molecular_lambda_, rho_, Cp_gaz_);
            }
          else
            {
              if (large_eddy_simulation_formulation_ == Nom("none")
                  || (formulation_favre_ && (!turbulent_viscosity_))
                  || (formulation_velocity_ && (!turbulent_viscosity_) && (!turbulent_diffusivity_)))
                {
                  dt_diff_mu =  diff_qdm_negligeable_  ? 1.e20 : calculer_dtstab_diffusion_temperature_local(molecular_mu_, rho_, 1.);
                  dt_diff    =  dt_diff_mu;
                }
              else if (formulation_favre_ && turbulent_viscosity_)
                {
                  dt_diff_mu =  diff_qdm_negligeable_  ? 1.e20 : calculer_dtstab_diffusion_temperature_local_avec_turbulent_favre(flag_nu_anisotropic_, molecular_mu_, ref_turbulent_mu_xx, ref_turbulent_mu_xy, ref_turbulent_mu_xz, ref_turbulent_mu_yy, ref_turbulent_mu_yz, ref_turbulent_mu_zz, rho_, 1.);
                  dt_diff    =  dt_diff_mu;
                }
              else if (formulation_velocity_ && turbulent_viscosity_ && turbulent_diffusivity_)
                {
                  dt_diff_mu =  diff_qdm_negligeable_  ? 1.e20 : calculer_dtstab_diffusion_temperature_local_avec_turbulent_velocity(flag_nu_anisotropic_, molecular_mu_, ref_turbulent_mu_xx, ref_turbulent_mu_xy, ref_turbulent_mu_xz, ref_turbulent_mu_yy, ref_turbulent_mu_yz, ref_turbulent_mu_zz, rho_, 1.);
                  double dt_diff_kappa =  calculer_dtstab_diffusion_temperature_local_sans_rho(flag_kappa_anisotropic_, ref_turbulent_kappa_x, ref_turbulent_kappa_y, ref_turbulent_kappa_z, 1.);
                  dt_diff    =  min(dt_diff_mu, dt_diff_kappa);
                }
              else if (formulation_velocity_ && (!turbulent_viscosity_) && turbulent_diffusivity_)
                {
                  dt_diff_mu =  diff_qdm_negligeable_  ? 1.e20 : calculer_dtstab_diffusion_temperature_local(molecular_mu_, rho_, 1.);
                  double dt_diff_kappa =  calculer_dtstab_diffusion_temperature_local_sans_rho(flag_kappa_anisotropic_, ref_turbulent_kappa_x, ref_turbulent_kappa_y, ref_turbulent_kappa_z, 1.);
                  dt_diff    =  min(dt_diff_mu, dt_diff_kappa);
                }
              else if (formulation_velocity_ && turbulent_viscosity_ && (!turbulent_diffusivity_))
                {
                  dt_diff_mu =  diff_qdm_negligeable_  ? 1.e20 : calculer_dtstab_diffusion_temperature_local_avec_turbulent_velocity(flag_nu_anisotropic_, molecular_mu_, ref_turbulent_mu_xx, ref_turbulent_mu_xy, ref_turbulent_mu_xz, ref_turbulent_mu_yy, ref_turbulent_mu_yz, ref_turbulent_mu_zz, rho_, 1.);
                  dt_diff    =  dt_diff_mu;
                }
              else
                {
                  Cerr << "This should not happen." << finl;
                  Process::exit();
                }
            }

          statistiques().end_count(cnt_dtstab);
          tmp_size3[0] = dt_conv;
          tmp_size3[1] = dt_diff;
          tmp_size3[2] = dt_diff_mu;
          mp_min_for_each_item(tmp_size3);
          dt_conv = tmp_size3[0];
          dt_diff = tmp_size3[1];
          dt_diff_mu = tmp_size3[2];

          // ajouter ici methode de non progression du dt si trop faible. (stats plus stables => plus precise).
          const double dt_theorique = 1. / (1./dt_conv + 1./dt_diff);
          double dt_regule = min(dt_theorique * timestep_facsec_, timestep_max_);

          // FA 17/03/14 amelioration pour le pas de temps !
          if (tstep == 0 ) /* au cas ou on a un dt_start */
            timestep_ = min(dt_regule,dt_start_);
          else if ( (dt_regule < timestep_) || ((dt_regule-timestep_)/timestep_ > 0.05 )) // on change le pas de temps si on gagne 5%  ou si il faut descendre
            {
              timestep_ = dt_regule;
            }
          else // sinon on garde
            {
              // Nothing to do
            }

          Cout << "T= " << current_time_
               << " dtconv= " << dt_conv
               << " dtdiff_t= " << dt_diff
               << " dtdiff_v= " << dt_diff_mu;
          Cout << " theorique_dt= " << dt_theorique
               << " dt_limited " << dt_regule
               << " timestep= " << timestep_ << finl;

          // // F.A 3 /04/2014 changement de la source est maintenant calculee au debut de rk_sub_step
          // statistiques().begin_count(cnt_calctermeacc);
          // double v_moy, rho_v_moy;
          // calculer_v_et_rhov_moyen(velocity_[0], rho_, delta_z_local_, volume_total_domaine_, v_moy, rho_v_moy);
          // double derivee_acceleration;
          // if (Process::je_suis_maitre()) {
          //   // Mise a jour de l'acceleration
          //   parser_derivee_acceleration_.setVar("force", terme_source_acceleration_);
          //   parser_derivee_acceleration_.setVar("v_moyen", v_moy);
          //   parser_derivee_acceleration_.setVar("rho_v_moyen", rho_v_moy);
          //   derivee_acceleration = parser_derivee_acceleration_.eval();
          //   terme_source_acceleration_ += derivee_acceleration * timestep_;
          // }
          //

          // F.A 16/04/14 changement de source (encore ...)
          // calcul de la force de recirculation avant modification des tableaux

          statistiques().begin_count(cnt_calctermeacc);
          double debit_old = debit_actuel_;
          calculer_debit(velocity_[0], rho_, delta_z_local_, Lx_tot_, debit_actuel_);
          double old_t_bulk = actual_t_bulk_;
          compute_rho_t_bulk(rho_bulk_, actual_t_bulk_);
          if (Process::je_suis_maitre())
            {
              double acceleration_du_dt = 0.;
              if (mode_terme_source_impose_)
                {
                  terme_source_acceleration_ = terme_source_acceleration_constant_;
                }
              else
                {
                  double ecart_debit =  (debit_cible_ - debit_actuel_ );
                  double ecart_ancien = ( debit_actuel_ - debit_old );
                  acceleration_du_dt = ecart_debit - ecart_ancien;
                  acceleration_du_dt /= dump_factor_ * Ly_tot_ * Lz_tot_ * timestep_;

                  terme_source_acceleration_ += acceleration_du_dt;
                }

              if(flag_t_bulk_forced_)
                {
                  /* double derivative_flux = (sum_entering_flux_ - old_entering_hf)/Lz_tot_; */
                  double derivative_2nd_order_t = aim_t_bulk_ - 2 * actual_t_bulk_ + old_t_bulk;
                  double d = derivative_2nd_order_t / timestep_;
                  /* double fac = dump_factor_2_ *( rho_bulk_ * constante_specifique_gaz_ * d / (gamma_ - 1) - derivative_flux); */
                  double fac = dump_factor_2_ * rho_bulk_ * constante_specifique_gaz_ * d / (gamma_ - 1);
                  puit_ -= fac;
                  Cout << "Current heat sink: " << puit_
                       << " Current T_bulk: " << actual_t_bulk_
                       << " Aim T_bulk: " << aim_t_bulk_
                       << " Old T_bulk: " << old_t_bulk
                       << " factor " << fac
                       << " Timestep " << timestep_
                       << " rho bulk " << rho_bulk_
                       /* << " sum_fluxes: " << sum_entering_flux_ */
                       /* << " old_sum_fluxes: " << old_entering_hf */
                       << finl;
                }


              // la derivee_acceleration n'est connue que sur le maitre
              Cout << "T= " << current_time_
                   << " debit cible= " << debit_cible_
                   << " debit actuel= " << debit_actuel_
                   << " acceleration= " << terme_source_acceleration_
                   << " da/dt= " << acceleration_du_dt << finl;


            }
          envoyer_broadcast(terme_source_acceleration_, 0);
          envoyer_broadcast(puit_, 0);


          // // F.A 3 /04/2014 changement de la source
          // if (Process::je_suis_maitre()) {
          //   // la derivee_acceleration n'est connue que sur le maitre
          //   Cout << "T= " << current_time_
          //        << " Vx_moyen= " << v_moy
          //        << " rho_vx_moyen= " << rho_v_moy
          //        << " acceleration= " << terme_source_acceleration_
          //        << " da/dt= " << derivee_acceleration << finl;
          // }
          //

          statistiques().end_count(cnt_calctermeacc);

          //////////////////////////
          for (int rk_step = 0; rk_step < 3; rk_step++)
            {
              rk3_sub_step(rk_step, timestep_);
              if (postraiter_sous_pas_de_temps_ && (tstep % dt_sauvegarde_ == dt_sauvegarde_-1))
                {
                  statistiques().begin_count(postraitement_counter_);
                  posttraiter_champs_instantanes(lata_name, current_time_ + timestep_ * rk_step / 3);
                  statistiques().end_count(postraitement_counter_);
                }
            }
        }
      // on s'assure que la force est la meme entre chaque pas de temps.
      // if (Process::je_suis_maitre())
      //   envoyer_broadcast(terme_source_acceleration_, 0);

      if (current_time_ >= t_debut_statistiques_)
        {
          statistiques().begin_count(postraitement_counter_);
          statistiques().begin_count(cnt_updtstat); // timer
          statistiques_.update_stat(velocity_, pressure_, temperature_, rho_, molecular_mu_, molecular_lambda_,
                                    delta_z_local_pour_delta_,
                                    flag_nu_anisotropic_, turbulent_viscosity_, ref_turbulent_mu_xx, ref_turbulent_mu_xy, ref_turbulent_mu_xz, ref_turbulent_mu_yy, ref_turbulent_mu_yz, ref_turbulent_mu_zz,
                                    flag_kappa_anisotropic_, turbulent_diffusivity_, ref_turbulent_kappa_x, ref_turbulent_kappa_y, ref_turbulent_kappa_z,
                                    structural_uu_, structural_uu_tensor_,
                                    structural_uscalar_, structural_uscalar_vector_,
                                    formulation_favre_, formulation_velocity_,
                                    Cp_gaz_, P_thermodynamique_, timestep_); // stats standard
          if (statistiques_.check_converge())
            statistiques_.update_stat_k(velocity_, pressure_,rho_, molecular_mu_, P_thermodynamique_, terme_source_acceleration_, timestep_);

          statistiques().end_count(cnt_updtstat); // fin timer
          statistiques().begin_count(fourier_upstat); // timer
          if ( dt_post_spectral_ > 0 )
            {
              int un_sur_20 = tstep%20;
              if ( !un_sur_20 )
                partie_fourier_.update(velocity_, pressure_, rho_, molecular_mu_);
            }
          statistiques().end_count(fourier_upstat); // fin timer
          statistiques().end_count(postraitement_counter_);
        }

      if (lecture_post_instantanes_)
        {
          // do nothing
        }
      else
        {
          current_time_ += timestep_;
        }

      // verification du fichier stop
      stop = 0;
      if (check_stop_file_ != "??")
        {
          if (je_suis_maitre())
            {
              EFichier f;
              stop = f.ouvrir(check_stop_file_);
              if (stop)
                {
                  // file exists, check if it contains 1:
                  f >> stop;
                }
            }
          envoyer_broadcast(stop, 0);
        }
      if (tstep == nb_timesteps_ - 1)
        stop = 1;

      if (tstep % dt_sauvegarde_ == dt_sauvegarde_-1 || stop)
        {
          statistiques().begin_count(cnt_sauvegarde);
          sauvegarder_qc(nom_sauvegarde_);
          statistiques().end_count(cnt_sauvegarde);
        }

      if (tstep % dt_post_ == dt_post_-1 || stop)
        {
          statistiques().begin_count(postraitement_counter_);
          posttraiter_champs_instantanes(lata_name, current_time_);
          statistiques().end_count(postraitement_counter_);
        }
      if ((tstep % dt_post_spectral_ == dt_post_-1 || stop) && (dt_post_spectral_ != -1))
        {
          statistiques().begin_count(postraitement_counter_);
          Nom Nom_post = "Spectrale_";
          Nom_post += Nom(current_time_);
          partie_fourier_.postraiter(Nom_post);
          statistiques().end_count(postraitement_counter_);
        }
      //  calculer_moyennes_flux();

      statistiques().end_count(timestep_counter_);
      if (Process::je_suis_maitre())
        {
          Cerr << "tstep " << tstep
               << " currenttime " << current_time_
               << " cpu_time " << statistiques().last_time(timestep_counter_) << finl;
        }
    }

  delete kernel_;
}


// int calculer_k_pour_bord(const IJK_Field_double& temperature, const bool bord_kmax)
// {
//   const int kmin = temperature.get_splitting().get_offset_local(DIRECTION_K);
//   const int nktot = temperature.get_splitting().get_nb_items_global(IJK_Splitting::ELEM, DIRECTION_K);
//   int k;
//   // calcul l'indice k de la couche de mailles voisine du bord. Si je n'ai pas de bord, on met k = -1
//   if (!bord_kmax)
//     {
//       // on veut le bord "k_global = 0"
//       if (kmin == 0)
//         {
//           // ce bord est chez moi... et il est en k=0
//           k = 0;
//         }
//       else
//         {
//           // ce bord n'est pas chez moi
//           k = -1;
//         }
//     }
//   else
//     {
//       // on veut le bord kmax
//       if (kmin + temperature.nk() == nktot)
//         {
//           // ce bord est chez moi... et il est en k= truc...
//           k = temperature.nk() - 1;
//         }
//       else
//         {
//           k = -1;
//         }
//     }
//   return k;
// }

// valeur de retour: indice local du plan de temperature voisin utilise,
//  -1 si on n'a pas le bord sur ce processeur
// Calcule l'integrale sur chaque face du bord demande du flux de chaleur a travers la face
// positif si le flux va vers les k positifs.
int calculer_flux_thermique_bord(const IJK_Field_double& temperature,
                                 const double lambda_de_t_paroi,
                                 const int turbulent_diffusivity,
                                 const IJK_Field_double& lambda_turbulent,
                                 const int flag_lambda_anisotropic,
                                 const int structural_uscalar,
                                 const IJK_Field_double& structural_uscalar_z,
                                 const double T_paroi_impose,
                                 IJK_Field_local_double& flux_bord,
                                 const bool bord_kmax)
{
  const int kmin = temperature.get_splitting().get_offset_local(DIRECTION_K);
  int k = calculer_k_pour_bord(temperature, bord_kmax);
  if (k == -1)
    return k;

  // redimensionne flux_bord avec ni * nj:
  const int ni = temperature.ni(); // nombre d'element local sur ce processeur
  const int nj = temperature.nj();
  flux_bord.allocate(ni, nj, 1, 0);

  const IJK_Grid_Geometry& geometry = temperature.get_splitting().get_grid_geometry();
  const double delta_k = geometry.get_delta(DIRECTION_K)[k + kmin]; // k+kmin est l'indice global de la maille locale k
  const ArrOfDouble& coord_z = geometry.get_node_coordinates(DIRECTION_K);

  for (int j = 0; j < nj; j++)
    {
      for (int i = 0; i < ni; i++)
        {
          const int sens = bord_kmax ? -1 : 1;

          // si bord bas:  x_p=coord_z[p]; y_k=velocity_k(i,j,p)
          // si bord haut: x_p=coord_z[k+kmin+1-p]; y_k=velocity_k(i,j,k+1-p)
          const double xf_0 = coord_z[k+kmin+bord_kmax];
          const double xf_1 = coord_z[k+kmin+bord_kmax+sens*1];
          const double xf_2 = coord_z[k+kmin+bord_kmax+sens*2];

          const double x_p = xf_0;
          const double x_0 = 0.5*(xf_0+xf_1);
          const double x_1 = 0.5*(xf_1+xf_2);

          double T_p = T_paroi_impose;

          const double T_0 = temperature(i,j,k);
          const double T_1 = temperature(i,j,k+sens*1);
          double derivee_premiere = (x_1-x_p)*(x_0-x_p)/((x_0-x_p)-(x_1-x_p))*(T_1/((x_1-x_p)*(x_1-x_p))-T_0/((x_0-x_p)*(x_0-x_p))-T_p*(((x_0-x_p)*(x_0-x_p))-((x_1-x_p)*(x_1-x_p)))/(((x_0-x_p)*(x_0-x_p))*((x_1-x_p)*(x_1-x_p))));

          double l;
          double flux;
          if (turbulent_diffusivity && (!flag_lambda_anisotropic))
            {
              l = lambda_turbulent(i,j,k) + lambda_de_t_paroi;
            }
          else if (turbulent_diffusivity && flag_lambda_anisotropic)
            {
              l = delta_k*lambda_turbulent(i,j,k) + lambda_de_t_paroi;
            }
          else
            {
              l = lambda_de_t_paroi;
            }
          /*
          // Faut-il ajouter ce if ? Si oui ou ajouter s ?
                if (structural_uscalar) {
                    const double s = structural_uscalar_z(i,j,k);
                    flux = ((T_paroi_impose - t) * l + s) * facteur;
                } else {
                    // le flux est positif s'il va vers les k croissants
                    flux = (T_paroi_impose - t) * l * facteur;
                }
          */
          if (structural_uscalar)
            {
              const double s = structural_uscalar_z(i,j,k);
              flux = -( derivee_premiere * l + s ) * geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
            }
          else
            {
              flux = -( derivee_premiere * l ) * geometry.get_constant_delta(DIRECTION_I) * geometry.get_constant_delta(DIRECTION_J);
            }
          flux_bord(i,j,0) = flux;
        }
    }
  return k;
}


// p_thermo = pression thermodynamique
void DNS_QC_double::calcul_p_thermo_et_bilan(const IJK_Field_double& rho,
                                             IJK_Field_double& temperature,
                                             const int turbulent_diffusivity,
                                             const IJK_Field_double& lambda_turbulent,
                                             const int flag_lambda_anisotropic,
                                             const int structural_uscalar,
                                             const IJK_Field_double& structural_uscalar_z,
                                             const double P_th_initial,
                                             double& P_th_final,
                                             const double fractionnal_timestep,
                                             double& d_Pth_divise_par_gammamoins1) const
{
  IJK_Field_local_double flux_bord;

  const int imax = temperature.ni();
  const int jmax = temperature.nj();

  double P_th = P_th_initial;
  // Boucle point fixe:
  const int max_point_fixe = 3;
  double somme_flux_entrants;

  for (int point_fixe_iter = 0; point_fixe_iter < max_point_fixe; point_fixe_iter++)
    {
      // Calcul t_star a partir de rho et de la valeur courante de pthermo
      // boucle sur les deux plans conditions aux limites:
      double somme_flux_kmin = 0.;
      double somme_flux_kmax = 0.;
      for (int plan_cl = 0; plan_cl < 2; plan_cl++)
        {
          double somme_flux = 0.;
          // indice local du plan de mailles voisin de ce bord:
          int k = calculer_k_pour_bord(temperature, plan_cl);
          if (k > -1)
            {
              for (int j = 0; j < jmax; j++)
                {
                  for (int i = 0; i < imax; i++)
                    {
                      temperature(i,j,k) = P_th / (constante_specifique_gaz_ * rho(i,j,k));
                    }
                }
              const double lambda_de_t_paroi = (plan_cl ? lambda_de_t_paroi_kmax_ : lambda_de_t_paroi_kmin_);
              const double T_paroi_impose = (plan_cl ? T_paroi_impose_kmax_ : T_paroi_impose_kmin_);

              calculer_flux_thermique_bord(temperature, lambda_de_t_paroi,
                                           0, lambda_turbulent, flag_lambda_anisotropic,
                                           0, structural_uscalar_z,
                                           T_paroi_impose, flux_bord, plan_cl);
              // integrale du flux:
              for (int j = 0; j < jmax; j++)
                {
                  for (int i = 0; i < imax; i++)
                    {
                      somme_flux += flux_bord(i,j,0);
                    }
                }
            }
          // le flux est positif si la chaleur va vers les k croissants:
          if (plan_cl)
            somme_flux_kmax = somme_flux;
          else
            somme_flux_kmin = somme_flux;
        }
      // somme sur tous les processeurs:
      somme_flux_kmin = mp_sum(somme_flux_kmin);
      somme_flux_kmax = mp_sum(somme_flux_kmax);
      // Somme des flux entrants dans le domaine:
      if ( diff_temp_negligeable_ && !(formulation_favre_ && turbulent_diffusivity) && !(formulation_favre_ && structural_uscalar) )  // si diffusion negligable
        somme_flux_entrants =0;
      else
        somme_flux_entrants = somme_flux_kmin - somme_flux_kmax;
      P_th = P_th_initial / (1. - (gamma_ - 1) * fractionnal_timestep / P_th * ( somme_flux_entrants / volume_total_domaine_ - puit_ ));

      assert_parallel(P_th);
      if (point_fixe_iter== max_point_fixe-1)
        {
          Cout << "calcul_p_thermo iter " << point_fixe_iter << " flux k=0: " << somme_flux_kmin
               << " flux k=kmax: " << somme_flux_kmax << " P_th: " << P_th << " ";
        }
    }
  P_th_final = P_th;

  d_Pth_divise_par_gammamoins1 = somme_flux_entrants / volume_total_domaine_ - puit_;
  Cout << "dP_th/dt= " << d_Pth_divise_par_gammamoins1 * (gamma_ - 1) << finl;
}

static void force_2d(IJK_Field_double& f)
{
  const int ni = f.ni();
  const int nj = f.nj();
  const int nk = f.nk();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 1; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              f(i, j, k) = f(i, 0, k);
            }
        }
    }
}

template<class T>
void DNS_QC_double::calculer_convection_vitesse(FixedVector<IJK_Field_double, 3>& rho_v,
                                                FixedVector<IJK_Field_double, 3>& velocity,
                                                const ArrOfDouble_with_ghost& delta_z,
                                                const double facteur_delta_x,
                                                const double facteur_delta_y,
                                                const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                                T& kernel,
                                                FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                                FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                                FixedVector<IJK_Field_double, 3>& d_velocity_tmp,
                                                FixedVector<IJK_Field_double, 3>& d_velocity,
                                                IJK_Field_double& u_div_rho_u)
{
  int ghost_size_filter;
  int ghost_size_d_velocity_tmp;
  if (flag_convection_qdm_sans_divergence_)
    {
      if (flag_filtrage_convection_qdm_)
        {
          rho_v[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          rho_v[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          rho_v[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          compute_divergence_times_constant(rho_v[0], rho_v[1], rho_v[2], 1., u_div_rho_u);
          u_div_rho_u.echange_espace_virtuel(1);

          multiplier_champ_rho_face_i(false, u_div_rho_u, 1., velocity[0], d_velocity_tmp[0]);
          multiplier_champ_rho_face_j(false, u_div_rho_u, 1., velocity[1], d_velocity_tmp[1]);
          multiplier_champ_rho_face_k(false, u_div_rho_u, 1., 0., velocity[2], d_velocity_tmp[2]);

          ghost_size_filter = 1 + kernel->ghost_size();
          ghost_size_d_velocity_tmp = max((int)2, ghost_size_filter);
          d_velocity_tmp[0].echange_espace_virtuel(ghost_size_d_velocity_tmp);
          d_velocity_tmp[1].echange_espace_virtuel(ghost_size_d_velocity_tmp);
          d_velocity_tmp[2].echange_espace_virtuel(ghost_size_d_velocity_tmp);

          const int flag_add = 1;
          filtrer_champ_elem(flag_add, d_velocity_tmp[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[0]);
          filtrer_champ_elem(flag_add, d_velocity_tmp[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[1]);
          filtrer_champ_face(flag_add, d_velocity_tmp[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[2]);
        }
      else
        {
          rho_v[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          rho_v[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          rho_v[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          compute_divergence_times_constant(rho_v[0], rho_v[1], rho_v[2], 1., u_div_rho_u);
          u_div_rho_u.echange_espace_virtuel(1);

          multiplier_champ_rho_face_i(true, u_div_rho_u, 1., velocity[0], d_velocity[0]);
          multiplier_champ_rho_face_j(true, u_div_rho_u, 1., velocity[1], d_velocity[1]);
          multiplier_champ_rho_face_k(true, u_div_rho_u, 1., 0., velocity[2], d_velocity[2]);
        }
    }
  else
    {
      if (flag_filtrage_convection_qdm_)
        {
          if (convection_velocity_amont_)
            {
              velocity_convection_op_amont_.calculer(rho_v[0], rho_v[1], rho_v[2],
                                                     velocity[0], velocity[1], velocity[2],
                                                     d_velocity[0], d_velocity[1], d_velocity[2]);

            }
          else if (convection_velocity_quicksharp_)
            {
              velocity_convection_op_quicksharp_.calculer(rho_v[0], rho_v[1], rho_v[2],
                                                          velocity[0], velocity[1], velocity[2],
                                                          d_velocity[0], d_velocity[1], d_velocity[2]);
            }
          else if(convection_velocity_centre2_)
            {
              velocity_convection_op_centre_2_.ajouter_avec_u_div_rhou(rho_v[0], rho_v[1], rho_v[2],
                                                                       velocity[0], velocity[1], velocity[2],
                                                                       d_velocity[0], d_velocity[1], d_velocity[2],
                                                                       u_div_rho_u);
            }
          else
            {
              // Schema centre4 (utilise par defaut)
              velocity_convection_op_.ajouter_avec_u_div_rhou(rho_v[0], rho_v[1], rho_v[2],
                                                              velocity[0], velocity[1], velocity[2],
                                                              d_velocity[0], d_velocity[1], d_velocity[2],
                                                              u_div_rho_u);
            }

          ghost_size_filter = 1 + kernel->ghost_size();
          ghost_size_d_velocity_tmp = max((int)2, ghost_size_filter);
          d_velocity_tmp[0].echange_espace_virtuel(ghost_size_d_velocity_tmp);
          d_velocity_tmp[1].echange_espace_virtuel(ghost_size_d_velocity_tmp);
          d_velocity_tmp[2].echange_espace_virtuel(ghost_size_d_velocity_tmp);

          const int flag_add = 1;
          filtrer_champ_elem(flag_add, d_velocity_tmp[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[0]);
          filtrer_champ_elem(flag_add, d_velocity_tmp[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[1]);
          filtrer_champ_face(flag_add, d_velocity_tmp[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[2]);
        }
      else
        {
          if (convection_velocity_amont_)
            {
              velocity_convection_op_amont_.calculer(rho_v[0], rho_v[1], rho_v[2],
                                                     velocity[0], velocity[1], velocity[2],
                                                     d_velocity[0], d_velocity[1], d_velocity[2]);
            }
          else if (convection_velocity_quicksharp_)
            {
              velocity_convection_op_quicksharp_.calculer(rho_v[0], rho_v[1], rho_v[2],
                                                          velocity[0], velocity[1], velocity[2],
                                                          d_velocity[0], d_velocity[1], d_velocity[2]);
            }
          else if (convection_velocity_centre2_)
            {
              velocity_convection_op_centre_2_.ajouter_avec_u_div_rhou(rho_v[0], rho_v[1], rho_v[2],
                                                                       velocity[0], velocity[1], velocity[2],
                                                                       d_velocity[0], d_velocity[1], d_velocity[2],
                                                                       u_div_rho_u);
            }
          else
            {
              // Schema centre4 (utilise par defaut)
              velocity_convection_op_.ajouter_avec_u_div_rhou(rho_v[0], rho_v[1], rho_v[2],
                                                              velocity[0], velocity[1], velocity[2],
                                                              d_velocity[0], d_velocity[1], d_velocity[2],
                                                              u_div_rho_u);
            }
        }
    }
}

template<class T>
void DNS_QC_double::calculer_turbulent_diffusion_vitesse(FixedVector<IJK_Field_double, 3>& velocity,
                                                         const IJK_Field_double& turbulent_mu_xx,
                                                         const IJK_Field_double& turbulent_mu_xy,
                                                         const IJK_Field_double& turbulent_mu_xz,
                                                         const IJK_Field_double& turbulent_mu_yy,
                                                         const IJK_Field_double& turbulent_mu_yz,
                                                         const IJK_Field_double& turbulent_mu_zz,
                                                         const ArrOfDouble_with_ghost& delta_z,
                                                         const double facteur_delta_x,
                                                         const double facteur_delta_y,
                                                         const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                                         T& kernel,
                                                         FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                                         FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                                         FixedVector<IJK_Field_double, 3>& d_velocity_tmp,
                                                         FixedVector<IJK_Field_double, 3>& d_velocity)
{
  const IJK_Field_double& turbulent_mu_yx = turbulent_mu_xy;
  const IJK_Field_double& turbulent_mu_zx = turbulent_mu_xz;
  const IJK_Field_double& turbulent_mu_zy = turbulent_mu_yz;
  int ghost_size_filter;
  int ghost_size_d_velocity_tmp;
  if (flag_filtrage_turbulent_diffusion_qdm_)
    {
      if (type_velocity_turbulent_diffusion_ == Nom("simple"))
        {
          velocity_turbulent_diffusion_op_simple_.calculer(velocity[0], velocity[1], velocity[2],
                                                           turbulent_mu_xx,
                                                           turbulent_mu_xy,
                                                           turbulent_mu_xz,
                                                           turbulent_mu_yx,
                                                           turbulent_mu_yy,
                                                           turbulent_mu_yz,
                                                           turbulent_mu_zx,
                                                           turbulent_mu_zy,
                                                           turbulent_mu_zz,
                                                           d_velocity_tmp[0], d_velocity_tmp[1], d_velocity_tmp[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("simple_with_transpose"))
        {
          velocity_turbulent_diffusion_op_simple_with_transpose_.calculer(velocity[0], velocity[1], velocity[2],
                                                                          turbulent_mu_xx,
                                                                          turbulent_mu_xy,
                                                                          turbulent_mu_xz,
                                                                          turbulent_mu_yx,
                                                                          turbulent_mu_yy,
                                                                          turbulent_mu_yz,
                                                                          turbulent_mu_zx,
                                                                          turbulent_mu_zy,
                                                                          turbulent_mu_zz,
                                                                          d_velocity_tmp[0], d_velocity_tmp[1], d_velocity_tmp[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("full"))
        {
          velocity[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          velocity[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          velocity[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          compute_divergence(velocity[0], velocity[1], velocity[2], divergence_);
          divergence_.echange_espace_virtuel(1);
          velocity_turbulent_diffusion_op_full_.calculer(velocity[0], velocity[1], velocity[2],
                                                         turbulent_mu_xx,
                                                         turbulent_mu_xy,
                                                         turbulent_mu_xz,
                                                         turbulent_mu_yx,
                                                         turbulent_mu_yy,
                                                         turbulent_mu_yz,
                                                         turbulent_mu_zx,
                                                         turbulent_mu_zy,
                                                         turbulent_mu_zz,
                                                         divergence_,
                                                         d_velocity_tmp[0], d_velocity_tmp[1], d_velocity_tmp[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("simple_anisotropic"))
        {
          velocity_turbulent_diffusion_op_simple_anisotropic_.calculer(velocity[0], velocity[1], velocity[2],
                                                                       turbulent_mu_xx,
                                                                       turbulent_mu_xy,
                                                                       turbulent_mu_xz,
                                                                       turbulent_mu_yx,
                                                                       turbulent_mu_yy,
                                                                       turbulent_mu_yz,
                                                                       turbulent_mu_zx,
                                                                       turbulent_mu_zy,
                                                                       turbulent_mu_zz,
                                                                       d_velocity_tmp[0], d_velocity_tmp[1], d_velocity_tmp[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("simple_with_transpose_anisotropic"))
        {
          velocity_turbulent_diffusion_op_simple_with_transpose_anisotropic_.calculer(velocity[0], velocity[1], velocity[2],
                                                                                      turbulent_mu_xx,
                                                                                      turbulent_mu_xy,
                                                                                      turbulent_mu_xz,
                                                                                      turbulent_mu_yx,
                                                                                      turbulent_mu_yy,
                                                                                      turbulent_mu_yz,
                                                                                      turbulent_mu_zx,
                                                                                      turbulent_mu_zy,
                                                                                      turbulent_mu_zz,
                                                                                      d_velocity_tmp[0], d_velocity_tmp[1], d_velocity_tmp[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("full_anisotropic"))
        {
          velocity[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          velocity[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          velocity[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          compute_divergence(velocity[0], velocity[1], velocity[2], divergence_);
          divergence_.echange_espace_virtuel(1);
          velocity_turbulent_diffusion_op_full_anisotropic_.calculer(velocity[0], velocity[1], velocity[2],
                                                                     turbulent_mu_xx,
                                                                     turbulent_mu_xy,
                                                                     turbulent_mu_xz,
                                                                     turbulent_mu_yx,
                                                                     turbulent_mu_yy,
                                                                     turbulent_mu_yz,
                                                                     turbulent_mu_zx,
                                                                     turbulent_mu_zy,
                                                                     turbulent_mu_zz,
                                                                     divergence_,
                                                                     d_velocity_tmp[0], d_velocity_tmp[1], d_velocity_tmp[2]);
        }
      else
        {
          Cerr << "Unknown velocity turbulent diffusion operator! " << finl;
          Process::exit();
        }

      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_d_velocity_tmp = max((int)2, ghost_size_filter);
      d_velocity_tmp[0].echange_espace_virtuel(ghost_size_d_velocity_tmp);
      d_velocity_tmp[1].echange_espace_virtuel(ghost_size_d_velocity_tmp);
      d_velocity_tmp[2].echange_espace_virtuel(ghost_size_d_velocity_tmp);

      const int flag_add = 1;
      filtrer_champ_elem(flag_add, d_velocity_tmp[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[0]);
      filtrer_champ_elem(flag_add, d_velocity_tmp[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[1]);
      filtrer_champ_face(flag_add, d_velocity_tmp[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[2]);
    }
  else
    {
      if (type_velocity_turbulent_diffusion_ == Nom("simple"))
        {
          velocity_turbulent_diffusion_op_simple_.ajouter(velocity[0], velocity[1], velocity[2],
                                                          turbulent_mu_xx,
                                                          turbulent_mu_xy,
                                                          turbulent_mu_xz,
                                                          turbulent_mu_yx,
                                                          turbulent_mu_yy,
                                                          turbulent_mu_yz,
                                                          turbulent_mu_zx,
                                                          turbulent_mu_zy,
                                                          turbulent_mu_zz,
                                                          d_velocity[0], d_velocity[1], d_velocity[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("simple_with_transpose"))
        {
          velocity_turbulent_diffusion_op_simple_with_transpose_.ajouter(velocity[0], velocity[1], velocity[2],
                                                                         turbulent_mu_xx,
                                                                         turbulent_mu_xy,
                                                                         turbulent_mu_xz,
                                                                         turbulent_mu_yx,
                                                                         turbulent_mu_yy,
                                                                         turbulent_mu_yz,
                                                                         turbulent_mu_zx,
                                                                         turbulent_mu_zy,
                                                                         turbulent_mu_zz,
                                                                         d_velocity[0], d_velocity[1], d_velocity[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("full"))
        {
          velocity[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          velocity[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          velocity[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          compute_divergence(velocity[0], velocity[1], velocity[2], divergence_);
          divergence_.echange_espace_virtuel(1);
          velocity_turbulent_diffusion_op_full_.ajouter(velocity[0], velocity[1], velocity[2],
                                                        turbulent_mu_xx,
                                                        turbulent_mu_xy,
                                                        turbulent_mu_xz,
                                                        turbulent_mu_yx,
                                                        turbulent_mu_yy,
                                                        turbulent_mu_yz,
                                                        turbulent_mu_zx,
                                                        turbulent_mu_zy,
                                                        turbulent_mu_zz,
                                                        divergence_,
                                                        d_velocity[0], d_velocity[1], d_velocity[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("simple_anisotropic"))
        {
          velocity_turbulent_diffusion_op_simple_anisotropic_.ajouter(velocity[0], velocity[1], velocity[2],
                                                                      turbulent_mu_xx,
                                                                      turbulent_mu_xy,
                                                                      turbulent_mu_xz,
                                                                      turbulent_mu_yx,
                                                                      turbulent_mu_yy,
                                                                      turbulent_mu_yz,
                                                                      turbulent_mu_zx,
                                                                      turbulent_mu_zy,
                                                                      turbulent_mu_zz,
                                                                      d_velocity[0], d_velocity[1], d_velocity[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("simple_with_transpose_anisotropic"))
        {
          velocity_turbulent_diffusion_op_simple_with_transpose_anisotropic_.ajouter(velocity[0], velocity[1], velocity[2],
                                                                                     turbulent_mu_xx,
                                                                                     turbulent_mu_xy,
                                                                                     turbulent_mu_xz,
                                                                                     turbulent_mu_yx,
                                                                                     turbulent_mu_yy,
                                                                                     turbulent_mu_yz,
                                                                                     turbulent_mu_zx,
                                                                                     turbulent_mu_zy,
                                                                                     turbulent_mu_zz,
                                                                                     d_velocity[0], d_velocity[1], d_velocity[2]);
        }
      else if (type_velocity_turbulent_diffusion_ == Nom("full_anisotropic"))
        {
          velocity[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          velocity[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          velocity[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          compute_divergence(velocity[0], velocity[1], velocity[2], divergence_);
          divergence_.echange_espace_virtuel(1);
          velocity_turbulent_diffusion_op_full_anisotropic_.ajouter(velocity[0], velocity[1], velocity[2],
                                                                    turbulent_mu_xx,
                                                                    turbulent_mu_xy,
                                                                    turbulent_mu_xz,
                                                                    turbulent_mu_yx,
                                                                    turbulent_mu_yy,
                                                                    turbulent_mu_yz,
                                                                    turbulent_mu_zx,
                                                                    turbulent_mu_zy,
                                                                    turbulent_mu_zz,
                                                                    divergence_,
                                                                    d_velocity[0], d_velocity[1], d_velocity[2]);
        }
      else
        {
          Cerr << "Unknown velocity turbulent diffusion operator! " << finl;
          Process::exit();
        }
    }
}

template<class T>
void DNS_QC_double::calculer_structural_diffusion_vitesse(FixedVector<IJK_Field_double, 3>& velocity,
                                                          const FixedVector<IJK_Field_double, 6>& structural_uu_tensor,
                                                          const ArrOfDouble_with_ghost& delta_z,
                                                          const double facteur_delta_x,
                                                          const double facteur_delta_y,
                                                          const ArrOfDouble_with_ghost& delta_z_pour_delta,
                                                          T& kernel,
                                                          FixedVector<IJK_Field_local_double, 18>& tmp_b,
                                                          FixedVector<IJK_Field_local_double, 18>& tmp_a,
                                                          FixedVector<IJK_Field_double, 3>& d_velocity_tmp,
                                                          FixedVector<IJK_Field_double, 3>& d_velocity)
{
  const IJK_Field_double& structural_uu_xx = structural_uu_tensor[0];
  const IJK_Field_double& structural_uu_xy = structural_uu_tensor[1];
  const IJK_Field_double& structural_uu_xz = structural_uu_tensor[2];
  const IJK_Field_double& structural_uu_yx = structural_uu_tensor[1];
  const IJK_Field_double& structural_uu_yy = structural_uu_tensor[3];
  const IJK_Field_double& structural_uu_yz = structural_uu_tensor[4];
  const IJK_Field_double& structural_uu_zx = structural_uu_tensor[2];
  const IJK_Field_double& structural_uu_zy = structural_uu_tensor[4];
  const IJK_Field_double& structural_uu_zz = structural_uu_tensor[5];
  int ghost_size_filter;
  int ghost_size_d_velocity_tmp;
  if (flag_filtrage_structural_diffusion_qdm_)
    {
      velocity_turbulent_diffusion_op_structural_.calculer(velocity[0], velocity[1], velocity[2],
                                                           structural_uu_xx,
                                                           structural_uu_xy,
                                                           structural_uu_xz,
                                                           structural_uu_yx,
                                                           structural_uu_yy,
                                                           structural_uu_yz,
                                                           structural_uu_zx,
                                                           structural_uu_zy,
                                                           structural_uu_zz,
                                                           d_velocity_tmp[0], d_velocity_tmp[1], d_velocity_tmp[2]);

      ghost_size_filter = 1 + kernel->ghost_size();
      ghost_size_d_velocity_tmp = max((int)2, ghost_size_filter);
      d_velocity_tmp[0].echange_espace_virtuel(ghost_size_d_velocity_tmp);
      d_velocity_tmp[1].echange_espace_virtuel(ghost_size_d_velocity_tmp);
      d_velocity_tmp[2].echange_espace_virtuel(ghost_size_d_velocity_tmp);

      const int flag_add = 1;
      filtrer_champ_elem(flag_add, d_velocity_tmp[0], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[0]);
      filtrer_champ_elem(flag_add, d_velocity_tmp[1], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[1]);
      filtrer_champ_face(flag_add, d_velocity_tmp[2], delta_z, facteur_delta_x, facteur_delta_y, delta_z_pour_delta, kernel, tmp_b, tmp_a, d_velocity[2]);
    }
  else
    {
      velocity_turbulent_diffusion_op_structural_.ajouter(velocity[0], velocity[1], velocity[2],
                                                          structural_uu_xx,
                                                          structural_uu_xy,
                                                          structural_uu_xz,
                                                          structural_uu_yx,
                                                          structural_uu_yy,
                                                          structural_uu_yz,
                                                          structural_uu_zx,
                                                          structural_uu_zy,
                                                          structural_uu_zz,
                                                          d_velocity[0], d_velocity[1], d_velocity[2]);
    }
}

void DNS_QC_double::calculer_diffusion_scalar(const IJK_Field_double& rho,
                                              const IJK_Field_double& turbulent_kappa_x,
                                              const IJK_Field_double& turbulent_kappa_y,
                                              const IJK_Field_double& turbulent_kappa_z,
                                              IJK_Field_double& d_rho,
                                              const IJK_Field_local_double& boundary_flux_kmin,
                                              const IJK_Field_local_double& boundary_flux_kmax)
{
  if (type_scalar_turbulent_diffusion_ == Nom("normal"))
    {
      operateur_diffusion_turbulent_scalar_.ajouter(rho,
                                                    turbulent_kappa_x,
                                                    turbulent_kappa_y,
                                                    turbulent_kappa_z,
                                                    d_rho,
                                                    boundary_flux_kmin, boundary_flux_kmax);
    }
  else if (type_scalar_turbulent_diffusion_ == Nom("anisotropic"))
    {
      operateur_diffusion_turbulent_scalar_anisotropic_.ajouter(rho,
                                                                turbulent_kappa_x,
                                                                turbulent_kappa_y,
                                                                turbulent_kappa_z,
                                                                d_rho,
                                                                boundary_flux_kmin, boundary_flux_kmax);
    }
  else
    {
      Cerr << "Unknown scalar turbulent diffusion operator! " << finl;
      Process::exit();
    }
}

// Perform one sub-step of rk3 for QC algorithm, called 3 times per time step.
// rk_step = 0, 1 or 2
// total_timestep = not the fractionnal timestep !
void DNS_QC_double::rk3_sub_step(const int rk_step, const double total_timestep)
{
  static Stat_Counter_Id cnt_conv_rho = statistiques().new_counter(1, "qc convection rho");
  static Stat_Counter_Id cnt_calcul_temp_rho_rhov = statistiques().new_counter(1, "qc calcul pthermo,t,rho_v,etc");
  static Stat_Counter_Id cnt_diff_v = statistiques().new_counter(1, "qc diffusion vitesse");
  static Stat_Counter_Id cnt_diff_t = statistiques().new_counter(1, "qc diffusion turbulente vitesse");
  static Stat_Counter_Id cnt_diff_struct = statistiques().new_counter(1, "qc diffusion vitesse modele structurel");
  static Stat_Counter_Id cnt_conv_v = statistiques().new_counter(1, "qc convection vitesse");
  static Stat_Counter_Id cnt_mettre_a_jour_v = statistiques().new_counter(1, "qc rk3 update");
  static Stat_Counter_Id cnt_diff_temp = statistiques().new_counter(1, "qc diffusion temperature");
  static Stat_Counter_Id cnt_kappa_t = statistiques().new_counter(1, "qc turbulent diffusivity");
  static Stat_Counter_Id cnt_kappa_struct = statistiques().new_counter(1, "qc turbulent diffusivity structurel");
  static Stat_Counter_Id cnt_prepare_rhs = statistiques().new_counter(1, "qc prepare rhs");
  static Stat_Counter_Id cnt_resoudre_systeme_poisson = statistiques().new_counter(1, "qc resoudre systeme");
  static Stat_Counter_Id cnt_ajouter_grad_p = statistiques().new_counter(1, "qc ajouter grad p");

  const double fractionnal_timestep = compute_fractionnal_timestep_rk3(total_timestep, rk_step);

  if (calcul_2d_)
    {
      velocity_[1].data() = 0;
      force_2d(rho_);
      force_2d(velocity_[0]);
      force_2d(velocity_[2]);
    }

  assert(rk_step>=0 && rk_step<3);
  // F.A deplace a la fin du sous pas de temps, pour le premier dt du premier rk_step c est fait pour post traitrer les champs initiaux.
  // rho_.echange_espace_virtuel(2);
  // velocity_[0].echange_espace_virtuel(1);
  // velocity_[1].echange_espace_virtuel(1);
  // velocity_[2].echange_espace_virtuel(1);

  //
  // Convection du champ de masse volumique
  //
  // L'operateur de convection calcule drho/dt integre sur le volume de controle:
  statistiques().begin_count(cnt_conv_rho);
  if (conv_rho_negligeable_) // cas ou la convection de rho serait negligeable
    {
      d_rho_.data()=0;
    }
  else
    {
      if (convection_rho_centre2_)
        {
          rho_convection_op_centre2_.calculer(rho_, velocity_[0], velocity_[1], velocity_[2], d_rho_);
        }
      else if (convection_rho_amont_)
        {
          rho_convection_op_amont_.calculer(rho_, rho_, rho_, velocity_[0], velocity_[1], velocity_[2], d_rho_, d_rho_, d_rho_);
        }
      else if (convection_rho_centre4_)
        {
          rho_convection_op_centre4_.calculer(rho_, rho_, rho_, velocity_[0], velocity_[1], velocity_[2], d_rho_, d_rho_, d_rho_);
        }
      else
        {
          // Schema quick (utilise par defaut)
          rho_convection_op_.calculer(rho_, velocity_[0], velocity_[1], velocity_[2], d_rho_);
        }
    }
  statistiques().end_count(cnt_conv_rho);

  DebogIJK::verifier("op_conv(rho)", d_rho_);

  if (formulation_velocity_)
    {
      if (turbulent_diffusivity_)
        {
          statistiques().begin_count(cnt_kappa_t);

          if (flag_kappa_vectorial_)
            {
              calculer_turbulent_mu_vector(flag_kappa_anisotropic_,
                                           turbulent_diffusivity_model_, turbulent_diffusivity_model_constant_,
                                           turbulent_diffusivity_vector_coefficients_,  variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_, pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_kappa_filtre_,
                                           turbulent_kappa_filtre_vector_,
                                           turbulent_kappa_vector_, splitting_);
            }
          else
            {
              calculer_turbulent_mu_scalar(flag_kappa_anisotropic_,
                                           turbulent_diffusivity_model_, turbulent_diffusivity_model_constant_, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_, pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_kappa_filtre_, turbulent_kappa_filtre_,
                                           turbulent_kappa_, splitting_);
            }
          statistiques().end_count(cnt_kappa_t);
        }

      if (structural_uscalar_)
        {
          statistiques().begin_count(cnt_kappa_struct);

          calculer_structural_uscalar(structural_uscalar_model_, structural_uscalar_model_constant_,
                                      structural_uscalar_vector_coefficients_,
                                      rho_,
                                      velocity_, velocity_filtre_,
                                      rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                      delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                      facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                      kernel_, tmp_b_, tmp_a_, structural_uscalar_tmp_vector_,
                                      flag_structural_uscalar_filtre_, structural_uscalar_filtre_vector_,
                                      structural_uscalar_vector_);

          statistiques().end_count(cnt_kappa_struct);
        }

      if (flag_kappa_vectorial_)
        {
          modification_modele_dynamic_uscalar_vector(flag_kappa_anisotropic_,
                                                     turbulent_diffusivity_dynamic_type_, structural_uscalar_dynamic_type_,
                                                     velocity_, velocity_filtre_,
                                                     rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                                     delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                     facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                     kernel_, tmp_b_, tmp_a_,
                                                     constante_modele_, ml_,
                                                     turbulent_diffusivity_, turbulent_kappa_vector_, turbulent_kappa_filtre_vector_,
                                                     structural_uscalar_, structural_uscalar_vector_, structural_uscalar_filtre_vector_);
        }
      else
        {
          modification_modele_dynamic_uscalar_scalar(flag_kappa_anisotropic_,
                                                     turbulent_diffusivity_dynamic_type_, structural_uscalar_dynamic_type_,
                                                     velocity_, velocity_filtre_,
                                                     rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                                     delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                     facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                     kernel_, tmp_b_, tmp_a_,
                                                     constante_modele_, ml_,
                                                     turbulent_diffusivity_, turbulent_kappa_, turbulent_kappa_filtre_,
                                                     structural_uscalar_, structural_uscalar_vector_, structural_uscalar_filtre_vector_);
        }

      if (turbulent_diffusivity_)
        {
          statistiques().begin_count(cnt_kappa_t);
          if (flag_kappa_vectorial_)
            {
              turbulent_kappa_vector_[0].echange_espace_virtuel(1);
              turbulent_kappa_vector_[1].echange_espace_virtuel(1);
              turbulent_kappa_vector_[2].echange_espace_virtuel(1);

              // on impose le flux au bord a zero dans tous les cas
              calculer_flux_thermique_bord(rho_, 0.,
                                           0, turbulent_kappa_vector_[2], flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           rho_paroi_impose_kmin_, boundary_flux_kmin_, 0 /* boundary kmin */);
              calculer_flux_thermique_bord(rho_, 0.,
                                           0, turbulent_kappa_vector_[2], flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           rho_paroi_impose_kmax_, boundary_flux_kmax_, 1 /* boundary kmax */);

              rho_.echange_espace_virtuel(1);
              turbulent_kappa_vector_[0].echange_espace_virtuel(1);
              turbulent_kappa_vector_[1].echange_espace_virtuel(1);
              turbulent_kappa_vector_[2].echange_espace_virtuel(1);

              calculer_diffusion_scalar(rho_,
                                        turbulent_kappa_vector_[0],
                                        turbulent_kappa_vector_[1],
                                        turbulent_kappa_vector_[2],
                                        d_rho_,
                                        boundary_flux_kmin_, boundary_flux_kmax_);
            }
          else
            {
              turbulent_kappa_.echange_espace_virtuel(1);

              // on impose le flux au bord a zero dans tous les cas
              calculer_flux_thermique_bord(rho_, 0.,
                                           0, turbulent_kappa_, flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           rho_paroi_impose_kmin_, boundary_flux_kmin_, 0 /* boundary kmin */);
              calculer_flux_thermique_bord(rho_, 0.,
                                           0, turbulent_kappa_, flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           rho_paroi_impose_kmax_, boundary_flux_kmax_, 1 /* boundary kmax */);

              rho_.echange_espace_virtuel(1);
              turbulent_kappa_.echange_espace_virtuel(1);

              calculer_diffusion_scalar(rho_,
                                        turbulent_kappa_,
                                        turbulent_kappa_,
                                        turbulent_kappa_,
                                        d_rho_,
                                        boundary_flux_kmin_, boundary_flux_kmax_);
            }
          statistiques().end_count(cnt_kappa_t);

          DebogIJK::verifier("rho", rho_);
          DebogIJK::verifier("op_conv(rho)_kappa", d_rho_);
        }

      if (structural_uscalar_)
        {
          statistiques().begin_count(cnt_kappa_struct);

          structural_uscalar_vector_[0].echange_espace_virtuel(1);
          structural_uscalar_vector_[1].echange_espace_virtuel(1);
          structural_uscalar_vector_[2].echange_espace_virtuel(1);

          // on impose le flux au bord a zero dans tous les cas
          calculer_flux_thermique_bord(rho_, 0.,
                                       0, turbulent_kappa_, flag_kappa_anisotropic_,
                                       0, structural_uscalar_vector_[2],
                                       rho_paroi_impose_kmin_, boundary_flux_kmin_, 0 /* boundary kmin */);
          calculer_flux_thermique_bord(rho_, 0.,
                                       0, turbulent_kappa_, flag_kappa_anisotropic_,
                                       0, structural_uscalar_vector_[2],
                                       rho_paroi_impose_kmax_, boundary_flux_kmax_, 1 /* boundary kmax */);

          rho_.echange_espace_virtuel(1);
          structural_uscalar_vector_[0].echange_espace_virtuel(1);
          structural_uscalar_vector_[1].echange_espace_virtuel(1);
          structural_uscalar_vector_[2].echange_espace_virtuel(1);

          DebogIJK::verifier("rho", rho_);
          DebogIJK::verifier("structural_uscalar_x", structural_uscalar_vector_[0]);
          DebogIJK::verifier("structural_uscalar_y", structural_uscalar_vector_[1]);
          DebogIJK::verifier("structural_uscalar_z", structural_uscalar_vector_[2]);

          operateur_diffusion_temperature_structural_.ajouter(rho_,
                                                              structural_uscalar_vector_[0],
                                                              structural_uscalar_vector_[1],
                                                              structural_uscalar_vector_[2],
                                                              d_rho_,
                                                              boundary_flux_kmin_, boundary_flux_kmax_);

          statistiques().end_count(cnt_kappa_struct);

          DebogIJK::verifier("op_conv(rho)_struct", d_rho_);
        }
    }

  {
    const int kmax = d_rho_.nk();
    for (int k = 0; k < kmax; k++)
      {
        // division par le volume de controle
        mass_solver_scalar(d_rho_, delta_z_local_, k);
        // calcul de rho_ au sous-pas de temps suivant, et mise a jour F_rho_ (valeur intermediaire de l'algo rk3)
        runge_kutta3_update(d_rho_, RK3_F_rho_, rho_, rk_step, k, total_timestep);
      }
  }
  DebogIJK::verifier("rk3 update rho", rho_);

  // On a besoin d'une epaisseur de 1 sur rho (pour calculer rho aux faces et ensuite
  //  pour calculer l'epaisseur 1 de la temperature et de mu/lambda.
  rho_.echange_espace_virtuel(1);

  statistiques().begin_count(cnt_calcul_temp_rho_rhov);
  // Resoudre EDO Pthermo: calcul de P_th_final et de la temperature en fonction de rho_
  const double P_th_initial = P_thermodynamique_;
  double P_th_final;
  double d_Pth_divise_par_gammamoins1;

  if (formulation_favre_)
    {
      if (turbulent_diffusivity_)
        {
          statistiques().begin_count(cnt_kappa_t);

          if (flag_kappa_vectorial_)
            {
              calculer_turbulent_mu_vector(flag_kappa_anisotropic_,
                                           turbulent_diffusivity_model_, turbulent_diffusivity_model_constant_,
                                           turbulent_diffusivity_vector_coefficients_, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_kappa_filtre_,
                                           turbulent_kappa_filtre_vector_,
                                           turbulent_kappa_vector_, splitting_);
            }
          else
            {
              calculer_turbulent_mu_scalar(flag_kappa_anisotropic_,
                                           turbulent_diffusivity_model_, turbulent_diffusivity_model_constant_, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_kappa_filtre_, turbulent_kappa_filtre_,
                                           turbulent_kappa_, splitting_);
            }
          statistiques().end_count(cnt_kappa_t);
        }

      if (structural_uscalar_)
        {
          statistiques().begin_count(cnt_kappa_struct);

          calculer_structural_uscalar(structural_uscalar_model_, structural_uscalar_model_constant_,
                                      structural_uscalar_vector_coefficients_,
                                      rho_,
                                      velocity_, velocity_filtre_,
                                      temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                      delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                      facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                      kernel_, tmp_b_, tmp_a_, structural_uscalar_tmp_vector_,
                                      flag_structural_uscalar_filtre_, structural_uscalar_filtre_vector_,
                                      structural_uscalar_vector_);

          statistiques().end_count(cnt_kappa_struct);
        }

      if (flag_kappa_vectorial_)
        {
          modification_modele_dynamic_uscalar_vector(flag_kappa_anisotropic_,
                                                     turbulent_diffusivity_dynamic_type_, structural_uscalar_dynamic_type_,
                                                     velocity_, velocity_filtre_,
                                                     temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                                     delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                     facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                     kernel_, tmp_b_, tmp_a_,
                                                     constante_modele_, ml_,
                                                     turbulent_diffusivity_, turbulent_kappa_vector_, turbulent_kappa_filtre_vector_,
                                                     structural_uscalar_, structural_uscalar_vector_, structural_uscalar_filtre_vector_);
        }
      else
        {
          modification_modele_dynamic_uscalar_scalar(flag_kappa_anisotropic_,
                                                     turbulent_diffusivity_dynamic_type_, structural_uscalar_dynamic_type_,
                                                     velocity_, velocity_filtre_,
                                                     temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                                     delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                     facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                     kernel_, tmp_b_, tmp_a_,
                                                     constante_modele_, ml_,
                                                     turbulent_diffusivity_, turbulent_kappa_, turbulent_kappa_filtre_,
                                                     structural_uscalar_, structural_uscalar_vector_, structural_uscalar_filtre_vector_);
        }

      if (turbulent_diffusivity_)
        {
          statistiques().begin_count(cnt_kappa_t);

          if (flag_kappa_vectorial_)
            {
              // we multiply by rho cp = rho (r/Pth) (Pth CP/r):
              //  * r/Pth: because we use grad T instead of grad 1/rho
              //  * Pth Cp/r: because div_lambda_grad_T will be multiplied by r/Cp/Pth
              //              but this should not
              //  * rho: required in the favre formulation
              multiplier_champ(rho_, Cp_gaz_, turbulent_kappa_vector_[0]);
              multiplier_champ(rho_, Cp_gaz_, turbulent_kappa_vector_[1]);
              multiplier_champ(rho_, Cp_gaz_, turbulent_kappa_vector_[2]);

              turbulent_kappa_vector_[0].echange_espace_virtuel(1);
              turbulent_kappa_vector_[1].echange_espace_virtuel(1);
              turbulent_kappa_vector_[2].echange_espace_virtuel(1);
            }
          else
            {
              // we multiply by rho cp = rho (r/Pth) (Pth CP/r):
              //  * r/Pth: because we use grad T instead of grad 1/rho
              //  * Pth Cp/r: because div_lambda_grad_T will be multiplied by r/Cp/Pth
              //              but this should not
              //  * rho: required in the favre formulation
              multiplier_champ(rho_, Cp_gaz_, turbulent_kappa_);

              turbulent_kappa_.echange_espace_virtuel(1);
            }
          statistiques().end_count(cnt_kappa_t);
        }

      if (structural_uscalar_)
        {
          statistiques().begin_count(cnt_kappa_struct);

          // we multiply by rho cp = rho (r/Pth) (Pth CP/r):
          //  * r/Pth: because we use grad T instead of grad 1/rho
          //  * Pth Cp/r: because div_lambda_grad_T will be multiplied by r/Cp/Pth
          //              but this should not
          //  * rho: required in the favre formulation
          multiplier_champ_rho_face_i(rho_, Cp_gaz_, structural_uscalar_vector_[0]);
          multiplier_champ_rho_face_j(rho_, Cp_gaz_, structural_uscalar_vector_[1]);
          multiplier_champ_rho_face_k(rho_, Cp_gaz_, rho_paroi_impose_kmin_, structural_uscalar_vector_[2]);

          structural_uscalar_vector_[0].echange_espace_virtuel(1);
          structural_uscalar_vector_[1].echange_espace_virtuel(1);
          structural_uscalar_vector_[2].echange_espace_virtuel(1);

          statistiques().end_count(cnt_kappa_struct);
        }
    }

  if (flag_kappa_vectorial_)
    {
      calcul_p_thermo_et_bilan(rho_, temperature_, turbulent_diffusivity_, turbulent_kappa_vector_[2], flag_kappa_anisotropic_, structural_uscalar_, structural_uscalar_vector_[2], P_th_initial, P_th_final, fractionnal_timestep, d_Pth_divise_par_gammamoins1);
    }
  else
    {
      calcul_p_thermo_et_bilan(rho_, temperature_, turbulent_diffusivity_, turbulent_kappa_, flag_kappa_anisotropic_, structural_uscalar_, structural_uscalar_vector_[2], P_th_initial, P_th_final, fractionnal_timestep, d_Pth_divise_par_gammamoins1);
    }
  P_thermodynamique_ = P_th_final;

  const double Pth_sur_R = P_thermodynamique_ / constante_specifique_gaz_;
  rho_paroi_impose_kmin_ = Pth_sur_R / T_paroi_impose_kmin_;
  rho_paroi_impose_kmax_ = Pth_sur_R / T_paroi_impose_kmax_;

  // Mise a jour de temperature, mu et lambda en fonction de rho, sur une epaisseur de 1
  calculer_temperature_mu_lambda_air(P_th_final, constante_specifique_gaz_,
                                     rho_, temperature_, molecular_mu_, molecular_lambda_, 1 /* epaisseur de joint */);

  DebogIJK::verifier("temperature", temperature_);

  // derivee_en_temps_impl_p1 Navier Stokes
  calculer_rho_v(rho_, velocity_, rho_v_);
  DebogIJK::verifier("rho_v_X", rho_v_[0]);
  DebogIJK::verifier("rho_v_Y", rho_v_[1]);
  DebogIJK::verifier("rho_v_Z", rho_v_[2]);
  statistiques().end_count(cnt_calcul_temp_rho_rhov);

  // Pour l'operateur de convection, on a besoin d'une epaisseur 2 sur la quantite convectee:
  rho_v_[0].echange_espace_virtuel(2);
  rho_v_[1].echange_espace_virtuel(2);
  rho_v_[2].echange_espace_virtuel(2);
  // on a besoin de l'epaisseur 1 sur mu et lambda

  // Le coefficient de diffusion est mu => resultat homogene a d(rho_v)/dt
  statistiques().begin_count(cnt_diff_v);
  if(diff_qdm_negligeable_) // si la diffusion est negligable
    {
      d_velocity_[0].data()=0;
      d_velocity_[1].data()=0;
      d_velocity_[2].data()=0;
    }
  else
    {
      if (type_velocity_diffusion_ == Nom("simple"))
        {
          velocity_diffusion_op_simple_.calculer(velocity_[0], velocity_[1], velocity_[2],
                                                 molecular_mu_,
                                                 d_velocity_[0], d_velocity_[1], d_velocity_[2]);
        }
      else if (type_velocity_diffusion_ == Nom("simple_with_transpose"))
        {
          velocity_diffusion_op_simple_with_transpose_.calculer(velocity_[0], velocity_[1], velocity_[2],
                                                                molecular_mu_,
                                                                d_velocity_[0], d_velocity_[1], d_velocity_[2]);
        }
      else if (type_velocity_diffusion_ == Nom("full"))
        {
          velocity_[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
          velocity_[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
          velocity_[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/
          compute_divergence(velocity_[0], velocity_[1], velocity_[2], divergence_);
          divergence_.echange_espace_virtuel(1);
          velocity_diffusion_op_full_.calculer(velocity_[0], velocity_[1], velocity_[2],
                                               molecular_mu_,
                                               divergence_,
                                               d_velocity_[0], d_velocity_[1], d_velocity_[2]);
        }
      else if (type_velocity_diffusion_ == Nom("none"))
        {
          d_velocity_[0].data()=0;
          d_velocity_[1].data()=0;
          d_velocity_[2].data()=0;
        }
      else
        {
          Cerr << "Unknown velocity diffusion operator! " << finl;
          Process::exit();
        }
    }
  statistiques().end_count(cnt_diff_v);
  DebogIJK::verifier("op_diff(velocity)X", d_velocity_[0]);
  DebogIJK::verifier("op_diff(velocity)Y", d_velocity_[1]);
  DebogIJK::verifier("op_diff(velocity)Z", d_velocity_[2]);

  if (formulation_favre_)
    {
      if (turbulent_viscosity_)
        {
          statistiques().begin_count(cnt_diff_t);

          if (flag_nu_tensorial_)
            {
              calculer_turbulent_mu_tensor(flag_nu_anisotropic_,
                                           turbulent_viscosity_model_, turbulent_viscosity_model_constant_,
                                           turbulent_viscosity_tensor_coefficients_, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_mu_filtre_,
                                           turbulent_mu_filtre_tensor_,
                                           turbulent_mu_tensor_, splitting_);
            }
          else
            {
              calculer_turbulent_mu_scalar(flag_nu_anisotropic_,
                                           turbulent_viscosity_model_, turbulent_viscosity_model_constant_, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_,  Re_tau_fr_, Re_tau_ch_, pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_mu_filtre_, turbulent_mu_filtre_,
                                           turbulent_mu_, splitting_);
            }
          statistiques().end_count(cnt_diff_t);
        }

      if (structural_uu_)
        {
          statistiques().begin_count(cnt_diff_struct);
          calculer_structural_uu(structural_uu_model_, structural_uu_model_constant_,
                                 structural_uu_tensor_coefficients_,
                                 rho_,
                                 velocity_, velocity_filtre_,
                                 delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                 facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                 kernel_, tmp_b_, tmp_a_, structural_uu_tmp_tensor_,
                                 flag_structural_uu_filtre_, structural_uu_filtre_tensor_,
                                 structural_uu_tensor_);
          statistiques().end_count(cnt_diff_struct);
        }

      if (flag_nu_tensorial_)
        {
          modification_modele_dynamic_uu_tensor(flag_nu_anisotropic_,
                                                turbulent_viscosity_dynamic_type_, structural_uu_dynamic_type_,
                                                velocity_, velocity_filtre_,
                                                temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                                delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                kernel_, tmp_b_, tmp_a_,
                                                constante_modele_, ml_,
                                                turbulent_viscosity_, turbulent_mu_tensor_, turbulent_mu_filtre_tensor_,
                                                structural_uu_, structural_uu_tensor_, structural_uu_filtre_tensor_);
        }
      else
        {
          modification_modele_dynamic_uu_scalar(flag_nu_anisotropic_,
                                                turbulent_viscosity_dynamic_type_, structural_uu_dynamic_type_,
                                                velocity_, velocity_filtre_,
                                                temperature_, temperature_filtre_, T_paroi_impose_kmin_, T_paroi_impose_kmax_,
                                                delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                kernel_, tmp_b_, tmp_a_,
                                                constante_modele_, ml_,
                                                turbulent_viscosity_, turbulent_mu_, turbulent_mu_filtre_,
                                                structural_uu_, structural_uu_tensor_, structural_uu_filtre_tensor_);
        }

      if (turbulent_viscosity_)
        {
          statistiques().begin_count(cnt_diff_t);

          if (flag_nu_tensorial_)
            {
              multiplier_champ(rho_, turbulent_mu_tensor_[0]);
              multiplier_champ(rho_, turbulent_mu_tensor_[1]);
              multiplier_champ(rho_, turbulent_mu_tensor_[2]);
              multiplier_champ(rho_, turbulent_mu_tensor_[3]);
              multiplier_champ(rho_, turbulent_mu_tensor_[4]);
              multiplier_champ(rho_, turbulent_mu_tensor_[5]);

              turbulent_mu_tensor_[0].echange_espace_virtuel(1);
              turbulent_mu_tensor_[1].echange_espace_virtuel(1);
              turbulent_mu_tensor_[2].echange_espace_virtuel(1);
              turbulent_mu_tensor_[3].echange_espace_virtuel(1);
              turbulent_mu_tensor_[4].echange_espace_virtuel(1);
              turbulent_mu_tensor_[5].echange_espace_virtuel(1);

              calculer_turbulent_diffusion_vitesse(velocity_,
                                                   turbulent_mu_tensor_[0],
                                                   turbulent_mu_tensor_[1],
                                                   turbulent_mu_tensor_[2],
                                                   turbulent_mu_tensor_[3],
                                                   turbulent_mu_tensor_[4],
                                                   turbulent_mu_tensor_[5],
                                                   delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                   kernel_, tmp_b_, tmp_a_,
                                                   d_velocity_tmp_,
                                                   d_velocity_);
            }
          else
            {
              multiplier_champ(rho_, turbulent_mu_);

              turbulent_mu_.echange_espace_virtuel(1);

              calculer_turbulent_diffusion_vitesse(velocity_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                   kernel_, tmp_b_, tmp_a_,
                                                   d_velocity_tmp_,
                                                   d_velocity_);
            }
          statistiques().end_count(cnt_diff_t);

          DebogIJK::verifier("op_difft(velocity)X", d_velocity_[0]);
          DebogIJK::verifier("op_difft(velocity)Y", d_velocity_[1]);
          DebogIJK::verifier("op_difft(velocity)Z", d_velocity_[2]);
        }

      if (structural_uu_)
        {
          statistiques().begin_count(cnt_diff_struct);

          multiplier_champ(rho_, structural_uu_tensor_[0]);
          multiplier_champ_rho_arete_ij(rho_, structural_uu_tensor_[1]);
          multiplier_champ_rho_arete_ik(rho_, rho_paroi_impose_kmin_, structural_uu_tensor_[2]);
          multiplier_champ(rho_, structural_uu_tensor_[3]);
          multiplier_champ_rho_arete_jk(rho_, rho_paroi_impose_kmin_, structural_uu_tensor_[4]);
          multiplier_champ(rho_, structural_uu_tensor_[5]);

          structural_uu_tensor_[0].echange_espace_virtuel(1);
          structural_uu_tensor_[1].echange_espace_virtuel(1);
          structural_uu_tensor_[2].echange_espace_virtuel(1);
          structural_uu_tensor_[3].echange_espace_virtuel(1);
          structural_uu_tensor_[4].echange_espace_virtuel(1);
          structural_uu_tensor_[5].echange_espace_virtuel(1);

          calculer_structural_diffusion_vitesse(velocity_,
                                                structural_uu_tensor_,
                                                delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                kernel_, tmp_b_, tmp_a_,
                                                d_velocity_tmp_,
                                                d_velocity_);

          statistiques().end_count(cnt_diff_struct);

          DebogIJK::verifier("op_diffstruct(velocity)X", d_velocity_[0]);
          DebogIJK::verifier("op_diffstruct(velocity)Y", d_velocity_[1]);
          DebogIJK::verifier("op_diffstruct(velocity)Z", d_velocity_[2]);
        }
    }

  // On transporte le champ rho_v par le champ v
  statistiques().begin_count(cnt_conv_v);
  if (conv_qdm_negligeable_ || flag_convection_qdm_sans_rho_)
    {
      u_div_rho_u_.data()=0;
    }
  else // sinon on calcule rho u div(u) comme div(rho u u) - u div(rho u)
    {
      calculer_convection_vitesse(rho_v_,
                                  velocity_,
                                  delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                  kernel_, tmp_b_, tmp_a_,
                                  d_velocity_tmp_,
                                  d_velocity_,
                                  u_div_rho_u_);
    }
  statistiques().end_count(cnt_conv_v);
  DebogIJK::verifier("op_conv(velocity)X", d_velocity_[0]);
  DebogIJK::verifier("op_conv(velocity)Y", d_velocity_[1]);
  DebogIJK::verifier("op_conv(velocity)Z", d_velocity_[2]);


  if(flag_oscillating_boundary)
    {
      int dir = DIRECTION_J;
      IJK_Field_double& dv = d_velocity_[dir];
      const int ni = dv.ni();
      const int nj = dv.nj();
      const int kmax = d_velocity_[dir].nk();
      const double volumic_force =  frequency_oscillating_boundary_ * amplitude_oscillating_boundary_ * cos(current_time_ * frequency_oscillating_boundary_);
      for(int k = 0; k < kmax; k++)
        {
          for (int j = 0; j < nj; j++)
            {
              for (int i = 0; i < ni; i++)
                {
                  const double volume = get_channel_control_volume(dv, k, delta_z_local_);
                  const double force = volumic_force * volume;
                  dv(i, j, k) += force;
                }
            }
        }
      /* Cout << "Volumic force at edge: " << volumic_force * get_channel_control_volume(dv, 0, delta_z_local_) */
      /* << " Volumic force at middle: " << volumic_force * get_channel_control_volume(dv, kmax/2, delta_z_local_) */
      /* << finl; */
    }




  // force_oscillating_boundary(
  //     velocity_[DIRECTION_J],
  //     current_time_,
  //     amplitude_oscillating_boundary_2_,
  //     frequency_oscillating_boundary_2_
  // );
  statistiques().begin_count(cnt_mettre_a_jour_v);
  for (int dir = 0; dir < 3; dir++)
    {
      const int kmax = d_velocity_[dir].nk();
      for (int k = 0; k < kmax; k++)
        {
          // Le terme source est ajouter avant solveur masse, donc homogene a des Newton:
          if (dir == DIRECTION_I)
            {
              const double force_volumique = terme_source_acceleration_;
              IJK_Field_double& dvx = d_velocity_[dir];
              const double volume = get_channel_control_volume(dvx, k, delta_z_local_);
              const double f = force_volumique * volume;
              const int ni = dvx.ni();
              const int nj = dvx.nj();

              for (int j = 0; j < nj; j++)
                {
                  for (int i = 0; i < ni; i++)
                    {
                      dvx(i,j,k) += f;
                    }
                }
            }

          density_solver_with_rho(d_velocity_[dir], rho_, delta_z_local_, k);
        }
    }

  // On transporte le champ v par le champ v
  statistiques().begin_count(cnt_conv_v);
  if (flag_convection_qdm_sans_rho_ && (!conv_qdm_negligeable_))
    {
      // on calcule u div(u) comme div(u u) - u div(u)
      calculer_convection_vitesse(velocity_,
                                  velocity_,
                                  delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                  kernel_, tmp_b_, tmp_a_,
                                  d_velocity_tmp_,
                                  d_velocity_,
                                  u_div_rho_u_);
    }
  statistiques().end_count(cnt_conv_v);
  DebogIJK::verifier("op_conv(velocity)X", d_velocity_[0]);
  DebogIJK::verifier("op_conv(velocity)Y", d_velocity_[1]);
  DebogIJK::verifier("op_conv(velocity)Z", d_velocity_[2]);

  if (formulation_velocity_)
    {
      if (turbulent_viscosity_)
        {
          statistiques().begin_count(cnt_diff_t);

          if (flag_nu_tensorial_)
            {
              calculer_turbulent_mu_tensor(flag_nu_anisotropic_,
                                           turbulent_viscosity_model_, turbulent_viscosity_model_constant_,
                                           turbulent_viscosity_tensor_coefficients_, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_mu_filtre_,
                                           turbulent_mu_filtre_tensor_,
                                           turbulent_mu_tensor_, splitting_);
            }
          else
            {
              calculer_turbulent_mu_scalar(flag_nu_anisotropic_,
                                           turbulent_viscosity_model_, turbulent_viscosity_model_constant_, variation_cste_modele_fonctionnel_, smoothing_center_fr_, smoothing_factor_fr_, Re_tau_fr_, Re_tau_ch_,  pond_fr_, pond_ch_, center_constant_, Lz_tot_,
                                           velocity_, velocity_filtre_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                           delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                           facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                           kernel_, tmp_b_, tmp_a_,
                                           flag_turbulent_mu_filtre_, turbulent_mu_filtre_,
                                           turbulent_mu_, splitting_);
            }

          statistiques().end_count(cnt_diff_t);
        }

      if (structural_uu_)
        {
          statistiques().begin_count(cnt_diff_struct);
          Cerr << "avant calculer_structural_uu" << finl;
          calculer_structural_uu(structural_uu_model_, structural_uu_model_constant_,
                                 structural_uu_tensor_coefficients_,
                                 rho_,
                                 velocity_, velocity_filtre_,
                                 delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                 facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                 kernel_, tmp_b_, tmp_a_, structural_uu_tmp_tensor_,
                                 flag_structural_uu_filtre_, structural_uu_filtre_tensor_,
                                 structural_uu_tensor_);

          statistiques().end_count(cnt_diff_struct);
        }

      if (flag_nu_tensorial_)
        {
          modification_modele_dynamic_uu_tensor(flag_nu_anisotropic_,
                                                turbulent_viscosity_dynamic_type_, structural_uu_dynamic_type_,
                                                velocity_, velocity_filtre_,
                                                rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                                delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                kernel_, tmp_b_, tmp_a_,
                                                constante_modele_, ml_,
                                                turbulent_viscosity_, turbulent_mu_tensor_, turbulent_mu_filtre_tensor_,
                                                structural_uu_, structural_uu_tensor_, structural_uu_filtre_tensor_);
        }
      else
        {
          modification_modele_dynamic_uu_scalar(flag_nu_anisotropic_,
                                                turbulent_viscosity_dynamic_type_, structural_uu_dynamic_type_,
                                                velocity_, velocity_filtre_,
                                                rho_, rho_filtre_, rho_paroi_impose_kmin_, rho_paroi_impose_kmax_,
                                                delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                facteur_delta_filtre_x_, facteur_delta_filtre_y_, delta_z_local_pour_delta_filtre_,
                                                kernel_, tmp_b_, tmp_a_,
                                                constante_modele_, ml_,
                                                turbulent_viscosity_, turbulent_mu_, turbulent_mu_filtre_,
                                                structural_uu_, structural_uu_tensor_, structural_uu_filtre_tensor_);
        }

      if (turbulent_viscosity_)
        {
          statistiques().begin_count(cnt_diff_t);

          if (flag_nu_tensorial_)
            {
              turbulent_mu_tensor_[0].echange_espace_virtuel(1);
              turbulent_mu_tensor_[1].echange_espace_virtuel(1);
              turbulent_mu_tensor_[2].echange_espace_virtuel(1);
              turbulent_mu_tensor_[3].echange_espace_virtuel(1);
              turbulent_mu_tensor_[4].echange_espace_virtuel(1);
              turbulent_mu_tensor_[5].echange_espace_virtuel(1);

              calculer_turbulent_diffusion_vitesse(velocity_,
                                                   turbulent_mu_tensor_[0],
                                                   turbulent_mu_tensor_[1],
                                                   turbulent_mu_tensor_[2],
                                                   turbulent_mu_tensor_[3],
                                                   turbulent_mu_tensor_[4],
                                                   turbulent_mu_tensor_[5],
                                                   delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                   kernel_, tmp_b_, tmp_a_,
                                                   d_velocity_tmp_,
                                                   d_velocity_);
            }
          else
            {
              turbulent_mu_.echange_espace_virtuel(1);

              calculer_turbulent_diffusion_vitesse(velocity_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   turbulent_mu_,
                                                   delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                   kernel_, tmp_b_, tmp_a_,
                                                   d_velocity_tmp_,
                                                   d_velocity_);
            }

          statistiques().end_count(cnt_diff_t);

          DebogIJK::verifier("op_difft(velocity)X", d_velocity_[0]);
          DebogIJK::verifier("op_difft(velocity)Y", d_velocity_[1]);
          DebogIJK::verifier("op_difft(velocity)Z", d_velocity_[2]);
        }

      if (structural_uu_)
        {
          statistiques().begin_count(cnt_diff_struct);

          structural_uu_tensor_[0].echange_espace_virtuel(1);
          structural_uu_tensor_[1].echange_espace_virtuel(1);
          structural_uu_tensor_[2].echange_espace_virtuel(1);
          structural_uu_tensor_[3].echange_espace_virtuel(1);
          structural_uu_tensor_[4].echange_espace_virtuel(1);
          structural_uu_tensor_[5].echange_espace_virtuel(1);

          calculer_structural_diffusion_vitesse(velocity_,
                                                structural_uu_tensor_,
                                                delta_z_local_, facteur_delta_x_, facteur_delta_y_, delta_z_local_pour_delta_,
                                                kernel_, tmp_b_, tmp_a_,
                                                d_velocity_tmp_,
                                                d_velocity_);

          statistiques().end_count(cnt_diff_struct);

          DebogIJK::verifier("op_diffstruct(velocity)X", d_velocity_[0]);
          DebogIJK::verifier("op_diffstruct(velocity)Y", d_velocity_[1]);
          DebogIJK::verifier("op_diffstruct(velocity)Z", d_velocity_[2]);
        }
    }

  // diffusion creates non zero velocity at walls, restore zero
  force_zero_normal_velocity_on_walls(d_velocity_[2]);

  for (int dir = 0; dir < 3; dir++)
    {
      const int kmax = d_velocity_[dir].nk();
      for (int k = 0; k < kmax; k++)
        {
          mass_solver_scalar(d_velocity_[dir], delta_z_local_, k);

          runge_kutta3_update(d_velocity_[dir], RK3_F_velocity_[dir], velocity_[dir], rk_step, k, total_timestep);
        }
    }
  statistiques().end_count(cnt_mettre_a_jour_v);
  DebogIJK::verifier("rk3update(velocity)X", velocity_[0]);
  DebogIJK::verifier("rk3update(velocity)Y", velocity_[1]);
  DebogIJK::verifier("rk3update(velocity)Z", velocity_[2]);
  // --------------------------------------
  // Projection sur div(u)=quelque chose
  // --------------------------------------
  // vpoint -= gradP / rho
  // secmem = - R / (Cp * Pth) * (div((lambda*grad(T))*volume+source)  - dpth2 / gamma_moins1 * volume)


  // B secmem =  [ ( ( R/(Cp*Pth)*(tab_W(i) + dpth2/gamma_moins_1*volumes(i) )) - div(v)) / dt - div(vpoint) ]

  // trouver P tel que div ( v + dt * vpoint - dt * 1/rho*grad(P) ) = (R/(cp*Pth)...)
  // equivalent a
  //   trouver P tel que div ( dt * 1/rho*gradP ) = div( v + dt * vpoint - (R/(cp*Pth)...) )
  // A  trouver P tel que div ( 1/rho*gradP ) = div( v/dt +  vpoint - (R/(cp*Pth)...)/dt )

  // calculer ici div(lambda*grad(T))*volume)

  calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmin_,
                               0, turbulent_kappa_, flag_kappa_anisotropic_,
                               0, structural_uscalar_vector_[2],
                               T_paroi_impose_kmin_, boundary_flux_kmin_, 0 /* boundary kmin */);
  calculer_flux_thermique_bord(temperature_, lambda_de_t_paroi_kmax_,
                               0, turbulent_kappa_, flag_kappa_anisotropic_,
                               0, structural_uscalar_vector_[2],
                               T_paroi_impose_kmax_, boundary_flux_kmax_, 1 /* boundary kmax */);

  temperature_.echange_espace_virtuel(1);
  molecular_lambda_.echange_espace_virtuel(1);
  DebogIJK::verifier("temp", temperature_);
  DebogIJK::verifier("lambda", molecular_lambda_);

  statistiques().begin_count(cnt_diff_temp);
  if (diff_temp_negligeable_) // si diffusion negligeable
    {
      div_lambda_grad_T_volume_.data()=0;
    }
  else
    {
      operateur_diffusion_temperature_.calculer(temperature_,
                                                molecular_lambda_,
                                                div_lambda_grad_T_volume_,
                                                boundary_flux_kmin_, boundary_flux_kmax_);
    }
  statistiques().end_count(cnt_diff_temp);

  DebogIJK::verifier("div_lambda_grad_T_volume", div_lambda_grad_T_volume_);

  if (formulation_favre_)
    {
      if (turbulent_diffusivity_)
        {
          statistiques().begin_count(cnt_kappa_t);
          if (flag_kappa_vectorial_)
            {
              calculer_flux_thermique_bord(temperature_, 0.,
                                           0, turbulent_kappa_vector_[2], flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           T_paroi_impose_kmin_, boundary_flux_kmin_, 0 /* boundary kmin */);
              calculer_flux_thermique_bord(temperature_, 0.,
                                           0, turbulent_kappa_vector_[2], flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           T_paroi_impose_kmax_, boundary_flux_kmax_, 1 /* boundary kmax */);

              temperature_.echange_espace_virtuel(1);
              turbulent_kappa_vector_[0].echange_espace_virtuel(1);
              turbulent_kappa_vector_[1].echange_espace_virtuel(1);
              turbulent_kappa_vector_[2].echange_espace_virtuel(1);

              calculer_diffusion_scalar(temperature_,
                                        turbulent_kappa_vector_[0],
                                        turbulent_kappa_vector_[1],
                                        turbulent_kappa_vector_[2],
                                        div_lambda_grad_T_volume_,
                                        boundary_flux_kmin_, boundary_flux_kmax_);
            }
          else
            {
              calculer_flux_thermique_bord(temperature_, 0.,
                                           0, turbulent_kappa_, flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           T_paroi_impose_kmin_, boundary_flux_kmin_, 0 /* boundary kmin */);
              calculer_flux_thermique_bord(temperature_, 0.,
                                           0, turbulent_kappa_, flag_kappa_anisotropic_,
                                           0, structural_uscalar_vector_[2],
                                           T_paroi_impose_kmax_, boundary_flux_kmax_, 1 /* boundary kmax */);

              temperature_.echange_espace_virtuel(1);
              turbulent_kappa_.echange_espace_virtuel(1);

              calculer_diffusion_scalar(temperature_,
                                        turbulent_kappa_,
                                        turbulent_kappa_,
                                        turbulent_kappa_,
                                        div_lambda_grad_T_volume_,
                                        boundary_flux_kmin_, boundary_flux_kmax_);
            }
          statistiques().end_count(cnt_kappa_t);

          DebogIJK::verifier("temp", temperature_);
          DebogIJK::verifier("div_lambda_grad_T_volume_kappa", div_lambda_grad_T_volume_);
        }

      if (structural_uscalar_)
        {
          statistiques().begin_count(cnt_kappa_struct);

          calculer_flux_thermique_bord(temperature_, 0.,
                                       0, turbulent_kappa_, flag_kappa_anisotropic_,
                                       0, structural_uscalar_vector_[2],
                                       T_paroi_impose_kmin_, boundary_flux_kmin_, 0 /* boundary kmin */);
          calculer_flux_thermique_bord(temperature_, 0.,
                                       0, turbulent_kappa_, flag_kappa_anisotropic_,
                                       0, structural_uscalar_vector_[2],
                                       T_paroi_impose_kmax_, boundary_flux_kmax_, 1 /* boundary kmax */);

          temperature_.echange_espace_virtuel(1);
          structural_uscalar_vector_[0].echange_espace_virtuel(1);
          structural_uscalar_vector_[1].echange_espace_virtuel(1);
          structural_uscalar_vector_[2].echange_espace_virtuel(1);

          DebogIJK::verifier("temp", temperature_);
          DebogIJK::verifier("structural_uscalar_x", structural_uscalar_vector_[0]);
          DebogIJK::verifier("structural_uscalar_y", structural_uscalar_vector_[1]);
          DebogIJK::verifier("structural_uscalar_z", structural_uscalar_vector_[2]);

          operateur_diffusion_temperature_structural_.ajouter(temperature_,
                                                              structural_uscalar_vector_[0],
                                                              structural_uscalar_vector_[1],
                                                              structural_uscalar_vector_[2],
                                                              div_lambda_grad_T_volume_,
                                                              boundary_flux_kmin_, boundary_flux_kmax_);

          statistiques().end_count(cnt_kappa_struct);

          DebogIJK::verifier("div_lambda_grad_T_volume_kappa_struct", div_lambda_grad_T_volume_);
        }
    }

  // Pour calculer la divergence, il faut les valeurs sur les faces de "droite" des elements, qui sont
  // des faces virtuelles a l'extremite du domaine, donc mettre a jour les faces de droite
  // A ete modifie lors de l'ajout de la source !!!!!
  velocity_[0].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_I*/
  velocity_[1].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_J*/
  velocity_[2].echange_espace_virtuel(1); /*, IJK_Field_double::EXCHANGE_GET_AT_RIGHT_K*/

  statistiques().begin_count(cnt_prepare_rhs);

  compute_divergence_times_constant(velocity_[0], velocity_[1], velocity_[2],
                                    - 1. / fractionnal_timestep, pressure_rhs_);
  DebogIJK::verifier("-div(u)/dt", pressure_rhs_);
  {
// Modif Martin 28-05-2019
    const int ni = rho_.ni();
    const int nj = rho_.nj();
    const int nk = rho_.nk();
    const double R_divise_par_Cp_Pth = constante_specifique_gaz_ / (Cp_gaz_ * P_th_final);
    const IJK_Grid_Geometry& geom = rho_.get_splitting().get_grid_geometry();
// On a besoin de vx et rho pour la pondÃ©ration du terme source par la vitesse
    IJK_Field_double& vx = velocity_[0];
    IJK_Field_double& rho = rho_;
    const double delta_x = geom.get_constant_delta(0);
    const double delta_y = geom.get_constant_delta(1);
    for (int k = 0; k < nk; k++)
      {
        const double delta_z = delta_z_local_[k];
        const double volume_maille = delta_x * delta_y * delta_z;
        const double facteur = d_Pth_divise_par_gammamoins1 * volume_maille;
// S il n y a pas d'ecoulement dans le canal, la source n'est pas ponderee par la vitesse
        if (debit_actuel_ < 0.0000001)
          {
            for (int j = 0; j < nj; j++)
              {
                for (int i = 0; i < ni; i++)
                  {
                    pressure_rhs_(i,j,k) += R_divise_par_Cp_Pth * (div_lambda_grad_T_volume_(i,j,k) - facteur - puit_ * volume_maille) / fractionnal_timestep;
                  }
              }
          }
        else
          {
            for (int j = 0; j < nj; j++)
              {
                for (int i = 0; i < ni; i++)
                  {
                    pressure_rhs_(i,j,k) += R_divise_par_Cp_Pth * (div_lambda_grad_T_volume_(i,j,k) - facteur - puit_ * rho(i,j,k) * vx(i,j,k) * Ly_tot_ * Lz_tot_ / debit_actuel_ * volume_maille) / fractionnal_timestep;
                  }
              }
          }
      }
  }
  DebogIJK::verifier("pressure_rhs", pressure_rhs_);
  statistiques().end_count(cnt_prepare_rhs);

  // Appel au solveur multigrille:
  if (!disable_solveur_poisson_)
    {
      statistiques().begin_count(cnt_resoudre_systeme_poisson);
      poisson_solver_.set_rho(rho_);
      poisson_solver_.resoudre_systeme_IJK(pressure_rhs_, pressure_);
      statistiques().end_count(cnt_resoudre_systeme_poisson);
      DebogIJK::verifier("pressure", pressure_);

// pressure gradient requires the "left" value in all directions:
      pressure_.echange_espace_virtuel(1 /*, IJK_Field_double::EXCHANGE_GET_AT_LEFT_IJK*/);
      statistiques().begin_count(cnt_ajouter_grad_p);
      add_gradient_times_constant_over_rho(pressure_, rho_, -fractionnal_timestep,
                                           velocity_[0], velocity_[1], velocity_[2]);
      statistiques().end_count(cnt_ajouter_grad_p);
    }

  // DD,2016-01-20: decalage de la pression en prenant comme reference la moyenne de la pression sur le volume.
  double reference_pression = 0;
  fixer_reference_pression(reference_pression);
  Cout << "reference_pression= " << reference_pression << finl;
  translation_pression(reference_pression);

  // Mise a jour des stats en prevision du prochain rk_step et du post_traitement. !!! (positionement optimum !!)
  rho_.echange_espace_virtuel(2);
  velocity_[0].echange_espace_virtuel(2);
  velocity_[1].echange_espace_virtuel(2);
  velocity_[2].echange_espace_virtuel(2);
}

void DNS_QC_double::fixer_reference_pression(double& reference_pression)
{
  const int ni = pressure_.ni();
  const int nj = pressure_.nj();
  const int nk = pressure_.nk();
  const double dx = pressure_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = pressure_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);

  double integrale_pression = 0.;
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              integrale_pression += pressure_(i,j,k) * dx * dy * delta_z_local_[k];
            }
        }
    }
  ArrOfDouble tmp(1);
  tmp[0] = integrale_pression;
  mp_sum_for_each_item(tmp);

  reference_pression = tmp[0] / (Lx_tot_*Ly_tot_*Lz_tot_);
}
void DNS_QC_double::translation_pression(const double& reference_pression)
{
  const int ni = pressure_.ni();
  const int nj = pressure_.nj();
  const int nk = pressure_.nk();
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              pressure_(i,j,k) -= reference_pression;
            }
        }
    }
}

void DNS_QC_double::compute_rho_bulk(double& rho_bulk)
{
  const int ni = rho_.ni();
  const int nj = rho_.nj();
  const int nk = rho_.nk();
  const double dx  = rho_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy  = rho_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);

  double integral_rho = 0.;
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              integral_rho += rho_(i,j,k) * dx * dy * delta_z_local_[k];
            }
        }
    }
  ArrOfDouble tmp(1);
  tmp[0] = integral_rho;
  mp_sum_for_each_item(tmp);

  rho_bulk = tmp[0] / (Lx_tot_*Ly_tot_*Lz_tot_);
}


void DNS_QC_double::compute_t_bulk(double& t_bulk)
{
  const int ni = temperature_.ni();
  const int nj = temperature_.nj();
  const int nk = temperature_.nk();
  const double dx = temperature_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = temperature_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);

  double integral_t = 0.;
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              integral_t += temperature_(i,j,k) * dx * dy * delta_z_local_[k];
            }
        }
    }
  ArrOfDouble tmp(1);
  tmp[0] = integral_t;
  mp_sum_for_each_item(tmp);

  t_bulk = tmp[0] / (Lx_tot_*Ly_tot_*Lz_tot_);
}

void DNS_QC_double::compute_rho_t_bulk(double& rho_bulk, double& t_bulk)
{
  const int ni = rho_.ni();
  const int nj = rho_.nj();
  const int nk = rho_.nk();
  const double dx = rho_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_I);
  const double dy = rho_.get_splitting().get_grid_geometry().get_constant_delta(DIRECTION_J);

  double integral_t   = 0.;
  double integral_rho = 0.;
  for (int k = 0; k < nk; k++)
    {
      for (int j = 0; j < nj; j++)
        {
          for (int i = 0; i < ni; i++)
            {
              integral_t   += temperature_(i,j,k) * dx * dy * delta_z_local_[k];
              integral_rho += rho_(i,j,k)         * dx * dy * delta_z_local_[k];
            }
        }
    }
  ArrOfDouble tmp(1);
  tmp[0] = integral_t;
  mp_sum_for_each_item(tmp);

  t_bulk    = tmp[0] / (Lx_tot_*Ly_tot_*Lz_tot_);

  ArrOfDouble tmp2(1);
  tmp2[0] = integral_rho;
  mp_sum_for_each_item(tmp2);

  rho_bulk = tmp2[0] / (Lx_tot_*Ly_tot_*Lz_tot_);
}
