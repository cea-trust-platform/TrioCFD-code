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
// File:        Calcul_Production_K_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Convection_Diffusion_Concentration.h>
#include <Calcul_Production_K_VEF.h>
#include <TRUSTTab.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <TRUSTVect.h>
#include <TRUSTTrav.h>
#include <Champ_P1NC.h>
#include <Periodique.h>
#include <Champ_Uniforme.h>
#include <Champ_Don_base.h>

/*! @brief Compute the production term for the turbulent kinetic energy.
 *
 * The total production term writes
 * \f$ P = -\frac{2}{3}k\text{div}(U) + 2\nu_t S_{ij}S_{ij}
 * Being in a incompressible flow, the first part is not computed.
 *
 * Like the TKE and epsilon, P is discretised on face centers.
 *
 * @param[in] const Domaine_VEF& domaine_VEF,
 * @param[in] const Domaine_Cl_VEF& zcl_VEF,
 * @param[in] DoubleTab& prodK,
 * @param[in] const DoubleTab& K_eps,
 * @param[in] const DoubleTab& vit,
 * @param[in] const DoubleTab& visco_turb,
 * @param[in] const int& interpol_visco,
 * @param[in] const double& limiteur
 * @param[out]
 * @return DoubleTab& prodK
 */
DoubleTab& Calcul_Production_K_VEF::calculer_terme_production_K(
  const Domaine_VEF& domaine_VEF,
  const Domaine_Cl_VEF& zcl_VEF,
  DoubleTab& prodK,
  const DoubleTab& K_eps,
  const DoubleTab& vit,
  const DoubleTab& visco_turb,
  const int& interpol_visco,
  const double& limiteur) const
{
  prodK = 0;

  // Compute the velocity gradient
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const int dimension = Objet_U::dimension;
  DoubleTab gradient_elem(nb_elem_tot, dimension, dimension);
  gradient_elem = 0.;
  Champ_P1NC::calcul_gradient(vit, gradient_elem, zcl_VEF);

  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleVect& volumes = domaine_VEF.volumes();

  // Boundary conditions
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int first_boundary_face = le_bord.num_premiere_face();
      const int last_boundary_face = first_boundary_face + le_bord.nb_faces();

      if (sub_type(Periodique, la_cl.valeur()))
        loop_for_internal_or_periodic_faces(prodK, gradient_elem, visco_turb, volumes,
                                            face_voisins, first_boundary_face, last_boundary_face,
                                            interpol_visco, limiteur);
      else
        loop_for_non_periodic_boundaries(prodK, gradient_elem, visco_turb, volumes,
                                         face_voisins, first_boundary_face, last_boundary_face,
                                         interpol_visco, limiteur);
    }

  // Internal faces
  const int premiere_face_int = domaine_VEF.premiere_face_int();
  const int nb_faces_ = domaine_VEF.nb_faces();
  loop_for_internal_or_periodic_faces(prodK, gradient_elem, visco_turb, volumes,
                                      face_voisins, premiere_face_int, nb_faces_,
                                      interpol_visco, limiteur);
  return prodK;
}

/*! @brief Compute production term on non periodic boundary faces
 *
 * Using the velocity gradient, the loop computes the production term written
 * \f$P = 2 \nu_t S_{ij}S_{ij}\f$
 * with \f$S_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i})\f$.
 *
 * @param[in] DoubleTab& prodK
 * @param[in] const DoubleTab& gradient_elem,
 * @param[in] const DoubleTab& visco_turb,
 * @param[in] const DoubleVect& volumes,
 * @param[in] const IntTab& face_voisins,
 * @param[in] const int nfaceinit,
 * @param[in] const int nfaceend,
 * @param[in] const int interpol_visco,
 * @param[in] const double limiteur
 * @return
 */
void Calcul_Production_K_VEF::loop_for_non_periodic_boundaries(DoubleTab& prodK,
                                                               const DoubleTab& gradient_elem,
                                                               const DoubleTab& visco_turb,
                                                               const DoubleVect& volumes,
                                                               const IntTab& face_voisins,
                                                               const int nfaceinit,
                                                               const int nfaceend,
                                                               const int interpol_visco,
                                                               const double limiteur
                                                              ) const
{
  for (int fac = nfaceinit; fac < nfaceend; fac++)
    {
      const int poly1 = face_voisins(fac, 0);
      const double visco_face = visco_turb(poly1);

      const double du_dx = gradient_elem(poly1, 0, 0);
      const double du_dy = gradient_elem(poly1, 0, 1);
      const double dv_dx = gradient_elem(poly1, 1, 0);
      const double dv_dy = gradient_elem(poly1, 1, 1);

      // Determination du terme de production
      prodK(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy)
                    + ((du_dy + dv_dx)*(du_dy + dv_dx)))*visco_face;

      if (Objet_U::dimension == 3)
        {
          const double du_dz = gradient_elem(poly1, 0, 2);
          const double dv_dz = gradient_elem(poly1, 1, 2);
          const double dw_dx = gradient_elem(poly1, 2, 0);
          const double dw_dy = gradient_elem(poly1, 2, 1);
          const double dw_dz = gradient_elem(poly1, 2, 2);

          // Determination du terme de production

          prodK(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz)
                        + ((du_dy + dv_dx) * (du_dy + dv_dx)
                           + (du_dz + dw_dx) * (du_dz + dw_dx)
                           + (dw_dy + dv_dz) * (dw_dy + dv_dz)))*visco_face;
        }
    }
}

/*! @brief Compute production term on internal and periodic boundary faces
 *
 * Using the velocity gradient, the loop computes the production term written
 * \f$P = 2 \nu_t S_{ij}S_{ij}\f$
 * with \f$S_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i})\f$.
 *
 * @param[in] DoubleTab& prodK
 * @param[in] const DoubleTab& gradient_elem,
 * @param[in] const DoubleTab& visco_turb,
 * @param[in] const DoubleVect& volumes,
 * @param[in] const IntTab& face_voisins,
 * @param[in] const int nfaceinit,
 * @param[in] const int nfaceend,
 * @param[in] const int interpol_visco,
 * @param[in] const double limiteur
 * @return
 */
void Calcul_Production_K_VEF::loop_for_internal_or_periodic_faces(DoubleTab& prodK,
                                                                  const DoubleTab& gradient_elem,
                                                                  const DoubleTab& visco_turb,
                                                                  const DoubleVect& volumes,
                                                                  const IntTab& face_voisins,
                                                                  const int nfaceinit,
                                                                  const int nfaceend,
                                                                  const int interpol_visco,
                                                                  const double limiteur
                                                                 ) const
{
  for (int fac = nfaceinit; fac < nfaceend; fac++)
    {
      const int poly1 = face_voisins(fac, 0);
      const int poly2 = face_voisins(fac, 1);
      const double a = volumes(poly1)/(volumes(poly1) + volumes(poly2));
      const double b = volumes(poly2)/(volumes(poly1) + volumes(poly2));

      const double visco_face = get_turbulent_viscosity(visco_turb, volumes,
                                                        interpol_visco, poly1, poly2,
                                                        limiteur);

      const double du_dx = a*gradient_elem(poly1, 0, 0) + b*gradient_elem(poly2, 0, 0);
      const double du_dy = a*gradient_elem(poly1, 0, 1) + b*gradient_elem(poly2, 0, 1);
      const double dv_dx = a*gradient_elem(poly1, 1, 0) + b*gradient_elem(poly2, 1, 0);
      const double dv_dy = a*gradient_elem(poly1, 1, 1) + b*gradient_elem(poly2, 1, 1);

      // Determination du terme de production

      prodK(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy)
                    + ((du_dy + dv_dx)*(du_dy + dv_dx)))*visco_face;

      if (Objet_U::dimension == 3)
        {
          const double du_dz = a*gradient_elem(poly1, 0, 2) + b*gradient_elem(poly2, 0, 2);
          const double dv_dz = a*gradient_elem(poly1, 1, 2) + b*gradient_elem(poly2, 1, 2);
          const double dw_dx = a*gradient_elem(poly1, 2, 0) + b*gradient_elem(poly2, 2, 0);
          const double dw_dy = a*gradient_elem(poly1, 2, 1) + b*gradient_elem(poly2, 2, 1);
          const double dw_dz = a*gradient_elem(poly1, 2, 2) + b*gradient_elem(poly2, 2, 2);

          // Determination du terme de production
          prodK(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz)
                        + ((du_dy + dv_dx)*(du_dy + dv_dx)
                           + (du_dz + dw_dx)*(du_dz + dw_dx)
                           + (dw_dy + dv_dz)*(dw_dy + dv_dz)))*visco_face;
        }
    }

}

DoubleTab& Calcul_Production_K_VEF::
calculer_terme_production_K_BiK(const Domaine_VEF& domaine_VEF,const Domaine_Cl_VEF& zcl_VEF,
                                DoubleTab& P,const DoubleTab& K,const DoubleTab& eps,
                                const DoubleTab& vit,const DoubleTab& visco_turb, const int& interpol_visco, const double& limiteur) const
{
  // P est discretise comme K et Eps i.e au centre des faces
  //
  // P(elem) = -(2/3)*k(i)*div_U(i) + nu_t(i) * F(u,v,w)
  //
  //                          2          2          2
  //    avec F(u,v,w) = 2[(du/dx)  + (dv/dy)  + (dw/dz) ] +
  //
  //                               2               2               2
  //                  (du/dy+dv/dx) + (du/dz+dw/dx) + (dw/dy+dv/dz)
  //
  // Rqs: On se place dans le cadre incompressible donc on neglige
  //      le terme (2/3)*k(i)*div_U(i)

  P= 0;

  // Calcul de F(u,v,w):
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleVect& volumes = domaine_VEF.volumes();
  const int dimension = Objet_U::dimension;

  DoubleTab gradient_elem(nb_elem_tot, dimension, dimension);
  gradient_elem = 0.;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //                        <
  // calcul des gradients;  < [ Ujp*np/vol(j) ]
  //                         j
  ////////////////////////////////////////////////////////////////////////////////////////////////

  Champ_P1NC::calcul_gradient(vit, gradient_elem, zcl_VEF);

  // Calcul des du/dx dv/dy et des derivees croisees sur les faces de chaque elements dans le cas 2D

  // Boucle sur les bords pour traiter les faces de bord
  // en distinguant le cas periodicite
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int poly1 = face_voisins(fac, 0);
              const int poly2 = face_voisins(fac, 1);
              const double a = volumes(poly1)/(volumes(poly1) + volumes(poly2));
              const double b = volumes(poly2)/(volumes(poly1) + volumes(poly2));

              const double visco_face = get_turbulent_viscosity(visco_turb, volumes,
                                                                interpol_visco, poly1, poly2,
                                                                limiteur);

              const double du_dx = a*gradient_elem(poly1, 0, 0) + b*gradient_elem(poly2, 0, 0);
              const double du_dy = a*gradient_elem(poly1, 0, 1) + b*gradient_elem(poly2, 0, 1);
              const double dv_dx = a*gradient_elem(poly1, 1, 0) + b*gradient_elem(poly2, 1, 0);
              const double dv_dy = a*gradient_elem(poly1, 1, 1) + b*gradient_elem(poly2, 1, 1);

              P(fac) = visco_face*( 2*(du_dx*du_dx + dv_dy*dv_dy)
                                    + ((du_dy+dv_dx)*(du_dy+dv_dx) ) );

              if (dimension == 3)
                {
                  const double du_dz = a*gradient_elem(poly1, 0, 2) + b*gradient_elem(poly2, 0, 2);
                  const double dv_dz = a*gradient_elem(poly1, 1, 2) + b*gradient_elem(poly2, 1, 2);
                  const double dw_dx = a*gradient_elem(poly1, 2, 0) + b*gradient_elem(poly2, 2, 0);
                  const double dw_dy = a*gradient_elem(poly1, 2, 1) + b*gradient_elem(poly2, 2, 1);
                  const double dw_dz = a*gradient_elem(poly1, 2, 2) + b*gradient_elem(poly2, 2, 2);

                  // Determination du terme de production
                  P(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz)
                            + ((du_dy + dv_dx)*(du_dy + dv_dx)
                               + (du_dz + dw_dx)*(du_dz + dw_dx)
                               + (dw_dy + dv_dz)*(dw_dy + dv_dz)))*visco_face;
                }
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int poly1 = face_voisins(fac,0);
              const double visco_face = visco_turb(poly1);

              const double du_dx = gradient_elem(poly1, 0, 0);
              const double du_dy = gradient_elem(poly1, 0, 1);
              const double dv_dx = gradient_elem(poly1, 1, 0);
              const double dv_dy = gradient_elem(poly1, 1, 1);

              P(fac) = visco_face*( 2*(du_dx*du_dx + dv_dy*dv_dy) + ((du_dy+dv_dx)*(du_dy+dv_dx)));

              if (dimension == 3)
                {
                  const double du_dz = gradient_elem(poly1, 0, 2);
                  const double dv_dz = gradient_elem(poly1, 1, 2);
                  const double dw_dx = gradient_elem(poly1, 2, 0);
                  const double dw_dy = gradient_elem(poly1, 2, 1);
                  const double dw_dz = gradient_elem(poly1, 2, 2);

                  P(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz)
                            + ( (du_dy + dv_dx)*(du_dy + dv_dx)
                                + (du_dz + dw_dx)*(du_dz + dw_dx)
                                + (dw_dy + dv_dz)*(dw_dy + dv_dz)))*visco_face;

                }
            }
        }
    }

  // Traitement des faces internes
  const int premiere_face_int = domaine_VEF.premiere_face_int();
  const int nb_faces_ = domaine_VEF.nb_faces();

  for (int fac = premiere_face_int; fac < nb_faces_; fac++)
    {
      const int poly1 = face_voisins(fac, 0);
      const int poly2 = face_voisins(fac, 1);
      const double a = volumes(poly1)/(volumes(poly1) + volumes(poly2));
      const double b = volumes(poly2)/(volumes(poly1) + volumes(poly2));

      const double visco_face = get_turbulent_viscosity(visco_turb, volumes,
                                                        interpol_visco, poly1, poly2,
                                                        limiteur);

      const double du_dx = a*gradient_elem(poly1, 0, 0) + b*gradient_elem(poly2, 0, 0);
      const double du_dy = a*gradient_elem(poly1, 0, 1) + b*gradient_elem(poly2, 0, 1);
      const double dv_dx = a*gradient_elem(poly1, 1, 0) + b*gradient_elem(poly2, 1, 0);
      const double dv_dy = a*gradient_elem(poly1, 1, 1) + b*gradient_elem(poly2, 1, 1);

      P(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy)
                + ((du_dy + dv_dx)*(du_dy + dv_dx)))*visco_face;
      if (dimension == 3)
        {
          const double du_dz = a*gradient_elem(poly1, 0, 2) + b*gradient_elem(poly2, 0, 2);
          const double dv_dz = a*gradient_elem(poly1, 1, 2) + b*gradient_elem(poly2, 1, 2);
          const double dw_dx = a*gradient_elem(poly1, 2, 0) + b*gradient_elem(poly2, 2, 0);
          const double dw_dy = a*gradient_elem(poly1, 2, 1) + b*gradient_elem(poly2, 2, 1);
          const double dw_dz = a*gradient_elem(poly1, 2, 2) + b*gradient_elem(poly2, 2, 2);

          P(fac) = (2*(du_dx*du_dx + dv_dy*dv_dy + dw_dz*dw_dz)
                    + ((du_dy + dv_dx)*(du_dy + dv_dx)
                       + (du_dz + dw_dx)*(du_dz + dw_dx)
                       + (dw_dy + dv_dz)*(dw_dy + dv_dz)))*visco_face;
        }
    }
  return P;
}



DoubleTab& Calcul_Production_K_VEF::calculer_terme_production_K_EASM(const Domaine_VEF& domaine_VEF,
                                                                     const Domaine_Cl_VEF& zcl_VEF,
                                                                     DoubleTab& P,
                                                                     const DoubleTab& K_eps,
                                                                     const DoubleTab& gradient_elem,
                                                                     const DoubleTab& visco_turb,
                                                                     const DoubleTab& Re,
                                                                     const int& interpol_visco,
                                                                     const double& limiteur) const
{
  //
  if (interpol_visco != 0 && interpol_visco != 1 && interpol_visco != 2)
    {
      Cerr << "interpol_visco must be equal to [0, 1, 2] and it is " << interpol_visco << ".\n";
      Process::exit();
    }

  // P : Production
  P = 0;

  // Compute velocity tensor on faces for
  const int nb_faces_ = domaine_VEF.nb_faces();
  const int dimension = Objet_U::dimension;

  DoubleTab gradient_face(nb_faces_, dimension, dimension);
  calcul_tenseur_face(gradient_face, gradient_elem, domaine_VEF, zcl_VEF);

  DoubleTab Re_face(nb_faces_, dimension, dimension);
  calcul_tenseur_face(Re_face, Re, domaine_VEF, zcl_VEF);

  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleVect& volumes = domaine_VEF.volumes();

  // Loop on boundary conditions
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int poly1 = face_voisins(fac, 0);
              const int poly2 = face_voisins(fac, 1);
              const double visco_face = get_turbulent_viscosity(visco_turb, volumes,
                                                                interpol_visco, poly1, poly2,
                                                                limiteur);
              compute_production_term_EASM(fac, visco_face, Re_face, gradient_face, P);
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int poly1 = face_voisins(fac, 0);
              const double visco_face = visco_turb(poly1);
              compute_production_term_EASM(fac, visco_face, Re_face, gradient_face, P);
            }
        }
    }

  // Loop on internal faces
  const int first_internal_face = domaine_VEF.premiere_face_int();
  for (int fac = first_internal_face; fac < nb_faces_; fac++)
    {
      const int poly1 = face_voisins(fac, 0);
      const int poly2 = face_voisins(fac, 1);
      const double visco_face = get_turbulent_viscosity(visco_turb, volumes,
                                                        interpol_visco, poly1, poly2,
                                                        limiteur);
      compute_production_term_EASM(fac, visco_face, Re_face, gradient_face, P);
    }
  return P;
}

void Calcul_Production_K_VEF::compute_production_term_EASM(const int face,
                                                           const double visco_face,
                                                           const DoubleTab& Re_face,
                                                           const DoubleTab& gradient_face,
                                                           DoubleTab& P) const
{
  for (int i = 0; i < Objet_U::dimension; i++)
    for (int j = 0; j < Objet_U::dimension; j++)
      P(face) += Re_face(face, i, j) * gradient_face(face, i, j);
  P(face) *= visco_face;
}

/*! @brief Get the turbulent viscosity depending on the interpolation choice.
 *
 * if type_interpo == 0: arithmetic interpolation (instable)
 * if type_interpo == 1: harmonic mean (used for the realisable k-epsilon)
 * if type_interpo == 2: weighted harmonic mean (used for the realisable k-epsilon)
 *
 * @param[in] const DoubleTab& visco_turb
 * @param[in] const DoubleVect& volumes
 * @param[in] const int type_interpo
 * @param[in] const int poly1
 * @param[in] const int poly2
 * @param[in] const double limiteur
 * @param[out]
 * @return double
 */
double Calcul_Production_K_VEF::get_turbulent_viscosity(const DoubleTab& visco_turb,
                                                        const DoubleVect& volumes,
                                                        const int type_interpo,
                                                        const int poly1,
                                                        const int poly2,
                                                        const double limiteur) const
{
  switch(type_interpo)
    {
    case 0: // cas initial, toutes les interpolation arithmetiques ponderees sont fortement instables
      return 0.5*(visco_turb(poly1) + visco_turb(poly2));

    case 1: //Moyenne harmonique (uniquement utilisee dans le cas du keps realisable)
      if (visco_turb(poly1) > 1.e-10 && visco_turb(poly2) > 1.e-10)
        return limiteur/(1./visco_turb(poly1) + 1./visco_turb(poly2));
      else
        return limiteur*0.5*(visco_turb(poly1) + visco_turb(poly2));

    case 2: //Moyenne harmonique ponderee pour garantir la continuite du tenseur des contraintes a la face (uniquement utilisee dans le cas du keps realisable)
      if (visco_turb(poly1) > 1.e-10 && visco_turb(poly2) > 1.e-10)
        {
          const double a = volumes(poly1)/(volumes(poly1) + volumes(poly2));
          const double b = volumes(poly2)/(volumes(poly1) + volumes(poly2));
          return limiteur*(visco_turb(poly1)*visco_turb(poly2))/(b*visco_turb(poly1) + a*visco_turb(poly2));
        }
      else
        return limiteur*0.5*(visco_turb(poly1) + visco_turb(poly2));

    default:
      Cerr << "Interpolation type should be 0, 1 or 2, not " << type_interpo << "\n";
      Process::exit();
      return 1;
    }
  return 0;
}

// Utilisée dans le bas Reynolds anisotherme
DoubleTab& Calcul_Production_K_VEF::calculer_terme_destruction_K_gen(
  const Domaine_VEF& domaine_VEF,
  const Domaine_Cl_VEF& zcl_VEF,
  DoubleTab& G,
  const DoubleTab& inconnue, // scalar, concentration, temperature, etc.
  const DoubleTab& alpha_turb,
  const Champ_Don_base& ch_beta,
  const DoubleVect& gravite,
  int nb_consti) const
{
  // G est discretise comme K et Eps i.e au centre des faces
  // G(face) = beta alpha_t(face) G . gradT(face)

  if (nb_consti < 0)
    Process::exit("nb_consti is negative");

  const int nb_compo = sub_type(Champ_Uniforme, ch_beta) ? 0 : ch_beta.nb_comp();

  if (nb_compo < 0)
    Process::exit("nb_compo is negative");

  const int dimension = Objet_U::dimension;

  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleVect& volumes = domaine_VEF.volumes();
  const DoubleTab& tab_beta = ch_beta.valeurs();

  G = 0;

  if (nb_consti == 0 || nb_consti == 1)
    {
      // Compute the gradient of the equation unknwon
      const int nb_elem_tot = domaine_VEF.nb_elem_tot();
      DoubleTrav gradient_elem(nb_elem_tot, dimension);
      gradient_elem = 0;
      Champ_P1NC::calcul_gradient(inconnue, gradient_elem, zcl_VEF);

      // Compute u_theta
      const int nb_faces_tot = domaine_VEF.nb_faces_tot();
      DoubleTrav u_theta(nb_faces_tot, dimension);
      u_theta = 0;

      if (nb_compo == 0)
        compute_utheta_nbConsti_le_1_nbCompo_eq_0(domaine_VEF, zcl_VEF, face_voisins,
                                                  volumes, tab_beta, alpha_turb, gradient_elem,
                                                  u_theta);
      else if (nb_compo == 1)
        compute_utheta_nbConsti_le_1_nbCompo_eq_1(domaine_VEF, zcl_VEF, face_voisins,
                                                  volumes, tab_beta, alpha_turb, gradient_elem,
                                                  u_theta);

      else if (nb_compo > 1)
        compute_utheta_nbConsti_le_1_nbCompo_gt_1(domaine_VEF, zcl_VEF, face_voisins,
                                                  volumes, tab_beta, alpha_turb, gradient_elem,
                                                  u_theta);


      // Calcul de gravite . u_teta
      const int nb_faces_ = domaine_VEF.nb_faces();
      for (int fac = 0; fac < nb_faces_; fac++)
        {
          G[fac] = gravite(0)*u_theta(fac, 0) + gravite(1)*u_theta(fac, 1);
          if (dimension == 3)
            G[fac] += gravite(2)*u_theta(fac, 2);
        }
    }
  else if (nb_consti > 1)
    {
      // Compute the gradient of the equation unknwon
      const int nb_elem_tot = domaine_VEF.nb_elem_tot();
      DoubleTrav gradient_elem(nb_elem_tot, nb_consti, dimension);
      gradient_elem = 0;
      Champ_P1NC::calcul_gradient(inconnue, gradient_elem, zcl_VEF);

      // Compute u_theta
      const int nb_faces_tot = domaine_VEF.nb_faces_tot();
      DoubleTrav u_theta(nb_faces_tot, nb_consti, dimension);
      u_theta = 0;

      if (nb_compo == 0)
        compute_utheta_nbConsti_gt_1_nbCompo_eq_0(domaine_VEF, zcl_VEF, face_voisins,
                                                  volumes, tab_beta, alpha_turb, gradient_elem,
                                                  nb_consti, u_theta);
      else if (nb_compo == 1)
        compute_utheta_nbConsti_gt_1_nbCompo_eq_1(domaine_VEF, zcl_VEF, face_voisins,
                                                  volumes, tab_beta, alpha_turb, gradient_elem,
                                                  nb_consti, u_theta);
      else if (nb_compo > 1)
        compute_utheta_nbConsti_gt_1_nbCompo_gt_1(domaine_VEF, zcl_VEF, face_voisins,
                                                  volumes, tab_beta, alpha_turb, gradient_elem,
                                                  nb_consti, u_theta);

      // Calcul de gravite . u_teta
      const int nb_faces_ = domaine_VEF.nb_faces();
      for (int fac = 0; fac < nb_faces_; fac++)
        for (int k = 0; k < nb_consti; k++)
          {
            G[fac] = gravite(0)*u_theta(fac, k, 0) + gravite(1)*u_theta(fac, k, 1);
            if (dimension == 3)
              G[fac] += gravite(2)*u_theta(fac, k, 2);
          }
    }

  return G;
}

void Calcul_Production_K_VEF::compute_utheta_nbConsti_le_1_nbCompo_eq_0(const Domaine_VEF& domaine_VEF,
                                                                        const Domaine_Cl_VEF& zcl_VEF,
                                                                        const IntTab& face_voisins,
                                                                        const DoubleVect& volumes,
                                                                        const DoubleTab& tab_beta,
                                                                        const DoubleTab& alpha_turb,
                                                                        const DoubleTrav& gradient_elem,
                                                                        DoubleTrav& u_theta) const
{
  // we treat the boundaries
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl   = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb ; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              const int elem2 = face_voisins(fac, 1);
              const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
              const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
              for (int i = 0; i < Objet_U::dimension; i++)
                u_theta(fac, i) = a*tab_beta(0, 0)*alpha_turb(elem1)*gradient_elem(elem1, i)
                                  + b*tab_beta(0, 0)*alpha_turb(elem2)*gradient_elem(elem2, i);
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              for (int i = 0; i < Objet_U::dimension; i++)
                u_theta(fac, i) = tab_beta(0, 0)*alpha_turb(elem1)*gradient_elem(elem1, i);
            }
        }
    }

  // we treat the internal faces
  const int nb_faces_tot = domaine_VEF.nb_faces_tot();
  for (int fac = 0; fac < nb_faces_tot; fac++)
    {
      const int elem1 = face_voisins(fac, 0);
      const int elem2 = face_voisins(fac, 1);
      if ((elem1 >= 0) && (elem2 >= 0))
        {
          double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
          double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
          for (int i = 0; i < Objet_U::dimension; i++)
            u_theta(fac, i) = a*tab_beta(0, 0)*alpha_turb(elem1)*gradient_elem(elem1, i)
                              + b*tab_beta(0, 0)*alpha_turb(elem2)*gradient_elem(elem2, i);
        }
    }
}

void Calcul_Production_K_VEF::compute_utheta_nbConsti_le_1_nbCompo_eq_1(const Domaine_VEF& domaine_VEF,
                                                                        const Domaine_Cl_VEF& zcl_VEF,
                                                                        const IntTab& face_voisins,
                                                                        const DoubleVect& volumes,
                                                                        const DoubleTab& tab_beta,
                                                                        const DoubleTab& alpha_turb,
                                                                        const DoubleTrav& gradient_elem,
                                                                        DoubleTrav& u_theta) const
{
  // we treat the boundaries
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              const int elem2 = face_voisins(fac, 1);
              const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
              const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
              for (int i = 0; i < Objet_U::dimension; i++)
                u_theta(fac, i) = a*tab_beta(elem1)*alpha_turb(elem1)*gradient_elem(elem1, i)
                                  + b*tab_beta(elem2)*alpha_turb(elem2)*gradient_elem(elem2, i);
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              for (int i = 0; i < Objet_U::dimension; i++)
                u_theta(fac, i) = tab_beta(elem1)*alpha_turb(elem1)*gradient_elem(elem1, i);
            }
        }
    }
  // we treat the internal faces
  const int nb_faces_tot = domaine_VEF.nb_faces_tot();
  for (int fac = 0; fac < nb_faces_tot; fac++)
    {
      const int elem1 = face_voisins(fac, 0);
      const int elem2 = face_voisins(fac, 1);
      if ((elem1 >= 0) && (elem2 >= 0))
        {
          const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
          const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
          for (int i = 0; i < Objet_U::dimension; i++)
            u_theta(fac, i) = a*tab_beta(elem1)*alpha_turb(elem1)*gradient_elem(elem1, i)
                              + b*tab_beta(elem2)*alpha_turb(elem2)*gradient_elem(elem2, i);
        }
    }
}

void Calcul_Production_K_VEF::compute_utheta_nbConsti_le_1_nbCompo_gt_1(const Domaine_VEF& domaine_VEF,
                                                                        const Domaine_Cl_VEF& zcl_VEF,
                                                                        const IntTab& face_voisins,
                                                                        const DoubleVect& volumes,
                                                                        const DoubleTab& tab_beta,
                                                                        const DoubleTab& alpha_turb,
                                                                        const DoubleTrav& gradient_elem,
                                                                        DoubleTrav& u_theta) const
{
  // we treat the boundaries
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              const int elem2 = face_voisins(fac, 1);
              const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
              const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
              for (int i = 0; i < Objet_U::dimension; i++)
                u_theta(fac, i) = a*tab_beta(elem1, 0)*alpha_turb(elem1)*gradient_elem(elem1, i)
                                  + b*tab_beta(elem2, 0)*alpha_turb(elem2)*gradient_elem(elem2, i);

            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              for (int i = 0; i < Objet_U::dimension; i++)
                u_theta(fac, i) = tab_beta(elem1, 0)*alpha_turb(elem1)*gradient_elem(elem1, i);
            }
        }
    }
  // we treat the internal faces
  const int nb_faces_tot = domaine_VEF.nb_faces_tot();
  for (int fac = 0; fac < nb_faces_tot; fac++)
    {
      const int elem1 = face_voisins(fac, 0);
      const int elem2 = face_voisins(fac, 1);
      if ((elem1 >= 0) && (elem2 >= 0))
        {
          const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
          const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
          for (int i = 0; i < Objet_U::dimension; i++)
            u_theta(fac, i) = a*tab_beta(elem1, 0)*alpha_turb(elem1)*gradient_elem(elem1, i)
                              + b*tab_beta(elem2, 0)*alpha_turb(elem2)*gradient_elem(elem2, i);
        }
    }
}

void Calcul_Production_K_VEF::compute_utheta_nbConsti_gt_1_nbCompo_eq_0(const Domaine_VEF& domaine_VEF,
                                                                        const Domaine_Cl_VEF& zcl_VEF,
                                                                        const IntTab& face_voisins,
                                                                        const DoubleVect& volumes,
                                                                        const DoubleTab& tab_beta,
                                                                        const DoubleTab& alpha_turb,
                                                                        const DoubleTrav& gradient_elem,
                                                                        const int nb_consti,
                                                                        DoubleTrav& u_theta) const
{
  // we treat the boundaries
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              const int elem2 = face_voisins(fac, 1);
              const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
              const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
              for (int i = 0; i < nb_consti; i++)
                for (int k = 0; k < Objet_U::dimension; k++)
                  u_theta(fac, i, k) = a*tab_beta(0, 0)*alpha_turb(elem1)*gradient_elem(elem1, i, k)
                                       + b*tab_beta(0, 0)*alpha_turb(elem2)*gradient_elem(elem2, i, k);
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              for (int i = 0; i < nb_consti; i++)
                for (int k = 0; k < Objet_U::dimension; k++)
                  u_theta(fac, i, k) = tab_beta(0, 0)*alpha_turb(elem1)*gradient_elem(elem1, i, k);
            }
        }
    }

  // we treat the internal faces
  const int nb_faces_tot = domaine_VEF.nb_faces_tot();
  for (int fac = 0; fac < nb_faces_tot; fac++)
    {
      const int elem1 = face_voisins(fac, 0);
      const int elem2 = face_voisins(fac, 1);
      if ((elem1 >= 0) && (elem2 >= 0))
        {
          const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
          const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
          for (int i = 0; i < nb_consti; i++)
            for (int k = 0; k < Objet_U::dimension; k++)
              u_theta(fac, i, k) = a*tab_beta(0, 0)*alpha_turb(elem1)*gradient_elem(elem1, i, k)
                                   + b*tab_beta(0, 0)*alpha_turb(elem2)*gradient_elem(elem2, i, k);
        }
    }
}

void Calcul_Production_K_VEF::compute_utheta_nbConsti_gt_1_nbCompo_eq_1(const Domaine_VEF& domaine_VEF,
                                                                        const Domaine_Cl_VEF& zcl_VEF,
                                                                        const IntTab& face_voisins,
                                                                        const DoubleVect& volumes,
                                                                        const DoubleTab& tab_beta,
                                                                        const DoubleTab& alpha_turb,
                                                                        const DoubleTrav& gradient_elem,
                                                                        const int nb_consti,
                                                                        DoubleTrav& u_theta) const
{
  // we treat the boundaries
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              const int elem2 = face_voisins(fac, 1);
              const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
              const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
              for (int i = 0; i < nb_consti; i++)
                for (int k = 0; k < Objet_U::dimension; k++)
                  u_theta(fac, i, k) = a*tab_beta(elem1)*alpha_turb(elem1)*gradient_elem(elem1, i, k)
                                       + b*tab_beta(elem2)*alpha_turb(elem2)*gradient_elem(elem2, i, k);
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              for (int i = 0; i < nb_consti; i++)
                for (int k = 0; k < Objet_U::dimension; k++)
                  u_theta(fac, i, k) = tab_beta(elem1)*alpha_turb(elem1)*gradient_elem(elem1, i, k);
            }
        }
    }
  // we treat the internal faces
  const int nb_faces_tot = domaine_VEF.nb_faces_tot();
  for (int fac = 0; fac < nb_faces_tot; fac++)
    {
      const int elem1 = face_voisins(fac, 0);
      const int elem2 = face_voisins(fac, 1);
      if ((elem1 >= 0) && (elem2 >= 0))
        {
          const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
          const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
          for (int i = 0; i < nb_consti; i++)
            for (int k = 0; k < Objet_U::dimension; k++)
              u_theta(fac, i, k) = a*tab_beta(elem1)*alpha_turb(elem1)*gradient_elem(elem1, i, k)
                                   + b*tab_beta(elem2)*alpha_turb(elem2)*gradient_elem(elem2, i, k);
        }
    }

}

void Calcul_Production_K_VEF::compute_utheta_nbConsti_gt_1_nbCompo_gt_1(const Domaine_VEF& domaine_VEF,
                                                                        const Domaine_Cl_VEF& zcl_VEF,
                                                                        const IntTab& face_voisins,
                                                                        const DoubleVect& volumes,
                                                                        const DoubleTab& tab_beta,
                                                                        const DoubleTab& alpha_turb,
                                                                        const DoubleTrav& gradient_elem,
                                                                        const int nb_consti,
                                                                        DoubleTrav& u_theta) const
{
  // we treat the boundaries
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl   = zcl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              const int elem2 = face_voisins(fac, 1);
              const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
              const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
              for (int i = 0; i < nb_consti; i++)
                for (int k = 0; k < Objet_U::dimension; k++)
                  u_theta(fac, i, k) = a*tab_beta(elem1, 0)*alpha_turb(elem1)*gradient_elem(elem1, i, k)
                                       + b*tab_beta(elem2, 0)*alpha_turb(elem2)*gradient_elem(elem2, i, k);
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int elem1 = face_voisins(fac, 0);
              for (int i = 0; i < nb_consti; i++)
                for (int k = 0; k < Objet_U::dimension; k++)
                  u_theta(fac, i, k) = tab_beta(elem1, 0)*alpha_turb(elem1)*gradient_elem(elem1, i, k);
            }
        }
    }
  // we treat the internal faces
  const int nb_faces_tot = domaine_VEF.nb_faces_tot();
  for (int fac = 0; fac < nb_faces_tot; fac++)
    {
      const int elem1 = face_voisins(fac, 0);
      const int elem2 = face_voisins(fac, 1);
      if ((elem1 >= 0) && (elem2 >= 0))
        {
          const double a = volumes(elem1)/(volumes(elem1) + volumes(elem2));
          const double b = volumes(elem2)/(volumes(elem1) + volumes(elem2));
          for (int i = 0; i < nb_consti; i++)
            for (int k = 0; k < Objet_U::dimension; k++)
              u_theta(fac, i, k) = a*tab_beta(elem1, 0)*alpha_turb(elem1)*gradient_elem(elem1, i, k)
                                   + b*tab_beta(elem2, 0)*alpha_turb(elem2)*gradient_elem(elem2, i, k);
        }
    }
}

// Calcul d'un tenseur aux faces a partir d'un tenseur aux elements
DoubleTab& Calcul_Production_K_VEF::calcul_tenseur_face(DoubleTab& Tenseur_face,
                                                        const DoubleTab& Tenseur_elem,
                                                        const Domaine_VEF& domaine_VEF,
                                                        const Domaine_Cl_VEF& domaine_Cl_VEF) const
{

  const int dimension = Objet_U::dimension;
  const IntTab& face_voisins = domaine_VEF.face_voisins();

  const Conds_lim& les_cl = domaine_Cl_VEF.les_conditions_limites();
  const int nb_cl = les_cl.size();
  const DoubleVect& volumes = domaine_VEF.volumes();

  for (int n_bord = 0; n_bord < nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int poly1 = face_voisins(fac, 0);
              const int poly2 = face_voisins(fac, 1);
              const double invVol = 1/(volumes(poly1) + volumes(poly2));
              const double a = volumes(poly1)*invVol;
              const double b = volumes(poly2)*invVol;
              for (int i = 0; i < dimension; i++)
                for (int j = 0; j < dimension; j++)
                  Tenseur_face(fac, i, j) = a*Tenseur_elem(poly1, i, j)
                                            + b*Tenseur_elem(poly2, i, j);
            }
        }
      else
        {
          for (int fac = ndeb; fac < nfin; fac++)
            {
              const int poly1 = face_voisins(fac, 0);
              for (int i = 0; i < dimension; i++)
                for (int j = 0; j < dimension; j++)
                  Tenseur_face(fac, i, j) = Tenseur_elem(poly1, i, j);
            }
        }
    }

  const int n0 = domaine_VEF.premiere_face_int();
  const int nb_faces = domaine_VEF.nb_faces();
  for (int fac = n0; fac < nb_faces; fac++)
    {
      const int poly1 = face_voisins(fac, 0);
      const int poly2 = face_voisins(fac, 1);
      const double invVol = 1/(volumes(poly1) + volumes(poly2));
      const double a = volumes(poly1)*invVol;
      const double b = volumes(poly2)*invVol;
      for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
          Tenseur_face(fac, i, j) = a*Tenseur_elem(poly1, i, j) + b*Tenseur_elem(poly2, i, j);
    }

  return Tenseur_face;
}
