/****************************************************************************
* Copyright (c) 2023, CEA
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

#ifndef Iterateur_VDF_Elem_FT_TCL_TPP_included
#define Iterateur_VDF_Elem_FT_TCL_TPP_included

#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Triple_Line_Model_FT_Disc.h>
#include <Probleme_FT_Disc_gen.h>

/*
 * XXX XXX XXX
 * Elie Saikali : NOTA BENE : fichier surcharger dans trio pour le FT, TCL model
 */

#define TCL_MODEL 1

template<class _TYPE_> template<typename Type_Double>
inline void Iterateur_VDF_Elem<_TYPE_>::fill_flux_tables_(const int face, const int ncomp, const double coeff, const Type_Double& flux, DoubleTab& resu) const
{
  DoubleTab& flux_bords = op_base->flux_bords();
  const int elem1 = elem(face, 0), elem2 = elem(face, 1);

#if TCL_MODEL
  const Equation_base& eq = op_base->equation();
  if (sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
    {
      Convection_Diffusion_Temperature_FT_Disc& eq_typee = ref_cast_non_const(Convection_Diffusion_Temperature_FT_Disc, eq);
      Transport_Interfaces_FT_Disc& eq_interface = eq_typee.eq_interface();
      const DoubleTab& indicatrice = eq_interface.get_update_indicatrice().valeurs();
      if (elem1 > -1)
        for (int k = 0; k < ncomp; k++)
          {
            resu(elem1, k) += coeff * flux[k];
            flux_bords(face, k) += coeff * flux[k];
            if ((indicatrice(elem1) != 0.) && (indicatrice(elem1) != 1.))
              {
                Cerr << "echange_externe_impose face-no= " << face << " elem1= " << elem1 << " flux= " << flux(k) << finl;
                // Energy is not resolved in elem0 (mixed-cell)
                eq_typee.mixed_elems().append_array(elem1);
                eq_typee.lost_fluxes().append_array(flux(k));
              }

          }
      if (elem2 > -1)
        for (int k = 0; k < ncomp; k++)
          {
            resu(elem2, k) -= coeff * flux[k];
            flux_bords(face, k) -= coeff * flux[k];

            if ((indicatrice(elem2) != 0.) && (indicatrice(elem2) != 1.))
              {
                Cerr << "echange_externe_impose face-no= " << face << " elem2= " << elem2 << " flux= " << flux(k) << finl;
                // Energy is not resolved in elem1 (mixed-cell)
                eq_typee.mixed_elems().append_array(elem2);
                eq_typee.lost_fluxes().append_array(-flux(k)); // see convention above, the flux should be taken with "-" for elems 1!
              }
          }
    }
  else
    {
      if (elem1 > -1)
        for (int k = 0; k < ncomp; k++)
          {
            resu(elem1, k) += coeff * flux[k];
            flux_bords(face, k) += coeff * flux[k];
          }
      if (elem2 > -1)
        for (int k = 0; k < ncomp; k++)
          {
            resu(elem2, k) -= coeff * flux[k];
            flux_bords(face, k) -= coeff * flux[k];
          }
    }

#else
  if (elem1 > -1)
    for (int k = 0; k < ncomp; k++)
      {
        resu(elem1, k) += coeff * flux[k];
        flux_bords(face, k) += coeff * flux[k];
      }
  if (elem2 > -1)
    for (int k = 0; k < ncomp; k++)
      {
        resu(elem2, k) -= coeff * flux[k];
        flux_bords(face, k) -= coeff * flux[k];
      }
#endif
}

template<class _TYPE_> template<typename Type_Double>
bool Iterateur_VDF_Elem<_TYPE_>::ajouter_blocs_bords_echange_ext_FT_TCL(const Echange_externe_impose& cl, const int ndeb, const int nfin, const int num_cl, const int N, const Front_VF& frontiere_dis,
                                                                        matrices_t mats, DoubleTab& resu, const tabs_t& semi_impl) const
{
  const DoubleTab& donnee = semi_impl.count(nom_ch_inco_) ? semi_impl.at(nom_ch_inco_) : le_champ_convecte_ou_inc->valeurs();
  Type_Double flux(N), aii(N), ajj(N), aef(N);
  int boundary_index = -1;
  if (le_dom.valeur().front_VF(num_cl).le_nom() == frontiere_dis.le_nom())
    boundary_index = num_cl;

  int e, Mv = le_ch_v.non_nul() ? le_ch_v->valeurs().line_size() : N;

#if TCL_MODEL
  const IntTab& face_voisins = le_dom->face_voisins();
  const Probleme_base& pb_gen = op_base->equation().probleme();
  const Equation_base& eq = op_base->equation();
  const Probleme_FT_Disc_gen *pbft = dynamic_cast<const Probleme_FT_Disc_gen*>(&pb_gen);
  const Triple_Line_Model_FT_Disc *tcl = pbft ? &pbft->tcl() : nullptr;
  const int k = 0; // component of the variable (which is a scalar)
#endif

  for (int face = ndeb; face < nfin; face++)
    {
      const int local_face = le_dom.valeur().front_VF(boundary_index).num_local_face(face);
      flux_evaluateur.flux_face(donnee, boundary_index, face, local_face, cl, ndeb, flux);
#if TCL_MODEL
      if (sub_type(Convection_Diffusion_Temperature_FT_Disc, eq) && tcl && tcl->is_activated())
        {
          // The model contribution has been computed before, and stored.
          // And then, it can be simply accessed in const mode here.
          // Cerr << "TCL is activated in flux iterator, we modify the contribution" << finl;
          // get the CL contributions:
          // HACK: dirty code, Where is my IT guy?
          // I would need to move these reference to tables before the loop "for" to avoid making them
          // so many times but I don't know how to proceed as "tcl" is only defined if sub_types are true.
          // Or is it? Maybe it works in any case if any pb can have an empty tcl model as an embedded object...
          const ArrOfInt& elems_with_CL_contrib = tcl->elems();
          const ArrOfInt& boundary_faces = tcl->boundary_faces();
          const ArrOfDouble& Q_from_CL = tcl->Q();
          //tcl.compute_TCL_fluxes_in_all_boundary_cells(elems_with_CL_contrib, mpoint_from_CL, Q_from_CL);
          const double sign = (face_voisins(face, 0) == -1) ? -1. : 1.;
          const int nb_contact_line_contribution = Q_from_CL.size_array();
          int nb_contrib = 0;
          for (int idx = 0; idx < nb_contact_line_contribution; idx++)
            {
              const int elemi = elems_with_CL_contrib[idx];
              const int facei = boundary_faces[idx];
              const int elem_bord_with_facei = face_voisins(facei, 0)+face_voisins(facei, 1) +1;
              if (facei == face)
                {
                  // The corresponding contribution should be assigned to the face
                  nb_contrib++;
                  const double TCL_wall_flux = Q_from_CL[idx];
                  // val should be : -rho*Cp * flux(W)
                  // probably because the whole energy equation is written with rhoCp somewhere...
                  // and the sign should be negative for incoming flux (towards the fluid) by convention.
                  const double val = sign*TCL_wall_flux;
                  Cerr << "GB face: " << face << " former-flux = " << flux[k];
                  if (nb_contrib == 1)
                    {
                      flux[k] = val; // To erase the standard calculation without the model and replace it
                      //             by part of the model
                    }
                  else
                    {
                      // If it's not the first contribution, it should be '+='
                      flux[k] += val;
                    }
                  Cerr << " new-flux = " << flux[k] << " [contrib #" << nb_contrib << "]" <<  finl;
                  if (elem_bord_with_facei != elemi)
                    {
                      // In this case, the flux is meso-region from a cell not adjacent to the wall.
                      // So what? Does it really matter?
                      // I believe that the fact that we take the interfacial flux from elemi and apply it at the wall
                      // face that is a boundary to elem_of_facei is not really an issue.
                    }
                }
            }
        }
#endif
      fill_flux_tables_(face, N, 1.0 /* coeff */, flux, resu);
    }

  Matrice_Morse *m_vit = (mats.count("vitesse") && is_convective_op()) ? mats.at("vitesse") : nullptr, *mat = (!is_pb_multiphase() && mats.count(nom_ch_inco_)) ? mats.at(nom_ch_inco_) : nullptr;
  VectorDeriv d_cc;
  fill_derivee_cc(mats, semi_impl, d_cc);

  //derivees : vitesse
  if (m_vit)
    {
      DoubleTab val_b = use_base_val_b_ ? le_champ_convecte_ou_inc->Champ_base::valeur_aux_bords() : le_champ_convecte_ou_inc->valeur_aux_bords();
      for (int face = ndeb; face < nfin; face++)
        {
          const int local_face = le_dom.valeur().front_VF(boundary_index).num_local_face(face);
          flux_evaluateur.coeffs_face_bloc_vitesse(donnee, val_b, boundary_index, face, local_face, cl, ndeb, aef);

          for (int i = 0; i < 2; i++)
            if ((e = elem(face, i)) >= 0)
              for (int n = 0, m = 0; n < N; n++, m += (Mv > 1)) (*m_vit)(N * e + n, Mv * face + m) += (i ? -1.0 : 1.0) * aef(n);
        }
    }

  //derivees : champ convecte
  if (mat || d_cc.size() > 0)
    for (int face = ndeb; face < nfin; face++)
      {
        const int local_face = le_dom.valeur().front_VF(boundary_index).num_local_face(face);
        flux_evaluateur.coeffs_face(boundary_index, face, local_face, ndeb, cl, aii, ajj);
        fill_coeffs_matrices(face, aii, ajj, mat, d_cc); // XXX : Attention Yannick pour d_cc c'est pas tout a fait comme avant ... N et M ...
      }

  return true; /* XXX ATTENTION : true dans trio ... */
}

#endif /* Iterateur_VDF_Elem_FT_TCL_TPP_included */
