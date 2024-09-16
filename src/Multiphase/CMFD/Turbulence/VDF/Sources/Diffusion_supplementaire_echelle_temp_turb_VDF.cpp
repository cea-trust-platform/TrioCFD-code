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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Diffusion_supplementaire_echelle_temp_turb_VDF.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/VDF
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Diffusion_supplementaire_echelle_temp_turb_VDF.h>

#include <Echelle_temporelle_turbulente.h>
#include <Op_Dift_Multiphase_VDF_Elem.h>
#include <Convection_Diffusion_std.h>
#include <Champ_Inc_P0_base.h>
#include <Domaine_Cl_VDF.h>
#include <Pb_Multiphase.h>
#include <Domaine_VDF.h>
#include <Dirichlet.h>

#include <cmath>
#include <vector>

Implemente_instanciable(Diffusion_supplementaire_echelle_temp_turb_VDF,"Diffusion_supplementaire_lin_echelle_temp_turb_VDF_P0_VDF|Diffusion_supplementaire_echelle_temp_turb_VDF_P0_VDF", Source_Diffusion_supplementaire_echelle_temp_turb);

Sortie& Diffusion_supplementaire_echelle_temp_turb_VDF::printOn(Sortie& os) const {  return Source_Diffusion_supplementaire_echelle_temp_turb::printOn(os);}

Entree& Diffusion_supplementaire_echelle_temp_turb_VDF::readOn(Entree& is) { return Source_Diffusion_supplementaire_echelle_temp_turb::readOn(is);}

void Diffusion_supplementaire_echelle_temp_turb_VDF::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VDF&                    domaine = ref_cast(Domaine_VDF, equation().domaine_dis());
  const Domaine_Cl_VDF&                     zcl = ref_cast(Domaine_Cl_VDF, equation().domaine_Cl_dis());
  const Champ_Inc_P0_base&                  tau = ref_cast(Champ_Inc_P0_base, equation().inconnue().valeur());
  const DoubleTab& tab_tau = semi_impl.count("tau") ? semi_impl.at("tau") : tau.passe(),
                   &nu_turb = ref_cast(Op_Dift_Multiphase_VDF_Elem, equation().operateur(0).l_op_base()).get_diffusivite_turbulente(),
                    &nu_visc  = ref_cast(Convection_Diffusion_std, equation()).diffusivite_pour_transport()->passe();

  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv();

  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  int N = tab_tau.dimension(1), ne = domaine.nb_elem(), ne_tot = domaine.nb_elem_tot(), nf_tot = domaine.nb_faces_tot(), D = dimension ;

  /* Calcul de grad de racine de tau aux faces */

  DoubleTrav grad_f_sqrt_tau(nf_tot, N);
  DoubleTrav sqrt_tau(ne_tot, N);
  for (int e = 0 ; e<ne_tot ; e++)
    for (int n=0 ; n< N ; n++) sqrt_tau(e,n) = std::sqrt(std::max(0., tab_tau(e,n)));
  assert ( tab_tau.dimension_tot(0) == ne_tot );
  const Convection_diffusion_turbulence_multiphase& eq_tau = ref_cast(Convection_diffusion_turbulence_multiphase, tau.equation());
  const Operateur_Grad& Op_Grad_tau = eq_tau.operateur_gradient_inconnue();
  // compute grad(sqrt(tau)) at faces ; Les CL sont des Dirichlet nulles, des symetries ou des Dirichlet non nulles => on devra corriger les dirichlet non nulles
  Op_Grad_tau.calculer(sqrt_tau,grad_f_sqrt_tau);

  // Boucle sur les bords pour corriger la contribution des cl de dirichlet non nulles au gradient
  for (int n_bord = 0; n_bord<domaine.nb_front_Cl(); n_bord++)
    for (int k=0; k<N; k++)
      {
        const Cond_lim& la_cl = zcl.les_conditions_limites(n_bord);
        const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
        const int ndeb = le_bord.num_premiere_face(), nfin = ndeb + le_bord.nb_faces();

        if (sub_type(Dirichlet,la_cl.valeur())) // Cas CL Dirichlet
          for (int num_face = ndeb, num_face_cl=0; num_face < nfin; num_face++, num_face_cl++)
            {
              int n0 = domaine.face_voisins(num_face,0);
              if (n0<0) n0 = domaine.face_voisins(num_face,1);
              int ori = domaine.orientation(num_face);
              grad_f_sqrt_tau(num_face, k) += -           ref_cast(Dirichlet, la_cl.valeur()).val_imp(num_face_cl, k)/(xp(n0,ori)- xv(num_face,ori));
              grad_f_sqrt_tau(num_face, k) -= - std::sqrt(ref_cast(Dirichlet, la_cl.valeur()).val_imp(num_face_cl, k))/(xp(n0,ori)- xv(num_face,ori));
            }
      }

  /* On interpole les champs aux elements */
  DoubleTab grad_e_sqrt_tau(ne_tot, N, D);
  face_to_elem(domaine, grad_f_sqrt_tau, grad_e_sqrt_tau);

  // On remplit le second membre

  for(int e = 0 ; e < ne ; e++)
    for (int n=0 ; n<N ; n++)
      {
        // Second membre
        double fac = pe(e)*ve(e) ;
        double secmem_en = 0 ;

        for (int d = 0; d<D; d++) secmem_en += grad_e_sqrt_tau(e, n, d)*grad_e_sqrt_tau(e, n, d);
        secmem(e, n) += fac * -8  * (nu_visc(e, n)+nu_turb(e,n)) * secmem_en;
      }
}

void Diffusion_supplementaire_echelle_temp_turb_VDF::face_to_elem(const Domaine_VF& domaine, const DoubleTab& tab_faces,DoubleTab& tab_elems) const
{
  const IntTab& elem_faces = domaine.elem_faces();
  const DoubleTab& nf = domaine.face_normales();
  const DoubleVect& fs = domaine.face_surfaces();
  const int nb_face_elem = elem_faces.line_size(), ne_tot= domaine.nb_elem_tot(), N = tab_faces.line_size(), D = dimension;
  assert (tab_elems.dimension_tot(0) == ne_tot && tab_faces.dimension_tot(0) == domaine.nb_faces_tot());

  tab_elems = 0.;
  for (int ele=0; ele<ne_tot; ele++)
    for (int n=0; n<N; n++)
      for (int d=0; d<D; d++)
        for (int s=0; s<nb_face_elem; s++)
          {
            int f = elem_faces(ele, s);
            if (fs[f]>1.e-12) tab_elems(ele, n, d) += tab_faces(f,n)*nf(f,d)/fs[f];
          }

  tab_elems *= 0.5;
}
