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
// File:        Source_Diffusion_croisee_echelle_temp_taux_diss_turb.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Diffusion_croisee_echelle_temp_taux_diss_turb_VDF.h>

#include <Convection_diffusion_turbulence_multiphase.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Champ_Inc_P0_base.h>
#include <Echange_impose_base.h>
#include <Domaine_VF.h>
#include <Domaine_Cl_PolyMAC.h>
#include <Pb_Multiphase.h>
#include <Dirichlet.h>

Implemente_instanciable(Diffusion_croisee_echelle_temp_taux_diss_turb_VDF,"Diffusion_croisee_echelle_temp_taux_diss_turb_VDF_P0_VDF", Source_Diffusion_croisee_echelle_temp_taux_diss_turb);

Sortie& Diffusion_croisee_echelle_temp_taux_diss_turb_VDF::printOn(Sortie& os) const {  return Source_Diffusion_croisee_echelle_temp_taux_diss_turb::printOn(os);}

Entree& Diffusion_croisee_echelle_temp_taux_diss_turb_VDF::readOn(Entree& is) { return Source_Diffusion_croisee_echelle_temp_taux_diss_turb::readOn(is); }

void Diffusion_croisee_echelle_temp_taux_diss_turb_VDF::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VF& domaine = ref_cast(Domaine_VF, equation().domaine_dis().valeur());
  const Champ_Inc_P0_base& ch_k = ref_cast(Champ_Inc_P0_base, equation().probleme().get_champ("k"));	// Champ k
  const DoubleTab& k_passe = ch_k.passe();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  const Champ_Inc_P0_base& ch_diss	= ref_cast(Champ_Inc_P0_base, equation().inconnue().valeur()); 		// Champ tau ou omega
  const DoubleTab& diss_passe = ch_diss.passe();
  const DoubleTab& diss = ch_diss.valeurs();

  const int nf_tot = domaine.nb_faces_tot(), D = dimension, ne = domaine.nb_elem(), ne_tot = domaine.nb_elem_tot() ;
  const int N = diss_passe.line_size();

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") abort();

//  Cerr <<"lol !! " ;

  assert(N == 1 || k_passe.line_size() == 1); // si Ntau > 1 il vaut mieux iterer sur les id_composites des phases turbulentes decrites par un modele k-tau dans le calcul de grad_f_k et dans le remplissage des matrices

  /* Calcul de grad de tau ou de omega aux faces */

  DoubleTab grad_f_diss(nf_tot, N);
  assert ( diss_passe.dimension_tot(0) == ne_tot );
  const Convection_diffusion_turbulence_multiphase& eq_diss = ref_cast(Convection_diffusion_turbulence_multiphase,equation());
  const Operateur_Grad& Op_Grad_diss = eq_diss.operateur_gradient_inconnue();
  Op_Grad_diss.calculer(diss_passe,grad_f_diss); // compute grad(diss) at faces

  /* Calcul de grad de k aux faces */

  DoubleTab grad_f_k(nf_tot, N);
  assert ( k_passe.dimension_tot(0) == ne_tot );
  const Convection_diffusion_turbulence_multiphase& eq_k = ref_cast(Convection_diffusion_turbulence_multiphase, ch_k.equation());
  const Operateur_Grad& Op_Grad_k = eq_k.operateur_gradient_inconnue();
  Op_Grad_k.calculer(k_passe,grad_f_k); // compute grad(diss) at faces

  /* On interpole les champs aux elements */
  DoubleTab grad_e_diss(ne_tot, N, D);
  DoubleTab grad_e_k(ne_tot, N, D);
  face_to_elem(domaine, grad_f_diss, grad_e_diss);
  face_to_elem(domaine, grad_f_k   , grad_e_k);

  /* Calcul de grad(tau/omega).(grad k) */

  DoubleTab grad_f_diss_dot_grad_f_k(ne, N);
  for (int n = 0; n < N; n++)
    for (int e = 0; e < ne; e++)
      {
        grad_f_diss_dot_grad_f_k(e, n) = 0;
        for (int d = 0 ; d < D ; d++) grad_f_diss_dot_grad_f_k(e, n) += grad_e_diss(e,n,d) * grad_e_k(e,n,d); // produit scalaire
      }

  /* remplissage des matrices et du second membre */

  Matrice_Morse *M = matrices.count(ch_diss.le_nom().getString()) ? matrices.at(ch_diss.le_nom().getString()) : nullptr;

  int e, n;

  for ( e = 0; e < ne; e++)
    for(n = 0; n<N ; n++)
      {
        if (Type_diss == "tau")
          {
            double secmem_en = pe(e) * ve(e) * sigma_d * diss(e, n) * std::min(grad_f_diss_dot_grad_f_k(e, n), 0.);
            secmem(e, n) += secmem_en;
            if (!(M==nullptr))     (*M)(N * e + n, N * e + n)       -= pe(e) * ve(e) * sigma_d * std::min(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee en tau
          }
        else if (Type_diss == "omega")
          {
            if (diss(e,n)>1.e-8) // Else everything = 0
              {
                double dp = std::max(diss_passe(e, n), 1.e-6);
                secmem(e, n) += pe(e) * ve(e) * sigma_d / dp*(2-diss(e, n)/dp)* std::max(grad_f_diss_dot_grad_f_k(e, n), 0.) ;
                if (!(M==nullptr))     (*M)(N * e + n, N * e + n)       -= pe(e) * ve(e) * sigma_d * (-1/(dp*dp)) * std::max(grad_f_diss_dot_grad_f_k(e, n), 0.); // derivee en omega
              }
          }
      }
}

void Diffusion_croisee_echelle_temp_taux_diss_turb_VDF::face_to_elem(const Domaine_VF& domaine, const DoubleTab& tab_faces,DoubleTab& tab_elems) const
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

