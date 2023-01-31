/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Source_BIF_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/Multiphase/CMFD/src/Turbulence/PolyMAC/Sources
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_BIF_PolyMAC_P0.h>

#include <Domaine_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Matrix_tools.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_PolyMAC_P0.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_base.h>
#include <Viscosite_turbulente_multiple.h>
#include <TRUSTTab_parts.h>


Implemente_instanciable(Source_BIF_PolyMAC_P0,"Source_BIF_Face_PolyMAC_P0", Source_base);

Sortie& Source_BIF_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Source_BIF_PolyMAC_P0::readOn(Entree& is)
{
  if (!sub_type(Viscosite_turbulente_multiple, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, ref_cast(Navier_Stokes_std, equation().probleme().equation(0)).operateur(0).l_op_base()).correlation().valeur())) Process::exit(que_suis_je() + " : the turbulence correlation must be multiple");

  return is;
}

void Source_BIF_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
// empty : no derivative to add in the blocks
}

void Source_BIF_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0&                      domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const Op_Diff_Turbulent_PolyMAC_P0_Face& Op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base());
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();

  const DoubleTab&                       vf_dir = domaine.volumes_entrelaces_dir(), &xp = domaine.xp(), &xv = domaine.xv();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes(), &fs = domaine.face_surfaces();
  const DoubleTab& normales_f = domaine.face_normales();
  const IntTab& voisins_f = domaine.face_voisins(), &e_f = domaine.elem_faces(), &f_e = domaine.face_voisins();

  const Viscosite_turbulente_multiple&    visc_turb = ref_cast(Viscosite_turbulente_multiple, Op_diff.correlation().valeur());

  int N = pb.get_champ("vitesse").valeurs().dimension(1), D = dimension, nf_tot = domaine.nb_faces_tot(), nf = domaine.nb_faces(), ne_tot = domaine.nb_elem_tot() ;

  // On recupere les tensions de reynolds des termes de BIF
  DoubleTrav Rij(0, N, D, D);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), Rij); //Necessary to compare size in reynolds_stress()
  visc_turb.reynolds_stress_BIF(Rij);

  DoubleTrav grad_Rij(0, N, D, D);
  MD_Vector_tools::creer_tableau_distribue(eq_qdm.vitesse()->valeurs().get_md_vector(), grad_Rij); //Necessary to exchange virtual elements after calculation of the gradient at the faces()

  const Champ_Elem_PolyMAC_P0& ch_alpha = ref_cast(Champ_Elem_PolyMAC_P0, equation().probleme().get_champ("alpha"));	// Champ alpha qui servira à obtenir les coeffs du gradient ; normalement toujours des CAL de Neumann ; terme source qui n'apparait qu'en multiphase
  ch_alpha.init_grad(0); // Initialisation des tables fgrad_d, fgrad_e, fgrad_w qui dependent de la discretisation et du type de conditions aux limites --> pas de mises a jour necessaires
  IntTab& fg_d = ch_alpha.fgrad_d, fg_e = ch_alpha.fgrad_e; // Tables utilisees dans domaine_PolyMAC_P0::fgrad pour le calcul du gradient
  DoubleTab fg_w = ch_alpha.fgrad_w;

  // On calcule le gradient de Rij aux faces
  for (int n = 0; n < N; n++)
    for (int f = 0; f < nf_tot; f++)
      for (int d_i = 0; d_i <D ; d_i++)
        for (int d_j = 0; d_j <D ; d_j++)
          {
            grad_Rij(f, n, d_i, d_j) = 0;
            // grad_Rij(face, phase, x-coord, y-coord)   or
            // grad_Rij(nb_face_tot + dimension*element*gradient_component, phase, x-coord, y-coord)
            for (int j = fg_d(f); j < fg_d(f+1) ; j++)
              {
                int e = fg_e(j);
                int f_bord;
                if ( (f_bord = e - ne_tot) < 0) //contribution d'un element
                  grad_Rij(f, n, d_i, d_j) += fg_w(j) * Rij(e, n, d_i, d_j);
                else if ( (ch_alpha.fcl()(f_bord, 0) == 1) || (ch_alpha.fcl()(f_bord, 0) == 2)
                          || (ch_alpha.fcl()(f_bord, 0) == 3) || (ch_alpha.fcl()(f_bord, 0) == 6))
                  {
                    Process::exit("You must have a neumann limit condition on alpha for RIJ_BIF to work !");
                  }
              }
          }

  // On interpole le gradient de Rij aux elements
  for (int n = 0; n < N; n++)
    for (int e = 0; e < ne_tot; e++)
      for (int d_d=0 ; d_d<D ; d_d++) // on derive / d_d
        for (int d_i = 0; d_i <D ; d_i++)
          for (int d_j = 0; d_j <D ; d_j++)
            {
              grad_Rij(nf_tot + D *e + d_d, n, d_i,  d_j) = 0;
              for (int j = 0, f; j < e_f.dimension(1) && (f = e_f(e, j)) >= 0; j++)
                grad_Rij(nf_tot + D *e + d_d, n, d_i,  d_j) += (e == f_e(f, 0) ? 1 : -1) * fs(f) * (xv(f, d_d) - xp(e, d_d)) / ve(e) * grad_Rij(f, n, d_i,  d_j);
            }

  // On calcule le second membre aux elements

  for(int e = 0 ; e < ne_tot ; e++)
    for(int n = 0; n<N ; n++)
      for (int d_i = 0; d_i < D; d_i++)
        {
          double secmem_en = 0;
          for (int d_j = 0; d_j < D; d_j++)
            secmem_en += grad_Rij(nf_tot + D *e + d_j, n, d_i,  d_j) ;
          secmem_en *= (-1) * pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) ; // For us, Rij = < u_i u_j >, therefore *(-1)
          secmem(nf_tot + e*D + d_i, n) += secmem_en;
        }

  // On calcule le second membre aux faces

  int e;

  for(int f = 0 ; f < nf ; f++)
    for(int n = 0; n<N ; n++)
      for (int i = 0; i < 2 && (e = voisins_f(f, i)) >= 0; i++)
        {
          DoubleTrav secmem_en(3); // Contains the vector of the divergence of R_ij at the face
          secmem_en = 0;
          for (int d_i = 0; d_i < D; d_i++)
            for (int d_j = 0; d_j < D; d_j++)
              secmem_en(d_i) += grad_Rij(nf_tot + D *e + d_j, n, d_i,  d_j) ;
          for (int d_i = 0; d_i < D; d_i++)
            secmem_en(d_i) *= (-1) * pe(e) * vf_dir(f, i) * tab_alp(e, n) * tab_rho(e, n);// For us, Rij = < u_i u_j >, therefore *(-1)
          double flux_face = domaine.dot(&normales_f(f, 0), &secmem_en(0));
          secmem(f, n) += flux_face;
        }

}
