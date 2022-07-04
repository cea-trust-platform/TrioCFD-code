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
// File:        Injection_QDM_nulle_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Injection_QDM_nulle_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Champ_Face_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Zone_PolyMAC_P0.h>
#include <Conds_lim.h>
#include <Neumann_paroi.h>
#include <Neumann_val_ext.h>
#include <Matrix_tools.h>
#include <Masse_ajoutee_base.h>
#include <Array_tools.h>
#include <math.h>

Implemente_instanciable(Injection_QDM_nulle_PolyMAC_P0, "Injection_QDM_nulle_Face_PolyMAC_P0", Source_base);

Sortie& Injection_QDM_nulle_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Injection_QDM_nulle_PolyMAC_P0::readOn(Entree& is)
{
  return is;
}

void Injection_QDM_nulle_PolyMAC_P0::completer() // We must wait for all readOn's to be sure that the bubble dispersion and lift correlations are created
{
}

void Injection_QDM_nulle_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const // 100 % explicit ?
{
}

void Injection_QDM_nulle_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_Face_PolyMAC_P0& ch = ref_cast(Champ_Face_PolyMAC_P0, equation().inconnue().valeur());
  const Champ_Elem_PolyMAC_P0& cha= ref_cast(Champ_Elem_PolyMAC_P0, equation().probleme().equation(1).inconnue().valeur()); // volume fraction
  const Zone_PolyMAC_P0&     zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const Conds_lim&           clsa = cha.zone_Cl_dis().les_conditions_limites();

  const IntTab&  fcl = ch.fcl(),
                 &fcla = cha.fcl(),
                  &f_e = zone.face_voisins(),
                   &e_f = zone.elem_faces();
  const DoubleVect& vf = zone.volumes_entrelaces(),
                    &fs = zone.face_surfaces();
  const DoubleTab& vf_dir = zone.volumes_entrelaces_dir();

  const DoubleTab& vit = ch.valeurs(),
                   &rho   = equation().milieu().masse_volumique().passe(), // passe car qdm
                    &alpha = cha.passe();

  Matrice_Morse *mat = matrices.count(ch.le_nom().getString()) ? matrices.at(ch.le_nom().getString()) : NULL; // Derivee locale/QDM

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : NULL;
  const Masse_ajoutee_base *corr = pbm && pbm->has_correlation("masse_ajoutee") ? &ref_cast(Masse_ajoutee_base, pbm->get_correlation("masse_ajoutee").valeur()) : NULL;

  int N = vit.line_size(), D = dimension, nf_tot = zone.nb_faces_tot(), nf = zone.nb_faces();

  for (int f = 0 ; f< nf_tot ; f ++)
    if (fcla(f, 0) == 4) // Neumann_paroi
      {
        int e = f_e(f, 0) ;

        DoubleTrav f_a_masse(N, N) ;
        DoubleTrav f_a(1, N);
        for (int n = 0 ; n<N ; n++)
          {
            f_a_masse(n, n) = ref_cast(Neumann_paroi, clsa[fcla(f, 1)].valeur()).flux_impose(fcla(f, 2), n) * rho(e, n) ;   // Pas de porosite
            f_a(0, n)       = ref_cast(Neumann_paroi, clsa[fcla(f, 1)].valeur()).flux_impose(fcla(f, 2), n) ;   // Pas de porosite
          }
        corr->ajouter_inj( &f_a(0,0) , &alpha(e, 0),   &rho(e, 0)    ,   f_a_masse   );

        for (int d = 0; d<D ; d++)
          for (int n = 0 ; n<N ; n++)
            for (int m = 0 ; m<N ; m++)
              {
                secmem(nf_tot + D * e + d, n) -= fs(f) * f_a_masse(n, m)  * vit( nf_tot + D * e + d, m);  // Force along -vit
                if (mat)
                  (*mat)( N * (nf_tot + D*e + d) + n , N * (nf_tot + D*e + d) + m ) += fs(f) * f_a_masse(n, m) ;
              }

        for (int i=0 ; i < e_f.line_size() ; i++)
          {
            int f2 = e_f(e, i);
            if ( (f2 >= 0) && (f2<nf) && (fcl(f2, 0) < 2)  ) // Si pas face de bord
              {
                int c = ( e == f_e(f, 0) ) ? 0 : 1 ;
                for (int n = 0 ; n<N ; n++)
                  for (int m = 0 ; m<N ; m++)
                    {
                      secmem(f2, n) -= fs(f) * f_a_masse(n, m) * vf_dir(f2, c) / vf(f2)  * vit( f2, m);
                      if (mat)
                        (*mat)( N * f2 + n , N * f2 + m ) += fs(f) * f_a_masse(n, m) * vf_dir(f2, c) / vf(f2);
                    }
              }
          }
      }
}