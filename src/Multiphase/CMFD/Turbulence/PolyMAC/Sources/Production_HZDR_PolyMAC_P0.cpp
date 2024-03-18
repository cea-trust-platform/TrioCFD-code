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
// File:        Production_HZDR_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_HZDR_PolyMAC_P0.h>

#include <Domaine_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>


Implemente_instanciable(Production_HZDR_PolyMAC_P0,"Production_HZDR_Elem_PolyMAC_P0", Source_base);

Sortie& Production_HZDR_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_HZDR_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("constante_gravitation", &g_);
  param.ajouter("C_k", &C_k);
  param.lire_avec_accolades_depuis(is);

  Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  return is;
}

void Production_HZDR_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
// empty : no derivative to add in the blocks
}

void Production_HZDR_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  // le terme de production s'écrit : C_k * (3/4 C_drag/d_b * alpha * rho_l * |u_r|**3) / (alpha_l * rho_l)
  // Le coefficient de trainée C_drag sera calculé avec le modèle de Tomiyama (codé en dur)

  const Domaine_PolyMAC_P0&             domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis().valeur());
  const DoubleTab&                      tab_rho = equation().probleme().get_champ("masse_volumique").passe();
  const DoubleTab&                      tab_alp = equation().probleme().get_champ("alpha").passe();
  const DoubleTab&                          vit = equation().probleme().get_champ("vitesse").passe();
  const DoubleTab&                         diam = equation().probleme().get_champ("diametre_bulles").valeurs();
  const DoubleTab&                           nu = equation().probleme().get_champ("viscosite_cinematique").passe();


  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();
  int N = ref_cast(Pb_Multiphase, equation().probleme()).nb_phases(), ne = domaine.nb_elem(), nf_tot = domaine.nb_faces_tot(), D = dimension ;

  // On récupère la tension superficielle sigma
  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());
  const DoubleTab& press = equation().probleme().get_champ("pression").passe();
  const DoubleTab& temp  = equation().probleme().get_champ("temperature").passe();
  const int nb_max_sat =  N * (N-1) /2; // oui !! suite arithmetique !!
  DoubleTrav Sigma_tab(ne,nb_max_sat);
  int Np = press.line_size();
  for (int k = 0; k < N; k++)
    {
      for (int l = k + 1; l < N; l++)
        {
          Interface_base& sat = milc.get_interface(k,l);
          const int ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // Et oui ! matrice triang sup !
          for (int i = 0 ; i<ne ; i++) Sigma_tab(i,ind_trav) = sat.sigma(temp(i,k),press(i,k * (Np > 1))) ;
        }
    }

  // On calcule le second membre aux elements : Prod_HZDR = C_k * (3/4 C_drag/d_b * alpha * rho_l * |u_r|**3) / (alpha_l * rho_l)
  for(int e = 0 ; e < ne ; e++)
    for (int k = 0 ; k<N ; k++)
      if (k!=n_l) // n_l est l'indice de la phase continue/liquide (n_l=0)
        {
          double u_r = 0;
          for (int d = 0; d < D; d++) u_r += (vit(nf_tot + D*e+d, k) - vit(nf_tot + D*e+d, n_l))*(vit(nf_tot + D*e+d, k) - vit(nf_tot + D*e+d, n_l)); // relative speed = gas speed - liquid speed
          u_r = std::sqrt(u_r);
          double Reb = diam(e,k)*u_r/nu(e,n_l); // Reynolds bulle
          int ind_trav = (k>n_l) ? (n_l*(N-1)-(n_l-1)*(n_l)/2) + (k-n_l-1) : (k*(N-1)-(k-1)*(k)/2) + (n_l-k-1);
          double Eo = g_ * std::abs(tab_rho(e, n_l)-tab_rho(e, k)) * diam(e, k)*diam(e, k)/Sigma_tab(e, ind_trav);
          double Cd = (u_r!=0) ? std::max( std::min( 16./Reb*(1+0.15*std::pow(Reb, 0.687)) , 48./Reb )   , 8.*Eo/(3.*(Eo+4.))) : 0; // Coef de trainee - remarque: si u_r=0 alors pas de trainée, pas de BIA donc production=0
          secmem(e, n_l) += ve(e) * pe(e) * C_k * (3./4.) * Cd / diam(e,k) * tab_alp(e, k) / tab_alp(e,n_l) * std::pow(u_r, 3);
        }
}
