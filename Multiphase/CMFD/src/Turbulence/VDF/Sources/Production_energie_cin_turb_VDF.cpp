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
// File:        Production_energie_cin_turb_VDF.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_energie_cin_turb_VDF.h>

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Viscosite_turbulente_base.h>
#include <Navier_Stokes_std.h>
#include <Domaine_Cl_VDF.h>
#include <Champ_Face_VDF.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>

Implemente_instanciable(Production_energie_cin_turb_VDF,"Production_energie_cin_turb_VDF_P0_VDF", Source_Production_energie_cin_turb);
// XD Production_energie_cin_turb source_base Production_energie_cin_turb 1 Production source term for the TKE equation


Sortie& Production_energie_cin_turb_VDF::printOn(Sortie& os) const {return Source_Production_energie_cin_turb::printOn(os);}
Entree& Production_energie_cin_turb_VDF::readOn(Entree& is) { return Source_Production_energie_cin_turb::readOn(is);}

void Production_energie_cin_turb_VDF::completer() 
{
  const Navier_Stokes_std&     eq_qdm 	= ref_cast(Navier_Stokes_std, equation().probleme().equation(0));
  if (ref_cast(Operateur_Diff_base, eq_qdm.operateur(0).l_op_base()).correlation_viscosite_turbulente()==nullptr) Process::exit(que_suis_je() + " : the momentum diffusion must be turbulent !");
}


void Production_energie_cin_turb_VDF::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VF&                   domaine = ref_cast(Domaine_VF, equation().domaine_dis().valeur());
  const Probleme_base&                     pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&             eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const Viscosite_turbulente_base&  visc_turb = ref_cast(Viscosite_turbulente_base, (*ref_cast(Operateur_Diff_base, eq_qdm.operateur(0).l_op_base()).correlation_viscosite_turbulente()).valeur());
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if      sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }


  int Nph = pb.get_champ("vitesse").valeurs().dimension(1), nb_elem = domaine.nb_elem(), D = dimension, ne_tot = domaine.nb_elem_tot() ;
  int N = equation().inconnue()->valeurs().line_size(),
      Np = equation().probleme().get_champ("pression").valeurs().line_size(),
      Nt = equation().probleme().get_champ("temperature").valeurs().line_size(),
      Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;
  int e, n, mp;

  const Champ_Face_VDF& ch_vit = ref_cast(Champ_Face_VDF, eq_qdm.inconnue().valeur());
  DoubleTrav    tab_grad(ne_tot, D, D);//, tab_vit_liq(nf_tot);  
//  boucle remplit  tab_vit_liq;
  ch_vit.calcul_duidxj(ch_vit.passe(), tab_grad, ref_cast(Domaine_Cl_VDF, ch_vit.domaine_Cl_dis().valeur()));

  double limiter_ = visc_turb.limiteur();
  double nut_l = -10000., fac;

  const DoubleTab& tab_rho = equation().probleme().get_champ("masse_volumique").passe(),
                   &tab_alp = equation().probleme().get_champ("alpha").passe(),
                    &nu =  equation().probleme().get_champ("viscosite_cinematique").passe(),
                     &k = equation().probleme().get_champ("k").valeurs(),
                      *diss = equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).valeurs() : nullptr,
                       *pdiss = equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).passe() : nullptr;

  const Champ_base&   ch_alpha_rho = sub_type(Pb_Multiphase,equation().probleme()) ? ref_cast(Pb_Multiphase,equation().probleme()).eq_masse.champ_conserve() : equation().milieu().masse_volumique().valeur();
  const DoubleTab&       alpha_rho = ch_alpha_rho.valeurs();
  const tabs_t&      der_alpha_rho = ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees

  if (Type_diss == "")
    {
      DoubleTrav nut(0, Nph);
      MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), nut); //Necessary to compare size in eddy_viscosity()
      visc_turb.eddy_viscosity(nut);

      Matrice_Morse *mat = matrices.count("k") ? matrices.at("k") : nullptr;

      for( e = 0 ; e < nb_elem ; e++)
        for( n = 0; n<N ; n++)
          {
            double secmem_en = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                secmem_en += ( tab_grad( e, d_U, d_X) + tab_grad( e, d_X, d_U) ) * tab_grad( e, d_U, d_X) ;
            secmem_en *= pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) * nut(e, n) ;

            secmem(e, n) += std::max(secmem_en, 0.);//  secmem(e, n) += fac_(e, n, 0) * std::max(secmem_en, 0.);

            if (mat) (*mat)(N * e + n, N * e + n) -= 0.;//fac_(e, n, 1) * std::max(secmem_en, 0.);
          }
    }

  else
    {
      for( e = 0 ; e < nb_elem ; e++)
        for(n = 0, mp = 0; n<N ; n++, mp += (Np > 1))
          {
            double grad_grad = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                grad_grad += ( tab_grad( e, d_U, d_X) + tab_grad( e, d_X, d_U) ) * tab_grad( e, d_U, d_X) ;

            fac = std::max(grad_grad, 0.) * pe(e) * ve(e) ;

            if      (Type_diss == "tau")   nut_l =                         std::max(k(e, n) * (*diss)(e, n), limiter_ * nu(e, n)) ;
            else if (Type_diss == "omega") nut_l = ( ((*pdiss)(e,n) > 0.) ? std::max(k(e, n) / (*pdiss)(e, n)*(2-(*diss)(e, n)/(*pdiss)(e, n)), limiter_ * nu(e, n)) : limiter_ * nu(e, n) );
            else Process::exit(que_suis_je() + " : ajouter_blocs : probleme !!!") ;

            secmem(e, n) += fac * alpha_rho(e, n) * nut_l ;
            for (auto &&i_m : matrices)
              {
                Matrice_Morse& mat = *i_m.second;
                if (i_m.first == "alpha") 	    mat(N * e + n, Na * e + n) -= fac * nut_l * (der_alpha_rho.count("alpha") ?       der_alpha_rho.at("alpha")(e, n) : 0 );			      // derivee par rapport au taux de vide
                if (i_m.first == "temperature") mat(N * e + n, Nt * e + n) -= fac * nut_l * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, n) : 0 );// derivee par rapport a la temperature
                if (i_m.first == "pression")    mat(N * e + n, Np * e + mp)-= fac * nut_l * (der_alpha_rho.count("pression") ?    der_alpha_rho.at("pression")(e, mp) : 0 );		  // derivee par rapport a la pression
              }

            if ( (Type_diss == "tau") && ((k(e, n)*(*diss)(e, n)) > (limiter_*nu(e, n))) )
              for (auto &&i_m : matrices)
                {
                  Matrice_Morse& mat = *i_m.second;
                  if (i_m.first == "k")         mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) *(*diss)(e,n);
                  if (i_m.first == "tau")       mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) * k(e,n);
                }
            if ( (Type_diss == "omega") && ((*pdiss)(e,n) > 0.) && ((k(e, n)/(*pdiss)(e, n)*(2-(*diss)(e, n)/(*pdiss)(e, n)) > (limiter_*nu(e, n))) ) )
              for (auto &&i_m : matrices)
                {
                  Matrice_Morse& mat = *i_m.second;
                  if (i_m.first == "k")         mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) *      1./(*pdiss)(e, n)*(2-(*diss)(e, n)/(*pdiss)(e, n)) ;
                  if (i_m.first == "omega")     mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) * -k(e,n)/((*pdiss)(e, n)*(*pdiss)(e, n)) ;
                }
          }
    }
}
