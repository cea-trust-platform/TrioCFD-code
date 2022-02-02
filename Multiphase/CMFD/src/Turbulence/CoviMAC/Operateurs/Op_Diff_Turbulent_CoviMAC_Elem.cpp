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
// File:        Op_Diff_CoviMAC_Elem.cpp
// Directory:   $TRUST_ROOT/src/CoviMAC/Operateurs
// Version:     1
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_Turbulent_CoviMAC_Elem.h>
#include <Op_Diff_Turbulent_CoviMAC_Face.h>
#include <Pb_Multiphase.h>
#include <Viscosite_turbulente_base.h>
#include <Transport_turbulent_base.h>
#include <Champ_P0_CoviMAC.h>
#include <CoviMAC_discretisation.h>

Implemente_instanciable( Op_Diff_Turbulent_CoviMAC_Elem, "Op_Diff_Turbulent_CoviMAC_Elem|Op_Diff_Turbulente_CoviMAC_Elem", Op_Diff_CoviMAC_Elem ) ;

Sortie& Op_Diff_Turbulent_CoviMAC_Elem::printOn( Sortie& os ) const
{
  Op_Diff_CoviMAC_base::printOn( os );
  return os;
}

Entree& Op_Diff_Turbulent_CoviMAC_Elem::readOn( Entree& is )
{
  //lecture de la correlation de viscosite turbulente
  corr.typer_lire(equation().probleme(), "transport_turbulent", is);

  champs_compris_.ajoute_nom_compris(nom_secmem_diff_ = Nom("Diffusion_turbulente_") + equation().que_suis_je());

  return is;
}

void Op_Diff_Turbulent_CoviMAC_Elem::creer_champ(const Motcle& motlu)
{
  if ((nom_secmem_diff_==motlu) && !(secmem_diff_.non_nul()))
    {
      int N = ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur()).valeurs().dimension(1) ;
      Noms noms(N), unites(N);
      for (int n = 0 ; n<N ; n++) {unites[0] = "m2/s"; noms[0] = nom_secmem_diff_;}
      Motcle typeChamp = "champ_elem" ;
      ref_cast(CoviMAC_discretisation, equation().discretisation()).discretiser_champ(typeChamp, equation().zone_dis(), scalaire, noms , unites, N, 0, secmem_diff_);
      champs_compris_.ajoute_champ(secmem_diff_);
    }
  Op_Diff_CoviMAC_Elem::creer_champ(motlu);
}

void Op_Diff_Turbulent_CoviMAC_Elem::mettre_a_jour(double temps)
{
  if (secmem_diff_.non_nul())
    {
      const Zone_CoviMAC&                      zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
      const Champ_P0_CoviMAC&                   tau = ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur());
      const DoubleTab&                      tab_tau = tau.valeurs();

      DoubleTrav sec_m(tab_tau); //residus
      matrices_t mat_m; //derivees vides
      tabs_t semi_impl;
      semi_impl[equation().inconnue().le_nom().getString()] = equation().inconnue().passe();
      ajouter_blocs(mat_m, sec_m, semi_impl); 

      int N = tab_tau.dimension(1), ne = zone.nb_elem() ;

      for(int e = 0 ; e < ne ; e++) for (int n=0 ; n<N ; n++)
          secmem_diff_(e,n) = sec_m(e, n) ;
    }
  Op_Diff_CoviMAC_Elem::mettre_a_jour(temps);
}


void Op_Diff_Turbulent_CoviMAC_Elem::completer()
{
  Op_Diff_CoviMAC_Elem::completer();
  //si la correlation a besoin du gradient de u, on doit le creer maintenant
  if (corr.non_nul() && ref_cast(Transport_turbulent_base, corr.valeur()).gradu_required())
    equation().probleme().creer_champ("gradient_vitesse");
}

void Op_Diff_Turbulent_CoviMAC_Elem::modifier_nu(DoubleTab& mu) const
{
  if (!corr.non_nul()) return; //rien a faire
  const Operateur_base& op_qdm = equation().probleme().equation(0).operateur(0).l_op_base();
  if (!sub_type(Op_Diff_Turbulent_CoviMAC_Face, op_qdm))
    Process::exit(que_suis_je() + ": no turbulent momentum diffusion found!");
  const Correlation& corr_visc = ref_cast(Op_Diff_Turbulent_CoviMAC_Face, op_qdm).corr;
  if (!corr.non_nul() || !sub_type(Viscosite_turbulente_base, corr_visc.valeur()))
    Process::exit(que_suis_je() + ": no turbulent viscosity correlation found!");
  //un "simple" appel a la correlation!
  ref_cast(Transport_turbulent_base, corr.valeur()).modifier_nu(ref_cast(Convection_Diffusion_std, equation()), ref_cast(Viscosite_turbulente_base, corr_visc.valeur()), mu);
}