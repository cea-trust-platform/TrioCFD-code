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

#include <Source_injection_masse_base.h>
#include <Pb_Multiphase.h>
#include <Discret_Thyd.h>
#include <Domaine_dis.h>

Implemente_base(Source_injection_masse_base, "Source_injection_masse_base", Sources_Multiphase_base);

Sortie& Source_injection_masse_base::printOn(Sortie& os) const {return os;}

Entree& Source_injection_masse_base::readOn(Entree& is)
{
  Cerr << "Lecture du Champ de masse injectee" << finl;

  Champ_Don flux_masse_lu_;
  Motcle type;
  is >> type;

  flux_masse_lu_.typer(type);
  Champ_Don_base& ch_flux_masse_lu_ = ref_cast(Champ_Don_base,flux_masse_lu_.valeur());
  is >> ch_flux_masse_lu_;
  const int nb_comp = ch_flux_masse_lu_.nb_comp();

  equation().probleme().discretisation().discretiser_champ("champ_elem", equation().domaine_dis(), "pp", "1",nb_comp,0., flux_masse_);
  flux_masse_lu_->fixer_nb_comp(nb_comp);
  if (ch_flux_masse_lu_.le_nom()=="anonyme") ch_flux_masse_lu_.nommer("Flux_masse_injectee");

  for (int n = 0; n < nb_comp; n++) flux_masse_lu_->fixer_nom_compo(n, ch_flux_masse_lu_.le_nom() + (nb_comp > 1 ? Nom(n) :""));
  for (int n = 0; n < nb_comp; n++) flux_masse_->fixer_nom_compo(n, ch_flux_masse_lu_.le_nom() + (nb_comp > 1 ? Nom(n) :""));
  equation().discretisation().nommer_completer_champ_physique(equation().domaine_dis(),ch_flux_masse_lu_.le_nom(),"1/s",flux_masse_lu_,equation().probleme());
  equation().discretisation().nommer_completer_champ_physique(equation().domaine_dis(),ch_flux_masse_lu_.le_nom(),"1/s",flux_masse_,equation().probleme());
  flux_masse_.valeur().valeurs() = 0;
  flux_masse_.valeur().affecter(flux_masse_lu_);

  const Pb_Multiphase& pb = ref_cast(Pb_Multiphase, equation().probleme());
  int N = pb.nb_phases();
  if (N != flux_masse_.valeurs().line_size()) Process::exit(que_suis_je() + " : you must input as many fluxes as there are phases !!");

  return is;
}

void Source_injection_masse_base::mettre_a_jour(double temps) { flux_masse_->mettre_a_jour(temps); }

void Source_injection_masse_base::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VF& domaine = ref_cast(Domaine_VF, equation().domaine_dis().valeur());
  const Pb_Multiphase& pb = ref_cast(Pb_Multiphase, equation().probleme());
  const DoubleTab& rho = pb.get_champ("masse_volumique").passe(),
            &tab_inj = flux_masse_.valeurs();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();
  int cI = (tab_inj.dimension(0)==1), cR = (rho.dimension(0)==1), N =rho.line_size();
  int e, k;

  for (e = 0; e < domaine.nb_elem(); e++)
    for (k = 0; k < N; k++)
      secmem(e, k) += pe[e] * ve[e] * tab_inj(!cI*e, k) * rho(!cR, k);

}