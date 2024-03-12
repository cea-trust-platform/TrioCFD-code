/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Navier_Stokes_QC_impl.h
// Directory:   $TRIO_U_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Navier_Stokes_QC_impl_included
#define Navier_Stokes_QC_impl_included

#include <Champ_Inc.h>
#include <TRUST_Ref.h>

class Navier_Stokes_std;
class DoubleTab;
class Entree;
class Sortie;
class Motcle;
class Fluide_Incompressible;
class Matrice;
class Assembleur;
class Matrice_Morse;

/*! @brief classe Navier_Stokes_QC_impl Cette classe porte les termes de l'equation de la dynamique
 *
 *     pour un fluide sans modelisation de la turbulence.
 *     On suppose l'hypothese de fluide quasi compressible.
 *     Sous ces hypotheses, on utilise la forme suivante des equations de
 *     Navier_Stokes:
 *        DU/dt = div(terme visqueux) - gradP/rho + sources/rho
 *        div U = W
 *     avec DU/dt : derivee particulaire de la vitesse
 *          rho   : masse volumique
 *     Rq : l'implementation de la classe permet bien sur de negliger
 *          certains termes de l'equation (le terme visqueux, le terme
 *          convectif, tel ou tel terme source).
 *     L'inconnue est le champ de vitesse.
 *
 * @sa Equation_base Pb_Thermohydraulique_QC Navier_Stokes_std
 */
class Navier_Stokes_QC_impl
{
public :
  Entree& lire_impl(Navier_Stokes_std&,const Motcle&, Entree&);
  void completer_impl(Navier_Stokes_std&);
  int impr_impl(const Navier_Stokes_std&,Sortie& os) const;
  int preparer_calcul_impl(Navier_Stokes_std&);
  void projeter_impl(Navier_Stokes_std&);
  //inline Champ_Inc& rho_la_vitesse();

  DoubleTab& derivee_en_temps_inco_impl(Navier_Stokes_std&,DoubleTab&,  Fluide_Incompressible& le_fluide,Matrice& matrice_pression_,Assembleur& assembleur_pression_ );
  DoubleTab& derivee_en_temps_inco_diffusion_impl(Navier_Stokes_std&,DoubleTab&,  Fluide_Incompressible& le_fluide,Matrice& matrice_pression_,Assembleur& assembleur_pression_ );
  void assembler_impl(const Navier_Stokes_std&, Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) ;
  void assembler_avec_inertie_impl(const Navier_Stokes_std&, Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) ;

  /*
    int comprend_champ_impl(Navier_Stokes_std&,const Motcle& ) const;
    int a_pour_Champ_Inc_impl(Navier_Stokes_std&,const Motcle&, REF(Champ_base)& ) const;
    int a_pour_Champ_Fonc_impl(Navier_Stokes_std&,const Motcle&, REF(Champ_base)& ) const;
  */
  void discretiser_impl(Navier_Stokes_std&);
  DoubleTab& derivee_en_temps_impl_p1(Navier_Stokes_std&,DoubleTab&,  Fluide_Incompressible& le_fluide,Matrice& matrice_pression_,Assembleur& assembleur_pression_, int is_implicite);
  DoubleTab& derivee_en_temps_impl_projection(Navier_Stokes_std&,DoubleTab&,  Fluide_Incompressible& le_fluide,Matrice& matrice_pression_,Assembleur& assembleur_pression_, int is_implicite);
protected:

  DoubleTab& derivee_en_temps_inco(Navier_Stokes_std&,DoubleTab&,  Fluide_Incompressible& le_fluide,Matrice& matrice_pression_,Assembleur& assembleur_pression_, int is_implicite);

  Champ_Inc rho_la_vitesse_;
  DoubleTab tab_W;

  IntVect orientation_VDF_;
};

DoubleTab& rho_vitesse_impl(const DoubleTab& tab_rho,const DoubleTab& vitesse,DoubleTab& rhovitesse) ;


#endif
