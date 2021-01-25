/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Cond_Lim_rayo_semi_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Cond_Lim_rayo_semi_transp.h>
#include <Modele_rayo_semi_transp.h>
#include <Flux_radiatif_base.h>
#include <Frontiere_dis_base.h>
#include <Symetrie.h>

// Description:
//    Associe une equation a l'objet.
//    Affecte le membre Cond_Lim_rayo_semi_transp::mon_modele avec l'objet
//    passe en parametre.
// Precondition:
// Parametre: Modele_rayo_semi_transp& modele
//    Signification: le modele auquel on veut s'associer
//    Valeurs par defaut:
//    Contraintes: reference constante
//    Acces: entree
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Cond_Lim_rayo_semi_transp::associer_modele(const Modele_rayo_semi_transp& un_modele)
{
  mon_modele=un_modele;
}


void Cond_Lim_rayo_semi_transp::recherche_emissivite_et_A()
{
  // On recupere le modele
  int num_cl_rayo=0;

  Equation_rayonnement_base& eq_rayo = modele().eq_rayo();
  Conds_lim& les_cl_rayo = eq_rayo.zone_Cl_dis().les_conditions_limites();

  int test_nom=0;
  for(num_cl_rayo = 0; num_cl_rayo<les_cl_rayo.size(); num_cl_rayo++)
    {
      Cond_lim& la_cl_rayo = eq_rayo.zone_Cl_dis().les_conditions_limites(num_cl_rayo);
      Nom nom_cl_rayo = la_cl_rayo.frontiere_dis().le_nom();
      Nom nom_cl_temp = la_cl().frontiere_dis().le_nom();
      if(nom_cl_temp == nom_cl_rayo)
        {
          test_nom = 1;
          if(sub_type(Flux_radiatif_base,la_cl_rayo.valeur()))
            {
              const Flux_radiatif_base& la_cl_rayon = ref_cast(Flux_radiatif_base,la_cl_rayo.valeur());
              emissivite_ = la_cl_rayon.emissivite();
              emissivite_.associer_fr_dis_base(la_cl_rayon.frontiere_dis());
              A_ = la_cl_rayon.A();
            }
          else if (sub_type(Symetrie,la_cl_rayo.valeur()))
            {
              ;
            }
          else
            {
              Cerr<<"Erreur dans Frontiere_ouverte_rayo_semi_transp::recherche_emissivite_et_A()"<<finl;
              Cerr<<"Les conditions a utiliser pour l'equation de rayonnement "<<finl;
              Cerr<<"doivent forcement etre du type Flux_radiatif_base ou symetrie"<<finl;
              Process::exit();
            }
        }
    }
  if (test_nom == 0)
    {
      Cerr<<"Erreur dans Frontiere_ouverte_rayo_semi_transp::recherche_emissivite_et_A()"<<finl;
      Cerr<<"Probleme de compatibilite entre les conditions de l'equation de"<<finl;
      Cerr<<"rayonnement et l'equation de temperature"<<finl;
      Process::exit();
    }
}
