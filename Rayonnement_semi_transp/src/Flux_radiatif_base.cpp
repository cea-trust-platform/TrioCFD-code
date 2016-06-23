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
// File:        Flux_radiatif_base.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#include <Flux_radiatif_base.h>
#include <Motcle.h>
#include <Front_VF.h>

Implemente_base(Flux_radiatif_base,"Flux_radiatif_base",Neumann_paroi);

// Description:
//    Imprime le type de l'equation sur un flot de sortie.
// Precondition:
// Parametre: Sortie& s
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Flux_radiatif_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}


// Description:
//    Appel Neumann_paroi::readOn(Entree& is)
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception: solveur pression non defini dans jeu de donnees
// Effets de bord:
// Postcondition:
Entree& Flux_radiatif_base::readOn(Entree& is)
{
  Nom nom_pb, nom_bord;
  Motcle nom_champ;

  Motcle motlu;
  Motcles les_motcles(2);
  {
    les_motcles[0] = "emissivite";
    les_motcles[1] = "A";
  }

  int ind = 0;
  while (ind < 2)
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);

      switch(rang)
        {
        case 0:
          {
            is >> emissivite_;
            break;
          }
        case 1:
          {
            is >> A_;
            break;
          }
        default:
          {
            Cerr << "Erreur a la lecture de la condition aux limites de type "<<finl;
            Cerr << "Flux_radiatif_base " << finl;
            Cerr << "On attendait " << les_motcles << "a la place de " <<  motlu << finl;
            exit();
          }
        }
      ind++;
    }

  champ_front().typer("Champ_front_fonc");
  champ_front()->fixer_nb_comp(1);

  flux_radiatif().typer("Champ_front_fonc");
  flux_radiatif()->fixer_nb_comp(1);

  return is;
}



// Description:
//    cette methode permet de typer les champs_front :
//        - le_champ_front
//        - flux_radiatif_
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Flux_radiatif_base::completer()
{
  Neumann_paroi::completer();
  // On type le champ_front flux_radiatif_ qui est associe a la condition a la limite
  const Front_VF& front_vf=ref_cast(Front_VF, le_champ_front.frontiere_dis());
  int nb_comp = 1;

  flux_radiatif()->nommer(front_vf.le_nom());
  DoubleTab& tab_flux = flux_radiatif().valeurs();
  tab_flux.resize(front_vf.nb_faces(),nb_comp);

  champ_front()->nommer(front_vf.le_nom());
  DoubleTab& tab= champ_front().valeurs();
  tab.resize(front_vf.nb_faces(),nb_comp);
  emissivite_.associer_fr_dis_base(front_vf);
}


// Description:
//   Renvoie la valeur de flux imposes a la paroi radiative
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
double Flux_radiatif_base::flux_impose(int i) const
{
  if (le_champ_front.valeurs().size()==1)
    return le_champ_front(0,0);
  else if (le_champ_front.valeurs().dimension(1)==1)
    return le_champ_front(i,0);
  else
    Cerr << "Flux_radiatif_base::flux_impose erreur" << finl;
  exit();
  return 0.;
}


// Description:
//   Renvoie la valeur de flux imposes a la paroi radiative
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
double Flux_radiatif_base::flux_impose(int i,int j) const
{
  if (le_champ_front.valeurs().dimension(0)==1)
    return le_champ_front(0,j);
  else
    return le_champ_front(i,j);
}
