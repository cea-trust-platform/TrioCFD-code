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
// File:        Champ_front_zoom.cpp
// Directory:   $TRUST_ROOT/src/Zoom/Kernel
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#include <Champ_front_zoom.h>
#include <Interprete.h>
#include <Equation_base.h>
#include <Pb_MG.h>

Implemente_instanciable(Champ_front_zoom,"Champ_front_zoom",Ch_front_var_instationnaire_dep);


// Description:
//    NE FAIT RIEN
// Precondition:
// Parametre: Sortie& os
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Champ_front_zoom::printOn(Sortie& os) const
{
  return os;
}

// Description:
//    Lit le nom d'un probleme multigrille, le nom d'un probleme base
//    qui est le probleme fin et
//    la nature des conditions limites a partir d'un flot d'entree.
//    Cree ensuite le champ de frontiere correspondant.
//    Format:
//      Champ_front_zoom nom_pbMG nom_pbF nature_cond_lim
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entre
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Champ_front_zoom::readOn(Entree& is)
{
  //Cerr << "Dans Champ_front_zoom::readOn(Entree& is) " << finl;
  Nom nom_pbMG;
  Nom nom_pbG;
  Nom nom_pbF;
  Nom bord;
  Motcle nom_inco;
  Nom nature;
  is >>  nom_pbMG >> nom_pbG >>  nom_pbF >> bord >> nom_inco ;
  creer(nom_pbMG,nom_pbG,nom_pbF, bord,nom_inco);
  return is;
}

// Description:
//    Cree l'objet Champ_front_zoom representant le prolongement
//    de l'inconnue grossiere sur le bord fin a partir des noms:
//         - du probleme multigrille portant le probleme a 2 grilles relatif
//           au pbG et au pbF,portant lui-meme les connectivites
//           necessaires au prolongement
//         - du probleme fin
//         - du nom de l'inconnue
//         - de la nature des conditions limites (dirichlet ou neumann)
// Precondition:
// Parametre: Nom& nom_pbMG
//    Signification: le nom du probleme MG auquel appartient le pb_2G contenant
//                   les problemes fin et grossier, ainsi que le tableau
//                   de connectivites
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree
// Parametre: Nom& nom_pbF
//    Signification: le nom du probleme fin pour lequel on veut
//                   prendre le prolongement de l'inconnue grossiere
//                   sur sa frontiere
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree
// Parametre: Motcle& nom_inco
//    Signification: le nom de l'inconnue
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree
// Parametre: Nom& nature
//    Signification: flag : dirichlet ou neumann
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Champ_front_zoom::creer(Nom& nom_pbMG,
                             Nom& nom_pbG, Nom& nom_pbF, Nom& nom_front,
                             Motcle& nom_inco)
{
  //  Cerr<<"Creation du champ_front_zoom"<<finl;
  Objet_U& ob1=Interprete::objet(nom_pbMG);
  Objet_U& ob2=Interprete::objet(nom_pbF);
  Objet_U& ob3=Interprete::objet(nom_pbG);
  nom_bord = nom_front;

  if(!sub_type(Pb_MG, ob1) )
    {
      Cerr << " Error in the type of problem " << nom_pbMG << finl;
      exit();
    }
  else if(!sub_type(Probleme_base, ob2) )
    {
      Cerr << " Error in the type of problem " << nom_pbF << finl;
      exit();
    }
  else if(!sub_type(Probleme_base, ob3) )
    {
      Cerr << " Error in the type of problem " << nom_pbG << finl;
      exit();
    }
  else
    {
      //reference au probleme multigrille
      Pb_MG& pbMG = ref_cast(Pb_MG, ob1);
      le_pb_MG_ = pbMG;
      //reference au probleme fin de base
      Probleme_base& pbF = ref_cast(Probleme_base, ob2);
      le_pb_ext_ = pbF;
      Probleme_base& pbG = ref_cast(Probleme_base, ob3);
      le_pb_courant_ = pbG;
      //nature des conditions aux limites
      //dirichlet_ou_neumann = nature;
      //reference au champ inconnue
      REF(Champ_base) rch;
      rch =  pbF.get_champ(nom_inco);
      if (sub_type(Champ_Inc_base,rch.valeur()))
        {
          l_inconnueF=ref_cast(Champ_Inc_base, rch.valeur()) ;
          fixer_nb_comp(rch.valeur().nb_comp());
        }
      else
        {
          Cerr << pbF.le_nom() << "does not have unknown field name" << nom_inco << finl;
          exit();
        }
      rch = pbG.get_champ(nom_inco);
      if (sub_type(Champ_Inc_base,rch.valeur()))
        {
          l_inconnue=ref_cast(Champ_Inc_base, rch.valeur()) ;
          fixer_nb_comp(rch.valeur().nb_comp());
        }
      else
        {
          Cerr << pbG.le_nom() << "does not have unknown field name" << nom_inco << finl;
          exit();
        }
    }

}




// Description:
//    Non code
// Precondition:
// Parametre: Champ_front_base& ch
//    Signification:
//    Valeurs par defaut:
//    Contraintes: reference constante
//    Acces: NON ACCEDE
// Retour: Champ_front_base&
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Champ_front_base& Champ_front_zoom::affecter_(const Champ_front_base& ch)
{
  return *this;
}


int Champ_front_zoom::initialiser(double temps, const Champ_Inc_base& inco)
{

  Ch_front_var_instationnaire_dep::initialiser(temps, inco);

  // Dimensionnement des tableaux
  int nb_faces=frontiere_dis().frontiere().nb_faces();

  // Initialiser valeurs a des valeurs raisonnables mais
  // independantes de l'exterieur !

  int ndeb=frontiere_dis().frontiere().num_premiere_face();
  for (int fac_front=0; fac_front<nb_faces; fac_front++)
    {
      int fac_glob=fac_front+ndeb;
      valeurs()(fac_front,0)=l_inconnue->valeurs()(fac_glob);
    }
  return 1;
}

// Description:
//    Mise a jour en temps du champ
//    Cette methode effectue seulement le calcul de la derivee en
//    temps de la CL grace a la methode GPoint : le champ est deja mis
//    a jour par ailleurs (a proscrire!)
// Precondition:
// Parametre: double
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
void Champ_front_zoom::mettre_a_jour(double temps)
{
  DoubleTab& tab=valeurs_au_temps(temps);
  tab.echange_espace_virtuel();
}


// Description:
//    Renvoie le probleme fin base sur lequel on va prolonger.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Nom&
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Probleme_base& Champ_front_zoom::le_pb_fin()
{
  return le_pb_ext_.valeur();
}



// Description:
//    Renvoie le probleme exterieur a celui concerne par ce champ_front.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Nom&
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Probleme_base& Champ_front_zoom::le_pb_exterieur()
{
  return le_pb_ext_.valeur();
}


// Description:
//    Renvoie le probleme  base contenant ce champ_front.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Nom&
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Probleme_base& Champ_front_zoom::le_pb_courant()
{
  return le_pb_courant_.valeur();
}


// Description:
//    Renvoie le probleme Multigrille sur lequel s applique le champ.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Nom&
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Pb_MG& Champ_front_zoom::le_pb_MG()
{
  return le_pb_MG_.valeur();
}



/**
 * Renvoie l'inconnue du probleme exterieur
 */
const Champ_Inc_base& Champ_front_zoom::inconnue() const
{
  return l_inconnueF.valeur();
}



// Description:
//    Renvoie la zone discretisee associee a l'equation
//    qui porte le champ inconnue dont on prend la trace.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Zone_dis_base&
//    Signification: la zone discretisee associee a l'equation
//                   qui porte le champ inconnue dont on prend la trace
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Zone_dis_base& Champ_front_zoom::zone_dis() const
{
  return inconnue().zone_dis_base();
}


/**
 * Renvoie le milieu du probleme exterieur
 */
const Milieu_base& Champ_front_zoom::milieu() const
{
  return le_pb_ext_->milieu();
}


// Description:
//    Renvoie l'equation associee au probleme exterieur.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Equation_base&
//    Signification: l'equation associee a l'inconnue dont on prend
//                   la trace
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Equation_base& Champ_front_zoom::equation() const
{
  return inconnue().equation();
}
// Description:
//    Renvoie la zone des conditions au limites discretisees
//    portee par l'equation qui porte le champ inconnue
//    dont on prend la trace
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Zone_Cl_dis_base&
//    Signification: la zone des conditions au limites discretisees
//                   portee par l'equation qui porte le champ inconnue
//                   dont on prend la trace
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Zone_Cl_dis_base& Champ_front_zoom::zone_Cl_dis() const
{
  return equation().zone_Cl_dis().valeur();
}

// Description:
//    Renvoie la frontiere discretisee correspondante
//    au domaine sur lequel prend la trace.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Frontiere_dis_base&
//    Signification: frontiere discretisee correspondante
//                   au domaine sur lequel prend la trace
//    Contraintes: reference constante
// Exception: frontiere du nom specifie introuvable
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Frontiere_dis_base& Champ_front_zoom::front_dis_exterieure() const
{
  //Cerr << "Dans Champ_front_zoom frontiere dis " << finl;
  const Zone_Cl_dis_base& zcl=zone_Cl_dis();
  int n=zcl.nb_cond_lim();
  for(int i=0; i<n ; i++)
    {
      const Frontiere_dis_base& fr_dis=zcl.les_conditions_limites(i).frontiere_dis();
      if (fr_dis.le_nom() == nom_bord)
        return fr_dis;
    }
  Cerr << "Error in Champ_front_calc : we did not find a boundary with name : "
       << nom_bord << finl;
  exit();
  return zcl.les_conditions_limites(0).frontiere_dis();
}
