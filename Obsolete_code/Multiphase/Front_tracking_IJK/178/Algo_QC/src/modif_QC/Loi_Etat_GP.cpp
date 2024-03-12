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
// File:        Loi_Etat_GP.cpp
// Directory:   $TRIO_U_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/25
//
//////////////////////////////////////////////////////////////////////////////

#include <Loi_Etat_GP.h>
#include <Motcle.h>
#include <Fluide_Quasi_Compressible.h>
#include <Champ_Uniforme.h>
#include <Zone_VF.h>
#include <Champ_Fonc_Tabule.h>
#include <Equation_base.h>
#include <Probleme_base.h>

Implemente_instanciable_sans_constructeur(Loi_Etat_GP,"Loi_Etat_Gaz_Parfait",Loi_Etat_base);


Loi_Etat_GP::Loi_Etat_GP()
{
  Cp_ = -1;
  R_ = -1;
}

// Description:
//    Imprime la loi sur un flot de sortie.
// Precondition:
// Parametre: Sortie& os
//    Signification: le flot de sortie pour l'impression
//    Valeurs par defaut:
//    Contraintes:
//    Acces: sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord: le flot de sortie est modifie
// Postcondition: la methode ne modifie pas l'objet
Sortie& Loi_Etat_GP::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

// Description:
//    Lecture d'une loi sur un flot d'entree.
// Precondition:
// Parametre: Entree& is
//    Signification: le flot d'entree pour la lecture des parametres
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception: accolade ouvrante attendue
// Effets de bord:
// Postcondition: l'objet est construit avec les parametres lus
Entree& Loi_Etat_GP::readOn(Entree& is)
{
  double Cv_ = -1;
  double gamma_ = -1;

  Motcle accferme="}";
  Motcle accouverte="{";

  Motcle motlu;
  is >> motlu;
  Cerr<<"Lecture de la loi d'etat Gaz Parfait"<<finl;
  if (motlu != accouverte)
    {
      Cerr<<" On attendait "<<accouverte<<" au lieu de "<<motlu<<finl;
      abort();
    }
  Motcles les_mots(7);
  {
    les_mots[0] = "Cp";
    les_mots[1] = "capacite_calorifique_pression_constante";
    les_mots[2] = "Cv";
    les_mots[3] = "capacite_calorifique_volume_constant";
    les_mots[4] = "gamma";
    les_mots[5] = "Prandtl";
    les_mots[6] = "rho_constant_pour_debug";
  }
  is >> motlu;
  while(motlu != accferme )
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :
        case 1 :
          {
            is>>Cp_;
            break;
          }
        case 2 :
        case 3 :
          {
            is>>Cv_;
            break;
          }
        case 4 :
          {
            is>>gamma_;
            if (gamma_<0)
              {
                gamma_ = -gamma_;
                debug=1;
              }
            else
              debug=0;
            break;
          }
        case 5 :
          {
            is>>Pr_;
            break;
          }
        case 6 :
          {
            is>>rho_constant_pour_debug_;
            break;
          }
        default :
          {
            Cerr<<"Une loi d'etat "<<que_suis_je()<<" n'a pas la propriete "<<motlu<<finl;
            Cerr<<"On attendait un mot dans :"<<finl<<les_mots<<finl;
            abort();
          }
        }
      is >> motlu;
    }

  if (Pr_==-1)
    {
      Cerr<<"ERREUR : on attendait la definition du nombre de Prandtl (constante)"<<finl;
      abort();
    }
  if (Cp_==-1)
    {
      Cerr<<"ERREUR : on attendait la definition du Cp (constante en gaz parfaits)"<<finl;
      abort();
    }
  if (Cv_!=-1)
    {
      R_ = Cp_ - Cv_;
    }
  else if (gamma_!=-1)
    {
      R_ = Cp_ *(1.-1./gamma_);
    }
  else
    {
      Cerr<<"ERREUR : on attendait la definition du Cv (constante en gaz parfaits)"<<finl;
      Cerr<<"ou de la constante gamma (constante en gaz parfaits)"<<finl;
      abort();
    }
  return is;
}

// Description:
//    Renvoie le type de fluide associe.
// Precondition:
// Parametre: Sortie& os
//    Signification: le flot de sortie pour l'impression
//    Valeurs par defaut:
//    Contraintes:
//    Acces: sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord: le flot de sortie est modifie
// Postcondition: la methode ne modifie pas l'objet
const Nom Loi_Etat_GP::type_fluide() const
{
  return "Gaz_Parfait";
}

// Description:
//    Associe le fluide a la loi d'etat
// Precondition:
// Parametre: Fluide_Quasi_Compressible& fl
//    Signification: le fluide associe
//    Valeurs par defaut:
//    Contraintes: reference constante
//    Acces: lecture
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Loi_Etat_GP::associer_fluide(const Fluide_Quasi_Compressible& fl)
{
  le_fluide = fl;
  le_fluide->set_Cp(Cp_);
}


void Loi_Etat_GP::initialiser_inco_ch()
{
  Cerr << "Loi_Etat_GP::initialiser_inco_ch()" << finl;
  Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,le_fluide.valeur());//modif YB 28/08/09
  const double Pth = fluide_QC.pression_th();//modif YB 28/08/09

  DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  int i,ntot=tab_rho.size_totale();
  DoubleTab& tab_rhoI =  le_fluide->inco_chaleur().valeurs();
  if (1)
    //if (le_fluide->inco_chaleur()->equation().probleme().reprise_effectuee()==0)
    {

      tab_rho_n.copy(tab_rho);
      tab_rho_np1.copy(tab_rho);

      le_fluide->ch_temperature().valeurs()=tab_rhoI;
      for (i=0 ; i<ntot ; i++)
        {
          tab_rho_np1[i]=Pth/(R_*tab_rhoI[i]);
        }

      tab_rho= tab_rho_np1;
      tab_rho_n=tab_rho;
      tab_rhoI=tab_rho_np1;
    }
  else
    {
      tab_rho=tab_rhoI;
      tab_rho_n=tab_rho;
      tab_rho_np1=tab_rho;
    }
  Cerr << "5 fin Loi_Etat_GP::initialiser_inco_ch()" << finl;
  mettre_a_jour(0.);
  tab_rho.echange_espace_virtuel();
  tab_rho_np1.echange_espace_virtuel();
}


// Description:
//    Initialise la loi d'etat : calcul Pth
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
void Loi_Etat_GP::initialiser()
{
  const DoubleTab& tab_rho = le_fluide->inco_chaleur().valeurs();			//modif YB 28/08/09

  Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,le_fluide.valeur());//modif YB 28/08/09
  const double Pth = fluide_QC.pression_th();//modif YB 28/08/09

  // const DoubleTab& tab_Temp = le_fluide->temperature();					modif YB 28/08/09
  int i;
  //int n = tab_Temp.dimension(0);
  int ntot=tab_rho.size_totale();

  DoubleTab& tab_T = temperature_.valeurs();
  assert(tab_T.size()==tab_rho.dimension(0));
//tab_T.resize(n);

  tab_rho_n=tab_rho;
//tab_rho_np1=tab_Temp;
//  tab_rho_n.resize(n);
// tab_rho_np1.resize(n);
  for (i=0 ; i<ntot ; i++)
    {
      tab_rho_n[i] = tab_rho[i];
      tab_T[i] = Pth/(R_*tab_rho[i]);
    }
//   Cerr<<"Loi_Etat_GP::initialiser tab_Temp="<<tab_Temp[0]<<"  tab_T="<<tab_T[0]<<finl;
}

// Description:
//    Remplit le tableau de la temperature : T=temp+273.15
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
void Loi_Etat_GP::remplir_T()
{
  //const DoubleTab& tab_rho = le_fluide->inco_chaleur().valeurs();
  int i;
  //  int n = tab_Temp.dimension(0);
  int ntot=tab_rho_np1.size_totale();//modif YB 28/08/09
  DoubleTab& tab_T = temperature_.valeurs();
  Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,le_fluide.valeur());//modif YB 28/08/09
  const double Pth = fluide_QC.pression_th();//modif YB 28/08/09

  for (i=0 ; i<ntot ; i++)
    {
      tab_T[i] = Pth/(R_*tab_rho_np1[i]);		//modif YB 28/08/09
    }
//   Cerr<<"---Loi_Etat_GP::remplir_T present="<<le_fluide->inco_chaleur().valeurs()(0)<<" passe="<<le_fluide->inco_chaleur().passe()(0)<<" futur="<<le_fluide->inco_chaleur().futur()(0)<<"  local T="<<tab_T(0)<<finl;
//Cerr<<"YB - 'Loi_Etat_GP' - Temperature = "<<tab_T[0]<<" "<<tab_T[ntot-1]<<finl;
//Cerr<<"YB - 'Loi_Etat_GP' - masse volumique inco_ch = "<<tab_rho[0]<<" "<<tab_rho[ntot-1]<<finl;
//Cerr<<"YB - 'Loi_Etat_GP' - masse volumique n = "<<tab_rho_n[0]<<" "<<tab_rho_n[ntot-1]<<finl;
//Cerr<<"YB - 'Loi_Etat_GP' - masse volumique  np1 = "<<tab_rho_np1[0]<<" "<<tab_rho_np1[ntot-1]<<finl;

}

// Description:
//    Calcule le Cp
//    NE FAIT RIEN : le Cp est constant
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
void Loi_Etat_GP::calculer_Cp()
{
}

// Description:
//    Calcule la conductivite
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
void Loi_Etat_GP::calculer_lambda()
{
  const Champ_Don&  mu          = le_fluide->viscosite_dynamique();
  const DoubleTab&  tab_mu      = mu.valeurs();
  Champ_Don&        lambda      = le_fluide->conductivite();
  DoubleTab&        tab_lambda  = lambda.valeurs();


  int i, n=tab_lambda.size();

  //La conductivite est soit un champ uniforme soit calculee a partir de la viscosite dynamique et du Pr
  //  Cerr << " Cp = " << Cp_ << finl;
  if (sub_type(Champ_Fonc_Tabule,lambda.valeur()))
    {
      lambda.valeur().mettre_a_jour(temperature_.valeur().temps());
      return;
    }
  if (!sub_type(Champ_Uniforme,lambda.valeur()))
    {
      //    DoubleTab& tab_T = temperature_.valeurs();
      if (sub_type(Champ_Uniforme,mu.valeur()))
        {
          double mu0 = tab_mu(0,0);
          for (i=0 ; i<n ; i++)
            {
              tab_lambda[i] = mu0 * Cp_ / Pr_;
            }
        }
      else
        {
          double rap=Cp_/Pr_;
          for (i=0 ; i<n ; i++)
            {
              tab_lambda(i) = tab_mu(i) * rap;
              //Corrections pour calculer lambda en codant la loi de Sutherland

              //tab_lambda[i]=3.406e-2*((400.+259.2)/(tab_T[i]+259.2))*pow(tab_T[i]/400.,1.5);
              //tab_lambda[i]=0.13883215*((1300.+6325.5108)/(tab_T[i]+6325.5108))*pow(tab_T[i]/1300.,1.5);
            }
        }
    }
  else
    {
      if (sub_type(Champ_Uniforme,mu.valeur()))
        {
          tab_lambda(0,0) = mu(0,0) * Cp_ / Pr_;
        }
      else
        {
          Cerr << finl ;
          Cerr << "If lambda is of type Champ_Uniform, mu must also be of type Champ_Uniforme !"<< finl;
          Cerr << "If needed, you can use a Champ_Fonc_Fonction or a Champ_Fonc_Tabule with a constant for mu."<< finl ;
          Cerr << finl ;
          exit();
        }
    }


  tab_lambda.echange_espace_virtuel();
  //Cerr<<"YB - 'Loi_Etat_GP' - lambda = "<<min(tab_lambda)<<" "<<max(tab_lambda)<<finl;
}

// Description:
//    Calcule la diffusivite
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
void Loi_Etat_GP::calculer_alpha()
{
  const Champ_Don& lambda = le_fluide->conductivite();
  const DoubleTab& tab_lambda = lambda.valeurs();
  Champ_Don& alpha=le_fluide->diffusivite();
  DoubleTab& tab_alpha = le_fluide->diffusivite().valeurs();
  const DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  int isVDF=0;
  if (alpha.valeur().que_suis_je()=="Champ_Fonc_P0_VDF") isVDF=1;
  int i, n=tab_alpha.size();
  if (isVDF)
    {
      if (sub_type(Champ_Uniforme,lambda.valeur()))
        {
          double lambda0 = tab_lambda(0,0);
          for (i=0 ; i<n ; i++)
            {
              tab_alpha[i] = lambda0 / (tab_rho[i]*Cp_);
            }
        }
      else
        {
          for (i=0 ; i<n ; i++)
            {
              tab_alpha[i] = tab_lambda(i) / (tab_rho[i]*Cp_);
            }
        }
    }
  else
    {
      //const IntTab& elem_faces=ref_cast(Zone_VF,ref_cast(Champ_Fonc_P0_VEF,alpha.valeur()).zone_dis_base()).elem_faces();
      const IntTab& elem_faces=ref_cast(Zone_VF,le_fluide->vitesse().zone_dis_base()).elem_faces();
      double rhoelem;
      int nfe=elem_faces.dimension(1),face;
      if (sub_type(Champ_Uniforme,lambda.valeur()))
        {
          double lambda0 = tab_lambda(0,0);
          for (i=0 ; i<n ; i++)
            {
              rhoelem=0;
              for (face=0; face<nfe; face++) rhoelem+=tab_rho(elem_faces(i,face));
              rhoelem/=nfe;
              tab_alpha[i] = lambda0 / (rhoelem*Cp_);
            }
        }
      else
        {
          for (i=0 ; i<n ; i++)
            {
              rhoelem=0;
              for (face=0; face<nfe; face++) rhoelem+=tab_rho(elem_faces(i,face));
              rhoelem/=nfe;
              tab_alpha[i] = tab_lambda(i) / (rhoelem*Cp_);
            }
        }
    }
  tab_alpha.echange_espace_virtuel();
  //Cerr<<"YB - 'Loi_Etat_GP' - alpha = "<<min(tab_alpha)<<" "<<max(tab_alpha)<<finl;
}

// Description:
//    Recalcule la masse volumique
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
double Loi_Etat_GP::calculer_masse_volumique(double P, double T) const
{
  if (rho_constant_pour_debug_.non_nul())
    {
      return rho_constant_pour_debug_(0,0);
    }
  if (inf_ou_egal(T,0))
    {
      Cerr << finl << "Temperature T must be defined in Kelvin." << finl;
      Cerr << "Check your data file." << finl;
      exit();
    }
  return P/(R_*T);
}

// Description:
//    Calcule la pression avec la temperature et la masse volumique
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
double Loi_Etat_GP::inverser_Pth(double T, double rho)
{
  return rho*R_*T;
}
void Loi_Etat_GP::calculer_masse_volumique()
{
  Loi_Etat_base::calculer_masse_volumique();
}
