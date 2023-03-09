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
// File:        Loi_Etat_base.cpp
// Directory:   $TRIO_U_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/29
//
//////////////////////////////////////////////////////////////////////////////

#include <Loi_Etat_base.h>
#include <Champ_Uniforme.h>
#include <Fluide_Quasi_Compressible.h>
#include <Champ_Fonc_Tabule.h>
#include <Domaine_VF.h>
#include <Debog.h>

Implemente_base_sans_constructeur(Loi_Etat_base,"Loi_Etat_base",Objet_U);

Loi_Etat_base::Loi_Etat_base()
{
  /*
    Noms& nom=champs_compris_.liste_noms_compris();
    nom.dimensionner(1);
    nom[0]="temperature";
  */
  Pr_=-1.;
  debug=-1;
}

/*! @brief Imprime la loi sur un flot de sortie.
 *
 * @param (Sortie& os) le flot de sortie pour l'impression
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Loi_Etat_base::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

/*! @brief Lecture d'une loi sur un flot d'entree.
 *
 * @param (Entree& is) le flot d'entree pour la lecture des parametres
 * @return (Entree&) le flot d'entree modifie
 * @throws accolade ouvrante attendue
 */
Entree& Loi_Etat_base::readOn(Entree& is)
{
  return is;
}

/*! @brief Associe le fluide a la loi d'etat
 *
 * @param (Fluide_Quasi_Compressible& fl) le fluide associe
 */
void Loi_Etat_base::associer_fluide(const Fluide_Quasi_Compressible& fl)
{
  le_fluide = fl;
}

/*! @brief Renvoie le champ de le temperature
 *
 */
const Champ_Don& Loi_Etat_base::ch_temperature() const
{
  return temperature_;
}
Champ_Don& Loi_Etat_base::ch_temperature()
{
  return temperature_;
}

/*! @brief Initialise l'inconnue de l'equation de chaleur : ne fai rien
 *
 */
void Loi_Etat_base::initialiser_inco_ch()
{
  DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();

  tab_rho_n=tab_rho;
  tab_rho_np1=tab_rho;
  calculer_masse_volumique();
  mettre_a_jour(0.);
  tab_rho.echange_espace_virtuel();
  tab_rho_np1.echange_espace_virtuel();

}

/*! @brief Prepare le fluide au calcul.
 *
 */
void Loi_Etat_base::preparer_calcul()
{
  remplir_T();

  calculer_masse_volumique();
  mettre_a_jour(0.);
}

/*! @brief Met a jour la loi d'etat
 *
 * @param (double temps) le temps de calcul
 */
void Loi_Etat_base::mettre_a_jour(double temps)
{
  //remplissage de rho avec rho(n+1)
  DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  int i, n=tab_rho.size_totale();
  for (i=0 ; i<n ; i++)
    {
      // tab_rho_nm1[i]=tab_rho_n[i];
      tab_rho_n[i] = tab_rho_np1[i];
      tab_rho[i] = tab_rho_np1[i];
    }
}

/*! @brief Calcule la viscosite
 *
 */
void Loi_Etat_base::calculer_mu()
{
  Champ_Don& mu = le_fluide->viscosite_dynamique();
  if (!sub_type(Champ_Uniforme,mu.valeur()))
    {
      if (sub_type(Champ_Fonc_Tabule,mu.valeur()))
        {
#if 0
          const DoubleTab& T = le_fluide->get_champ("temperature_qc").valeurs();
          DoubleTab& tab=mu.valeurs();
          int n=tab.dimension(0);
          //Cerr<< "ll "<<tab.nb_dim()<<finl;exit();
          for (int i=0; i<n; i++)
            {
              double x=T(i);
              if (1)
                tab(i)=1.461e-6*x*sqrt(x)/(x+111.);

              else
                tab(i)=((((2.55137e-14)*x-5.7633e-11)*x+7.49077e-08)*x+4.98671e-07);


            }
#else
          mu.mettre_a_jour(temperature_.valeur().temps());
#endif
        }
      else
        {

          const DoubleTab& T = le_fluide->get_champ("temperature_qc").valeurs();
          DoubleTab& tab=mu.valeurs();
          int n=tab.dimension(0);
          //Cerr<< "ll "<<tab.nb_dim()<<finl;exit();
          for (int i=0; i<n; i++)
            {
              double val=T(i);
              tab(i)=1.461e-6*val*sqrt(val)/(val+111.);
            }
          tab.echange_espace_virtuel();
          /* Cerr<<"The viscosity field mu of type "<<mu.valeur().que_suis_je()<<" is not recognized.";
           exit();
           */
        }
    }
}

/*! @brief Calcule la conductivite
 *
 */
void Loi_Etat_base::calculer_lambda()
{
  const Champ_Don& mu = le_fluide->viscosite_dynamique();
  const DoubleTab& tab_Cp = le_fluide->capacite_calorifique().valeurs();
  const DoubleTab& tab_mu = mu.valeurs();
  Champ_Don& lambda = le_fluide->conductivite();
  DoubleTab& tab_lambda = lambda.valeurs();

  int i, n=tab_lambda.size();
  //La conductivite est soit un champ uniforme soit calculee a partir de la vis  //cosite dynamique et du Pr

  if (!sub_type(Champ_Uniforme,lambda.valeur()))
    {

      if (sub_type(Champ_Uniforme,mu.valeur()))
        {
          double mu0 = tab_mu(0,0);
          for (i=0 ; i<n ; i++)
            {
              tab_lambda[i] = mu0 * tab_Cp[i] / Pr_;
            }
        }
      else
        {
          for (i=0 ; i<n ; i++)
            {
              tab_lambda[i] = tab_mu[i] * tab_Cp[i] / Pr_;
            }
        }

    }
  tab_lambda.echange_espace_virtuel();
}

/*! @brief Calcule la viscosite cinematique
 *
 */
void Loi_Etat_base::calculer_nu()
{
  const Champ_Don& mu = le_fluide->viscosite_dynamique();
  const DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  const DoubleTab& tab_mu = mu.valeurs();
  Champ_Don& nu= le_fluide->viscosite_cinematique ();
  DoubleTab& tab_nu = nu.valeurs();
  //Cerr<<le_fluide->viscosite_cinematique().valeur().que_suis_je()<<finl;
  int isVDF=0;
  if (nu.valeur().que_suis_je()=="Champ_Fonc_P0_VDF") isVDF=1;
  int i, n=tab_nu.size();
  if (isVDF)
    {
      if (sub_type(Champ_Uniforme,mu.valeur()))
        {
          double mu0 = tab_mu(0,0);
          for (i=0 ; i<n ; i++)
            {
              tab_nu[i] = mu0 / tab_rho[i];
            }
        }
      else
        {
          for (i=0 ; i<n ; i++)
            {
              tab_nu[i] = tab_mu[i] / tab_rho[i];
            }
        }
    }
  else
    {
      //const IntTab& elem_faces=ref_cast(Domaine_VF,ref_cast(Champ_Fonc_P0_VEF,nu.valeur()).domaine_dis_base()).elem_faces();
      const IntTab& elem_faces=ref_cast(Domaine_VF,le_fluide->vitesse().domaine_dis_base()).elem_faces();
      double rhoelem;
      int nfe=elem_faces.dimension(1),face;
      if (sub_type(Champ_Uniforme,mu.valeur()))
        {
          double mu0 = tab_mu(0,0);
          for (i=0 ; i<n ; i++)
            {
              rhoelem=0;
              for (face=0; face<nfe; face++) rhoelem+=tab_rho(elem_faces(i,face));
              rhoelem/=nfe;
              tab_nu[i] = mu0 /rhoelem;
            }
        }
      else
        {
          for (i=0 ; i<n ; i++)
            {
              rhoelem=0;
              for (face=0; face<nfe; face++) rhoelem+=tab_rho(elem_faces(i,face));
              rhoelem/=nfe;
              tab_nu[i] = tab_mu[i] /rhoelem;
            }
        }

    }
  tab_nu.echange_espace_virtuel();
  Debog::verifier("Loi_Etat_base::calculer_nu tab_nu",tab_nu);
}

/*! @brief Calcule la diffusivite
 *
 */
void Loi_Etat_base::calculer_alpha()
{
  const Champ_Don& lambda = le_fluide->conductivite();
  const DoubleTab& tab_lambda = lambda.valeurs();
  const DoubleTab& tab_Cp = le_fluide->capacite_calorifique().valeurs();
  const DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  DoubleTab& tab_alpha = le_fluide->diffusivite().valeurs();

  int i, n=tab_alpha.size();
  if (sub_type(Champ_Uniforme,lambda.valeur()))
    {
      double lambda0 = tab_lambda(0,0);
      for (i=0 ; i<n ; i++)
        {
          tab_alpha[i] = lambda0 / (tab_rho[i]*tab_Cp[i]);
        }
    }
  else
    {
      for (i=0 ; i<n ; i++)
        {
          tab_alpha[i] = tab_lambda(i) / (tab_rho[i]*tab_Cp[i]);
        }
    }
  tab_alpha.echange_espace_virtuel();
  Debog::verifier("Loi_Etat_base::calculer_alpha alpha",tab_alpha);
}


/*! @brief Ne fait rien Surcharge dans Loi_Etat_Melange_GP
 *
 */
void Loi_Etat_base::calculer_mu_sur_Sc()
{

}
/*! @brief Recalcule la masse volumique
 *
 */
void Loi_Etat_base::calculer_masse_volumique()
{
  DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  //int som, n=tab_rho.size();


  Debog::verifier("Loi_Etat_base::calculer_masse_volumique, tab_rho 0",tab_rho);
  Debog::verifier("Loi_Etat_base::calculer_masse_volumique, debug 0",debug);

  tab_rho_np1 = le_fluide->inco_chaleur().valeurs();
  //for (som=0 ; som<n ; som++) {
  tab_rho=tab_rho_np1; 		//le_fluide.masse_volumique renvoie rho n+1		//[som] = (tab_rho_n[som] + tab_rho_np1[som])/2.;
  //}
  //Cerr<<"'Loi_Etat_base :: calculer_masse_volumique()'  ---> tab_rho ="<<min(tab_rho)<<"  "<<max(tab_rho)<<finl;
  tab_rho.echange_espace_virtuel();
  tab_rho_np1.echange_espace_virtuel();
  le_fluide->calculer_rho_face(tab_rho_np1); // l'argument de cette fonction est sans effet ! AT 04 09 09
  Debog::verifier("Loi_Etat_base::calculer_masse_volumique, tab_rho_n 1",tab_rho_n);
  Debog::verifier("Loi_Etat_base::calculer_masse_volumique, tab_rho_np1 1",tab_rho_np1);
  Debog::verifier("Loi_Etat_base::calculer_masse_volumique, tab_rho 1",tab_rho);
}

/*! @brief Cas gaz parfait : ne fait rien Cas gaz Reel : doit recalculer l'enthalpie a partir de la pression et la temperature
 *
 */
double Loi_Etat_base::calculer_H(double Pth_, double T_) const
{
  return T_;
}

double Loi_Etat_base::Drho_DP(double,double) const
{
  Cerr<<"Drho_DP doit etre code dans la classe fille "<<que_suis_je()<<" pour etre utilisee"<<finl;
  abort();
  return 0;
}
double Loi_Etat_base::Drho_DT(double,double) const
{
  Cerr<<"Drho_DT doit etre code dans la classe fille "<<que_suis_je()<<" pour etre utilisee"<<finl;
  abort();
  return 0;
}
double Loi_Etat_base::De_DP(double,double) const
{
  Cerr<<"De_DP doit etre code dans la classe fille "<<que_suis_je()<<" pour etre utilisee"<<finl;
  abort();
  return 0;
}
double Loi_Etat_base::De_DT(double,double) const
{
  Cerr<<"De_DT doit etre code dans la classe fille "<<que_suis_je()<<" pour etre utilisee"<<finl;
  abort();
  return 0;
}


void Loi_Etat_base::creer_champ(const Motcle& motlu)
{
}

const Champ_base& Loi_Etat_base::get_champ(const Motcle& nom) const
{
  return champs_compris_.get_champ(nom);
}
void Loi_Etat_base::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  if (opt==DESCRIPTION)
    Cerr<<"Loi_Etat_base : "<<champs_compris_.liste_noms_compris()<<finl;
  else
    nom.add(champs_compris_.liste_noms_compris());

}

void Loi_Etat_base::abortTimeStep()
{
  int i, n=tab_rho_n.size_totale();
  for (i=0 ; i<n ; i++)
    {
      tab_rho_np1[i] = tab_rho_n[i];
    }
}
