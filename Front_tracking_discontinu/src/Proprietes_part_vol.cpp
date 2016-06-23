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
// File:        Proprietes_part_vol.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Proprietes_part_vol.h>
#include <EFichier.h>
#include <Motcle.h>
#include <DoubleVect.h>
#include <ArrOfIntFT.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Proprietes_part_vol,"Proprietes_part_vol",Objet_U);

Proprietes_part_vol::Proprietes_part_vol()
{
  int dim=Objet_U::dimension;
  nb_particules_ = 0;
  vitesse_p_.resize(0,dim);
  deltat_v_.resize(0,dim);
  temperature_p_.resize(0,1);
  masse_volumique_p_.resize(0,1);
  diametre_p_.resize(0,1);
}

Sortie& Proprietes_part_vol::printOn(Sortie& os) const
{
  assert(0);
  exit();
  return os;
}

// Description:
//  Lecture des proprietes sur un flot d'entree.
//  Format de lecture :
// -cas 1 lecture dans un fichier (motcle fichier)
// -cas 2 specification des valeurs uniformes (motcle distribution)
//  Appel a lire_distribution()

Entree& Proprietes_part_vol::readOn(Entree& is)
{
  Cerr << "Reading of data for the particles properties" << finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  const LIST(Nom)& params_lu = param. get_list_mots_lus();
  if (params_lu.contient(Motcle("fichier")) && params_lu.contient(Motcle("distribution")))
    {
      Cerr<<"Both keyword fichier and distribution have been readen."<<finl;
      Cerr<<"Only one of these keywords can be selected."<<finl;
      Cerr<<"Please, see the User_Manual_TRUST for more details."<<finl;
      exit();
    }
  return is;
}

void Proprietes_part_vol::set_param(Param& param)
{
  param.ajouter_non_std("fichier",(this));
  param.ajouter_non_std("distribution",(this));
}

int Proprietes_part_vol::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="fichier")
    {
      Nom nomfic;
      is >> nomfic;
      EFichier fic;
      Cerr << "Reading the file " << nomfic << finl;
      if (!fic.ouvrir(nomfic))
        {
          Cerr << " Error while opening the file "<< nomfic << finl;
          exit();
        }
      fic >> vitesse_p_ >> temperature_p_ >> masse_volumique_p_>> diametre_p_;
      nb_particules_ = vitesse_p_.dimension(0);
      return 1;
    }
  else if (mot=="distribution")
    {
      lire_distribution(is);
      return 1;
    }
  else
    {
      Cerr << mot << " is not a keyword understood by " << que_suis_je() << " in lire_motcle_non_standard"<< finl;
      exit();
      return -1;
    }
  return 1;
}

//Dimensionnement de delta_v_ et initialisation a 0
//Dimensionnement de volume_p_ et initialisation a :
// -cas 2D : pi*dp*dp/4
// -cas 3D : pi*dp*dp*dp/6
void Proprietes_part_vol::completer()
{
  //Dimensionnement de deltat_v_
  //Valeurs initialisees a zero car pas de sens en premiere iteration
  //pour le terme source de masse ajoutee
  int dim =Objet_U::dimension;
  deltat_v_.resize(nb_particules_,dim);
  deltat_v_=0.;

  //Calcul du volume des particules
  double pi = M_PI;
  double dp;
  volume_p_.resize(nb_particules_,1);

  if (dim==2)
    {
      for (int i=0; i<nb_particules_; i++)
        {
          dp = diametre_p_(i,0);
          volume_p_(i,0) = (pi/4.)*dp*dp;
        }
    }
  else if (dim==3)
    {
      for (int i=0; i<nb_particules_; i++)
        {
          dp = diametre_p_(i,0);
          volume_p_(i,0) = (pi/6.)*dp*dp*dp;
        }
    }
  else
    {
      Cerr<<"Method Proprietes_part_vol::completer()"<<finl;
      Cerr<<"The case dimension "<<dim<<" is not treated."<<finl;
      exit();
    }
}

//Methode qui permet a un objet exterieur de modifier la valeur de nb_particules_
void Proprietes_part_vol::fixer_nb_particules(const int& nb_part)
{
  nb_particules_ = nb_part;
}

//Nettoyage par suppression des particules qui ne sont plus utilisees
//La particule i est supprimee si som_utilise(i)=0
void Proprietes_part_vol::nettoyer(const ArrOfInt& som_utilises)
{
  // Compteur de sommets apres suppression
  int n = 0;
  for (int i = 0; i<nb_particules_; i++)
    {
      if (som_utilises[i])
        {
          for (int dim=0; dim<dimension; dim++)
            {
              vitesse_p_(n,dim) = vitesse_p_(i,dim);
              deltat_v_(n,dim) = deltat_v_(i,dim);
            }
          temperature_p_(n,0) = temperature_p_(i,0);
          masse_volumique_p_(n,0) = masse_volumique_p_(i,0);
          diametre_p_(n,0) = diametre_p_(i,0);
          volume_p_(n,0) = volume_p_(i,0);
          n++;
        }
    }

  vitesse_p_.resize(n,dimension);
  temperature_p_.resize(n,1);
  masse_volumique_p_.resize(n,1);
  diametre_p_.resize(n,1);
  volume_p_.resize(n,1);
  nb_particules_ = n;
}

//Ajout des proprietes contenues dans proprietes_tmp a celles de l objet courant
void Proprietes_part_vol::ajouter_proprietes(const Proprietes_part_vol& proprietes_tmp)
{
  const DoubleTab& vitesse_par = proprietes_tmp.vitesse_particules();
  const DoubleTab& delta_vit = proprietes_tmp.delta_v();
  const DoubleTab& temperature_par = proprietes_tmp.temperature_particules();
  const DoubleTab& masse_vol_par = proprietes_tmp.masse_vol_particules();
  const DoubleTab& diametre_par = proprietes_tmp.diametre_particules();
  const DoubleTab& volume_par = proprietes_tmp.volume_particules();

  const int nb_part_tmp = proprietes_tmp.nb_particules();
  const int nb_part_tot = nb_particules_ + nb_part_tmp;
  int par, par_tmp;
  int k;

  int dim1 = vitesse_p_.dimension(1);
  vitesse_p_.resize(nb_part_tot,dim1);
  deltat_v_.resize(nb_part_tot,dim1);
  temperature_p_.resize(nb_part_tot,1);
  masse_volumique_p_.resize(nb_part_tot,1);
  diametre_p_.resize(nb_part_tot,1);
  volume_p_.resize(nb_part_tot,1);

  for (par_tmp=0 ; par_tmp<nb_part_tmp ; par_tmp++)
    {
      par = par_tmp + nb_particules_;
      for (k=0 ; k<dim1 ; k++)
        {
          vitesse_p_(par,k) = vitesse_par(par_tmp,k);
          deltat_v_(par,k) = delta_vit(par_tmp,k);
        }

      temperature_p_(par,0) = temperature_par(par_tmp,0);
      masse_volumique_p_(par,0) = masse_vol_par(par_tmp,0);
      diametre_p_(par,0) = diametre_par(par_tmp,0);
      volume_p_(par,0) = volume_par(par_tmp,0);
    }

  nb_particules_ = nb_part_tot;
}


//  Lecture des valeurs (uniformes) de la vitesse, de la temperature,
//  de la densite, du diametre
//  Format de lecture :
//  "distribution"   "nb_particules"    val_nb_particules
//                   "vitesse"          val_1 val_2 [val_3]
//                   "temperaturee      val_temp
//                   "masse_volumique"  val_rho
//                   "diametre"         val_diam
void Proprietes_part_vol::lire_distribution(Entree& is)
{
  int dim = Objet_U::dimension;
  ArrOfDouble vitesse(dim);
  double temperature;
  double masse_volumique;
  double diametre;

  Param param_distrib(que_suis_je());
  param_distrib.ajouter("nb_particules",&nb_particules_,Param::REQUIRED);
  param_distrib.ajouter_arr_size_predefinie("vitesse",&vitesse,Param::REQUIRED);
  param_distrib.ajouter("temperature",&temperature,Param::REQUIRED);
  param_distrib.ajouter("masse_volumique",&masse_volumique,Param::REQUIRED);
  param_distrib.ajouter("diametre",&diametre,Param::REQUIRED);
  param_distrib.lire_avec_accolades_depuis(is);

  vitesse_p_.resize(nb_particules_,dim);
  temperature_p_.resize(nb_particules_,1);
  masse_volumique_p_.resize(nb_particules_,1);
  diametre_p_.resize(nb_particules_,1);

  for (int j=0; j<dim; j++)
    {
      for (int i=0; i<nb_particules_; i++)
        vitesse_p_(i,j) = vitesse(j);
    }

  for (int i=0; i<nb_particules_; i++)
    {
      temperature_p_(i,0) = temperature;
      masse_volumique_p_(i,0) = masse_volumique;
      diametre_p_(i,0) = diametre;
    }
}


