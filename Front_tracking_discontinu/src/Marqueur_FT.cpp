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
// File:        Marqueur_FT.cpp
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////
#include <Marqueur_FT.h>
#include <Probleme_base.h>
#include <Zone_VF.h>
#include <Transport_Interfaces_FT_Disc.h>

Implemente_instanciable(Marqueur_FT,"Marqueur_FT",Marqueur_Lagrange_base);

Sortie& Marqueur_FT::printOn(Sortie& os) const
{

  return os;
}

//Lecture d un marqueur
// ensemble_points         : On lit les informations necessaires pour localiser les particules
//                          -soit lues dans un fichier
//                          -soit generees dans des sous zones
// t_debut_integration  : temps de debut d integration des trajectoires (fixe a t_init par defaut)

Entree& Marqueur_FT::readOn(Entree& is)
{
  //temporaire
  /*
    Motcles les_mots(2);
    {
    les_mots[0] = "ensemble_points";
    les_mots[1] = "t_debut_integration";
    }

    Motcle motlu, accolade_fermee="}", accolade_ouverte="{";
    is >> motlu;
    if(motlu!=accolade_ouverte)
    {
    Cerr << "On attendait une { a la lecture d'une " << que_suis_je() << finl;
    Cerr << "et non : " << motlu << finl;
    exit();
    }
    is >> motlu;

    while (motlu != accolade_fermee)
    {
    int rang=les_mots.search(motlu);
    switch(rang){
    case 0 :
    {
    is>>ensemble_points_;
    break;
    }
    case 1 :
    {
    is>>t_debut_integr_;
    break;
    }
    }
    is >> motlu;
    }
  */
  //temporaire
  exit();
  return is;
}

//Surcharge de Marqueur_Lagrange_base::discretiser()
//Les coordonnees des points a suivre sont :
//  -soit lues dans un fichier
//  -soit generees dans des sous zones

void Marqueur_FT::discretiser(const Probleme_base& pb, const  Discretisation_base& dis)
{
  //temporaire
  /*
    Marqueur_Lagrange_base::discretiser(pb,dis);
    Maillage_FT_Disc& maillage = ref_cast(Maillage_FT_Disc,ensemble_points());
    const DoubleTab & sommets = maillage.sommets() ;

    DoubleTab sommets_tmp;
    //Premier cas : les sommets ont ete lus dans un fichier
    //On copie sommets dans un tableau temporaire sommets_tmp
    if (maillage.nb_marqs_par_sz().size()==0)
    {
    int dim0=sommets.dimension(0);
    int dim1=sommets.dimension(1);
    sommets_tmp.resize(dim0,dim1);

    for (int i=0;i<dim0;i++)
    for (int j=0;j<dim1;j++)
    sommets_tmp(i,j)=sommets(i,j);
    }
    else
    //Deuxieme cas : les sommets sont crees dans des sous zones
    //On remplit le tableau temporaire
    maillage.generer_marqueurs_sz(sommets_tmp);

    maillage.construire_points(sommets_tmp);
    calculer_valeurs_champs();
  */

  //temporaire
  exit();
}

//Calcul des valeurs du champ postraitable densite_particules_
//Estimation du nombre de partcules par maille
void Marqueur_FT::calculer_valeurs_champs()
{
  //temporaire
  /*
    DoubleTab& val_densite = densite_particules_.valeurs();
    const Maillage_FT_Disc& ens_points = ref_cast(Maillage_FT_Disc,ensemble_points());
    const ArrOfInt & sommet_elem = ens_points.sommet_elem();
    val_densite = 0.;
    int nb_som = ens_points.nb_sommets();
    int elem;

    for (int som=0; som<nb_som; som++)
    {
    elem = sommet_elem[som];
    if (elem!=-1)
    val_densite(elem) += 1.;
    else
    {
    if (je_suis_maitre())
    Cerr<<"Une particule n est pas localisee dans le domaine"<<finl;
    exit();
    }
    }
  */

  //temporaire
  exit();
}
