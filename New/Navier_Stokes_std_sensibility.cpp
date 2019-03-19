/****************************************************************************
* Copyright (c) 2019, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Navier_Stokes_std_sensibility.cpp
// Directory : $BALTIK_COUPLAGE_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////

#include <Navier_Stokes_std_sensibility.h>
#include <Probleme_base.h>
#include <Param.h>
#include <Debog.h>
#include <LecFicDiffuse.h>
#include <communications.h>
#include <Interprete.h>

Implemente_instanciable( Navier_Stokes_std_sensibility, "Navier_Stokes_standard_sensibility", Navier_Stokes_std) ;

Sortie& Navier_Stokes_std_sensibility::printOn( Sortie& os ) const
{
  return Navier_Stokes_std::printOn( os );
}

Entree& Navier_Stokes_std_sensibility::readOn( Entree& is )
{
  Navier_Stokes_std::readOn(is);

  return is;
}

void Navier_Stokes_std_sensibility::set_param(Param& param)
{

  Navier_Stokes_std::set_param(param);
  param.ajouter_non_std("state",(this),Param::REQUIRED);

}


int Navier_Stokes_std_sensibility::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="state")
    {
      int lu_info_evaluateur = 0;
      Cerr << "Reading and typing of the state : " << finl;
      Motcle motlu, accolade_fermee="}", accolade_ouverte="{";
      is >> motlu;
      if(motlu!=accolade_ouverte)
        {
          Cerr << "We expected a { while reading of " << que_suis_je() << finl;
          Cerr << "and not : " << motlu << finl;
          exit();
        }

      is >> motlu;
      Cerr<<"word read="<<motlu<<finl;
      if(motlu=="pb_champ_evaluateur")
        {
          is >> nom_pb_etat >> nom_inco_etat;
          lu_info_evaluateur = 1;
          if( nom_pb_etat  ==   accolade_fermee  || nom_pb_etat  ==  accolade_ouverte   ||  nom_inco_etat  ==   accolade_fermee  || nom_inco_etat  ==  accolade_ouverte )
            {
              Cerr<<"We expected the name of a problem and a fluid fild while reading of "<<motlu<< finl;
              Cerr << "and not : " <<accolade_ouverte << " or "<<  accolade_fermee << finl;
              exit();
            }
        }
      else
        {
          Cerr<<"Navier_Stokes_std_sensibility::lire_motcle_non_standard: keyword "<<motlu<<" is not recognized."<<finl;
          Cerr<<"The recognized keywords are :"<<"pb_champ_evaluateur "<<finl;
          exit();
        }
      is>>motlu;
      if(motlu != accolade_fermee)
        {
          Cerr << "We expected a } while reading of " << que_suis_je() << finl;
          Cerr << "and not : " << motlu << finl;
          exit();
        }
      if (lu_info_evaluateur!=1)
        {
          Cerr<<"Keyword pb_champ_evaluateur must be specified with associated data"<<finl;
          Cerr<<"when Navier_Stokes_std_sensibility is used."<<finl;
          exit();
        }
      associer_champ_evaluateur( nom_pb_etat, nom_inco_etat);
      return 1;
    }
  else
    return Navier_Stokes_std::lire_motcle_non_standard(mot, is);
}

//Initialisation de la reference au champ evaluateur (l_inconnue_etat)
void Navier_Stokes_std_sensibility::associer_champ_evaluateur(const Nom& un_nom_pb_etat,const Motcle& un_nom_inco_etat)
{
  Cerr <<"Navier_Stokes_std_sensibility: recoup state from "<<un_nom_pb_etat<< finl;
  if(probleme().le_nom() ==un_nom_pb_etat )
    {
      Cerr <<"Navier_Stokes_std_sensibility: we expect here the nom of the state problem and not the name of the sensibility problem"<< finl;
      exit();
    }
  Objet_U& ob= Interprete::objet(un_nom_pb_etat);
  REF(Probleme_base) pb;
  REF(Champ_base) rch;

  if(sub_type(Probleme_base,ob))
    {
      pb = ref_cast(Probleme_base,ob);
    }
  else
    {
      Cerr <<"No problem named "<< un_nom_pb_etat <<" has been found."<< finl;
      exit();
    }
  rch = pb->get_champ(un_nom_inco_etat);

  if(sub_type(Champ_Inc_base,rch.valeur()))
    l_inconnue_etat = ref_cast(Champ_Inc_base,rch.valeur()) ;
  else
    {
      Cerr<<pb->le_nom()<<" has no unknown field named "<<un_nom_inco_etat<<finl;
      exit();
    }
  Cerr <<"  l_inconnue_etat.que_suis_je() = "<<l_inconnue_etat->valeurs().que_suis_je()<<finl; //DoubleTab
  Cerr <<"  l_inconnue_etat nb comp = "<<l_inconnue_etat->nb_comp()<<finl;// nb_composants vitesse (2,3) ou 1 pour la pression
  Cerr <<"  l_inconnue_etat nom = "<<l_inconnue_etat->le_nom()<<finl;//vitesse ou pression
  Cerr<<" innconue valeurs  =  "<<l_inconnue_etat->valeurs()<<finl;// les valeurs
}

void Navier_Stokes_std_sensibility::mettre_a_jour(double temps)
{
  // Mise a jour de la classe mere (on tourne la roue).
  Navier_Stokes_std::mettre_a_jour(temps);
  associer_champ_evaluateur( nom_pb_etat, nom_inco_etat);
}


