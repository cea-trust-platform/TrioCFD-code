/****************************************************************************
* Copyright (c) 2020, CEA
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
// Directory : $Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Navier_Stokes_std_sensibility.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Param.h>
#include <Debog.h>
#include <LecFicDiffuse.h>
#include <communications.h>
#include <Interprete.h>

Implemente_instanciable_sans_constructeur_ni_destructeur( Navier_Stokes_std_sensibility, "Navier_Stokes_standard_sensibility", Navier_Stokes_std) ;
// XD Navier_Stokes_standard_sensibility navier_stokes_standard Navier_Stokes_standard_sensibility -1 Resolution of Navier-Stokes sensitivity problem
Navier_Stokes_std_sensibility::Navier_Stokes_std_sensibility() :  poly_chaos(0)
{

}
Navier_Stokes_std_sensibility::~Navier_Stokes_std_sensibility()
{

}
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
  param.ajouter_non_std("state",(this),Param::REQUIRED); // XD_ADD_P bloc_lecture Block to indicate the state problem. Between the braces, you must specify the key word 'pb_champ_evaluateur' then the name of the state problem and the velocity unknown  NL2 Example:  state { pb_champ_evaluateur pb_state  velocity }
  param.ajouter_non_std("uncertain_variable",(this),Param::REQUIRED); // XD_ADD_P bloc_lecture Block to indicate the name of the uncertain variable. Between the braces, you must specify the name of the unknown variable. Choice between velocity and mu.  NL2 Example: uncertain_variable { velocity }
  param.ajouter_non_std("polynomial_chaos",(this),Param::OPTIONAL); // XD_ADD_P floattant It is the method that we will use to study the sensitivity of the Navier Stokes equation: NL2 if poly_chaos=0, the sensitivity will be treated by the standard sentivity method. If different than 0, it will be treated by the polynomial chaos method

  /* if (schema_temps().diffusion_implicite())
     {
       Cerr<<"diffusion implicite forbidden within Navier_Stokes_std_sensibility  "<<finl;
       exit();
     }*/
  if( schema_temps().que_suis_je() != "Schema_euler_explicite" )
    {
      Cerr<<"Time  scheme: "<<schema_temps().que_suis_je() <<finl;
      Cerr<<"Only  Scheme_euler_explicit time scheme within Navier_Stokes_std_sensibility  "<<finl;
      exit();
    }

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
          is >> name_state_pb >>name_state_field;
          lu_info_evaluateur = 1;
          if( (name_state_pb  ==   accolade_fermee)  || (name_state_pb  ==  accolade_ouverte)   || (name_state_field  ==   accolade_fermee)
              || (name_state_field  ==  accolade_ouverte) )
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
      associate_evaluator_field( name_state_pb,name_state_field);
      return 1;
    }
  else  if (mot=="uncertain_variable")
    {

      Cerr << "Reading and typing of the  uncertain variable: " << finl;
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
      uncertain_var =  motlu;
      is>>motlu;
      if(motlu != accolade_fermee)
        {
          Cerr << "We expected a } while reading of " << que_suis_je() << finl;
          Cerr << "and not : " << motlu << finl;
          exit();
        }
      return 1;
    }
  else if (mot=="polynomial_chaos")
    {
      Cerr << "Reading and typing of the  option polynomial chaos: " << finl;
      double value;
      is >> value;
      poly_chaos = value;
      return 1;
    }
  else
    return Navier_Stokes_std::lire_motcle_non_standard(mot, is);
}

//Initialization of the reference to the evaluator field (state_field)
void Navier_Stokes_std_sensibility::associate_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field)
{
  Cerr <<"Navier_Stokes_std_sensibility: recoup state from "<<one_name_state_pb<< finl;
  if(probleme().le_nom() ==one_name_state_pb )
    {
      Cerr <<"Navier_Stokes_std_sensibility: we expect here the nom of the state problem and not the name of the sensibility problem"<< finl;
      exit();
    }
  Objet_U& ob= Interprete::objet(one_name_state_pb);
  OBS_PTR(Probleme_base) pb;
  OBS_PTR(Champ_base) rch;

  if(sub_type(Probleme_base,ob))
    {
      pb = ref_cast(Probleme_base,ob);
    }
  else
    {
      Cerr <<"No problem named "<< one_name_state_pb <<" has been found."<< finl;
      exit();
    }
  rch = pb->get_champ(one_name_state_field);

  if(sub_type(Champ_Inc_base,rch.valeur()))
    state_field = ref_cast(Champ_Inc_base,rch.valeur()) ;
  else
    {
      Cerr<<pb->le_nom()<<" has no unknown field named "<<one_name_state_field<<finl;
      exit();
    }
}
//Update the reference to the evaluator field (state_field)
void Navier_Stokes_std_sensibility::update_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field)
{
  Cerr <<"Navier_Stokes_std_sensibility: update  state from "<<one_name_state_pb<< finl;

  Objet_U& ob= Interprete::objet(one_name_state_pb);
  OBS_PTR(Probleme_base) pb;
  OBS_PTR(Champ_base) rch;

  pb = ref_cast(Probleme_base,ob);
  rch = pb->get_champ(one_name_state_field);
  state_field = ref_cast(Champ_Inc_base,rch.valeur()) ;

// Cerr <<"  state_field.que_suis_je() = "<<state_field->valeurs().que_suis_je()<<finl; //DoubleTab
  //Cerr <<"  state_field nb comp = "<<state_field->nb_comp()<<finl;// nb_composants vitesse (2,3) ou 1 pour la pression
  //Cerr <<"  state_field nom = "<<state_field->le_nom()<<finl;//vitesse ou pression
  //Cerr<<" innconue valeurs  =  "<<state_field->valeurs()<<finl;// les valeurs
}
void Navier_Stokes_std_sensibility::mettre_a_jour(double temps)
{
  // Mise a jour de la classe mere (on tourne la roue).
  Navier_Stokes_std::mettre_a_jour(temps);
  update_evaluator_field( name_state_pb,name_state_field);
}
const DoubleTab& Navier_Stokes_std_sensibility::get_state_field() const
{
  return state_field->valeurs();
}
const Champ_Inc_base& Navier_Stokes_std_sensibility::get_state() const
{
  return state_field;
}

const Motcle& Navier_Stokes_std_sensibility::get_uncertain_variable_name() const
{
  return uncertain_var;
}

const double& Navier_Stokes_std_sensibility::get_poly_chaos_value() const
{

  return poly_chaos;
}
