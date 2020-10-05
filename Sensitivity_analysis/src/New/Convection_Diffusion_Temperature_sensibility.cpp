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
// File      : Convection_Diffusion_Temperature_sensibility.cpp
// Directory : $BALTIK_SENSIBILITY_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////

#include <Convection_Diffusion_Temperature_sensibility.h>
#include <Probleme_base.h>
#include <Schema_Temps.h>
#include <Param.h>
#include <Debog.h>
#include <LecFicDiffuse.h>
#include <communications.h>
#include <Interprete.h>

Implemente_instanciable( Convection_Diffusion_Temperature_sensibility, "Convection_Diffusion_Temperature_sensibility", Convection_Diffusion_Temperature ) ;
// XD Convection_Diffusion_Temperature_sensibility convection_diffusion_temperature  Convection_Diffusion_Temperature_sensibility -1 Energy sensitivity equation (temperature diffusion convection)
//XD attr velocity_state interprete velocity_state 0 Block to indicate the state problem. Between the braces, you must specify the key word 'pb_champ_evaluateur' then the name of the state problem and the velocity unknown  NL2 Example:  velocity_state { pb_champ_evaluateur pb_state  velocity }
//XD attr temperature_state interprete temperature_state 0 Block to indicate the state problem. Between the braces, you must specify the key word 'pb_champ_evaluateur' then the name of the state problem and the temperature unknown  NL2 Example:  velocity_state { pb_champ_evaluateur pb_state  temperature }
//XD attr uncertain_variable interprete uncertain_variable 0 Block to indicate the name of the uncertain variable. Between the braces, you must specify the name of the unknown variable (choice between: temperature, beta_th, boussinesq_temperature, Cp and lambda .  NL2 Example: uncertain_variable { temperature }
// XD  attr convection_sensibility convection_deriv sensibility 0 Choice between: amont and muscl   NL2 Example: convection {  Sensibility { amont } }


Sortie& Convection_Diffusion_Temperature_sensibility::printOn( Sortie& os ) const
{
  Convection_Diffusion_Temperature::printOn( os );
  return os;
}

Entree& Convection_Diffusion_Temperature_sensibility::readOn( Entree& is )
{
  Convection_Diffusion_Temperature::readOn( is );
  return is;
}
void Convection_Diffusion_Temperature_sensibility::set_param(Param& param)
{

  Convection_Diffusion_Temperature::set_param(param);
  param.ajouter_non_std("velocity_state",(this),Param::REQUIRED);
  param.ajouter_non_std("temperature_state",(this),Param::REQUIRED);
  param.ajouter_non_std("uncertain_variable",(this),Param::REQUIRED);

}


int Convection_Diffusion_Temperature_sensibility::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="velocity_state")
    {
      int lu_info_evaluateur = 0;
      Cerr << "Reading and typing of the velocity state : " << finl;
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
          is >> name_state_pb >>name_velocity_state_field;
          lu_info_evaluateur = 1;
          if( (name_state_pb  ==   accolade_fermee)  || (name_state_pb  ==  accolade_ouverte)   || (name_velocity_state_field  ==   accolade_fermee)
              || (name_velocity_state_field  ==  accolade_ouverte) )
            {
              Cerr<<"We expected the name of a problem and a fluid fild while reading of "<<motlu<< finl;
              Cerr << "and not : " <<accolade_ouverte << " or "<<  accolade_fermee << finl;
              exit();
            }
        }
      else
        {
          Cerr<<"Convection_Diffusion_Temperature_sensibility::lire_motcle_non_standard: keyword "<<motlu<<" is not recognized."<<finl;
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
          Cerr<<"when Convection_Diffusion_Temperature_sensibility is used."<<finl;
          exit();
        }
      const Motcle velocity;
      associate_evaluator_field( name_state_pb,name_velocity_state_field);
      return 1;
    }
  else if (mot=="temperature_state")
    {
      int lu_info_evaluateur = 0;
      Cerr << "Reading and typing of the temperature state : " << finl;
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
          is >> name_state_pb >>name_temperature_state_field;
          lu_info_evaluateur = 1;
          if( (name_state_pb  ==   accolade_fermee)  || (name_state_pb  ==  accolade_ouverte)   || (name_temperature_state_field  ==   accolade_fermee)
              || (name_temperature_state_field  ==  accolade_ouverte) )
            {
              Cerr<<"We expected the name of a problem and a fluid fild while reading of "<<motlu<< finl;
              Cerr << "and not : " <<accolade_ouverte << " or "<<  accolade_fermee << finl;
              exit();
            }
        }
      else
        {
          Cerr<<"Convection_Diffusion_Temperature_sensibility::lire_motcle_non_standard: keyword "<<motlu<<" is not recognized."<<finl;
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
          Cerr<<"when Convection_Diffusion_Temperature_sensibility is used."<<finl;
          exit();
        }
      associate_evaluator_field( name_state_pb,name_temperature_state_field);
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
  else
    return Convection_Diffusion_Temperature::lire_motcle_non_standard(mot, is);
}

//Initialization of the reference to the evaluator field (state_field)
void Convection_Diffusion_Temperature_sensibility::associate_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field)
{
  Cerr <<"Convection_Diffusion_Temperature_sensibility::get_state_field(): recoup state from "<<one_name_state_pb<< finl;
  if(probleme().le_nom() ==one_name_state_pb )
    {
      Cerr <<"Convection_Diffusion_Temperature_sensibility: we expect here the nom of the state problem and not the name of the sensibility problem"<< finl;
      exit();
    }
  Objet_U& ob= Interprete::objet(one_name_state_pb);
  REF(Probleme_base) pb;
  REF(Champ_base) rch;

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
    {
      if(one_name_state_field==name_velocity_state_field)
        velocity_state_field = ref_cast(Champ_Inc_base,rch.valeur()) ;
      else
        temperature_state_field = ref_cast(Champ_Inc_base,rch.valeur()) ;
    }
  else
    {
      Cerr<<pb->le_nom()<<" has no unknown field named "<<one_name_state_field<<finl;
      exit();
    }
}
//Update the reference to the evaluator field (state_field)
void Convection_Diffusion_Temperature_sensibility::update_evaluator_field(const Nom& one_name_state_pb,const Motcle& one_name_state_field)
{
  Cerr <<"Convection_Diffusion_Temperature_sensibility: update  state from "<<one_name_state_pb<< finl;

  Objet_U& ob= Interprete::objet(one_name_state_pb);
  REF(Probleme_base) pb;
  REF(Champ_base) rch;

  pb = ref_cast(Probleme_base,ob);
  rch = pb->get_champ(one_name_state_field);
  if(one_name_state_field==name_velocity_state_field)
    velocity_state_field = ref_cast(Champ_Inc_base,rch.valeur()) ;
  else
    temperature_state_field = ref_cast(Champ_Inc_base,rch.valeur()) ;
}
void Convection_Diffusion_Temperature_sensibility::mettre_a_jour(double temps)
{
  // Mise a jour de la classe mere (on tourne la roue).
  Convection_Diffusion_Temperature::mettre_a_jour(temps);
  update_evaluator_field( name_state_pb,name_velocity_state_field);
  update_evaluator_field( name_state_pb,name_temperature_state_field);

  /*  Cout <<"velocity_state_field.que_suis_je() = "<< velocity_state_field->valeurs().que_suis_je()<<finl; //DoubleTab
    Cout <<"velocity_state_field nb comp = "<< velocity_state_field->nb_comp()<<finl;// nb_composants vitesse (2,3) ou 1 pour la pression
    Cout <<"velocity_state_field nom = "<< velocity_state_field->le_nom()<<finl;//vitesse ou pression
    Cout<<"innconue valeurs  =  "<< velocity_state_field->valeurs()<<finl;// les valeurs
    getchar();

    Cout <<"temperature_state_field.que_suis_je() = "<< temperature_state_field->valeurs().que_suis_je()<<finl; //DoubleTab
    Cout <<"temperature_state_field nb comp = "<< temperature_state_field->nb_comp()<<finl;// nb_composants 1 pour la temperature
    Cout <<"temperature_state_field nom = "<< temperature_state_field->le_nom()<<finl;//temperature
    Cout<<"innconue valeurs  =  "<< temperature_state_field->valeurs()<<finl;// les valeurs
    getchar();*/

}
const DoubleTab& Convection_Diffusion_Temperature_sensibility::get_velocity_state_field() const
{
  return velocity_state_field->valeurs();
}
const DoubleTab& Convection_Diffusion_Temperature_sensibility::get_temperature_state_field() const
{
  return temperature_state_field->valeurs();
}

const Champ_Inc_base& Convection_Diffusion_Temperature_sensibility::get_velocity_state() const
{
  return velocity_state_field;
}
const Champ_Inc_base& Convection_Diffusion_Temperature_sensibility::get_temperature_state() const
{
  return temperature_state_field;
}

const Motcle& Convection_Diffusion_Temperature_sensibility::get_uncertain_variable_name() const
{
  return uncertain_var;
}



