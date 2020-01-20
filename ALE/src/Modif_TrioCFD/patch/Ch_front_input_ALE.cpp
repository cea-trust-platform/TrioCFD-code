/****************************************************************************
 * Copyright (c) 2015 - 2018, CEA
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
// File:        Ch_front_input_ALE.cpp
// Directory:   New class
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Ch_front_input_ALE.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Frontiere_dis_base.h>
#include <communications.h>
#include <Convert_ICoCoTrioField.h>
#include <MD_Vector_tools.h>
#include <Roue.h>

Implemente_instanciable_sans_constructeur(Ch_front_input_ALE,"Ch_front_input_ALE",Ch_front_input);
// XD Ch_front_input_ALE  front_field_base  Ch_front_input_ALE  -1 Class to define a boundary condition on a moving boundary of a mesh (only for the Arbitrary Lagrangian-Eulerian  framework ) .  NL2  Example:  Ch_front_input_ALE { nb_comp 3 nom VITESSE_IN_ALE probleme pb initial_value 3 1. 0. 0. }

Ch_front_input_ALE::Ch_front_input_ALE():
        alreadyInit_(false)
{
}

Entree& Ch_front_input_ALE::readOn(Entree& is)
{
  return Ch_front_input::readOn(is);
}

Sortie& Ch_front_input_ALE::printOn(Sortie& os) const
{
  return Ch_front_input::printOn(os);
}

void Ch_front_input_ALE::getTemplate(TrioField& afield) const
{
  buildSommetsFaces(); //From Ch_front_input in order to update mesh
  Ch_front_input::getTemplate(afield);
}

void Ch_front_input_ALE::setValue(const TrioField& afield)
{

  int NrOfFacesIn_les_valeurs=les_valeurs[1].valeurs().size()/nb_comp();

  for(int i=0; i<NrOfFacesIn_les_valeurs; i++)
    {
      for(int j=0; j<afield._nb_field_components; j++)
	{
	  les_valeurs[1].valeurs()(i,j)=afield._field[afield._nb_field_components*i+j]; //Values to faces.
	}
    }
}

int Ch_front_input_ALE::initialiser(double temps, const Champ_Inc_base& inco)
{
  if (alreadyInit_)
    return 1;
  alreadyInit_ = true;
  if (!Ch_front_input::initialiser(temps,inco))
    return 0;

  //ALE part:
  mettre_a_jour(temps);

  return 1;
}

void Ch_front_input_ALE::mettre_a_jour(double temps)
{
  Cerr << "Champ_front_input_ALE ::mettre_a_jour" << finl;

  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Faces& faces=front.faces();
  int nbsf=faces.nb_som_faces();
  int i,j,k;

  DoubleTab& tab=valeurs_au_temps(temps);
  remplir_vit_som_bord_ALE(temps);
  tab=0.;

  for( i=0; i<nb_faces; i++)
    {
      for( j=0; j<nb_comp(); j++)
        {
          for( k=0; k<nbsf; k++)
            {
              tab(i,j)+=vit_som_bord_ALE(faces.sommet(i,k),j);
            }
        }
    }
  tab/=nbsf;
  tab.echange_espace_virtuel();
}

void Ch_front_input_ALE::remplir_vit_som_bord_ALE(double tps)
{
  Cerr<<"Champ_front_input_ALE::remplir_vit_som_bord_ALE"<<finl;

  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Zone& zone=front.zone();
  const Faces& faces=front.faces();
  const Domaine& domaine=zone.domaine();
  int nbsf=faces.nb_som_faces();
  int i,j,k;
  int nb_som_tot=domaine.nb_som_tot();
  vit_som_bord_ALE.resize(nb_som_tot,nb_comp());
  /*if (vit_som_bord_ALE.dimension(0) != domaine.nb_som())
    {
      vit_som_bord_ALE.resize(domaine.nb_som(),nb_comp());
      const MD_Vector& md = zone.domaine().md_vector_sommets();
      MD_Vector_tools::creer_tableau_distribue(md, vit_som_bord_ALE);
    }*/


  DoubleTab DenominatorValues(nb_som_tot); //To hold number of faces connected to each node.
  vit_som_bord_ALE=0.;
  double InitialTimeValue=mon_pb->schema_temps().temps_init();

  if(tps==InitialTimeValue){ //Initial values.
      for( i=0; i<nb_faces; i++)
        {
          for( k=0; k<nbsf; k++)
            {
              for( j=0; j<nb_comp(); j++)
                {
                  if (initial_value_.size_array()==0)
                    vit_som_bord_ALE(faces.sommet(i,k),j)=0; // By default initializes to zero.
                  else
                    vit_som_bord_ALE(faces.sommet(i,k),j)=initial_value_(j);
                }
            }
        }
  }
  else{ //Nodes will have the average value from the faces connected to it.
      for( i=0; i<nb_faces; i++)
        {
          for( k=0; k<nbsf; k++)
            {
        	  DenominatorValues(faces.sommet(i,k))=DenominatorValues(faces.sommet(i,k))+1.0;
              for( j=0; j<nb_comp(); j++)
                {
                  double val =les_valeurs[1].valeurs()(i,j);
                  vit_som_bord_ALE(faces.sommet(i,k),j)+=val; //Sum up velocity values from the connected faces.
                }
            }
        }

      for( i=0; i<nb_som_tot; i++)
        {
    	  if(DenominatorValues(i)>=1.e-12){
    		  for( j=0; j<nb_comp(); j++)
    		                  {
    			  	  			vit_som_bord_ALE(i,j)=vit_som_bord_ALE(i,j)/DenominatorValues(i); //Divide by the number of connected faces.
    		                  }
    	  }
        }
  }
  //vit_som_bord_ALE.echange_espace_virtuel();
}

