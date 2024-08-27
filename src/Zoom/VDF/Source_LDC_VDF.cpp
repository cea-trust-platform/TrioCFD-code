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
// File:        Source_LDC_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_LDC_VDF.h>
#include <Domaine_VDF.h>
#include <Operateur.h>
#include <Schema_Temps_base.h>

Implemente_instanciable(Source_LDC_VDF,"Source_LDC_VDF_P0_VDF",Source_Correction_Deficitaire);


//// printOn
//

Sortie& Source_LDC_VDF::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Source_LDC_VDF::readOn(Entree& s )
{
  return s ;
}


void Source_LDC_VDF::associer_domaines(const Domaine_dis& domaine_dis,
                                       const Domaine_Cl_dis& domaine_cl_dis)
{
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_dis.valeur());
  // const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF,domaine_cl_dis.valeur());

  // const Domaine& mon_dom = zvdf.domaine();
  Champ_Inc_base& inco = equation().inconnue();
  int nb_compo = inco.nb_comp();
  //Cerr<<" nb_compo dans Source_LDC_VDF::associer_domaines = "<<nb_compo<<finl;

  if(nb_compo == 1)
    {
      int nb_elem = zvdf.nb_elem();
      la_correction.resize(nb_elem);
      le_residu.resize(nb_elem);
    }
  else
    {
      int nb_face = zvdf.nb_faces();
      la_correction.resize(nb_face);
      le_residu.resize(nb_face);
    }
  // Cout << "Dans Source_LDC_VDF Init la_correction = " << nb_elem << " " << la_correction << finl;

}


DoubleTab& Source_LDC_VDF::calculer_residu(Connectivites_IndGros& connect, Restriction_base& restriction, Equation_base& eq_fine)
{
  Equation_base& eq = equation();
  DoubleTab& passe = eq.inconnue()->passe();
  double dtgros = eq.schema_temps().pas_de_temps();
  DoubleTab& present = eq.inconnue()->valeurs();
  DoubleTab& present_fin = eq_fine.inconnue()->valeurs();
  DoubleTab& passe_fin = eq_fine.inconnue()->passe();
  double dt = eq_fine.schema_temps().pas_de_temps();
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF, eq.domaine_dis().valeur());
  const Domaine_VDF& le_dom_fine = ref_cast(Domaine_VDF, eq_fine.domaine_dis().valeur());
  const DoubleVect& volumes = le_dom.volumes();
  const DoubleVect& volumes_fin = le_dom_fine.volumes();
  const IntVect& indice_gros = connect.indice_gros();
  int i, nb_compo;

  le_residu = 0;

  // On calcule la somme des operateurs fins puis leur restriction
  DoubleTab op_fin, op_fin_restreint;

  //TENONS COMPTE DU NOMBRE DE COMPOSANTES DE L'INCONNUE

  /* --> JR */
  nb_compo = eq_fine.inconnue()->nb_comp();
  /* <-- JR */

  if(nb_compo == 1)
    {
      op_fin.resize(le_dom_fine.nb_elem());
      op_fin_restreint.resize(le_dom.nb_elem());
    }
  else
    {
      op_fin.resize(le_dom_fine.nb_faces());
      op_fin_restreint.resize(le_dom.nb_faces());
    }
  //////////////////


  //Cerr<<"somme des operateurs fins"<<finl;
  for(i=0; i<eq_fine.nombre_d_operateurs(); i++)
    {
      eq_fine.operateur(i).ajouter(op_fin) ;
    }

  for (i=0; i<present_fin.size(); i++)
    {
      op_fin(i)-=(present_fin(i)-passe_fin(i))*volumes_fin(i)/dt;
    }


  //restriction
  restriction.restreindre(le_dom, le_dom_fine,connect.connectivites_elemF_elemG(),op_fin_restreint , op_fin,1);
  la_correction = present;
  restriction.restreindre(le_dom, le_dom_fine,connect.connectivites_elemF_elemG(),la_correction , present_fin,1);

  Cout << " ope fin restreint = " << op_fin_restreint << finl;


  //Calcul maintenant de la somme des operateurs de la restriction

  Cout << "Dans Source_LDC_VDF on ajoute une contrib dans le calcul du residu : la_correct = "  << la_correction << finl;

  for(i=0; i<eq.nombre_d_operateurs(); i++)
    {
      Cout << " ope = " << eq.operateur(i).ajouter(la_correction,le_residu) << finl;
    }




  //le_residu+=op_fin_restreint;

  //ajout de la derivee en temps
  for (i=0; i<present.size(); i++)
    {
      if (indice_gros(i)==-1) le_residu(i) = 0;
      else le_residu(i)+= - (la_correction(i)-passe(i))*volumes(i)/dtgros - op_fin_restreint(i);
      //else le_residu(i)-=op_fin_restreint(i);
    }



  le_residu *= -1.;

  Cout << " Dans Source_LDC_VDF le_residu = " << le_residu << finl;


  //le_residu = 0;



  return le_residu;
}

