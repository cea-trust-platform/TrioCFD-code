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
// File:        Source_LDC_VDF_NS.cpp
// Directory:   $TRUST_ROOT/src/Zoom/VDF
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_LDC_VDF_NS.h>
#include <Zone_VDF.h>
#include <Operateur.h>
#include <Equation_base.h>

Implemente_instanciable(Source_LDC_VDF_NS,"Source_LDC_VDF_Face",Source_Correction_Deficitaire);


//// printOn
//

Sortie& Source_LDC_VDF_NS::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Source_LDC_VDF_NS::readOn(Entree& s )
{
  return s ;
}


void Source_LDC_VDF_NS::associer_zones(const Zone_dis& zone_dis,
                                       const Zone_Cl_dis& zone_cl_dis)
{
  const Zone_VDF& zvdf = ref_cast(Zone_VDF,zone_dis.valeur());
  // const Zone_Cl_VDF& zclvdf = ref_cast(Zone_Cl_VDF,zone_cl_dis.valeur());

  // const Zone& ma_zone = zvdf.zone();
  Champ_Inc_base& inco = equation().inconnue();
  int nb_compo = inco.nb_comp();
  Cerr<<" nb_compo dans Source_LDC_VDF_NS::associer_zones = "<<nb_compo<<finl;

  int nb_face = zvdf.nb_faces();
  la_correction.resize(nb_face);
  le_residu.resize(nb_face);

  // Cout << "Dans Source_LDC_VDF_NS Init la_correction = " << nb_elem << " " << la_correction << finl;

}


DoubleTab& Source_LDC_VDF_NS::calculer_residu(Connectivites_IndGros& connect, Restriction_base& restriction, Equation_base& eq_fine)
{
  //Cerr<<"debut de Source_LDC_VDF_NS::calculer_residu"<<finl;
  Equation_base& eq = equation();
  // Equation_base& eq_fine = pb_fin.equation(0); // 1 seule equation pour l instant !!
  DoubleTab& present = eq.inconnue().valeurs();
  DoubleTab& present_fin = eq_fine.inconnue().valeurs();
  const Zone_VDF& la_zone = ref_cast(Zone_VDF, eq.zone_dis().valeur());
  const Zone_VDF& la_zone_fine = ref_cast(Zone_VDF, eq_fine.zone_dis().valeur());
  const IntVect& indice_gros = connect.indice_gros();
  int i;
  //AJOUT !
  Champ_Inc_base& incoG = eq.inconnue();
  int nb_compo = incoG.nb_comp();
  //int nbC;
  const IntTab& face_voisG = la_zone.face_voisins();
  int num_elemG;

  le_residu = 0;



  //On calcule le second membre grossier
  /* LIST_CURSEUR(Source) curseur(eq.sources());
     while(curseur){
     if (!sub_type(Source_Correction_Deficitaire,curseur.valeur().valeur()))
     {
     curseur->ajouter(le_residu);
     }
     ++curseur;
     }
  */

  //Cout << "source = " << le_residu << finl;


  // On calcule la somme des operateurs fins puis leur restriction
  DoubleTab op_fin, op_fin_restreint;


  op_fin.resize(la_zone_fine.nb_faces());
  op_fin_restreint.resize(la_zone.nb_faces());
  //////////////////


  //Cerr<<"somme des operateurs fins"<<finl;
  for(i=0; i<eq_fine.nombre_d_operateurs(); i++)
    {
      eq_fine.operateur(i).ajouter(op_fin) ;
    }
  /*
        for (i=0;i<present_fin.size();i++)
        {
        op_fin(i)-=(present_fin(i)-passe_fin(i))*volumes_fin(i)/dt;
        }

  */


  //restriction
  //Cerr<<"restriction de la somme des operateurs fins"<<finl;
  restriction.restreindre(la_zone, la_zone_fine,connect.connectivites_elemF_elemG(),op_fin_restreint , op_fin,nb_compo);
  la_correction = present;


  restriction.restreindre(la_zone, la_zone_fine,connect.connectivites_elemF_elemG(),la_correction , present_fin,nb_compo);

  //Cout << " ope fin restreint = " << op_fin_restreint << finl;


  //Calcul maintenant de la somme des operateurs de la restriction

  //Cerr<<"somme des operateurs de la restriction fine"<<finl;
  for(i=0; i<eq.nombre_d_operateurs(); i++)
    {
      Cout << "Dans Source_LDC_VDF_NS on ajoute une contrib dans le calcul du residu : la_correct = "  << la_correction << finl;

      Cout << " ope = " << eq.operateur(i).ajouter(la_correction,le_residu) << finl;
    }

  //AFFICHAGE !!
  //  Cerr<<"AFFICHAGE apres laplacien de la restriction"<<finl;
  //  for (i=0;i<present.size();i++)
  //    {
  //      Cerr<<"residu = "<<le_residu(i)<<finl;
  //    }


  //le_residu+=op_fin_restreint;

  //ajout de la derivee en temps
  //Cerr<<"ajout de la derivee en temps"<<finl;
  //ATTENTION : ON TRAITE PLUSIEURS FOIS LA MEME FACE
  //MAIS C'EST SANS CONSEQUENCE "MATHEMATIQUE ET PHYSIQUE"!!
  for (i=0; i<present.size(); i++)
    {
      num_elemG = face_voisG(i,0);
      if(num_elemG != -1)
        {
          if (indice_gros(num_elemG)==-1)
            le_residu(i) = 0;
          //        else le_residu(i)-=(la_correction(i)-present(i))*volumes(i)/dt;
          else
            le_residu(i)-=op_fin_restreint(i);
        }
      else
        {
          num_elemG = face_voisG(i,1);
          if(num_elemG != -1)
            {
              if (indice_gros(num_elemG)==-1)
                le_residu(i) = 0;
              //        else le_residu(i)-=(la_correction(i)-present(i))*volumes(i)/dt;
              else
                le_residu(i)-=op_fin_restreint(i);
            }
        }
    }


  //AFFICHAGE !!
  //  Cerr<<"AFFICHAGE apres ajout de la restriction du laplacien"<<finl;
  //  for (i=0;i<present.size();i++)
  //    {
  //      Cerr<<"residu = "<<le_residu(i)<<finl;
  //    }



  //Cerr<<"le_residu = - le_residu"<<finl;
  le_residu *= -1.;

  Cout << " Dans Source_LDC_VDF_NS le_residu = " << le_residu << finl;


  //  le_residu = 0;


  return le_residu;
}

