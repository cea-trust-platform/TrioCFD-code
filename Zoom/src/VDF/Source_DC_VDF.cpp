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
// File:        Source_DC_VDF.cpp
// Directory:   $TRUST_ROOT/src/Zoom/VDF
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_DC_VDF.h>
#include <Zone_VDF.h>
#include <Operateur.h>
#include <Equation_base.h>
#include <Schema_Temps_base.h>

Implemente_instanciable(Source_DC_VDF,"Source_DC_VDF_P0_VDF",Source_Correction_Deficitaire);


//// printOn
//

Sortie& Source_DC_VDF::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Source_DC_VDF::readOn(Entree& s )
{
  return s ;
}


void Source_DC_VDF::associer_zones(const Zone_dis& zone_dis,
                                   const Zone_Cl_dis& zone_cl_dis)
{
  const Zone_VDF& zvdf = ref_cast(Zone_VDF,zone_dis.valeur());
  // const Zone_Cl_VDF& zclvdf = ref_cast(Zone_Cl_VDF,zone_cl_dis.valeur());

  // const Zone& ma_zone = zvdf.zone();


  int nb_elem = zvdf.nb_elem();
  la_correction.resize(nb_elem);
  le_residu.resize(nb_elem);
  Cerr << "Dans Source_DC_VDF Init la_correction = " << nb_elem << " " << la_correction << finl;

}


DoubleTab& Source_DC_VDF::calculer_residu(Connectivites_base& connect, LIST(Prolongement)& P, Equation_base& eqG, Nom& type_eq)
{
  Cerr<<"debut de Source_DC_VDF::calculer_residu"<<finl;
  le_residu = 0;
  int i;

  /* ---> JR  */
  Cerr << "Source_DC_VDF::calculer_residu --------------> JRRRRRRRRRRRR" << finl;
  Cerr << "type_eq = " << type_eq << finl;
  /* <--- JR */

  if(type_eq == "cond")
    {
      Prolongement& P1 = P(0);

      Equation_base& eqF = equation();
      Zone_VDF& la_zone = ref_cast(Zone_VDF, eqG.zone_dis().valeur());
      Zone_VDF& la_zone_fine = ref_cast(Zone_VDF, eqF.zone_dis().valeur());
      DoubleTab& present = eqG.inconnue().valeurs();
      DoubleTab& passe = eqG.inconnue().passe();
      //int tailleG = present.dimension(0);
      DoubleTab& presentF = eqF.inconnue().valeurs();
      //DoubleTab& passeF = eqF.inconnue().passe();
      int tailleF = presentF.dimension(0);
      const DoubleVect& volumes_fin = la_zone_fine.volumes();
      //const DoubleVect& volumes_gros = la_zone.volumes();
      double dtGros = eqG.schema_temps().pas_de_temps();
      //double dtGros = eqG.schema_temps().pas_de_temps();

      DoubleTab presentG_prol;
      presentG_prol.resize(tailleF);
      presentG_prol = 0.;

      DoubleTab passeG_prol;
      passeG_prol.resize(tailleF);
      passeG_prol = 0.;


      const Bord front_fictive;


      Cerr<<"tailleF = "<<tailleF<<finl;
      Cerr<<"presentG = "<<present<<finl;
      //Prolongement de l'inconnue grossiere present et passe
      P1.prolonger(la_zone, la_zone_fine,front_fictive, connect.connectivites_elemF_elemG(), present, presentG_prol, 1);
      Cerr<<"presentG_prol = "<<presentG_prol<<finl;



      Cerr<<"passeG = "<<passe<<finl;
      P1.prolonger(la_zone, la_zone_fine,front_fictive, connect.connectivites_elemF_elemG(), passe, passeG_prol, 1);
      Cerr<<"passeG_prol = "<<passeG_prol<<finl;


      //somme des operateurs fins appliques a l'inconnue grossiere prolongee
      //dans le residu
      for(i=0; i<eqF.nombre_d_operateurs(); i++)
        {
          //Cerr << "Dans Source_DC_VDF on ajoute une contrib dans le calcul du residu : incoG_prol = "  <<incoG_prol  << finl;

          Cout << " ope = " << eqF.operateur(i).ajouter(presentG_prol,le_residu) << finl;
        }


      Cerr<<"voici le residu apres somme des operateurs fins"<<le_residu<<finl;




      DoubleTab toto(le_residu);
      toto = 0.;
      //somme des operateurs fins appliques a l'inconnue grossiere prolongee
      //dans le residu
      for(i=0; i<eqG.nombre_d_operateurs(); i++)
        {
          //Cout << "Dans Source_DC_VDF on ajoute une contrib dans le calcul du residu : incoG_prol = "  <<incoG_prol  << finl;

          Cout << " ope = " << eqG.operateur(i).ajouter(present,toto) << finl;
        }


      //Cerr<<"Laplacien de T  G : "<<toto<<finl;





      //ajout de (-derivee du prolongement)
      for (i=0; i<tailleF; i++)
        {
          le_residu(i) -= (presentG_prol(i)-passeG_prol(i))*volumes_fin(i)/dtGros;
          Cerr<<"num_faceF = "<<i<<finl;
          Cerr<<"passeG_prol(i) = "<<passeG_prol(i)<<finl;
          Cerr<<"presentG_prol(i) = "<<presentG_prol(i)<<finl;
          Cerr<<"volumes_fin = "<<volumes_fin(i)<<finl;
          Cerr<<"dtGros = "<<dtGros<<finl;

          //else le_residu(i)-=op_fin_restreint(i);
        }

      Cerr<<"le_residu apres ajout de (-derivee du prolongement)"<<le_residu<<finl;
      //le_residu = - le_residu;

    }


  else if(type_eq == "conv_diff")
    {
      Cerr<<"Je suis bien un pb de conduction et donc je sais calculer le residu !!"<<finl;
      Prolongement& P1 = P(2);

      Equation_base& eqF = equation();
      Zone_VDF& la_zone = ref_cast(Zone_VDF, eqG.zone_dis().valeur());
      Zone_VDF& la_zone_fine = ref_cast(Zone_VDF, eqF.zone_dis().valeur());
      DoubleTab& present = eqG.inconnue().valeurs();
      DoubleTab& passe = eqG.inconnue().passe();
      //int tailleG = present.dimension(0);
      DoubleTab& presentF = eqF.inconnue().valeurs();
      //DoubleTab& passeF = eqF.inconnue().passe();
      int tailleF = presentF.dimension(0);
      const DoubleVect& volumes_fin = la_zone_fine.volumes();
      //const DoubleVect& volumes_gros = la_zone.volumes();
      double dtGros = eqG.schema_temps().pas_de_temps();
      //double dtGros = eqG.schema_temps().pas_de_temps();

      DoubleTab presentG_prol;
      presentG_prol.resize(tailleF);
      presentG_prol = 0.;

      DoubleTab passeG_prol;
      passeG_prol.resize(tailleF);
      passeG_prol = 0.;


      const Bord front_fictive;


      //   Cerr<<"tailleF = "<<tailleF<<finl;
      //       Cerr<<"presentG = "<<present<<finl;
      //Prolongement de l'inconnue grossiere present et passe
      P1.prolonger(la_zone, la_zone_fine,front_fictive, connect.connectivites_elemF_elemG(), present, presentG_prol, 1);
      //  Cerr<<"presentG_prol = "<<presentG_prol<<finl;



      //       Cerr<<"passeG = "<<passe<<finl;
      P1.prolonger(la_zone, la_zone_fine,front_fictive, connect.connectivites_elemF_elemG(), passe, passeG_prol, 1);
      // Cerr<<"passeG_prol = "<<passeG_prol<<finl;

      //somme des operateurs fins appliques a l'inconnue grossiere prolongee
      //dans le residu
      for(i=0; i<eqF.nombre_d_operateurs(); i++)
        {
          //Cout << "Dans Source_DC_VDF on ajoute une contrib dans le calcul du residu : incoG_prol = "  <<incoG_prol  << finl;

          Cout << " ope = " << eqF.operateur(i).ajouter(presentG_prol,le_residu) << finl;
        }


      //   Cerr<<"voici le residu avant ajout de (-derivee du prolongement)"<<le_residu<<finl;




      //ajout de (-derivee du prolongement)
      for (i=0; i<tailleF; i++)
        {
          le_residu(i) =  le_residu(i) - (presentG_prol(i)-passeG_prol(i))*volumes_fin(i)/dtGros;
          //else le_residu(i)-=op_fin_restreint(i);
        }

      //le_residu = - le_residu;
    }



  else
    {
      Cerr<<"Ce type de probleme n'a pas ete traite !!"<<finl;
      exit();
    }

  return le_residu;
}

