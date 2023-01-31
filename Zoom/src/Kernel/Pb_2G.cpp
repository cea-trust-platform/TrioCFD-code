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
// File:        Pb_2G.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#include <Pb_2G.h>
#include <Pb_MG.h>
#include <Champ_front_zoom.h>
#include <Schema_Temps_base.h>
#include <Domaine_VF.h>
#include <Equation_base.h>




Implemente_instanciable(Pb_2G,"Pb_2G",Objet_U);

//// printOn
//

Sortie& Pb_2G::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Pb_2G::readOn(Entree& s )
{
  return s ;
}


/*! @brief Calcul des connectivites entre faces fines et faces grossieres, et entre elements fins et elements grossiers
 *
 */
void Pb_2G::calculer_connectivites_2G()
{
  //Cerr<<"debut de Pb_2G::calculer_connectivites_2G"<<finl;

  int i, j;
  Probleme_base& pbF = pb_Fin();
  Probleme_base& pb_G = pbG();

  Domaine_dis& domF = pbF.domaine_dis();
  Domaine_dis_base& domaine_disF = domF;
  Domaine_VF& domaine_VFF = ref_cast(Domaine_VF, domaine_disF);

  Domaine& domaineG = pb_G.domaine();
  Domaine_dis& domG = pb_G.domaine_dis();
  Domaine_dis_base& domaine_disG = domG;
  Domaine_VF& domaine_VFG = ref_cast(Domaine_VF, domaine_disG);


  IntVect& connect_ff = connectivites().connectivites_faceF_faceG();
  IntVect& connect_ee = connectivites().connectivites_elemF_elemG();

  connectivites_ff_ee.calculer_connectivites(domaine_VFF, domaine_VFG, domaineG);
  if (nb_prolongement_>0)
    {
      /*if (sub_type(Prolongement_elem_face_FMG,mon_prolongement_.valeur()))
        for(i=0; i<nb_prolongement_; i++)
        mon_prolongement_(i).calculer(domaine_VFF, domaine_VFG, connect_ee);
        else */
      for(i=0; i<nb_prolongement_; i++)
        {
          mon_prolongement_(i).calculer(domaine_VFF, domaine_VFG, connect_ff);
        }
    }

  if (nb_restriction_>0)
    {

      for(j=0; j<nb_restriction_; j++)
        {
          Restriction_base& rest = ma_restriction_(j).valeur();
          /*if(sub_type(Restriction_elemF_faceG_1, rest))
            rest.calculer(domaine_VFG, domaine_VFF, connect_ff);
            else */
          rest.calculer(domaine_VFG, domaine_VFF, connect_ee);
        }
    }
  //Cerr<<"fin de Pb_2G::calculer_connectivites_2G"<<finl;
}


/*! @brief Prolongement de l'inconnue grossiere sur la frontiere fine pour inconnue TEMPERATURE ET PRESSION
 *
 *     en VDF pour T du domG, en VDF et VEF pour p du domG
 *     remarque : pas de couplage!!
 *     car valeur
 *
 */
void Pb_2G::prolonger_2G(IntVect& connect, DoubleTab& valG, int nb_compo,
                         const Frontiere& frontF,
                         DoubleTab& tab, int num_prolongement)
{
  //Cerr<<"debut de Pb_2G::prolonger_2G"<<finl;


  //AJOUT
  //ON INITIALISE TAB A ZERO
  tab = 0.;
  //////////////////////////
  //////////////////////////


  //Tableau de connectivites necessaire au prolongement
  //  Connectivites& connections = connectivites();
  //   IntVect& connect = connections.connectivites_faceF_faceG();
  //Tableau des distances


  //probleme fin
  Probleme_base& pbF = pb_Fin();
  Domaine_dis_base& domainef = pbF.domaine_dis();
  Domaine_VF& domaine_VFF = ref_cast(Domaine_VF, domainef);
  //probleme grossier
  Probleme_base& pb_G = pbG();
  Domaine_dis_base& domaineg = pb_G.domaine_dis();
  Domaine_VF& domaine_VFG = ref_cast(Domaine_VF, domaineg);

  //Cout << "Inco Grossiere = " << incoG.valeurs() << finl;

  mon_prolongement_(num_prolongement).prolonger(domaine_VFG, domaine_VFF, frontF, connect, valG, tab, nb_compo);

  //Cerr<<"fin de Pb_2G::prolonger_2G"<<finl;
}



/*! @brief Correction du probleme grossier!! restriction de l'inconnue fine sur le domaine grossier
 *
 *     pour inconnue TEMPERATURE ET PRESSION
 *     en VDF-VDF pour T du domG, en VDF-VDF, VDF-VEF, VEF-VEF
 *     et VEF-VDF pour p du domG
 *     remarque : pas de couplage!!
 *     car valeur
 *
 */
void Pb_2G::restreindre_2G(IntVect& connect, DoubleTab& val_incoG, const DoubleTab& val_incoF, int nb_comp_incoG, int num_rest)
{
  //probleme fin
  Probleme_base& pbF = pb_Fin();
  Domaine_dis_base& domainef = pbF.domaine_dis();
  Domaine_VF& domaine_VFF = ref_cast(Domaine_VF, domainef);
  //probleme grossier
  Probleme_base& pb_G = pbG();
  Domaine_dis_base& domaineg = pb_G.domaine_dis();
  Domaine_VF& domaine_VFG = ref_cast(Domaine_VF, domaineg);
  //Tableau de connectivites necessaire a la restriction
  //   Connectivites& connections = connectivites();
  //   IntVect& connect = connections.connectivites_elemF_elemG();
  //les inconnues  //A VOIR!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //   Champ_Inc_base& incoG = pbG.equation(0).inconnue();
  //   Champ_Inc_base& incoF = pbF.equation(0).inconnue();
  //   Cerr<<"nom de l'incoG : "<<incoG.le_nom()<<finl;
  //   Cerr<<"nom de l'incoF : "<<incoF.le_nom()<<finl;

  ma_restriction_(num_rest).restreindre(domaine_VFG, domaine_VFF,connect, val_incoG, val_incoF,nb_comp_incoG);
}





//POUR LA RESTRICTION_ELEMF_ELEMG_1
void Pb_2G::restreindre_2G(IntVect& connect, DoubleTab& val_incoG, const DoubleTab& val_incoF, int nb_comp_incoG, int num_prem_face_frontG, int num_rest)
{
  //probleme fin
  Probleme_base& pbF = pb_Fin();
  Domaine_dis_base& domainef = pbF.domaine_dis();
  Domaine_VF& domaine_VFF = ref_cast(Domaine_VF, domainef);
  //probleme grossier
  Probleme_base& pb_G = pbG();
  Domaine_dis_base& domaineg = pb_G.domaine_dis();
  Domaine_VF& domaine_VFG = ref_cast(Domaine_VF, domaineg);

  ma_restriction_(num_rest).restreindre(domaine_VFG, domaine_VFF,connect, val_incoG, val_incoF,nb_comp_incoG, num_prem_face_frontG);
}






//preparation des calculs pour son pb fin
void Pb_2G::preparer_calcul_pbFin()
{
  pb_Fin().preparer_calcul();
}


//postraitement du probleme fin
void Pb_2G::postraiter_pbFin()
{
  pb_Fin().postraiter();
}


//prolongement de la solution grossiere
//et resolution du probleme fin
void Pb_2G::prolonger_et_resoudre_pbFin()
{

  Probleme_base& pbF=pb_Fin();
  Schema_Temps_base& schF = pbF.schema_temps();
  schF.initTimeStep(schF.pas_de_temps());
  pbF.solveTimeStep();
  //mise a jour du probleme fin en temps
  double tempsF=schF.temps_courant()+schF.pas_de_temps();
  pbF.mettre_a_jour(tempsF);
  schF.mettre_a_jour();

  //impression des resultats
  schF.imprimer(Cout);
  pbF.imprimer(Cout);
  Cerr<<"end of Pb_fin::prolonger_et_resoudre_pbF"<<finl;

}

Probleme_base& Pb_2G::pbG()
{
  assert(pb_MG.non_nul());
  return pb_MG->pbG_MG();
}

void Pb_2G::associer_probleme_fin(Pb_MG& mg, int index)
{
  pb_MG=mg;
  index_pb_fin=index;
}

Probleme_base& Pb_2G::pb_Fin()
{
  assert(pb_MG.non_nul());
  return ref_cast(Probleme_base,pb_MG->probleme(index_pb_fin));
}

const Probleme_base& Pb_2G::pb_Fin() const
{
  assert(pb_MG.non_nul());
  return ref_cast(Probleme_base,pb_MG->probleme(index_pb_fin));
}


