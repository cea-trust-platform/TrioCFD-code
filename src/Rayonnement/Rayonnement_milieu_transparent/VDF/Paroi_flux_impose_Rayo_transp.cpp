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
// File:        Paroi_flux_impose_Rayo_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src/VDF
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_flux_impose_Rayo_transp.h>
#include <Domaine_VDF.h>

#include <Equation_base.h>
#include <Milieu_base.h>
#include <Champ_Uniforme.h>
#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
//
// printOn et readOn

Implemente_instanciable(Paroi_flux_impose_Rayo_transp,"Paroi_flux_impose_Rayo_transp",Paroi_Rayo_transp);
Sortie& Paroi_flux_impose_Rayo_transp::printOn(Sortie& s ) const
{
  return s;
}

Entree& Paroi_flux_impose_Rayo_transp::readOn(Entree& is )
{
  return Neumann_paroi::readOn(is);
}

void Paroi_flux_impose_Rayo_transp::mettre_a_jour(double temps)
{
  Neumann_paroi::mettre_a_jour(temps);
  calculer_Teta_i();
}
void Paroi_flux_impose_Rayo_transp::completer()
{
  preparer_surface(frontiere_dis(),domaine_Cl_dis());
  is_VDF_=0;

  if (sub_type(Domaine_VDF,domaine_Cl_dis().domaine_dis()))
    {
      is_VDF_=1;
      domaine_VDF = ref_cast(Domaine_VDF,domaine_Cl_dis().domaine_dis());

      // initialisation de Teta_i
      const Domaine_VDF& zvdf = domaine_VDF.valeur();

      const DoubleTab& T_f = mon_dom_cl_dis->equation().inconnue().valeurs();
      const Front_VF& la_frontiere_VF = ref_cast(Front_VF,frontiere_dis());
      int ndeb = la_frontiere_VF.num_premiere_face();
      int nb_faces_bord = la_frontiere_VF.nb_faces();

      const IntTab& face_voisins = zvdf.face_voisins();



      for (int numfa=0; numfa<nb_faces_bord; numfa++)
        {
          int elem = face_voisins(numfa+ndeb,0);
          if (elem<0)
            elem = face_voisins(numfa+ndeb,1);
          Teta_i(numfa) =T_f(elem);
        }
    }
  else
    {
      // initialisation de Teta_i en VEF
      const DoubleTab& T_p = mon_dom_cl_dis->equation().inconnue().valeurs();
      const Front_VF& la_frontiere_VF = ref_cast(Front_VF,frontiere_dis());
      int ndeb = la_frontiere_VF.num_premiere_face();
      int nb_faces_bord = la_frontiere_VF.nb_faces();

      for (int numfa=0; numfa<nb_faces_bord; numfa++)
        {
          Teta_i[numfa]=T_p(numfa+ndeb);
        }
    }
  //  Cerr<<"ici completer "<<Teta_i<<finl;

}

void Paroi_flux_impose_Rayo_transp::calculer_Teta_i()
{
  if (is_VDF_)
    calculer_Teta_i_VDF();
  else
    calculer_Teta_i_VEF();
}


void Paroi_flux_impose_Rayo_transp::calculer_Teta_i_VDF()
{
  const Domaine_VDF& le_dom_vdf = domaine_VDF.valeur();
  const Milieu_base& le_milieu = mon_dom_cl_dis->equation().milieu();
  ////const Champ_Uniforme& Lambda = ref_cast(Champ_Uniforme,le_milieu.conductivite().valeur());
  ////double d_Lambda = Lambda(0,0);
  const DoubleTab& T_f = mon_dom_cl_dis->equation().inconnue().valeurs();
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF,frontiere_dis());
  int ndeb = la_frontiere_VF.num_premiere_face();
  int nb_faces_bord = la_frontiere_VF.nb_faces();
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_Cl_dis().domaine_dis());
  const IntTab& face_voisins = zvdf.face_voisins();
  int is_rho_unif=0;
  int is_conduc_unif=0;
  int is_Cp_unif=0;

  double d_rho=0;
  double d_Lambda=0;
  double d_Cp =0;

  const DoubleTab& rho=le_milieu.masse_volumique()->valeurs();
  const DoubleTab& Lambda = le_milieu.conductivite()->valeurs();
  const DoubleTab& Cp = le_milieu.capacite_calorifique()->valeurs();

  if (sub_type(Champ_Uniforme,le_milieu.masse_volumique().valeur()))
    {
      is_rho_unif=1;

      d_rho= rho(0,0);
    }
  if (sub_type(Champ_Uniforme,le_milieu.conductivite().valeur()))
    {
      is_conduc_unif=1;
      d_Lambda= Lambda(0,0);
    }

  if (sub_type(Champ_Uniforme,le_milieu.capacite_calorifique().valeur()))
    {
      is_Cp_unif=1;
      d_Cp = Cp(0,0);
    }

  Schema_Temps_base& sch = mon_dom_cl_dis->equation().probleme().schema_temps();
  double dt= sch.pas_de_temps() ;

  int is_relax=1;
  if (le_modele_rayo->relaxation()==0)
    is_relax=0;

  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    {
      // QUI QU'A BU ?????
      // T_f (numfa)!!!! balaise Tf(elem) plus judicieux!!!
      int elem = face_voisins(numfa+ndeb,0);
      if (elem<0)
        elem = face_voisins(numfa+ndeb,1);

      if (!is_conduc_unif)
        d_Lambda = Lambda(elem);
      double omega;
      double e = le_dom_vdf.dist_norm_bord(numfa+ndeb);
      if (is_relax)
        {

          if (!is_rho_unif)
            d_rho=rho(elem);
          if (!is_Cp_unif)
            d_Cp = Cp(elem);

          omega=d_Lambda*dt/(d_Lambda*dt+e*e*d_rho*d_Cp);
        }
      else
        omega=1.;
      // BUG
      //elem=numfa;
      //attention omega pirate
      //omega=1;
      //      Cerr<<"omega = "<<omega<<finl;
      double flux_radia=le_modele_rayo->flux_radiatif(numfa+ndeb);

      if (le_champ_front->valeurs().size()==1)
        {
          Teta_i(numfa) = omega*((le_champ_front->valeurs()(0,0)-flux_radia)
                                 /(d_Lambda/e) + T_f(elem))+(1-omega)*Teta_i(numfa);
        }
      else if (le_champ_front->valeurs().dimension(1)==1)
        {
//      MR : T_f(elem) a la place de T_f(elem,0) !!!
          Teta_i(numfa) = omega*((le_champ_front->valeurs()(numfa,0)-flux_radia)
                                 /(d_Lambda/e) + T_f(elem))+(1-omega)*Teta_i(numfa);
        }
      else
        {
          Cerr << "Paroi_flux_impose_Rayo_transp::calculer_Teta_i() erreur" << finl;
          exit();
        }
    }
  // Impression:
  if (zvdf.domaine().bords_a_imprimer().contient(la_frontiere_VF.le_nom()) && sch.limpr())
    {
      Cout << "Impression des temperatures de paroi sur la frontiere " << la_frontiere_VF.le_nom() << " :" << finl;
      Cout << "---------------------------------------------------------------------" << finl;
      for (int numfa=0; numfa<nb_faces_bord; numfa++)
        Cout << "T("<<numfa<<") : "<< Teta_i(numfa) << " K." << finl;
    }
}

void Paroi_flux_impose_Rayo_transp::calculer_Teta_i_VEF()
{
  const DoubleTab& T_p = mon_dom_cl_dis->equation().inconnue().valeurs();
  double Temp;
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF,frontiere_dis());
  int ndeb = la_frontiere_VF.num_premiere_face();
  int nb_faces_bord = la_frontiere_VF.nb_faces();
  int is_relax=1;
  if (le_modele_rayo->relaxation()==0)
    is_relax=0;
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    {
      Temp=T_p(numfa+ndeb);
      //Cerr<< " Verif "<< Temp<<" face "<<numfa+ndeb<<finl;
      // omega devra etre saisi par l'utilisateur
      double omega = 0.8;
      if (is_relax==0)
        omega=1.;
      //  omega=0;
      //Teta_i[numfa] = omega*(Temp/nb_faces_bord) + (1-omega)*Teta_i[numfa];
      //Cerr<<" Teta_i "<<numfa+ndeb << " "<<Teta_i[numfa] << " "<<Temp<<finl;
      Teta_i[numfa] = omega*(Temp) + (1-omega)*Teta_i[numfa];

    }
}
