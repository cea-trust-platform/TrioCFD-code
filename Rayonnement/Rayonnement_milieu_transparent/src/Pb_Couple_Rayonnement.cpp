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
// File:        Pb_Couple_Rayonnement.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src
// Version:     /main/27
//
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////

// Description de Pb_Couple_Rayonnement:Classe heritant de Pb_Couple
// Precondition :
// Parametre :
//     Signification :
//     Valeurs par defaut :
//     Contraintes :
//     Entree :
//     Entree/Sortie :
//     Sortie :
// Retour :
//     Signification :
//     Contraintes :
// Exception :
// Effets de bord :
// Postcondition :
//

#include <Pb_Couple_Rayonnement.h>
#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Paroi_Rayo_transp.h>
#include <Frontiere_Ouverte_temperature_imposee_Rayo_transp.h>
#include <Probleme_base.h>
#include <Echange_contact_Rayo_transp_VDF.h>
#include <Echange_externe_impose_rayo_transp.h>
#include <Temperature_imposee_paroi_rayo_transp.h>
#include <Frontiere_Ouverte_Rayo_transp.h>
#include <verif_cast.h>
#include <Probleme_base.h>

Implemente_instanciable(Pb_Couple_Rayonnement,"Pb_Couple_Rayonnement",Probleme_Couple);

Entree& Pb_Couple_Rayonnement::readOn(Entree& is)
{
  return is;
}

Sortie& Pb_Couple_Rayonnement::printOn(Sortie& os) const
{
  return Probleme_Couple::printOn(os);
}

void Pb_Couple_Rayonnement::initialize( )
{
  completer();
  Probleme_Couple::initialize();
  le_modele_rayo().preparer_calcul();
}

int Pb_Couple_Rayonnement::associer_(Objet_U& ob)
{
  int set_type_rayo=0;
  if (Probleme_Couple::associer_(ob))
    {
      set_type_rayo=1;

    }
  else
    {
      if( sub_type(Modele_Rayonnement_base, ob))
        {
          set_type_rayo=1;
          Cerr << "association du modele au pbc" << finl;
          le_modele_rayo_associe(ref_cast(Modele_Rayonnement_base, ob));
        }
      else
        return 0;
    }
  if (set_type_rayo==1)
    {
      int nb_pb_fluide=0;
      for(int l=0; l< nb_problemes(); l++)
        {

          Probleme_base& pb = ref_cast(Probleme_base,probleme(l));

          if (sub_type(Fluide_base,pb.milieu()))
            nb_pb_fluide++;
        }

      for(int l=0; l< nb_problemes(); l++)
        {

          Probleme_base& pb = ref_cast(Probleme_base,probleme(l));

          if (sub_type(Fluide_base,pb.milieu()))
            {
              Fluide_base& fluide=ref_cast(Fluide_base,pb.milieu());
              if (nb_pb_fluide==1)
                fluide.fixer_type_rayo();
              else
                fluide.reset_type_rayo();
            }
        }

      return 1;
    }
  else
    return 0;
}

void Pb_Couple_Rayonnement::le_modele_rayo_associe(const Modele_Rayonnement_base& un_modele_de_rayonnement)
{
  le_modele_de_rayo = un_modele_de_rayonnement;
}

void Pb_Couple_Rayonnement::associer_cl_base(const Cond_lim_base& les_cl_)
{
  les_cl = les_cl_;
}

int Pb_Couple_Rayonnement::postraiter(int force)
{
  int ok=Probleme_Couple::postraiter(force);
  if (!ok)
    return 0;

  // Impression en plus du modele de rayonnement
  const Modele_Rayonnement_Milieu_Transparent& mod_rayo  =
    ref_cast(Modele_Rayonnement_Milieu_Transparent,le_modele_de_rayo.valeur());
  if (mod_rayo.processeur_rayonnant() != -1)
    if (schema_temps().limpr())
      {
        Cout << "Impression des flux radiatifs sur les bords de rayonnement" << finl;
        Cout << "----------------------------------------------------------------"<< finl;
        le_modele_de_rayo->imprimer_flux_radiatifs(Cout);
      }
  return 1;
}

void Pb_Couple_Rayonnement::validateTimeStep()
{
  Probleme_Couple::validateTimeStep();
  le_modele_rayo().mettre_a_jour(presentTime());
}


int is_la_cl_rayo (const Cond_lim_base& la_cl,Cond_Lim_Rayo*& la_cl_rayo)
{
  //if  (sub_type(Paroi_Rayo_transp,la_cl)|| sub_type(Frontiere_Ouverte_temperature_imposee_Rayo_transp,la_cl))

  if (sub_type(Paroi_Rayo_transp,la_cl))
    {
      la_cl_rayo =&( verif_cast(Cond_Lim_Rayo&, ref_cast(Paroi_Rayo_transp,la_cl)));
      return 1;
    }
  if sub_type(Frontiere_Ouverte_temperature_imposee_Rayo_transp,la_cl)
    {
      la_cl_rayo =&( verif_cast(Cond_Lim_Rayo&,  ref_cast(Frontiere_Ouverte_temperature_imposee_Rayo_transp,la_cl)));
      return 1;
    }
  if sub_type( Echange_contact_Rayo_transp_VDF,la_cl)
    {
      la_cl_rayo =&( verif_cast(Cond_Lim_Rayo&,  ref_cast( Echange_contact_Rayo_transp_VDF,la_cl)));
      return 1;
    }
  if sub_type( Echange_externe_impose_rayo_transp,la_cl)
    {
      la_cl_rayo =&( verif_cast(Cond_Lim_Rayo&,  ref_cast( Echange_externe_impose_rayo_transp,la_cl)));
      return 1;
    }
  if sub_type( Temperature_imposee_paroi_rayo_transp,la_cl)
    {
      la_cl_rayo =&( verif_cast(Cond_Lim_Rayo&,  ref_cast( Temperature_imposee_paroi_rayo_transp,la_cl)));
      return 1;
    }
  if sub_type( Frontiere_Ouverte_Rayo_transp,la_cl)
    {
      la_cl_rayo =&( verif_cast(Cond_Lim_Rayo&,  ref_cast(Frontiere_Ouverte_Rayo_transp ,la_cl)));
      return 1;
    }
  return 0;
}

void Pb_Couple_Rayonnement::completer()
{
  le_modele_de_rayo->discretiser(ref_cast(Probleme_base,probleme(0)).discretisation(), ref_cast(Probleme_base,probleme(0)).domaine());
  int nb_pb_fluide=0;
  int clrayo=0;
  Modele_Rayonnement_Milieu_Transparent& mod_rayo  =
    ref_cast(Modele_Rayonnement_Milieu_Transparent,le_modele_de_rayo.valeur());

  int compte_nb_bords_rayo =0;
  int l;
  int pb_fluide=-1;
  int is_pb_nom_existe=0;
  for(l=0; l< nb_problemes(); l++)
    {
      Probleme_base& le_pb = ref_cast(Probleme_base,probleme(l));
      if (sub_type(Fluide_base,le_pb.milieu()))
        {
          nb_pb_fluide++;
          pb_fluide=l;
          if (le_pb.le_nom()==mod_rayo.nom_pb_rayonnant())
            is_pb_nom_existe=1;
        }
    }
  if (nb_pb_fluide>1)
    {
      if (mod_rayo.nom_pb_rayonnant()=="non_donne")
        {
          Cerr<<"On ne sait traiter qu'un seul pb fluide"<<finl;
          Cerr<<" a moins d'indiquer le nom_pb_rayonnant au modele de rayonnement"<<finl;
          //  is_pb_fluide=0;
          exit();
        }
      else if (is_pb_nom_existe==0)
        {
          Cerr<<"On a au moins deux problemes fluides, le nom du pb rayonnant indique est "<<mod_rayo.nom_pb_rayonnant()<<" et ne corrrespond pas a un probleme fluide existant"<<finl;
          exit();
        }
    }
  else
    mod_rayo.nom_pb_rayonnant()= probleme(pb_fluide).le_nom();

  for(l=0; l< nb_problemes(); l++)
    {
      Probleme_base& le_pb = ref_cast(Probleme_base,probleme(l));
      //          Probleme_base& le_pb = probleme(l);
      int is_pb_fluide=( mod_rayo.nom_pb_rayonnant()== le_pb.le_nom());

      if (is_pb_fluide==1)
        {
          // normalement deja fait
          assert(le_pb.milieu().is_rayo_transp()==1);
          //          ref_cast(Fluide_base,le_pb.milieu()).fixer_type_rayo();
          Cerr << "Le probleme rayonnant trouve est : " << le_pb.le_nom() << finl;
        }
      for (int j=0; j<le_pb.nombre_d_equations(); j++)
        {
          Zone_Cl_dis& la_zcl = le_pb.equation(j).zone_Cl_dis();
          for (int num_cl=0; num_cl<la_zcl.nb_cond_lim(); num_cl++)
            {
              Cond_lim_base& la_cl = la_zcl.les_conditions_limites(num_cl).valeur();

              Cond_Lim_Rayo* la_cl_rayo;
              if (is_la_cl_rayo(la_cl,la_cl_rayo))
                {
                  ((*la_cl_rayo)).associer_modele_rayo(mod_rayo);

                  // on associe la cl liee au pb fluide
                  if (is_pb_fluide)
                    {
                      int ok=0;
                      for (int i =0; i <mod_rayo.nb_faces_totales(); i++)
                        {


                          if (mod_rayo.face_rayonnante(i).nom_bord_rayo()==la_zcl.les_conditions_limites(num_cl).frontiere_dis().le_nom())
                            //if (la_cl.frontiere_dis().frontiere().nb_faces()!=0)
                            {
                              if (mod_rayo.face_rayonnante(i).emissivite()!=-1)
                                ok=1;
                              //Cerr<< mod_rayo.face_rayonnante(i).nom_bord_rayo()<<" associe a "<<la_zcl.les_conditions_limites(num_cl).frontiere_dis().le_nom()<<finl;
                              clrayo++;
                              mod_rayo.face_rayonnante(i).ensembles_faces_bord(0).associer_les_cl(la_cl);
                              compte_nb_bords_rayo += 1;
                            }
                        }
                      if (ok==0)
                        {
                          Cerr<<"La condition limite de nom "<<la_zcl.les_conditions_limites(num_cl).frontiere_dis().le_nom()<<" est definie comme rayonnante, mais n'est pas dans la liste des faces rayonnantes ou son emissivite vaut -1"<<finl;
                          exit();
                        }
                    }
                }
            }
        }
      if (is_pb_fluide)
        {
          for (int i=0; i <mod_rayo.nb_faces_totales(); i++)
            {
              if (!mod_rayo.face_rayonnante(i).ensembles_faces_bord(0).is_ok()&& (mod_rayo.face_rayonnante(i).emissivite()!=-1) )
                {
                  Cerr<<"Le bord " << mod_rayo.face_rayonnante(i).nom_bord_rayo_lu()<<" n'a pas ete asssocie a une condition limite rayonnante."<<finl;
                  Cerr<<"Soit vous mettez une condition limite rayonnante pour "<<mod_rayo.face_rayonnante(i).nom_bord_rayo()<<finl;
                  Cerr<<"Soit vous affectez une emissivite de -1 a ce bord." << finl;
                  Cerr<<finl;
                }
            }
        }
    }
  //assert(nb_pb_fluide==1);

  if (compte_nb_bords_rayo!=mod_rayo.nb_faces_rayonnantes()) abort();
  if (nproc() == 1)
    {
      mod_rayo.associer_processeur_rayonnant(me());
    }
  else
    {
      //      if (compte_nb_bords_rayo == mod_rayo.nb_faces_rayonnantes())
      //         {
      //           mod_rayo.associer_processeur_rayonnant(me());
      //           Cerr<<"le processeur rayonnant est "<<me()<<finl;
      //         }
      //       else
      {
        if (compte_nb_bords_rayo!=0)
          {

            //Cerr<<me() <<nom_rayos<<finl;
            LIST(Nom) collectnoms;
            for (int i=0; i<mod_rayo.nb_faces_rayonnantes(); i++)
              {
                if (mod_rayo.face_rayonnante(i).ensembles_faces_bord(0).nb_faces_bord()!=0) collectnoms.add(mod_rayo.face_rayonnante(i).nom_bord_rayo_lu());
              }
            Cerr<<me() << collectnoms<<finl;
            // on verifie que l'on a bien tous les noms
            // pour verifier a la fin;
            /*
              for (int n=0;n<nb_proc()-1;n++)

              abort();
            */
            if (me()==0) mod_rayo.associer_processeur_rayonnant(me());
            else
              mod_rayo.associer_processeur_rayonnant(-1);




          }
        else
          {
            //tout est ok
            Cerr<<"On redimenssionne le tableau de faces de bord"<<finl;
            Cerr<<"compte_nb_bords_rayo = "<<compte_nb_bords_rayo<<finl;
            Cerr<<"mod_rayo.nb_faces_rayonnantes() = "<<mod_rayo.nb_faces_rayonnantes()<<finl;
            mod_rayo.associer_processeur_rayonnant(-1);
          }
      }
    }

  //  Cerr << "Pb_Couple_Rayonnement::completer() Fin" << finl;
}

void Pb_Couple_Rayonnement::discretiser(const Discretisation_base& dis)
{
  //   for(int i=0; i< nb_problemes(); i++)
  //     probleme(i).discretiser(dis);
  Probleme_Couple::discretiser(dis);


}
