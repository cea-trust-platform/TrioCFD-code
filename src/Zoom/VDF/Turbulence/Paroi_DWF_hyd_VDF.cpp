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
// File:        Paroi_DWF_hyd_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Turbulence
// Version:     /main/23
//
//////////////////////////////////////////////////////////////////////////////


#include <Paroi_DWF_hyd_VDF.h>
#include <Navier_Stokes_std.h>
#include <Modele_turbulence_hyd_base.h>
#include <Fluide_Incompressible.h>
#include <Champ_Uniforme.h>
#include <Domaine_Cl_VDF.h>
#include <Source_DWF_VDF.h>
#include <EChaine.h>
#include <Connectivites_DWF.h>
#include <Schema_Implicite_base.h>
#include <Champ_front_zoom.h>
#include <Interprete.h>
#include <Convection_Diffusion_Temperature.h>

Implemente_instanciable(Paroi_DWF_hyd_VDF,"Dynamic_WF_VDF",Paroi_std_hyd_VDF);

///////////////////////////////////////////////////////////////////////////
//On conserve la classe car possible developpement dans le futur.
//Actuellement pas utilisee et ajout de commentaires pour
//compilation suite a integration de developpements Couplage_Multi_Echelles
///////////////////////////////////////////////////////////////////////////

//#################################################################
Sortie& Paroi_DWF_hyd_VDF::printOn(Sortie& s ) const
//#################################################################
{
  return s << que_suis_je() << " " << le_nom();
}

//#################################################################
Entree& Paroi_DWF_hyd_VDF::readOn(Entree& is )
//#################################################################
// WEC : could be entirely rewritten...
{
  Nom fine_grid;

  CHT=0;

  is >> fine_grid;
  const Objet_U& ob_pbfin=Interprete::objet(fine_grid); // le probleme fin thermohydraulique (PAS le PB COUPLE !)
  pb_fin=ref_cast(Probleme_base, ob_pbfin);

  Motcle mot_lu;
  is >> mot_lu;
  if(mot_lu=="CHT")
    {
      is >> CHT;
      if((CHT!=1)&&(CHT!=0))
        {
          Cerr << "Check the value (0 or 1) behind key-word CHT..." << finl;
          Cerr << "Trio will now stop." << finl;
          exit();
        }
    }

  // WEC : search for closing paren
  while (mot_lu!="}")
    is >> mot_lu;

  return is ;
}


//#################################################################
int Paroi_DWF_hyd_VDF::init_lois_paroi()
//#################################################################
{

  /**
   * Definition des liens entre probleme a resoudre et probleme fin
   * (servant a la loi de paroi dynamique)
   */

  thermo=0;

  Probleme_base& pb_thhyd = mon_modele_turb_hyd->equation().probleme();

  //Nom nom_pbmg("export Pb_MG pbMG");
  //EChaine ch_pbmg(nom_pbmg);
  //interpreter(ch_pbmg);

  Nom pbMGName("pbMG");
  const Objet_U& ob_pbmg=Interprete::objet(pbMGName);
  pbMG=ref_cast(Pb_MG, ob_pbmg);

  pbMG.associer_pbMG_pbGglobal_(pb_thhyd); // Probleme grossier  = le probleme a resoudre
  pbMG.associer_pbMG_pbFin_(pb_fin.valeur());

  // On n'utilise pas d'algorithme MG donc il faut calculer les connectivites "a la main"
  // sans passer par  pbMG.calculer_connectivites_MG()
  Pb_2G& pb2G = pbMG.pb_2G(0); // On a un seul probleme "2 grilles"
  pb2G.typer_Connectivites("Connectivites_DWF");
  pb2G.set_nb_prol(3);
  pb2G.set_nb_rest(0);
  OWN_PTR(Prolongement_base) P1;
  P1.typer("Prolongement_elem_elem_DWF"); // Pour la pression ou la temperature
  pb2G.mon_prolongement().add(P1);

  OWN_PTR(Prolongement_base) P2;
  P2.typer("Prolongement_face_face_DWF"); // Pour la vitesse
  pb2G.mon_prolongement().add(P2);

  OWN_PTR(Prolongement_base) P3;
  P3.typer("Prolongement_face_face_FMG"); // Pour les termes sources : prolongement partout
  pb2G.mon_prolongement().add(P3);

  pb2G.calculer_connectivites_2G();


  // ALEX : 10 juin 2003
  // On regarde pour chaque equation (QDM et Temperature) du probleme
  // s'il y a des CL de nom "Interface" et si c'est le cas,
  // on doit retyper maintenant le champ de la CL "Interface" en Champ_front_zoom
  for(int i=0; i<pb_fin->nombre_d_equations(); i++)
    {
      Conds_lim& les_cl = pb_fin->equation(i).domaine_Cl_dis().les_conditions_limites();
      for (int icl = 0; icl<les_cl.size(); icl++)
        {
          const Frontiere_dis_base& cl = les_cl[icl]->frontiere_dis();

          if (sub_type(Navier_Stokes_std,pb_fin->equation(i)))
            {
              const RefObjU& modele_turbulence = pb_fin->equation(i).get_modele(TURBULENCE);
              if (modele_turbulence.non_nul() && sub_type(Modele_turbulence_hyd_base,modele_turbulence.valeur()) && (cl.le_nom() == "Interface"))
                {
                  Nom champ("Champ_front_zoom pbMG ");
                  champ+=pb_fin->le_nom();
                  champ+=Nom(" ");
                  champ+=pb_thhyd.le_nom();
                  champ+=Nom(" ");
                  champ+=Nom(mon_modele_turb_hyd->equation().domaine_Cl_dis().les_conditions_limites(0)->frontiere_dis().le_nom());
                  champ+=Nom(" vitesse ");
                  EChaine chp(champ);
                  chp >> les_cl[icl]->champ_front();
                  les_cl[icl]->associer_fr_dis_base(cl);
                }
            }

          // ######################
          // # THERMIQUE S'IL Y A #
          // ######################
          else if(sub_type(Convection_Diffusion_Temperature,pb_fin->equation(i)) && (cl.le_nom() == "Interface"))
            {
              thermo=1;
              Nom champ("Champ_front_zoom pbMG ");
              champ+=pb_fin->le_nom();
              champ+=Nom(" ");
              champ+=pb_thhyd.le_nom();
              champ+=Nom(" ");
              // marchait avant
              //          champ+=Nom(modele_turbulence().equation().domaine_Cl_dis().les_conditions_limites(0)->frontiere_dis().le_nom());
              champ+=Nom(mon_modele_turb_hyd->equation().domaine_Cl_dis().les_conditions_limites(0)->frontiere_dis().le_nom());
              champ+=Nom(" temperature ");
              EChaine chp(champ);
              chp >> les_cl[icl]->champ_front();
              les_cl[icl]->associer_fr_dis_base(cl);
            }
          // #################
          // # FIN THERMIQUE #
          // #################

        }
    }



  pbMG.initialiser_champ_front_zoom(); // dimensionnement et mise a zero des champs frontieres

  //////////////////////////////////////////////////////////////////////////
  //Commente suite a integration des developpements Couplage_Multi_Echelles
  //////////////////////////////////////////////////////////////////////////
  ///if(CHT)
  ////{
  ////Probleme_Couple& pb_couple = pb_fin->probleme_couple();
  ////pb_couple.preparer_calcul(); // on prepare le probleme couple
  ////pb_couple.postraiter();
  ////}
  ////else pb_fin->preparer_calcul(); // on prepare le probleme fin

  Paroi_std_hyd_VDF::init_lois_paroi();

  return 1;
}


// calculer_hyd pour le k-epsilon
//#################################################################
int Paroi_DWF_hyd_VDF::calculer_hyd(DoubleTab& tab_k_eps)
//#################################################################
{
  Cerr << " Paroi_DWF_hyd_VDF::calculer_hyd(DoubleTab& tab_k_eps) " << finl;
  Cerr <<  "on ne doit pas entrer dans cette methode" << finl;
  Cerr << " car elle est definie uniquement pour la LES " << finl ;
  return 1 ;
}




//#################################################################
int Paroi_DWF_hyd_VDF::calculer_hyd(DoubleTab& tab_nu_t,DoubleTab& tab_k)
//#################################################################
{

  Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const int nb_face = domaine_VDF.nb_faces();

  const double temps = eqn_hydr.schema_temps().temps_courant();

  // la vitesse resolue par l'equation hydr
  const DoubleTab& vit = eqn_hydr.inconnue().valeurs();
  const int nb_compo = eqn_hydr.inconnue().nb_comp();

  // Variables sur le probleme fin
  Domaine_VDF& domaine_fine = ref_cast(Domaine_VDF,pb_fin->domaine_dis());


  // Calcul des termes sources grossiers
  DoubleTab sources_grossiers(nb_face);
  eqn_hydr.sources().calculer(sources_grossiers);


  // We must divide the source terms by the coarse embedded volumes.
  // In a few lines, we will multiply them by the embedded volumes
  // of the fine mesh.
  const DoubleVect& volentre_G = domaine_VDF.volumes_entrelaces();
  for(int i=0; i<nb_face; i++) sources_grossiers(i)/=volentre_G(i);

  // Remplissage des termes sources fins
  Equation_base& eqF = pb_fin->equation(0);
  Sources& les_sources =  eqF.sources();
  Connectivites_DWF& connect = ref_cast(Connectivites_DWF, pbMG.pb_2G(0).connectivites().valeur());
  for (int nb_source = 0; nb_source < les_sources.size(); nb_source++)
    {
      if (sub_type(Source_DWF_VDF,les_sources(nb_source).valeur()))
        {

          Source_DWF_VDF& mon_source = ref_cast(Source_DWF_VDF,les_sources(nb_source).valeur());
          DoubleTab& val  = mon_source.getValeurs();
          Prolongement_base& P = pbMG.pb_2G(0).mon_prolongement(2);
          const Bord front_fictive;
          P.prolonger(domaine_VDF, domaine_fine,front_fictive,connect.connectivites_elemF_elemG(),sources_grossiers , val ,1);

          // We now multiply the source terms by the embedded volumes
          // of the fine mesh.
          const int nb_face_F = domaine_fine.nb_faces();
          const DoubleVect& volentre_F = domaine_fine.volumes_entrelaces();
          for(int i=0; i<nb_face_F; i++) val(i)*=volentre_F(i);
          break;
        }
    }

  // Pas de temps grossier :
  double dt_gros = eqn_hydr.schema_temps().pas_de_temps();


  // Remplissage de la CL sur le bord de nom "Interface" du probleme fin
  Conds_lim& les_cl = pb_fin->equation(0).domaine_Cl_dis().les_conditions_limites();
  Prolongement_base& P_cl = pbMG.pb_2G(0).mon_prolongement(1);

  int icl,i_de_la_cl=0;

  for (icl = 0; icl<les_cl.size(); icl++)
    {
      const Frontiere_dis_base& cl = les_cl[icl]->frontiere_dis();
      if (cl.le_nom() == "Interface")
        {
          i_de_la_cl = icl;
        }
    }


  // On prolonge la vitesse sur ce bord
  const Frontiere_dis_base& cl = les_cl[i_de_la_cl]->frontiere_dis();
  DoubleTab& val  = les_cl[i_de_la_cl]->champ_front().valeurs();
  DoubleTab pente(val);

  P_cl.prolonger(domaine_VDF,domaine_fine,cl.frontiere(),connect.connectivites_elemF_elemG(),vit, pente,nb_compo);

  if(dt_gros>-1e-6)
    {
      pente -=val;
      pente/=dt_gros;
    }
  else
    {
      // Le pas de temps est tres petit : ca veut peut-etre dire
      // que c'est la premiere iteration et que val n'est pas encore initialise.
      // ... Donc on le fait !
      P_cl.prolonger(domaine_VDF,domaine_fine,cl.frontiere(),connect.connectivites_elemF_elemG(),vit, val,nb_compo);
      // Et on met pente a 0.
      pente = 0;
    }
  ref_cast(Champ_front_var_instationnaire,les_cl[i_de_la_cl]->champ_front()).set_derivee_en_temps(pente);

  // #############
  // # THERMIQUE #
  // #############
  if(!thermo)
    {
      Cerr << "DWF only works for THERMO-hydraulic problems" << finl;
      Cerr << "Trio will now stop." << finl;
      exit();
    }

  Equation_base& eqn_NRJ_F = pb_fin->equation(1);
  const Probleme_base& pb_gros = mon_modele_turb_hyd->equation().probleme();
  const DoubleTab& tempiotte_grosse = pb_gros.equation(1).inconnue().valeurs();

  Conds_lim& les_cl_th = eqn_NRJ_F.domaine_Cl_dis().les_conditions_limites();
  Prolongement_base& P_cl_th = pbMG.pb_2G(0).mon_prolongement(0);

  i_de_la_cl=0;

  for (icl = 0; icl<les_cl_th.size(); icl++)
    {
      const Frontiere_dis_base& cl_th = les_cl_th[icl]->frontiere_dis();
      if (cl_th.le_nom() == "Interface")
        {
          i_de_la_cl = icl;
        }
    }
  const Frontiere_dis_base& cl_th = les_cl_th[i_de_la_cl]->frontiere_dis();
  DoubleTab& val_th  = les_cl_th[i_de_la_cl]->champ_front().valeurs();

  DoubleTab pente_th(val_th);

  P_cl_th.prolonger(domaine_VDF,domaine_fine,cl_th.frontiere(),connect.connectivites_elemF_elemG(),tempiotte_grosse, pente_th,1);

  if(dt_gros>1e-6)
    {
      pente_th -=val_th;
      pente_th/=dt_gros;
    }
  else
    {
      // Le pas de temps est tres petit : ca veut peut-etre dire
      // que c'est la premiere iteration et que val_th n'est pas encore initialise.
      // ... Donc on le fait !
      P_cl_th.prolonger(domaine_VDF,domaine_fine,cl_th.frontiere(),connect.connectivites_elemF_elemG(),tempiotte_grosse, val_th,1);
      // Et on met pente a 0.
      pente_th = 0;
    }
  // #################
  // # FIN THERMIQUE #
  // #################

  Schema_Temps_base& sch = eqF.schema_temps();

  double temps_courant = sch.temps_courant();
  Probleme_base& pb_base = pbMG.pb_2G(0).pb_Fin();

  //******************************************************
  // On monitore le maximum de divergence
  //  Navier_Stokes_std& eqF_typee = ref_cast(Navier_Stokes_std, eqF);
  //  SFichier fic("Max_divU.dat",ios::app);
  //  Solveur_Null& solv_press_type = ref_cast(Solveur_Null,eqF_typee.solveur_pression().valeur());
  //  const double le_seuil_fixe = solv_press_type.crit_seuil;
  //******************************************************


  // On fait une sauvegarde s'il y a besoin
  // et de tous les problemes en meme temps pour etre sur de les avoir dans le



  ////////////////////////////////////////////////////////////////////////
  //Commente suite a integration des developpements Couplage_Multi_Echelles
  ////////////////////////////////////////////////////////////////////////
  if(sch.lsauv())
    {
      if(CHT)
        {
          // Sauvegarde probleme Fin
          ////Probleme_Couple& pb_coupleF = pb_fin->probleme_couple();
          ////pb_coupleF.sauver();
          // Sauvegarde probleme Grossier
          ////Probleme_base& pb_grosse = ref_cast(Probleme_base,mon_modele_turb_hyd->equation().probleme());
          ////Probleme_Couple& pb_coupleG = pb_grosse.probleme_couple();
          ////pb_coupleG.sauver();
        }
      else
        {
          // Sauvegarde probleme Fin
          pb_base.sauver();
          // Sauvegarde probleme Grossier
          Probleme_base& pb_grosse = ref_cast_non_const(Probleme_base,mon_modele_turb_hyd->equation().probleme());
          pb_grosse.sauver();
        }
    }





  //************************************************
  // On commence l'avancee en temps du probleme fin.
  //************************************************

  while (temps_courant < temps)
    {

      /******************************************************/
      // On monitore le maximum de divergence
      // const DoubleTab& div =eqF_typee.div().valeurs();
      //       for (int ii=0;ii<div.size_array();ii++)
      //         if (std::fabs(div(ii))>1e-3) Cerr<<" Pb en "<<ii<<" : "<<div(ii)<<finl;
      //      double divU = max_abs(eqF_typee.div().valeurs());
      //      fic << temps_courant << "     " << divU << finl;
      /******************************************************/
      //      if(divU>le_seuil_fixe) eqF_typee.projeter();
      double dt;
      if(CHT)
        {
          ////Probleme_Couple& pb_couple = pb_fin->probleme_couple();
          ////dt=pb_couple.calculer_pas_de_temps();
        }
      else dt=pb_base.calculer_pas_de_temps();


      // On met a jour le pas de temps du schema en temps du probleme fin
      // corrige pour le cas de l'implicite pour que le dt soit multiplie par le facsec.
      sch.corriger_dt_calcule(dt);

      // Cependant, si le pas de temps depasse le temps du probleme grossier,
      // on limite la valeur du dt fin pour arriver au temps grossier.
      ////if(temps_courant+dt> temps) sch.changer_pas_de_temps(dt=temps-temps_courant);


      // On construit la pente pour une condition aux limites qui varie lineairement.
      // La vitesse :
      val.ajoute_sans_ech_esp_virt(dt, pente, VECT_ALL_ITEMS);

      // ######################
      // # THERMIQUE S'IL Y A #
      // ######################
      // La temperature :
      val_th.ajoute_sans_ech_esp_virt(dt, pente_th, VECT_ALL_ITEMS);
      // #################
      // # FIN THERMIQUE #
      // #################

      //////////////////////////////////////////////////////////////////////////
      //Commente suite a integration des developpements Couplage_Multi_Echelles
      //////////////////////////////////////////////////////////////////////////

      ////sch.preparer_pas_temps();


      ////if(CHT)
      ////{
      ////Probleme_Couple& pb_couple = pb_fin->probleme_couple();
      ////sch.faire_un_pas_de_temps_pb_couple(pb_couple);
      ////}
      ////else sch.faire_un_pas_de_temps_pb_base(pb_base);

      // Mise a jour du probleme fin en temps.
      ////double pas_de_temps_ = sch.pas_de_temps();
      ////double temps_pour_mise_a_jour = pas_de_temps_+temps_courant;

      // Pour les schemas implicites, la mise a jour des pbs est comprise dans la resolution
      if(!sub_type(Schema_Implicite_base, sch))
        {
          ////if(CHT)
          ////{
          ////Probleme_Couple& pb_couple = pb_fin->probleme_couple();
          ////pb_couple.mettre_a_jour(temps_pour_mise_a_jour);
          ////}
          ////else pb_base.mettre_a_jour(temps_pour_mise_a_jour);
        }


      sch.mettre_a_jour();
      if(CHT) sch.imprimer(Cout,pb_fin);
      else sch.imprimer(Cout,pb_base);

      temps_courant = sch.temps_courant();

    } // Fin des iterations pour aller au temps.



  // On remet dans val (qui est l'adresse memoire de la condition aux limites a l'interface)
  // la bonne valeur du prolongement au temps (n+1) afin de repartir sur de bonnes
  // bases a la prochaine iteration pour etre sur de ne pas accumuler
  // les erreurs numeriques tout au long du calcul.
  P_cl.prolonger(domaine_VDF,domaine_fine,cl.frontiere(),connect.connectivites_elemF_elemG(),vit, val,nb_compo);

  // Pareil pour la temperature :
  P_cl_th.prolonger(domaine_VDF,domaine_fine,cl_th.frontiere(),connect.connectivites_elemF_elemG(),tempiotte_grosse, val_th,1);


  //////////////////////////////////////////////////////////////////////////
  //Commente suite a integration des developpements Couplage_Multi_Echelles
  //////////////////////////////////////////////////////////////////////////

  ////if(CHT)
  ////{
  ////Probleme_Couple& pb_couple = pb_fin->probleme_couple();
  ////pb_couple.traiter_postraitement();
  ////}
  ////else pb_base.traiter_postraitement();


  // 06/05/2003
  // on appelle directement la loi de paroi standard :
  // c'est elle qui va donner le frottement et non pas le probleme fin.
  Paroi_std_hyd_VDF::calculer_hyd(tab_nu_t,tab_k);

  return 1;
}


