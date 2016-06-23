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
// File:        Echange_contact_VDF_Zoom_fin.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Cond_Lim
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_VDF_Zoom_fin.h>
#include <Connectivites_faces_couple.h>
#include <Conduction.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Zone_VDF.h>

Implemente_instanciable(Echange_contact_VDF_Zoom_fin,"Paroi_Echange_contact_VDF_Zoom_fin",Echange_contact_VDF_Zoom_base);

Sortie& Echange_contact_VDF_Zoom_fin::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Echange_contact_VDF_Zoom_fin::readOn(Entree& s )
{
  return Echange_contact_VDF_Zoom_base::readOn(s);
}

void Echange_contact_VDF_Zoom_fin::mettre_a_jour(double temps)
{
  int num_faceG;
  if(sub_type(Champ_front_zoom, T_ext().valeur()))
    {
      Champ_front_zoom& ch = ref_cast(Champ_front_zoom, T_ext().valeur());
      DoubleTab& Text_valeurs = T_ext().valeur().valeurs_au_temps(temps);

      int elemG;
      REF(Pb_2G) le_pb2G;
      Pb_MG& pbMG = ch.le_pb_MG();

      int indice_pb=pbMG.indice_probleme(ch.le_pb_courant().le_nom());
      le_pb2G=pbMG.pb_2G(indice_pb);

      Connectivites_faces_couple& connections =  ref_cast(Connectivites_faces_couple,le_pb2G->connectivites().valeur());
      IntVect& connect = connections.connectivites_faceF_faceG();


      const Frontiere_dis_base& front=ch.frontiere_dis();
      const Probleme_base& pbG = ch.le_pb_exterieur();
      const Zone_dis_base& zone_disG = pbG.domaine_dis().zone_dis(0);
      const Zone_VDF& zvdfG = ref_cast(Zone_VDF, zone_disG);
      const IntTab& face_voisinsG = zvdfG.face_voisins();
      ///////////////////////////
      ///////////////////////////

      const Frontiere_dis_base& front_ext=ch.front_dis_exterieure();
      const Front_VF& front_vf_ext=ref_cast(Front_VF, front_ext);

      const Front_VF& front_vf=ref_cast(Front_VF, front);
      //const IntTab& face_voisins = zvdf_2.face_voisins();
      //A VOIR !!!!!!!!!!!!!!!!!!!!!!!!
      //ch.inconnue() sort l'inconnue du probleme grossier
      const Milieu_base& le_milieu = ch.inconnue().equation().milieu();
      int nb_comp = le_milieu.conductivite()->nb_comp();
      int i;



      h_imp_.typer("Champ_front_fonc");
      h_imp_->fixer_nb_comp(nb_comp);
      DoubleTab& tab= h_imp_->valeurs();



      DoubleVect e;

      Zone_dis_base& zone_dis1 = zone_Cl_dis().zone_dis().valeur();
      Nom nom_racc1=frontiere_dis().frontiere().le_nom();

      //MODIF
      //ATTENTION : POUR LE MOMENT, ON NE TRAITE PAS LES FRONTIERES DE TYPE "Raccord_distant_homogene"
      //DONC ON N'A PAS BESOIN DU nom_bord_oppose()
      //DE TOUTE FACON, ON NE L'A PAS DANS Camp_front_zoom
      //Nom nom_racc2=ch.nom_bord_oppose();




      //Cerr << "Domaine : "<<zone_dis1.zone().domaine().le_nom()<<finl;
      //Cerr << "On traite le raccord de nom " << nom_racc1 << finl;
      int nb_faces_raccord1 = zone_dis1.zone().raccord(nom_racc1).valeur().nb_faces();
      tab.resize(nb_faces_raccord1,nb_comp);
      //Cerr<<"nb_faces_raccord1 = "<<nb_faces_raccord1<<finl;
      //Cerr << "Nb faces raccord 1: "<<nb_faces_raccord1<<finl;



      if (zone_dis1.zone().raccord(nom_racc1).valeur().que_suis_je() =="Raccord_distant_homogene")
        {
          //POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!
          Cerr<<"POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!"<<finl;
          exit();

          // if (sub_type(Convection_Diffusion_Temperature_Turbulent,ch.equation()))
          //         {
          //           const Convection_Diffusion_Temperature_Turbulent& eq_turb
          //             = ref_cast(Convection_Diffusion_Temperature_Turbulent,ch.equation());
          //           if (eq_turb.modele_turbulence().loi_paroi_non_nulle())
          //             dequiv=1;

          //           if (dequiv)
          //             {
          //               const Turbulence_paroi& loipar =  eq_turb.modele_turbulence().loi_paroi();

          //               if (sub_type(Paroi_std_scal_hyd_VDF,loipar.valeur()))
          //                 {
          //                   Paroi_std_scal_hyd_VDF& paroi_vdf = ref_cast(Paroi_std_scal_hyd_VDF,loipar.valeur());
          //                   trace_face_raccord_distant(front_vf,paroi_vdf.tab_d_equiv(),e);
          //                 }
          //               else if (sub_type(Paroi_2couches_scal_VDF,loipar.valeur()))
          //                 {
          //                   Paroi_2couches_scal_VDF& paroi_vdf = ref_cast(Paroi_2couches_scal_VDF,loipar.valeur());
          //                   trace_face_raccord_distant(front_vf,paroi_vdf.tab_d_equiv(),e);
          //                 }
          //             }
          //           else // turbulence sans loi de paroi
          //             {
          //               DoubleVect dist;
          //               trace_face_raccord_distant(front_vf,zvdf_2.dist_norm_bord(dist,nom_racc2),e);
          //             }
          //         }
          //       else  // cas laminaire
          //         {
          //           DoubleVect dist;
          //           trace_face_raccord_distant(front_vf,zvdf_2.dist_norm_bord(dist,nom_racc2),e);
          //         }

          //       // Calcul de tab = 1/(e/lambda + 1/h_paroi)

          //       if(!sub_type(Champ_Uniforme,le_milieu.conductivite().valeur()))
          //         {
          //           //Cerr << " Raccord distant homogene et conductivite non uniforme " << finl;
          //           DoubleTab lambda;
          //           trace_raccord_distant(front_vf,zvdf_2,le_milieu.conductivite().valeurs(),lambda);
          //           if (lambda.nb_dim() == 1)
          //             {
          //               assert(nb_comp==1);
          //               for (int face=0; face<nb_faces_raccord1; face++)
          //                 {
          //                   assert(lambda(face)!=0.);
          //                   assert(h_paroi!=0.);
          //                   for(i=0; i<nb_comp; i++)
          //                     tab(face,i) = 1./(e(face)/lambda(face)+1./h_paroi);
          //                 }
          //             }
          //           else
          //             for (int face=0; face<nb_faces_raccord1; face++)
          //               for(i=0; i<nb_comp; i++){
          //                 assert(lambda(face,i)!=0.);
          //                 assert(h_paroi!=0.);
          //                 tab(face,i) = 1./(e(face)/lambda(face,i)+1./h_paroi);
          //               }
          //         }
          //       else  // la conductivite est un Champ uniforme
          //         {
          //           //Cerr << "cas d'une conductivite uniforme " << finl;
          //           const DoubleTab& lambda = le_milieu.conductivite().valeurs();
          //           for (int face=0; face<nb_faces_raccord1; face++)
          //             for(i=0; i<nb_comp; i++){
          //               assert(lambda(0,i)!=0.);
          //               assert(h_paroi!=0.);
          //               tab(face,i) = 1./(e(face)/lambda(0,i)+1./h_paroi);
          //             }
          //         }
        }        // fin du cas Raccord_distant_homogene
      else // Raccord_local_homogene
        {
          const int ndeb = front_vf.num_premiere_face();
          const int nfin = ndeb + front_vf.nb_faces();
          const int ind_prem_face_bord = zvdfG.premiere_face_bord();
          const int ndebG = front_vf_ext.num_premiere_face();
          const int nfinG = ndebG + front_vf_ext.nb_faces();

          e.resize(front_vf_ext.nb_faces());
          /////////////////
          /////////////////

          const Equation_base& eq_turb = ch.inconnue().equation();
          const RefObjU& mod = eq_turb.get_modele(TURBULENCE);

          if ((sub_type(Convection_Diffusion_Temperature,eq_turb)) && (mod.non_nul()))
            {

              //find the associated boundary
              int boundary_index=-1;
              int nb_boundaries=zvdfG.zone().nb_front_Cl();
              for (int n_bord=0; n_bord<nb_boundaries; n_bord++)
                {
                  if (zvdfG.front_VF(n_bord).le_nom() == front_vf.le_nom())
                    boundary_index=n_bord;
                }

              const Modele_turbulence_scal_base& le_mod_turb_th = ref_cast(Modele_turbulence_scal_base,mod.valeur());
              const Turbulence_paroi_scal& loipar =  le_mod_turb_th.loi_paroi();
              if (sub_type(Turbulence_paroi_scal_base,loipar.valeur()))
                {
                  const Turbulence_paroi_scal_base& paroi_vdf = ref_cast(Turbulence_paroi_scal_base,loipar.valeur());
                  for (int face=ndeb; face<nfin; face++)
                    {
                      num_faceG = connect(face);
                      // e(num_faceG-ndebG) =  paroi_vdf.d_equiv(num_faceG-ind_prem_face_bord);
                      int global_face=num_faceG-ind_prem_face_bord;
                      int local_face=zvdfG.front_VF(boundary_index).num_local_face(global_face);
                      e(global_face)=paroi_vdf.equivalent_distance(boundary_index,local_face);
                    }
                }
              //MODIF A VOIR !!!!!!!!!!
              //on calcule plusieurs fois e(num_faceG) pour la meme face !!
              else // turbulence sans loi de paroi
                {
                  for (int face=ndeb; face<nfin; face++)
                    {
                      num_faceG = connect(face);
                      e(num_faceG-ndebG) = zvdfG.dist_norm_bord(num_faceG);
                    }
                }
            }
          else  // cas laminaire
            {
              //Cerr<<"cas laminaire"<<finl;
              for (int face=ndebG; face<nfinG; face++)
                {
                  //        Cerr << " Face rempli " << face << finl;
                  e(face-ndebG) = zvdfG.dist_norm_bord(face);
                }
            }

          // Calcul de tab = 1/(e/lambda + 1/h_paroi)
          if(!sub_type(Champ_Uniforme,le_milieu.conductivite().valeur()))
            {
              //Cerr << "raccord local homogene et conductivite non uniforme" << finl;


              //ON DOIT RECUPERER LE lambda DU MILIEU GROSSIER !!
              const DoubleTab& lambda = le_milieu.conductivite().valeurs();


              assert(h_paroi!=0.);
              if (lambda.nb_dim() == 1)
                {
                  assert(nb_comp==1);
                  for (int face=ndeb; face<nfin; face++)
                    {
                      //on recupere la face grossiere correspondante
                      num_faceG = connect(face);
                      elemG =  face_voisinsG(num_faceG,0);
                      if(elemG == -1)
                        elemG =  face_voisinsG(num_faceG,1);

                      assert(lambda(elemG)!=0.);
                      for(i=0; i<nb_comp; i++)
                        tab(face-ndeb,i) = 1./(e(num_faceG-ndebG)/lambda(elemG)+1./h_paroi);
                    }
                }
              else
                for (int face=ndeb; face<nfin; face++)
                  {
                    //on recupere la face grossiere correspondante
                    num_faceG = connect(face);
                    elemG =  face_voisinsG(num_faceG,0);
                    if(elemG == -1)
                      elemG =  face_voisinsG(num_faceG,1);

                    for(i=0; i<nb_comp; i++)
                      {
                        assert(lambda(elemG,i)!=0.);
                        tab(face-ndeb,i) = 1./(e(num_faceG-ndebG)/lambda(elemG,i)+1./h_paroi);
                      }
                  }
            }
          else  // la conductivite est un Champ uniforme
            {
              //Cerr << "raccord local homogene et conductivite uniforme " << finl;
              assert(h_paroi!=0.);
              for (int face=ndeb; face<nfin; face++)
                {
                  num_faceG = connect(face);
                  //        Cerr << " Face connect " << num_faceG << finl;
                  // Cerr<<"num_faceG = "<<num_faceG<<finl;
                  //    Cerr << "Dans Contact Zoom fin h = " << 1./(e(num_faceG)/le_milieu.conductivite()(0,0)+1./h_paroi) << finl;
                  for(i=0; i<nb_comp; i++)
                    {
                      assert(le_milieu.conductivite()(0,i)!=0.);
                      tab(face-ndeb,i) = 1./(e(num_faceG-ndebG)/le_milieu.conductivite()(0,i)+1./h_paroi);
                    }
                }
            }

          //ON RECUPERE LA TEMPERATURE DU PROBLEME GROSSIER

          Text_valeurs = 0.;

          /*      const int nb_eq_totG = pbG.nombre_d_equations();
                  int trouve_eq = 0;
                  int nb_eqG = 0;

                  while((trouve_eq == 0) && (nb_eqG<nb_eq_totG))
          */
          {
            //          const Equation_base& eqG = pbG.equation(nb_eqG);
            const Equation_base& eqG = ch.equation();
            if (sub_type(Conduction, eqG) || sub_type(Convection_Diffusion_std, eqG))
              {
                //  trouve_eq = 1;
                const Champ_Inc_base& incoG = eqG.inconnue();
                const DoubleTab& valG = incoG.valeurs();
                //ON PROLONGE LA TEMPERATURE GROSSIERE
                //POUR CHAQUE EQUATION FINE
                for (int face=ndeb; face<nfin; face++)
                  {
                    num_faceG = connect(face);
                    assert(num_faceG != -1);
                    elemG =  face_voisinsG(num_faceG,0);
                    if(elemG == -1)
                      elemG =  face_voisinsG(num_faceG,1);
                    Text_valeurs(face-ndeb,0) = valG(elemG);
                    /*                      Cerr << " Text fin = " << Text_valeurs(face-ndeb,0)<< finl;
                                            Cerr << " Text fin faceG = " <<  num_faceG << finl;
                                            Cerr << " Text fin elemG = " <<  elemG << finl;
                    */
                  }

              }
            else
              {
                Cerr << " L'equation " << eqG << " n'a pas de champs inconnu temperature !!" << finl;
                exit();
              }
          }

        } // fin du cas Raccord_local_homogene
    }

  else
    {
      Cerr<<"Ce type de champ a la frontiere ne peux pas etre utilise"<<finl;
      Cerr<<"Il doit avoir acces au pb_MG et aux connectivites !!"<<finl;
      exit();
    }
  Echange_impose_base::mettre_a_jour(temps);
}



