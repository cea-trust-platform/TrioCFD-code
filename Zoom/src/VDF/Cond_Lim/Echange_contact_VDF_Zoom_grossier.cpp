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
// File:        Echange_contact_VDF_Zoom_grossier.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Cond_Lim
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_VDF_Zoom_grossier.h>
#include <Conduction.h>
#include <Connectivites_faces_couple.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Zone_VDF.h>
#include <Ref_Milieu_base.h>

Implemente_instanciable(Echange_contact_VDF_Zoom_grossier,"Paroi_Echange_contact_VDF_Zoom_grossier",Echange_contact_VDF_Zoom_base);

Sortie& Echange_contact_VDF_Zoom_grossier::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Echange_contact_VDF_Zoom_grossier::readOn(Entree& s )
{
  return Echange_contact_VDF_Zoom_base::readOn(s);
}

void Echange_contact_VDF_Zoom_grossier::mettre_a_jour(double temps)
{
  int num_faceG;
  if(sub_type(Champ_front_zoom, T_ext().valeur()))
    {
      //   Cout << "Dans Echange_contact_VDF_Zoom_grossier Text = " << T_ext().valeur().valeurs() << finl;
      Champ_front_zoom& ch = ref_cast(Champ_front_zoom, T_ext().valeur());
      DoubleTab& Text_valeurs = T_ext().valeur().valeurs_au_temps(temps);
      int elemF;
      REF(Pb_2G) le_pb2G;
      Pb_MG& pbMG = ch.le_pb_MG();

      //POUR LE MOMENT, CONSIDERONS QUE LE PB GROSSIER A UN SEUL PB FIN !!!!



      const Probleme_base& pbF = ch.le_pb_exterieur();
      int indice_pb=pbMG.indice_probleme(pbF.le_nom());
      le_pb2G=pbMG.pb_2G(indice_pb);

      const Zone_dis_base& zone_disF = pbF.domaine_dis();
      const Zone_VDF& zvdfF = ref_cast(Zone_VDF, zone_disF);
      const DoubleVect& surfacesF = zvdfF.face_surfaces();
      const IntTab& face_voisinsF = zvdfF.face_voisins();


      const Connectivites_faces_couple& connections =  ref_cast(Connectivites_faces_couple,le_pb2G->connectivites().valeur());
      const IntVect& connect = connections.connectivites_faceF_faceG();

      const Frontiere_dis_base& front_ext=ch.front_dis_exterieure();

      const Frontiere_dis_base& front=ch.frontiere_dis();

      int nb_face_frontG;

      //    Cerr << "Dans Echange_contact_VDF_Zoom_grossier frontier fine = " << front.le_nom() << finl;

      ///////////////////////////
      ///////////////////////////
      int face;


      //MODIF
      //const Zone_VDF& zvdf_2=ref_cast(Zone_VDF, ch.zone_dis());
      const Front_VF& front_vf_ext=ref_cast(Front_VF, front_ext);
      const Front_VF& front_vf=ref_cast(Front_VF, front);
      //const IntTab& face_voisins = zvdf_2.face_voisins();



      //ON VEUT LE MILIEU FIN !!!!!!!!!
      //const Milieu_base& le_milieu = ch.inconnue().equation().milieu();
      REF(Champ_base) ch_tampon;
      REF(Milieu_base) le_milieu;
      le_milieu = pbF.milieu();

      const int nb_comp = le_milieu->conductivite()->nb_comp();
      int i;



      h_imp_.typer("Champ_front_fonc");
      h_imp_->fixer_nb_comp(nb_comp);
      DoubleTab& tab= h_imp_->valeurs();



      DoubleVect e;

      const Zone_dis_base& zone_dis1 = zone_Cl_dis().zone_dis().valeur();
      const Nom nom_racc1=frontiere_dis().frontiere().le_nom();





      //Cerr << "Zone : "<<zone_dis1.zone().domaine().le_nom()<<finl;
      //Cerr << "On traite le raccord de nom " << nom_racc1 << finl;
      const int nb_faces_front_gros = front_vf.nb_faces();
      const int num_prem_faceG = front_vf.num_premiere_face();
      //Cerr<<" num_prem_faceG = "<<num_prem_faceG<<finl;
      DoubleVect cumulSurfaces(nb_faces_front_gros);
      cumulSurfaces = 0.;
      DoubleVect cumulHiSi(nb_faces_front_gros);
      cumulHiSi = 0.;
      DoubleVect valHiSi(front_vf_ext.nb_faces());
      valHiSi = 0.;
      tab.resize(nb_faces_front_gros,nb_comp);
      tab = 0.;

      //      Cerr<<"Nb faces grossieres = "<<nb_faces_front_gros << finl;
      //  Cerr << "Nb faces fines : "<<front_vf_ext.nb_faces() << finl;



      if (zone_dis1.zone().raccord(nom_racc1).valeur().que_suis_je() =="Raccord_distant_homogene")
        {
          //POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!
          Cerr<<"POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!"<<finl;
          exit();
          //MODIF
          //ATTENTION : POUR LE MOMENT, ON NE TRAITE PAS LES FRONTIERES DE TYPE "Raccord_distant_homogene"
          //DONC ON N'A PAS BESOIN DU nom_bord_oppose()
          //DE TOUTE FACON, ON NE L'A PAS DANS Camp_front_zoom
          //Nom nom_racc2=ch.nom_bord_oppose();

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
          //               for (int face=0; face<nb_faces_front_gros; face++)
          //                 {
          //                   assert(lambda(face)!=0.);
          //                   assert(h_paroi!=0.);
          //                   for(i=0; i<nb_comp; i++)
          //                     tab(face,i) = 1./(e(face)/lambda(face)+1./h_paroi);
          //                 }
          //             }
          //           else
          //             for (int face=0; face<nb_faces_front_gros; face++)
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
          //           for (int face=0; face<nb_faces_front_gros; face++)
          //             for(i=0; i<nb_comp; i++){
          //               assert(lambda(0,i)!=0.);
          //               assert(h_paroi!=0.);
          //               tab(face,i) = 1./(e(face)/lambda(0,i)+1./h_paroi);
          //             }
          //         }
        }        // fin du cas Raccord_distant_homogene
      else // Raccord_local_homogene
        {

          //MODIF!!!!!!!!!!
          const int ndeb = front_vf_ext.num_premiere_face();
          const int nfin = ndeb + front_vf_ext.nb_faces();
          const int ind_prem_face_bord = zvdfF.premiere_face_bord();




          //MODIF !!!
          //e.resize(front_vf.nb_faces());
          //TROP GRAND !!!!!!!!
          e.resize(front_vf_ext.nb_faces());
          /////////////////
          /////////////////

          const Equation_base& eq_turb = ch.inconnue().equation();
          const RefObjU& mod = eq_turb.get_modele(TURBULENCE);

          if ((sub_type(Convection_Diffusion_Temperature,eq_turb)) && (mod.non_nul()))
            {

              //find the associated boundary
              int boundary_index=-1;
              int nb_boundaries=zvdfF.zone().nb_front_Cl();
              for (int n_bord=0; n_bord<nb_boundaries; n_bord++)
                {
                  if (zvdfF.front_VF(n_bord).le_nom() == front_vf.le_nom())
                    boundary_index=n_bord;
                }

              const Modele_turbulence_scal_base& le_mod_turb_th = ref_cast(Modele_turbulence_scal_base,mod.valeur());
              const Turbulence_paroi_scal& loipar =  le_mod_turb_th.loi_paroi();
              if (sub_type(Turbulence_paroi_scal_base,loipar.valeur()))
                {
                  {
                    const Turbulence_paroi_scal_base& paroi_vdf = ref_cast(Turbulence_paroi_scal_base,loipar.valeur());
                    for (int faceb=ndeb; faceb<nfin; faceb++)
                      {
                        num_faceG = connect(faceb);
                        if(num_faceG != -1)
                          {
                            //e(faceb-ndeb) =  paroi_vdf.d_equiv(faceb-ind_prem_face_bord);
                            int global_face=faceb-ind_prem_face_bord;
                            int local_face=zvdfF.front_VF(boundary_index).num_local_face(global_face);
                            e(local_face)=paroi_vdf.equivalent_distance(boundary_index,local_face);
                          }
                        else
                          {
                            Cerr << "Une face fine de l'interface ne correspond a aucune face grossiere !!!" << finl;
                            exit();
                          }
                      }
                  }
                }
              //MODIF A VOIR !!!!!!!!!!
              else // turbulence sans loi de paroi
                {
                  Cerr<<"turbulence sans loi de paroi"<<finl;
                  for (face=ndeb; face<nfin; face++)
                    {
                      num_faceG = connect(face);
                      if(num_faceG != -1)
                        {
                          e(face-ndeb) = zvdfF.dist_norm_bord(face);
                          //   Cerr<<"num_faceF = "<<face<<finl;
                          //                           Cerr<<"e(face) = "<<e(face)<<finl;
                        }
                      else
                        {
                          Cerr << "Une face fine de l'interface ne correspond a aucune face grossiere !!!" << finl;
                          exit();
                        }
                    }
                }
            }
          else  // cas laminaire
            {
              //Cerr<<"cas laminaire"<<finl;
              for (face=ndeb; face<nfin; face++)
                {
                  num_faceG = connect(face);
                  if(num_faceG != -1)
                    {
                      e(face-ndeb) = zvdfF.dist_norm_bord(face);
                      //   Cerr<<"num_faceF = "<<face<<finl;
                      //                       Cerr<<"e(face) = "<<e(face)<<finl;
                    }
                  else
                    {
                      Cerr << "Une face fine de l'interface ne correspond a aucune face grossiere !!!" << finl;
                      exit();
                    }
                }
            }

          // Calcul de tab = 1/(e/lambda + 1/h_paroi)
          if(!sub_type(Champ_Uniforme,le_milieu->conductivite().valeur()))
            {
              //Cerr << "raccord local homogene et conductivite non uniforme" << finl;
              const DoubleTab& lambda = le_milieu->conductivite().valeurs();
              assert(h_paroi!=0.);

              if (lambda.nb_dim() == 1)
                {
                  assert(nb_comp==1);
                  for (face=ndeb; face<nfin; face++)
                    {
                      //on recupere la face grossiere correspondante
                      num_faceG = connect(face);
                      elemF =  face_voisinsF(face,0);
                      if(elemF == -1)
                        elemF =  face_voisinsF(face,1);

                      // int elem = face_voisins(face,0);
                      //                   if (elem == -1)
                      //                     elem = face_voisins(face,1);
                      assert(lambda(elemF)!=0.);
                      cumulSurfaces(num_faceG-num_prem_faceG) += surfacesF(face);
                      valHiSi(face-ndeb) = surfacesF(face)/(e(face-ndeb)/lambda(elemF)+1./h_paroi);
                      cumulHiSi(num_faceG-num_prem_faceG) += valHiSi(face-ndeb);
                      for(i=0; i<nb_comp; i++)
                        {
                          tab(num_faceG-num_prem_faceG,i) = tab(num_faceG-num_prem_faceG,i) + valHiSi(face);
                        }
                    }
                  for(nb_face_frontG = 0; nb_face_frontG<nb_faces_front_gros; nb_face_frontG++)
                    {
                      for(i=0; i<nb_comp; i++)
                        if (cumulSurfaces(nb_face_frontG)>0.) tab(nb_face_frontG,i) = tab(nb_face_frontG,i)/ cumulSurfaces(nb_face_frontG);
                    }
                }
              else
                for (face=ndeb; face<nfin; face++)
                  {
                    //on recupere la face grossiere correspondante
                    num_faceG = connect(face);
                    elemF =  face_voisinsF(face,0);
                    if(elemF == -1)
                      elemF =  face_voisinsF(face,1);

                    // int elem = face_voisins(face,0);
                    //                 if (elem == -1)
                    //                   elem = face_voisins(face,1);
                    cumulSurfaces(num_faceG-num_prem_faceG) += surfacesF(face);
                    valHiSi(face-ndeb) = surfacesF(face)/(e(face-ndeb)/lambda(elemF,0)+1./h_paroi);
                    cumulHiSi(num_faceG-num_prem_faceG) += valHiSi(face-ndeb);
                    for(i=0; i<nb_comp; i++)
                      {
                        assert(lambda(elemF,i)!=0.);
                        tab(num_faceG-num_prem_faceG,i) = tab(num_faceG-num_prem_faceG,i) + surfacesF(face)/(e(face-ndeb)/lambda(elemF,i)+1./h_paroi);
                      }
                  }
              for(nb_face_frontG = 0; nb_face_frontG<nb_faces_front_gros; nb_face_frontG++)
                {
                  for(i=0; i<nb_comp; i++)
                    if (cumulSurfaces(nb_face_frontG)>0.) tab(nb_face_frontG,i) = tab(nb_face_frontG,i)/ cumulSurfaces(nb_face_frontG);
                }
            }
          else  // la conductivite est un Champ uniforme
            {
              assert(h_paroi!=0.);
              const DoubleTab& lambda = le_milieu->conductivite().valeurs();
              for(i=0; i<nb_comp; i++)
                assert(lambda(0,i)!=0.); // juste des asserts

              //Cerr << "raccord local homogene et conductivite uniforme " << finl;
              for (face=ndeb; face<nfin; face++)
                {
                  num_faceG = connect(face);
                  // Cerr<<"faceF = "<<face<<finl;
                  //                   Cerr<<"num_faceG = "<<num_faceG<<finl;
                  //                   Cerr<<"num_faceG-num_prem_faceG = "<<num_faceG-num_prem_faceG<<finl;
                  // if (num_faceG != -1)
                  assert(num_faceG != -1);
                  {
                    //               Cerr<<"num_faceG = "<<num_faceG<<finl;
                    //                       Cerr<<"num_faceG-num_prem_faceG = "<<num_faceG-num_prem_faceG<<finl;
                    cumulSurfaces(num_faceG-num_prem_faceG) += surfacesF(face);
                    valHiSi(face-ndeb) = surfacesF(face)/(e(face-ndeb)/lambda(0,0)+1./h_paroi);
                    cumulHiSi(num_faceG-num_prem_faceG) += valHiSi(face-ndeb);
                    //   Cerr << "Dans Contact Zoom gros h = " << 1./(e(face-ndeb)/lambda(0,0)+1./h_paroi) << finl;
                    //   Cerr << "Dans Contact Zoom gros valHiSi(face-ndeb) = " << valHiSi(face-ndeb) << finl;

                    for(i=0; i<nb_comp; i++)
                      {
                        //         Cerr<<"faceF : "<<face<<finl;
                        //                         Cerr<<"e(face) = "<<e(face)<<finl;
                        //                         Cerr<<"le_milieu->conductivite()(0,i) = "<<le_milieu->conductivite()(0,i)<<finl;
                        tab(num_faceG-num_prem_faceG,i) = tab(num_faceG-num_prem_faceG,i) +surfacesF(face)/(e(face-ndeb)/lambda(0,i)+1./h_paroi) ;
                      }
                  }
                }
              for(nb_face_frontG = 0; nb_face_frontG<nb_faces_front_gros; nb_face_frontG++)
                {
                  for(i=0; i<nb_comp; i++)
                    if (cumulSurfaces(nb_face_frontG)>0.) tab(nb_face_frontG,i) = tab(nb_face_frontG,i)/ cumulSurfaces(nb_face_frontG);
                }
            }


          //Calcul de Text = (Somme_i h_i*S_i*T_i)/(Somme_i h_i*S_i)
          assert(nb_comp == 1); // on sait pas faire pour l'instant sinon!!
          if (nb_comp != 1)
            {
              Cerr << "On ne sait pas traiter un pb couple 'zoom' lorsque nb_composante de la conductivite est plus grand que 1 " << finl;
              exit();
            }
          //ON RECUPERE LA TEMPERATURE FINE
          Text_valeurs = 0.;

          /*    const int nb_eqF = pbF.nombre_d_equations();
                int trouve_eq = 0;
                int nb_equationF = 0;
                while((trouve_eq == 0) && (nb_equationF<nb_eqF))*/
          {
            const  Equation_base& eqF = ch.equation();
            if (sub_type(Conduction, eqF) || sub_type(Convection_Diffusion_std, eqF))
              {
                //   trouve_eq = 1;
                const Champ_Inc_base& incoF = eqF.inconnue();
                const DoubleTab& valF = incoF.valeurs();
                for (face=ndeb; face<nfin; face++)
                  {
                    num_faceG = connect(face);
                    assert(num_faceG != -1);
                    elemF =  face_voisinsF(face,0);
                    if(elemF == -1)
                      elemF =  face_voisinsF(face,1);
                    Text_valeurs(num_faceG-num_prem_faceG,0) += valF(elemF)*valHiSi(face-ndeb);
                  }


                for(nb_face_frontG = 0; nb_face_frontG<nb_faces_front_gros; nb_face_frontG++)
                  {
                    //    Cout << "cumulHiSi(nb_face_frontG)" << cumulHiSi(nb_face_frontG) << finl;
                    if (cumulHiSi(nb_face_frontG)>0.) Text_valeurs(nb_face_frontG,0) =  Text_valeurs(nb_face_frontG,0) /(cumulHiSi(nb_face_frontG));
                    //  Cerr << " Text gros = " << Text_valeurs(nb_face_frontG,0)<< finl;
                  }
              }
            else
              {
                Cerr << " L'equation " << eqF << " n'a pas de champs inconnu temperature !!" << finl;
                exit();
              }

          }

          //Cout << "Text valeurs " << Text_valeurs << finl;
          //Cout << "size " <<nb_faces_front_gros  << finl;






        } // fin du cas Raccord_local_homogene
      /*        }
                else
                {
                Cerr<<"On ne peut pas faire de couplage avec ce type de connectivites !!!"<<finl;
                exit();
                }*/
    }
  else
    {
      Cerr<<"Ce type de champ a la frontiere ne peut pas etre utilise"<<finl;
      Cerr<<"Il doit avoir acces au pb_MG et aux connectivites !!"<<finl;
      exit();
    }

  Echange_impose_base::mettre_a_jour(temps);
  //  Cerr<<"fin de Echange_contact_VDF_Zoom_grossier::mettre_a_jour"<<finl;
}



