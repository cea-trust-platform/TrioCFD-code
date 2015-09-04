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
// File:        Echange_contact_VEF_VDF_Zoom.cpp
// Directory:   $TRUST_ROOT/src/Zoom/VDF/Cond_Lim
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_VEF_VDF_Zoom.h>
#include <Echange_contact_VDF_VEF_Zoom.h>
#include <Conduction.h>
#include <Connectivites_faces_couple.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Zone_VDF.h>

Implemente_instanciable(Echange_contact_VEF_VDF_Zoom,"Contact_VEF_VDF",Temperature_imposee_paroi);









Sortie& Echange_contact_VEF_VDF_Zoom::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Echange_contact_VEF_VDF_Zoom::readOn(Entree& s )
{
  return Temperature_imposee_paroi::readOn(s);
}

void Echange_contact_VEF_VDF_Zoom::mettre_a_jour(double temps)
{
  if(sub_type(Champ_front_zoom, le_champ_front.valeur()))
    {
      Champ_front_zoom& ch = ref_cast(Champ_front_zoom, le_champ_front.valeur());
      DoubleTab& Text_valeurs = ch.valeurs_au_temps(temps);
      int elemF;
      int  elemG;
      REF(Pb_2G) le_pb2G;
      Pb_MG& pbMG = ch.le_pb_MG();
      const double h_paroi=1.e30;

      const Probleme_base& pbG = ch.le_pb_exterieur();
      const Probleme_base& pbF = ch.le_pb_courant();
      int indice_pb=pbMG.indice_probleme(pbF.le_nom());
      le_pb2G=pbMG.pb_2G(indice_pb);

      const Zone_dis_base& zone_disG = pbG.domaine_dis().zone_dis(0);
      const Zone_VDF& zvdf = ref_cast(Zone_VDF, zone_disG); // ce doit etre en principe une zone VDF !!!
      const IntTab& face_voisinsG = zvdf.face_voisins();
      const DoubleVect& surfacesVDF = zvdf.face_surfaces();

      const Zone_dis_base& zone_disF = pbF.domaine_dis().zone_dis(0);
      const Zone_VEF& zvef = ref_cast(Zone_VEF, zone_disF); // ce doit etre en principe une zone VEF !!!
      const IntTab& face_voisinsF = zvef.face_voisins();


      const Connectivites_faces_couple& connections =  ref_cast(Connectivites_faces_couple,le_pb2G->connectivites().valeur());
      const IntVect& connect = connections.connectivites_faceF_faceG();

      const Frontiere_dis_base& front_ext=ch.front_dis_exterieure();

      const Frontiere_dis_base& front=ch.frontiere_dis();

      int face;

      const Front_VF& front_vf_ext=ref_cast(Front_VF, front_ext);
      const Front_VF& front_vf=ref_cast(Front_VF, front);
      //const IntTab& face_voisins = zvdf_2.face_voisins();



      REF(Milieu_base) le_milieu;
      le_milieu = pbF.milieu();
      REF(Milieu_base) le_milieu_vdf;
      le_milieu_vdf = pbG.milieu();

      const int nb_comp_vef = le_milieu->conductivite()->nb_comp();
      const int nb_comp_vdf = le_milieu_vdf->conductivite()->nb_comp();
      int i;


      const Zone_dis_base& zone_dis1 = zone_Cl_dis().zone_dis().valeur();
      const Nom nom_racc1=frontiere_dis().frontiere().le_nom();



      const int nb_faces_front_vdf = front_vf_ext.nb_faces();
      const int num_prem_face_vdf = front_vf_ext.num_premiere_face();
      const int nb_faces_front_vef = front_vf.nb_faces();
      const int num_prem_face_vef = front_vf.num_premiere_face();


      //      assert(nb_faces_front_vdf == nb_faces_front_vef);


      DoubleTab h_vdf;
      DoubleTab h_vef;
      DoubleVect e;

      e.resize(nb_faces_front_vdf);

      h_vef.resize(nb_faces_front_vef,nb_comp_vef);
      h_vef = 0.;

      h_vdf.resize(nb_faces_front_vdf,nb_comp_vdf);
      h_vdf = 0.;





      if (zone_dis1.zone().raccord(nom_racc1).valeur().que_suis_je() =="Raccord_distant_homogene")
        {
          //POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!
          Cerr<<"POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!"<<finl;
          exit();
        }   // fin du cas Raccord_distant_homogene
      else // Raccord_local_homogene
        {

          const int ndeb = num_prem_face_vef;
          const int nfin = ndeb + nb_faces_front_vef;
          int num_faceG;

          const Equation_base& eq_turb = ch.inconnue().equation();
          const RefObjU& mod = eq_turb.get_modele(TURBULENCE);

          if ((sub_type(Convection_Diffusion_Temperature,eq_turb)) && (mod.non_nul()))
            {

              //find the associated boundary
              int boundary_index=-1;
              int nb_boundaries=zvdf.zone().nb_front_Cl();
              for (int n_bord=0; n_bord<nb_boundaries; n_bord++)
                {
                  if (zvdf.front_VF(n_bord).le_nom() == front_vf.le_nom())
                    boundary_index=n_bord;
                }

              const Modele_turbulence_scal_base& le_mod_turb_th = ref_cast(Modele_turbulence_scal_base,mod.valeur());
              const Turbulence_paroi_scal& loipar =  le_mod_turb_th.loi_paroi();
              if (sub_type(Turbulence_paroi_scal_base,loipar.valeur()))
                {
                  const Turbulence_paroi_scal_base& paroi_vdf = ref_cast(Turbulence_paroi_scal_base,loipar.valeur());
                  for (face=ndeb; face<nfin; face++)
                    {
                      num_faceG = connect(face);
                      //e(num_faceG-num_prem_face_vdf) =  paroi_vdf.d_equiv(num_faceG-num_prem_face_vdf);
                      int global_face=num_faceG-num_prem_face_vdf;
                      int local_face=zvdf.front_VF(boundary_index).num_local_face(global_face);
                      e(global_face)=paroi_vdf.equivalent_distance(boundary_index,local_face);
                    }
                }
              else // turbulence sans loi de paroi
                {
                  const int nfinG = num_prem_face_vdf+nb_faces_front_vdf;
                  for (face=num_prem_face_vdf; face<nfinG; face++)
                    {
                      e(face-num_prem_face_vdf) = zvdf.dist_norm_bord(face);
                    }
                }

            }
          else  // cas laminaire
            {
              const int nfinG = num_prem_face_vdf+nb_faces_front_vdf;
              //Cerr<<"cas laminaire"<<finl;
              for (face=num_prem_face_vdf; face<nfinG; face++)
                {
                  e(face-num_prem_face_vdf) = zvdf.dist_norm_bord(face);
                }
            }


          // Calcul de tab = 1/(e/lambda + 1/h_paroi)
          if(!sub_type(Champ_Uniforme,le_milieu->conductivite().valeur()))
            {
              //Cerr << "raccord local homogene et conductivite non uniforme" << finl;
              const DoubleTab& lambda = le_milieu->conductivite().valeurs();
              const DoubleTab& lambda_vdf = le_milieu_vdf->conductivite().valeurs();
              assert(h_paroi!=0.);

              if (lambda.nb_dim() == 1)
                {
                  assert(nb_comp_vef==1);
                  for (face=ndeb; face<nfin; face++)
                    {
                      //on recupere la face grossiere correspondante
                      num_faceG = connect(face);
                      elemF =  face_voisinsF(face,0);
                      if(elemF == -1)
                        elemF =  face_voisinsF(face,1);

                      elemG =  face_voisinsG(num_faceG,0);
                      if(elemG == -1)
                        elemG =  face_voisinsG(num_faceG,1);

                      assert(lambda(elemF)!=0.);
                      for(i=0; i<nb_comp_vef; i++)
                        {
                          h_vef(face-ndeb,i) = 1./(1./pdt_scalSqrt(zvef,face,face,elemF,dimension,lambda(elemF))+1./h_paroi) ;
                        }
                      for(i=0; i<nb_comp_vdf; i++)
                        {
                          h_vdf(num_faceG-num_prem_face_vdf,i) = 1./(e(num_faceG-num_prem_face_vdf)/lambda_vdf(elemG)+1./h_paroi) ;
                        }
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
                    elemG =  face_voisinsG(num_faceG,0);
                    if(elemG == -1)
                      elemG =  face_voisinsG(num_faceG,1);

                    for(i=0; i<nb_comp_vef; i++)
                      {
                        h_vef(face-ndeb,i) = 1./(1./pdt_scalSqrt(zvef,face,face,elemF,dimension,lambda(elemF,i))+1./h_paroi) ;
                      }
                    for(i=0; i<nb_comp_vdf; i++)
                      {
                        assert(lambda_vdf(elemF,i)!=0.);
                        h_vdf(num_faceG-num_prem_face_vdf,i) = 1./(e(num_faceG-num_prem_face_vdf)/lambda_vdf(elemG,i)+1./h_paroi) ;
                      }
                  }
            }
          else  // la conductivite est un Champ uniforme
            {
              assert(h_paroi!=0.);
              const DoubleTab& lambda = le_milieu->conductivite().valeurs();
              const DoubleTab& lambda_vdf = le_milieu_vdf->conductivite().valeurs();
              for(i=0; i<nb_comp_vef; i++)
                assert(lambda(0,i)!=0.); // juste des asserts

              for (face=ndeb; face<nfin; face++)
                {
                  num_faceG = connect(face);
                  elemF =  face_voisinsF(face,0);
                  if(elemF == -1)
                    elemF =  face_voisinsF(face,1);
                  elemG =  face_voisinsG(num_faceG,0);
                  if(elemG == -1)
                    elemG =  face_voisinsG(num_faceG,1);

                  assert(num_faceG != -1);

                  for(i=0; i<nb_comp_vef; i++)
                    {
                      h_vef(face-ndeb,i) = 1./(1./pdt_scalSqrt(zvef,face,face,elemF,dimension,lambda(0,i))+1./h_paroi);
                    }
                  for(i=0; i<nb_comp_vdf; i++)
                    {
                      h_vdf(num_faceG-num_prem_face_vdf,i) = 1./(e(num_faceG-num_prem_face_vdf)/lambda_vdf(0,i)+1./h_paroi) ;
                    }
                }
            }


          // Calcul de Text
          assert(nb_comp_vef == 1); // on sait pas faire pour l'instant sinon!!
          assert(nb_comp_vdf == 1);
          if (nb_comp_vef != 1 || nb_comp_vdf != 1)
            {
              Cerr << "On ne sait pas traiter un pb couple 'zoom' lorsque nb_composante de la conductivite est plus grand que 1 " << finl;
              exit();
            }



          //ON RECUPERE LA TEMPERATURE FINE
          {
            int face2, face3, face4;
            const Equation_base& eqG = ch.equation();
            const int nb_eqF = pbF.nombre_d_equations();
            int nb_equationF = 0;
            int trouve_eq = 0;

            while((trouve_eq == 0) && (nb_equationF<nb_eqF))
              {
                const  Equation_base& eqF = pbF.equation(nb_equationF);

                if (sub_type(Conduction, eqF) || sub_type(Convection_Diffusion_std, eqF))
                  {
                    const Champ_Inc_base& incoG = eqG.inconnue();
                    const DoubleTab& T_VDF = incoG.valeurs();
                    const DoubleTab& T_VEF = eqF.inconnue().valeurs();
                    DoubleVect cumulHiSi(nb_faces_front_vdf);
                    DoubleVect Tparoi(nb_faces_front_vdf);

                    cumulHiSi = 0.;
                    Tparoi = 0.;
                    trouve_eq = 1;


                    //ON PROLONGE LA TEMPERATURE GROSSIERE
                    //POUR CHAQUE EQUATION FINE
                    for (face=ndeb; face<nfin; face++)
                      {
                        double tempo = 0.;
                        num_faceG = connect(face);
                        assert(num_faceG != -1);
                        elemF =  face_voisinsF(face,0);
                        if(elemF == -1)
                          elemF =  face_voisinsF(face,1);


                        // on recupere les faces 2, 3 et 4 de l'element contenant "face"
                        face2 = zvef.elem_faces(elemF,0);
                        face3 = zvef.elem_faces(elemF,1);
                        if (face2 == face)
                          face2 = zvef.elem_faces(elemF,dimension);
                        else if (face3 == face)
                          face3 = zvef.elem_faces(elemF,dimension);

                        cumulHiSi(num_faceG-num_prem_face_vdf)+=surfacesVEF(zvef,face,dimension)*h_vef(face-ndeb,0);
                        tempo += pdt_scal(zvef,face,face2,elemF,dimension,1.)*T_VEF(face2) + pdt_scal(zvef,face,face3,elemF,dimension,1.)*T_VEF(face3);

                        // Cerr << "tempo= " << tempo << finl;
                        if(dimension == 3)
                          {
                            face4 = zvef.elem_faces(elemF,2);
                            if (face4 == face)
                              face4 = zvef.elem_faces(elemF,dimension);
                            tempo += pdt_scal(zvef,face,face4,elemF,dimension,1.)*T_VEF(face4);
                          }
                        //       Cerr << "face " << face << " " << face2 << " " << face3 << " " << face4 << finl;
                        //  Cerr << "tempo 2 = " << tempo << finl;
                        /*                     Cerr << "psc1  = " << le_pdt_scal(zvef,face,face2,elemF,dimension,1.) << finl;
                                               Cerr << "psc2 = " << le_pdt_scal(zvef,face,face3,elemF,dimension,1.) << finl;
                                               Cerr << "psc3 = " << le_pdt_scal(zvef,face,face4,elemF,dimension,1.) << finl;
                        */                     tempo *= -surfacesVEF(zvef,face,dimension)*h_vef(face-ndeb,0)/pdt_scal(zvef,face,face,elemF,dimension,1.);
                        Tparoi(num_faceG-num_prem_face_vdf)+=tempo;
                      }

                    for (face=0; face<nb_faces_front_vdf; face++)
                      {
                        // Cerr << "Tparoi= " << Tparoi(face) << finl;
                        const int face_glob = face+num_prem_face_vdf;
                        elemG =  face_voisinsG(face_glob,0);
                        if(elemG == -1)
                          elemG =  face_voisinsG(face_glob,1);

                        const double hS = surfacesVDF(face_glob)*h_vdf(face,0);
                        Tparoi(face) += (hS*T_VDF(elemG));
                        Tparoi(face) /= (hS+cumulHiSi(face));
                      }
                    for (face=ndeb; face<nfin; face++)
                      {
                        num_faceG = connect(face);
                        Text_valeurs(face-ndeb,0) = Tparoi(num_faceG-num_prem_face_vdf);
                        // Cerr << "Text = " << Text_valeurs(face-ndeb,0) << finl;
                      }
                    double flux_tot=0;
                    for (face=0; face<nb_faces_front_vdf; face++)
                      {
                        // Cerr << "Tparoi= " << Tparoi(face) << finl;
                        const int face_glob = face+num_prem_face_vdf;
                        elemG =  face_voisinsG(face_glob,0);
                        if(elemG == -1)
                          elemG =  face_voisinsG(face_glob,1);
                        flux_tot+=h_vdf(face,0)*surfacesVDF(face_glob)*(Tparoi(face)-T_VDF(elemG));
                      }
                    // Cerr << "flux tot = " << flux_tot << finl;
                  }
                nb_equationF++;
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

}



