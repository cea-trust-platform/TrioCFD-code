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
// File:        Echange_contact_VDF_VEF_Zoom.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Cond_Lim
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_VDF_VEF_Zoom.h>
#include <Conduction.h>
#include <Connectivites_faces_couple.h>
#include <Convection_Diffusion_std.h>
#include <Milieu_base.h>
#include <TRUST_Ref.h>


Implemente_instanciable(Echange_contact_VDF_VEF_Zoom,"Contact_VDF_VEF",Echange_contact_VDF_Zoom_base);

// Fonctions hors classe pour le pdt scal n_i.n_j
// Attention : la normal a la face num_face est suppose sortante
// tandis que celle a la face num2 est reorientee dans la methode
// afin d etre sortante.
double pdt_scal(const Domaine_VEF& le_dom,int num_face,int num2,
                int num_elem,int dimension,double diffu)
{
  double pscal;
  const IntTab& face_voisinsF = le_dom.face_voisins();

  pscal = le_dom.face_normales(num_face,0)*le_dom.face_normales(num2,0)
          + le_dom.face_normales(num_face,1)*le_dom.face_normales(num2,1);
  if (dimension == 3)
    pscal += le_dom.face_normales(num_face,2)*le_dom.face_normales(num2,2);

  //if (!(face_voisinsF(num2,0) == num_elem)) pscal = -pscal;
  if ( (face_voisinsF(num_face,0) == face_voisinsF(num2,0)) ||
       (face_voisinsF(num_face,1) == face_voisinsF(num2,1))) pscal = -pscal;

  return (pscal*diffu)/le_dom.volumes(num_elem);
}


// Attention : la normal a la face num_face est suppose sortante
// tandis que celle a la face num2 est reorientee dans la methode
// afin d etre sortante.
double pdt_scalSqrt(const Domaine_VEF& le_dom,int num_face,int num2,
                    int num_elem,int dimension,double diffu)
{
  double pscal;
  //const IntTab& face_voisinsF = le_dom.face_voisins();

  pscal = le_dom.face_normales(num_face,0)*le_dom.face_normales(num2,0)
          + le_dom.face_normales(num_face,1)*le_dom.face_normales(num2,1);
  if (dimension == 3)
    pscal += le_dom.face_normales(num_face,2)*le_dom.face_normales(num2,2);
  pscal = sqrt(std::fabs(pscal));
  //if (!(face_voisinsF(num2,0) == num_elem)) pscal = -pscal;
  //  if ( (face_voisinsF(num_face,0) == face_voisinsF(num2,0)) ||
  //     (face_voisinsF(num_face,1) == face_voisinsF(num2,1))) pscal = -pscal;

  return (pscal*diffu)/le_dom.volumes(num_elem);
}

double surfacesVEF(const Domaine_VEF& le_dom,int num_face,int dimension)
{
  double pscal;

  pscal = le_dom.face_normales(num_face,0)*le_dom.face_normales(num_face,0)
          + le_dom.face_normales(num_face,1)*le_dom.face_normales(num_face,1);
  if (dimension == 3)
    pscal += le_dom.face_normales(num_face,2)*le_dom.face_normales(num_face,2);
  pscal = sqrt(pscal);

  return pscal;
}


int Echange_contact_VDF_VEF_Zoom::compatible_avec_discr(const Discretisation_base& discr) const
{
  return 1;
}


Sortie& Echange_contact_VDF_VEF_Zoom::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Echange_contact_VDF_VEF_Zoom::readOn(Entree& s )
{

  s >> le_champ_front;

  h_paroi=1.e30;
  // GF on le type plus tard
  //h_imp_.typer("Champ_front_fonc");

  return s ;
}


/**
 * VEF toujours le probleme fin !!
 */
void Echange_contact_VDF_VEF_Zoom::mettre_a_jour(double temps)
{
  int num_faceG;
  if(sub_type(Champ_front_zoom, T_ext()))
    {
      //   Cout << "Dans Echange_contact_VDF_VEF_Zoom Text = " << T_ext().valeurs() << finl;
      Champ_front_zoom& ch = ref_cast(Champ_front_zoom, T_ext());
      DoubleTab& Text_valeurs = T_ext().valeurs_au_temps(temps);
      int elemF;
      OBS_PTR(Pb_2G) le_pb2G;
      Pb_MG& pbMG = ch.le_pb_MG();

      const Probleme_base& pbF = ch.le_pb_exterieur();
      int indice_pb=pbMG.indice_probleme(pbF.le_nom());
      le_pb2G=pbMG.pb_2G(indice_pb);

      const Domaine_dis_base& domaine_disF = pbF.domaine_dis();
      const Domaine_VEF& zvef = ref_cast(Domaine_VEF, domaine_disF); // ce doit etre en principe un domaine VEF !!!
      const IntTab& face_voisinsF = zvef.face_voisins();


      const Connectivites_faces_couple& connections = ref_cast(Connectivites_faces_couple,le_pb2G->connectivites().valeur());
      const IntVect& connect = connections.connectivites_faceF_faceG();

      const Frontiere_dis_base& front_ext=ch.front_dis_exterieure();

      const Frontiere_dis_base& front=ch.frontiere_dis();

      int face;

      const Front_VF& front_vf_ext=ref_cast(Front_VF, front_ext);
      const Front_VF& front_vf=ref_cast(Front_VF, front);

      OBS_PTR(Milieu_base) le_milieu;
      le_milieu = pbF.milieu();

      const int nb_comp = le_milieu->conductivite().nb_comp();
      int i=0;


      h_imp_.typer("Champ_front_fonc");
      h_imp_->fixer_nb_comp(nb_comp);
      DoubleTab& tab= h_imp_->valeurs();


      const Domaine_dis_base& domaine_dis1 = domaine_Cl_dis().domaine_dis();
      const Nom nom_racc1=frontiere_dis().frontiere().le_nom();



      const int nb_faces_front_vdf = front_vf.nb_faces();
      const int num_prem_face_vdf = front_vf.num_premiere_face();



      DoubleVect cumulSurfaces(nb_faces_front_vdf);
      cumulSurfaces = 0.;
      DoubleVect cumulHiSi(nb_faces_front_vdf);
      cumulHiSi = 0.;
      DoubleVect valHiSi(front_vf_ext.nb_faces());
      valHiSi = 0.;


      tab.resize(nb_faces_front_vdf,nb_comp);
      tab = 0.;


      if (domaine_dis1.domaine().raccord(nom_racc1)->que_suis_je() =="Raccord_distant_homogene")
        {
          //POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!
          Cerr<<"POUR LE MOMENT ON NE TRAITE PAS CE CAS !!!!!!!!!"<<finl;
          exit();
        }        // fin du cas Raccord_distant_homogene
      else // Raccord_local_homogene
        {
          const int ndeb = front_vf_ext.num_premiere_face();
          const int nfin = ndeb + front_vf_ext.nb_faces();

          // Calcul de tab = 1/(e/lambda + 1/h_paroi)
          if(!sub_type(Champ_Uniforme,le_milieu->conductivite()))
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

                      assert(lambda(elemF)!=0.);
                      cumulSurfaces(num_faceG-num_prem_face_vdf) += surfacesVEF(zvef,face,dimension);
                      valHiSi(face-ndeb) = surfacesVEF(zvef,face,dimension)/(1./pdt_scalSqrt(zvef,face,face,elemF,dimension,lambda(elemF))+1./h_paroi);
                      cumulHiSi(num_faceG-num_prem_face_vdf) += valHiSi(face-ndeb);
                      for(i=0; i<nb_comp; i++)
                        {
                          tab(num_faceG-num_prem_face_vdf,i) += valHiSi(face-ndeb) ;
                        }
                    }
                  for(face = 0; face<nb_faces_front_vdf; face++)
                    {
                      for(i=0; i<nb_comp; i++)
                        if (cumulSurfaces(face)>0.) tab(face,i) = tab(face,i)/ cumulSurfaces(face);
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

                    cumulSurfaces(num_faceG-num_prem_face_vdf) += surfacesVEF(zvef,face,dimension);
                    valHiSi(face-ndeb) = surfacesVEF(zvef,face,dimension)/(1./pdt_scalSqrt(zvef,face,face,elemF,dimension,lambda(elemF,0))+1./h_paroi);
                    cumulHiSi(num_faceG-num_prem_face_vdf) += valHiSi(face-ndeb);
                    for(i=0; i<nb_comp; i++)
                      {
                        assert(lambda(elemF,i)!=0.);
                        tab(num_faceG-num_prem_face_vdf,i) += surfacesVEF(zvef,face,dimension)/(1./pdt_scalSqrt(zvef,face,face,elemF,dimension,lambda(elemF,i))+1./h_paroi);
                      }
                  }
              for(face = 0; face<nb_faces_front_vdf; face++)
                {
                  for(i=0; i<nb_comp; i++)
                    if (cumulSurfaces(face)>0.) tab(face,i) = tab(face,i)/ cumulSurfaces(face);
                }
            }
          else  // la conductivite est un OWN_PTR(Champ_base) uniforme
            {
              assert(h_paroi!=0.);
              const DoubleTab& lambda = le_milieu->conductivite().valeurs();
              for(i=0; i<nb_comp; i++)
                assert(lambda(0,i)!=0.); // juste des asserts

              //Cerr << "raccord local homogene et conductivite uniforme " << finl;
              for (face=ndeb; face<nfin; face++)
                {
                  num_faceG = connect(face);
                  elemF =  face_voisinsF(face,0);
                  if(elemF == -1)
                    elemF =  face_voisinsF(face,1);
                  assert(num_faceG != -1);
                  cumulSurfaces(num_faceG-num_prem_face_vdf) += surfacesVEF(zvef,face,dimension);
                  valHiSi(face-ndeb) = surfacesVEF(zvef,face,dimension)/(1./pdt_scalSqrt(zvef,face,face,elemF,dimension,lambda(0,0))+1./h_paroi);
                  cumulHiSi(num_faceG-num_prem_face_vdf) += valHiSi(face-ndeb);

                  for(i=0; i<nb_comp; i++)
                    {
                      tab(num_faceG-num_prem_face_vdf,i) += valHiSi(face-ndeb) ;
                    }
                }
              for(face = 0; face<nb_faces_front_vdf; face++)
                {
                  for(i=0; i<nb_comp; i++)
                    if (cumulSurfaces(face)>0.) tab(face,i) = tab(face,i)/ cumulSurfaces(face);
                }

            }


          //Calcul de Text = (Somme_i h_i*S_i*T_exti)/(Somme_i h_i*S_i)
          assert(nb_comp == 1); // on sait pas faire pour l'instant sinon!!
          if (nb_comp != 1)
            {
              Cerr << "On ne sait pas traiter un pb couple 'zoom' lorsque nb_composante de la conductivite est plus grand que 1 " << finl;
              exit();
            }


          //ON RECUPERE LA TEMPERATURE FINE
          Text_valeurs = 0.;
          const int nb_eqF = pbF.nombre_d_equations();
          int trouve_eq = 0;
          int nb_equationF = 0;
          int face2, face3, face4;
          while((trouve_eq == 0) && (nb_equationF<nb_eqF))
            {
              const  Equation_base& eqF = pbF.equation(nb_equationF);
              if (sub_type(Conduction, eqF) || sub_type(Convection_Diffusion_std, eqF))
                {
                  trouve_eq = 1;
                  const Champ_Inc_base& incoF = eqF.inconnue();
                  const DoubleTab& valF = incoF.valeurs();
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


                      tempo += pdt_scal(zvef,face,face2,elemF,dimension,1.)*valF(face2) + pdt_scal(zvef,face,face3,elemF,dimension,1.)*valF(face3)  ;

                      if(dimension == 3)
                        {
                          face4 = zvef.elem_faces(elemF,2);
                          if (face4 == face)
                            face4 = zvef.elem_faces(elemF,dimension);
                          tempo += pdt_scal(zvef,face,face4,elemF,dimension,1.)*valF(face4);
                        }
                      tempo *=-valHiSi(face-ndeb)/pdt_scal(zvef,face,face,elemF,dimension,1.);
                      Text_valeurs(num_faceG-num_prem_face_vdf,0) += tempo;

                    }
                  for(face = 0; face<nb_faces_front_vdf; face++)
                    {
                      if (cumulHiSi(face)>0.) Text_valeurs(face,0) =  Text_valeurs(face,0) /(cumulHiSi(face));
                      // Cerr << "Text = " << Text_valeurs(face,0) << finl;
                    }

                }
              else
                nb_equationF++;
            }


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
  //  Cerr<<"fin de Echange_contact_VDF_VEF_Zoom::mettre_a_jour"<<finl;
}



