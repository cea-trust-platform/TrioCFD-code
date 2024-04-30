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
// File:        Prolongement_face_face_lin.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_face_face_lin.h>

Implemente_instanciable(Prolongement_face_face_lin,"Prolongement_face_face_lin",Prolongement_base);

//// printOn
//

Sortie& Prolongement_face_face_lin::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_face_face_lin::readOn(Entree& s )
{
  return s ;
}




/*! @brief Prolongement de l'inconnue grossiere sur la frontiere fine pour inconnue vitesse en VDF 2D uniquement
 *
 *     prolongement lineaire pour les composantes normales
 *     et moyenne pour les composantes tangentielles.
 *     remarque : pas de couplage!!
 *     car valeur
 *
 */
void Prolongement_face_face_lin::prolonger(Domaine_VF& domaine_VFG, Domaine_VF& domaine_VFF,
                                           const Frontiere& frontF,
                                           IntVect& connect,
                                           const DoubleTab& valG, DoubleTab& tab,
                                           int nb_compo)
{

  int nbFacesF;
  int num_faceG;
  tab = 0.;
  int comp;

  const int prem_face_bord_fin  =  frontF.num_premiere_face();

  int nb_faces_front_fine = frontF.nb_faces();

  for (nbFacesF=0; nbFacesF<nb_faces_front_fine; nbFacesF++)
    {
      int num_face=prem_face_bord_fin+nbFacesF;
      //on recupere le numero de la face G correspondante
      num_faceG = connect(num_face);

      /* --->  partie non modifiee : A voir... */

      if (nb_compo == 1)
        {
          tab(nbFacesF, 0) = valG(num_faceG);
        }
      else if(sub_type(Domaine_VEF, domaine_VFG))
        {
          for (comp=0; comp<nb_compo; comp++)
            {
              tab(nbFacesF, comp) = valG(num_faceG, comp);
            }
        }
      else
        /* <--- fin de partie non modifieee */

        if(sub_type(Domaine_VDF, domaine_VFG) && sub_type(Domaine_VDF, domaine_VFF))
          {

            Cerr << "debut modif prolongement" << finl;

            Domaine_VDF& domaine_vdfG = ref_cast(Domaine_VDF, domaine_VFG);
            const IntVect& oriG = domaine_vdfG.orientation();
            const IntTab& face_voisG = domaine_vdfG.face_voisins();

            const Domaine& domaineg = domaine_vdfG.domaine();
            int nb_faces_elemG = domaineg.type_elem()->nb_faces();
            IntTab& elem_facesG =  domaine_vdfG.elem_faces();

            const DoubleTab& xvG = domaine_vdfG.xv();

            Domaine_VDF& domaine_vdfF = ref_cast(Domaine_VDF, domaine_VFF);
            const DoubleTab& xvF = domaine_vdfF.xv();

            int num_face_gros;
            int nbFacesG;
            int num_elemG;
            int num_elemG1 = -1;
            int num_elemG2 = -1;
            int modulodim;
            int num_faceG_interne=-1;

            int ori_faceG = oriG(num_faceG);

            /* modif de cette valeur */
            //              tab(nbFacesF, ori_faceG) = valG(num_faceG);
            /* modif ici : 2D CDF uniquemant!!! */

            /* pour la vitesse normale aux faces : */

            num_elemG = face_voisG(num_faceG,0);
            if(num_elemG == -1)
              num_elemG = face_voisG(num_faceG,1);

            /* recherche du numero interne de la face grossiere dans l'element grossier */
            for (nbFacesG=0; nbFacesG<nb_faces_elemG; nbFacesG++)
              {
                num_face_gros = elem_facesG(num_elemG, nbFacesG);
                if ( num_face_gros == num_faceG )
                  {
                    num_faceG_interne = nbFacesG;
                    break;
                  }
              }

            Cerr << "num face interne = " << num_faceG_interne << finl;

            /* recherche des faces d'orientation normales pour recuperer les valeurs dans les mailles voisines */
            modulodim = (ori_faceG + 1)%(Objet_U::dimension); // pour le 2D uniquement

            num_face_gros = elem_facesG(num_elemG,modulodim);
            /* recuperation de la valeur grossiere dans les mailles voisines */
            /* premiere maille voisine trouvee ie premier passage */
            Cerr << "premier passage" << finl;
            num_elemG1 = face_voisG(num_face_gros,0);
            if (num_elemG1 == num_elemG)
              num_elemG1 = face_voisG(num_face_gros,1);
            if (num_elemG1 == -1)
              num_elemG1 = num_elemG;

            num_face_gros = elem_facesG(num_elemG,modulodim+Objet_U::dimension);
            /* c'est la 2ieme face */
            Cerr << "deuxieme  passage" << finl;
            num_elemG2 = face_voisG(num_face_gros,1);
            if (num_elemG2 == num_elemG)
              num_elemG2 = face_voisG(num_face_gros,0);
            if (num_elemG2 == -1)
              num_elemG2 = num_elemG;

            Cerr << "num fG "<< num_elemG1 << " " << num_elemG2 << finl;

            // ici on a les deux faces voisines => on calcul la valeur a mettre dans le prolongement


            double diffcoord1,diffcoord2,diffvalG;

            diffcoord1 = xvF(num_face,modulodim)-xvG(elem_facesG(num_elemG1,num_faceG_interne),modulodim);
            diffvalG = valG(elem_facesG(num_elemG2,num_faceG_interne))-valG(elem_facesG(num_elemG1,num_faceG_interne));
            diffcoord2 = xvG(elem_facesG(num_elemG2,num_faceG_interne),modulodim)- xvG(elem_facesG(num_elemG1,num_faceG_interne),modulodim);

            tab(nbFacesF, ori_faceG) = diffcoord1*diffvalG/diffcoord2 + valG(elem_facesG(num_elemG1,num_faceG_interne));

            Cerr << diffcoord2 << " tab(" << nbFacesF <<","<< ori_faceG << ") = " << tab(nbFacesF, ori_faceG) << finl;
            Cerr << "fin modif prolongement" << finl;

            /* pour la vitesse tangentielle aux faces : */
            //              int modulodim;
            for (comp=1 ; comp < nb_compo ; comp++)
              {

                modulodim = (ori_faceG + comp)%(Objet_U::dimension);
                num_elemG = face_voisG(num_faceG,0);
                int nb_faces_gros = 0;
                if(num_elemG != -1)
                  {
                    //pour chaque faceG
                    for(nbFacesG=0; nbFacesG<nb_faces_elemG; nbFacesG++)
                      {
                        num_face_gros = elem_facesG(num_elemG, nbFacesG);
                        if(oriG(num_face_gros) == modulodim)
                          {
                            tab(nbFacesF, modulodim) = tab(nbFacesF, modulodim)+ valG(num_face_gros);
                            nb_faces_gros = nb_faces_gros+1;
                          }
                      }
                  }
                num_elemG = face_voisG(num_faceG,1);
                if(num_elemG != -1)
                  {
                    //pour chaque faceG
                    for(nbFacesG=0; nbFacesG<nb_faces_elemG; nbFacesG++)
                      {
                        num_face_gros = elem_facesG(num_elemG, nbFacesG);
                        if(oriG(num_face_gros) == modulodim)
                          {
                            tab(nbFacesF, modulodim) = tab(nbFacesF, modulodim) + valG(num_face_gros);
                            nb_faces_gros = nb_faces_gros+1;
                          }
                      }
                  }
                tab(nbFacesF, modulodim) = tab(nbFacesF, modulodim) / nb_faces_gros;
              } // du for modulodim
          } // du if VDF
    } // du nb_face


  //Cerr<<"fin de Prolongement_face_face_lin::prolonger"<<finl;
}




//NE FAIT RIEN
void Prolongement_face_face_lin::calculer(Domaine_VF& domainef,
                                          Domaine_VF& domaineg,
                                          IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}
