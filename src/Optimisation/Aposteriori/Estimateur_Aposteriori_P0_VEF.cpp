/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Estimateur_Aposteriori_P0_VEF.h>
#include <Navier_Stokes_Aposteriori.h>
#include <Terme_Source_Qdm_VEF_Face.h>
#include <Champ_P1_isoP1Bulle.h>
#include <Champ_P1NC.h>
#include <Champ_Don.h>
#include <Source.h>

Implemente_instanciable( Estimateur_Aposteriori_P0_VEF, "Estimateur_Aposteriori_P0_VEF", Champ_Fonc_P0_VEF ) ;

Sortie& Estimateur_Aposteriori_P0_VEF::printOn( Sortie& s ) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Estimateur_Aposteriori_P0_VEF::readOn( Entree& s ) { return s; }

void Estimateur_Aposteriori_P0_VEF::associer_champ(const Champ_P1NC& la_vitesse, const Champ_P1_isoP1Bulle& la_pression, const Champ_Don& la_viscosite_cinematique, const Domaine_Cl_dis_base& le_dom_Cl_dis_base)
{
  le_dom_Cl_VEF  = ref_cast(Domaine_Cl_VEF, le_dom_Cl_dis_base);
  vitesse_= la_vitesse;
  pression_p1isop1b_ = la_pression;
  viscosite_cinematique_ = la_viscosite_cinematique;
}

void Estimateur_Aposteriori_P0_VEF::mettre_a_jour(double tps)
{
  const Domaine_Cl_VEF& domaine_cl_VEF = le_dom_Cl_VEF.valeur();
  const Domaine_VEF& domaine_VEF = domaine_cl_VEF.domaine_vef();
  const Domaine_VEF& domaine_VEF_p = ref_cast(Domaine_VEF, le_dom_VF.valeur());
  const Domaine& domaine = domaine_VEF_p.domaine();

  const DoubleTab& face_normales = domaine_VEF.face_normales();
  const DoubleVect& face_surfaces = domaine_VEF.face_surfaces();
  const int nb_faces = domaine_VEF.nb_faces_tot();
  DoubleVect coeff_faces(nb_faces);
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  int premiere_face_int = domaine_VEF.premiere_face_int();
  const IntTab& som_elem=domaine.les_elems();

  double sautdirectionnel;
  double estimtot;
  double contri_elem_i;
  double signe;
  double estimloc;

  int dim = Objet_U::dimension;
  int fac,elem1,elem2;
  int i,j,ielem;
  int som_opp;

  int nps=domaine_VEF_p.numero_premier_sommet();
  int nb_elem = domaine_VEF.nb_elem();
  int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const DoubleVect& inverse_volumes = domaine_VEF.inverse_volumes();
  DoubleTab la_vitesse = vitesse_->valeurs();
  DoubleTab la_pression = pression_p1isop1b_->valeurs();

  DoubleTab le_terme_source(vitesse_->valeurs());


  Navier_Stokes_Aposteriori& eq = ref_cast(Navier_Stokes_Aposteriori,vitesse_->equation());
  //const Sources& sources_eq =  eq.sources();

  OWN_PTR(Champ_base) espace_stockage;
  const Champ_base& ch_bs = eq.get_champ_source().get_champ(espace_stockage);
  const DoubleTab& src = ch_bs.valeurs();

  //Cerr << src << finl;

  le_terme_source = src;
  DoubleTab gradient_elem(nb_elem_tot,dim,dim);
  gradient_elem=0.;
  const DoubleTab tabnu=viscosite_cinematique_->valeur().valeurs();
  int sz=tabnu.size();

  DoubleVect visc(nb_elem);
  if ( sz == 1)
    {
      visc = tabnu(0,0);
    }
  else
    {
      visc = tabnu(0);
    }
  coeff_faces = 1.0; // Si dim = 2
  if (dim == 3)
    {
      for (fac=0; fac<nb_faces; fac++)
        {
          coeff_faces(fac) = pow(face_surfaces(fac),0.5);
        }
    }


  DoubleTab& estim = valeurs(); // Estimateur
  estim = 0.;

// Residus aux elements (le Laplacien de vitesse et le gradient de pression P0 sont nuls; il reste le terme source, le gradient de pression P1 et le terme de convection)

  // Calcul du terme de convection "u . grad u"
  // Calcul d'une vitesse moyenne par element, comme moyenne des (dim + 1) vitesses à ses faces; sert à calculer "u" dans "u . grad u"

  DoubleTab vite_moye(nb_elem,dim);
  vite_moye=0.;
  DoubleTab convec(nb_elem,dim);
  convec=0.;

  for (fac=0; fac<nb_faces; fac++) // calcul de vite_moye comme moyenne de la vitesse
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      if (elem1>=0)
        {
          for (i=0; i<dim; i++)
            {
              vite_moye(elem1,i) += la_vitesse(fac,i);
            }
        }
      if (elem2>=0)
        {
          for (i=0; i<dim; i++)
            {
              vite_moye(elem2,i) += la_vitesse(fac,i);
            }
        }
    }
  vite_moye /= (dim+1.0);

  // Calcul du gradient des vitesses; sert à calculer "grad u" dans "u . grad u"
  Champ_P1NC::calcul_gradient(la_vitesse, gradient_elem, le_dom_Cl_VEF.valeur());


  // Calcul de "u . grad u"

  for (ielem=0; ielem<nb_elem; ielem++)
    {
      for (i=0; i<dim; i++)
        {
          for (j=0; j<dim; j++)
            {
              convec(ielem,i) += vite_moye(ielem,j) * gradient_elem(ielem,i,j);
            }
        }
    }
  // Calcul du gradient de pression P1, si applicable
  DoubleTab grad_pression(nb_elem,dim);
  grad_pression = 0.;
  if (domaine_VEF_p.get_alphaS()) // On a une pression P1
    {
      for (fac=0; fac<nb_faces; fac++)
        {
          elem1=face_voisins(fac,0);
          elem2=face_voisins(fac,1);

          if (elem1>=0)
            {
              signe=1.;
              som_opp = le_dom_VF->get_num_fac_loc(fac, 0);
              som_opp = som_elem(elem1, som_opp);
              for (j=0; j<dim; j++)
                {
                  grad_pression(elem1,j) -= signe * la_pression(nps+som_opp) * face_normales(fac,j);
                }
            }

          if (elem2>=0)
            {
              signe=-1.;
              som_opp = le_dom_VF->get_num_fac_loc(fac, 1);
              som_opp = som_elem(elem2, som_opp);
              for (j=0; j<dim; j++)
                {
                  grad_pression(elem2,j) -= signe * la_pression(nps+som_opp) * face_normales(fac,j);
                }
            }
        }
      for (ielem=0; ielem<nb_elem; ielem++)
        {
          for (j=0; j<dim; j++)
            {
              grad_pression(ielem,j) *= inverse_volumes(ielem)/dim;
            }
        }
    }
  // Calcul du terme source par element, comme moyenne des sources aux faces
  DoubleTab source_moye(nb_elem,dim);
  source_moye=0.;
  for (fac=0; fac<nb_faces; fac++) // calcul de vite_moye comme moyenne de la vitesse
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      if (elem1>=0)
        {
          for (i=0; i<dim; i++)
            {
              source_moye(elem1,i) += le_terme_source(fac,i);
            }
        }
      if (elem2>=0)
        {
          for (i=0; i<dim; i++)
            {
              source_moye(elem2,i) += le_terme_source(fac,i);
            }
        }
    }
  source_moye /= (dim+1.0);

// Calcul des Residus aux elements
  for (ielem=0; ielem<nb_elem; ielem++)
    {
      estimloc = 0.;
      for (i=0; i<dim; i++)
        {
          contri_elem_i = ( source_moye(ielem,i) - grad_pression(ielem,i) - convec(ielem,i) );
          estimloc += contri_elem_i * contri_elem_i;
        }
      estim(ielem) += estimloc / pow( inverse_volumes(ielem), (1.0 + 2.0/dim) );
    }

// Fin Residus aux elements







//  Sauts au travers des faces - partie dans la direction de la normale (seules les faces internes contribuent)
  if (domaine_VEF_p.get_alphaE()) // Si on a une pression P0 elle est discontinue et donc contribue aux sauts au travers des faces
    {
      for (fac=premiere_face_int; fac<nb_faces; fac++)
        {
          elem1=face_voisins(fac,0);
          elem2=face_voisins(fac,1);
          if ( (elem1>=0) && (elem2>=0) )
            {
              for (i=0; i<dim; i++)
                {
                  sautdirectionnel = 0.;
                  for (j=0; j<dim; j++)
                    {
                      sautdirectionnel += ( visc(elem1)*gradient_elem(elem1,i,j) - visc(elem2)*gradient_elem(elem2,i,j) ) * face_normales(fac,j);
                    }
                  sautdirectionnel -= ( la_pression(elem1) - la_pression(elem2) ) * face_normales(fac,i);
                  sautdirectionnel *= sautdirectionnel;
                  sautdirectionnel /= coeff_faces(fac);
                  estim(elem1) += sautdirectionnel/2.0;
                  estim(elem2) += sautdirectionnel/2.0;
                }
            }
        }
    }
  else // Sinon la pression est continue et ne contribue donc pas aux sauts
    {
      for (fac=premiere_face_int; fac<nb_faces; fac++)
        {
          elem1=face_voisins(fac,0);
          elem2=face_voisins(fac,1);
          if ( (elem1>=0) && (elem2>=0) )
            {
              for (i=0; i<dim; i++)
                {
                  sautdirectionnel = 0.;
                  for (j=0; j<dim; j++)
                    {
                      sautdirectionnel += ( visc(elem1)*gradient_elem(elem1,i,j) - visc(elem2)*gradient_elem(elem2,i,j) ) * face_normales(fac,j);
                    }
                  sautdirectionnel *= sautdirectionnel;
                  sautdirectionnel /= coeff_faces(fac);
                  estim(elem1) += sautdirectionnel/2.0;
                  estim(elem2) += sautdirectionnel/2.0;
                }
            }
        }
    }

//  Sauts au travers des faces - partie dans la direction tangente à la face (1 contribution si dim=2 et 3 contributions si dim=3)

  if (dim == 3)
    {
      for (fac=premiere_face_int; fac<nb_faces; fac++) // faces internes
        {
          elem1=face_voisins(fac,0);
          elem2=face_voisins(fac,1);
          if ( (elem1>=0) && (elem2>=0) )
            {
              for (i=0; i<dim; i++)
                {
                  sautdirectionnel  = ( visc(elem1)*gradient_elem(elem1,i,1) - visc(elem2)*gradient_elem(elem2,i,1) ) * face_normales(fac,2) - ( visc(elem1)*gradient_elem(elem1,i,2) - visc(elem2)*gradient_elem(elem2,i,2) ) * face_normales(fac,1); // composante 1 du produit vectoriel : [g(u_i)]_2 * n_3 - [g(u_i)]_3 * n_2
                  sautdirectionnel *= sautdirectionnel;
                  sautdirectionnel /= coeff_faces(fac);
                  estim(elem1) += sautdirectionnel/2.0;
                  estim(elem2) += sautdirectionnel/2.0;

                  sautdirectionnel  = ( visc(elem1)*gradient_elem(elem1,i,2) - visc(elem2)*gradient_elem(elem2,i,2) ) * face_normales(fac,0) - ( visc(elem1)*gradient_elem(elem1,i,0) - visc(elem2)*gradient_elem(elem2,i,0) ) * face_normales(fac,2); // composante 2 du produit vectoriel : [g(u_i)]_3 * n_1 - [g(u_i)]_1 * n_3
                  sautdirectionnel *= sautdirectionnel;
                  sautdirectionnel /= coeff_faces(fac);
                  estim(elem1) += sautdirectionnel/2.0;
                  estim(elem2) += sautdirectionnel/2.0;

                  sautdirectionnel  = ( visc(elem1)*gradient_elem(elem1,i,0) - visc(elem2)*gradient_elem(elem2,i,0) ) * face_normales(fac,1) - ( visc(elem1)*gradient_elem(elem1,i,1) - visc(elem2)*gradient_elem(elem2,i,1) ) * face_normales(fac,0); // composante 3 du produit vectoriel : [g(u_i)]_1 * n_2 - [g(u_i)]_2 * n_1
                  sautdirectionnel *= sautdirectionnel;
                  sautdirectionnel /= coeff_faces(fac);
                  estim(elem1) += sautdirectionnel/2.0;
                  estim(elem2) += sautdirectionnel/2.0;
                }
            }
        }

      /*      for (fac=0; fac<premiere_face_int; fac++) // faces du bord (non implemente en attendant un traitement complet des differentes conditions aux limites possibles)
              {
                elem1=face_voisins(fac,0);
                elem2=face_voisins(fac,1);
                if ( elem1 > 0)
                  {
                    elem = elem1;
                  }
                if ( elem2 > 0)
                  {
                    elem = elem2;
                  }
                for (i=0; i<dim; i++)
                  {
                    sautdirectionnel  = visc(elem)*gradient_elem(elem,i,1) * face_normales(fac,2) - visc(elem)*gradient_elem(elem,i,2) * face_normales(fac,1); // composante 1 du produit vectoriel : [g(u_i)]_2 * n_3 - [g(u_i)]_3 * n_2
                    sautdirectionnel *= sautdirectionnel;
                    sautdirectionnel /= coeff_faces(fac);
                    estim(elem) += sautdirectionnel;

                    sautdirectionnel  = visc(elem)*gradient_elem(elem,i,2) * face_normales(fac,0) - visc(elem)*gradient_elem(elem,i,0) * face_normales(fac,2); // composante 2 du produit vectoriel : [g(u_i)]_3 * n_1 - [g(u_i)]_1 * n_3
                    sautdirectionnel *= sautdirectionnel;
                    sautdirectionnel /= coeff_faces(fac);
                    estim(elem) += sautdirectionnel;

                    sautdirectionnel  = visc(elem)*gradient_elem(elem,i,0) * face_normales(fac,1) - visc(elem)*gradient_elem(elem,i,1) * face_normales(fac,0); // composante 3 du produit vectoriel : [g(u_i)]_1 * n_2 - [g(u_i)]_2 * n_1
                    sautdirectionnel *= sautdirectionnel;
                    sautdirectionnel /= coeff_faces(fac);
                    estim(elem) += sautdirectionnel;
                  }
              } */
    }



  if (dim == 2)
    {
      for (fac=premiere_face_int; fac<nb_faces; fac++) // faces internes
        {
          elem1=face_voisins(fac,0);
          elem2=face_voisins(fac,1);
          if ( (elem1>=0) && (elem2>=0) )
            {
              for (i=0; i<dim; i++)
                {
                  sautdirectionnel = 0.;
                  sautdirectionnel  = ( visc(elem1)*gradient_elem(elem1,i,0) - visc(elem2)*gradient_elem(elem2,i,0) ) * face_normales(fac,1) - ( visc(elem1)*gradient_elem(elem1,i,1) - visc(elem2)*gradient_elem(elem2,i,1) ) * face_normales(fac,0);
                  sautdirectionnel *= sautdirectionnel;
                  sautdirectionnel /= coeff_faces(fac);
                  estim(elem1) += sautdirectionnel/2.0;
                  estim(elem2) += sautdirectionnel/2.0;
                }
            }
        }

      /*  for (fac=0; fac<premiere_face_int; fac++) // faces du bord (non implemente en attendant un traitement complet des differentes conditions aux limites possibles)
          {
            elem1=face_voisins(fac,0);
            elem2=face_voisins(fac,1);
            cout << "17 " << endl;
            if ( elem1 > 0)
              {
                ielem = elem1;
              }
            if ( elem2 > 0)
              {
                ielem = elem2;
              }
            for (i=0; i<dim; i++)
              {
                cout << "18 ielem" << ielem << endl;
                sautdirectionnel  = visc(ielem)*gradient_elem(ielem,i,0) * face_normales(fac,1) - visc(ielem)*gradient_elem(ielem,i,1) * face_normales(fac,0);
                sautdirectionnel *= sautdirectionnel;
                sautdirectionnel /= coeff_faces(fac);
                estim(ielem) += sautdirectionnel;
              }
          } */
    }


  estimtot = 0.;

  for (ielem=0; ielem<nb_elem; ielem++)
    {
      estimtot += estim(ielem);
      estim(ielem) = sqrt(estim(ielem));
    }
  estimtot = sqrt(estimtot);
  cout << "Estimateur_Aposteriori_P0_VEF::mettre_a_jour - estimateur = " << estimtot << endl;

  changer_temps(tps);
  Champ_Fonc_base::mettre_a_jour(tps);
}
