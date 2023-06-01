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
// File:        Modele_Jones_Launder_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Jones_Launder_VDF.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <TRUSTTrav.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>
#include <Dirichlet_paroi_defilante.h>
#include <Echange_externe_impose.h>
#include <Periodique.h>
#include <Symetrie.h>
#include <Neumann.h>
#include <Neumann_homogene.h>
#include <Champ_Face_VDF.h>
#include <Champ_Uniforme.h>
#include <Milieu_base.h>

Implemente_instanciable(Modele_Jones_Launder_VDF,"Modele_Jones_Launder_VDF",Modele_Fonc_Bas_Reynolds_Base);
// XD Jones_Launder modele_fonction_bas_reynolds_base Jones_Launder -1 Model described in ' Jones, W. P. and Launder, B. E. (1972), The prediction of laminarization with a two-equation model of turbulence, Int. J. of Heat and Mass transfer, Vol. 15, pp. 301-314.'
// printOn et readOn

Sortie& Modele_Jones_Launder_VDF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_Jones_Launder_VDF::readOn(Entree& is )
{
  Motcle motlu, accolade_fermee="}", accolade_ouverte="{";
  is >> motlu;
  if (motlu==accolade_ouverte)
    {
      is >> motlu;
      if (motlu != accolade_fermee)
        {
          Cerr << "Erreur a la lecture du Modele fonc bas reynolds Jones et Launder" << finl;
          Cerr << "On attendait } a la place de " << motlu << finl;
          exit();
        }
    }
  else
    {
      Cerr << "Erreur a la lecture du Modele fonc bas reynolds Jones et Launder" << finl;
      Cerr << "On attendait { a la place de " << motlu << finl;
      exit();
    }
  return is;
}

Entree& Modele_Jones_Launder_VDF::lire(const Motcle& , Entree& is)
{
  return is;
}
///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////

void  Modele_Jones_Launder_VDF::associer(const Domaine_dis& domaine_dis,
                                         const Domaine_Cl_dis& domaine_Cl_dis)
{
  //  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  //  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis.valeur());
}



DoubleTab& Modele_Jones_Launder_VDF::Calcul_D(DoubleTab& D,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                              const DoubleTab& vitesse,const DoubleTab& K_eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);

  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis.valeur());
  D = 0;
  //  return D;
  //  const DoubleVect& volumes = le_dom.volumes();
  const DoubleVect& porosite_surf = domaine_Cl_dis->equation().milieu().porosite_face();
  const DoubleVect& volume_entrelaces = le_dom.volumes_entrelaces();
  //  int nb_elem = le_dom.nb_elem();
  int nb_elem_tot = le_dom.nb_elem_tot();
  const Domaine& domaine=le_dom.domaine();

  int nb_faces_elem = domaine.nb_faces_elem();
  IntTrav numfa(nb_faces_elem);
  double coef;
  //  const IntTab& elem_faces = le_dom.elem_faces();
  const IntTab& face_voisins = le_dom.face_voisins();
  int nb_faces = le_dom.nb_faces();

  double gradk;
  int num_face,poly1,poly2,ori, ndeb, nfin;

  // Calcul de Gradient de racine de K.
  if (mp_min_vect(K_eps_Bas_Re)<0)
    {
      Cerr << "Il y'a des valeurs negatives dans les valeurs de K" << finl;
      Cerr << "dans Modele_Jones_Launder_VDF::Calcul_D" << finl;
      Cerr << "On arrete le calcul." << finl;
      exit();
    }
  // Boucle sur les bords pour traiter les conditions aux limites
  for (int n_bord=0; n_bord<le_dom.nb_front_Cl(); n_bord++)

    {
      const Cond_lim& la_cl = le_dom_Cl.les_conditions_limites(n_bord);

      if ( sub_type(Dirichlet,la_cl.valeur())||
           sub_type(Dirichlet_homogene,la_cl.valeur()) ||
           sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ||
           sub_type(Echange_externe_impose,la_cl.valeur())
         )


        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              gradk = 0;
              poly1 =  face_voisins(num_face,0);
              if (poly1 != -1)
                {
                  // coef = 0.5;
                  coef = volume_entrelaces(num_face)*porosite_surf(num_face)*0.5;
                  gradk = ( - sqrt(K_eps_Bas_Re(poly1,0)))/le_dom.dist_norm_bord(num_face);
                  if (!is_visco_const)
                    visco=tab_visco[poly1];
                  D[poly1] += 2*visco*(gradk*gradk)*coef;
                }
              else
                {
                  poly2 = face_voisins(num_face,1);
                  // coef = 0.5;
                  coef = volume_entrelaces(num_face)*porosite_surf(num_face)*0.5;
                  gradk = ((sqrt(K_eps_Bas_Re(poly2,0)) ))/le_dom.dist_norm_bord(num_face);
                  //
                  if (!is_visco_const)
                    visco=tab_visco[poly2];
                  D[poly2] += 2*visco*(gradk*gradk)*coef;
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              gradk = 0;
              poly1 = face_voisins(num_face,0);
              poly2 = face_voisins(num_face,1);
              ori = le_dom.orientation(num_face);
              // coef = 0.5;

              coef = volume_entrelaces(num_face)*porosite_surf(num_face);


              gradk =  (sqrt(K_eps_Bas_Re(poly2,0))-sqrt(K_eps_Bas_Re(poly1,0)))/le_dom.dist_elem_period(poly1,poly2,ori);
              if (!is_visco_const)
                visco=tab_visco[poly1];
              D[poly1] += 2*visco*(gradk*gradk)*coef;
              if (!is_visco_const)
                visco=tab_visco[poly2];
              D[poly2] += 2*visco*(gradk*gradk)*coef;
            }

        }
      else if (sub_type(Symetrie,la_cl.valeur()))
        ;
      else if ( (sub_type(Neumann,la_cl.valeur()))
                ||
                (sub_type(Neumann_homogene,la_cl.valeur()))
              )
        {
          // do nothing
          ;
        }
      else
        {
          Cerr<<la_cl.valeur().que_suis_je()<< "not implemented in calculer_D"<<finl;
          exit();
        }
    }

  // Traitement des faces internes
  for (num_face=le_dom.premiere_face_int(); num_face<nb_faces; num_face++)
    {
      poly1 = face_voisins(num_face,0);
      poly2 = face_voisins(num_face,1);
      ori = le_dom.orientation(num_face);
      // coef = 0.5;

      coef = volume_entrelaces(num_face)*porosite_surf(num_face);

      gradk =  (sqrt(K_eps_Bas_Re(poly2,0))-sqrt(K_eps_Bas_Re(poly1,0)))/(le_dom.xp(poly2,ori)- le_dom.xp(poly1,ori));
      //      Cerr<<" ici "<< num_face<< " "<<gradk*gradk/K_eps_Bas_Re(poly2,0)<<" K "<<K_eps_Bas_Re(poly2,0)/K_eps_Bas_Re(0,0)<<finl;
      if (num_face==-396)
        Cerr << "K_eps_Bas_Re(poly2,0)=" << K_eps_Bas_Re(poly2,0)/K_eps_Bas_Re(0,0) << " K_eps_Bas_Re(poly1,0)=" << K_eps_Bas_Re(poly1,0)/K_eps_Bas_Re(0,0) << " test "<<sqrt(  K_eps_Bas_Re(poly2,0))/sqrt( K_eps_Bas_Re(0,0))<<" "<<sqrt(  K_eps_Bas_Re(poly1,0))/sqrt( K_eps_Bas_Re(0,0))<< " "<<(sqrt(K_eps_Bas_Re(poly2,0))-sqrt(K_eps_Bas_Re(poly1,0)))/sqrt( K_eps_Bas_Re(0,0))<<" "<< K_eps_Bas_Re(0,0)<<finl;
      if (!is_visco_const)
        visco=tab_visco[poly1];
      D[poly1] += 2*visco*(gradk*gradk)*coef;
      if (!is_visco_const)
        visco=tab_visco[poly2];
      D[poly2] += 2*visco*(gradk*gradk)*coef;
    }
  // GF on a nb_face_elem contributions par elem ....
  // mais 2 contributions dans chaque direction

  // provisoire
  D/=2;
  const DoubleVect& volumes= le_dom.volumes();
  for (int i=0; i<nb_elem_tot; i++)
    D(i)/=volumes(i);
  //Cerr<<D.mp_min_vect()<<" DDDDDDDDDDDDDD "<<D.mp_max_vect()<<finl;
  //D=0;
  return D;
  // D abord sans le 1/3 2/3
}

// void Modele_Jones_Launder_VDF::associer_domaines(const Domaine_dis& domaine_dis,
//                                                         const Domaine_Cl_dis& domaine_Cl_dis)
// {
//   le_dom_VDF = ref_cast(Domaine_VDF, domaine_dis.valeur());
//   le_dom_Cl_VDF = ref_cast(Domaine_Cl_VDF, domaine_Cl_dis.valeur());
// }


DoubleTab& Modele_Jones_Launder_VDF::Calcul_E(DoubleTab& E,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& vit,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis.valeur());
  E = 0;
  //return E;
  // provisoire
  if (1)
    {
      const IntTab& face_voisins = le_dom.face_voisins();
      const IntTab& elem_faces = le_dom.elem_faces();

      const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eq_hydraulique->inconnue().valeur());
      int nb_elem_tot=le_dom.nb_elem_tot();
      assert (vitesse.valeurs().line_size() == 1);
      DoubleTab gij(nb_elem_tot,dimension,dimension, vitesse.valeurs().line_size());
      // Rque methode non const Pourquoi ?

      ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj(vitesse.valeurs(),gij,ref_cast(Domaine_Cl_VDF,eq_hydraulique->domaine_Cl_dis().valeur()));
      DoubleTab der_seconde(dimension,dimension,dimension);
      // der_seconde (i,j,k) d2 ui dxj dxk
      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          // calcul de la derivee seconde
          ArrOfInt num(dimension*2);
          for (int dir=0; dir<dimension; dir++)
            for (int sens=0; sens<2; sens++)
              {
                int f=elem_faces(elem,dir+dimension*sens);
                int num0=face_voisins(f,sens);
                if (num0 == -1)
                  num0=elem;
                num[dir+dimension*sens]=num0;
              }
          for (int i=0; i<dimension; i++)
            for (int j=0; j<dimension; j++)
              for (int k=0; k<dimension; k++)
                {
                  der_seconde(i,j,k)=0.5*(gij(num[k+dimension],i,j,0)-gij(num[k],i,j,0))/le_dom.dim_elem(elem,k);
                }

          if (!is_visco_const)
            visco=tab_visco[elem];
          if (dimension==2)
            {
              double val2=der_seconde(0,1,1)*der_seconde(0,1,1); //d2u/dy2 ^2
              double val3=der_seconde(1,0,0)*der_seconde(1,0,0); //d2v/dx2  ^2

              E[elem]= 2*visco*visco_turb(elem)*(val2+val3);
            }
          else
            {
              double val2 = der_seconde(1,0,0)+der_seconde(2,0,0); //d2v_dx2(fac)+d2w_dx2(fac);
              val2 *= val2;
              double val3 = der_seconde(0,1,1)+der_seconde(2,1,1); //d2u_dy2(fac)+d2w_dy2(fac);
              val3 *= val3;
              double val4 = der_seconde(0,2,2)+der_seconde(1,2,2); //d2u_dz2(fac)+d2v_dz2(fac);
              val4 *= val4;
              E[elem] += 2*visco*visco_turb(elem)*(val2+val3+val4);
            }
        }
      // E/=2.;
      //Cerr<<E.mp_min_vect()<<" E "<<E.mp_max_vect()<<finl;
      return E;

    }


  int nb_faces_tot = le_dom.nb_faces_tot();
  const IntTab& face_voisins = le_dom.face_voisins();
  int n_bord,poly1, poly2,fac=-1;

  DoubleTab derivee_premiere(0, Objet_U::dimension, Objet_U::dimension);
  le_dom.creer_tableau_faces(derivee_premiere);
  DoubleTab derivee_seconde(derivee_premiere);

  calcul_derivees_premieres_croisees(derivee_premiere,domaine_dis,domaine_Cl_dis,vit );
  calcul_derivees_secondes_croisees(derivee_seconde,domaine_dis,domaine_Cl_dis,derivee_premiere);

  // traitement par faces pour avoir la contribution par element de D
  double val2,val3,val4;
  for (n_bord=0; n_bord<le_dom.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();
      for (fac=ndeb; fac<nfin; fac++)
        {
          poly1 = face_voisins(fac,0);
          if (poly1 != -1)
            {
              if (Objet_U::dimension == 2 )
                {
                  val2 = derivee_seconde(fac,0,1)*derivee_seconde(fac,0,1);
                  val3 = derivee_seconde(fac,1,0)*derivee_seconde(fac,1,0);
                  if (!is_visco_const)
                    visco=tab_visco[poly1];
                  E[poly1]+= 2*visco*visco_turb(poly1)*(val2+val3);
                }
              else
                {
                  val2 = derivee_seconde(fac,1,0)+derivee_seconde(fac,2,0); //d2v_dx2(fac)+d2w_dx2(fac);
                  val2 *= val2;
                  val3 = derivee_seconde(fac,0,1)+derivee_seconde(fac,2,1); //d2u_dy2(fac)+d2w_dy2(fac);
                  val3 *= val3;
                  val4 = derivee_seconde(fac,0,2)+derivee_seconde(fac,1,2); //d2u_dz2(fac)+d2v_dz2(fac);
                  val4 *= val4;
                  if (!is_visco_const)
                    visco=tab_visco[poly1];
                  E[poly1] += 2*visco*visco_turb(poly1)*(val2+val3+val4);
                }
            }
          else
            {
              poly2 = face_voisins(fac,1);
              if (Objet_U::dimension == 2 )
                {
                  val2 = derivee_seconde(fac,0,1)*derivee_seconde(fac,0,1);
                  val3 = derivee_seconde(fac,1,0)*derivee_seconde(fac,1,0);
                  if (!is_visco_const)
                    visco=tab_visco[poly2];
                  E[poly2]+= 2*visco*visco_turb(poly2)*(val2+val3);
                }
              else
                {
                  val2 = derivee_seconde(fac,1,0)+derivee_seconde(fac,2,0); //d2v_dx2(fac)+d2w_dx2(fac);
                  val2 *= val2;
                  val3 = derivee_seconde(fac,0,1)+derivee_seconde(fac,2,1); //d2u_dy2(fac)+d2w_dy2(fac);
                  val3 *= val3;
                  val4 = derivee_seconde(fac,0,2)+derivee_seconde(fac,1,2); //d2u_dz2(fac)+d2v_dz2(fac);
                  val4 *= val4;
                  if (!is_visco_const)
                    visco=tab_visco[poly2];
                  E[poly2] += 2*visco*visco_turb(poly2)*(val2+val3+val4);
                }
            }
        }
    }

  // Traitement des faces internes
  for (; fac<nb_faces_tot; fac++)
    {
      poly1=face_voisins(fac,0);
      poly2=face_voisins(fac,1);
      if (Objet_U::dimension == 2 )
        {
          val2 = derivee_seconde(fac,0,1)*derivee_seconde(fac,0,1);
          val3 = derivee_seconde(fac,1,0)*derivee_seconde(fac,1,0);

          if (poly1>=0)
            {
              if (!is_visco_const)
                visco=tab_visco[poly1];
              E[poly1]+= 2*visco*visco_turb(poly1)*(val2+val3);
            }
          if (poly2>=0)
            {
              if (!is_visco_const)
                visco=tab_visco[poly2];
              E[poly2]+= 2*visco*visco_turb(poly2)*(val2+val3);
            }
        }
      else
        {
          val2 = derivee_seconde(fac,1,0)+derivee_seconde(fac,2,0); //d2v_dx2(fac)+d2w_dx2(fac);
          val2 *= val2;
          val3 = derivee_seconde(fac,0,1)+derivee_seconde(fac,2,1); //d2u_dy2(fac)+d2w_dy2(fac);
          val3 *= val3;
          val4 = derivee_seconde(fac,0,2)+derivee_seconde(fac,1,2); //d2u_dz2(fac)+d2v_dz2(fac);
          val4 *= val4;

          if (poly1>=0)
            {
              if (!is_visco_const)
                visco=tab_visco[poly1];
              E[poly1] += 2*visco*visco_turb(poly1)*(val2+val3+val4);
            }
          if (poly2>=0)
            {
              if (!is_visco_const)
                visco=tab_visco[poly1];
              E[poly2] += 2*visco*visco_turb(poly2)*(val2+val3+val4);
            }
        }
    }
  return E;
}

DoubleTab& Modele_Jones_Launder_VDF::calcul_derivees_premieres_croisees(DoubleTab& derivee_premiere, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& vit ) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis.valeur());
  //  Cerr<<le_dom_Cl.equation().le_nom()<<finl;exit();
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eq_hydraulique->inconnue().valeur());
  const Domaine_Cl_VDF& zcl_hydro = ref_cast(Domaine_Cl_VDF,eq_hydraulique->domaine_Cl_dis().valeur());
  //  int nb_faces = le_dom.nb_faces();
  const IntTab& Qdm = le_dom.Qdm();
  const IntVect& orientation = le_dom.orientation();
  IntTrav num(4);
  int i,num_arete;
  double pond1,pond2;

  int signe,ori0,ori1;
  // DEBUT CALCUL DES DERIVEES PREMIERES

  if (Objet_U::dimension == 2)
    {
      DoubleVect coef(2);
      int ndeb,nfin;
      // Boucle sur les aretes internes  pour le calcul
      // des moyennes des derivees croisees

      ndeb = le_dom.premiere_arete_interne();
      nfin = le_dom.nb_aretes_internes() + ndeb;

      for (num_arete =ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_premiere(num[0],ori0,ori1) +=  coef[0];
          derivee_premiere(num[1],ori0,ori1) +=  coef[0];
          derivee_premiere(num[2],ori1,ori0) +=  coef[1];
          derivee_premiere(num[3],ori1,ori0) +=  coef[1];
        }

      // Boucle sur les aretes_mixte
      ndeb = le_dom.premiere_arete_mixte();
      nfin = ndeb + le_dom.nb_aretes_mixtes();

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_premiere(num[0],ori0,ori1) +=  coef[0];
          derivee_premiere(num[1],ori0,ori1) +=  coef[0];
          derivee_premiere(num[2],ori1,ori0) +=  coef[1];
          derivee_premiere(num[3],ori1,ori0) +=  coef[1];
        }

      //*******************************
      //Prise en compte des CL
      //*******************************

      //*******************************
      //On parcourt les aretes bords
      //*******************************

      ndeb = le_dom.premiere_arete_bord();
      nfin = ndeb + le_dom.nb_aretes_bord();
      int n_type;

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          n_type=le_dom_Cl.type_arete_bord(num_arete-ndeb);

          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          if (n_type == TypeAreteBordVDF::PERIO_PERIO)
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              coef[0] =  0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] =  0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);
              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
              derivee_premiere(num[3],ori1,ori0) +=  coef[1];
            }
          else
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              signe = num[3];

              pond1 = le_dom.face_normales(num[0],ori0);
              pond2 = le_dom.face_normales(num[1],ori0);
              double tps=vitesse.temps();
              double vit_imp = pond2*Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[0],ori1,zcl_hydro)+
                               Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[1],ori1,zcl_hydro)*pond1; // val tangentielle
              vit_imp /= pond1+pond2;
              //               double vit_imp = 0.5*(vitesse.val_imp_face_bord(num[0],ori1)+
              //                                     vitesse.val_imp_face_bord(num[1],ori1));                // val tangentielle
              coef[0] =  0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              ////coef[1] =  0.5*(vit_imp-vit(num[2]))/le_dom.dist_norm_bord(num[2])*signe;
              coef[1] =  0.5*(vit_imp-vit(num[2]))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
            }
        }


      //*******************************
      //On parcourt les aretes coins
      //*******************************

      ndeb = le_dom.premiere_arete_coin();
      nfin = ndeb + le_dom.nb_aretes_coin();

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          n_type=le_dom_Cl.type_arete_coin(num_arete-ndeb);
          //***************************************
          // Traitement des aretes coin perio-perio
          //***************************************

          if (n_type == TypeAreteCoinVDF::PERIO_PERIO) // arete de type periodicite-periodicite
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              coef[0] =  0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] =  0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);
              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
              derivee_premiere(num[3],ori1,ori0) +=  coef[1];
            }

          //***************************************
          // Traitement des aretes coin perio-paroi
          //***************************************
          else if (n_type == TypeAreteCoinVDF::PERIO_PAROI) // arete de type periodicite-paroi
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              signe = num[3];

              pond1 = le_dom.face_normales(num[0],ori0);
              pond2 = le_dom.face_normales(num[1],ori0);
              double tps=vitesse.temps();
              double vit_imp = pond2*Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[0],ori1,zcl_hydro)+
                               Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[1],ori1,zcl_hydro)*pond1;
              // val tangentielle
              vit_imp /= pond1+pond2;
              //               double vit_imp = 0.5*(vitesse.val_imp_face_bord(num[0],ori1)+
              //                                     vitesse.val_imp_face_bord(num[1],ori1));                // val tangentielle
              coef[0] =  0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] =  0.5*(vit_imp-vit(num[2]))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
            }
          else
            {
              //Cerr << "Oh!!! Un coin!!!" << finl;
            }

        }

      //     Cerr << "Il y a : " << compt_coin << " coins!!" << finl;
    }

  else if (Objet_U::dimension == 3)
    {
      DoubleVect coef(3);

      // Boucle sur les aretes internes  pour le calcul
      // des moyennes des derivees croisees
      int ndeb = le_dom.premiere_arete_interne();
      int nfin = le_dom.nb_aretes_internes() + ndeb;

      for (num_arete =ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_premiere(num[0],ori0,ori1) +=  coef[0];
          derivee_premiere(num[1],ori0,ori1) +=  coef[0];
          derivee_premiere(num[2],ori1,ori0) +=  coef[1];
          derivee_premiere(num[3],ori1,ori0) +=  coef[1];
        }
      // Boucle sur les aretes_mixte
      ndeb = le_dom.premiere_arete_mixte();
      nfin = ndeb + le_dom.nb_aretes_mixtes();

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_premiere(num[0],ori0,ori1) +=  coef[0];
          derivee_premiere(num[1],ori0,ori1) +=  coef[0];
          derivee_premiere(num[2],ori1,ori0) +=  coef[1];
          derivee_premiere(num[3],ori1,ori0) +=  coef[1];
        }

      //*******************************
      //Prise en compte des CL
      //*******************************

      //*******************************
      //On parcourt les aretes bords
      //*******************************

      ndeb = le_dom.premiere_arete_bord();
      nfin = ndeb + le_dom.nb_aretes_bord();
      int n_type;

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          n_type=le_dom_Cl.type_arete_bord(num_arete-ndeb);
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          if (n_type == TypeAreteBordVDF::PERIO_PERIO)
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              if (num[0]==num[1] || num[2]==num[3])
                {
                  Cerr << "2 bords avec une condition limite de periodicite ne sont separees que d'une maille !" << finl;
                  Cerr << "Cela n'est pas valide pour le Modele Jones Launder et son calcul de derivees croisees..." << finl;
                  exit();
                }
              coef[0] = 0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] = 0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);

              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
              derivee_premiere(num[3],ori1,ori0) +=  coef[1];
            }
          else
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              signe = num[3];

              //               double vit_imp = 0.5*(vitesse.val_imp_face_bord(num[0],ori1)+
              //                                     vitesse.val_imp_face_bord(num[1],ori1));                // val tangentielle
              pond1 = le_dom.face_normales(num[0],ori0);
              pond2 = le_dom.face_normales(num[1],ori0);
              double tps=vitesse.temps();
              double vit_imp = pond2*Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[0],ori1,zcl_hydro)+
                               Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[1],ori1,zcl_hydro)*pond1; // val tangentielle
              // val tangentielle
              vit_imp /= pond1+pond2;

              coef[0] =  0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] =  0.5*(vit_imp-vit(num[2]))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
            }
        }

      //*******************************
      //On parcourt les aretes coins
      //*******************************

      ndeb = le_dom.premiere_arete_coin();
      nfin = ndeb + le_dom.nb_aretes_coin();

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          n_type=le_dom_Cl.type_arete_coin(num_arete-ndeb);

          //***************************************
          // Traitement des aretes coin perio-perio
          //***************************************

          if (n_type == TypeAreteCoinVDF::PERIO_PERIO) // arete de type periodicite-periodicite
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              coef[0] = 0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] = 0.5*(vit(num[3])-vit(num[2]))/le_dom.dist_face(num[2],num[3],ori0);

              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
              derivee_premiere(num[3],ori1,ori0) +=  coef[1];
            }

          //***************************************
          // Traitement des aretes coin perio-paroi
          //***************************************
          else if (n_type == TypeAreteCoinVDF::PERIO_PAROI) // arete de type periodicite-paroi
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              signe = num[3];

              pond1 = le_dom.face_normales(num[0],ori0);
              pond2 = le_dom.face_normales(num[1],ori0);
              double tps=vitesse.temps();
              double vit_imp = pond2*Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[0],ori1,zcl_hydro)+
                               Champ_Face_get_val_imp_face_bord_sym(vitesse.valeurs(),tps,num[1],ori1,zcl_hydro)*pond1;

              // val tangentielle
              vit_imp /= pond1+pond2;
              //               double vit_imp = 0.5*(vitesse.val_imp_face_bord(num[0],ori1)+
              //                                     vitesse.val_imp_face_bord(num[1],ori1));                // val tangentielle

              coef[0] =  0.5*(vit(num[1])-vit(num[0]))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] =  0.5*(vit_imp-vit(num[2]))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_premiere(num[0],ori0,ori1) +=  coef[0];
              derivee_premiere(num[1],ori0,ori1) +=  coef[0];
              derivee_premiere(num[2],ori1,ori0) +=  coef[1];
            }
          else
            {
              //Cerr << "Oh!!! Un coin!!!" << finl;
            }
        }
      // Cerr << "Il y a : " << compt_coin << " coins!!" << finl;
    }
  derivee_premiere.echange_espace_virtuel();
  return derivee_premiere;
}

DoubleTab& Modele_Jones_Launder_VDF::calcul_derivees_secondes_croisees(DoubleTab& derivee_seconde, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& derivee_premiere ) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis.valeur());
  //  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eq_hydraulique->inconnue().valeur());

  //  int nb_faces = le_dom.nb_faces();
  const IntTab& Qdm = le_dom.Qdm();
  const IntVect& orientation = le_dom.orientation();
  //  const int nb_cond_lim = le_dom_Cl.nb_cond_lim();
  IntTrav num(4);
  int i,num_arete,ori0,ori1;

  int ndeb=-10000,nfin=100000;//,signe;
  // DEBUT CALCUL DES DERIVEES SECONDES

  if (Objet_U::dimension == 2)
    {
      DoubleVect coef(2);

      ndeb = le_dom.premiere_arete_interne();
      nfin = le_dom.nb_aretes_internes() + ndeb;
      // Boucle sur les aretes internes  pour le calcul
      // des moyennes des derivees croisees

      for (num_arete =ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_seconde(num[0],ori0,ori1) +=  coef[0];
          derivee_seconde(num[1],ori0,ori1) +=  coef[0];
          derivee_seconde(num[2],ori1,ori0) +=  coef[1];
          derivee_seconde(num[3],ori1,ori0) +=  coef[1];
        }

      // Boucle sur les aretes_mixte
      ndeb = le_dom.premiere_arete_mixte();
      nfin = ndeb + le_dom.nb_aretes_mixtes();

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_seconde(num[0],ori0,ori1) +=  coef[0];
          derivee_seconde(num[1],ori0,ori1) +=  coef[0];
          derivee_seconde(num[2],ori1,ori0) +=  coef[1];
          derivee_seconde(num[3],ori1,ori0) +=  coef[1];
        }

      //*******************************
      //Prise en compte des CL
      //*******************************

      //*******************************
      //On parcourt les aretes bords
      //*******************************
      ndeb = le_dom.premiere_arete_bord();
      nfin = ndeb + le_dom.nb_aretes_bord();
      int n_type;

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          n_type=le_dom_Cl.type_arete_bord(num_arete-ndeb);

          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          if (n_type == TypeAreteBordVDF::PERIO_PERIO)
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              coef[0] =  0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] =  0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);
              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              derivee_seconde(num[2],ori1,ori0) +=  coef[1];
              derivee_seconde(num[3],ori1,ori0) +=  coef[1];
            }
          else
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              // signe = num[3];

              // ///???????????????
              //               double  val_deriv_prem = 0.5*(vit.val_imp_face_bord(num[0],ori1)+
              //                                     vit.val_imp_face_bord(num[1],ori1));                // val tangentielle
              // ///???????????????
              coef[0] =  0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              //coef[1] =  0.5*(val_deriv_prem-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              derivee_seconde(num[2],ori1,ori0) =  2.*derivee_seconde(num[2],ori1,ori0);
            }
        }

      //*******************************
      //On parcourt les aretes coins
      //*******************************

      ndeb = le_dom.premiere_arete_coin();
      nfin = ndeb + le_dom.nb_aretes_coin();

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          n_type=le_dom_Cl.type_arete_coin(num_arete-ndeb);
          //***************************************
          // Traitement des aretes coin perio-perio
          //***************************************

          if (n_type == TypeAreteCoinVDF::PERIO_PERIO) // arete de type periodicite-periodicite
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              coef[0] =  0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] =  0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);
              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              derivee_seconde(num[2],ori1,ori0) +=  coef[1];
              derivee_seconde(num[3],ori1,ori0) +=  coef[1];
            }

          //***************************************
          // Traitement des aretes coin perio-paroi
          //***************************************
          else if (n_type == TypeAreteCoinVDF::PERIO_PAROI) // arete de type periodicite-paroi
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              //signe = num[3];

              // ///???????????????
              //               double  val_deriv_prem = 0.5*(vit.val_imp_face_bord(num[0],ori1)+
              //                                     vit.val_imp_face_bord(num[1],ori1));                // val tangentielle
              // ///???????????????
              coef[0] =  0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              //coef[1] =  0.5*(val_deriv_prem-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              //              derivee_seconde(num[2],ori1,ori0) +=  coef[1];
              derivee_seconde(num[2],ori1,ori0) =  2.*derivee_seconde(num[2],ori1,ori0);
            }
          else
            {
              //Cerr << "Oh!!! Un coin!!!" << finl;
            }

        }
      //  Cerr << "Il y a : " << compt_coin << " coins!!" << finl;
    }

  else if (Objet_U::dimension == 3)
    {
      //   Cerr << " tableau des derivee_premiere " << derivee_premiere << finl;
      DoubleVect coef(3);

      // Boucle sur les aretes internes  pour le calcul
      // des moyennes des derivees croisees
      ndeb = le_dom.premiere_arete_interne();
      nfin = le_dom.nb_aretes_internes() + ndeb;

      for (num_arete =ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_seconde(num[0],ori0,ori1) +=  coef[0];
          derivee_seconde(num[1],ori0,ori1) +=  coef[0];
          derivee_seconde(num[2],ori1,ori0) +=  coef[1];
          derivee_seconde(num[3],ori1,ori0) +=  coef[1];
        }

      // Boucle sur les aretes_mixte
      ndeb = le_dom.premiere_arete_mixte();
      nfin = ndeb + le_dom.nb_aretes_mixtes();

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          ori0 = orientation(num[0]);
          ori1 = orientation(num[2]);
          coef[0] = 0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
          coef[1] = 0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);

          derivee_seconde(num[0],ori0,ori1) +=  coef[0];
          derivee_seconde(num[1],ori0,ori1) +=  coef[0];
          derivee_seconde(num[2],ori1,ori0) +=  coef[1];
          derivee_seconde(num[3],ori1,ori0) +=  coef[1];
        }

      //*******************************
      //Prise en compte des CL
      //*******************************

      //*******************************
      //On parcourt les aretes bords
      //*******************************

      ndeb = le_dom.premiere_arete_bord();
      nfin = ndeb + le_dom.nb_aretes_bord();
      int n_type;

      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          n_type=le_dom_Cl.type_arete_bord(num_arete-ndeb);
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          if (n_type == TypeAreteBordVDF::PERIO_PERIO)
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              coef[0] = 0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] = 0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);

              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              derivee_seconde(num[2],ori1,ori0) +=  coef[1];
              derivee_seconde(num[3],ori1,ori0) +=  coef[1];
            }
          else
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              //signe = num[3];

              // ///???????????????
              //               double  val_deriv_prem = 0.5*(vit.val_imp_face_bord(num[0],ori1)+
              //                                     vit.val_imp_face_bord(num[1],ori1));                // val tangentielle
              // ///???????????????

              coef[0] =  0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              //coef[1] =  0.5*(val_deriv_prem-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              //              derivee_seconde(num[2],ori1,ori0) +=  coef[1];
              derivee_seconde(num[2],ori1,ori0) =  2.*derivee_seconde(num[2],ori1,ori0);
            }
        }

      //*******************************
      //On parcourt les aretes coins
      //*******************************

      ndeb = le_dom.premiere_arete_coin();
      nfin = ndeb + le_dom.nb_aretes_coin();


      for (num_arete=ndeb; num_arete<nfin; num_arete++)
        {
          for (i=0; i<4; i++)
            num[i] = Qdm(num_arete,i);

          n_type=le_dom_Cl.type_arete_coin(num_arete-ndeb);
          //***************************************
          // Traitement des aretes coin perio-perio
          //***************************************
          if (n_type == TypeAreteCoinVDF::PERIO_PERIO) // arete de type periodicite-periodicite
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              coef[0] = 0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              coef[1] = 0.5*(derivee_premiere(num[3],ori1,ori0)-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_face(num[2],num[3],ori0);

              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              derivee_seconde(num[2],ori1,ori0) +=  coef[1];
              derivee_seconde(num[3],ori1,ori0) +=  coef[1];
            }

          //***************************************
          // Traitement des aretes coin perio-paroi
          //***************************************
          else if (n_type == TypeAreteCoinVDF::PERIO_PAROI) // arete de type periodicite-paroi
            {
              ori0 = orientation(num[0]);
              ori1 = orientation(num[2]);
              //signe = num[3];

              // ///???????????????
              //               double  val_deriv_prem = 0.5*(vit.val_imp_face_bord(num[0],ori1)+
              //                                     vit.val_imp_face_bord(num[1],ori1));                // val tangentielle
              // ///???????????????

              coef[0] =  0.5*(derivee_premiere(num[1],ori0,ori1)-derivee_premiere(num[0],ori0,ori1))/le_dom.dist_face(num[0],num[1],ori1);
              //coef[1] =  0.5*(val_deriv_prem-derivee_premiere(num[2],ori1,ori0))/le_dom.dist_norm_bord(num[0])*signe;
              derivee_seconde(num[0],ori0,ori1) +=  coef[0];
              derivee_seconde(num[1],ori0,ori1) +=  coef[0];
              derivee_seconde(num[2],ori1,ori0) =  2.*derivee_seconde(num[2],ori1,ori0);
              //              derivee_seconde(num[2],ori1,ori0) +=  coef[1];
            }
          else
            {
              //Cerr << "Oh!!! Un coin!!!" << finl;
            }
        }
      //  Cerr << "Il y a : " << compt_coin << " coins!!" << finl;
    }
  // Cerr << "derivee_seconde " << derivee_seconde << finl;
  derivee_seconde.echange_espace_virtuel();
  return derivee_seconde;
}

DoubleTab& Modele_Jones_Launder_VDF::Calcul_F1( DoubleTab& F1, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  int nb_elem = le_dom.nb_elem();
  for (int elem=0; elem <nb_elem; elem ++ )
    F1[elem] = 1.;
  return F1;
}

DoubleTab& Modele_Jones_Launder_VDF::Calcul_F2( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis& domaine_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  int nb_elem = le_dom.nb_elem();
  double Re;
  int elem;

  for (elem=0; elem< nb_elem ; elem++)
    {
      if (!is_visco_const)
        visco=tab_visco[elem];
      if (visco>0)
        {
          Re = K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0)/(K_eps_Bas_Re(elem,1)+DMINFLOAT)/visco;
          F2[elem] = 1. - (0.3*exp(-1.*carre(Re)));
        }
      else
        F2[elem] = 1.;
    }
  // F2=1;
//  Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}
/*
  DoubleTab& Modele_Jones_Launder_VDF::Calcul_F2( DoubleTab& F2, DoubleTab& D, const Domaine_dis& domaine_dis,const DoubleTab& K_eps_Bas_Re, const DoubleTab& tab_visco ) const
  {
  exit();
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  int nb_elem = le_dom.nb_elem();
  DoubleTab Re(nb_elem);
  int elem;

  for (elem=0; elem< nb_elem ; elem++) {
  if (K_eps_Bas_Re(elem,1)>0) {
  Re(elem) = (K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0))/(tab_visco(elem)*K_eps_Bas_Re(elem,1));
  F2[elem] = 1. - (0.3*exp(-1*carre(Re(elem))));
  } else {
  F2[elem] = 1.;
  }
  }
  return F2;
  }
*/

DoubleTab&  Modele_Jones_Launder_VDF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis& domaine_dis,const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  Fmu = 0;
  int nb_elem = le_dom.nb_elem();
  double Re;
  int elem;
  //  Cerr << " Calc Fmu " << finl;
  for (elem=0; elem< nb_elem ; elem++)
    {
      if (!is_visco_const)
        visco=tab_visco[elem];
      if (visco>0)
        {
          Re = K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0)/(K_eps_Bas_Re(elem,1)+DMINFLOAT)/visco;
          Fmu[elem] = exp(-2.5/(1.+Re/50.));
        }
      else
        Fmu[elem] = 1;
      // provisoire
      //        Fmu[elem] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)));
    }
  //Cerr<<Fmu.mp_min_vect()<<" fmuuuuuuuuuuuuuuu " <<Fmu.mp_max_vect()<<finl;
  /*
    Cerr<<K_eps_Bas_Re(0,0)<<" ke "<<K_eps_Bas_Re(0,1)<<finl;
    Cerr<<"re "<< K_eps_Bas_Re(0,0)*K_eps_Bas_Re(0,0)/K_eps_Bas_Re(0,1)/visco;
    Cerr<<"visco "<<visco<<finl;
  */
  // Fmu=1;
  return Fmu;
}
/*
  DoubleTab&  Modele_Jones_Launder_VDF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis& domaine_dis,const DoubleTab& K_eps_Bas_Re,const DoubleTab& tab_visco ) const
  {
  exit();
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  double Re;
  Fmu = 0;
  int nb_elem = le_dom.nb_elem();
  int elem;
  for (elem=0; elem< nb_elem ; elem++) {
  //     Fmu[elem] = exp(-2.5/(1.+K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0)/(tab_visco(elem)*K_eps_Bas_Re(elem,1))));
  Re = K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0)/K_eps_Bas_Re(elem,1)/tab_visco(elem);
  Fmu[elem] = exp(-2.5/(1.+Re/50.));
  Fmu[elem] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)));

  }
  return Fmu;
  }
*/
void  Modele_Jones_Launder_VDF::mettre_a_jour(double temps)
{
  ;
}

DoubleTab&  Modele_Jones_Launder_VDF::Calcul_Fmu_BiK( DoubleTab& Fmu,const Domaine_dis& domaine_dis,const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  Fmu = 0;
  int nb_elem = le_dom.nb_elem();
  double Re;
  int elem;
  //  Cerr << " Calc Fmu " << finl;
  for (elem=0; elem< nb_elem ; elem++)
    {
      if (!is_visco_const)
        visco=tab_visco[elem];
      if (visco>0)
        {
          Re = K_Bas_Re(elem)*K_Bas_Re(elem)/(eps_Bas_Re(elem)+DMINFLOAT)/visco;
          Fmu[elem] = exp(-2.5/(1.+Re/50.));
        }
      else
        Fmu[elem] = 1;
      // provisoire
      //        Fmu[elem] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)));
    }
  //Cerr<<Fmu.mp_min_vect()<<" fmuuuuuuuuuuuuuuu " <<Fmu.mp_max_vect()<<finl;
  /*
    Cerr<<K_eps_Bas_Re(0,0)<<" ke "<<K_eps_Bas_Re(0,1)<<finl;
    Cerr<<"re "<< K_eps_Bas_Re(0,0)*K_eps_Bas_Re(0,0)/K_eps_Bas_Re(0,1)/visco;
    Cerr<<"visco "<<visco<<finl;
  */
  // Fmu=1;
  return Fmu;
}


DoubleTab& Modele_Jones_Launder_VDF::Calcul_F2_BiK( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis& domaine_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  int nb_elem = le_dom.nb_elem();
  double Re;
  int elem;

  for (elem=0; elem< nb_elem ; elem++)
    {
      if (!is_visco_const)
        visco=tab_visco[elem];
      if (visco>0)
        {
          Re = K_Bas_Re(elem)*K_Bas_Re(elem)/(eps_Bas_Re(elem)+DMINFLOAT)/visco;
          F2[elem] = 1. - (0.3*exp(-1.*carre(Re)));
        }
      else
        F2[elem] = 1.;
    }
  // F2=1;
//  Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}



DoubleTab& Modele_Jones_Launder_VDF::Calcul_F1_BiK( DoubleTab& F1, const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_Bas_Re, const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  int nb_elem = le_dom.nb_elem();
  for (int elem=0; elem <nb_elem; elem ++ )
    F1[elem] = 1.;
  return F1;
}


DoubleTab& Modele_Jones_Launder_VDF::Calcul_E_BiK(DoubleTab& E,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& vit,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  return Calcul_E( E, domaine_dis, domaine_Cl_dis, vit, K_Bas_Re, ch_visco, visco_turb );
}


DoubleTab& Modele_Jones_Launder_VDF::Calcul_D_BiK(DoubleTab& D,const Domaine_dis& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                  const DoubleTab& vitesse,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);

  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis.valeur());
  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis.valeur());
  D = 0;
  //  return D;
  //  const DoubleVect& volumes = le_dom.volumes();
  const DoubleVect& porosite_surf = domaine_Cl_dis->equation().milieu().porosite_face();
  const DoubleVect& volume_entrelaces = le_dom.volumes_entrelaces();
  //  int nb_elem = le_dom.nb_elem();
  int nb_elem_tot = le_dom.nb_elem_tot();
  const Domaine& domaine=le_dom.domaine();

  int nb_faces_elem = domaine.nb_faces_elem();
  IntTrav numfa(nb_faces_elem);
  double coef;
  //  const IntTab& elem_faces = le_dom.elem_faces();
  const IntTab& face_voisins = le_dom.face_voisins();
  int nb_faces = le_dom.nb_faces();

  double gradk;
  int num_face,poly1,poly2,ori, ndeb, nfin;

  // Calcul de Gradient de racine de K.
  if (mp_min_vect(K_Bas_Re)<0)
    {
      Cerr << "Il y'a des valeurs negatives dans les valeurs de K" << finl;
      Cerr << "dans Modele_Jones_Launder_VDF::Calcul_D" << finl;
      Cerr << "On arrete le calcul." << finl;
      exit();
    }
  // Boucle sur les bords pour traiter les conditions aux limites
  for (int n_bord=0; n_bord<le_dom.nb_front_Cl(); n_bord++)

    {
      const Cond_lim& la_cl = le_dom_Cl.les_conditions_limites(n_bord);

      if ( sub_type(Dirichlet,la_cl.valeur())||
           sub_type(Dirichlet_homogene,la_cl.valeur()) ||
           sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ||
           sub_type(Echange_externe_impose,la_cl.valeur())
         )


        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              gradk = 0;
              poly1 =  face_voisins(num_face,0);
              if (poly1 != -1)
                {
                  // coef = 0.5;
                  coef = volume_entrelaces(num_face)*porosite_surf(num_face)*0.5;
                  gradk = ( - sqrt(K_Bas_Re(poly1)))/le_dom.dist_norm_bord(num_face);
                  if (!is_visco_const)
                    visco=tab_visco[poly1];
                  D[poly1] += 2*visco*(gradk*gradk)*coef;
                }
              else
                {
                  poly2 = face_voisins(num_face,1);
                  // coef = 0.5;
                  coef = volume_entrelaces(num_face)*porosite_surf(num_face)*0.5;
                  gradk = ((sqrt(K_Bas_Re(poly2)) ))/le_dom.dist_norm_bord(num_face);
                  //
                  if (!is_visco_const)
                    visco=tab_visco[poly2];
                  D[poly2] += 2*visco*(gradk*gradk)*coef;
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              gradk = 0;
              poly1 = face_voisins(num_face,0);
              poly2 = face_voisins(num_face,1);
              ori = le_dom.orientation(num_face);
              // coef = 0.5;

              coef = volume_entrelaces(num_face)*porosite_surf(num_face);


              gradk =  (sqrt(K_Bas_Re(poly2))-sqrt(K_Bas_Re(poly1)))/le_dom.dist_elem_period(poly1,poly2,ori);
              if (!is_visco_const)
                visco=tab_visco[poly1];
              D[poly1] += 2*visco*(gradk*gradk)*coef;
              if (!is_visco_const)
                visco=tab_visco[poly2];
              D[poly2] += 2*visco*(gradk*gradk)*coef;
            }

        }
      else if (sub_type(Symetrie,la_cl.valeur()))
        ;
      else if ( (sub_type(Neumann,la_cl.valeur()))
                ||
                (sub_type(Neumann_homogene,la_cl.valeur()))
              )
        {
          // do nothing
          ;
        }
      else
        {
          Cerr<<la_cl.valeur().que_suis_je()<< "not implemented in calculer_D"<<finl;
          exit();
        }
    }

  // Traitement des faces internes
  for (num_face=le_dom.premiere_face_int(); num_face<nb_faces; num_face++)
    {
      poly1 = face_voisins(num_face,0);
      poly2 = face_voisins(num_face,1);
      ori = le_dom.orientation(num_face);
      // coef = 0.5;

      coef = volume_entrelaces(num_face)*porosite_surf(num_face);

      gradk =  (sqrt(K_Bas_Re(poly2))-sqrt(K_Bas_Re(poly1)))/(le_dom.xp(poly2,ori)- le_dom.xp(poly1,ori));
      //      Cerr<<" ici "<< num_face<< " "<<gradk*gradk/K_eps_Bas_Re(poly2,0)<<" K "<<K_eps_Bas_Re(poly2,0)/K_eps_Bas_Re(0,0)<<finl;
      if (num_face==-396)
        Cerr << "K_eps_Bas_Re(poly2,0)=" << K_Bas_Re(poly2)/K_Bas_Re(0) << " K_eps_Bas_Re(poly1,0)=" << K_Bas_Re(poly1)/K_Bas_Re(0) << " test "<<sqrt(  K_Bas_Re(poly2,0))/sqrt( K_Bas_Re(0))<<" "<<sqrt(  K_Bas_Re(poly1))/sqrt( K_Bas_Re(0))<< " "<<(sqrt(K_Bas_Re(poly2))-sqrt(K_Bas_Re(poly1)))/sqrt( K_Bas_Re(0))<<" "<< K_Bas_Re(0)<<finl;
      if (!is_visco_const)
        visco=tab_visco[poly1];
      D[poly1] += 2*visco*(gradk*gradk)*coef;
      if (!is_visco_const)
        visco=tab_visco[poly2];
      D[poly2] += 2*visco*(gradk*gradk)*coef;
    }
  // GF on a nb_face_elem contributions par elem ....
  // mais 2 contributions dans chaque direction

  // provisoire
  D/=2;
  const DoubleVect& volumes= le_dom.volumes();
  for (int i=0; i<nb_elem_tot; i++)
    D(i)/=volumes(i);
  //Cerr<<D.mp_min_vect()<<" DDDDDDDDDDDDDD "<<D.mp_max_vect()<<finl;
  //D=0;
  return D;
  // D abord sans le 1/3 2/3
}

