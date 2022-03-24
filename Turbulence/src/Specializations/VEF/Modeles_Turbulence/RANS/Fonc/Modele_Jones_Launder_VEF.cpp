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
// File:        Modele_Jones_Launder_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Jones_Launder_VEF.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Champ_Uniforme.h>
#include <Scatter.h>
#include <Champ_P1NC.h>
#include <Champ_P0_VEF.h>

Implemente_instanciable(Modele_Jones_Launder_VEF,"Modele_Jones_Launder_VEF",Modele_Fonc_Bas_Reynolds_Base);

// printOn et readOn

Sortie& Modele_Jones_Launder_VEF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_Jones_Launder_VEF::readOn(Entree& is )
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

Entree& Modele_Jones_Launder_VEF::lire(const Motcle& , Entree& is)
{
  return is;
}
///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////

void  Modele_Jones_Launder_VEF::associer(const Zone_dis& zone_dis,
                                         const Zone_Cl_dis& zone_Cl_dis)
{
  //  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  //  const Zone_Cl_VEF& la_zone_Cl = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());
}

DoubleTab& Modele_Jones_Launder_VEF::Calcul_D(DoubleTab& D,const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,
                                              const DoubleTab& vitesse,const DoubleTab& K_eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  //  const Zone_Cl_VEF& la_zone_Cl = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());
  const DoubleVect& volumes = la_zone.volumes();
  //  int nb_faces = la_zone.nb_faces();
  int nb_faces_tot = la_zone.nb_faces();
  //  int nb_elem = la_zone.nb_elem();
  int nb_elem_tot = la_zone.nb_elem_tot();
  //  int nb_elem_tot = la_zone.nb_elem_tot();
  const Zone& zone=la_zone.zone();
  const IntTab& elem_faces = la_zone.elem_faces();
  const int nb_faces_elem = zone.nb_faces_elem();
  const DoubleTab& face_normales = la_zone.face_normales();
  const DoubleVect& vol_ent = la_zone.volumes_entrelaces();

  D = 0;
  //return D;
  int num_elem,i,face_loc,face_glob,num_face;
  double deriv,Vol;

  // Algo :
  //   * Boucle sur les elements
  //      * boucle locale dans l'element sur les faces -> calcul du gradient de racine de k
  //      * distribution de la valeur de l'integrale sur les faces de l'element
  //  Cela doit etre OK!!!
  // Rq : k et epsilon sont definis comme la vitesse P1NC.

  // ET le // ????? nb_elem ou nb_elem_tot
  // Pbls des C.L. en paroi???? non car on impose epsilon-D = 0!!
  for (num_elem=0; num_elem<nb_elem_tot; num_elem++)
    {
      Vol = volumes(num_elem);
      for(i=0; i<dimension; i++)
        {
          deriv = 0;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              deriv += sqrt(K_eps_Bas_Re(face_glob,0))*face_normales(face_glob,i)*la_zone.oriente_normale(face_glob,num_elem);
            }
          deriv /= Vol;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              /*              Cerr<<"face_glob = "<<face_glob<<finl;
                              Cerr<<"K_eps_Bas_Re(face_glob,0) = "<<K_eps_Bas_Re(face_glob,0)<<finl;
                              Cerr<<"deriv = "<<deriv<<", Vol = "<<Vol<<", nb_faces_elem = "<<nb_faces_elem<<finl;*/
              if (!is_visco_const)
                visco=tab_visco(num_elem);
              // on divise par nb_faces_elem pour dire que chaque elem
              // contribue pour 1/3 de son volume en 2D.
              D(face_glob) += visco*deriv*deriv*Vol/nb_faces_elem;
            }
        }
    }
  // on divise par vol_ent car en 2d c 2/3 du volume. Donc en tout en 2D
  // chaque face recoit en tout 2 * (1/3) * vol / (vol * (2/3))=1 !
  for (num_face=0; num_face<nb_faces_tot; num_face++)
    D(num_face) *= 2./vol_ent(num_face);

  // RQ : il faut diviser par le volume entrelace car ce que nous calculons c est l integrale en V
  // Par contre dans les termes sources on multiplie par le volume entrelace daonc pour rester coherent avec le reste....
  //D=0;
  //  D*=0.5;
  return D;
}


// Fonction utile visc
// Pour le calcul de la derivee seconde
// 1/|v|*(Si,elem)ind_der*(Sj,elem)ind_der
/* pas appele
   static double viscA(const Zone_VEF& la_zone,int num_face,int num2,int dimension,
   int num_elem, int ind_der)
   {
   double pscal;

   pscal = la_zone.face_normales(num_face,ind_der)*la_zone.face_normales(num2,ind_der);

   if ( (la_zone.face_voisins(num_face,0) == la_zone.face_voisins(num2,0)) ||
   (la_zone.face_voisins(num_face,1) == la_zone.face_voisins(num2,1)))
   return -pscal/la_zone.volumes(num_elem);
   else
   return pscal/la_zone.volumes(num_elem);
   }
*/

DoubleTab& Modele_Jones_Launder_VEF::Calcul_E(DoubleTab& E,const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  //const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());
  const Zone_VEF& zone_VEF =  ref_cast(Zone_VEF,zone_dis.valeur());
  E = 0;
  //return E;
  //  const IntTab& elem_faces = zone_VEF.elem_faces();
  const IntTab& face_voisins = zone_VEF.face_voisins();
  const DoubleVect& volumes = zone_VEF.volumes();
  //const DoubleTab& face_normales = zone_VEF.face_normales();
  //  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  //  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  //  const Zone& zone = zone_VEF.zone();
  const int nb_faces = zone_VEF.nb_faces();
  //  int nb_faces_tot = zone_VEF.nb_faces_tot();
  //  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem = zone_VEF.nb_elem();
  //  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  //int premiere_face_int = zone_VEF.premiere_face_int();
  int premiere_face_int = zone_VEF.premiere_face_int();

  int i,j;
  int elem,elem0,elem2;
  //int n_bord;//alpha,beta,elem0;
  //  int nb_faces_elem = zone_VEF.zone().nb_faces_elem();

  //  double val;

  const int ncomp_ch_transporte = transporte.line_size();

  DoubleTab gradient(0, ncomp_ch_transporte, dimension);
  zone_VEF.creer_tableau_faces(gradient);

  DoubleTab deriv_seconde_de_gradient_elem(0, ncomp_ch_transporte, dimension, dimension);
  zone_VEF.zone().creer_tableau_elements(deriv_seconde_de_gradient_elem);
  deriv_seconde_de_gradient_elem=0.;

  DoubleTab champ_face(0,dimension);
  zone_VEF.creer_tableau_faces(champ_face);

  // Calcul de la derivee premiere de U
  // stokage d'un gradient_elem par element.
  DoubleTab gradient_elem(0, ncomp_ch_transporte, dimension);
  zone_VEF.zone().creer_tableau_elements(gradient_elem);

  DoubleTab champ_elem(0, ncomp_ch_transporte, dimension);
  zone_VEF.zone().creer_tableau_elements(champ_elem);

  DoubleTab deriv_seconde_de_gradient(0, ncomp_ch_transporte, dimension, dimension);
  zone_VEF.creer_tableau_faces(deriv_seconde_de_gradient);
  deriv_seconde_de_gradient=0.;

  int fac=-1,elem1,comp;
  // Rque methode non const Pourquoi ?

  const Champ_P1NC& vitesse = ref_cast(Champ_P1NC,eq_hydraulique->inconnue().valeur());
  ref_cast_non_const(Champ_P1NC,vitesse).calcul_gradient(transporte,gradient_elem,ref_cast(Zone_Cl_VEF,eq_hydraulique->zone_Cl_dis().valeur()));

  gradient_elem.echange_espace_virtuel();
  // On a les gradient_elem par elements
  // Boucle sur les faces
  //
  for (fac=0; fac< premiere_face_int; fac++)
    for (comp=0; comp<ncomp_ch_transporte; comp++)
      for (i=0; i<dimension; i++)
        {
          elem1=face_voisins(fac,0);
          if (elem1 != -1)
            gradient(fac, comp, i) = gradient_elem(elem1, comp, i);
          else
            {
              elem1=face_voisins(fac,1);
              gradient(fac, comp, i) = gradient_elem(elem1, comp, i);
            }
        }
// fin du for faces

  for (; fac<nb_faces; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      double vol1=volumes(elem1);
      double vol2=volumes(elem2);
      double voltot=vol1+vol2;
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (i=0; i<dimension; i++)
          {
            double grad1=gradient_elem(elem1, comp, i);
            double grad2=gradient_elem(elem2, comp, i);
            gradient(fac, comp, i) =  (vol1*grad1 + vol2*grad2)/voltot;
          }
    } // fin du for faces

  gradient.echange_espace_virtuel();
//on a les gradients sur les faces.


  // maintenant on calcul la derivee seconde
  // stokage de la derivee seconde par element.
  for (i=0; i<dimension; i++)
    {
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (fac=0; fac<nb_faces; fac++)
          champ_face(fac,comp)=gradient(fac, comp, i);

      ref_cast_non_const(Champ_P1NC,vitesse).calcul_gradient(champ_face,champ_elem,ref_cast(Zone_Cl_VEF,eq_hydraulique->zone_Cl_dis().valeur()));

      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for ( elem=0; elem<nb_elem; elem++)
          for (j=0; j<dimension; j++)
            deriv_seconde_de_gradient_elem(elem,comp,i,j) = champ_elem(elem,comp,j);

    }

  deriv_seconde_de_gradient_elem.echange_espace_virtuel();
  // On a les gradient_elem par elements

  // Boucle sur les faces
  //

  for (fac=0; fac< premiere_face_int; fac++)
    {
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (j=0; j<dimension; j++)
          for (i=0; i<dimension; i++)
            {
              elem1=face_voisins(fac,0);
              if (elem1 != -1)
                deriv_seconde_de_gradient(fac, comp, i, j) = deriv_seconde_de_gradient_elem(elem1, comp, i, j);
              else
                {
                  elem1=face_voisins(fac,1);
                  deriv_seconde_de_gradient(fac, comp, i, j) = deriv_seconde_de_gradient_elem(elem1, comp, i, j);
                }
            }

    } // fin du for faces

  for (; fac<nb_faces; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      double vol1=volumes(elem1);
      double vol2=volumes(elem2);
      double voltot=vol1+vol2;
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (j=0; j<dimension; j++)
          for (i=0; i<dimension; i++)
            {
              double grad1=deriv_seconde_de_gradient_elem(elem1, comp, i, j);
              double grad2=deriv_seconde_de_gradient_elem(elem2, comp, i, j);
              deriv_seconde_de_gradient(fac, comp, i, j) = (vol1*grad1 + vol2*grad2)/voltot;
            }
    } // fin du for faces

  // Calcul de E
  double deriv,nuturb;
  for (fac=0; fac<nb_faces; fac++)
    {
      deriv = 0;
      if (dimension == 2)
        {
          double val2 = deriv_seconde_de_gradient(fac,0,1,1)*deriv_seconde_de_gradient(fac,0,1,1);
          double val3 = deriv_seconde_de_gradient(fac,1,0,0)*deriv_seconde_de_gradient(fac,1,0,0);
          deriv = val2 + val3;
          //deriv = gradient(fac, 0, 1);
        }
      else
        {
          double val2 = deriv_seconde_de_gradient(fac,1,0,0)+deriv_seconde_de_gradient(fac,2,0,0); //d2v_dx2(fac)+d2w_dx2(fac);
          val2 *= val2;
          double val3 = deriv_seconde_de_gradient(fac,0,1,1)+deriv_seconde_de_gradient(fac,2,1,1); //d2u_dy2(fac)+d2w_dy2(fac);
          val3 *= val3;
          double val4 = deriv_seconde_de_gradient(fac,0,2,2)+deriv_seconde_de_gradient(fac,1,2,2); //d2u_dz2(fac)+d2v_dz2(fac);
          val4 *= val4;
          deriv = val2 + val3 + val4;

        }
      elem0 = face_voisins(fac,0);
      elem1 = face_voisins(fac,1);

      if (elem1!=-1)
        {
          nuturb = visco_turb(elem0)*volumes(elem0)+visco_turb(elem1)*volumes(elem1);
          nuturb /= volumes(elem0) + volumes(elem1);
          if (!is_visco_const)
            visco=(tab_visco(elem0)+tab_visco(elem1))/2.;
        }
      else
        {
          nuturb =  visco_turb(elem0);
          if (!is_visco_const)
            visco=(tab_visco(elem0));
        }

      E(fac) = 2.*visco*deriv*nuturb;

    }
  //E=0;
//  Cerr<<E.mp_min_vect()<<" EEEEEEEEEEEEEEEEEEEEEEEE "<<E.mp_max_vect()<<finl;

  return E;
}

DoubleTab& Modele_Jones_Launder_VEF::Calcul_F1( DoubleTab& F1, const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  for (int num_face=0; num_face <nb_faces; num_face ++ )
    F1[num_face] = 1.;
  return F1;
}

DoubleTab& Modele_Jones_Launder_VEF::Calcul_F2( DoubleTab& F2, DoubleTab& Deb, const Zone_dis& zone_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face;
  double Re;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (K_eps_Bas_Re(num_face,1)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = la_zone.face_voisins(num_face,0);
              int elem1 = la_zone.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
                  visco /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);

            }
          if (visco>0)
            {
              Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*K_eps_Bas_Re(num_face,1));
              F2[num_face] = 1. - (0.3*exp(-Re*Re));
            }
        }
      else
        {
          F2[num_face] = 1.;
        }
    }
  //Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}
/*
  DoubleTab& Modele_Jones_Launder_VEF::Calcul_F2( DoubleTab& F2, DoubleTab& D, const Zone_dis& zone_dis,const DoubleTab& K_eps_Bas_Re, const DoubleTab& tab_visco ) const
  {
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face,elem0,elem1;
  double Re,nulam;

  for (num_face=0; num_face<nb_faces  ; num_face++)
  {
  elem0 = la_zone.face_voisins(num_face,0);
  elem1 = la_zone.face_voisins(num_face,1);
  if (elem1!=-1)
  {
  nulam = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
  nulam /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
  }
  else
  nulam =  tab_visco(elem0);

  if (K_eps_Bas_Re(num_face,1)>0)
  {
  Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(nulam*K_eps_Bas_Re(num_face,1));
  F2[num_face] = 1. - (0.3*exp(-Re*Re));
  }
  else
  {
  F2[num_face] = 1.;
  }
  }
  return F2;
  }

*/
DoubleTab&  Modele_Jones_Launder_VEF::Calcul_Fmu( DoubleTab& Fmu,const Zone_dis& zone_dis,const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face;
  double Re;
  Fmu = 0;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (1) //   if (K_eps_Bas_Re(num_face,1)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = la_zone.face_voisins(num_face,0);
              int elem1 = la_zone.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
                  visco /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);
            }
          Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*(K_eps_Bas_Re(num_face,1)+DMINFLOAT));
          Fmu[num_face] = exp(-2.5/(1.+Re/50));
        }
      else
        {
          Fmu[num_face] = 0.;
        }
    }

  return Fmu;
}

/*
  DoubleTab&  Modele_Jones_Launder_VEF::Calcul_Fmu( DoubleTab& Fmu,const Zone_dis& zone_dis,const DoubleTab& K_eps_Bas_Re,const DoubleTab& tab_visco ) const
  {
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face,elem0,elem1;
  double Re,nulam;
  Fmu = 0;

  for (num_face=0; num_face<nb_faces  ; num_face++)
  {
  elem0 = la_zone.face_voisins(num_face,0);
  elem1 = la_zone.face_voisins(num_face,1);
  if (elem1!=-1)
  {
  nulam = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
  nulam /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
  }
  else
  nulam =  tab_visco(elem0);

  if (K_eps_Bas_Re(num_face,1)>0)
  {
  Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(nulam*K_eps_Bas_Re(num_face,1));
  Fmu[num_face] = exp(-2.5/(1.+Re/50.));
  }
  else
  {
  Fmu[num_face] = 0.;
  }
  }

  return Fmu;
  }
*/
void  Modele_Jones_Launder_VEF::mettre_a_jour(double temps)
{
  ;
}

DoubleTab&  Modele_Jones_Launder_VEF::Calcul_Fmu_BiK( DoubleTab& Fmu,const Zone_dis& zone_dis,const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face;
  double Re;
  Fmu = 0;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (1) //   if (eps_Bas_Re(num_face)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = la_zone.face_voisins(num_face,0);
              int elem1 = la_zone.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
                  visco /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);
            }
          Re = (K_Bas_Re(num_face)*K_Bas_Re(num_face))/(visco*(eps_Bas_Re(num_face)+DMINFLOAT));
          Fmu[num_face] = exp(-2.5/(1.+Re/50));
        }
      else
        {
          Fmu[num_face] = 0.;
        }
    }

  return Fmu;
}


DoubleTab& Modele_Jones_Launder_VEF::Calcul_F2_BiK( DoubleTab& F2, DoubleTab& Deb, const Zone_dis& zone_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  int num_face;
  double Re;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (eps_Bas_Re(num_face)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = la_zone.face_voisins(num_face,0);
              int elem1 = la_zone.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*la_zone.volumes(elem0)+tab_visco(elem1)*la_zone.volumes(elem1);
                  visco /= la_zone.volumes(elem0) + la_zone.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);

            }
          if (visco>0)
            {
              Re = (K_Bas_Re(num_face)*K_Bas_Re(num_face))/(visco*eps_Bas_Re(num_face));
              F2[num_face] = 1. - (0.3*exp(-Re*Re));
            }
        }
      else
        {
          F2[num_face] = 1.;
        }
    }
  //Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}

DoubleTab& Modele_Jones_Launder_VEF::Calcul_F1_BiK( DoubleTab& F1, const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& P, const DoubleTab& K_Bas_Re, const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  int nb_faces = la_zone.nb_faces();
  for (int num_face=0; num_face <nb_faces; num_face ++ )
    F1[num_face] = 1.;
  return F1;
}


DoubleTab& Modele_Jones_Launder_VEF::Calcul_E_BiK(DoubleTab& E,const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  return Calcul_E( E, zone_dis, zone_Cl_dis, transporte, K_Bas_Re, ch_visco, visco_turb );
}

DoubleTab& Modele_Jones_Launder_VEF::Calcul_D_BiK(DoubleTab& D,const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,
                                                  const DoubleTab& vitesse,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Zone_VEF& la_zone = ref_cast(Zone_VEF,zone_dis.valeur());
  //  const Zone_Cl_VEF& la_zone_Cl = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());
  const DoubleVect& volumes = la_zone.volumes();
  //  int nb_faces = la_zone.nb_faces();
  int nb_faces_tot = la_zone.nb_faces();
  //  int nb_elem = la_zone.nb_elem();
  int nb_elem_tot = la_zone.nb_elem_tot();
  //  int nb_elem_tot = la_zone.nb_elem_tot();
  const Zone& zone=la_zone.zone();
  const IntTab& elem_faces = la_zone.elem_faces();
  const int nb_faces_elem = zone.nb_faces_elem();
  const DoubleTab& face_normales = la_zone.face_normales();
  const DoubleVect& vol_ent = la_zone.volumes_entrelaces();

  D = 0;
  //return D;
  int num_elem,i,face_loc,face_glob,num_face;
  double deriv,Vol;

  // Algo :
  //   * Boucle sur les elements
  //      * boucle locale dans l'element sur les faces -> calcul du gradient de racine de k
  //      * distribution de la valeur de l'integrale sur les faces de l'element
  //  Cela doit etre OK!!!
  // Rq : k et epsilon sont definis comme la vitesse P1NC.

  // ET le // ????? nb_elem ou nb_elem_tot
  // Pbls des C.L. en paroi???? non car on impose epsilon-D = 0!!
  for (num_elem=0; num_elem<nb_elem_tot; num_elem++)
    {
      Vol = volumes(num_elem);
      for(i=0; i<dimension; i++)
        {
          deriv = 0;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              deriv += sqrt(K_Bas_Re(face_glob))*face_normales(face_glob,i)*la_zone.oriente_normale(face_glob,num_elem);
            }
          deriv /= Vol;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              /*              Cerr<<"face_glob = "<<face_glob<<finl;
                              Cerr<<"K_eps_Bas_Re(face_glob,0) = "<<K_eps_Bas_Re(face_glob,0)<<finl;
                              Cerr<<"deriv = "<<deriv<<", Vol = "<<Vol<<", nb_faces_elem = "<<nb_faces_elem<<finl;*/
              if (!is_visco_const)
                visco=tab_visco(num_elem);
              // on divise par nb_faces_elem pour dire que chaque elem
              // contribue pour 1/3 de son volume en 2D.
              D(face_glob) += visco*deriv*deriv*Vol/nb_faces_elem;
            }
        }
    }
  // on divise par vol_ent car en 2d c 2/3 du volume. Donc en tout en 2D
  // chaque face recoit en tout 2 * (1/3) * vol / (vol * (2/3))=1 !
  for (num_face=0; num_face<nb_faces_tot; num_face++)
    D(num_face) *= 2./vol_ent(num_face);

  // RQ : il faut diviser par le volume entrelace car ce que nous calculons c est l integrale en V
  // Par contre dans les termes sources on multiplie par le volume entrelace daonc pour rester coherent avec le reste....
  //D=0;
  //  D*=0.5;
  return D;
}



