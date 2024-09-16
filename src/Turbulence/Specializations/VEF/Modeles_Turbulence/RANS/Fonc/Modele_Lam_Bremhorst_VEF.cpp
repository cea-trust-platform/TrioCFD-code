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
// File:        Modele_Lam_Bremhorst_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////
#include <Modele_Lam_Bremhorst_VEF.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Periodique.h>
#include <Champ_Uniforme.h>
#include <Scatter.h>
#include <Champ_P1NC.h>
#include <Champ_P0_VEF.h>
#include <Discretisation_base.h>
#include <Check_espace_virtuel.h>
#include <LecFicDiffuse.h>
#include <EcritureLectureSpecial.h>

Implemente_instanciable(Modele_Lam_Bremhorst_VEF,"Modele_Lam_Bremhorst_VEF",Modele_Fonc_Bas_Reynolds_Base);
// XD Lam_Bremhorst modele_fonction_bas_reynolds_base Lam_Bremhorst -1 Model described in ' C.K.G.Lam and K.Bremhorst, A modified form of the k- epsilon model for predicting wall turbulence, ASME J. Fluids Engng., Vol.103, p456, (1981)'. Only in VEF.
// printOn et readOn

Sortie& Modele_Lam_Bremhorst_VEF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_Lam_Bremhorst_VEF::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  if (is_Reynolds_stress_isotrope_ == 0)
    is_Cmu_constant_ = 0;
  return is;
}

void Modele_Lam_Bremhorst_VEF::set_param(Param& param)
{

  is_Cmu_constant_ = 1;
  is_Reynolds_stress_isotrope_ = 1;
  param.ajouter("fichier_distance_paroi",&nom_fic);    // XD_ADD_P chaine refer to distance_paroi keyword
  param.ajouter("Reynolds_stress_isotrope",&is_Reynolds_stress_isotrope_); // XD_ADD_P int keyword for isotropic Reynolds stress
}

int Modele_Lam_Bremhorst_VEF::Calcul_is_Reynolds_stress_isotrope() const
{
  return is_Reynolds_stress_isotrope_;
}

int Modele_Lam_Bremhorst_VEF::Calcul_is_Cmu_constant() const
{
  return is_Cmu_constant_;
}

Entree& Modele_Lam_Bremhorst_VEF::lire(const Motcle& , Entree& is)
{
  return is;
}
///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////

void  Modele_Lam_Bremhorst_VEF::associer(const Domaine_dis_base& domaine_dis,
                                         const Domaine_Cl_dis& domaine_Cl_dis)
{
  le_dom_VEF = ref_cast(Domaine_VEF,domaine_dis);
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());
}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_D(DoubleTab& D,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                              const DoubleTab& vitesse,const DoubleTab& K_eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  D = 0;
  return D;
}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_E(DoubleTab& E,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  E = 0;
  return E;
}

/* DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_F1( DoubleTab& F1, const Domaine_dis_base& domaine_dis) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  for (int num_face=0; num_face <nb_faces; num_face ++ )
    F1[num_face] = 1.;
  return F1;
}
*/
DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_F1( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{

  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());
  const DoubleTab& wall_length = BR_wall_length_->valeurs();
  DoubleTab wall_length_face(0);
  le_dom.creer_tableau_faces(wall_length_face);
  DoubleTab Pderive(0);
  le_dom.creer_tableau_faces(Pderive);
  DoubleTab Fmu_loc(0);
  le_dom.creer_tableau_faces(Fmu_loc);
  int nb_faces = le_dom.nb_faces();
  const Conds_lim& les_cl = domaine_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const IntTab& face_voisins = le_dom.face_voisins();
  int num_face;
  double Rey,Re;
  /*
  	for (num_face=0; num_face <nb_faces; num_face ++ )
  		F1[num_face] = 1.;
  	return F1;
  */
  // Calcul de la distance a la paroi aux faces
  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              int elem2 = face_voisins(num_face,1);
              wall_length_face(num_face) = 0.5*wall_length(elem1) + 0.5*wall_length(elem2);
            }
        }
      /* else if (sub_type(Pa,la_cl.valeur()))
      {
        for (int num_face=ndeb; num_face<nfin; num_face++)
        {
      	  int elem1 = face_voisins(num_face,0);
      	  int elem2 = face_voisins(num_face,1);
      			  wall_length_face(num_face) = 0.;
        }
      }*/
      else
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              wall_length_face(num_face) = wall_length(elem1);
            }
        }
    }
  int n0 = le_dom.premiere_face_int();
  for (num_face=n0; num_face<nb_faces; num_face++)
    {
      int elem0 = le_dom.face_voisins(num_face,0);
      int elem1 = le_dom.face_voisins(num_face,1);
      wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
    }
  // Calcul de la distance a la paroi aux faces
  /*    for (num_face=0; num_face< le_dom.premiere_face_int(); num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  if (elem0 != -1)
    		  wall_length_face(num_face) = wall_length(elem0);
    	  else
    	  {
    		  elem0 = le_dom.face_voisins(num_face,1);
    		  wall_length_face(num_face) = wall_length(elem0);
    	  }
      }

      for (; num_face<nb_faces; num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  int elem1 = le_dom.face_voisins(num_face,1);
    	  wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
      }
  */
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      if (visco>BR_EPS && K_eps_Bas_Re(num_face,1)>BR_EPS && K_eps_Bas_Re(num_face,0))
        {
          Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*K_eps_Bas_Re(num_face,1));
          Rey = wall_length_face(num_face)*sqrt(K_eps_Bas_Re(num_face,0))/visco;
          Fmu_loc(num_face) = (1. - exp(-0.0165*Rey-BR_EPS))*(1. - exp(-0.0165*Rey-BR_EPS))*(1.+20.5/(Re+BR_EPS));
        }
      else
        Fmu_loc(num_face) = 1.;

      if (Fmu_loc(num_face) > 1. + 1.e-8)
        {
          Fmu_loc(num_face) = 1.;
          //Cerr <<  " On force Fmu_loc a 1 " << finl;
          //exit();
        }
    }
//  Cerr<< Fmu_loc.mp_min_vect()<< " Fmu_loc "<<Fmu_loc.mp_max_vect()<<finl;

  // Calcul de F1
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      if (Fmu_loc(num_face) > BR_EPS && (visco>BR_EPS && K_eps_Bas_Re(num_face,1)>BR_EPS && K_eps_Bas_Re(num_face,0)))
        F1[num_face] = 1.+ (0.05/(Fmu_loc(num_face)+BR_EPS))*(0.05/(Fmu_loc(num_face)+BR_EPS))*(0.05/(Fmu_loc(num_face)+BR_EPS));
      else
        F1[num_face] = 1.;
    }
//  Cerr<< F1.mp_min_vect()<< " F1 "<<F1.mp_max_vect()<<finl;
  return F1;
}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_F2( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis_base& domaine_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  int num_face;
  double Re;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (!is_visco_const)
        {
          int elem0 = le_dom.face_voisins(num_face,0);
          int elem1 = le_dom.face_voisins(num_face,1);
          if (elem1!=-1)
            {
              visco = tab_visco(elem0)*le_dom.volumes(elem0)+tab_visco(elem1)*le_dom.volumes(elem1);
              visco /= le_dom.volumes(elem0) + le_dom.volumes(elem1);
            }
          else
            visco =  tab_visco(elem0);

        }
      if (visco>BR_EPS && K_eps_Bas_Re(num_face,1)>BR_EPS)
        {
          Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*K_eps_Bas_Re(num_face,1));
          F2[num_face] = 1. - (exp(-Re*Re));
        }
      else
        F2[num_face] = 1.;

    }
//  Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}

DoubleTab&  Modele_Lam_Bremhorst_VEF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());
  const DoubleTab& wall_length = BR_wall_length_->valeurs();
  DoubleTab wall_length_face(0);
  le_dom.creer_tableau_faces(wall_length_face);
  DoubleTab Pderive(0);
  le_dom.creer_tableau_faces(Pderive);
  int nb_faces = le_dom.nb_faces();
  const Conds_lim& les_cl = domaine_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const IntTab& face_voisins = le_dom.face_voisins();
  int num_face;
  double Rey,Re;
  /*
  	for (num_face=0; num_face <nb_faces; num_face ++ )
  		F1[num_face] = 1.;
  	return F1;
  */
  // Calcul de la distance a la paroi aux faces
  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              int elem2 = face_voisins(num_face,1);
              wall_length_face(num_face) = 0.5*wall_length(elem1) + 0.5*wall_length(elem2);
            }
        }
      /* else if (sub_type(Pa,la_cl.valeur()))
      {
        for (int num_face=ndeb; num_face<nfin; num_face++)
        {
      	  int elem1 = face_voisins(num_face,0);
      	  int elem2 = face_voisins(num_face,1);
      			  wall_length_face(num_face) = 0.;
        }
      }*/
      else
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              wall_length_face(num_face) = wall_length(elem1);
            }
        }
    }
  int n0 = le_dom.premiere_face_int();
  for (num_face=n0; num_face<nb_faces; num_face++)
    {
      int elem0 = le_dom.face_voisins(num_face,0);
      int elem1 = le_dom.face_voisins(num_face,1);
      wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
    }
  // Calcul de la distance a la paroi aux faces
  /*    for (num_face=0; num_face< le_dom.premiere_face_int(); num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  if (elem0 != -1)
    		  wall_length_face(num_face) = wall_length(elem0);
    	  else
    	  {
    		  elem0 = le_dom.face_voisins(num_face,1);
    		  wall_length_face(num_face) = wall_length(elem0);
    	  }
      }

      for (; num_face<nb_faces; num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  int elem1 = le_dom.face_voisins(num_face,1);
    	  wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
      }
  */
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      if (visco>BR_EPS && K_eps_Bas_Re(num_face,1)>BR_EPS && K_eps_Bas_Re(num_face,0))
        {
          Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*K_eps_Bas_Re(num_face,1));
          Rey = wall_length_face(num_face)*sqrt(K_eps_Bas_Re(num_face,0))/visco;
          Fmu(num_face) = (1. - exp(-0.0165*Rey-BR_EPS))*(1. - exp(-0.0165*Rey-BR_EPS))*(1.+20.5/(Re+BR_EPS));

        }
      else
        {
          //Cerr <<  " visco " << visco << " Eps " << K_eps_Bas_Re(num_face,1) << " K " << K_eps_Bas_Re(num_face,0) <<finl;
          Fmu(num_face) = 1.;
        }
      if (Fmu(num_face) > 1. + 1.e-8)
        {
          //Cerr <<  " Fmu=" << Fmu(num_face) << " Eps=" << K_eps_Bas_Re(num_face,1) << " K=" << K_eps_Bas_Re(num_face,0) << " Re=" << Re << " Rey=" << Rey <<finl;
          Fmu(num_face) = 1.;
          //Cerr <<  " On force Fmu a 1 " << finl;
          //exit();
        }
    }
//  Cerr<< Fmu.mp_min_vect()<< " Fmu "<<Fmu.mp_max_vect()<<finl;

  return Fmu;
}

void  Modele_Lam_Bremhorst_VEF::mettre_a_jour(double temps)
{
  ;
}

// Initialisation d'une matrice aux elements
void Modele_Lam_Bremhorst_VEF::init_tenseur_elem(DoubleTab& Tenseur, const Domaine_VEF& domaine_VEF, const int ndim)  const
{
  if(!Tenseur.get_md_vector().non_nul())
    {
      if (ndim==1)
        Tenseur.resize(0, Objet_U::dimension);
      else if (ndim==2)
        Tenseur.resize(0, Objet_U::dimension, Objet_U::dimension);
      domaine_VEF.domaine().creer_tableau_elements(Tenseur);
    }
  Tenseur = 0.;
}


// Initialisation d'une matrice aux faces
void  Modele_Lam_Bremhorst_VEF::init_tenseur_face(DoubleTab& Tenseur,const Domaine_VEF& domaine_VEF, const int ndim) const
{
  if(!Tenseur.get_md_vector().non_nul())
    {
      if (ndim==1)
        Tenseur.resize(0, Objet_U::dimension);
      else if (ndim==2)
        Tenseur.resize(0, Objet_U::dimension, Objet_U::dimension);
      domaine_VEF.creer_tableau_faces(Tenseur);
    }
  Tenseur = 0.;
}

// Calcul d'un tenseur aux faces a partir d'un tenseur aux elements
DoubleTab& Modele_Lam_Bremhorst_VEF::calcul_tenseur_face(DoubleTab& Tenseur_face, const DoubleTab& Tenseur_elem,
                                                         const Domaine_VEF& domaine_VEF, const Domaine_Cl_VEF& domaine_Cl_VEF) const
{
  assert_espace_virtuel_vect(Tenseur_elem);
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  int nb_faces = domaine_VEF.nb_faces();

  const Conds_lim& les_cl = domaine_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const DoubleVect& volumes = domaine_VEF.volumes();

  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (int fac=ndeb; fac<nfin; fac++)
            {
              int poly1 = face_voisins(fac,0);
              int poly2 = face_voisins(fac,1);
              double a=volumes(poly1)/(volumes(poly1)+volumes(poly2));
              double b=volumes(poly2)/(volumes(poly1)+volumes(poly2));
              for (int i=0; i<dimension; i++)
                for (int j=0; j<dimension; j++)
                  Tenseur_face(fac,i,j) = a*Tenseur_elem(poly1,i,j) + b*Tenseur_elem(poly2,i,j);
            }
        }
      else
        {
          for (int fac=ndeb; fac<nfin; fac++)
            {
              int poly1 = face_voisins(fac,0);
              for (int i=0; i<dimension; i++)
                for (int j=0; j<dimension; j++)
                  Tenseur_face(fac,i,j) = Tenseur_elem(poly1,i,j);
            }
        }
    }
  int n0 = domaine_VEF.premiere_face_int();
  for (int fac = n0; fac<nb_faces; fac++)
    {
      int poly1 = face_voisins(fac,0);
      int poly2 = face_voisins(fac,1);
      double a=volumes(poly1)/(volumes(poly1)+volumes(poly2));
      double b=volumes(poly2)/(volumes(poly1)+volumes(poly2));
      for (int i=0; i<dimension; i++)
        for (int j=0; j<dimension; j++)
          Tenseur_face(fac,i,j) = a*Tenseur_elem(poly1,i,j) + b*Tenseur_elem(poly2,i,j);
    }

  return Tenseur_face;
}

bool Modele_Lam_Bremhorst_VEF::calcul_tenseur_Re(const DoubleTab& nu_turb, const DoubleTab& G, DoubleTab& Re) const
{
  const Domaine_dis_base& domaine_dis = mon_equation->domaine_dis();
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  int nelem = G.dimension(0);

  DoubleTab S, R;
  init_tenseur_elem(S,domaine_VEF,2);
  init_tenseur_elem(R,domaine_VEF,2);

  for (int elem=0; elem<nelem; elem++)
    for (int i=0; i<Objet_U::dimension; i++)
      for (int j=0; j<Objet_U::dimension; j++)
        {
          S(elem,i,j) = G(elem,i,j) + G(elem,j,i);
          R(elem,i,j) = G(elem,i,j) - G(elem,j,i);
        }

  DoubleTab Sn = calcul_norme_elem(domaine_VEF,S);

  Re = calcul_tenseur_Re_elem(mon_equation->discretisation(),domaine_dis,G,S,R,Sn,mon_equation->inconnue());

  return true;
}

bool Modele_Lam_Bremhorst_VEF::calcul_tenseur_Re_BiK(const DoubleTab& nu_turb, const DoubleTab& G, DoubleTab& Re) const
{
  const Domaine_dis_base& domaine_dis = mon_equation->domaine_dis();
  const Domaine_VEF&       domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  int nelem = G.dimension(0);

  DoubleTab S, R;
  init_tenseur_elem(S,domaine_VEF,2);
  init_tenseur_elem(R,domaine_VEF,2);

  for (int elem=0; elem<nelem; elem++)
    for (int i=0; i<Objet_U::dimension; i++)
      for (int j=0; j<Objet_U::dimension; j++)
        {
          S(elem,i,j) = G(elem,i,j) + G(elem,j,i);
          R(elem,i,j) = G(elem,i,j) - G(elem,j,i);
        }

  DoubleTab Sn = calcul_norme_elem(domaine_VEF,S);

  Re = calcul_tenseur_Re_elem_BiK(mon_equation->discretisation(),domaine_dis,G,S,R,Sn,mon_equation->inconnue(),ma_seconde_equation->inconnue());

  return true;
}

DoubleTab Modele_Lam_Bremhorst_VEF::calcul_tenseur_Re_elem(const Discretisation_base& dis,
                                                           const Domaine_dis_base& domaine_dis,
                                                           const DoubleTab& G, const DoubleTab& S,
                                                           const DoubleTab& R, const DoubleTab& Sn,
                                                           const Champ_base& K_Eps) const
{
  int nelem = G.dimension(0);

  Champ_Fonc K_Eps_elem;
  Noms noms(2), unites(2);
  noms[0]="K";
  noms[1]="eps";
  unites[0]="m2/s2";
  unites[1]="m2/s3";
  const Motcle wtype= "champ_elem";
  dis.discretiser_champ(wtype,domaine_dis,multi_scalaire,noms,unites,2,0.,K_Eps_elem);
  K_Eps_elem->affecter(K_Eps);
  DoubleTab tab_K_Eps  = K_Eps_elem->valeurs();

  DoubleTab ReNL;
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  init_tenseur_elem(ReNL,domaine_VEF,2);

  ReNL = 0.;
  for (int elem=0; elem<nelem; elem++)
    {
      double kseps;
      if (K_Eps.valeurs()(elem,1) <= BR_EPS)
        {
          kseps = K_Eps.valeurs()(elem,0)/BR_EPS;
        }
      else
        {
          kseps = tab_K_Eps(elem,0)/tab_K_Eps(elem,1);
        }
      double Cmu = (2./3.)/(A1+sqrt(0.5)*Sn(elem)*kseps);
      double C1 = CNL1/((CNL4+CNL5*pow(sqrt(0.5)*Sn(elem)*kseps,3.))*Cmu)* kseps;
      double C2 = CNL2/((CNL4+CNL5*pow(sqrt(0.5)*Sn(elem)*kseps,3.))*Cmu)* kseps;
      double C3 = CNL3/((CNL4+CNL5*pow(sqrt(0.5)*Sn(elem)*kseps,3.))*Cmu)* kseps;
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            ReNL(elem,i,j) = S(elem,i,j);
            for (int k=0; k<Objet_U::dimension; k++)
              ReNL(elem,i,j) -= C1* S(elem,i,k)*S(elem,k,j) +
                                C2*(R(elem,i,k)*S(elem,k,j) + R(elem,j,k)*S(elem,k,i)) +
                                C3 *R(elem,i,k)*R(elem,j,k);
            if (i==j)
              {
                for (int k=0; k<Objet_U::dimension; k++)
                  for (int l=0; l<Objet_U::dimension; l++)
                    ReNL(elem,i,j) += 1./3.*(C1*S(elem,k,l)*S(elem,k,l)+C3*R(elem,k,l)*R(elem,k,l));
              }
          }
    }

  return ReNL;
}

DoubleTab Modele_Lam_Bremhorst_VEF::calcul_tenseur_Re_elem_BiK(const Discretisation_base& dis,
                                                               const Domaine_dis_base& domaine_dis,
                                                               const DoubleTab& G, const DoubleTab& S,
                                                               const DoubleTab& R, const DoubleTab& Sn,
                                                               const Champ_base& K, const Champ_base& Eps) const
{
  int nelem = G.dimension(0);

  Champ_Fonc K_elem;
  Champ_Fonc Eps_elem;
  Noms noms(1), unites(1), noms2(1), unites2(1);
  noms[0]="K";
  noms2[0]="eps";
  unites[0]="m2/s2";
  unites2[0]="m2/s3";
  const Motcle wtype= "champ_elem";
  dis.discretiser_champ(wtype,domaine_dis,multi_scalaire,noms,unites,1,0.,K_elem);
  dis.discretiser_champ(wtype,domaine_dis,multi_scalaire,noms2,unites2,1,0.,Eps_elem);

  K_elem->affecter(K);
  Eps_elem->affecter(Eps);
  DoubleTab tab_K    = K_elem->valeurs();
  DoubleTab tab_Eps  = Eps_elem->valeurs();

  DoubleTab ReNL;
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  init_tenseur_elem(ReNL,domaine_VEF,2);

  ReNL = 0.;
  for (int elem=0; elem<nelem; elem++)
    {
      double kseps;
      if (Eps.valeurs()(elem,0) <= BR_EPS)
        {
          kseps = K.valeurs()(elem,0)/BR_EPS;
        }
      else
        {
          kseps = tab_K(elem,0)/tab_Eps(elem,0);
        }
      double Cmu = (2./3.)/(A1+sqrt(0.5)*Sn(elem)*kseps);
      double C1 = CNL1/((CNL4+CNL5*pow(sqrt(0.5)*Sn(elem)*kseps,3.))*Cmu)* kseps;
      double C2 = CNL2/((CNL4+CNL5*pow(sqrt(0.5)*Sn(elem)*kseps,3.))*Cmu)* kseps;
      double C3 = CNL3/((CNL4+CNL5*pow(sqrt(0.5)*Sn(elem)*kseps,3.))*Cmu)* kseps;
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            ReNL(elem,i,j) = S(elem,i,j);
            for (int k=0; k<Objet_U::dimension; k++)
              ReNL(elem,i,j) -= C1* S(elem,i,k)*S(elem,k,j) +
                                C2*(R(elem,i,k)*S(elem,k,j) + R(elem,j,k)*S(elem,k,i)) +
                                C3 *R(elem,i,k)*R(elem,j,k);
            if (i==j)
              {
                for (int k=0; k<Objet_U::dimension; k++)
                  for (int l=0; l<Objet_U::dimension; l++)
                    ReNL(elem,i,j) += 1./3.*(C1*S(elem,k,l)*S(elem,k,l)+C3*R(elem,k,l)*R(elem,k,l));
              }
          }
    }

  return ReNL;
}


DoubleTab Modele_Lam_Bremhorst_VEF::calcul_tenseur_Re_shih(const Discretisation_base& dis,
                                                           const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                           const DoubleTab& G,
                                                           const Champ_base& K_Eps) const
{

  const Domaine_VEF&       domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  int nelem = G.dimension(0);

  DoubleTab S, R;
  init_tenseur_elem(S,domaine_VEF,2);
  init_tenseur_elem(R,domaine_VEF,2);

  for (int elem=0; elem<nelem; elem++)
    for (int i=0; i<Objet_U::dimension; i++)
      for (int j=0; j<Objet_U::dimension; j++)
        {
          S(elem,i,j) = 0.5*(G(elem,i,j) + G(elem,j,i));
          R(elem,i,j) = 0.5*(G(elem,i,j) - G(elem,j,i));
          if (i==j)
            for (int k=0; k<Objet_U::dimension; k++)
              {
                S(elem,i,j) -= 1/3 * G(elem,k,k);
              }
        }

  DoubleTab Sn = calcul_norme_elem(domaine_VEF,S);

  DoubleTab Rn = calcul_norme_elem(domaine_VEF,R);

  DoubleTab ReNL = calcul_tenseur_Re_elem_shih(dis,domaine_dis,G,S,R,Sn,Rn,K_Eps);

  return ReNL;
}

DoubleTab Modele_Lam_Bremhorst_VEF::calcul_tenseur_Re_elem_shih(const Discretisation_base& dis,
                                                                const Domaine_dis_base& domaine_dis,
                                                                const DoubleTab& G, const DoubleTab& S,
                                                                const DoubleTab& R, const DoubleTab& Sn,const DoubleTab& Rn,
                                                                const Champ_base& K_Eps) const
{
  int nelem = G.dimension(0);

  Champ_Fonc K_Eps_elem;
  Noms noms(2), unites(2);
  noms[0]="K";
  noms[1]="eps";
  unites[0]="m2/s2";
  unites[1]="m2/s3";
  const Motcle wtype= "champ_elem";
  dis.discretiser_champ(wtype,domaine_dis,multi_scalaire,noms,unites,2,0.,K_Eps_elem);
  K_Eps_elem->affecter(K_Eps);
  DoubleTab tab_K_Eps  = K_Eps_elem->valeurs();

  DoubleTab ReNL;
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  init_tenseur_elem(ReNL,domaine_VEF,2);

  ReNL = 0.;
  for (int elem=0; elem<nelem; elem++)
    {
      double kseps;
      if (K_Eps.valeurs()(elem,1) <= BR_EPS)
        {
          kseps = K_Eps.valeurs()(elem,0)/BR_EPS;
        }
      else
        {
          kseps = tab_K_Eps(elem,0)/tab_K_Eps(elem,1);
        }
      double Cmu = (1.)/(6.5+sqrt(pow(Sn(elem),2)+pow(Rn(elem),2))*1.8*kseps);
      double C2 = sqrt(1-pow(3*Cmu*Sn(elem)*kseps,2))/(1+6.*Sn(elem)*kseps*Rn(elem)*kseps);
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            ReNL(elem,i,j) = 2*S(elem,i,j);
            for (int k=0; k<Objet_U::dimension; k++)
              ReNL(elem,i,j) -= 2*C2*kseps/Cmu*(-S(elem,i,k)*R(elem,k,j)+S(elem,k,j)*R(elem,i,k));
          }
    }

  return ReNL;
}

// Calcul de la norme d'un tenseur aux elements
DoubleTab Modele_Lam_Bremhorst_VEF::calcul_norme_elem(const Domaine_VEF& domaine_VEF,const DoubleTab Tenseur) const
{
  int dim_tens = 0;
  DoubleTab Tnorme;
  init_tenseur_elem(Tnorme,domaine_VEF,dim_tens);
  int nb_elems = Tnorme.dimension(0);
  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          Tnorme(num_elem) += Tenseur(num_elem,i,j)*Tenseur(num_elem,i,j);
      Tnorme(num_elem) = sqrt(Tnorme(num_elem));
    }
  return Tnorme;
}

void Modele_Lam_Bremhorst_VEF::lire_distance_paroi( )
{

  // PQ : 25/02/04 recuperation de la distance a la paroi dans Wall_length.xyz

  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  DoubleTab& wall_length = BR_wall_length_->valeurs();
  wall_length=-1.;
//  return;

  LecFicDiffuse fic;
  fic.set_bin(1);
  if(!fic.ouvrir(nom_fic))
    {
      Cerr << " File " <<nom_fic<< " doesn't exist. To generate it, please, refer to html.doc (Distance_paroi) " << finl;
      exit();
    }

  Noms nom_paroi;

  fic >> nom_paroi;

  if(je_suis_maitre())
    {

      Cerr << "Recall : " <<nom_fic<< " obtained from boundaries nammed : "<< finl;
      for (int b=0; b<nom_paroi.size(); b++)
        {
          Cerr << nom_paroi[b]<< finl;
          //test pour s'assurer de la coherence de Wall_length.xyz avec le jeu de donnees :
          domaine_VEF.rang_frontiere(nom_paroi[b]);
        }
    }
  EcritureLectureSpecial::lecture_special(domaine_VEF, fic, wall_length);

}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_Cmu(DoubleTab& Cmu,
                                                const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN) const
{

  const Domaine_VEF&       domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());

  DoubleTab gradient_elem;
  init_tenseur_elem(gradient_elem,domaine_VEF,2);
  Champ_P1NC::calcul_gradient(vitesse,gradient_elem,domaine_Cl_VEF);

  DoubleTab S_elem;
  init_tenseur_elem(S_elem,domaine_VEF,2);
  int nelem = S_elem.dimension(0);
  for (int elem=0; elem<nelem; elem++)
    for (int i=0; i<dimension; i++)
      for (int j=0; j<dimension; j++)
        S_elem(elem,i,j) =  gradient_elem(elem,i,j) + gradient_elem(elem,j,i);
  S_elem.echange_espace_virtuel();

  DoubleTab S_face;
  init_tenseur_face(S_face,domaine_VEF,2);

  calcul_tenseur_face(S_face,S_elem,domaine_VEF,domaine_Cl_VEF);

  int nfaces = S_face.dimension(0);
  DoubleTab Snorme_face;
  domaine_VEF.creer_tableau_faces(Snorme_face);

  for (int face=0; face<nfaces; face++)
    {
      double somme = 0.;
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          somme += S_face(face,i,j)*S_face(face,i,j);
      Snorme_face(face) = sqrt(somme);
    }

  for (int face=0; face<nfaces; face++)
    {
      // Definition d'un Cmu minimum base sur EPS_MIN = 1e-10
      if (K_Eps(face,1) <= EPS_MIN)
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K_Eps(face,0)/BR_EPS);
      else
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K_Eps(face,0)/K_Eps(face,1));
    }

  return Cmu;
}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_Cmu_Paroi(DoubleTab& Cmu,
                                                      const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                      const DoubleTab& visco, const DoubleTab& visco_turb,
                                                      const DoubleTab& loi_paroi,const int idt,
                                                      const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN) const
{

  const Domaine_VEF&       domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());

  DoubleTab gradient_elem;
  init_tenseur_elem(gradient_elem,domaine_VEF,2);
  Champ_P1NC::calcul_gradient(vitesse,gradient_elem,domaine_Cl_VEF);

  if (idt>0)
    Champ_P1NC::calcul_duidxj_paroi(gradient_elem,visco,visco_turb,loi_paroi,domaine_Cl_VEF);
  DoubleTab S_elem;
  init_tenseur_elem(S_elem,domaine_VEF,2);
  int nelem = S_elem.dimension(0);
  for (int elem=0; elem<nelem; elem++)
    for (int i=0; i<dimension; i++)
      for (int j=0; j<dimension; j++)
        S_elem(elem,i,j) =  gradient_elem(elem,i,j) + gradient_elem(elem,j,i);
  S_elem.echange_espace_virtuel();

  DoubleTab S_face;
  init_tenseur_face(S_face,domaine_VEF,2);

  calcul_tenseur_face(S_face,S_elem,domaine_VEF,domaine_Cl_VEF);

  int nfaces = S_face.dimension(0);
  DoubleTab Snorme_face;
  domaine_VEF.creer_tableau_faces(Snorme_face);

  for (int face=0; face<nfaces; face++)
    {
      double somme = 0.;
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          somme += S_face(face,i,j)*S_face(face,i,j);
      Snorme_face(face) = sqrt(somme);
    }

  for (int face=0; face<nfaces; face++)
    {
      // Definition d'un Cmu minimum base sur EPS_MIN = 1e-10
      if (K_Eps(face,1) <= EPS_MIN)
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K_Eps(face,0)/BR_EPS);
      else
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K_Eps(face,0)/K_Eps(face,1));
    }

  // TMA nettoyage de restes de debug
  // Cerr<<Cmu.mp_min_vect()<<" Cmuuuuuuuuuuuuuuuuuuuu " <<Cmu.mp_max_vect()<<finl;
  return Cmu;
}


DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_Cmu_BiK(DoubleTab& Cmu,
                                                    const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                    const DoubleTab& vitesse, const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) const
{

  const Domaine_VEF&       domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());

  DoubleTab gradient_elem;
  init_tenseur_elem(gradient_elem,domaine_VEF,2);
  Champ_P1NC::calcul_gradient(vitesse,gradient_elem,domaine_Cl_VEF);

  DoubleTab S_elem;
  init_tenseur_elem(S_elem,domaine_VEF,2);
  int nelem = S_elem.dimension(0);
  for (int elem=0; elem<nelem; elem++)
    for (int i=0; i<dimension; i++)
      for (int j=0; j<dimension; j++)
        S_elem(elem,i,j) =  gradient_elem(elem,i,j) + gradient_elem(elem,j,i);
  S_elem.echange_espace_virtuel();

  DoubleTab S_face;
  init_tenseur_face(S_face,domaine_VEF,2);

  calcul_tenseur_face(S_face,S_elem,domaine_VEF,domaine_Cl_VEF);

  int nfaces = S_face.dimension(0);
  DoubleTab Snorme_face;
  domaine_VEF.creer_tableau_faces(Snorme_face);

  for (int face=0; face<nfaces; face++)
    {
      double somme = 0.;
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          somme += S_face(face,i,j)*S_face(face,i,j);
      Snorme_face(face) = sqrt(somme);
    }

  for (int face=0; face<nfaces; face++)
    {
      // Definition d'un Cmu minimum base sur EPS_MIN = 1e-10
      if (Eps(face) <= EPS_MIN)
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K(face)/BR_EPS);
      else
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K(face)/Eps(face));
    }

  return Cmu;
}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_Cmu_Paroi_BiK(DoubleTab& Cmu,
                                                          const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                          const DoubleTab& visco, const DoubleTab& visco_turb,
                                                          const DoubleTab& loi_paroi,const int idt,
                                                          const DoubleTab& vitesse, const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN) const
{

  const Domaine_VEF&       domaine_VEF = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());

  DoubleTab gradient_elem;
  init_tenseur_elem(gradient_elem,domaine_VEF,2);
  Champ_P1NC::calcul_gradient(vitesse,gradient_elem,domaine_Cl_VEF);

  if (idt>0)
    Champ_P1NC::calcul_duidxj_paroi(gradient_elem,visco,visco_turb,loi_paroi,domaine_Cl_VEF);
  DoubleTab S_elem;
  init_tenseur_elem(S_elem,domaine_VEF,2);
  int nelem = S_elem.dimension(0);
  for (int elem=0; elem<nelem; elem++)
    for (int i=0; i<dimension; i++)
      for (int j=0; j<dimension; j++)
        S_elem(elem,i,j) =  gradient_elem(elem,i,j) + gradient_elem(elem,j,i);
  S_elem.echange_espace_virtuel();

  DoubleTab S_face;
  init_tenseur_face(S_face,domaine_VEF,2);

  calcul_tenseur_face(S_face,S_elem,domaine_VEF,domaine_Cl_VEF);

  int nfaces = S_face.dimension(0);
  DoubleTab Snorme_face;
  domaine_VEF.creer_tableau_faces(Snorme_face);

  for (int face=0; face<nfaces; face++)
    {
      double somme = 0.;
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          somme += S_face(face,i,j)*S_face(face,i,j);
      Snorme_face(face) = sqrt(somme);
    }

  for (int face=0; face<nfaces; face++)
    {
      // Definition d'un Cmu minimum base sur EPS_MIN = 1e-10
      if (Eps(face) <= EPS_MIN)
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K(face)/BR_EPS);
      else
        Cmu[face] = (2./3.)/(A1+sqrt(0.5)*Snorme_face[face]*K(face)/Eps(face));
    }

  // TMA nettoyage de restes de debug
  // Cerr<<Cmu.mp_min_vect()<<" Cmuuuuuuuuuuuuuuuuuuuu " <<Cmu.mp_max_vect()<<finl;
  return Cmu;
}


DoubleTab&  Modele_Lam_Bremhorst_VEF::Calcul_Fmu_BiK( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis& domaine_Cl_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());
  const DoubleTab& wall_length = BR_wall_length_->valeurs();
  DoubleTab wall_length_face(0);
  le_dom.creer_tableau_faces(wall_length_face);
  DoubleTab Pderive(0);
  le_dom.creer_tableau_faces(Pderive);
  int nb_faces = le_dom.nb_faces();
  const Conds_lim& les_cl = domaine_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const IntTab& face_voisins = le_dom.face_voisins();
  int num_face;
  double Rey,Re;
  /*
  	for (num_face=0; num_face <nb_faces; num_face ++ )
  		F1[num_face] = 1.;
  	return F1;
  */
  // Calcul de la distance a la paroi aux faces
  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              int elem2 = face_voisins(num_face,1);
              wall_length_face(num_face) = 0.5*wall_length(elem1) + 0.5*wall_length(elem2);
            }
        }
      /* else if (sub_type(Pa,la_cl.valeur()))
      {
        for (int num_face=ndeb; num_face<nfin; num_face++)
        {
      	  int elem1 = face_voisins(num_face,0);
      	  int elem2 = face_voisins(num_face,1);
      			  wall_length_face(num_face) = 0.;
        }
      }*/
      else
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              wall_length_face(num_face) = wall_length(elem1);
            }
        }
    }
  int n0 = le_dom.premiere_face_int();
  for (num_face=n0; num_face<nb_faces; num_face++)
    {
      int elem0 = le_dom.face_voisins(num_face,0);
      int elem1 = le_dom.face_voisins(num_face,1);
      wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
    }
  // Calcul de la distance a la paroi aux faces
  /*    for (num_face=0; num_face< le_dom.premiere_face_int(); num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  if (elem0 != -1)
    		  wall_length_face(num_face) = wall_length(elem0);
    	  else
    	  {
    		  elem0 = le_dom.face_voisins(num_face,1);
    		  wall_length_face(num_face) = wall_length(elem0);
    	  }
      }

      for (; num_face<nb_faces; num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  int elem1 = le_dom.face_voisins(num_face,1);
    	  wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
      }
  */
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      if (visco>BR_EPS && eps_Bas_Re(num_face)>BR_EPS && K_Bas_Re(num_face))
        {
          Re = (K_Bas_Re(num_face)*K_Bas_Re(num_face))/(visco*eps_Bas_Re(num_face));
          Rey = wall_length_face(num_face)*sqrt(K_Bas_Re(num_face))/visco;
          Fmu(num_face) = (1. - exp(-0.0165*Rey-BR_EPS))*(1. - exp(-0.0165*Rey-BR_EPS))*(1.+20.5/(Re+BR_EPS));

        }
      else
        {
          //Cerr <<  " visco " << visco << " Eps " << K_eps_Bas_Re(num_face,1) << " K " << K_eps_Bas_Re(num_face,0) <<finl;
          Fmu(num_face) = 1.;
        }
      if (Fmu(num_face) > 1. + 1.e-8)
        {
          //Cerr <<  " Fmu=" << Fmu(num_face) << " Eps=" << K_eps_Bas_Re(num_face,1) << " K=" << K_eps_Bas_Re(num_face,0) << " Re=" << Re << " Rey=" << Rey <<finl;
          Fmu(num_face) = 1.;
          //Cerr <<  " On force Fmu a 1 " << finl;
          //exit();
        }
    }
  Cerr<< Fmu.mp_min_vect()<< " Fmu "<<Fmu.mp_max_vect()<<finl;

  return Fmu;
}


DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_F2_BiK( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis_base& domaine_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  int num_face;
  double Re;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (!is_visco_const)
        {
          int elem0 = le_dom.face_voisins(num_face,0);
          int elem1 = le_dom.face_voisins(num_face,1);
          if (elem1!=-1)
            {
              visco = tab_visco(elem0)*le_dom.volumes(elem0)+tab_visco(elem1)*le_dom.volumes(elem1);
              visco /= le_dom.volumes(elem0) + le_dom.volumes(elem1);
            }
          else
            visco =  tab_visco(elem0);

        }
      if (visco>BR_EPS && eps_Bas_Re(num_face)>BR_EPS)
        {
          Re = (K_Bas_Re(num_face)*K_Bas_Re(num_face))/(visco*eps_Bas_Re(num_face));
          F2[num_face] = 1. - (exp(-Re*Re));
        }
      else
        F2[num_face] = 1.;

    }
  Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_F1_BiK( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_Bas_Re, const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const
{

  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis.valeur());
  const DoubleTab& wall_length = BR_wall_length_->valeurs();
  DoubleTab wall_length_face(0);
  le_dom.creer_tableau_faces(wall_length_face);
  DoubleTab Pderive(0);
  le_dom.creer_tableau_faces(Pderive);
  DoubleTab Fmu_loc(0);
  le_dom.creer_tableau_faces(Fmu_loc);
  int nb_faces = le_dom.nb_faces();
  const Conds_lim& les_cl = domaine_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const IntTab& face_voisins = le_dom.face_voisins();
  int num_face;
  double Rey,Re;
  /*
  	for (num_face=0; num_face <nb_faces; num_face ++ )
  		F1[num_face] = 1.;
  	return F1;
  */
  // Calcul de la distance a la paroi aux faces
  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              int elem2 = face_voisins(num_face,1);
              wall_length_face(num_face) = 0.5*wall_length(elem1) + 0.5*wall_length(elem2);
            }
        }
      /* else if (sub_type(Pa,la_cl.valeur()))
      {
        for (int num_face=ndeb; num_face<nfin; num_face++)
        {
      	  int elem1 = face_voisins(num_face,0);
      	  int elem2 = face_voisins(num_face,1);
      			  wall_length_face(num_face) = 0.;
        }
      }*/
      else
        {
          for (num_face=ndeb; num_face<nfin; num_face++)
            {
              int elem1 = face_voisins(num_face,0);
              wall_length_face(num_face) = wall_length(elem1);
            }
        }
    }
  int n0 = le_dom.premiere_face_int();
  for (num_face=n0; num_face<nb_faces; num_face++)
    {
      int elem0 = le_dom.face_voisins(num_face,0);
      int elem1 = le_dom.face_voisins(num_face,1);
      wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
    }
  // Calcul de la distance a la paroi aux faces
  /*    for (num_face=0; num_face< le_dom.premiere_face_int(); num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  if (elem0 != -1)
    		  wall_length_face(num_face) = wall_length(elem0);
    	  else
    	  {
    		  elem0 = le_dom.face_voisins(num_face,1);
    		  wall_length_face(num_face) = wall_length(elem0);
    	  }
      }

      for (; num_face<nb_faces; num_face++)
      {
    	  int elem0 = le_dom.face_voisins(num_face,0);
    	  int elem1 = le_dom.face_voisins(num_face,1);
    	  wall_length_face(num_face) = 0.5*wall_length(elem0)+0.5*wall_length(elem1);
      }
  */
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      if (visco>BR_EPS && eps_Bas_Re(num_face)>BR_EPS && K_Bas_Re(num_face))
        {
          Re = (K_Bas_Re(num_face)*K_Bas_Re(num_face))/(visco*eps_Bas_Re(num_face));
          Rey = wall_length_face(num_face)*sqrt(K_Bas_Re(num_face))/visco;
          Fmu_loc(num_face) = (1. - exp(-0.0165*Rey-BR_EPS))*(1. - exp(-0.0165*Rey-BR_EPS))*(1.+20.5/(Re+BR_EPS));
        }
      else
        Fmu_loc(num_face) = 1.;

      if (Fmu_loc(num_face) > 1. + 1.e-8)
        {
          Fmu_loc(num_face) = 1.;
          //Cerr <<  " On force Fmu_loc a 1 " << finl;
          //exit();
        }
    }
  Cerr<< Fmu_loc.mp_min_vect()<< " Fmu_loc "<<Fmu_loc.mp_max_vect()<<finl;

  // Calcul de F1
  for (num_face=0; num_face <nb_faces; num_face ++ )
    {
      if (Fmu_loc(num_face) > BR_EPS && (visco>BR_EPS && eps_Bas_Re(num_face)>BR_EPS && K_Bas_Re(num_face)))
        F1[num_face] = 1.+ (0.05/(Fmu_loc(num_face)+BR_EPS))*(0.05/(Fmu_loc(num_face)+BR_EPS))*(0.05/(Fmu_loc(num_face)+BR_EPS));
      else
        F1[num_face] = 1.;
    }
  Cerr<< F1.mp_min_vect()<< " F1 "<<F1.mp_max_vect()<<finl;
  return F1;
}



DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_E_BiK(DoubleTab& E,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  E = 0;
  return E;
}

DoubleTab& Modele_Lam_Bremhorst_VEF::Calcul_D_BiK(DoubleTab& D,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis,
                                                  const DoubleTab& vitesse,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  D = 0;
  return D;
}









