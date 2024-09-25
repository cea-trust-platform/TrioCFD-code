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
// File:        Implicite_steady.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Schema_Euler_Implicite_Stationnaire/src/New
// Version:     /main/29
//
//////////////////////////////////////////////////////////////////////////////

#include <Implicit_steady.h>
#include <Schema_Euler_Implicite_Stationnaire.h>
#include <Domaine_VF.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Domaine_VDF.h>
#include <Navier_Stokes_std.h>
#include <EChaine.h>
#include <Debog.h>
#include <Matrice_Bloc.h>
#include <Assembleur_base.h>
#include <Discretisation_base.h>
#include <Statistiques.h>
#include <Schema_Temps_base.h>
#include <TRUSTTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Dirichlet.h>
#include <Probleme_base.h>
#include <assert.h>

Implemente_instanciable(Implicit_steady,"Implicit_steady",Simpler);
// XD implicit_steady implicite  implicit_steady 1  this is the implicit solver using a dual time step. Remark: this solver can be used only with the Implicit_Euler_Steady_Scheme time scheme.

Sortie& Implicit_steady::printOn(Sortie& os ) const
{
  return Simpler::printOn(os);
}

Entree& Implicit_steady::readOn(Entree& is )
{
  Simpler::readOn(is);
  return is;
}




void test_impose_bound_cond(Equation_base& eqn,DoubleTab& current2,const char * msg,int flag)
{
  return;
  DoubleTab& present = eqn.inconnue().futur();
  DoubleTab sauv(present);
  const Schema_Temps_base& sch = eqn.probleme().schema_temps();
  eqn.domaine_Cl_dis().imposer_cond_lim(eqn.inconnue(),sch.temps_courant()+sch.pas_de_temps());
  present -= sauv;
  // BM, je remplace max_abs par mp_pax_abs: du coup la methode doit etre appelee simultanement par tous les procs.
  double ecart_max=mp_max_abs_vect(present);
  Cout<<msg <<" "<<ecart_max<<finl;

  if ((ecart_max>1e-10))
    abort();
  present = sauv;
}


//This function slightly modifies the "Piso :: iterate_NS with avancement_crank_=0" function in order to take into account a variable time step
//Use of a variable time step with the IMPLICITE algorithm
//replaces the global time step with a local time step
void Implicit_steady::iterer_NS(Equation_base& eqn,DoubleTab& current,DoubleTab& pression, double dt,
                                Matrice_Morse& matrice,double seuil_resol,DoubleTrav& secmem,int nb_iter,int& converge, int& ok)
{
  Schema_Temps_base& sch = eqn.probleme().schema_temps();
  Parametre_implicite& param_eqn = get_and_set_parametre_implicite(eqn);
  SolveurSys& le_solveur_ = param_eqn.solveur();
  converge = 1;
  Navier_Stokes_std& eqnNS = ref_cast(Navier_Stokes_std,eqn);
  if(!sub_type(Schema_Euler_Implicite_Stationnaire,sch))
    {
      Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
      Cerr << "Implicit_steady solveur can be used only with the Implicit_Euler_Steady_Scheme or Schema_Euler_Implicite_Stationnaire scheme!" << finl;
      Cerr << "Please, contact TRUST support." << finl;
      exit();
    }


  //DoubleTab& passe = eqnNS.inconnue().passe(2);
  DoubleTab present(current);
  DoubleTrav gradP(current);
  DoubleTrav correction_en_pression(pression);
  DoubleTrav resu(current);

  Operateur_Grad& gradient = eqnNS.operateur_gradient();
  Operateur_Div& divergence = eqnNS.operateur_divergence();

  //Construction de matrice et resu
  //matrice = A[Un] = M/dt_locaux + CONV +DIFF
  //resu = M/dt_locaux*Un  -BtPn + Sv
  gradient.calculer(pression,gradP);
  resu -= gradP;

  eqnNS.assembler_avec_inertie(matrice,current,resu);
  le_solveur_->reinit();
  Debog::verifier("Implicit_steady::iterer_NS resu apres assembler_avec_inertie",resu);

  //Definition de matrice_en_pression_2
  Matrice& matrice_en_pression_2 = eqnNS.matrice_pression();
  SolveurSys& solveur_pression_ = eqnNS.solveur_pression();

  //Etape predicteur
  Cout << "Compute U* :" << finl;
  //Resolution du systeme A[Un]U* = -BtPn + Sv + Ss
  //current = U*
  le_solveur_.resoudre_systeme(matrice,resu,current);
  test_impose_bound_cond(eqn,current,"apres resolution ",0);
  Debog::verifier("Implicit_steady::iterer_NS current apres CL",current);
  /*
  #ifndef NDEBUG
    // Check periodic BC on U*
    test_periodic_solution(eqnNS, current);
  #endif
  */
  current.echange_espace_virtuel();

  //Calcul de secmem = BU*
  if (eqn.probleme().is_dilatable())
    {
      Cerr<<"Steady option is not compatible with the quasi/weakly compressible models !"<<finl;
      Cerr << "Please, contact TRUST support." << finl;
      Process::exit();
    }
  else
    {
      divergence.calculer(current,secmem);
    }
  secmem *= -1;
  secmem.echange_espace_virtuel();

  //Etape de correction
  Cout << "Solving mass equation :" << finl;
  //Description du cas implicite
  //Resolution du systeme (B*dt_locaux*M-1Bt)P' = Bu*
  //correction_en_pression = P'

  DoubleVect dt_locaux = sch.pas_de_temps_locaux();

  //Calcul du M/dt_locaux
  DoubleVect m_dt(dt_locaux);
  const Nom discr=eqnNS.discretisation().que_suis_je();
  if (discr=="VDF")
    calcul_mat_masse_diviser_par_dt_vdf(eqnNS, m_dt, dt_locaux);
  else
    calcul_mat_masse_diviser_par_dt_vef(eqnNS, m_dt, dt_locaux);

  //DoubleVect  m_dt(dt_locaux);
  //eqnNS.solv_masse().get_masse_divide_by_local_dt(m_dt, dt_locaux, 0);


  //Construction de matrice_en_pression_2 = B*dt_locaux*M-1Bt
  eqnNS.assembleur_pression()->assembler_mat(matrice_en_pression_2,m_dt,1,1);
  solveur_pression_->reinit();
  //Resolution du systeme (B*dt_locaux*M-1Bt)P' = Bu*
  solveur_pression_.resoudre_systeme(matrice_en_pression_2.valeur(),
                                     secmem,correction_en_pression);
  correction_en_pression.echange_espace_virtuel();
  //Calcul de M^-1BtP'=gradP
  gradient->multvect(correction_en_pression,gradP);
  eqn.solv_masse().appliquer(gradP);
  //Calcul de Un+1 = U* -dt_locaux*gradP -M*dt_locaux*deltaU/dt
  //dt = pas de temps global
  //deltaU = Un+1 -Un
  int size=gradP.size();

  //dt_locaux_taille_vitesse pas de temps local: Vect de la meme taille que le vecteur vitesse
  //construit a partir de dt_face: pas de temps par face
  //dt_locaux_taille_vitesse par face a la meme valeurs dans les directions : x,y,z
  DoubleVect  dt_locaux_taille_vitesse(current);
  if (discr=="VDF")
    {
      dt_locaux_taille_vitesse = dt_locaux;
    }
  else
    {
      int nb_dim = eqnNS.inconnue().dimension;
      assert(dt_locaux.size()*nb_dim == dt_locaux_taille_vitesse.size());
      int j = 0;
      for(int i=0; i<dt_locaux_taille_vitesse.size(); i++)
        {
          if(i != 0 && i%nb_dim == 0) j++;
          dt_locaux_taille_vitesse[i] = dt_locaux[j];
        }
      dt_locaux_taille_vitesse.echange_espace_virtuel();
    }
  DoubleVect  dt_locaux_masse (dt_locaux_taille_vitesse);
  //Calcul du M*dt_locaux
  eqnNS.solv_masse().get_masse_dt_local(dt_locaux_masse, dt_locaux_taille_vitesse, 0);


  for(int i=0; i<size; i++)
    {
      present.addr()[i] *=dt_locaux_masse[i]/dt; //M*dt_locaux*Un/dt
      gradP.addr()[i] *=dt_locaux_taille_vitesse[i];  //dt_locaux*gradP
    }


  present.echange_espace_virtuel();
  gradP.echange_espace_virtuel();
  current -= gradP;
  // ajustement pour prendre en compte le terme M*dt_locaux*deltaU/dt
  current += present;
  for(int i=0; i<size; i++)
    {
      current.addr()[i] *= dt/(dt+dt_locaux_masse[i]);
    }
  test_impose_bound_cond(eqn,current,"apres resolution ",0);
  //current=Un+1
  current.echange_espace_virtuel();
  divergence.calculer(current,secmem);


  //Calcul de Pn+1 = Pn + P'
  pression += correction_en_pression;
  Debog::verifier("Implicit_steady::iterer_NS pression fin", pression);
  eqnNS.assembleur_pression()->modifier_solution(pression);
  pression.echange_espace_virtuel();
  return;
}
//Calcule M/dt_locaux
void Implicit_steady::calcul_mat_masse_diviser_par_dt_vef(Navier_Stokes_std& eqnNS, DoubleVect& m_dt, DoubleVect& dt_locaux)
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,eqnNS.domaine_dis());
  const DoubleVect& volumes_entrelaces_ref=le_dom.volumes_entrelaces();
  const Domaine_Cl_VEF& le_dom_cl = ref_cast(Domaine_Cl_VEF,eqnNS.domaine_Cl_dis());
  const DoubleVect& volumes_entrelaces_cl=le_dom_cl.volumes_entrelaces_Cl();
  DoubleVect volumes_entrelaces(volumes_entrelaces_ref);
  int size_cl=volumes_entrelaces_cl.size();

  for (int f=0; f<size_cl; f++)
    if (volumes_entrelaces_cl(f)!=0)
      {
        volumes_entrelaces(f)=volumes_entrelaces_cl(f);
      }
  volumes_entrelaces.echange_espace_virtuel();

  int size=volumes_entrelaces.size_totale();
  // Si rho n'est pas constant
  const DoubleVect& masse_volumique = eqnNS.fluide().masse_volumique().valeurs();
  if(masse_volumique.size_totale()==size)
    {
      for (int face=0; face<size; face++)
        {
          m_dt(face) = volumes_entrelaces(face)*masse_volumique(face) / dt_locaux(face);
        }
    }
  else
    {
      for (int face=0; face<size; face++)
        {
          m_dt(face) = volumes_entrelaces(face)/ dt_locaux(face);
        }
    }
  m_dt.echange_espace_virtuel();

}

void Implicit_steady::calcul_mat_masse_diviser_par_dt_vdf(Navier_Stokes_std& eqnNS, DoubleVect& m_dt, DoubleVect& dt_locaux)
{
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,eqnNS.domaine_dis());
  const DoubleVect& volumes_entrelaces=le_dom.volumes_entrelaces();
  int size=volumes_entrelaces.size_totale();
  // Si rho n'est pas constant
  const DoubleVect& masse_volumique = eqnNS.fluide().masse_volumique().valeurs();
  if(masse_volumique.size_totale()==size)
    {
      for (int face=0; face<size; face++)
        {
          m_dt(face) = volumes_entrelaces(face)*masse_volumique(face) / dt_locaux(face);
        }
    }
  else
    {
      for (int face=0; face<size; face++)
        {
          m_dt(face) = volumes_entrelaces(face)/ dt_locaux(face);
        }
    }
  m_dt.echange_espace_virtuel();

}

/*//Fonction utile uniquement pour le Debug
void Implicit_steady::test_periodic_solution(Navier_Stokes_std& eqnNS, DoubleTab& current) const
{

  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,eqnNS.domaine_dis());
  const Domaine_Cl_VEF& le_dom_cl = ref_cast(Domaine_Cl_VEF,eqnNS.domaine_Cl_dis());
  int nb_comp=current.dimension(1);
  for (int n_bord=0; n_bord<le_dom.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = le_dom_cl.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int nb_faces_bord=le_bord.nb_faces();
          for (int ind_face=0; ind_face<nb_faces_bord; ind_face++)
            {
              int ind_face_associee = la_cl_perio.face_associee(ind_face);
              int face = le_bord.num_face(ind_face);
              int face_associee = le_bord.num_face(ind_face_associee);
              for (int comp=0; comp<nb_comp; comp++)
                {
                  if (!est_egal(current(face, comp),current(face_associee, comp),1.e-4))
                    {
                      Cerr << "Periodic boundary condition is not correct in Implicite_steady::test_periodic_solution" << finl;
                      Cerr << "vit1("<<face<< ","<<comp<<")=" << current(face, comp) << finl;
                      Cerr << "vit2("<<face_associee<<","<<comp<<")=" << current(face_associee,comp) << finl;
                      Cerr << "Delta=" << current(face,comp)-current(face_associee,comp) << finl;
                      exit();
                    }
                }
            }
        }
    }

}*/
