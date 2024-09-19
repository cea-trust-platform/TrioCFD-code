/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Paroi_std_hyd_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_std_hyd_VEF.h>
#include <Champ_Q1NC.h>
#include <Champ_P1NC.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Periodique.h>
#include <Static_Int_Lists.h>
#include <Debog.h>
#include <TRUSTList.h>
#include <EcrFicPartage.h>
#include <Modele_turbulence_hyd_Longueur_Melange_base.h>
#include <Neumann_sortie_libre.h>
#include <Modele_turbulence_hyd_combinaison.h>
#include <Param.h>
#include <Paroi_rugueuse.h>
#include <SFichier.h>
#include <Paroi_decalee_Robin.h>
#include <Schema_Temps_base.h>
#include <communications.h>

// methode_calcul_face_keps_impose_ value:
// 0: std avant
// 1: toutes les faces accrochees
// 2: comme avant mais toutes les faces des elts accroches.
// 3: comme avant sans test si bord...
// 4: que les faces des elements Dirichlet accroches
// 5: 4 avec test de distance pour //
// -1: 1 en 2D, 5 en 3D

Implemente_instanciable_sans_constructeur(Paroi_std_hyd_VEF,"loi_standard_hydr_VEF",Paroi_hyd_base_VEF);
Implemente_instanciable(Loi_expert_hydr_VEF,"Loi_expert_hydr_VEF",Paroi_std_hyd_VEF);

Paroi_std_hyd_VEF::Paroi_std_hyd_VEF(): methode_calcul_face_keps_impose_(-1),is_u_star_impose_ (0) {}

Sortie& Paroi_std_hyd_VEF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Sortie& Loi_expert_hydr_VEF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}
//// readOn
//

Entree& Paroi_std_hyd_VEF::readOn(Entree& s)
{
  return Paroi_hyd_base_VEF::readOn(s);
}

Entree& Loi_expert_hydr_VEF::readOn(Entree& s)
{

  Param param(que_suis_je());
  Paroi_std_hyd_VEF::set_param(param);

  param.ajouter("u_star_impose",&u_star_impose_);
  param.ajouter("methode_calcul_face_keps_impose",&methode_calcul_face_keps_impose_);

  param.dictionnaire("toutes_les_faces_accrochees",1);
  param.dictionnaire("que_les_faces_des_elts_dirichlet",4);
  param.dictionnaire("que_les_faces_des_elts_dirichlet_et_test_distance",5);
  param.lire_avec_accolades_depuis(s);

  Nom mot_test;
  mot_test = "u_star_impose";

  if (param.get_list_mots_lus().rang(Motcle(mot_test))!=-1)
    is_u_star_impose_ = 1;

  return s ;
}

void Paroi_std_hyd_VEF::set_param(Param& param)
{
  Paroi_hyd_base_VEF::set_param(param);
  Paroi_log_QDM::set_param(param);
}

/////////////////////////////////////////////////////////////////////
//
//  Implementation des fonctions de la classe Paroi_std_hyd_VEF
//
/////////////////////////////////////////////////////////////////////
bool compute_u_plus(const int ind_face,
                    const double u_plus_d_plus,
                    const double erugu,
                    const double Kappa,
                    double& r, int& iter, double& u_plus)
{
  // PQ : 05/04/07 : formulation continue de la loi de paroi
  // construite d'apres la loi de Reichardt
  if(u_plus_d_plus < 1.)
    {
      u_plus = sqrt(u_plus_d_plus); // pour eviter de faire tourner la procedure iterative
      return true;
    }
  u_plus = u_plus_d_plus/100.;
  r = 1.;
  const double seuil = 0.001;
  iter = 0;
  const int itmax  = 25;
  const double A = (1/Kappa)*log(erugu/Kappa) ; //  (=7.44, contre 7.8 dans la loi d'origine)
  //  permettant d'avoir en l'infini la loi de
  //  Reichardt se calant sur : u+ = (1/Kappa).ln(Erugu.y+)

  // Fixed point method for Reichardt equation
  while((iter++ < itmax) && (r > seuil))
    {
      const double d_plus  = u_plus_d_plus/u_plus;
      const double u_plus2 = (1./Kappa)*log(1. + Kappa*d_plus)
                             + A*(1. - exp(-d_plus/11.) - exp(-d_plus/3.)*d_plus/11.);
      u_plus  = 0.5*(u_plus + u_plus2);
      r       = std::fabs(u_plus - u_plus2)/u_plus;
    }
  return (iter < itmax);
}

KOKKOS_FUNCTION
double calculer_u_plus_kokkos(const int ind_face, const double u_plus_d_plus, const double erugu, const double Kappa, DoubleArrView seuil_LP, IntArrView iterations_LP)
{
  double r, u_plus;
  int iter;
  if (!compute_u_plus(ind_face, u_plus_d_plus, erugu, Kappa, r, iter, u_plus))
#ifdef KOKKOS
    Kokkos::abort("The iterative process of u* did not converge!");
#else
    Process::exit("The iterative process of u* did not converge!");
#endif
  seuil_LP(ind_face) = r;
  iterations_LP(ind_face) = iter;
  return u_plus;
}

KOKKOS_FUNCTION
int compute_k_eps(double& k, double& eps ,
                  const double yp, const double u_star,
                  const double d_visco, const double dist, const double Cmu, const double Kappa)
{
  // PQ : 05/04/07 : formulation continue de k et epsilon
  //  assurant le bon comportement asymptotique
  const double u_star_carre = u_star * u_star;
  k    = 0.07*yp*yp*(exp(-yp/9.));
  k += 1./sqrt(Cmu)*(1. - exp(-yp/20.))*(1. - exp(-yp/20.));  // k_plus
  k   *= u_star_carre;
  // PL: 50625=15^4 on evite d'utiliser pow car lent
  eps  = 1./(Kappa*pow(yp*yp*yp*yp+50625, 0.25));  // eps_plus
  eps *= u_star_carre*u_star_carre/d_visco;

  return 1;
}

double Paroi_std_hyd_VEF::calculer_u_plus(const int ind_face,
                                          const double u_plus_d_plus,
                                          const double erugu)
{
  double r, u_plus;
  int iter;
  if (!compute_u_plus(ind_face, u_plus_d_plus, erugu, Kappa_, r, iter, u_plus))
    erreur_non_convergence();
  seuil_LP_(ind_face) = r;
  iterations_LP_(ind_face) = iter;
  return u_plus;
}

int Paroi_std_hyd_VEF::init_lois_paroi()
{
  tab_u_star_.resize(le_dom_VEF->nb_faces_tot());
  tab_d_plus_.resize(le_dom_VEF->nb_faces_tot());
  uplus_.resize(le_dom_VEF->nb_faces_tot());
  if (!Cisaillement_paroi_.get_md_vector().non_nul())
    {
      Cisaillement_paroi_.resize(0, dimension);
      le_dom_VEF->creer_tableau_faces(Cisaillement_paroi_);
    }
  seuil_LP_.resize(le_dom_VEF->nb_faces_tot());
  iterations_LP_.resize(le_dom_VEF->nb_faces_tot());

  check_turbulence_model();

  return init_lois_paroi_hydraulique();
}

// Remplissage de la table

int Paroi_std_hyd_VEF::init_lois_paroi_hydraulique()
{
  Cmu_ = mon_modele_turb_hyd->get_Cmu();
  init_lois_paroi_hydraulique_();
  return 1;
}

// face_keps_imposee_ constitue une connectivite entre face connectee a un bord et face de bord
//-pour une face pas connectee a un bord portant une CL de Dirichlet ou n appartenant pas a ce type de bord : -2
//-pour une face appartenant a un bord portant une CL de Dirichlet : -1
//-pour une face connectee a un bord portant une CL de Dirichlet : numero_de_face de bord la plus proche

void remplir_face_keps_imposee_gen(int& flag_face_keps_imposee_,
                                   IntVect& face_keps_imposee_,
                                   int& methode_calcul_face_keps_impose_,
                                   const Domaine_VEF& domaine_VEF,
                                   const REF(Domaine_Cl_VEF) le_dom_Cl_VEF,
                                   int is_champ_P1NC)
{
  if (flag_face_keps_imposee_==0)
    {
      flag_face_keps_imposee_=1;
      Cerr<<"construction de fac_keps_imposee_ par la methode generale"<<finl;
      int nb_faces_tot=domaine_VEF.nb_faces_tot();

      // Dimensionnement et initialisation de face_keps_imposee_
      if (!face_keps_imposee_.get_md_vector().non_nul())
        domaine_VEF.creer_tableau_faces(face_keps_imposee_, RESIZE_OPTIONS::NOCOPY_NOINIT);
      face_keps_imposee_=-2;

      // Remplissage de face_bords_diri donnant les faces
      // de bord avec une CL de Dirichlet
      int compt=0;
      ArrOfInt face_bords_diri(nb_faces_tot);
      for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
        {
          const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
          if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) ||
              (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ) ||
              (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE"))
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int ndeb = 0;
              int nfin = le_bord.nb_faces_tot();
              // On traite les faces de bord reelles.
              for (int ind_face=ndeb; ind_face<nfin; ind_face++)
                {
                  face_bords_diri[compt++]= le_bord.num_face(ind_face);
                }
            }
        }
      face_bords_diri.resize_array(compt);

      // Remplissage du tableau is_sommet_sur_bord qui contient
      // les sommets appartenant a une face de bord avec une CL Dirichlet
      ArrOfInt is_sommet_sur_bord(domaine_VEF.nb_som_tot());
      const IntTab& face_sommets=domaine_VEF.face_sommets();
      int nb_som_face=face_sommets.line_size();
      for (int fac=0; fac<compt; fac++)
        {
          int face=face_bords_diri[fac];
          for (int som=0; som<nb_som_face; som++)
            {
              int sommet=face_sommets(face,som);
              is_sommet_sur_bord[sommet]++;
            }
        }

      //On exploite is_sommet_sur_bord pour fixer la taille des List que va repertorier som_face_bord
      //-som_face_bord contient nb_som_tot List
      //-chaque List est dimensionnee par le nombre de faces de bord (portant une CL de Dirichlet)
      // connectees au sommet considere
      //-une List donnee pour un sommet repertorie le numero des faces qui lui sont connectees

      Static_Int_Lists som_face_bord;
      som_face_bord.set_list_sizes(is_sommet_sur_bord);
      is_sommet_sur_bord=0;
      for (int fac=0; fac<compt; fac++)
        {
          int face=face_bords_diri[fac];
          for (int som=0; som<nb_som_face; som++)
            {
              int sommet=face_sommets(face,som);
              int n=(is_sommet_sur_bord[sommet])++;
              som_face_bord.set_value(sommet,n,face);
            }
        }

      //-on parcourt toutes les faces du domaine pour identifier celles qui sont connectees a une face de bord
      // portant une CL de Dirichlet
      //-pour une face connectee on determine le nombre de sommets appartenant au(x) bord(s) (ils sont listes dans test)
      //int nb_faces_testes=nb_faces_tot;
      ArrOfInt traite_face(nb_faces_tot);
      if (methode_calcul_face_keps_impose_==5)
        {
          traite_face=0;
          // on marque que les faces sur les elements acrroches..
          const IntTab& face_voisins = domaine_VEF.face_voisins();
          int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
          // on parcourt les bords dirichlet, on cherche les voisins...
          for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
            {
              const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
              if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) ||
                  (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ) ||
                  (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE"))
                {
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
                  const IntTab& elem_faces = domaine_VEF.elem_faces();
                  int ndeb = 0;
                  int nfin = le_bord.nb_faces_tot();
                  // On traite toutes les faces de bord
                  for (int ind_face=ndeb; ind_face<nfin; ind_face++)
                    {
                      int num_face = le_bord.num_face(ind_face);
                      int elem = face_voisins(num_face,0);
                      for (int i=0; i<nb_faces_elem; i++)
                        traite_face[elem_faces(elem,i)]=1;
                    }
                }
            }
        }
      else
        {
          traite_face=1;
          assert(methode_calcul_face_keps_impose_==1);
        }
      const DoubleTab& xv=domaine_VEF.xv();
      for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
        if (traite_face[ind_face])
          {
            // on regarde combien de sommets sont en commun avec le bord
            IntList test;
            for (int som=0; som<nb_som_face; som++)
              {
                int sommet=face_sommets(ind_face,som);
                if (is_sommet_sur_bord[sommet])
                  test.add(sommet);
              }
            int test_size=test.size();
            if (test_size==0) //||(test_size==nb_som_face))
              {
                ; // on a rien a faire pour ce type de face (touche pas le bord ou sont confondus avec une face de bord (et non cas tordu ou une face est constituee que de sommets de bord mais pas sur le meme bord...
              }
            else
              {
                // la face ind_face est connectee a un bord associe a une CL de Dirichlet
                // on va determiner l'ensemble des faces de bord auxquelles elles est connectee
                // -on selectionne une face de bord f_pos parmi celles de la liste des faces associees
                // au premier sommet de la liste test
                // -on recherche ensuite si cette face est connue dans chacune des listes de faces associees
                // avec les autres sommets de la liste test. Si c est le cas, cette face est retenue comme possible
                // pour etre connectee a ind_face

                IntList possible;
                for (int i=0; i<som_face_bord.get_list_size(test[0]); i++)
                  {
                    int f_pos=som_face_bord(test[0],i);
                    int ok=1;
                    for (int lis=1; lis<test_size; lis++)
                      {
                        int som=test[lis];
                        int ok2=0;
                        for (int k=0; k<som_face_bord.get_list_size(som); k++)
                          {
                            int f2=som_face_bord(som,k);
                            if (f2==f_pos) ok2=1;
                          }
                        if (ok2==0) ok=0;
                      }
                    if (ok==1) possible.add(f_pos);
                  }

                //Dans le cas ou la face traitee appartient a deux bords differents portant une CL de type Dirichlet
                //il se peut que tous les sommets connectes n aient pas une face de bord commune et par consequent
                //possible est vide
                //Dans cette situation, on construit la liste possible a partir de toutes les faces contenues
                //dans les listes identifiees par som_face_bord pour les sommets contenus dans test

                if ((possible.size()==0))
                  {
                    //Cerr<<"Attention choix par la distance pour la face "<<ind_face<<finl;
                    for (int k=0; k<test_size; k++)
                      for (int i=0; i<som_face_bord.get_list_size(test[k]); i++)
                        {
                          int f_pos=som_face_bord(test[k],i);

                          possible.add_if_not(f_pos);
                        }
                  }

                assert(possible.size());
                int combien_de_face_possible=possible.size();

                //On remplit maintenant face_keps_imposee
                // -soit il n y a qu une face possible et on retient celle-ci
                // -soit il y en a plusieurs et on retient celle qui satisfait les deux criteres suivants :
                // 1 : la distance normale entre ind_face et face_possible est au moins 10 % de la distance totale
                // 2 : Si plusieurs faces possibles satisfont le critere 1 on retient celle est a la plus petite distance
                // Dans le cas ou il reste plusieurs faces (qui sont donc a la meme distance de ind_face),
                // on recherche celle qui est a une distance minimum d un face fictive obtenue en translatant ind_face.
                // La translation est effectuee dans les trois directions jusqu a obtention d une face possible unique.
                // Ce tri doit assurer le parallelisme quant au choix de la face a connecter a ind_face.

                if (combien_de_face_possible==1)
                  face_keps_imposee_[ind_face]=possible[0];
                else
                  {
                    //On supprime les faces possibles qui ne respecte pas le critere de distance normale
                    int ind=-1;
                    int elem_suppr=0;
                    int pos=0;
                    int size_initiale = possible.size();
                    while (pos<size_initiale)
                      {
                        int new_pos = pos-elem_suppr;
                        int f_new_pos=possible[new_pos];
                        double dist=0;
                        double dist2=distance_face(ind_face,f_new_pos,domaine_VEF);
                        for (int dir=0; dir<Objet_U::dimension; dir++)
                          {
                            dist+=(xv(f_new_pos,dir)-xv(ind_face,dir))*(xv(f_new_pos,dir)-xv(ind_face,dir));
                          }
                        if (dist2<=0.1*dist)
                          {
                            possible.suppr(f_new_pos);
                            elem_suppr++;
                          }
                        pos++;
                      }

                    if (possible.size()!=0)
                      {
                        //Parmi les faces possibles restantes on cherche celle qui est a une distance minimale
                        //ce choix etant fait en assurant le parallelisme

                        int dimension = Objet_U::dimension;
                        DoubleTab coord_possible;
                        DoubleVect coord_ref;
                        int size_liste = possible.size();
                        coord_possible.resize(size_liste,dimension,RESIZE_OPTIONS::NOCOPY_NOINIT);
                        coord_ref.resize(dimension,RESIZE_OPTIONS::NOCOPY_NOINIT);
                        coord_ref(0) = xv(ind_face,0);
                        coord_ref(1) = xv(ind_face,1);
                        if (dimension>2)
                          coord_ref(2) = xv(ind_face,2);;

                        for (int i_face=0; i_face<size_liste; i_face++)
                          {
                            int f_pos = possible[i_face];
                            for (int dir=0; dir<dimension; dir++)
                              coord_possible(i_face,dir) = xv(f_pos,dir);
                          }
                        Domaine::identifie_item_unique(possible,coord_possible,coord_ref);
                        ind = possible[0];
                      }

                    // Si pas de face trouvee on considere que c est parce que le critere de normalite
                    // n est pas respecte. Dans ce cas la on ne connecte pas la face
                    if (ind==-1)
                      ind = -2;
                    face_keps_imposee_[ind_face] = ind;
                  }
              }
          }

      // on reparcourt les bords pour remettre les bon flags
      for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
        {
          const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = 0;
          int nfin = le_bord.nb_faces_tot();
          if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) ||
              (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ) ||
              (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE"))
            {
              for (int ind_face=ndeb; ind_face<nfin; ind_face++)
                {
                  int num_face=le_bord.num_face(ind_face);
                  face_keps_imposee_[num_face]=-1;
                }
            }
          else if (sub_type(Periodique,la_cl.valeur()) || sub_type(Neumann_sortie_libre,la_cl.valeur()) || sub_type(Symetrie,la_cl.valeur()))
            {
              //exit();
            }
          else
            {
              for (int ind_face=ndeb; ind_face<nfin; ind_face++)
                {
                  int num_face=le_bord.num_face(ind_face);
                  face_keps_imposee_[num_face]=-2;
                }
            }
        }
      Cerr<<"fin construction de fac_keps_imposee_ par la methode generale"<<finl;
    }
}

void remplir_face_keps_imposee(int& flag_face_keps_imposee_,
                               int methode_calcul_face_keps_impose_,
                               IntVect& face_keps_imposee_,
                               const Domaine_VEF& domaine_VEF,
                               const REF(Domaine_Cl_VEF) le_dom_Cl_VEF,
                               int is_champ_P1NC)
{
  // Deja calcule, on quitte
  if (flag_face_keps_imposee_ == 1)
    return;

  if (methode_calcul_face_keps_impose_ == -1)
    {
      if (Objet_U::dimension == 2)
        methode_calcul_face_keps_impose_ = 1;
      else
        {
          assert(Objet_U::dimension == 3);
          methode_calcul_face_keps_impose_ = 5 ;
          // GF pour avoir comme dans la 156 decommenter la ligne suivante
          // methode_calcul_face_keps_impose_=1 ;
        }
    }
  Cout <<"Calcul des faces ou keps est impose par la methode : ";
  if (methode_calcul_face_keps_impose_==1)
    Cout<<"toutes_les_faces_accrochees"<<finl;
  else  if (methode_calcul_face_keps_impose_==4)
    Cout<<"que_les_faces_des_elts_dirichlet"<<finl;
  else if (methode_calcul_face_keps_impose_==5)
    Cout<<"que_les_faces_des_elts_dirichlet_et_test_distance"<<finl;
  else
    {
      Cout <<" Cas non prevu "<<finl;
      // GF Les cas 0 2 3 sont appeles a disparaitre du code....
      Process::exit();
    }

  if (flag_face_keps_imposee_ == 0)
    {
      domaine_VEF.creer_tableau_faces(face_keps_imposee_, RESIZE_OPTIONS::NOCOPY_NOINIT);
      face_keps_imposee_ = -2;
    }
  if ((methode_calcul_face_keps_impose_==1) || (methode_calcul_face_keps_impose_==5))
    remplir_face_keps_imposee_gen(flag_face_keps_imposee_,
                                  face_keps_imposee_,
                                  methode_calcul_face_keps_impose_,
                                  domaine_VEF,
                                  le_dom_Cl_VEF,
                                  is_champ_P1NC);
  else
    {
      // methode_calcul_face_keps_impose_=4
      flag_face_keps_imposee_ = 1;
      const IntTab& face_voisins = domaine_VEF.face_voisins();
      // on parcourt les bords dirichlet, on cherche les voisins...
      for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
        {
          const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
          if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) ||
              (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ) ||
              (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE"))
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              const IntTab& elem_faces = domaine_VEF.elem_faces();
              int ndeb = 0;
              int nfin = le_bord.nb_faces_tot();
              // On traite toutes les faces de bord
              for (int ind_face=ndeb; ind_face<nfin; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  int elem=face_voisins(num_face,0);
                  const Domaine& domaine = domaine_VEF.domaine();
                  int nfac = domaine.nb_faces_elem();
                  int nf2;
                  ArrOfInt num(nfac);
                  for (nf2=0; nf2<nfac; nf2++)
                    {
                      num[nf2] = elem_faces(elem,nf2);
                    }
                  for (nf2=0; nf2<nfac-1; nf2++)
                    {
                      num[nf2] = elem_faces(elem,nf2);
                      if (num[nf2] == num_face)
                        {
                          num[nf2] = num[nfac-1];
                          num[nfac-1] = num_face;
                        }
                    }

                  for (int nf=0; nf<nfac; nf++)
                    {
                      int num0=num[nf];//
                      num0=elem_faces(elem,nf);

                      if (num0==num_face)
                        face_keps_imposee_[num0]=-1;
                      else
                        {
                          // Ce sont les faces directement connectes a num_face
                          //                          int face_asso=face_keps_imposee_[num2];
                          //          if (face_asso<0)
                          face_keps_imposee_[num0]=num_face;
                          // on cherche les faces accroches
                          int elem_voisin;
                          if (face_voisins(num0,0)!=elem)
                            elem_voisin=face_voisins(num0,0);
                          else
                            elem_voisin=face_voisins(num0,1);

                          if (elem_voisin!=-1)
                            {
                              if ((methode_calcul_face_keps_impose_==0)||(methode_calcul_face_keps_impose_==3))
                                {
                                  if (Objet_U::dimension==2)
                                    {
                                      const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
                                      int num1=elem_faces(elem_voisin,0);
                                      int num2=elem_faces(elem_voisin,1);
                                      if (num1==num0) num1=elem_faces(elem_voisin,2);
                                      else if (num2==num0) num2=elem_faces(elem_voisin,2);
                                      if (rang_elem_non_std(elem_voisin)==-1)

                                        {
                                          double dist1=distance_face(num_face,num1,domaine_VEF);
                                          double dist2=distance_face(num_face,num2,domaine_VEF);

                                          // Prend seulement le sommet le plus pres de la paroi :
                                          int numb;
                                          if (dist1<dist2)
                                            {
                                              numb=num1;
                                              //  dist=dist1;
                                            }
                                          else
                                            {
                                              numb=num2;
                                              //dist=dist2;
                                            }
                                          // On regarde si la face num0 appartient a un element
                                          // qui est en contact avec un bord
                                          int elem0 = face_voisins(numb,0);
                                          int elem1 = face_voisins(numb,1);
                                          int test_rang;
                                          if ((elem0 == elem_voisin) && (elem1>=0))
                                            test_rang = rang_elem_non_std(elem1);
                                          else if ((elem1 == elem_voisin) && (elem0>=0))
                                            test_rang = rang_elem_non_std(elem0);
                                          else
                                            test_rang = 1;

                                          if (test_rang<0)
                                            face_keps_imposee_[numb]=num_face;
                                        }
                                    }
                                  else if (Objet_U::dimension==3)
                                    {
                                      int num1=elem_faces(elem_voisin,0);
                                      int num2=elem_faces(elem_voisin,1);
                                      int num3=elem_faces(elem_voisin,2);
                                      if (num1==num0) num1=elem_faces(elem_voisin,3);
                                      else if (num2==num0) num2=elem_faces(elem_voisin,3);
                                      else if (num3==num0) num3=elem_faces(elem_voisin,3);
                                      // if (rang_elem_non_std(elem_voisin)==-1)
                                      int ndebint = domaine_VEF.premiere_face_int();
                                      if ((methode_calcul_face_keps_impose_==3)||(num1>=ndebint && num2>=ndebint && num3>=ndebint))
                                        {
                                          // elem_voisin n'est pas un element de bord
                                          double dist1=distance_face(num_face,num1,domaine_VEF);
                                          double dist2=distance_face(num_face,num2,domaine_VEF);
                                          double dist3=distance_face(num_face,num3,domaine_VEF);
                                          // Prend seulement les deux sommets les plus pres de la paroi :
                                          if (dist2<dist1)
                                            {
                                              int n=num1;
                                              double d=dist1;
                                              num1=num2;
                                              dist1=dist2;
                                              num2=n;
                                              dist2=d;
                                            }
                                          if (dist3<dist2)
                                            {
                                              int n=num2;
                                              double d=dist2;
                                              num2=num3;
                                              dist2=dist3;
                                              num3=n;
                                              dist3=d;
                                            }
                                          if (dist1 > 0)
                                            {
                                              face_keps_imposee_[num1] = num_face;
                                            }
                                          if (dist2 > 0)
                                            {
                                              face_keps_imposee_[num2] = num_face;
                                            }
                                        }
                                    }
                                }
                              else if (methode_calcul_face_keps_impose_ == 2)
                                for (int dir = 0; dir < nfac; dir++)
                                  {
                                    assert(methode_calcul_face_keps_impose_ == 2);
                                    const int num3 = elem_faces(elem_voisin, dir);
                                    if ((num3 > 0) && (num3 != num0))
                                      {
                                        int face_asso = face_keps_imposee_[num3];
                                        if (face_asso < 0)
                                          face_keps_imposee_[num3] = num_face;

                                        else
                                          {
                                            // face deja associe
                                            // quelle est la plus pres de face_asso et num_face ...  non fait pour l'instant
                                            face_keps_imposee_[num3] = num_face;
                                          }
                                      }

                                  }
                              else
                                {
                                  // on ne rajoute pas de faces
                                  assert(methode_calcul_face_keps_impose_ == 4);
                                }

                            }
                        }
                    }
                }
            }
        }

      // on reparcourt les bords pour remettre les bon flags
      for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
        {
          const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          const int ndeb = 0;
          const int nfin = le_bord.nb_faces_tot();
          if (sub_type(Dirichlet_paroi_fixe, la_cl.valeur())
              || (sub_type(Dirichlet_paroi_defilante, la_cl.valeur()))
              || (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE"))
            {
              for (int ind_face = ndeb; ind_face < nfin; ind_face++)
                {
                  const int num_face = le_bord.num_face(ind_face);
                  face_keps_imposee_[num_face] = -1;
                }
            }
          else if (sub_type(Periodique, la_cl.valeur())
                   || sub_type(Neumann_sortie_libre, la_cl.valeur())
                   || sub_type(Symetrie, la_cl.valeur()))
            {
              //exit();
            }
          else
            {
              for (int ind_face = ndeb; ind_face < nfin; ind_face++)
                {
                  const int num_face = le_bord.num_face(ind_face);
                  face_keps_imposee_[num_face] = -2;
                }
            }
        }
    }

  Debog::verifier_indices_items("Paroi_std_hyd_VEF::remplir_face_keps_imposee face_keps_imposee",
                                domaine_VEF.face_sommets().get_md_vector(),
                                face_keps_imposee_);
  // Il ne trouve pas les memes faces en sequentiel et parallele d'ou
  // des legers ecarts en sequentiel et parallele
}

int Paroi_std_hyd_VEF::calculer_hyd_BiK(DoubleTab& tab_k,DoubleTab& tab_eps)
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();
  const DoubleTab& vitesse = eqn_hydr.inconnue().valeurs();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  const Domaine& domaine = domaine_VEF.domaine();
  int nfac = domaine.nb_faces_elem();

  double visco=-1;
  int l_unif;
  if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
    {
      visco = std::max(tab_visco(0,0),DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;
  if ((!l_unif) && (tab_visco.local_min_vect()<DMINFLOAT))
    //   on ne doit pas changer tab_visco ici !
    {
      Cerr<<" visco <=0 ?"<<finl;
      exit();
    }
  // tab_visco+=DMINFLOAT;

  double norm_v=-1;
  double dist=-1,d_visco=-1;
  double u_plus_d_plus,u_plus,d_plus,u_star;
  ArrOfDouble val(dimension);
  Cisaillement_paroi_=0;

  int is_champ_Q1NC=sub_type(Champ_Q1NC,eqn_hydr.inconnue());
  remplir_face_keps_imposee( flag_face_keps_imposee_, methode_calcul_face_keps_impose_, face_keps_imposee_, domaine_VEF, le_dom_Cl_VEF, !is_champ_Q1NC);

  IntVect num(nfac);
  ArrOfDouble stock_erugu(domaine_VEF.nb_faces_tot());
  ArrOfInt is_defilante_face(domaine_VEF.nb_faces_tot());

  // Loop on boundaries
  int nb_bords=domaine_VEF.nb_front_Cl();
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);

      // Only Dirichlet conditions:
      //if (sub_type(Dirichlet,la_cl.valeur()) ||
      //    (sub_type(Dirichlet,la_cl.valeur())) ||
      if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) ||
          (sub_type(Dirichlet_paroi_defilante,la_cl.valeur())) ||
          (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE") )
        {
          int is_defilante=sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ;

          if(la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE")
            {
              is_defilante = (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE");
            }

          // Recuperation de la valeur Erugu
          double erugu=Erugu;
          if (sub_type(Paroi_rugueuse,la_cl.valeur()))
            erugu=ref_cast(Paroi_rugueuse,la_cl.valeur()).get_erugu();

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          const IntTab& elem_faces = domaine_VEF.elem_faces();

          // Loop on real faces
          int ndeb = 0;
          int nfin = le_bord.nb_faces_tot();
          for (int ind_face=ndeb; ind_face<nfin; ind_face++)
            {
              int num_face=le_bord.num_face(ind_face);
              int elem=face_voisins(num_face,0);

              is_defilante_face[num_face]=is_defilante;
              stock_erugu[num_face]=erugu;

              // on determine les faces de l'element qui sont autres que le num_face traite
              for (int nf2=0; nf2<nfac; nf2++)
                num[nf2] = elem_faces(elem,nf2);

              // Maintenant on place le num_face en fin de tableau
              for (int nf2=0; nf2<nfac-1; nf2++)
                {
                  num[nf2] = elem_faces(elem,nf2);
                  if (num[nf2] == num_face)
                    {
                      num[nf2] = num[nfac-1];
                      num[nfac-1] = num_face;
                    }
                }

              int nfc=0;
              // Boucle sur les faces :
              for (int nf=0; nf<nfac; nf++)
                {
                  if (num[nf]==num_face)
                    {
                      // Strategie pour les tetras :
                      // On impose k et eps a la paroi :
                      // approximation: d(k)/d(n) = 0 a la paroi
                      // c'est faux mais ca marche
                      tab_k(num[nf])  =0.;
                      tab_eps(num[nf])=0.;
                      int nk=0;

                      for (int k=0; k<nfac; k++)
                        //if ( (num[k] >= ndebint) && (k != nf))
                        if ( (face_keps_imposee_[num[k]]>-1) && (k != nf))
                          {

                            tab_k(num[nf])  += tab_k(num[k]);
                            tab_eps(num[nf])+= tab_eps(num[k]);
                            nk++;
                          }
                      if (nk != 0 )
                        {
                          tab_k(num[nf])  /=nk;
                          tab_eps(num[nf])/=nk;
                        }
                    }

                  // On verifie si num[nf] n'est pas une face de bord :
                  else if ( (face_keps_imposee_[num[nf]]>-1))//if (num[nf]>=ndebint)
                    {
                      nfc++;
                      norm_v=norm_vit_lp_k(vitesse,num[nf],num_face,domaine_VEF,val,is_defilante);
                      if (!is_champ_Q1NC)
                        dist=distance_face(num_face,num[nf],domaine_VEF);
                      else
                        dist=distance_face_elem(num_face,elem,domaine_VEF);

                      if (l_unif)
                        d_visco = visco;
                      else
                        d_visco = tab_visco(elem,0);

                      u_plus_d_plus = norm_v*dist/d_visco;
                      u_plus = calculer_u_plus(nf,u_plus_d_plus,erugu);

                      if (!is_u_star_impose_)
                        {
                          if(u_plus)
                            {
                              u_star = norm_v/u_plus ;
                              d_plus = u_plus_d_plus/u_plus ;
                            }
                          else
                            {
                              u_star = 0.;
                              d_plus = 0.;
                            }
                        }
                      else
                        {
                          u_star = u_star_impose_  ;
                          d_plus = 0.;
                        }
                      calculer_k_eps(tab_k(num[nf]),tab_eps(num[nf]),d_plus,u_star,d_visco,dist);

                      // Calcul de la contrainte tangentielle
                      for (int dir=0; dir<dimension; dir++)
                        Cisaillement_paroi_(num_face,dir) += u_star*u_star*val[dir];
                      // Fin de la strategie du calcul generalise de la loi de paroi
                    }
                }

              // A voir si juste :
              if (nfc != 0 )
                for (int dir=0; dir<dimension; dir++)
                  Cisaillement_paroi_(num_face,dir)/=nfc;

              double u_starbis=0;
              for (int dir=0; dir<dimension; dir++)
                u_starbis+=Cisaillement_paroi_(num_face,dir)*Cisaillement_paroi_(num_face,dir);

              u_starbis=sqrt(sqrt(u_starbis));
              tab_u_star_(num_face)=u_starbis;
              tab_d_plus_(num_face)=u_starbis*dist/d_visco;
              if (u_starbis != 0) uplus_(num_face)=norm_v/u_starbis;

            } // End loop on real faces

        } // End hlet conditions

      // Robin condition:
      else if (sub_type(Paroi_decalee_Robin,la_cl.valeur()))
        {
          // Recuperation de la valeur Erugu
          const double erugu = Erugu;

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          const Paroi_decalee_Robin& Paroi = ref_cast(Paroi_decalee_Robin,la_cl.valeur());
          const DoubleTab& normales = domaine_VEF.face_normales();
          const double delta = Paroi.get_delta();

          // Loop on real faces
          const int ndeb = 0;
          const int nfin = le_bord.nb_faces_tot();
          for (int ind_face=ndeb; ind_face<nfin; ind_face++)
            {
              const int num_face=le_bord.num_face(ind_face);
              const int elem=face_voisins(num_face,0);

              double psc=0, norm=0;
              norm_v=0;

              for(int comp=0; comp<dimension; comp++)
                {
                  psc += vitesse(num_face,comp)*normales(num_face,comp);
                  norm += normales(num_face,comp)*normales(num_face,comp);
                }
              // psc /= norm; // Fixed bug: Arithmetic exception
              if (std::fabs(norm)>=DMINFLOAT) psc/=norm;

              for(int comp=0; comp<dimension; comp++)
                {
                  val[comp]=vitesse(num_face,comp)-psc*normales(num_face,comp);
                  norm_v += val[comp]*val[comp];
                }

              norm_v = sqrt(norm_v);
              val /= norm_v;
              dist = delta;

              // Common to Dirichlet

              if (l_unif)
                d_visco = visco;
              else
                d_visco = tab_visco(elem,0);

              u_plus_d_plus = norm_v*dist/d_visco;

              u_plus = calculer_u_plus(ind_face,u_plus_d_plus,erugu);

              if (!is_u_star_impose_)
                {
                  if(u_plus)
                    {
                      u_star = norm_v/u_plus ;
                      d_plus = u_plus_d_plus/u_plus ;
                    }
                  else
                    {
                      u_star = 0.;
                      d_plus = 0.;
                    }
                }
              else
                {
                  u_star = u_star_impose_;
                  d_plus = 0.;
                }

              calculer_k_eps(tab_k(num_face),tab_eps(num_face),d_plus,u_star,d_visco,dist);

              // Calcul de la contrainte tangentielle
              for (int j=0; j<dimension; j++)
                Cisaillement_paroi_(num_face,j) = u_star*u_star*val[j];

              double u_starbis=0;
              for (int dir=0; dir<dimension; dir++)
                u_starbis+=Cisaillement_paroi_(num_face,dir)*Cisaillement_paroi_(num_face,dir);
              u_starbis=sqrt(sqrt(u_starbis));

              tab_u_star_(num_face)=u_starbis;
              tab_d_plus_(num_face)=u_starbis*dist/d_visco;
              if (u_starbis != 0) uplus_(num_face)=norm_v/u_starbis;

              // End common to Dirichlet

            } // End loop on real faces

        } // End Robin condition

    } // End loop on boundaries


#ifdef CONTROL
  if ((flag_face_keps_imposee_)&&(!is_champ_Q1NC))
    {
      if (0)
        {
          int nb_faces_tot=domaine_VEF.nb_faces_tot();
          DoubleVect toto(nb_faces_tot);
          const ArrOfInt& renum=Debog::renum_faces();
          for (int i=0; i<nb_faces_tot; i++)
            if ((face_keps_imposee_(i)>-1)&&(Debog::mode_db==1))
              {
                if (i>1000) Cerr<<me()<<" face_virt "<<i<<finl;
                toto(i)=renum(face_keps_imposee_(i));
              }
            else toto(i)=face_keps_imposee_(i);
          barrier();
          Debog::verifier("face_keps ",toto);
        }
      int tutu=0;
      ArrOfInt test;
      remplir_face_keps_imposee_gen( tutu, test, domaine_VEF,le_dom_Cl_VEF,!is_champ_Q1NC);
      test-=face_keps_imposee_;
      if (max(test)>0|| min(test)<0)
        {
          const DoubleTab& xv=domaine_VEF.xv();
          Cerr<<"TEST "<<finl;
          int compteur=0,compteur2=0;
          for (int i=0; i<test.size_array(); i++)
            {
              test(i)+=face_keps_imposee_(i);
              if (test(i)!=face_keps_imposee_(i))
                {

                  Cerr<<me()<<" face "<<i<<" : " <<face_keps_imposee_(i)<<" "<<test(i)<<" pos "<<xv(i,0)<<" "<<xv(i,1);
                  if (dimension==3) Cerr<<" "<<xv(i,2);

                  if (face_keps_imposee_(i)!=-2)
                    {
                      Cerr<<" pos trouve "<<xv(face_keps_imposee_(i),0)<<" "<<xv(face_keps_imposee_(i),1);
                      if (dimension==3) Cerr<<" "<<xv(face_keps_imposee_(i),2);
                    }
                  if (test(i)!=-2)
                    {
                      Cerr<<" pos test "<<xv(test(i),0)<<" "<<xv(test(i),1);
                      if (dimension==3) Cerr<<" "<<xv(test(i),2);
                    }
                  Cerr<<finl;
                }
              if (test(i)>-1) compteur++;
              if (face_keps_imposee_(i)>-1) compteur2++;
            }
          Cerr<<"compteurs "<<compteur2<<" "<<compteur<<finl;
          exit();
        }
    }
#endif

  // on recalcule partout ou c'est impose
  int nb_faces_tot=domaine_VEF.nb_faces_tot();
  for (int face=0; face<nb_faces_tot; face++)
    {
      int num_face=face_keps_imposee_[face];
      if (num_face>-1)
        {
          // int elem_voisin;
          int elem=face_voisins(num_face,0);
          //if (face_voisins(face,0)!=elem) elem_voisin=face_voisins(face,0);
          //else elem_voisin=face_voisins(face,1);
          //int elem_voisin=face_voisins(num_face,0);
          // ce n'est pas le bon voisin!!!!
          double distb;
          if (!is_champ_Q1NC)
            {
              distb=distance_face(num_face,face,domaine_VEF);
            }
          else
            {
              assert(sub_type(Champ_Q1NC,eqn_hydr.inconnue())) ;
              distb=distance_face_elem(num_face,elem,domaine_VEF);
            }

          norm_v=norm_vit_lp_k(vitesse,face,num_face,domaine_VEF,val,is_defilante_face[num_face]);
          if (l_unif)
            d_visco = visco;
          else
            d_visco = tab_visco(elem,0);

          u_plus_d_plus = norm_v*distb/d_visco;
          u_plus = calculer_u_plus(face,u_plus_d_plus,stock_erugu[num_face]);

          if(u_plus)
            {
              u_star = norm_v/u_plus ;
              d_plus = u_plus_d_plus/u_plus ;
            }
          else
            {
              u_star = 0.            ;
              d_plus = 0.                   ;
            }

          calculer_k_eps(tab_k(face),tab_eps(face),d_plus,u_star,d_visco,distb);
        }
    }

  Cisaillement_paroi_.echange_espace_virtuel();
  tab_k.echange_espace_virtuel();
  tab_eps.echange_espace_virtuel();
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd tab_k",tab_k);
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd tab_eps",tab_eps);
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd Cisaillement_paroi_",Cisaillement_paroi_);
  return 1;
}  // fin de calcul_hyd_BiK (K-eps bicephale)

int Paroi_std_hyd_VEF::calculer_hyd(DoubleTab& tab_2eq)
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();
  const DoubleTab& vitesse = eqn_hydr.inconnue().valeurs();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  const Domaine& domaine = domaine_VEF.domaine();

  double visco {-1};
  int l_unif {0};
  if (sub_type(Champ_Uniforme, ch_visco_cin.valeur()))
    {
      visco = std::max(tab_visco(0, 0), DMINFLOAT);
      l_unif = 1;
    }

  if ((!l_unif) && (tab_visco.local_min_vect() < DMINFLOAT))
    //   on ne doit pas changer tab_visco ici !
    {
      Cerr << " visco <=0 ?" << finl;
      exit();
    }
  // tab_visco+=DMINFLOAT;

  double norm_v {-1};
  double dist {-1}, d_visco {-1};
  double u_plus_d_plus {0}, u_plus {0}, d_plus {0}, u_star {0};
  ArrOfDouble val(dimension);
  Cisaillement_paroi_ = 0;

  int is_champ_Q1NC = sub_type(Champ_Q1NC, eqn_hydr.inconnue());
  remplir_face_keps_imposee(flag_face_keps_imposee_, methode_calcul_face_keps_impose_,
                            face_keps_imposee_, domaine_VEF,
                            le_dom_Cl_VEF, !is_champ_Q1NC);

  const int nfac = domaine.nb_faces_elem();
  IntVect num(nfac);
  ArrOfDouble stock_erugu(domaine_VEF.nb_faces_tot());
  ArrOfInt is_defilante_face(domaine_VEF.nb_faces_tot());

  // Loop on boundaries
  const int nb_bords = domaine_VEF.nb_front_Cl();
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);

      // Only Dirichlet conditions:
      if (sub_type(Dirichlet_paroi_fixe, la_cl.valeur())
          || (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()))
          || (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE"))
        {
          int is_defilante = sub_type(Dirichlet_paroi_defilante, la_cl.valeur()) ;

          if(la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE")
            is_defilante = (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE");

          // Recuperation de la valeur Erugu
          // double erugu = Erugu;
          // if (sub_type(Paroi_rugueuse, la_cl.valeur()))
          //   erugu = ref_cast(Paroi_rugueuse, la_cl.valeur()).get_erugu();

          const double erugu = sub_type(Paroi_rugueuse, la_cl.valeur())
                               ? ref_cast(Paroi_rugueuse, la_cl.valeur()).get_erugu()
                               : Erugu;

          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          const IntTab& elem_faces = domaine_VEF.elem_faces();

          // Loop on real faces
          const int ndeb = 0;
          const int nfin = le_bord.nb_faces_tot();
          for (int ind_face = ndeb; ind_face < nfin; ind_face++)
            {
              const int num_face = le_bord.num_face(ind_face);
              const int elem = face_voisins(num_face, 0);

              is_defilante_face[num_face] = is_defilante;
              stock_erugu[num_face] = erugu;

              // on determine les faces de l'element qui sont autres que le num_face traite
              for (int nf2 = 0; nf2 < nfac; nf2++)
                num[nf2] = elem_faces(elem, nf2);

              // Maintenant on place le num_face en fin de tableau
              for (int nf2 = 0; nf2 < nfac - 1; nf2++)
                {
                  num[nf2] = elem_faces(elem, nf2);
                  if (num[nf2] == num_face)
                    {
                      num[nf2] = num[nfac-1];
                      num[nfac-1] = num_face;
                    }
                }

              int nfc = 0;
              // Boucle sur les faces :
              for (int nf = 0; nf < nfac; nf++)
                {
                  // cAlan : split this function to avoid the ifs!
                  // cAlan : and remove this block from the loop as
                  // cAlan : num_face was placed at the end of the num tab
                  if (num[nf] == num_face)
                    {
                      // Strategie pour les tetras :
                      // On impose k et eps a la paroi :
                      // approximation: d(k)/d(n) = 0 a la paroi
                      // c'est faux mais ca marche
                      tab_2eq(num[nf], 0) = 0.;
                      if (turbulence_model_type == 1) // second equation is epsilon
                        tab_2eq(num[nf], 1) = 0.;
                      else if (turbulence_model_type == 2) // second equation is omega
                        {
                          // cAlan : c'est sale car on reprend la valeur de d_visco et de dist de
                          // l'itération précédente. A améliorer quand on optimisera cette fonction.
                          assert(dist > 0);
                          tab_2eq(num[nf], 1) = 6*d_visco/(BETA_OMEGA*dist*dist);
                        }
                      int nk = 0;

                      for (int k = 0; k < nfac; k++)
                        //if ( (num[k] >= ndebint) && (k != nf))
                        if ((face_keps_imposee_[num[k]] > -1) && (k != nf))
                          {

                            tab_2eq(num[nf], 0) += tab_2eq(num[k], 0);
                            tab_2eq(num[nf], 1) += tab_2eq(num[k], 1);
                            nk++;
                          }
                      if (nk != 0)
                        {
                          tab_2eq(num[nf], 0) /= nk;
                          tab_2eq(num[nf], 1) /= nk;
                        }
                    }

                  // On verifie si num[nf] n'est pas une face de bord :
                  else if ((face_keps_imposee_[num[nf]] > -1))//if (num[nf]>=ndebint)
                    {
                      nfc++;
                      norm_v = norm_vit_lp_k(vitesse,num[nf], num_face, domaine_VEF,
                                             val, is_defilante);

                      if (!is_champ_Q1NC)
                        dist = distance_face(num_face, num[nf],domaine_VEF);
                      else
                        dist = distance_face_elem(num_face, elem, domaine_VEF);

                      if (l_unif)
                        d_visco = visco;
                      else
                        d_visco = tab_visco(elem, 0);

                      u_plus_d_plus = norm_v*dist/d_visco;
                      u_plus = calculer_u_plus(nf, u_plus_d_plus, erugu);

                      if (!is_u_star_impose_)
                        {
                          if((bool)u_plus)
                            {
                              u_star = norm_v/u_plus ;
                              d_plus = u_plus_d_plus/u_plus ;
                            }
                          else
                            {
                              u_star = 0.;
                              d_plus = 0.;
                            }
                        }
                      else
                        {
                          u_star = u_star_impose_;
                          d_plus = 0.;
                        }
                      // calculer_k_eps(tab_k_eps(num[nf],0),tab_k_eps(num[nf],1),d_plus,u_star,d_visco,dist);

                      // Calcul de la contrainte tangentielle
                      for (int dir = 0; dir < dimension; dir++)
                        Cisaillement_paroi_(num_face, dir) += u_star*u_star*val[dir];
                      // Fin de la strategie du calcul generalise de la loi de paroi

                      // Evaluate the turbulent quantities
                      compute_turbulent_quantities(tab_2eq(num[nf], 0), tab_2eq(num[nf], 1),
                                                   d_plus, u_star, d_visco, dist);
                    }
                }

              // A voir si juste : // cAlan !! WTF ce commentaire
              if (nfc != 0)
                for (int dir = 0; dir < dimension; dir++)
                  Cisaillement_paroi_(num_face, dir) /= nfc;

              double u_starbis = 0;
              for (int dir = 0; dir < dimension; dir++)
                u_starbis += Cisaillement_paroi_(num_face, dir)*Cisaillement_paroi_(num_face, dir);

              u_starbis=sqrt(sqrt(u_starbis));
              tab_u_star_(num_face) = u_starbis;
              tab_d_plus_(num_face) = u_starbis*dist/d_visco;
              if (u_starbis != 0)
                uplus_(num_face) = norm_v/u_starbis;

            } // End loop on real faces

        } // End Dirichlet conditions

      // Robin condition:
      else if (sub_type(Paroi_decalee_Robin, la_cl.valeur()))
        {
          // Recuperation de la valeur Erugu
          double erugu = Erugu;

          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          const Paroi_decalee_Robin& Paroi = ref_cast(Paroi_decalee_Robin, la_cl.valeur());
          const DoubleTab& normales = domaine_VEF.face_normales();
          const double delta = Paroi.get_delta();

          // Loop on real faces
          const int ndeb = 0;
          const int nfin = le_bord.nb_faces_tot();
          for (int ind_face = ndeb; ind_face < nfin; ind_face++)
            {
              const int num_face = le_bord.num_face(ind_face);
              const int elem = face_voisins(num_face, 0);

              double psc = 0;
              double norm = 0;
              for(int comp = 0; comp < dimension; comp++)
                {
                  psc += vitesse(num_face, comp)*normales(num_face, comp);
                  norm += normales(num_face, comp)*normales(num_face, comp);
                }
              // psc /= norm; // Fixed bug: Arithmetic exception
              if (std::fabs(norm) >= DMINFLOAT)
                psc /= norm;

              norm_v = 0;
              for(int comp = 0; comp < dimension; comp++)
                {
                  val[comp] = vitesse(num_face, comp) - psc*normales(num_face, comp);
                  norm_v += val[comp]*val[comp];
                }

              norm_v = sqrt(norm_v);
              if (std::fabs(norm_v) >= DMINFLOAT)
                val /= norm_v;
              dist = delta;

              // Common to Dirichlet

              if (l_unif)
                d_visco = visco;
              else
                d_visco = tab_visco(elem, 0);

              u_plus_d_plus = norm_v*dist/d_visco;

              u_plus = calculer_u_plus(ind_face, u_plus_d_plus, erugu);

              if (!is_u_star_impose_)
                {
                  if((bool)u_plus)
                    {
                      u_star = norm_v/u_plus ;
                      d_plus = u_plus_d_plus/u_plus ;
                    }
                  else
                    {
                      u_star = 0.;
                      d_plus = 0.;
                    }
                }
              else
                {
                  u_star = u_star_impose_;
                  d_plus = 0.;
                }

              // calculer_k_eps(tab_2eq(num_face,0),tab_2eq(num_face,1),d_plus,u_star,d_visco,dist);

              // Calcul de la contrainte tangentielle
              for (int j = 0; j < dimension; j++)
                Cisaillement_paroi_(num_face, j) = u_star*u_star*val[j];

              // cAlan : pourquoi recalcule-t-on ustar à partir de tau_w ?
              double u_starbis = 0;
              for (int dir = 0; dir < dimension; dir++)
                u_starbis += Cisaillement_paroi_(num_face, dir)*Cisaillement_paroi_(num_face, dir);
              u_starbis = sqrt(sqrt(u_starbis));

              tab_u_star_(num_face) = u_starbis;
              tab_d_plus_(num_face) = u_starbis*dist/d_visco;
              if (u_starbis != 0)
                uplus_(num_face) = norm_v / u_starbis;

              // End common to Dirichlet

              // Evaluate the turbulent quantities
              compute_turbulent_quantities(tab_2eq(num_face, 0), tab_2eq(num_face, 1),
                                           d_plus, u_star, d_visco, dist);

            } // End loop on real faces

        } // End Robin condition

    } // End loop on boundaries


#ifdef CONTROL
  if ((flag_face_keps_imposee_)&&(!is_champ_Q1NC))
    {
      if (0)
        {
          int nb_faces_tot=domaine_VEF.nb_faces_tot();
          DoubleVect toto(nb_faces_tot);
          const ArrOfInt& renum=Debog::renum_faces();
          for (int i=0; i<nb_faces_tot; i++)
            if ((face_keps_imposee_(i)>-1)&&(Debog::mode_db==1))
              {
                if (i>1000) Cerr<<me()<<" face_virt "<<i<<finl;
                toto(i)=renum(face_keps_imposee_(i));
              }
            else toto(i)=face_keps_imposee_(i);
          barrier();
          Debog::verifier("face_keps ",toto);
        }
      int tutu=0;
      ArrOfInt test;
      remplir_face_keps_imposee_gen( tutu, test, domaine_VEF,le_dom_Cl_VEF,!is_champ_Q1NC);
      test-=face_keps_imposee_;
      if (std::max(test)>0|| std::min(test)<0)
        {
          const DoubleTab& xv=domaine_VEF.xv();
          Cerr<<"TEST "<<finl;
          int compteur=0,compteur2=0;
          for (int i=0; i<test.size_array(); i++)
            {
              test(i)+=face_keps_imposee_(i);
              if (test(i)!=face_keps_imposee_(i))
                {

                  Cerr<<me()<<" face "<<i<<" : " <<face_keps_imposee_(i)<<" "<<test(i)<<" pos "<<xv(i,0)<<" "<<xv(i,1);
                  if (dimension==3) Cerr<<" "<<xv(i,2);

                  if (face_keps_imposee_(i)!=-2)
                    {
                      Cerr<<" pos trouve "<<xv(face_keps_imposee_(i),0)<<" "<<xv(face_keps_imposee_(i),1);
                      if (dimension==3) Cerr<<" "<<xv(face_keps_imposee_(i),2);
                    }
                  if (test(i)!=-2)
                    {
                      Cerr<<" pos test "<<xv(test(i),0)<<" "<<xv(test(i),1);
                      if (dimension==3) Cerr<<" "<<xv(test(i),2);
                    }
                  Cerr<<finl;
                }
              if (test(i)>-1) compteur++;
              if (face_keps_imposee_(i)>-1) compteur2++;
            }
          Cerr<<"compteurs "<<compteur2<<" "<<compteur<<finl;
          exit();
        }
    }
#endif

  // on recalcule partout ou c'est impose
  const int nb_faces_tot = domaine_VEF.nb_faces_tot();
  for (int face = 0; face < nb_faces_tot; face++)
    {
      const int num_face = face_keps_imposee_[face];
      if (num_face > -1)
        {
          // int elem_voisin;
          const int elem = face_voisins(num_face, 0);
          //if (face_voisins(face,0)!=elem) elem_voisin=face_voisins(face,0);
          //else elem_voisin=face_voisins(face,1);
          //int elem_voisin=face_voisins(num_face,0);
          // ce n'est pas le bon voisin!!!!
          double distb;
          if (!is_champ_Q1NC)
            {
              distb = distance_face(num_face, face, domaine_VEF);
            }
          else
            {
              assert(sub_type(Champ_Q1NC, eqn_hydr.inconnue())) ;
              distb = distance_face_elem(num_face, elem, domaine_VEF);
            }

          norm_v = norm_vit_lp_k(vitesse, face, num_face, domaine_VEF,
                                 val, is_defilante_face[num_face]);

          if (l_unif)
            d_visco = visco;
          else
            d_visco = tab_visco(elem, 0);

          u_plus_d_plus = norm_v*distb/d_visco;
          u_plus = calculer_u_plus(face, u_plus_d_plus, stock_erugu[num_face]);

          if((bool)u_plus)
            {
              u_star = norm_v/u_plus ;
              d_plus = u_plus_d_plus/u_plus ;
            }
          else
            {
              u_star = 0.;
              d_plus = 0.;
            }

          // calculer_k_eps(tab_k_eps(face,0),tab_k_eps(face,1),d_plus,u_star,d_visco,distb);
          // Evaluate the turbulent quantities
          compute_turbulent_quantities(tab_2eq(face, 0), tab_2eq(face, 1),
                                       d_plus, u_star, d_visco, distb);
        }
    }

  Cisaillement_paroi_.echange_espace_virtuel();
  tab_2eq.echange_espace_virtuel();
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd tab_2eq", tab_2eq);
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd Cisaillement_paroi_", Cisaillement_paroi_);
  return 1;
}  // fin de calcul_hyd (K-eps)

void Paroi_std_hyd_VEF::compute_turbulent_quantities(double& tke, double& var2,
                                                     const double d_plus, const double u_star,
                                                     const double d_visco, const double dist)
{
  if (turbulence_model_type == 1) // K-Epsilon model
    calculer_k_eps(tke, var2, d_plus, u_star, d_visco, dist);
  else if (turbulence_model_type == 2) // K-Omega model
    compute_k_omega(tke, var2, d_plus, u_star, d_visco, dist);
  else
    Cerr << "The turbulence model is neither k-epsilon nor k-omega. Implementation required?\n";
  return;
}

KOKKOS_FUNCTION  // A reverser dans VEF/Domaine (?)
double norm_vit_lp(int dim, const double* v, int fac, CDoubleTabView face_normale, double* val)
{
  double r[3] {};
  double norme = 0;
  for(int i = 0; i < dim; i++)
    {
      r[i] = face_normale(fac, i);
      norme += r[i] * r[i];
    }
  norme = sqrt(norme);
  for(int i = 0; i < dim; i++)
    r[i] /= norme;

  double sum_carre=0;
  double psc=0;
  for (int i=0; i<dim; i++)
    {
      sum_carre += carre(v[i]);
      psc += v[i] * r[i];
    }
  double norm_vit = sqrt(std::fabs(sum_carre-carre(psc)));

  // val1,val2 val3 sont les vitesses tangentielles
  for(int i = 0; i < dim; i++)
    val[i] = (v[i] - psc*r[i])/(norm_vit + DMINFLOAT);

  return norm_vit;
}

double norm_vit_lp(const ArrOfDouble& vit, int face, const Domaine_VEF& domaine, ArrOfDouble& val)
{
  // A reverser dans VEF/Domaine (?)

  const DoubleTab& face_normale = domaine.face_normales();
  const int dim = Objet_U::dimension;
  ArrOfDouble r(dim);

  for(int i = 0; i < dim; i++)
    r[i] = face_normale(face, i);

  r /= norme_array(r);
  const double psc = dotproduct_array(r, vit);

  double norm_vit {0};
  if (dim == 3)
    norm_vit = vitesse_tangentielle(vit[0],vit[1],vit[2],r[0],r[1],r[2]);
  else
    norm_vit = vitesse_tangentielle(vit[0],vit[1],r[0],r[1]);

  for(int i = 0; i < dim; i++)
    val[i] = (vit[i] - psc*r[i])/(norm_vit + DMINFLOAT);

  return norm_vit;
}

// cAlan: unused function?
int Paroi_std_hyd_VEF::calculer_hyd(DoubleTab& tab_nu_t,DoubleTab& tab_k)
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  const Domaine& domaine = domaine_VEF.domaine();
  int nfac = domaine.nb_faces_elem();

  double visco0=-1;
  int l_unif;
  if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
    {
      visco0 = std::max(tab_visco(0,0),DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;
  if ((!l_unif) && (tab_visco.local_min_vect()<DMINFLOAT))
    //   on ne doit pas changer tab_visco ici !
    {
      Cerr<<" visco <=0 ?"<<finl;
      exit();
    }
  //tab_visco+=DMINFLOAT;


  Cisaillement_paroi_=0;

  double dist_corr=1.;
  double coef_vit=nfac;
  if (sub_type(Champ_P1NC,eqn_hydr.inconnue()))
    {
      dist_corr=double(dimension+1)/double(dimension);
      coef_vit=nfac-1;
    }

  bool LM   = sub_type(Modele_turbulence_hyd_Longueur_Melange_base,mon_modele_turb_hyd.valeur()); // Longueur de Melange
  bool COMB = sub_type(Modele_turbulence_hyd_combinaison,mon_modele_turb_hyd.valeur());  //Modele Combinaison (fonction analytique et (ou) dependance a des champs sources)
  int nombre_sources = COMB ? ref_cast(Modele_turbulence_hyd_combinaison,mon_modele_turb_hyd.valeur()).nombre_sources() : -1;

  // Loop on boundaries
  int nb_bords=domaine_VEF.nb_front_Cl();
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);

      bool dirichlet = sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) || sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) || la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE";
      bool robin = sub_type(Paroi_decalee_Robin,la_cl.valeur());
      // Only Dirichlet or Robin conditions:
      if (dirichlet || robin)
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          // Loop on real faces
          int ndeb = 0;
          int nfin = le_bord.nb_faces_tot();
          int dim = Objet_U::dimension;
          double Cmu = Cmu_;
          double Kappa = Kappa_;
          int is_u_star_impose = is_u_star_impose_;
          double u_star_impose = u_star_impose_;
          const double delta = robin ? ref_cast(Paroi_decalee_Robin,la_cl.valeur()).get_delta() : 0;
          double erugu=Erugu;
          if (sub_type(Paroi_rugueuse,la_cl.valeur()))
            erugu=ref_cast(Paroi_rugueuse,la_cl.valeur()).get_erugu();
          // Lecture
          CIntArrView le_bord_num_face = le_bord.num_face().view_ro();
          CIntTabView face_voisins = domaine_VEF.face_voisins().view_ro();
          CIntTabView elem_faces = domaine_VEF.elem_faces().view_ro();
          CDoubleTabView vitesse = eqn_hydr.inconnue().valeurs().view_ro();
          CDoubleTabView face_normales = domaine_VEF.face_normales().view_ro();
          CDoubleTabView xp = domaine_VEF.xp().view_ro();
          CDoubleTabView xv = domaine_VEF.xv().view_ro();
          CDoubleTabView visco = tab_visco.view_ro();
          // Ecriture
          DoubleArrView seuil_LP = seuil_LP_.view_wo();
          IntArrView iterations_LP = iterations_LP_.view_wo();
          DoubleArrView uplus = uplus_.view_wo();
          DoubleArrView tab_d_plus = tab_d_plus_.view_wo();
          DoubleArrView tab_u_star = tab_u_star_.view_wo();
          DoubleArrView nu_t = static_cast<ArrOfDouble&>(tab_nu_t).view_wo();
          DoubleArrView k = static_cast<ArrOfDouble&>(tab_k).view_wo();
          DoubleTabView Cisaillement_paroi = Cisaillement_paroi_.view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(ndeb, nfin), KOKKOS_LAMBDA (const int ind_face)
          {
            double val[3] {};
            double vit[3] {};
            int num_face = le_bord_num_face(ind_face);
            int elem = face_voisins(num_face,0);

            double dist=0, norm_v=0;
            if (dirichlet)
              {
                // Dirichlet
                for (int j = 0; j < dim; j++)
                  vit[j] = 0.;
                for (int i = 0; i < nfac; i++)
                  {
                    int face = elem_faces(elem, i);
                    for (int j = 0; j < dim; j++)
                      vit[j] += (vitesse(face, j) -
                                 vitesse(num_face, j)); // permet de soustraire la vitesse de glissement eventuelle
                  }
                for (int j = 0; j < dim; j++)
                  vit[j] /= coef_vit;

                dist = distance(dim, num_face, elem, xp, xv, face_normales);
                dist *= dist_corr; // pour passer du centre de gravite au milieu des faces en P1NC

                norm_v = norm_vit_lp(dim, vit, num_face, face_normales, val);
              }
            else
              {
                // Robin
                double psc=0, norm=0;
                for(int comp=0; comp<dim; comp++)
                  {
                    psc += vitesse(num_face,comp)*face_normales(num_face,comp);
                    norm += face_normales(num_face,comp)*face_normales(num_face,comp);
                  }
                // psc /= norm; // Fixed bug: Arithmetic exception
                if (std::fabs(norm)>=DMINFLOAT) psc/=norm;

                for(int comp=0; comp<dim; comp++)
                  {
                    val[comp]=vitesse(num_face,comp)-psc*face_normales(num_face,comp);
                    norm_v += val[comp]*val[comp];
                  }
                norm_v = sqrt(norm_v);
                for (int j = 0; j < dim; j++)
                  val[j] /= norm_v;
                dist = delta;
              }

            double d_visco = l_unif ? visco0 : visco(elem,0);

            double u_plus_d_plus = norm_v*dist/d_visco;

            double u_plus = calculer_u_plus_kokkos(ind_face,u_plus_d_plus,erugu,Kappa,seuil_LP,iterations_LP);

            double u_star, d_plus;
            if (is_u_star_impose)
              {
                u_star = u_star_impose;
                d_plus = 0.;
              }
            else
              {
                u_star = u_plus ? norm_v/u_plus : 0.;
                d_plus = u_plus ? u_plus_d_plus/u_plus : 0.;
              }
            double eps;
            compute_k_eps(k(elem),eps,d_plus,u_star,d_visco,dist, Cmu, Kappa);

            // Calcul de la contrainte tangentielle
            for (int j=0; j<dim; j++)
              Cisaillement_paroi(num_face,j) = u_star*u_star*val[j];

            // Remplissage des tableaux (dans le cas de Longueur de melange on laisse la viscosite telle quelle)
            if((!LM) && (!COMB)) nu_t(elem) = Cmu*k(elem)*k(elem)/(eps+DMINFLOAT);

            uplus(num_face) = u_plus;
            tab_d_plus(num_face) = d_plus;
            tab_u_star(num_face) = u_star;

            // Modification de nu_t (et par consequent lambda_t) pour exploiter la valeur de nu_t (lambda_t) en y=deq_lam.
            // La valeur de dist_corr n est valable que dans le cas particuler ou nu_t est fonction lineaire de y
            if (COMB && nombre_sources==0)
              nu_t(elem) *= dist_corr;

          }); // End loop on real faces
          end_gpu_timer(Objet_U::computeOnDevice, __KERNEL_NAME__);

        } // End BC conditions

    } // End loop on boundaries

  Cisaillement_paroi_.echange_espace_virtuel();
  tab_nu_t.echange_espace_virtuel();
  tab_k.echange_espace_virtuel();
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd k",tab_k);
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd tab_nu_t",tab_nu_t);
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd Cisaillement_paroi_",Cisaillement_paroi_);
  return 1;
}  // fin du calcul_hyd (nu-t)

void Paroi_std_hyd_VEF::compute_k(double& k, const double yp, const double u_star)
{
  // Original comment:
  // PQ : 05/04/07 : formulation continue de k et epsilon
  //  assurant le bon comportement asymptotique

  k = 0.07*yp*yp*exp(-yp/9.);
  k += 1./sqrt(Cmu_)*(1. - exp(-yp/20.))*(1. - exp(-yp/20.));
  k *= u_star*u_star; // = k+ * u_star^2
}

void Paroi_std_hyd_VEF::compute_epsilon(double& eps, const double yp,
                                        const double u_star, const double d_visco)
{
  double const u_star_squared = u_star*u_star;

  // 50625 = 15^4
  eps  = 1./(Kappa_*pow(yp*yp*yp*yp + 50625, 0.25)); // eps_plus
  eps *= u_star_squared*u_star_squared/d_visco;
}

void Paroi_std_hyd_VEF::compute_omega(double& omega, const double yp,
                                      const double u_star, const double d_visco,
                                      const double dist)
{
  // cAlan: à mutualiser avec le calcul du k_omega de Multiphase
  const double w_vis = 6.*d_visco/(BETA_OMEGA*dist*dist);
  const double w_log = u_star/(std::sqrt(BETA_K)*Kappa_*dist);
  const double w_1 = w_vis + w_log ;
  const double w_2 = std::pow(std::pow(w_vis, 1.2) + std::pow(w_log, 1.2), 1./1.2);
  const double blending = std::tanh(yp/10*yp/10*yp/10*yp/10);

  omega = blending*w_1 + (1 - blending)*w_2;
}

int Paroi_std_hyd_VEF::calculer_k_eps(double& k, double& eps ,
                                      const double yp, const double u_star,
                                      const double d_visco, const double dist)
{
  return compute_k_eps(k, eps, yp, u_star, d_visco, dist, Cmu_, Kappa_);
}


void Paroi_std_hyd_VEF::compute_k_epsilon(double& k, double& epsilon,
                                          double yp, double u_star,
                                          double d_visco, double dist)
{
  compute_k(k, yp, u_star);
  compute_epsilon(epsilon, yp, u_star, d_visco);
  return;
}

void Paroi_std_hyd_VEF::compute_k_omega(double& k, double& omega,
                                        double yp, double u_star,
                                        double d_visco, double dist)
{
  compute_k(k, yp, u_star);
  compute_omega(omega, yp, u_star, d_visco, dist);
  return;
}

void Paroi_std_hyd_VEF::imprimer_ustar(Sortie& os) const
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  int ndeb,nfin;
  double upmoy,dpmoy,utaumoy;
  double seuil_moy,iter_moy;
  double norme_L2=0.;

  upmoy=0.;
  dpmoy=0.;
  utaumoy=0.;
  seuil_moy=0.;
  iter_moy=0.;

  int compt=0;

  EcrFicPartage Ustar;
  ouvrir_fichier_partage(Ustar,"Ustar");

  for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet_paroi_fixe,la_cl.valeur())) ||
           (sub_type(Dirichlet_paroi_defilante,la_cl.valeur())) ||
           (sub_type(Paroi_decalee_Robin,la_cl.valeur()) ) ||
           (la_cl->que_suis_je() == "Frontiere_ouverte_vitesse_imposee_ALE"))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          if(je_suis_maitre())
            {
              Ustar << finl;
              Ustar << "Bord " << le_bord.le_nom() << finl;
              if (dimension == 2)
                {
                  Ustar << "-------------------------------------------------------------------------------------------";
                  Ustar << "--------------------------------------------------------------------------------------------" << finl;
                  Ustar << "\tFace a\t\t\t\t|\t\t\t\t\t\t\t\t\t| TAU=Nu.Grad(Ut) [m2/s2]" << finl;
                  Ustar << "----------------------------------------|--------------------------------------------------";
                  Ustar << "---------------------|----------------------------------------------------------------------" << finl;
                  Ustar << "X\t\t| Y\t\t\t| u+\t\t\t| d+\t\t\t| u*\t\t\t| ||TAU||_2\t\t| |TAUx|\t\t| |TAUy|" << finl;
                  Ustar << "----------------|-----------------------|-----------------------|-----------------------|--";
                  Ustar << "---------------------|-----------------------|-----------------------|----------------------" << finl;
                }
              if (dimension == 3)
                {
                  Ustar << "-----------------------------------------------------------------------------------------------------------------";
                  Ustar << "-----------------------------------------------------------------------------------------------------------------" << finl;
                  Ustar << "\tFace a\t\t\t\t\t\t\t|\t\t\t\t\t\t\t\t\t| TAU=Nu.Grad(Ut) [m2/s2]" << finl;
                  Ustar << "----------------------------------------------------------------|------------------------------------------------";
                  Ustar << "-----------------------|-----------------------------------------------------------------------------------------" << finl;
                  Ustar << "X\t\t| Y\t\t\t| Z\t\t\t| u+\t\t\t| d+\t\t\t| u*\t\t\t| ||TAU||_2\t\t| |TAUx|\t\t| |TAUy|\t\t| |TAUz|" << finl;
                  Ustar << "----------------|-----------------------|-----------------------|-----------------------|-----------------------|";
                  Ustar << "-----------------------|-----------------------|-----------------------|-----------------------|-----------------" << finl;
                }
            }
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              double x=domaine_VEF.xv(num_face,0);
              double y=domaine_VEF.xv(num_face,1);
              norme_L2= Cisaillement_paroi_(num_face,0)*Cisaillement_paroi_(num_face,0) + Cisaillement_paroi_(num_face,1)*Cisaillement_paroi_(num_face,1);
              if (dimension == 2)
                Ustar << x << "\t| " << y;
              if (dimension == 3)
                {
                  double z=domaine_VEF.xv(num_face,2);
                  Ustar << x << "\t| " << y << "\t| " << z;
                  norme_L2+= Cisaillement_paroi_(num_face,2)*Cisaillement_paroi_(num_face,2);
                }
              norme_L2=sqrt(norme_L2);
              Ustar << "\t| " << uplus_(num_face) << "\t| " << tab_d_plus(num_face) << "\t| " << tab_u_star(num_face);
              Ustar << "\t| " << norme_L2 << "\t| " << Cisaillement_paroi_(num_face,0) << "\t| " << Cisaillement_paroi_(num_face,1) ;
              if (dimension == 3)
                Ustar << "\t| " << Cisaillement_paroi_(num_face,2) << finl;
              else
                Ustar << finl;

              // PQ : 03/03 : Calcul des valeurs moyennes (en supposant maillage regulier)

              upmoy +=uplus_(num_face);
              dpmoy +=tab_d_plus(num_face);
              utaumoy +=tab_u_star(num_face);
              seuil_moy += seuil_LP_(num_face);
              iter_moy += iterations_LP_(num_face);
              compt +=1;
            }
          Ustar.syncfile();
        }
    }
  /* Reduce 6 mp_sum to 1 by using mp_sum_for_each_item:
  upmoy = mp_sum(upmoy);
  dpmoy = mp_sum(dpmoy);
  utaumoy = mp_sum(utaumoy);
  seuil_moy = mp_sum(seuil_moy);
  iter_moy = mp_sum(iter_moy);
  compt = mp_sum(compt);
  */

  ArrOfDouble array(6);
  array[0]=upmoy;
  array[1]=dpmoy;
  array[2]=utaumoy;
  array[3]=seuil_moy;
  array[4]=iter_moy;
  array[5]=compt;
  mp_sum_for_each_item(array);
  upmoy=array[0];
  dpmoy=array[1];
  utaumoy=array[2];
  seuil_moy=array[3];
  iter_moy=array[4];
  compt=(int)array[5];
  if (je_suis_maitre())
    {
      if (compt)
        {
          Ustar << finl;
          Ustar << "-------------------------------------------------------------" << finl;
          Ustar << "Calcul des valeurs moyennes (en supposant maillage regulier):" << finl;
          Ustar << "<u+>= " << upmoy/compt << " <d+>= " << dpmoy/compt << " <u*>= " << utaumoy/compt << " seuil_LP= " << seuil_moy/compt << " iterations_LP= " << iter_moy/compt << finl;
        }
      Ustar << finl << finl;
    }
  Ustar.syncfile();
}
