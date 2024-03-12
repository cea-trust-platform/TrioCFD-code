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
// File:        EDO_Pression_th_VDF_Gaz_Parfait.cpp
// Directory:   $TRIO_U_ROOT/src/ThHyd/Quasi_Compressible/VDF
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#include <EDO_Pression_th_VDF_Gaz_Parfait.h>
#include <Fluide_Quasi_Compressible.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Schema_Temps_base.h>
#include <Loi_Etat_GP.h>
//#include <Les_Cl.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps.h>
#include <DoubleTrav.h>
#include <Probleme_base.h>
#include <Op_Diff_VDF_base.h>

Implemente_instanciable(EDO_Pression_th_VDF_Gaz_Parfait,"EDO_Pression_th_VDF_Gaz_Parfait",EDO_Pression_th_VDF);


/*! @brief Imprime la loi sur un flot de sortie.
 *
 * @param (Sortie& os) le flot de sortie pour l'impression
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& EDO_Pression_th_VDF_Gaz_Parfait::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

/*! @brief Lecture d'une EDO sur un flot d'entree.
 *
 * @param (Entree& is) le flot d'entree pour la lecture des parametres
 * @return (Entree&) le flot d'entree modifie
 */
Entree& EDO_Pression_th_VDF_Gaz_Parfait::readOn(Entree& is)
{
  return is;
}


/*! @brief Resoud l'EDO
 *
 * @param (double Pth_n) La pression a l'etape precedente
 * @return (double) La nouvelle valeur de la pression
 */
double EDO_Pression_th_VDF_Gaz_Parfait::resoudre(double Pth_n)
{
  int traitPth = le_fluide_->getTraitementPth();
  if ((traitPth!=2)) // ajout d'une condition pour que Pth puisse utiliser le mode constant.
    {
      const Loi_Etat_GP& loi_ = ref_cast(Loi_Etat_GP,le_fluide_->loi_etat().valeur());
      double Cp=loi_.Cp();
      double R= loi_.R();
      double gamma_moins_1=R/(Cp-R);
      DoubleTrav tab_W(le_fluide_->inco_chaleur().valeurs());
      DoubleTrav S(tab_W);
      const Probleme_base& pb=  le_fluide_->inco_chaleur().valeur().equation().probleme();
      pb.equation(1).sources().ajouter(S);


      // a priori deja calcule...
      const Op_Diff_VDF_base& op=  ref_cast(Op_Diff_VDF_base,pb.equation(1).operateur(0).l_op_base());
      // pour mettre a jour les flux bords
      const DoubleTab& tab_rho_np1=loi_.rho_np1();
      double Pth_n_old=Pth_n;
      double dt = le_fluide_->vitesse()->equation().schema_temps().pas_de_temps();


      // recuperation vitesse
      const DoubleTab& vitesse= pb.equation(0).inconnue().valeurs();
      const Domaine_VDF& domainevdf=ref_cast(Domaine_VDF,pb.equation(1).domaine_dis().valeur());
      const DoubleVect& surface=domainevdf.face_surfaces();
      double int_divu=0;
      const IntTab& face_voisins=domainevdf.face_voisins();
      int nb_faces_bord=domainevdf.nb_faces_bord();
      for(int face=0; face<nb_faces_bord; face++)
        {

          if (face_voisins(face,0)==-1)
            {
              // pas de voisin a gauche
              int_divu-=vitesse(face)*surface(face);
            }
          else
            {
              if (face_voisins(face,1)==-1)
                int_divu+=vitesse(face)*surface(face);
            }
        }

      for (int boucle=0; boucle<3; boucle++)
        {
          DoubleTrav T_etoile(tab_rho_np1);
          int size_rho=T_etoile.dimension_tot(0);
          for (int i=0; i<size_rho; i++) T_etoile(i)=Pth_n/R/tab_rho_np1(i);
          //op.calculer(T_etoile,tab_W);
          op.calculer_flux_bord(T_etoile);
          //op.ajouter(le_fluide_.valeur().temperature(),tab_W);
          const DoubleTab& flux_bord = op.flux_bords();
          int size=flux_bord.dimension(0);
          double bilan=0;
          for (int f=0; f<size; f++)
            {
              bilan+=flux_bord(f,0);
            }
// bilan=mp_sum(bilan);

          double bilan_bis=0;
          double vtot=0;
          if (pb.equation(1).sources().size()>0)
            {
              int elem, nb_elem=le_dom->nb_elem();
// double bilan_inu=0;
              for (elem=0 ; elem<nb_elem ; elem++)
                {
                  double v = le_dom->volumes(elem);
                  vtot+=v;
                  //  bilan_inu+=tab_W(elem);
                  bilan_bis=+S(elem);
                }
              vtot=mp_sum(vtot);
              if (!est_egal(vtot,le_dom->domaine().volume_total())) Process::exit();
            }
          bilan_bis+=bilan;
          bilan_bis-=int_divu*Pth_n*Cp/R;

          bilan_bis=mp_sum(bilan_bis);
          vtot=le_dom->domaine().volume_total();

          //bilan*=R/Pth_n;
//  Cerr<<" IIII "<< Pth_n ;

          // Pth_n=Pth_n+dt*(gamma_moins_1/vtot)*bilan;

          Pth_n=Pth_n_old/(1-gamma_moins_1/vtot*dt*bilan_bis/Pth_n);

          //Pth_n=Pth_n+Pth_n*gamma_moins_1/vtot*dt*bilan/R;
          //Cerr <<" " <<boucle <<" new "<< Pth_n <<" dpth/dt "<<   (Pth_n-Pth_n_old)/dt ;
          //Cerr <<" dt "<<dt<<" bilan "<< bilan<<" bilan_bis " << bilan_bis<< " dbilan "<<bilan_bis-bilan<<" int_divu" << int_divu<<finl;
        }
      double temps=le_fluide_->vitesse()->equation().schema_temps().temps_courant();
      Cout << temps << " "<<(Pth_n-Pth_n_old)/dt<< " "<<temps+dt<<" suivi dpth"<<finl;

    }

  return Pth_n;


#if 0

  int traitPth = le_fluide_->getTraitementPth();
  const DoubleTab& tempnp1 = le_fluide_->inco_chaleur().valeurs();       //actuel
  if (traitPth==0)
    {
      Cerr<<"on choisit le traitement "<<finl;
      int n_bord ;
      traitPth=1;
      for (n_bord=0; n_bord<le_dom->nb_front_Cl(); n_bord++)
        {
          const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(n_bord);
          if (sub_type(Neumann_sortie_libre, la_cl.valeur()))
            {
              traitPth=2;
              return Pth_n;
            }
        }
    }

  double dt = le_fluide_->vitesse()->equation().schema_temps().pas_de_temps();
  if (traitPth==2)   //Cerr << "trait cst" << finl;
    {
      return Pth_n;
    }
  else if (traitPth == 1)
    {
      /*
      //Cerr << "trait cons" << finl;
      double Sco=0,v;
      int elem, nb_elem=le_dom->nb_elem();
      for (elem=0 ; elem<nb_elem ; elem++) {
      v = le_dom->volumes(elem);
      //V += v;
      //rho_moy += tab_rho(elem) *v;
      Sco += v/tempnp1(elem);
      }

      double Pth = loi_.R()*getM0()/Sco;
      //        Cerr<<"ICI " <<(Pth-Pth_n)/dt<<finl;
      return Pth;
      }
      else
      {
      */
      const DoubleTab& tempn = le_fluide_->inco_chaleur().passe();
      double cn1=0,cn=0,v;
      int elem, nb_elem=le_dom->nb_elem();
      for (elem=0 ; elem<nb_elem ; elem++)
        {
          v = le_dom->volumes(elem);
          //V += v;
          //rho_moy += tab_rho(elem) *v;
          cn1 += v/tempnp1(elem);
          cn+=v/tempn(elem);
        }
      int n_bord ;

      double cm=0;
      const IntTab& face_voisins=le_dom->face_voisins();
      const DoubleVect& Surface=le_dom->face_surfaces();
      // ce n'est pas la bonne vitesse mais on essaye
      const IntVect& orientation=le_dom->orientation();
      for (n_bord=0; n_bord<le_dom->nb_front_Cl(); n_bord++)
        {
          const Cond_lim_base& la_cl = le_dom_Cl->les_conditions_limites(n_bord).valeur();
          if (sub_type(Dirichlet, la_cl))
            {
              const Dirichlet& diri=ref_cast(Dirichlet,la_cl);
              const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
              int ndeb = la_front_dis.num_premiere_face();
              int nfin = ndeb + la_front_dis.nb_faces();
              for (int num_face=ndeb; num_face<nfin; num_face++)
                {
                  int n0 = face_voisins(num_face, 0);
                  double S=Surface(num_face);
                  if (n0 == -1)
                    {
                      n0 = face_voisins(num_face, 1);
                      S*=-1;
                    }
                  cm+=S/tempnp1(n0)*diri.val_imp(num_face-ndeb,orientation(num_face));
                }
            }

          else if (sub_type(Neumann_sortie_libre, la_cl))
            {
              Cerr<<la_cl.que_suis_je()<<" est incompatible avec le traitement conservation_masse pour l'instant"<<finl;
              abort();
            }
        }

      //exit();
      //Cerr<<"ici "<<cm;
      double cnt=mp_sum(cn);
      double cn1t=mp_sum(cn1);
      double cmt=mp_sum(cm);
      double Pt= Pth_n*cnt/cn1t/(1+dt/cn1t*cmt);
      return Pt;
    }
  Cerr<<" Utilisez traitement conservation_masse ou constant " <<finl;
  abort();
  // traitement edo ???
  int n_bord;
  for (n_bord=0; n_bord<le_dom->nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(n_bord);
      if (sub_type(Neumann_sortie_libre, la_cl.valeur()))
        return Pth_n;

    }
  Cerr<<" Le codage ici est faux !!!!"<<finl;
  Cerr<<" Pour l'instant on bloque"<<finl;
  assert(0);
  exit();
  const DoubleTab& tab_vit = ref_cast(Navier_Stokes_std,le_fluide_->vitesse()->equation()).vitesse();
  const DoubleTab& tempn = le_fluide_->inco_chaleur().passe();        //passe
  const DoubleTab& tab_rho = le_fluide_->masse_volumique().valeurs();    //actuel


  Cerr<<"---EDO : Tnp1="<<tempnp1(0)<<" Tn="<<tempn(0)<<finl;

  int elem, nb_elem=le_dom->nb_elem();
  double V = 0; //mesure du domaine
  double F = 0; //integrale du gradient de U
  double S = 0; //second membre

  //double Pth1 = le_fluide_->pression_th1();


  double v;


  const IntTab& elem_faces = le_dom->elem_faces();
  DoubleTrav divU(tab_vit.dimension(0));
  ref_cast(Navier_Stokes_std,le_fluide_->vitesse()->equation()).operateur_divergence().calculer(tab_vit,divU);
  DoubleTrav gradT(tab_vit.dimension(0));
  DoubleTrav Tstar(tab_vit.dimension(0));
  for (elem=0 ; elem<nb_elem ; elem++)
    {
      Tstar(elem) = .5*(tempn(elem)+tempnp1(elem));
    }
  calculer_grad(Tstar,gradT);
  DoubleTab u_gradT(nb_elem);
  int f1,f2,i;
  for (elem=0 ; elem<nb_elem ; elem++)
    {
      u_gradT(elem)=0;
      for (i=0 ; i<dimension ; i++)
        {
          f1 = elem_faces(elem,i);
          f2 = elem_faces(elem,i+dimension);
          //u_gradT(elem) += .25* (gradT(f1)+gradT(f2)) * (tab_vit(f1)+tab_vit(f2));
          u_gradT(elem) += .5* (gradT(f1)*tab_vit(f1) + gradT(f2)*tab_vit(f2));
        }
    }

  for (elem=0 ; elem<nb_elem ; elem++)
    {
      v = le_dom->volumes(elem);
      V += v;

      S += v * tab_rho(elem) * ((tempnp1(elem)-tempn(elem))/dt + u_gradT(elem));

    }

  //  int nb_cond_lim = le_dom->nb_front_Cl();
  int ndeb, nfin, face;
  double norm;
  DoubleVect norme(dimension);
  for (n_bord=0; n_bord<le_dom->nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(n_bord);
      const Front_VF& frontiere_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
      ndeb = frontiere_dis.num_premiere_face();
      nfin = ndeb + frontiere_dis.nb_faces();
      //if (sub_type(Neumann_sortie_libre, la_cl.valeur()) || sub_type(Dirichlet_entree_fluide, la_cl.valeur())) {
      for (face=ndeb ; face<nfin ; face++)
        {
          elem = le_dom->face_voisins(face,0);
          norm = le_dom->face_surfaces(face);
          if (elem==-1)
            {
              elem = le_dom->face_voisins(face,1);
              norm *= -1;
            }
          //calcul de F=Som(div(U))
          F += tab_vit(face) * norm;
          //}
        }
    }
  const Loi_Etat_GP& loi_ = ref_cast(Loi_Etat_GP,le_fluide_->loi_etat().valeur());
  S *= loi_.R();

  //   Cout<<"Volume="<<V<<finl;
  //   Cout<<"F="<<F<<finl;
  //   Cout<<"S="<<S<<finl;
  //   Cout<<"Vit_CN="<<tab_vit<<finl;
  //   Cout<<"temp="<<temp<<finl;
  //   Cout<<"tempn="<<tempn<<finl;


  double Pth = (1./(V/dt+F/2.)) * ((V/dt-F/2.)*Pth_n + S);
  //double Pth = (1./(V)) * ((-F)*Pth_n + S);


  Cerr<<"Pression thermo recalculee = "<<Pth<<finl;
  //   Cout<<"Pression thermo recalculee = "<<Pth<<finl;

  //return Pth1;
  return Pth;

  /*


  // Resolution par Euler Implicite => plus stable que Crank Nicholson mais ordre 1

  int traitPth = le_fluide_->getTraitementPth();

  if (traitPth==2) { //Cerr << "trait cst" << finl;
  return Pth_n;
  }

  int n_bord ;
  const DoubleTab& tab_vit_CN = le_fluide_->vitesse().valeurs();
  const DoubleTab& tempnp1 = le_fluide_->inco_chaleur().valeurs();       //actuel
  const DoubleTab& tempn = le_fluide_->inco_chaleur().passe();        //passe
  const DoubleTab& tab_rho = le_fluide_->masse_volumique().valeurs();    //actuel
  const Loi_Etat_GP& loi_ = ref_cast(Loi_Etat_GP,le_fluide_->loi_etat().valeur());


  int elem, nb_elem=le_dom->nb_elem();
  double V = 0; //mesure du domaine
  double F = 0; //integrale du gradient de U
  double S = 0; //second membre
  double rho_moy = 0; //second membre

  //double Pth1 = le_fluide_->pression_th1();


  double dt = le_fluide_->vitesse()->equation().schema_temps().pas_de_temps();
  double v;

  const IntTab& elem_faces = le_dom->elem_faces();
  double Sco=0,Pth;


  if (traitPth == 0)
  { //Cerr << "trait edo" << finl;
  for (elem=0 ; elem<nb_elem ; elem++) {
  v = le_dom->volumes(elem);
  V += v;

  rho_moy += tab_rho(elem) *v;
  S += v * tab_rho(elem) * ((tempnp1(elem)-tempn(elem))/dt)/tempnp1(elem);// + u_gradT(elem));
  // Sco += v/tempnp1(elem);
  }

  //double dp = (getM0()-rho_moy)*loi_.R()/Sco;

  S/=rho_moy;

  Pth = Pth_n/(1- S*dt) ;
  //Cerr<<"Pth "<<Pth<<finl;
  //Cerr << "dt ODE = " << 1/S << finl;

  }
  else if (traitPth == 1)
  { //Cerr << "trait cons" << finl;
  for (elem=0 ; elem<nb_elem ; elem++) {
  v = le_dom->volumes(elem);
  V += v;

  rho_moy += tab_rho(elem) *v;
  Sco += v/tempnp1(elem);
  }

  Pth = loi_.R()*getM0()/Sco;
  //Cerr<<"Pth "<<Pth<<finl;
  }

  Cerr<<"Pression thermo recalculee = "<<Pth<<finl;


  return Pth;

  */



#endif

  Cerr<<"YB - 'EDO_Pression_th_VDF_GP' - Pression thermo = "<<Pth_n<<finl;
  return Pth_n;
}



