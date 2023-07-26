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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Remaillage_FT_IJK.cpp
// Directory : $IJK_ROOT/src/FT/front_tracking
//
/////////////////////////////////////////////////////////////////////////////

#include <Remaillage_FT_IJK.h>
#include <Param.h>
#include <Maillage_FT_Disc.h>
#include <Comm_Group.h>
#include <stat_counters.h>

Implemente_instanciable_sans_constructeur(Remaillage_FT_IJK,"Remaillage_FT_IJK",Remaillage_FT) ;

Remaillage_FT_IJK::Remaillage_FT_IJK()
{
  nb_iter_barycentrage_ = 0; // Par defaut, pas de barycentrage
  nb_iter_bary_volume_seul_ = 0; // Pas de correction de volume
  nb_iter_remaillage_ = 1;
  facteur_longueur_ideale_ = -1; // Pas de redecoupage ou suppression d'arretes.
  //                                Sinon, la valeur recommandee en 3D est 1.5.
  critere_arete_ = 0.35;   // Je crois que c'est sa valeur par defaut en FTD. Je n'arrive pas a le retrouver.
  equilateral_ = 1;        // Par defaut, sur maillage etire, on fait des triangles equilateraux avec
  //                         un facteur_longueur qui vaut 1 si l'arrete a la meme taille que la diagonale.
  //                         Comportement different du FTD classique de trio (equilateral = 0) qui regarde
  //                         l'orientation de la facette pour ajuster sa taille a l'element eulerien.

}

Sortie& Remaillage_FT_IJK::printOn(Sortie& os) const
{
  //  Objet_U::printOn(os);
  os << "{\n"
     << "     pas_remaillage " << dt_remaillage_ << "\n"
     << "     nb_iter_barycentrage " << nb_iter_barycentrage_ << "\n"
     << "     relax_barycentrage " << relax_barycentrage_ << "\n";
  os     << "     critere_arete " << critere_arete_ << "\n"
         << "     seuil_dvolume_residuel " << seuil_dvolume_residuel_ << "\n";
  os     << "     nb_iter_correction_volume " << nb_iter_bary_volume_seul_ << "\n"
         << "     nb_iter_remaillage " << nb_iter_remaillage_ << "\n"
         << "     facteur_longueur_ideale " << facteur_longueur_ideale_ << "\n";
  os << "     equilateral " << equilateral_ << "\n"
     << "     lissage_courbure_coeff " << lissage_courbure_coeff_ << "\n"
     << "     lissage_courbure_iterations_systematique " << lissage_courbure_iterations_systematique_ << "\n"
     << "     lissage_courbure_iterations_si_remaillage " << lissage_courbure_iterations_si_remaillage_ << "\n";
  os << "   }\n" ;
  return os;
}
// XD remaillage_ft_ijk interprete nul 1 not_set
Entree& Remaillage_FT_IJK::readOn(Entree& is)
{
  Param p(que_suis_je());
  p.ajouter("pas_remaillage", &dt_remaillage_); // XD_ADD_P floattant not_set
  p.ajouter("nb_iter_barycentrage", &nb_iter_barycentrage_); // XD_ADD_P entier not_set
  p.ajouter("relax_barycentrage", &relax_barycentrage_); // XD_ADD_P floattant not_set
  p.ajouter("critere_arete", &critere_arete_); // XD_ADD_P floattant not_set
  p.ajouter("seuil_dvolume_residuel", &seuil_dvolume_residuel_); // XD_ADD_P floattant not_set
  p.ajouter("nb_iter_correction_volume",  &nb_iter_bary_volume_seul_); // XD_ADD_P entier not_set
  p.ajouter("nb_iter_remaillage", &nb_iter_remaillage_); // XD_ADD_P entier not_set
  p.ajouter("facteur_longueur_ideale", &facteur_longueur_ideale_); // XD_ADD_P floattant not_set
  p.ajouter("equilateral", &equilateral_); // XD_ADD_P entier not_set
  p.ajouter("lissage_courbure_coeff", &lissage_courbure_coeff_); // XD_ADD_P floattant not_set
  p.ajouter("lissage_courbure_iterations_systematique", &lissage_courbure_iterations_systematique_); // XD_ADD_P entier not_set
  p.ajouter("lissage_courbure_iterations_si_remaillage", &lissage_courbure_iterations_si_remaillage_); // XD_ADD_P entier not_set
  p.lire_avec_accolades_depuis(is);
  Cout << "Remaillage_FT_IJK::readOn : Les options lues sont : " << finl;
  p.print(Cout);

  if ((dt_remaillage_ > 0. ) && (facteur_longueur_ideale_<0.))
    {
      Cerr << "Erreur dans les parametres de Remaillage_FT_IJK."
           << " Il faut specifier la valeur de facteur_longueur_ideale si pas_remaillage>0."
           << " La valeur recommandee est : facteur_longueur_ideale 1.5 (en 3D)." << finl;
      Process::exit();
    }
  return is;
}

void Remaillage_FT_IJK::barycentrer_lisser_systematique_ijk(Maillage_FT_Disc& maillage,
                                                            ArrOfDouble& var_volume)
{
  static Stat_Counter_Id barycentre_lissage_sys_counter_ = statistiques().new_counter(2, "remaillage interf: bary/lissage systematiques");
  statistiques().begin_count(barycentre_lissage_sys_counter_);
  regulariser_maillage(maillage,
                       var_volume,
                       relax_barycentrage_,
                       lissage_courbure_coeff_,
                       nb_iter_barycentrage_,
                       lissage_courbure_iterations_systematique_,
                       nb_iter_bary_volume_seul_,
                       seuil_dvolume_residuel_);
  supprimer_facettes_bord(maillage);
  statistiques().end_count(barycentre_lissage_sys_counter_);


}

void Remaillage_FT_IJK::barycentrer_lisser_apres_remaillage(Maillage_FT_Disc& maillage, ArrOfDouble& var_volume)
{
  static Stat_Counter_Id barycentre_lissage_apres_counter_ = statistiques().new_counter(2, "remaillage local : bary/lissage apres remaillage");
  statistiques().begin_count(barycentre_lissage_apres_counter_);
  regulariser_maillage(maillage, var_volume,
                       relax_barycentrage_,
                       lissage_courbure_coeff_,
                       nb_iter_barycentrage_,
                       lissage_courbure_iterations_si_remaillage_,
                       nb_iter_bary_volume_seul_,
                       seuil_dvolume_residuel_);
  supprimer_facettes_bord(maillage);
  // Dans le doute, je laisse l'appel a nettoyer_maillage :
  nettoyer_maillage(maillage);
  statistiques().end_count(barycentre_lissage_apres_counter_);

}

// Surcharge de Remaillage_FT::diviser_grandes_aretes(Maillage_FT_Disc& maillage) const
// A la creation des facettes, il faut leur attribuer un numero de compo connexe dans compo_connex_facettes_.
int Remaillage_FT_IJK::diviser_grandes_aretes(Maillage_FT_IJK& maillage) const
{
  static Stat_Counter_Id sup_div_aretes_counter_ = statistiques().new_counter(2, "remaillage local : suppressions / divisions aretes");
  statistiques().begin_count(sup_div_aretes_counter_);
  static int compteur = 0;
  static int test_val = -1;

  Process::Journal()<<"Remaillage_FT_IJK::diviser_grandes_aretes "<<temps_<<"  nb_som="<<maillage.nb_sommets()<<finl;
  Process::Journal()<<" Compteur = " << compteur << finl;
  compteur++;
  if (compteur == test_val)
    {
      Process::Journal() << " STOP." << finl;
    }
  //  int res = 1;

  maillage.nettoyer_elements_virtuels();

  //tableaux de stockage
  IntTabFT tab_aretesMarquees;
  ArrOfIntFT tab_somD;
  DoubleTabFT tab_deplacement_somD;

  //on commence par marquer les grandes aretes
  marquer_aretes(maillage,
                 tab_aretesMarquees,
                 tab_somD,
                 tab_deplacement_somD,
                 1 /* marquage des aretes trop grandes */);
  // resultat =  tab_aretesMarquees : [ fa7 iarete pe som ]

  const int nb_aretes_divis = tab_aretesMarquees.dimension(0);


  int nb_facettes = maillage.nb_facettes();
  const int nb_facettes0 = nb_facettes;
  IntTab& facettes = maillage.facettes_;
  const int nb_som_par_facette = facettes.dimension(1);

  const ArrOfInt& sommet_num_owner = maillage.sommet_num_owner_;
  int fa7,iarete, isom,som,isom_s,som_s,isom_ss,som_ss, pe_somD,numOwner_somD,somD,somD_s,somD_ss;

  const int dimension3 = (dimension==3);
  const int nb_aretes_par_facette = (dimension3)?3:1;
  //tableau stockant le sommet servant a decouper l'arete (ou -1 si arete non divisee)
  IntTabFT tab_fa7Divis(nb_facettes,nb_aretes_par_facette);
  for (fa7=0 ; fa7<nb_facettes ; fa7++)
    {
      for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
        {
          tab_fa7Divis(fa7,isom) = -1;
        }
    }

  //on va balayer les aretes a memoriser l'ensemble des aretes a scinder
  for (iarete=0 ; iarete<nb_aretes_divis ; iarete++)
    {
      fa7 = tab_aretesMarquees(iarete,0);
      isom = tab_aretesMarquees(iarete,1);
      isom_s = (isom+1)%nb_som_par_facette;
      //sommet qui va rester dans fa7
      som = facettes(fa7,isom);
      //sommet qui va aller dans une nouvelle facette
      som_s = facettes(fa7,isom_s);
      //sommet a inserer dans l'arete
      pe_somD = tab_aretesMarquees(iarete,2);
      numOwner_somD = tab_aretesMarquees(iarete,3);

      if (pe_somD!=me())
        {
          //je ne connais pas l'indice du sommet a inserer dans me()
          maillage.convertir_numero_distant_local(maillage.desc_sommets(),sommet_num_owner,numOwner_somD,pe_somD,somD);
          assert(somD >= 0);
        }
      else
        {
          //je suis le proprietaire de somD : je connais donc le bon indice du sommet a inserer
          somD = numOwner_somD;
        }

      tab_fa7Divis(fa7,isom) = somD;
    }

  // Specifique IJK :
  ArrOfInt& compo_connexe_facettes = maillage.compo_connexe_facettes_non_const();

  //on va ensuite balayer les facettes et les scinder
  //la configuration depend du nb d'aretes a scinder par facette
  int nb_areteScinder, isom0=-1,isom1=-1;
  for (fa7=0 ; fa7<nb_facettes0 ; fa7++)
    {
      //on compte le nombre d'aretes a scinder
      nb_areteScinder = 0;
      for (isom=0 ; isom<nb_aretes_par_facette ; isom++)
        {
          if (tab_fa7Divis(fa7,isom)>=0)
            {
              if (nb_areteScinder==0)
                {
                  isom0 = isom;
                }
              else
                {
                  isom1 = isom;
                }
              nb_areteScinder++;
            }
        }
      if (nb_areteScinder==1)
        {
          //s'il n'y a qu'une arete a scinder
          somD = tab_fa7Divis(fa7,isom0);
          //on modifie la facettes d'origine
          isom_s = (isom0+1)%nb_som_par_facette;
          som_s = facettes(fa7,isom_s);
          facettes(fa7,isom_s) = somD;
          //on cree une nouvelle facette
          if (nb_facettes>=facettes.dimension(0))
            {
              //Redimensionnement
              facettes.resize(nb_facettes+10,facettes.dimension(1));
              compo_connexe_facettes.resize_array(nb_facettes+10);
            }
          isom_ss = (isom_s+1)%nb_som_par_facette;
          som_ss = facettes(fa7,isom_ss);
          compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
          if (dimension==2)
            {
              facettes(nb_facettes,isom0) = somD;
              facettes(nb_facettes,isom_s) = som_s;
            }
          else
            {
              facettes(nb_facettes,0) = som_ss;
              facettes(nb_facettes,1) = somD; //sommet insere
              facettes(nb_facettes,2) = som_s;
            }
          nb_facettes++;
        }
      else if (nb_areteScinder==2)
        {
          //si 2 aretes sont a scinder
          //on cree deux nouvelles facettes
          if (nb_facettes+2>=facettes.dimension(0))
            {
              //Redimensionnement
              facettes.resize(nb_facettes+12,facettes.dimension(1));
              compo_connexe_facettes.resize_array(nb_facettes+12);
            }
          //on positionne isom et isom_s tq elles soient les aretes a scinder
          if (isom1==(isom0+1)%nb_som_par_facette)
            {
              isom = isom0;
              isom_s = isom1;
            }
          else
            {
              isom = isom1;
              isom_s = isom0;
            }
          isom_ss = ((isom_s+1)%nb_som_par_facette);
          som = facettes(fa7,isom);
          som_s = facettes(fa7,isom_s);
          som_ss = facettes(fa7,isom_ss);
          somD = tab_fa7Divis(fa7,isom);
          somD_s = tab_fa7Divis(fa7,isom_s);
          somD_ss = tab_fa7Divis(fa7,isom_ss);
          //on modifie la facette existante
          facettes(fa7,0) = som;
          facettes(fa7,1) = somD;
          facettes(fa7,2) = som_ss;
          //premiere nouvelle facette
          compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = som_s;
          facettes(nb_facettes,2) = somD_s;
          nb_facettes++;
          //seconde nouvelle facette
          compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = somD_s;
          facettes(nb_facettes,2) = som_ss;
          nb_facettes++;
        }
      else if (nb_areteScinder==3)
        {
          //si toutes les aretes sont a scinder
          //on cree trois nouvelles facettes
          if (nb_facettes+3>=facettes.dimension(0))
            {
              //Redimensionnement
              facettes.resize(nb_facettes+13,facettes.dimension(1));
              compo_connexe_facettes.resize_array(nb_facettes+13);
            }
          isom = 0;
          isom_s = 1;
          isom_ss = 2;
          som = facettes(fa7,isom);
          som_s = facettes(fa7,isom_s);
          som_ss = facettes(fa7,isom_ss);
          somD = tab_fa7Divis(fa7,isom);
          somD_s = tab_fa7Divis(fa7,isom_s);
          somD_ss = tab_fa7Divis(fa7,isom_ss);
          //on modifie la facette existante
          facettes(fa7,0) = som;
          facettes(fa7,1) = somD;
          facettes(fa7,2) = somD_ss;
          //premiere nouvelle facette
          compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = som_s;
          facettes(nb_facettes,2) = somD_s;
          nb_facettes++;
          //seconde nouvelle facette
          compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
          facettes(nb_facettes,0) = somD_s;
          facettes(nb_facettes,1) = som_ss;
          facettes(nb_facettes,2) = somD_ss;
          nb_facettes++;
          //troisieme nouvelle facette
          compo_connexe_facettes[nb_facettes] = compo_connexe_facettes[fa7];
          facettes(nb_facettes,0) = somD;
          facettes(nb_facettes,1) = somD_s;
          facettes(nb_facettes,2) = somD_ss;
          nb_facettes++;
        }
    }
  //Redimensionnement
  facettes.resize(nb_facettes,facettes.dimension(1));
  compo_connexe_facettes.resize_array(nb_facettes);
  maillage.desc_facettes_.calcul_schema_comm(nb_facettes);
  maillage.corriger_proprietaires_facettes();

  ArrOfIntFT liste_sommets_sortis;
  ArrOfIntFT numero_face_sortie;
  maillage.deplacer_sommets(tab_somD,tab_deplacement_somD,liste_sommets_sortis,numero_face_sortie);
  maillage.corriger_proprietaires_facettes();

  int nb_aretes_divis_tot = Process::mp_sum(nb_aretes_divis);
  Process::Journal()<<"FIN Remaillage_FT::diviser_grandes_aretes "<<temps_<<"  nb_som="<<maillage.nb_sommets()
                    <<"  nb_aretes_divisees_on_proc="<< nb_aretes_divis;
  Process::Journal()<< "  nb_aretes_divisees_tot="<< nb_aretes_divis_tot<<finl;
  maillage.maillage_modifie(Maillage_FT_Disc::MINIMAL);

  statistiques().end_count(sup_div_aretes_counter_);

  return nb_aretes_divis_tot;
  //  return res;
}

void Remaillage_FT_IJK::remaillage_local_interface(double temps, Maillage_FT_IJK& maillage)
{
  static Stat_Counter_Id remaillage_loc_interf_counter_ = statistiques().new_counter(2, "remaillage interf: remaillage local");
  statistiques().begin_count(remaillage_loc_interf_counter_);
  temps_dernier_remaillage_ = temps_dernier_lissage_ = temps_ = temps;

  maillage.nettoyer_elements_virtuels();
  maillage.check_mesh();
  //boucle sur les remaillages
  int iter;
  ArrOfDoubleFT varVolume;

  for (iter = 0; iter < nb_iter_remaillage_; iter++)
    {
      const int nb_sommets = maillage.nb_sommets();
      varVolume.resize_array(nb_sommets);
      varVolume = 0.;
      variation_volume_ = 0.;
      // n = nombre de sommets supprimes
      const int n = supprimer_petites_aretes(maillage,varVolume);

      if (Comm_Group::check_enabled()) maillage.check_mesh();
      if (n > 0)
        {
          supprimer_doublons_facettes(maillage);
          if (Comm_Group::check_enabled()) maillage.check_mesh();
          // On a supprime les petites aretes en deplacant un noeud sur
          //  un autre, la variation de volume engendree a ete mise dans varVolume:
          // On barycentre et on lisse, notamment pour recuperer cette variation.
          barycentrer_lisser_apres_remaillage(maillage, varVolume);
          if (Comm_Group::check_enabled()) maillage.check_mesh();
        }


      // m est le nombre d'aretes divisees.
      const int m = diviser_grandes_aretes(maillage);

      if (Comm_Group::check_enabled()) maillage.check_mesh();
      if (m > 0)
        {
          varVolume.resize_array(maillage.nb_sommets());
          varVolume = 0.;
          barycentrer_lisser_apres_remaillage(maillage, varVolume);
          if (Comm_Group::check_enabled()) maillage.check_mesh();
        }
      if (Process::je_suis_maitre())
        Journal() << "remaillage_local_interface t= " << temps << " suppressions: " << n << " divisions: " << m << finl;

    }

  nettoyer_maillage(maillage);
  statistiques().end_count(remaillage_loc_interf_counter_);

}

Vecteur3 Remaillage_FT_IJK::get_delta_euler(const Maillage_FT_IJK& maillage) const
{
  //  const IJK_Splitting & s = maillage.ref_splitting().valeur();
  const IJK_Grid_Geometry& geom = maillage.ref_splitting().valeur().get_grid_geometry();
  //  const IJK_Grid_Geometry & geom = s.get_grid_geometry();
  Vecteur3 delta(0., 0., 0.);
  const int dim = Objet_U::dimension;
  for (int k = 0; k < dim; k++)
    delta[k] = geom.get_constant_delta(k);
  // Le carre de la diagonale des elements :
  //  diago_elem2_ = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
  return delta;
}

double Remaillage_FT_IJK::calculer_longueurIdeale2_arete(const Maillage_FT_Disc& maillage,
                                                         int som0,
                                                         double x, double y, double z) const
{
  double lgrId2 = 0.;
  if (facteur_longueur_ideale_ > 0.)
    {
      const Maillage_FT_IJK& maillage_ijk = ref_cast(Maillage_FT_IJK, maillage);
      const Vecteur3 delta_xv = get_delta_euler(maillage_ijk);
      const int dim = Objet_U::dimension;
      int k;
      if (equilateral_)
        {
          // On calcul la diagonale de l'element eulerien. Puis longueur au carre.
          for (k = 0; k < dim; k++)
            {
              lgrId2 += delta_xv[k]*delta_xv[k];
            }
        }
      else
        {
          const DoubleTab& sommets = maillage.sommets();
          const Vecteur3 xyz(x, y, z);
          Vecteur3 v(0., 0., 0.);
          double norme2 = 0.;
          for (k = 0; k < dim; k++)
            {
              v[k] = xyz[k] - sommets(som0, k);
              norme2 += v[k] * v[k];
            }
          if (norme2 == 0)
            {
              v[0] = 1.;
              v[1] = 1.;
              v[2] = dim==3 ? 1. : 0.;
              norme2 = dim;
            }
          double f = 1. / sqrt(norme2);
          norme2 = 0.;
          for (k = 0; k < dim; k++)
            {
              v[k] *= f * delta_xv[k];
              norme2 += v[k] * v[k];
            }
          lgrId2 = norme2;

        }
      lgrId2 *= facteur_longueur_ideale_ * facteur_longueur_ideale_;
    }
  else
    {
      Cerr << "Erreur Remaillage_FT_IJK::calculer_longueurIdeale2_arete" << finl;
      exit();
    }
  return lgrId2;
}
