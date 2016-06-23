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
// File:        Eq_rayo_semi_transp_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VDF
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#include <Eq_rayo_semi_transp_VDF.h>
#include <Fluide_Incompressible.h>
#include <Modele_rayo_semi_transp.h>
#include <Operateur_Diff_base.h>
#include <Champ_Uniforme.h>
#include <Champ_front_uniforme.h>
#include <Echange_contact_rayo_semi_transp_VDF.h>
#include <Echange_externe_impose_rayo_semi_transp.h>
#include <Echange_global_impose_rayo_semi_transp.h>
#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Flux_radiatif_VDF.h>
#include <Neumann_paroi_rayo_semi_transp_VDF.h>
#include <Symetrie.h>
#include <Ref_Champ_front.h>
#include <Zone_VDF.h>
#include <DoubleTrav.h>

Implemente_instanciable(Eq_rayo_semi_transp_VDF,"Eq_rayo_semi_transp_VDF",Equation_rayonnement_base);

// Description:
//    Imprime le type de l'equation sur un flot de sortie.
// Precondition:
// Parametre: Sortie& s
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Eq_rayo_semi_transp_VDF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}



// Description:
// Precondition:
// Parametre: Entree& s
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Eq_rayo_semi_transp_VDF::readOn(Entree& s )
{
  return Equation_rayonnement_base::readOn(s);
}



// Description:
//    Calcule le champ de l'irradiance en fonction des CLs
//    au temps donne en parametre et de la temperature au temps par defaut.
//    Met le resultat dans inconnue.valeurs()
//    Suppose que le schema en temps est Euler explicite...
// Precondition:
// Parametre:
// Retour:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
void Eq_rayo_semi_transp_VDF::resoudre(double temps)
{
  //  Cerr<<"Eq_rayo_semi_transp_VDF::resoudre : Debut"<<finl;
  const Zone_VF& zone_VF = ref_cast(Zone_VF, zone_dis().valeur());
  int nb_elem = zone_VF.nb_elem();
  const DoubleTab& kappa = fluide().kappa().valeurs();
  //calcul du second membre
  DoubleTrav secmem(inconnue().valeurs());
  Probleme_base& pb = Modele().probleme();
  double n,k;

  assert(pb.equation(1).inconnue().le_nom()=="temperature");
  const DoubleTab& temper = pb.equation(1).inconnue().valeurs();
  const DoubleTab& indice = fluide().indice().valeurs();
  const double& sigma = Modele().valeur_sigma();

  secmem = 0;
  int elem;
  for (elem=0; elem<nb_elem; elem++)
    {
      assert(fluide().indice().nb_comp());
      if(sub_type(Champ_Uniforme,fluide().indice().valeur()))
        n = indice(0,0);
      else
        n = indice(elem,0);

      assert(fluide().kappa().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide().kappa().valeur()))
        k = kappa(0,0);
      else
        k = kappa(elem,0);

      double vol = zone_VF.volumes(elem);
      double T = temper(elem);
      secmem(elem) += +4*n*n*sigma*pow(T,4)*k*vol;
      //  Cerr<<"T = "<<T<<", n = "<<n<<", sigma = "<<sigma<<", vol = "<<vol<<", secmem(elem) = "<<secmem(elem)<<finl;
    }

  // On met a jour les champs associes aux conditions aux limites
  // avant d'evaluer leur contribution dans la matrice de discretisation
  evaluer_cl_rayonnement(temps);
  terme_diffusif.valeur().contribuer_au_second_membre(secmem);

  if (solveur.valeur().que_suis_je() == "Solv_GCP")
    if (sub_type(Champ_Uniforme,fluide().kappa().valeur()))
      {
        Matrice matrice_tmp;
        dimensionner_Mat_Bloc_Morse_Sym(matrice_tmp);
        Mat_Morse_to_Mat_Bloc(matrice_tmp);
        solveur.resoudre_systeme(matrice_tmp.valeur(),secmem,irradiance_.valeurs());
      }
    else
      {
        Cerr<<finl;
        Cerr<<"Erreur dans Eq_rayo_semi_transp_VDF::resoudre()"<<finl;
        Cerr<<"Attention, on ne peut pas resoudre l'equation"<<finl;
        Cerr<<"de rayonnement semi transparent avec le solveur : "<<solveur.que_suis_je()<<finl;
        Cerr<<"car kappa n'est pas constant, donc, la_matrice"<<finl;
        Cerr<<"n'est pas symetrique"<<finl;
        exit();
      }
  else if (solveur.valeur().que_suis_je() == "Solv_Gmres")
    if (Process::nproc() == 1)
      solveur.resoudre_systeme(la_matrice,secmem,irradiance_.valeurs());
    else
      {
        /*        Matrice matrice_tmp;
                  dimensionner_Mat_Bloc_Morse(matrice_tmp);
                  Mat_Morse_to_Mat_Bloc(matrice_tmp);
                  solveur.resoudre_systeme(matrice_tmp.valeur(),secmem,irradiance_.valeurs(),irradiance_.valeur());*/
        Cerr<<finl;
        Cerr<<"Erreur dans Eq_rayo_semi_transp_VDF::resoudre()"<<finl;
        Cerr<<"Attention, on ne peut pas resoudre un probleme"<<finl;
        Cerr<<"de rayonnement semi transparent en parallele en"<<finl;
        Cerr<<"utilisant le solveur Solv_Gmres. Si vous traitez"<<finl;
        Cerr<<"un probleme avec kappa constant, vous contournerez"<<finl;
        Cerr<<"cette limitation en utilisant le solveur GCP avec"<<finl;
        Cerr<<"un preconditionnement GCP"<<finl;
        exit();
      }
  else
    {
      Cerr<<finl;
      Cerr<<"Erreur dans Eq_rayo_semi_transp_VDF::resoudre()"<<finl;
      Cerr<<"Attention, on ne peut pas utiliser le solveur : "<<solveur.que_suis_je()<<finl;
      Cerr<<"pour resoudre l'equation de rayonnement dans un "<<finl;
      Cerr<<"probleme de rayonnement semi transparent"<<finl;
      exit();
    }
}


void Eq_rayo_semi_transp_VDF::evaluer_cl_rayonnement(double temps)
{
  //  Cerr<<"Eq_rayo_semi_transp_VDF::evaluer_cl_rayonnement : Debut"<<finl;
  // Boucle sur les conditions aux limites de l'equation de rayonnement
  Conds_lim& les_cl_rayo = zone_Cl_dis().les_conditions_limites();

  // recherche des conditions aux limites associes au l'equation de temperature
  Equation_base& eq_temp = Modele().probleme().equation(1);
  assert(eq_temp.inconnue().le_nom()=="temperature");

  Conds_lim& les_cl_temp = eq_temp.zone_Cl_dis().les_conditions_limites();
  int num_cl_rayo=0;
  for(num_cl_rayo = 0; num_cl_rayo<les_cl_rayo.size(); num_cl_rayo++)
    {
      Cond_lim& la_cl_rayo =  zone_Cl_dis().les_conditions_limites(num_cl_rayo);
      if(sub_type(Flux_radiatif_VDF,la_cl_rayo.valeur()))
        {
          Flux_radiatif_VDF& la_cl_rayon = ref_cast(Flux_radiatif_VDF,la_cl_rayo.valeur());
          // Recherche des temperatures de bord pour cette frontiere
          Nom nom_cl_rayo = la_cl_rayo.frontiere_dis().le_nom();
          int num_cl_temp = 0;

          REF(Champ_front) Tb;
          int test_remplissage_Tb = 0;
          for(num_cl_temp = 0; num_cl_temp<les_cl_temp.size(); num_cl_temp++)
            {
              Cond_lim& la_cl_temp =  eq_temp.zone_Cl_dis().les_conditions_limites(num_cl_temp);
              Nom nom_cl_temp = la_cl_temp.frontiere_dis().le_nom();
              if(nom_cl_temp == nom_cl_rayo)
                {
                  if (sub_type(Neumann_paroi_rayo_semi_transp_VDF,la_cl_temp.valeur()))
                    {
                      Neumann_paroi_rayo_semi_transp_VDF& la_cl_temper
                        = ref_cast(Neumann_paroi_rayo_semi_transp_VDF,la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Echange_contact_rayo_semi_transp_VDF,la_cl_temp.valeur()))
                    {
                      Echange_contact_rayo_semi_transp_VDF& la_cl_temper
                        = ref_cast(Echange_contact_rayo_semi_transp_VDF,la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Echange_externe_impose_rayo_semi_transp,la_cl_temp.valeur()))
                    {
                      Echange_externe_impose_rayo_semi_transp& la_cl_temper
                        = ref_cast(Echange_externe_impose_rayo_semi_transp,la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl_temp.valeur()))
                    {
                      Frontiere_ouverte_temperature_imposee_rayo_semi_transp& la_cl_temper
                        = ref_cast(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Frontiere_ouverte_rayo_semi_transp,la_cl_temp.valeur()))
                    {
                      Frontiere_ouverte_rayo_semi_transp& la_cl_temper
                        = ref_cast(Frontiere_ouverte_rayo_semi_transp,la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Echange_global_impose_rayo_semi_transp,la_cl_temp.valeur()))
                    {
                      Echange_global_impose_rayo_semi_transp& la_cl_temper
                        = ref_cast(Echange_global_impose_rayo_semi_transp,la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else
                    {
                      Cerr <<"Coder pour les autres condition limites de l'equation de temperature 2 "<<finl;
                      Cerr<<"le cas d'une CL thermique  "<<la_cl_temp->que_suis_je()<<" n'est pas code"<<finl;
                      exit();
                    }
                }
            }
          if ( test_remplissage_Tb == 0)
            // On n'a pas remplie le tableau des temperatures de bord !!!!
            Cerr<<"On n'a pas remplie le tableau des temperatures de bord !!!!"<<finl;
          const Zone_VF& zvf = ref_cast(Zone_VF,zone_dis().valeur());
          la_cl_rayon.evaluer_cl_rayonnement(Tb.valeur(), fluide().kappa(), fluide().longueur_rayo(),
                                             fluide().indice(),zvf,Modele().valeur_sigma(),temps);
        }
      else if (sub_type(Symetrie,la_cl_rayo.valeur()))
        {
          ;
        }
      else
        {
          Cerr<<"La condition a la limite "<<la_cl_rayo.valeur().que_suis_je()<<" n'est pas connue"<<finl;
          Cerr<<" pour l'equation de rayonnement "<< finl;
          exit();
        }
    }
  //  Cerr<<"Eq_rayo_semi_transp_VDF::evaluer_cl_rayonnement : Fin"<<finl;
}


// Description:
//        modifie la matrice pour prendre en compte la presence de faces
//        rayonnantes au voisinage des elements de bord
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Eq_rayo_semi_transp_VDF::modifier_matrice()
{
  //  Cerr<<"Eq_rayo_semi_transp_VDF::modifier_matrice : Debut"<<finl;
  // On fait une boucle sur les conditions aux limites associees a l'equations
  Conds_lim& les_cl = zone_Cl_dis().les_conditions_limites();
  const Zone_VDF& zvdf = ref_cast(Zone_VDF,zone_dis().valeur());
  const IntTab& face_voisins=zvdf.face_voisins();
  const DoubleVect& face_surfaces = zvdf.face_surfaces();

  int num_cl=0;
  for(num_cl = 0; num_cl<les_cl.size(); num_cl++)
    {
      Cond_lim& la_cl = zone_Cl_dis().les_conditions_limites(num_cl);
      if (sub_type(Flux_radiatif_VDF,la_cl.valeur()))
        {
          Flux_radiatif_VDF& cl_radiatif = ref_cast(Flux_radiatif_VDF,la_cl.valeur());
          const DoubleTab& epsilon = cl_radiatif.emissivite().valeurs();
          const DoubleTab& long_rayo = fluide().longueur_rayo().valeurs();
          const DoubleTab& kappa = fluide().kappa().valeurs();
          double A = cl_radiatif.A();

          if (sub_type(Front_VF,la_cl.frontiere_dis()))
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              // Boucle sur les faces de le_bord
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              int face = ndeb;
              for(face = ndeb; face<nfin; face++)
                {
                  int elem = face_voisins(face,0);
                  if(elem == -1)
                    elem = face_voisins(face,1);

                  double eF = zvdf.dist_norm_bord(face);
                  double k,l_r;
                  assert(fluide().longueur_rayo().nb_comp() == 1);
                  assert(fluide().kappa().nb_comp() == 1);
                  if(sub_type(Champ_Uniforme,fluide().kappa().valeur()))
                    {
                      k = kappa(0,0);
                      l_r = long_rayo(0,0);
                    }
                  else
                    {
                      k = kappa(elem,0);
                      l_r = long_rayo(elem,0);
                    }
                  double epsi;
                  assert(cl_radiatif.emissivite().nb_comp() == 1);
                  if(sub_type(Champ_front_uniforme,cl_radiatif.emissivite().valeur()))
                    epsi = epsilon(0,0);
                  else
                    epsi = epsilon(face-ndeb,0);


                  //                  Cerr<<"face_surfaces(face) = "<<face_surfaces(face)<<finl;
                  double numer_coeff = l_r*face_surfaces(face);
                  double denum_coeff = 3*k*epsi;
                  denum_coeff = 1/denum_coeff;
                  denum_coeff *= A*(2-epsi);
                  denum_coeff = denum_coeff + eF;

                  double coeff = numer_coeff/denum_coeff;

                  // On rajoute ce coefficient sur la diagonale de la matrice de discretisation
                  if(epsi<DMINFLOAT)
                    ;
                  else
                    la_matrice(elem,elem) +=  coeff;
                }
            }
          else
            {
              Cerr<<"Erreur dans Eq_rayo_semi_transp_VDF::modifier_matrice()"<<finl;
              Cerr<<"la frontiere associee a la_cl ne derive pas de Front_VF"<<finl;
              exit();
            }
        }
      else if (sub_type(Symetrie,la_cl.valeur()))
        {
          ;
        }
      else
        {
          Cerr<<"La condition a la limite utilisee n'est pas connue pour l'equation"<<finl;
          Cerr<<" de rayonnement "<< finl;
          exit();
        }
    }
  //  Cerr<<"Eq_rayo_semi_transp_VDF::modifier_matrice : Fin"<<finl;
}


void Eq_rayo_semi_transp_VDF::completer()
{
  const Zone_dis_base& une_zone_dis = zone_dis().valeur();
  int n = une_zone_dis.nb_front_Cl();

  int ii;
  for (ii =0; ii<n; ii++)
    {
      const Frontiere_dis_base& la_fr_dis = une_zone_dis.frontiere_dis(ii);
      la_zone_Cl_dis.valeur().les_conditions_limites(ii)->associer_fr_dis_base(la_fr_dis);
    }

  Equation_rayonnement_base::completer();

  // On assemble la matrice une fois pour toute au debut du calcul
  // Attention, ceci n'est valable que si kappa est constant au cours
  // du temps
  assembler_matrice();
}

void Eq_rayo_semi_transp_VDF::assembler_matrice()
{
  const Zone_VF& zone_VF = ref_cast(Zone_VF, zone_dis().valeur());
  int nb_elem_tot = zone_VF.nb_elem_tot();

  const DoubleTab& irradi = irradiance_.valeurs();

  int i;

  la_matrice.clean();

  // Prise en compte de la partie div((1/3K)grad(irradiance)) dans la matrice
  // de discretisation.
  terme_diffusif.valeur().contribuer_a_avec(irradi,la_matrice);

  // Modification de la matrice pour prendre en compte le second membre en K*irradiance
  const DoubleTab& kappa = fluide().kappa().valeurs();

  Cerr<<"On verifie lors du calcul de la matrice de discretisation que "<<finl;
  Cerr<<"l'ordre de la matrice est bien egale au nombre d'elements"<<finl;
  if(la_matrice.ordre() != nb_elem_tot)
    exit();
  Cerr<<"Ordre de la matrice OK"<<finl;

  //  int nb_comp = fluide().kappa().nb_comp();
  double k;
  for (i=0; i<la_matrice.ordre(); i++)
    {
      assert(fluide().kappa().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide().kappa().valeur()))
        k = kappa(0,0);
      else
        k = kappa(i,0);

      double vol = zone_VF.volumes(i);

      la_matrice(i,i) = la_matrice(i,i) + k*vol;
    }

  // On modifie la matrice pour prendre en compte l'effet des parois
  // rayonnantes sur les elements de bord
  modifier_matrice();
}

void Eq_rayo_semi_transp_VDF::typer_op_grad()
{
  ;
}

/*
  void Eq_rayo_semi_transp_VDF::dimensionner_Mat_Bloc_Morse_Sym(Matrice& matrice_tmp)
  {
  const Zone_VDF& zone_VDF = ref_cast(Zone_VDF,zone_dis().valeur());
  const IntTab& face_voisins = zone_VDF.face_voisins();

  int elem1,elem2;
  int num_face,face;
  int n1 = nb_colonnes_tot();
  int n2 = nb_colonnes();

  matrice_tmp.typer("Matrice_Bloc");
  Matrice_Bloc& matrice=ref_cast(Matrice_Bloc, matrice_tmp.valeur());
  matrice.dimensionner(1,2);
  matrice(0,0).typer("Matrice_Morse_Sym");
  matrice(0,1).typer("Matrice_Morse");

  Matrice_Morse_Sym& MBrr =  ref_cast(Matrice_Morse_Sym,matrice(0,0).valeur());
  Matrice_Morse& MBrv =  ref_cast(Matrice_Morse,matrice(0,1).valeur());
  MBrr.dimensionner(n2,0);
  MBrv.dimensionner(n2,0);

  IntVect& tab1RR=MBrr.tab1_;
  IntVect& tab2RR=MBrr.tab2_;
  IntVect& tab1RV=MBrv.tab1_;
  IntVect& tab2RV=MBrv.tab2_;

  // On traite les faces internes: dimensionnement des matrices reelles et virtuelles.
  int ndeb = zone_VDF.premiere_face_int();
  int nfin = zone_VDF.nb_faces();
  //  int pourcent=0;
  //  int tmp;

  IntVect rang_voisinRR(n2);
  IntVect rang_voisinRV(n2);
  rang_voisinRR=1;
  rang_voisinRV=0;

  for (num_face=ndeb; num_face<nfin; num_face++)
  {
  elem1 = face_voisins(num_face,0);
  elem2 = face_voisins(num_face,1);
  if (elem1 > elem2)
  {
  if(elem1<n2)
  {
  (rang_voisinRR(elem2))++;
  }
  else
  {
  if(elem2<n2)
  {
  (rang_voisinRV(elem2))++;
  }
  }
  }
  else // elem2 >= elem1
  {
  if(elem2<n2)
  {
  (rang_voisinRR(elem1))++;
  }
  else
  {
  if(elem1<n2)
  {
  (rang_voisinRV(elem1))++;
  }
  }
  }
  }

  // Prise en compte des conditions de type periodicite
  int i;
  const Conds_lim& les_cl = zone_Cl_dis().les_conditions_limites();

  for (i=0; i<les_cl.size(); i++)
  {
  const Cond_lim& la_cl = les_cl[i];

  if (sub_type(Periodique,la_cl.valeur()))
  {
  const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
  const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
  int numdeb = la_front_dis.num_premiere_face();
  int nfaces = la_front_dis.nb_faces();
  int ind_face_global;
  IntVect fait(nfaces);
  fait = 0;
  for (face=0; face<nfaces; face++)
  {
  if (fait[face] == 0)
  {
  fait[face] = 1;
  fait[la_cl_perio.face_associee(face)] = 1;
  ind_face_global = face+numdeb;
  elem1 = face_voisins(ind_face_global,0);
  elem2 = face_voisins(ind_face_global,1);
  assert(elem1>=0);
  assert(elem2>=0);
  if (elem1 > elem2)
  {
  if(elem1<n2)
  {
  (rang_voisinRR(elem2))++;
  }
  else
  {
  if(elem2<n2)
  {
  (rang_voisinRV(elem2))++;
  }
  }
  }
  else // elem2 >= elem1
  {
  if(elem2<n2)
  {
  (rang_voisinRR(elem1))++;
  }
  else
  {
  if(elem1<n2)
  {
  (rang_voisinRV(elem1))++;
  }
  }
  }
  }
  }
  }
  }

  tab1RR(0)=1;
  tab1RV(0)=1;
  for(i=0; i<n2; i++)
  {
  tab1RR(i+1)=rang_voisinRR(i)+tab1RR(i);
  tab1RV(i+1)=rang_voisinRV(i)+tab1RV(i);
  }
  MBrr.dimensionner(n2,tab1RR(n2)-1);
  MBrv.dimensionner(n2,n1-n2,tab1RV(n2)-1);
  for(i=0; i<n2; i++)
  {
  tab2RR[tab1RR[i]-1]=i+1;
  rang_voisinRR[i]=tab1RR[i];
  rang_voisinRV[i]=tab1RV[i]-1;
  }


  for (num_face=ndeb; num_face<nfin; num_face++)
  {
  elem1 = face_voisins(num_face,0);
  elem2 = face_voisins(num_face,1);

  if (elem1 > elem2)
  {
  if(elem1<n2)
  {
  tab2RR[rang_voisinRR[elem2]]=elem1+1;
  rang_voisinRR[elem2]++;
  }
  else
  {
  if(elem2<n2)
  {
  tab2RV[rang_voisinRV[elem2]]=(elem1-n2)+1;
  rang_voisinRV[elem2]++;
  }
  }
  }
  else
  {
  if(elem2<n2)
  {
  tab2RR[rang_voisinRR[elem1]]=elem2+1;
  rang_voisinRR[elem1]++;
  }
  else
  {
  if(elem1<n2)
  {
  tab2RV[rang_voisinRV[elem1]]=(elem2-n2)+1;
  rang_voisinRV[elem1]++;
  }
  }
  }
  }
  Cerr << finl;
  // On traite les conditions aux limites
  for (i=0; i<les_cl.size(); i++)
  {
  // Le traitement depend du type de la condition aux limites :
  const Cond_lim& la_cl = les_cl[i];
  const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
  int ndeb = la_front_dis.num_premiere_face();
  int nfin = ndeb + la_front_dis.nb_faces();

  if (sub_type(Periodique,la_cl.valeur()) )
  {
  const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
  int ind_face_local;
  IntVect fait(nfin-ndeb);
  fait = 0;
  for (num_face=ndeb; num_face<nfin; num_face++)
  {
  ind_face_local = num_face - ndeb;
  if (fait[ind_face_local] == 0)
  {
  fait[ind_face_local] = 1;
  fait[la_cl_perio.face_associee(ind_face_local)] = 1;
  elem1 = face_voisins(num_face,0);
  elem2 = face_voisins(num_face,1);

  // diagonale :
  if (elem1 > elem2) {
  if(elem1<n2) {
  tab2RR[rang_voisinRR[elem2]]=elem1+1;
  rang_voisinRR[elem2]++;
  }
  else {
  if(elem2<n2) {
  tab2RV[rang_voisinRV[elem2]]=elem1+1;
  rang_voisinRV[elem2]++;
  }
  }
  }
  else {
  if(elem2<n2) {
  tab2RR[rang_voisinRR[elem1]]=elem2+1;
  rang_voisinRR[elem1]++;
  }
  else {
  if(elem1<n2) {
  tab2RV[rang_voisinRV[elem1]]=elem2+1;
  rang_voisinRV[elem1]++;
  }
  }
  }
  }
  }
  }
  }
  }*/


/*
  void Eq_rayo_semi_transp_VDF::dimensionner_Mat_Bloc_Morse(Matrice& matrice_tmp)
  {
  const Zone_VDF& zone_VDF = ref_cast(Zone_VDF,zone_dis().valeur());
  const IntTab& face_voisins = zone_VDF.face_voisins();

  int elem1,elem2;
  int num_face,face;
  int n1 = nb_colonnes_tot();
  int n2 = nb_colonnes();

  matrice_tmp.typer("Matrice_Bloc");
  Matrice_Bloc& matrice=ref_cast(Matrice_Bloc, matrice_tmp.valeur());
  matrice.dimensionner(1,2);
  matrice(0,0).typer("Matrice_Morse");
  matrice(0,1).typer("Matrice_Morse");

  Matrice_Morse& MBrr =  ref_cast(Matrice_Morse,matrice(0,0).valeur());
  Matrice_Morse& MBrv =  ref_cast(Matrice_Morse,matrice(0,1).valeur());
  MBrr.dimensionner(n2,0);
  MBrv.dimensionner(n2,0);

  IntVect& tab1RR=MBrr.tab1_;
  IntVect& tab2RR=MBrr.tab2_;
  IntVect& tab1RV=MBrv.tab1_;
  IntVect& tab2RV=MBrv.tab2_;

  // On traite les faces internes: dimensionnement des matrices reelles et virtuelles.
  int ndeb = zone_VDF.premiere_face_int();
  int nfin = zone_VDF.nb_faces();
  //  int pourcent=0;
  //  int tmp;

  IntVect rang_voisinRR(n2);
  IntVect rang_voisinRV(n2);
  rang_voisinRR=1;
  rang_voisinRV=0;

  for (num_face=ndeb; num_face<nfin; num_face++)
  {
  elem1 = face_voisins(num_face,0);
  elem2 = face_voisins(num_face,1);
  if ((elem1 < n2) && (elem2 < n2))
  {
  (rang_voisinRR(elem2))++;
  (rang_voisinRR(elem1))++;
  }
  else
  {
  if(elem2<n2)
  (rang_voisinRV(elem2))++;
  else
  (rang_voisinRV(elem1))++;
  }
  }

  // Prise en compte des conditions de type periodicite
  int i;
  const Conds_lim& les_cl = zone_Cl_dis().les_conditions_limites();

  for (i=0; i<les_cl.size(); i++)
  {
  const Cond_lim& la_cl = les_cl[i];

  if (sub_type(Periodique,la_cl.valeur()))
  {
  const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
  const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
  int numdeb = la_front_dis.num_premiere_face();
  int nfaces = la_front_dis.nb_faces();
  int ind_face_global;
  IntVect fait(nfaces);
  fait = 0;
  for (face=0; face<nfaces; face++)
  {
  if (fait[face] == 0)
  {
  fait[face] = 1;
  fait[la_cl_perio.face_associee(face)] = 1;
  ind_face_global = face+numdeb;
  elem1 = face_voisins(ind_face_global,0);
  elem2 = face_voisins(ind_face_global,1);
  assert(elem1>=0);
  assert(elem2>=0);
  if ((elem1 < n2) && (elem2 < n2))
  {
  (rang_voisinRR(elem2))++;
  (rang_voisinRR(elem1))++;
  }
  else
  {
  if(elem2<n2)
  (rang_voisinRV(elem2))++;
  else
  (rang_voisinRV(elem1))++;
  }
  }
  }
  }
  }

  tab1RR(0)=1;
  tab1RV(0)=1;
  for(i=0; i<n2; i++)
  {
  tab1RR(i+1)=rang_voisinRR(i)+tab1RR(i);
  tab1RV(i+1)=rang_voisinRV(i)+tab1RV(i);
  }
  MBrr.dimensionner(n2,tab1RR(n2)-1);
  MBrv.dimensionner(n2,n1-n2,tab1RV(n2)-1);
  for(i=0; i<n2; i++)
  {
  tab2RR[tab1RR[i]-1]=i+1;
  rang_voisinRR[i]=tab1RR[i];
  rang_voisinRV[i]=tab1RV[i]-1;
  }

  for (num_face=ndeb; num_face<nfin; num_face++)
  {
  elem1 = face_voisins(num_face,0);
  elem2 = face_voisins(num_face,1);
  if ((elem1 < n2) && (elem2 < n2))
  {
  tab2RR[rang_voisinRR[elem2]]=elem1+1;
  tab2RR[rang_voisinRR[elem1]]=elem2+1;
  rang_voisinRR[elem2]++;
  rang_voisinRR[elem1]++;
  }
  else
  {
  if(elem2<n2)
  {
  tab2RV[rang_voisinRV[elem2]]=(elem1-n2)+1;
  rang_voisinRV[elem2]++;
  }
  else
  {
  tab2RV[rang_voisinRV[elem1]]=(elem2-n2)+1;
  rang_voisinRV[elem1]++;
  }
  }
  }

  // On traite les conditions aux limites
  for (i=0; i<les_cl.size(); i++)
  {
  // Le traitement depend du type de la condition aux limites :
  const Cond_lim& la_cl = les_cl[i];
  const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
  int ndeb = la_front_dis.num_premiere_face();
  int nfin = ndeb + la_front_dis.nb_faces();

  if (sub_type(Periodique,la_cl.valeur()) )
  {
  const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
  int ind_face_local;
  IntVect fait(nfin-ndeb);
  fait = 0;
  for (num_face=ndeb; num_face<nfin; num_face++)
  {
  ind_face_local = num_face - ndeb;
  if (fait[ind_face_local] == 0)
  {
  fait[ind_face_local] = 1;
  fait[la_cl_perio.face_associee(ind_face_local)] = 1;
  elem1 = face_voisins(num_face,0);
  elem2 = face_voisins(num_face,1);

  // diagonale :
  if ((elem1 < n2) && (elem2 < n2))
  {
  tab2RR[rang_voisinRR[elem2]]=elem1+1;
  tab2RR[rang_voisinRR[elem1]]=elem2+1;
  rang_voisinRR[elem2]++;
  rang_voisinRR[elem1]++;
  }
  else
  {
  if(elem2<n2)
  {
  tab2RV[rang_voisinRV[elem2]]=(elem1-n2)+1;
  rang_voisinRV[elem2]++;
  }
  else
  {
  tab2RV[rang_voisinRV[elem1]]=(elem2-n2)+1;
  rang_voisinRV[elem1]++;
  }
  }
  }
  }
  }
  }
  }*/


int Eq_rayo_semi_transp_VDF::nb_colonnes_tot()
{
  const Zone_VF& zone_VF = ref_cast(Zone_VF,zone_dis().valeur());
  return zone_VF.zone().nb_elem_tot();
}

int Eq_rayo_semi_transp_VDF::nb_colonnes()
{
  const Zone_VF& zone_VF = ref_cast(Zone_VF,zone_dis().valeur());
  return zone_VF.zone().nb_elem();
}
