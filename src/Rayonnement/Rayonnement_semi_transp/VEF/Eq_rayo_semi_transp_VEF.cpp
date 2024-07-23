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
// File:        Eq_rayo_semi_transp_VEF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VEF
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#include <Eq_rayo_semi_transp_VEF.h>
#include <Flux_radiatif_VEF.h>
#include <Neumann_paroi_rayo_semi_transp_VEF.h>
#include <Temperature_imposee_paroi_rayo_semi_transp.h>
#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Champ_front_uniforme.h>
#include <Debog.h>
#include <Operateur_Diff_base.h>
#include <Modele_rayo_semi_transp.h>
#include <Fluide_Incompressible.h>
#include <Champ_Uniforme.h>
#include <Domaine_VEF.h>
#include <Symetrie.h>
#include <TRUST_Ref.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Eq_rayo_semi_transp_VEF,"Eq_rayo_semi_transp_VEF",Equation_rayonnement_base);

/*! @brief Imprime le type de l'equation sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Eq_rayo_semi_transp_VEF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}



/*! @brief
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Eq_rayo_semi_transp_VEF::readOn(Entree& s )
{
  return Equation_rayonnement_base::readOn(s);
}


/*! @brief Calcule le champ de l'irradiance en fonction des CLs au temps donne en parametre et de la temperature au temps par defaut.
 *
 *     Met le resultat dans inconnue.valeurs()
 *     Suppose que le schema en temps est Euler explicite...
 *
 */
void Eq_rayo_semi_transp_VEF::resoudre(double temps)
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, domaine_dis().valeur());
  const DoubleVect& volumes_entrelaces =  domaine_VF.volumes_entrelaces();
  int nb_faces = domaine_VF.nb_faces();

  //
  // Remarque : l'assemblage de la matrice est realise une fois pour toute au debut du
  //            calcul, ce qui n'est valable que pour les cas ou kappa ne depend pas du
  //            temps

  //calcul du second membre
  DoubleTrav secmem(inconnue().valeurs());
  Probleme_base& pb = Modele().probleme();
  double n;
  double k;

  assert(pb.equation(1).inconnue().le_nom()=="temperature");
  const DoubleTab& temper = pb.equation(1).inconnue()->valeurs();

  const DoubleTab& indice = fluide().indice().valeurs();
  const DoubleTab& kappa = fluide().kappa().valeurs();
  double sigma = Modele().valeur_sigma();

  secmem = 0;
  int face;
  for (face=0; face<nb_faces; face++)
    {
      assert(fluide().indice().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide().indice().valeur()))
        n = indice(0,0);
      else
        n = indice(face,0);

      assert(fluide().kappa().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide().kappa().valeur()))
        k = kappa(0,0);
      else
        k = kappa(face,0);

      double vol = volumes_entrelaces(face);
      double T = temper(face);
      secmem(face) += +4*n*n*sigma*pow(T,4)*k*vol;
    }

  // On met a jour les champs associes aux conditions aux limites
  // avant d'evaluer leur contribution dans la matrice de discretisation
  evaluer_cl_rayonnement(temps);
  terme_diffusif->contribuer_au_second_membre(secmem);
  secmem.echange_espace_virtuel();

  if (solveur->que_suis_je() == "Solv_GCP")
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
  else if (solveur->que_suis_je() == "Solv_Gmres")
    if (Process::nproc() == 1)
      solveur.resoudre_systeme(la_matrice,secmem,irradiance_.valeurs());
    else
      {
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


  Debog::verifier("Eq_rayo_semi_transp_VEF::resoudre irradiance",irradiance_.valeurs());
  //Cerr<<"irradiance.std::max() = "<<irradiance_.std::max()<<", irradiance.std::min() = "<<irradiance_.std::min()<<finl;
  //Cerr<<"Eq_rayo_semi_transp_VEF::resoudre : Fin"<<finl;
}




/*! @brief modifie la matrice pour prendre en compte la presence de faces rayonnantes au voisinage des elements de bord
 *
 * @return le flot d'entree modifie
 */
void Eq_rayo_semi_transp_VEF::modifier_matrice()
{
  // On fait une boucle sur les conditions aux limites associees a l'equations
  Conds_lim& les_cl = domaine_Cl_dis().les_conditions_limites();
  const Domaine_VEF& zvef = ref_cast(Domaine_VEF,domaine_dis().valeur());
  const IntTab& face_voisins=zvef.face_voisins();
  const DoubleTab& face_normales = zvef.face_normales();

  int num_cl=0;
  for(num_cl = 0; num_cl<les_cl.size(); num_cl++)
    {
      Cond_lim& la_cl = domaine_Cl_dis().les_conditions_limites(num_cl);
      if (sub_type(Flux_radiatif_VEF,la_cl.valeur()))
        {
          Flux_radiatif_VEF& cl_radiatif = ref_cast(Flux_radiatif_VEF,la_cl.valeur());
          const DoubleTab& epsilon = cl_radiatif.emissivite().valeurs();
          double A = cl_radiatif.A();

          if (sub_type(Front_VF,la_cl.frontiere_dis()))
            {
              Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              // Boucle sur les faces de le_bord
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              int face = ndeb;
              for(face = ndeb; face<nfin; face++)
                {
                  int elem = face_voisins(face,0);
                  if(elem == -1)
                    elem = face_voisins(face,1);


                  double epsi;
                  assert(cl_radiatif.emissivite().nb_comp() == 1);
                  if (sub_type(Champ_front_uniforme,cl_radiatif.emissivite().valeur()))
                    epsi = epsilon(0,0);
                  else
                    epsi = epsilon(face-ndeb,0);

                  double surface=0;
                  int i;
                  for (i=0; i<dimension; i++) surface += (face_normales(face,i) * face_normales(face,i));
                  surface = sqrt(surface);

                  double coeff = epsi*surface;
                  coeff /= A*(2-epsi);

                  // On rajoute ce coefficient sur la diagonale de la matrice de discretisation
                  la_matrice(face,face) +=  coeff;
                }
            }
          else
            {
              Cerr<<"Erreur dans Eq_rayo_semi_transp_VEF::modifier_matrice()"<<finl;
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
}


void Eq_rayo_semi_transp_VEF::evaluer_cl_rayonnement(double temps)
{
  // Boucle sur les conditions aux limites de l'equation de rayonnement
  Conds_lim& les_cl_rayo = domaine_Cl_dis().les_conditions_limites();

  // recherche des conditions aux limites associes au l'equation de temperature
  Probleme_base& pb = Modele().probleme();
  Equation_base& eq_temp = pb.equation(1);

  Conds_lim& les_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites();

  int num_cl_rayo=0;
  for(num_cl_rayo = 0; num_cl_rayo<les_cl_rayo.size(); num_cl_rayo++)
    {
      Cond_lim& la_cl_rayo = domaine_Cl_dis().les_conditions_limites(num_cl_rayo);
      if(sub_type(Flux_radiatif_VEF,la_cl_rayo.valeur()))
        {
          Flux_radiatif_VEF& la_cl_rayon = ref_cast(Flux_radiatif_VEF,la_cl_rayo.valeur());
          // Recherche des temperatures de bord pour cette frontiere
          Nom nom_cl_rayo = la_cl_rayo.frontiere_dis().le_nom();
          int num_cl_temp = 0;
          REF(Champ_front) Tb;

          int test_remplissage_Tb = 0;
          for(num_cl_temp = 0; num_cl_temp<les_cl_temp.size(); num_cl_temp++)
            {
              Cond_lim& la_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites(num_cl_temp);
              Nom nom_cl_temp = la_cl_temp.frontiere_dis().le_nom();
              if(nom_cl_temp == nom_cl_rayo)
                {
                  if (sub_type(Neumann_paroi_rayo_semi_transp_VEF,la_cl_temp.valeur()))
                    {
                      Neumann_paroi_rayo_semi_transp_VEF& la_cl_temper
                        = ref_cast(Neumann_paroi_rayo_semi_transp_VEF,la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Temperature_imposee_paroi_rayo_semi_transp,la_cl_temp.valeur()))
                    {
                      Temperature_imposee_paroi_rayo_semi_transp& la_cl_temper
                        = ref_cast(Temperature_imposee_paroi_rayo_semi_transp,la_cl_temp.valeur());
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
                  else
                    {
                      Cerr<<"Coder pour les autres condition limites de l'equation de temperature 1 "<<finl;
                      Cerr<<"La condition a la limite "<<la_cl_rayo.que_suis_je()<<" n'est pas connue"<<finl;
                      exit();
                    }
                }
            }
          if ( test_remplissage_Tb == 0)
            // On n'a pas remplie le tableau des temperatures de bord !!!!
            Cerr<<"On n'a pas remplie le tableau des temperatures de bord !!!!"<<finl;

          const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis().valeur());
          la_cl_rayon.evaluer_cl_rayonnement(Tb.valeur(),fluide().kappa(), fluide().longueur_rayo(),
                                             fluide().indice(),zvf,Modele().valeur_sigma(),temps);
        }
      else if (sub_type(Symetrie,la_cl_rayo.valeur()))
        {
          ;
        }
      else
        {
          Cerr<<"La condition a la limite "<<la_cl_rayo.que_suis_je()<<" n'est pas connue"<<finl;
          Cerr<<" pour l'equation de rayonnement "<< finl;
          exit();
        }
    }
}


void Eq_rayo_semi_transp_VEF::assembler_matrice()
{

  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, domaine_dis().valeur());
  const DoubleVect& volumes_entrelaces =  domaine_VF.volumes_entrelaces();
  //int nb_faces = domaine_VF.nb_faces();

  const DoubleTab& irradi = irradiance_.valeurs();

  la_matrice.clean();

  // Prise en compte de la partie div((1/3K)grad(irradiance)) dans la matrice
  // de discretisation.

  terme_diffusif->contribuer_a_avec(irradi,la_matrice);


  // Modification de la matrice pour prendre en compte le second membre en K*irradiance
  const DoubleTab& kappa = fluide().kappa().valeurs();

  Cerr<<"On verifie lors du calcul de la matrice de discretisation que "<<finl;
  Cerr<<"l'ordre de la matrice est bien egale au nombre de faces"<<finl;
  if(la_matrice.ordre() != domaine_VF.nb_faces_tot())
    exit();
  Cerr<<"Ordre de la matrice OK"<<finl;

  int i;
  double k;
  for (i=0; i<la_matrice.ordre(); i++)
    {
      assert(fluide().kappa().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide().kappa().valeur()))
        k = kappa(0,0);
      else
        k = kappa(i,0);

      double vol = volumes_entrelaces(i);
      la_matrice(i,i) = la_matrice(i,i) + k*vol;
    }


  // On modifie la matrice pour prendre en compte l'effet des parois
  // rayonnantes sur les elements de bord
  modifier_matrice();
}

void Eq_rayo_semi_transp_VEF::completer()
{
  const Domaine_dis_base& un_domaine_dis = domaine_dis().valeur();
  int n = un_domaine_dis.nb_front_Cl();

  int ii;
  for (ii =0; ii<n; ii++)
    {
      const Frontiere_dis_base& la_fr_dis = un_domaine_dis.frontiere_dis(ii);
      le_dom_Cl_dis->les_conditions_limites(ii)->associer_fr_dis_base(la_fr_dis);
    }

  Equation_rayonnement_base::completer();

  // On assemble la matrice une fois pour toute au debut du calcul
  // Attention, ceci n'est valable que si kappa est constant au cours
  // du temps
  assembler_matrice();
}


void Eq_rayo_semi_transp_VEF::typer_op_grad()
{
  ;
}


int Eq_rayo_semi_transp_VEF::nb_colonnes_tot()
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis().valeur());
  return domaine_VF.nb_faces_tot();
}

int Eq_rayo_semi_transp_VEF::nb_colonnes()
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,domaine_dis().valeur());
  return domaine_VF.nb_faces();
}
