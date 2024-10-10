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
// File:        Flux_radiatif_VEF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VEF
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#include <Flux_radiatif_VEF.h>
#include <Neumann_paroi_rayo_semi_transp_VEF.h>
#include <Temperature_imposee_paroi_rayo_semi_transp.h>
#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Eq_rayo_semi_transp_VEF.h>
#include <Champ_front_uniforme.h>
#include <Schema_Temps_base.h>
#include <Debog.h>
#include <Modele_rayo_semi_transp.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <Domaine_VEF.h>
#include <TRUST_Ref.h>

Implemente_instanciable(Flux_radiatif_VEF,"Flux_radiatif_VEF",Flux_radiatif_base);

/*! @brief Imprime le type de l'equation sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Flux_radiatif_VEF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}



/*! @brief
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Flux_radiatif_VEF::readOn(Entree& s )
{
  return Flux_radiatif_base::readOn(s);
}


/*! @brief
 *
 */
void Flux_radiatif_VEF::evaluer_cl_rayonnement(Champ_front_base& Tb, const Champ_Don_base&
                                               coeff_abs, const Champ_Don_base& longueur_rayo,
                                               const Champ_Don_base& indice,const Domaine_VF& zvf,
                                               const double sigma, double temps)
{
  const DoubleTab& n = indice.valeurs();
  const DoubleTab& epsilon = emissivite().valeurs();

  //  int nb_elem = zvf.nb_elem();
  const Front_VF& le_bord = ref_cast(Front_VF,frontiere_dis());
  //  const IntTab& face_voisins=zvf.face_voisins();

  // On dimensionne le DoubleTab associe a le_champ_front
  assert(le_champ_front->nb_comp() == 1);
  DoubleTab& tab = le_champ_front->valeurs_au_temps(temps);

  // Boucle sur les faces de le_bord
  int ndeb = le_bord.num_premiere_face();
  int nfin = ndeb + le_bord.nb_faces();
  int face = ndeb;
  for(face = ndeb; face<nfin; face++)
    {
      double epsi;
      assert(emissivite().nb_comp() == 1);
      if(sub_type(Champ_front_uniforme,emissivite()))
        epsi = epsilon(0,0);
      else
        epsi = epsilon(face-ndeb,0);

      double nn;
      assert(indice.nb_comp() == 1);
      if(sub_type(Champ_Uniforme,indice))
        nn = n(0,0);
      else
        nn = n(face,0);

      double T;
      assert(Tb.nb_comp() == 1);
      if(sub_type(Champ_front_uniforme,Tb))
        T = Tb.valeurs()(0,0);
      else
        T = Tb.valeurs_au_temps(temps)(face-ndeb,0);

      double numer_coeff = 4*nn*nn*sigma*pow(T,4)*epsi;
      double denum_coeff = A()*(2-epsi);

      tab(face-ndeb,0) = numer_coeff/denum_coeff;
    }
  tab.echange_espace_virtuel();
}


void Flux_radiatif_VEF::calculer_flux_radiatif(const Equation_base& eq_temp)
{
  //Cerr<<"Flux_radiatif_VEF::calculer_flux_radiatif : Debut"<<finl;
  // On doit recuperer la temperature de bord
  const Front_VF& le_bord = ref_cast(Front_VF,frontiere_dis());
  int nb_faces = le_bord.nb_faces();
  OBS_PTR(Champ_front_base) Tb;
  const Conds_lim& les_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites();
  int num_cl_temp = 0;

  int test_nom=0;
  for(num_cl_temp = 0; num_cl_temp<les_cl_temp.size(); num_cl_temp++)
    {
      const Cond_lim& la_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites(num_cl_temp);
      Nom nom_cl_temp = la_cl_temp->frontiere_dis().le_nom();
      if(nom_cl_temp == frontiere_dis().le_nom())
        {
          test_nom = 1;
          if (sub_type(Neumann_paroi_rayo_semi_transp_VEF,la_cl_temp.valeur()))
            {
              Neumann_paroi_rayo_semi_transp_VEF& la_cl_temper
                = ref_cast_non_const(Neumann_paroi_rayo_semi_transp_VEF,la_cl_temp.valeur());
              // WEC
              //               la_cl_temper.calculer_temperature_bord();
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Temperature_imposee_paroi_rayo_semi_transp,la_cl_temp.valeur()))
            {
              Temperature_imposee_paroi_rayo_semi_transp& la_cl_temper
                = ref_cast_non_const(Temperature_imposee_paroi_rayo_semi_transp,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl_temp.valeur()))
            {
              const Frontiere_ouverte_temperature_imposee_rayo_semi_transp& la_cl_temper
                = ref_cast(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Frontiere_ouverte_rayo_semi_transp,la_cl_temp.valeur()))
            {
              const Frontiere_ouverte_rayo_semi_transp& la_cl_temper
                = ref_cast(Frontiere_ouverte_rayo_semi_transp,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else
            {
              Cerr <<"Coder pour les autres condition limites de l'equation de temperature 1 "<<finl;
              exit();
            }
        }
    }
  if(test_nom == 0)
    {
      Cerr<<"Erreur : il n'y a pas de condition limite sur une frontiere portant le nom : "<<frontiere_dis().le_nom()<<finl;
      exit();
    }
  //  Tb contient les temperatures de bord
  // Calcul du flux radiatif
  DoubleTab& Flux = flux_radiatif().valeurs();
  Flux.resize(le_bord.nb_faces(),1);
  Eq_rayo_semi_transp_VEF& eq_rayo = ref_cast( Eq_rayo_semi_transp_VEF,domaine_Cl_dis().equation());
  Fluide_base& fluide = eq_rayo.fluide();
  DoubleTab& indice = fluide.indice().valeurs();
  DoubleTab& irradiance = eq_rayo.inconnue().valeurs();

  const Domaine_VEF& zvef = ref_cast(Domaine_VEF,domaine_Cl_dis().domaine_dis());
  //  const IntTab& face_voisins = zvef.face_voisins();
  const DoubleTab& face_normales = zvef.face_normales();

  double bilan_flux=0;
  // On fait une boucle sur les faces
  int ndeb = le_bord.num_premiere_face();
  int face=0;
  //Cerr << "On traite le bord de nom " << frontiere_dis().le_nom() << finl;
  for(face=0; face<nb_faces; face++)
    {
      double epsi;
      assert(emissivite().nb_comp() == 1);
      if (sub_type(Champ_front_uniforme,emissivite()))
        epsi = emissivite().valeurs()(0,0);
      else
        epsi = emissivite().valeurs()(face,0);
      double n;
      assert(fluide.indice().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide.indice()))
        n = indice(0,0);
      else
        n = indice(face+ndeb,0);
      double sigma = eq_rayo.Modele().valeur_sigma();
      double Tbord;
      assert(Tb->nb_comp() == 1);
      if(sub_type(Champ_front_uniforme,Tb.valeur()))
        Tbord = Tb->valeurs()(0,0);
      else
        Tbord = Tb->valeurs()(face,0);
      double irra = irradiance(face+ndeb);

      double denum = A()*(2-epsi);
      double numer = epsi*(irra - 4*n*n*sigma*pow(Tbord,4));
      Flux(face,0) = -numer/denum;
      //      Cerr << "face+ndeb,Tbord,irra" << " " << face+ndeb << " " << Tbord << " " << irra << finl;
      //      Cerr << "denum,numer " << denum << " " << numer << finl;
      //      Cerr << "Flux(face,0) " << " " << Flux(face,0) << finl;

      // Calcul des bilans
      double surface=0;
      int i;
      for (i=0; i<dimension; i++) surface += (face_normales(face+ndeb,i) * face_normales(face+ndeb,i));
      surface = sqrt(surface);
      bilan_flux += surface*Flux(face,0);
    }

  //Debog::verifier(" Flux_radiatif_VEF::calculer_flux_radiatif_irra",irradiance);
  //Debog::verifier_bord(" Flux_radiatif_VEF::calculer_flux_radiatif_Tb",Tb.valeurs()(0,0),ndeb);
  Debog::verifier_bord(" Flux_radiatif_VEF::calculer_flux_radiatif_Flux",Flux,ndeb);

  if (eq_rayo.schema_temps().limpr())
    Cout << "Flux radiatif sur le bord "<<le_bord.le_nom()<<" : "<<bilan_flux<<finl;
  //Cerr<<"Flux_radiatif_VEF::calculer_flux_radiatif : Fin"<<finl;
}
