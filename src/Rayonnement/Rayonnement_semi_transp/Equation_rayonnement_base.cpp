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
// File:        Equation_rayonnement_base.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/30
//
//////////////////////////////////////////////////////////////////////////////

#include <Equation_rayonnement_base.h>
#include <Discret_Thyd.h>
#include <Modele_rayo_semi_transp.h>
#include <Operateur_Diff_base.h>
#include <Schema_Temps_base.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <Matrice_Bloc.h>
#include <Matrice_Morse_Sym.h>
#include <Param.h>

Implemente_base_sans_constructeur(Equation_rayonnement_base,"Equation_rayonnement_base",Equation_base);

Equation_rayonnement_base::Equation_rayonnement_base()
{
  /*
    Noms& nom=champs_compris_.liste_noms_compris();
    nom.dimensionner(1);
    nom[0]="irradiance";
  */
}

bool Equation_rayonnement_base::initTimeStep(double dt)
{
  schema_temps().set_dt()=dt;
  return Equation_base::initTimeStep(dt);
}

bool Equation_rayonnement_base::solve()
{
  resoudre(schema_temps().temps_courant()+schema_temps().pas_de_temps());
  return true;
}


/*! @brief Imprime le type de l'equation sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Equation_rayonnement_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}


/*! @brief cf Equation_base::readOn(Entree& is)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 * @throws solveur pression non defini dans jeu de donnees
 */
Entree& Equation_rayonnement_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  // typage de l'operateur de diffusion
  if (sub_type(Fluide_base,fluide()))
    if (fluide().is_rayo_semi_transp())
      {
        if (fluide().longueur_rayo_is_discretised())
          {
            const Champ_Don& long_rayo = fluide().longueur_rayo();
            terme_diffusif.associer_diffusivite(long_rayo);
          }
        else
          {
            Cerr<<"Error while reading the Radiation equation."<<finl;
            Cerr<<"You may not have discretized the problem of type Pb_Couple_rayo_semi_transp."<<finl;
            Cerr<<"Remark : The discretization of a coupling problem carries out those of the involved coupled "<<finl;
            Cerr<<"         problems therefore it is not useful to keep the discretization of these last ones.  "<<finl;
            Cerr<<finl;
            exit();
          }
      }
    else
      {
        Cerr<<"Error : the radiative properties of the incompressible fluid have not"<<finl;
        Cerr<<"been defined while a semi transparent radiation model is used."<<finl;
        Cerr<<"The fields kappa and indice which respectively define the absoption coefficient"<<finl;
        Cerr<<"and the refraction index must be added to the luid properties."<<finl;
        exit();
      }
  else
    {
      Cerr << "Error while reading the Radiation equation." << finl;
      Cerr << "Your fluid is of type " <<fluide().que_suis_je()<< finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered "<< finl;
      Cerr << "with the semi transparent radiation model."<<finl;
      exit();
    }

  Nom type="Op_Diff_";
  Nom discr = discretisation().que_suis_je();
  // les operateurs C_D_Turb_T sont communs aux discretisations VEF et VEFP1B
  if(discr=="VEFPreP1B") discr = "VEF";
  type +=discr;

  Nom nb_inc;
  if(sub_type(Champ_Uniforme,terme_diffusif.diffusivite()))
    nb_inc = "_";
  else
    nb_inc = "_Mult_inco_";
  type+=nb_inc;

  Nom type_inco=inconnue()->que_suis_je();
  type+=(type_inco.suffix("Champ_"));

  if (axi)
    type+="_Axi";

  terme_diffusif.typer(type);
  terme_diffusif.l_op_base().associer_eqn(*this);


  if (sub_type(Fluide_base,fluide()))
    if (fluide().kappa().non_nul())
      {
        const Champ_Don& long_rayo = fluide().longueur_rayo();
        terme_diffusif.valeur().associer_diffusivite(long_rayo);
        terme_diffusif.completer();
        terme_diffusif.valeur().dimensionner(la_matrice);
      }
    else
      {
        Cerr<<"Error : the radiative properties of the incompressible fluid have not"<<finl;
        Cerr<<"been defined while a semi transparent radiation model is used."<<finl;
        Cerr<<"The fields kappa and indice which respectively define the absoption coefficient"<<finl;
        Cerr<<"and the refraction index must be added to the luid properties."<<finl;
        exit();
      }
  else
    {
      Cerr << "Error while reading the Radiation equation." << finl;
      Cerr << "Your fluid is of type " <<fluide().que_suis_je()<< finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered "<< finl;
      Cerr << "with the semi transparent radiation model."<<finl;
      exit();
    }

  return is;
}

void Equation_rayonnement_base::set_param(Param& param)
{
  param.ajouter_non_std("conditions_limites|boundary_conditions",(this),Param::REQUIRED);
  param.ajouter_non_std("solveur",(this),Param::REQUIRED);
}

int Equation_rayonnement_base::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  int retval = 1;
  if (mot=="conditions_limites|boundary_conditions")
    {
      lire_cl(is);
      verif_Cl();
    }
  else if (mot=="solveur")
    {
      Cerr << "Reading and typing of the radiation equation solver :" << finl;
      Nom nom_solveur("Solv_");
      Nom type_solv_sys;
      is >> type_solv_sys;
      nom_solveur+=type_solv_sys;
      Cerr<<"Name of the radiation equation solver : "<<nom_solveur<<finl;
      solveur.typer(nom_solveur);
      is >> solveur.valeur();
      solveur.nommer("solveur_irradiance");
    }
  else retval = -1;

  return retval;
}

/*! @brief Associe un milieu physique a l'equation
 *
 * @param (Milieu_base& un_milieu) le milieu physique a associer a l'equation
 */
void Equation_rayonnement_base::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_base,un_milieu))
    if (fluide().kappa().non_nul())
      {
        const Fluide_base& un_fluide = ref_cast(Fluide_base,un_milieu);
        associer_fluide(un_fluide);
        le_fluide = un_fluide;
      }
    else
      {
        Cerr<<"Error : the radiative properties of the incompressible fluid have not"<<finl;
        Cerr<<"been defined while a semi transparent radiation model is used."<<finl;
        Cerr<<"The fields kappa and indice which respectively define the absoption coefficient"<<finl;
        Cerr<<"and the refraction index must be added to the luid properties."<<finl;
        exit();
      }
  else
    {
      Cerr << "Error for the method Equation_rayonnement_base::associer_milieu_base" << finl;
      Cerr << "Your fluid is of type " <<un_milieu.que_suis_je()<< finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered "<< finl;
      Cerr << "with the semi transparent radiation model."<<finl;
      exit();
    }
}

/*! @brief Associe le modele de rayonnement a l'equation de rayonnement
 *
 * @param (Modele_rayo_semi_transp& un_modele) le modele de rayonnement associe a l'equation de rayonnement
 */
void Equation_rayonnement_base::associer_modele_rayonnement(const Modele_rayo_semi_transp& un_modele)
{
  le_modele = un_modele;
}

/*! @brief Renvoie le milieu physique de l'equation (le Fluide_base upcaste en Milieu_base)
 *
 * @return (Milieu_base&) le Fluide_base de l'equation upcaste en Milieu_base
 */
const Milieu_base& Equation_rayonnement_base::milieu() const
{
  if (!le_fluide.non_nul())
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}


/*! @brief Renvoie le milieu physique de l'equation (le Fluide_base upcaste en Milieu_base)
 *
 *     (version const)
 *
 * @return (Milieu_base&) le Fluide_base de l'equation upcaste en Milieu_base
 */
Milieu_base& Equation_rayonnement_base::milieu()
{
  if (!le_fluide.non_nul())
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

/*! @brief Renvoie le nombre d'operateurs de l'equation.
 *
 * Ici 1.
 *
 * @return (int) le nombre d'operateurs de l'equation
 */
int Equation_rayonnement_base::nombre_d_operateurs() const
{
  return 1;
}


/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      exit si i>0
 *     (version const)
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 1 operateur
 */
const Operateur& Equation_rayonnement_base::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    default :
      Cerr << "Error for Equation_rayonnement_base::operateur(int i)" << finl;
      Cerr << "Equation_rayonnement_base has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  return terme_diffusif;
}


/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      exit si i>0
 *     (version const)
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 1 operateur
 */
Operateur& Equation_rayonnement_base::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    default :
      Cerr << "Error for Equation_rayonnement_base::operateur(int i)" << finl;
      Cerr << "Equation_rayonnement_base has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  return terme_diffusif;
}


/*! @brief Renvoie l'irradiance (champ inconnue de l'equation de rayonnement) (version const)
 *
 * @return (Champ_Inc&) le champ inconnue representant l'irradience
 */
const Champ_Inc& Equation_rayonnement_base::inconnue() const
{
  return irradiance_;
}


/*! @brief Renvoie l'irradiance (champ inconnue de l'equation de rayonnement) (version const)
 *
 * @return (Champ_Inc&) le champ inconnue representant l'irradience
 */
Champ_Inc& Equation_rayonnement_base::inconnue()
{
  return irradiance_;
}


/*! @brief Dicretise l'equation.
 *
 */
void Equation_rayonnement_base::discretiser()
{
  //
  // Discretisation de l'equation de rayonnement
  //
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  Cerr <<"Radiation equation discretisation" << finl;
  Cerr <<"Do not matter with the fact that discretization of the temperature"<<finl;
  Cerr <<"is indicated indeed it is the irradiance wich is discretized."<<finl;
  dis.temperature(schema_temps(), domaine_dis(), irradiance_);

  // La methode temperature(schema_temps(), domaine_dis(), irradiance_) permet
  // d'associer tout la bonne discretisation au champ de l'irradiance
  // toutefois, il faut modifier certaines grandeurs telles que le nom ou l'unite
  // de l'irradiance.
  irradiance_.valeur().fixer_nb_valeurs_temporelles(1);
  irradiance_.valeur().nommer("irradiance");
  irradiance_.valeur().fixer_unite("w/m2");
  champs_compris_.ajoute_champ(irradiance_);
  //
  // Fin de discretisation de l'equation de rayonnement
  //

  Equation_base::discretiser();
}


/*! @brief Renvoie la discretisation associee a l'equation.
 *
 * @return (Discretisation_base&) a discretisation associee a l'equation
 * @throws pas de probleme associe
 */
const Discretisation_base& Equation_rayonnement_base::discretisation() const
{
  //  if(!le_modele->probleme().non_nul())
  //  {
  //    Cerr << "Erreur : " << que_suis_je()
  //           << "n'a pas ete associee a un probleme ! " << finl;
  //    exit();
  //  }
  return le_modele->discretisation();
}

void Equation_rayonnement_base::completer()
{
  Equation_base::completer();
  //  gradient.completer();
}


void Equation_rayonnement_base::associer_pb_base(const Probleme_base& pb)
{
  Equation_base::associer_pb_base(pb);
  gradient.associer_eqn(*this);
}


Operateur_Grad& Equation_rayonnement_base::operateur_gradient()
{
  return gradient;
}

const Operateur_Grad& Equation_rayonnement_base::operateur_gradient() const
{
  return gradient;
}


void Equation_rayonnement_base::Mat_Morse_to_Mat_Bloc(Matrice& matrice_tmp)
{
  int n1 = nb_colonnes_tot();
  int n2 = nb_colonnes();

  Matrice_Bloc& matrice= ref_cast(Matrice_Bloc,matrice_tmp.valeur());
  Matrice_Morse& MBrr =  ref_cast(Matrice_Morse,matrice.get_bloc(0,0).valeur());
  Matrice_Morse& MBrv =  ref_cast(Matrice_Morse,matrice.get_bloc(0,1).valeur());

  IntVect& tab1RR=MBrr.get_set_tab1();
  IntVect& tab2RR=MBrr.get_set_tab2();
  DoubleVect& coeffRR=MBrr.get_set_coeff();
  IntVect& tab1RV=MBrv.get_set_tab1();
  IntVect& tab2RV=MBrv.get_set_tab2();
  DoubleVect& coeffRV=MBrv.get_set_coeff();

  DoubleTab ligne_tmp(n1);
  for(int i=0; i<n2; i++)
    {
      int k;
      // On recopie le premier bloc de la matrice dans un tableau :
      //      ligne_tmp = 0;
      for ( k=la_matrice.get_tab1()(i)-1; k<la_matrice.get_tab1()(i+1)-1; k++)
        ligne_tmp(la_matrice.get_tab2()(k) - 1) = la_matrice.get_coeff()(k);

      // On complete la partie reelle de la matrice
      for ( k=tab1RR(i)-1; k<tab1RR(i+1)-1; k++)
        coeffRR[k] = ligne_tmp(tab2RR[k] - 1);

      // On complete la partie virtuelle
      for ( k=tab1RV(i)-1; k<tab1RV(i+1)-1; k++)
        coeffRV[k] = ligne_tmp(n2 + tab2RV[k] - 1);
    }


  /*  Cerr<<"Impression de la_matrice"<<finl;
      la_matrice.imprimer_formatte(Cerr);
      Cerr<<"Impression de MBrr"<<finl;
      MBrr.imprimer_formatte(Cerr);
      Cerr<<"Impression de MBrv"<<finl;
      MBrv.imprimer_formatte(Cerr);

      Debog::verifier_Mat_faces("avant resolution systeme : la_matrice",la_matrice);
      Debog::verifier_Mat_faces("avant resolution systeme : MBrr",MBrr);*/
}


void Equation_rayonnement_base::dimensionner_Mat_Bloc_Morse_Sym(Matrice& matrice_tmp)
{
  int n1 = nb_colonnes_tot();
  int n2 = nb_colonnes();
  int iligne;
  const IntVect& tab1=la_matrice.get_set_tab1();
  const IntVect& tab2=la_matrice.get_set_tab2();

  matrice_tmp.typer("Matrice_Bloc");
  Matrice_Bloc& matrice=ref_cast(Matrice_Bloc, matrice_tmp.valeur());
  matrice.dimensionner(1,2);
  matrice.get_bloc(0,0).typer("Matrice_Morse_Sym");
  matrice.get_bloc(0,1).typer("Matrice_Morse");

  Matrice_Morse_Sym& MBrr =  ref_cast(Matrice_Morse_Sym,matrice.get_bloc(0,0).valeur());
  Matrice_Morse& MBrv =  ref_cast(Matrice_Morse,matrice.get_bloc(0,1).valeur());
  MBrr.dimensionner(n2,0);
  MBrv.dimensionner(n2,0);

  IntVect& tab1RR=MBrr.get_set_tab1();
  IntVect& tab2RR=MBrr.get_set_tab2();
  IntVect& tab1RV=MBrv.get_set_tab1();
  IntVect& tab2RV=MBrv.get_set_tab2();

  IntVect compteur_MBrr(n2);
  IntVect compteur_MBrv(n2);
  compteur_MBrr=0;
  compteur_MBrv=0;

  // On parcours les lignes de la_matrice pour compter les elements
  // non nuls de chaque ligne
  int jcolonne;
  for (iligne=0; iligne<n2; iligne++)
    {
      int k;
      for ( k=tab1(iligne)-1; k<tab1(iligne+1)-1; k++)
        {
          jcolonne = tab2(k)-1;
          if (jcolonne < n2)
            {
              // l'element correspondant est dans la partie RR de la_matrice
              if ((jcolonne >= iligne) && (jcolonne < n2))
                {
                  // l'element correspondant est  situe au dessus de la diagonale de la_matrice
                  compteur_MBrr(iligne)++;
                }
            }
          else
            {
              // l'element correspondant est dans la partie RV de la_matrice
              compteur_MBrv(iligne)++;
            }
        }
    }


  // On remplie tab1RR et tab1RV
  tab1RR(0)=1;
  tab1RV(0)=1;
  for(int i=0; i<n2; i++)
    {
      tab1RR(i+1)=compteur_MBrr(i)+tab1RR(i);
      tab1RV(i+1)=compteur_MBrv(i)+tab1RV(i);
    }
  // On dimensionne tab2RR et tab2RV
  MBrr.dimensionner(n2,tab1RR(n2)-1);
  MBrv.dimensionner(n2,n1-n2,tab1RV(n2)-1);

  // On remplit tab2RR et tab2RV
  int compteurRR,compteurRV;
  for (iligne=0; iligne<n2; iligne++)
    {
      int k;
      compteurRR = tab1RR(iligne)-1;
      compteurRV = tab1RV(iligne)-1;
      for ( k=tab1(iligne)-1; k<tab1(iligne+1)-1; k++)
        {
          jcolonne = tab2(k)-1;
          if (jcolonne < n2)
            {
              // l'element correspondant est dans la partie RR de la_matrice
              if ((jcolonne >= iligne) && (jcolonne < n2))
                {
                  // l'element correspondant est  situe au dessus de la diagonale de la_matrice
                  tab2RR(compteurRR) = tab2(k);
                  compteurRR++;
                }
            }
          else
            {
              // l'element correspondant est dans la partie RV de la_matrice
              tab2RV(compteurRV) = tab2(k)-n2;
              compteurRV++;
            }
        }
    }
}
