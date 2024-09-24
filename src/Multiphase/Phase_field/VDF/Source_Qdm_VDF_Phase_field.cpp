/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Source_Qdm_VDF_Phase_field.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src/VDF
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Qdm_VDF_Phase_field.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Probleme_base.h>
#include <Milieu_base.h>
#include <Navier_Stokes_phase_field.h>
#include <Source_Con_Phase_field.h>
#include <SFichier.h>

Implemente_instanciable(Source_Qdm_VDF_Phase_field,"Source_Qdm_Phase_field_VDF_Face",Source_base);
// XD source_qdm_phase_field source_base source_qdm_phase_field 1 Keyword to define the capillary force into the Navier Stokes equation for the Phase Field problem.
// XD attr forme_du_terme_source entier forme_du_terme_source 0 Kind of the source term (1, 2, 3 or 4).


//// printOn
//

Sortie& Source_Qdm_VDF_Phase_field::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Source_Qdm_VDF_Phase_field::readOn(Entree& is )
{
  Cerr<<"Source_Qdm_VDF_Phase_field::readOn"<<finl;

  Motcle motlu;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="{")
    {
      Cerr<<"On attendait { dans Source_Qdm_VDF_Phase_field::readOn"<<finl;
      exit();
    }

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="forme_du_terme_source")
    {
      Cerr<<"On attendait forme_du_terme_source dans Source_Qdm_VDF_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> terme_source;

  // Le nouveau modele est incompressible
  compressible=0;

  is>>motlu;
  Cerr << motlu << finl;
  if (motlu!="}")
    {
      Cerr<<"On attendait }  dans Source_Qdm_VDF_Phase_field::readOn"<<finl;
      exit();
    }

  switch(terme_source)
    {
    case 1:
      {
        methode=&Source_Qdm_VDF_Phase_field::methode_1;
      };
      break;

    case 2:
      {
        methode=&Source_Qdm_VDF_Phase_field::methode_2;
      };
      break;

    case 3:
      {
        methode=&Source_Qdm_VDF_Phase_field::methode_3;
      };
      break;

    case 4:
      {
        methode=&Source_Qdm_VDF_Phase_field::methode_4;
      };
      break;

    default:
      Cerr<<"Le choix de la methode :"<<terme_source<<", n'est pas valide."<<finl;
      exit();
      break;
    }

  return is ;
}

void Source_Qdm_VDF_Phase_field::associer_pb(const Probleme_base& pb)
{
  le_probleme2=pb;
  Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));
  rho0=eq_ns.rho0();
  boussi_=eq_ns.get_boussi_();
  if (boussi_!=1 && boussi_!=0)
    {
      Cerr << "Erreur dans le choix du parametre boussi_" << finl;
      exit();
    }

  eq_ns.getset_terme_source_()=terme_source;
  eq_ns.getset_compressible_()=compressible;
}

void Source_Qdm_VDF_Phase_field::associer_domaines(const Domaine_dis_base& domaine_dis,
                                                   const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_VDF = ref_cast(Domaine_VDF, domaine_dis);
  le_dom_Cl_VDF = ref_cast(Domaine_Cl_VDF, domaine_Cl_dis);
}


DoubleTab& Source_Qdm_VDF_Phase_field::methode_1(DoubleTab& resu) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const DoubleVect& volumes = domaine_VDF.volumes();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const int nb_comp = eq_c.constituant().nb_constituants();

  double cface;
  int ndeb=domaine_VDF.premiere_face_int();
  int nbfaces=domaine_VDF.nb_faces();
  int el0,el1;
  double vol0,vol1;

  const Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));
  DoubleTab rhoPF=eq_ns.rho()->valeurs();
  double rho_face;

  // Forme en c*Grad(mutilde)
  //=============================

  // on calcule Grad(mutilde)
  //-----------------------------

  DoubleVect u_carre;
  DoubleTab mutilde_NS;
  Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  int type_systeme_naire = source_pf.get_type_systeme_naire();
  mutilde_NS.resize(eq_c.get_mutilde_().valeurs().dimension_tot(0), 1);

  const DoubleTab& mutilde=eq_c.get_mutilde_().valeurs();

  for (int i = 0; i < eq_c.get_mutilde_().valeurs().dimension_tot(0); i++) mutilde_NS(i,0) = mutilde(i);

  u_carre=source_pf.get_u_carre();

  if (eq_c.get_mutype_())
    {
      int taille=mutilde.size_totale();
      for (int i=0; i<taille; i++)
        {
          mutilde_NS(i)-=(0.5*u_carre(i))*source_pf.drhodc(i);
        }
    }
  // Dans le cas mutype_==1, on utilise mutilde_d dans CH, mais mutilde dans NS

  if (type_systeme_naire==0)
    {
      DoubleTab& grad_mutilde=ref_cast_non_const(DoubleTab, grad_mutilde_);
      if (grad_mutilde.size()==0) grad_mutilde=eq_ns.inconnue().valeurs();
      grad_mutilde=0.;
      const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
      opgrad.calculer(mutilde_NS,grad_mutilde);

      // on interpole c et on calcule la source
      //------------------------------------------------

      for (int fac=ndeb; fac<nbfaces; fac++)
        {
          el0=face_voisins(fac,0);
          el1=face_voisins(fac,1);
          vol0=volumes(el0);
          vol1=volumes(el1);

          cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);
          if (boussi_==1)
            {
              resu(fac) -= cface*grad_mutilde(fac) / rho0; // Cas approximation de Boussinesq
            }
          else if (boussi_==0)
            {
              rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);
              resu(fac) -= cface*grad_mutilde(fac) / rho_face;
            }
        }
    }
  else if (type_systeme_naire==1)
    {
      ///// Kim2012 terme source somme c.Grad(mutilde)

      DoubleTab& grad_mutilde=ref_cast_non_const(DoubleTab,  grad_mutilde_);
      if (grad_mutilde.size()==0) grad_mutilde=eq_ns.inconnue().valeurs();
      grad_mutilde=0.;
      DoubleTab temp_mutilde_NS(mutilde_NS.dimension_tot(0),1);
      const Operateur_Grad& opgrad=eq_ns.operateur_gradient();

      for (int j=0; j<nb_comp; j++)
        {
          grad_mutilde=0.;
          for (int i=0; i<temp_mutilde_NS.dimension(0); i++)
            {
              temp_mutilde_NS(i,0)=mutilde_NS(i,j);
            }
          opgrad.calculer(temp_mutilde_NS, grad_mutilde);

          // on interpole c et on calcule la source
          //------------------------------------------------

          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);
              cface=(vol0*c(el0,j)+vol1*c(el1,j))/(vol0+vol1);
              if (boussi_==1)
                {
                  resu(fac) -= cface*grad_mutilde(fac) / rho0; // Cas approximation de Boussinesq
                }
              else if (boussi_==0)
                {
                  rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);
                  resu(fac) -= cface*grad_mutilde(fac) / rho_face;
                }
            }
        }
      //resu.echange_espace_virtuel(); //ajoute mr264902
    }
  //===============================================
  return resu;

}

DoubleTab& Source_Qdm_VDF_Phase_field::methode_2(DoubleTab& resu) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const DoubleVect& volumes = domaine_VDF.volumes();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  const Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));

  double cface;
  int ndeb=domaine_VDF.premiere_face_int();
  int nbfaces=domaine_VDF.nb_faces();
  int el0,el1;
  double vol0,vol1;

  Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  int type_systeme_naire = source_pf.get_type_systeme_naire();

  // Forme en c*Grad(Div(alpha*rho*Grad(C)))
  //========================================

  // on calcule Grad(div_alpha_rho_gradC)
  //-------------------------------------
  if (type_systeme_naire==0)
    {
      const DoubleTab& div_alpha_rho_gradC=eq_c.get_div_alpha_rho_gradC();
      DoubleTab& grad_div_alpha_rho_gradC=ref_cast_non_const(DoubleTab,  grad_div_alpha_rho_gradC_);
      if (grad_div_alpha_rho_gradC.size()==0) grad_div_alpha_rho_gradC=eq_ns.inconnue().valeurs();
      grad_div_alpha_rho_gradC=0.;
      const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
      opgrad.calculer(div_alpha_rho_gradC, grad_div_alpha_rho_gradC);

      DoubleTab rhoPF=eq_ns.rho()->valeurs();
      double rho_face;

      // Division par rho0 necessaire dans le cas incompressible
      //--------------------------------------------------------
      if (boussi_==1)
        {
          if(compressible==0) grad_div_alpha_rho_gradC /= rho0; // Cas approximation de Boussinesq
        }
      else if (boussi_==0)
        {
          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);

              rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);
              grad_div_alpha_rho_gradC(fac) /= rho_face;
            }
        }

      // on interpole c et on calcule la source
      //------------------------------------------------

      for (int fac=ndeb; fac<nbfaces; fac++)
        {
          el0=face_voisins(fac,0);
          el1=face_voisins(fac,1);
          vol0=volumes(el0);
          vol1=volumes(el1);

          cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);
          if (boussi_==1)
            {
              resu(fac) += cface*grad_div_alpha_rho_gradC(fac) / rho0; // Cas approximation de Boussinesq
            }
          else if (boussi_==0)
            {
              rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);
              resu(fac) += cface*grad_div_alpha_rho_gradC(fac) / rho_face;
            }
        }
    }
  return resu;
}

DoubleTab& Source_Qdm_VDF_Phase_field::methode_3(DoubleTab& resu) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const DoubleVect& volumes = domaine_VDF.volumes();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  const Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));

  double cface;
  int ndeb=domaine_VDF.premiere_face_int();
  int nbfaces=domaine_VDF.nb_faces();
  int el0,el1;
  double vol0,vol1;

  Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  int type_systeme_naire = source_pf.get_type_systeme_naire();

  // Forme en c*Grad(Div(alpha*rho*Grad(C)))-alpha*rho*Grad((Grad(C))^2)/2
  //======================================================================

  // on calcule Grad(div_alpha_rho_gradC)
  //-------------------------------------
  if (type_systeme_naire==0)
    {
      const DoubleTab& div_alpha_rho_gradC=eq_c.get_div_alpha_rho_gradC();
      DoubleTab& grad_div_alpha_rho_gradC=ref_cast_non_const(DoubleTab,  grad_div_alpha_rho_gradC_);
      if (grad_div_alpha_rho_gradC.size()==0) grad_div_alpha_rho_gradC=eq_ns.inconnue().valeurs();
      grad_div_alpha_rho_gradC=0.;
      const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
      opgrad.calculer(div_alpha_rho_gradC, grad_div_alpha_rho_gradC);

      DoubleTab rhoPF=eq_ns.rho()->valeurs();
      double rho_face;

      // Division par rho0 necessaire dans le cas incompressible
      //--------------------------------------------------------
      if (boussi_==1)
        {
          if(compressible==0) grad_div_alpha_rho_gradC /= rho0; // Cas approximation de Boussinesq
        }
      else if (boussi_==0)
        {
          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);

              rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);
              grad_div_alpha_rho_gradC(fac) /= rho_face;
            }
        }

      // on calcule Grad(alpha_gradC_carre)
      //-----------------------------------

      const DoubleTab& alpha_gradC_carre=eq_c.get_alpha_gradC_carre();
      DoubleTab& grad_alpha_gradC_carre=ref_cast_non_const(DoubleTab,  grad_alpha_gradC_carre_);
      if (grad_alpha_gradC_carre.size()==0) grad_alpha_gradC_carre=eq_ns.inconnue().valeurs();
      grad_alpha_gradC_carre=0.;
      opgrad.calculer(alpha_gradC_carre,grad_alpha_gradC_carre);

      // on interpole c et on calcule la source
      //------------------------------------------------

      if((compressible==0 && boussi_==1) || boussi_==0)
        {
          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);

              cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);
              rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);

              if (boussi_==1)
                {
                  resu(fac) += (cface*grad_div_alpha_rho_gradC(fac) - grad_alpha_gradC_carre(fac)/2.) / rho0; // Incompressible - Cas approximation de Boussinesq
                }
              else if (boussi_==0)
                {
                  resu(fac) += (cface*grad_div_alpha_rho_gradC(fac) - grad_alpha_gradC_carre(fac)/2.) / rho_face;
                }
            }
        }
      else if (compressible==1 && boussi_==1)
        {
          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);

              cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);

              // Quasi-compressible
              //-------------------
              resu(fac) += (cface*grad_div_alpha_rho_gradC(fac) - rho0*grad_alpha_gradC_carre(fac)/2.); // Cas approximation de Boussinesq
            }
        }
    }
  return resu;
}

DoubleTab& Source_Qdm_VDF_Phase_field::methode_4(DoubleTab& resu) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const DoubleVect& volumes = domaine_VDF.volumes();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  const Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));

  //   double cface;
  int ndeb=domaine_VDF.premiere_face_int();
  int nbfaces=domaine_VDF.nb_faces();
  int el0,el1;
  double vol0,vol1;
  Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  int type_systeme_naire = source_pf.get_type_systeme_naire();

  // Forme en -Div(alpha*rho*Grad(C))*Grad(C)
  //=========================================

  // on calcule Grad(C)
  //-------------------
  if (type_systeme_naire==0)
    {
      DoubleTab& gradC=ref_cast_non_const(DoubleTab,  gradC_);
      if (gradC.size()==0) gradC=eq_ns.inconnue().valeurs();
      gradC=0.;
      const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
      opgrad.calculer(c,gradC);
      const DoubleTab& div_alpha_rho_gradC=eq_c.get_div_alpha_rho_gradC();
      double div_alpha_rho_gradCface;

      DoubleTab rhoPF=eq_ns.rho()->valeurs();
      double rho_face;

      if((compressible==0 && boussi_==1) || boussi_==0)
        {
          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);

              div_alpha_rho_gradCface=(vol0*div_alpha_rho_gradC(el0)+vol1*div_alpha_rho_gradC(el1))/(vol0+vol1);

              // Division par rho0 necessaire dans le cas incompressible
              //--------------------------------------------------------
              if (boussi_==1)
                {
                  div_alpha_rho_gradCface /= rho0;
                  resu(fac) -= div_alpha_rho_gradCface*gradC(fac) / rho0; // Cas approximation de Boussinesq
                }
              else if (boussi_==0)
                {
                  rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);
                  div_alpha_rho_gradCface /= rho_face;
                  resu(fac) -= div_alpha_rho_gradCface*gradC(fac) / rho_face;
                }

            }
        }
      else if (compressible==1 && boussi_==1)
        {
          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);

              div_alpha_rho_gradCface=(vol0*div_alpha_rho_gradC(el0)+vol1*div_alpha_rho_gradC(el1))/(vol0+vol1);

              resu(fac) -= div_alpha_rho_gradCface*gradC(fac) / rho0; // Cas approximation de Boussinesq

            }
        }
    }
  //=====================================================================

  return resu;
}

DoubleTab& Source_Qdm_VDF_Phase_field::ajouter(DoubleTab& resu) const
{
  //   const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  //   const IntTab& face_voisins = domaine_VDF.face_voisins();
  //   const DoubleVect& volumes = domaine_VDF.volumes();

  (this->*methode)(resu);

  //Cerr<<"max de la force "<<grad_mutilde.max_abs()<<finl;
  return resu;
}

DoubleTab& Source_Qdm_VDF_Phase_field::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

void Source_Qdm_VDF_Phase_field::mettre_a_jour(double temps)
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));

  const DoubleTab& inco=eq_ns.inconnue().valeurs();
  // GF a revoir car ne sert a rien en ce moment
  if (0)
    {
      double etot=0.;
      int nb_face=inco.dimension(0);

      for (int i=0; i<nb_face; i++)
        etot+=inco(i)*inco(i);
      Nom nom_fich=Objet_U::nom_du_cas();
      nom_fich+="_ecc";
      SFichier fic (nom_fich,ios::app);
      fic.precision(16);
      fic<<temps<<" "<<etot<<" Etot*volume si dx=dy=dz=cste!!! "<<finl;
      // Cout.precision(8);
    }
}

