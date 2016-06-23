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
// File:        Source_Con_Phase_field.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Phase_field/src/VDF
// Version:     /main/25
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Con_Phase_field.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Probleme_base.h>
#include <Milieu_base.h>
#include <Navier_Stokes_phase_field.h>
#include <Dimension.h>
#include <Source_Qdm_VDF_Phase_field.h>
#include <Check_espace_virtuel.h>
#include <Convection_Diffusion_Phase_field.h>

Implemente_instanciable(Source_Con_Phase_field,"Source_Con_Phase_field_VDF_P0_VDF",Source_Con_Phase_field_base);


//// printOn
//

Sortie& Source_Con_Phase_field::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Source_Con_Phase_field::readOn(Entree& is )
{
  Cerr<<"Source_Con_Phase_field::readOn"<<finl;

  Motcle motlu;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="{")
    {
      Cerr<<"On attendait { dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="Temps_d_affichage")
    {
      Cerr<<"On attendait Temps_d_affichage dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> tpsaff;
  if(tpsaff<0. || tpsaff >100.)
    {
      Cerr << "Le temps d'affichage doit etre compris entre 0 et 100 secondes." << finl;
      exit();
    }

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="alpha")
    {
      Cerr<<"On attendait alpha dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> alpha;

  is>>motlu;
  Cerr << motlu << finl;
  if (motlu!="beta")
    {
      Cerr<<"On attendait beta dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> beta;

  is >>motlu;
  Cerr << motlu << finl;
  if (motlu!="kappa")
    {
      Cerr<<"On attendait kappa dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> kappa;

  is >>motlu;
  Cerr << motlu << finl;
  if (motlu!="kappa_variable")
    {
      Cerr<<"On attendait kappa_variable dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    {
      Motcle temp_type_kappa_;
      is >> temp_type_kappa_;
      if(temp_type_kappa_=="oui")
        {
          type_kappa_=1;
        }
      else
        {
          type_kappa_=0;
        }
    }

  is >>motlu;
  Cerr << motlu << finl;
  if (motlu!="moyenne_de_kappa")
    {
      Cerr<<"On attendait moyenne_de_kappa dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    {
      Motcle temp_kappa_moy_;
      is >> temp_kappa_moy_;
      if(temp_kappa_moy_=="arithmetique")
        {
          kappa_moy_=0;
        }
      else if(temp_kappa_moy_=="harmonique")
        {
          kappa_moy_=1;
        }
      else
        {
          kappa_moy_=2;
        }
    }

  is >>motlu;
  Cerr << motlu << finl;
  if (motlu!="multiplicateur_de_kappa")
    {
      Cerr<<"On attendait multiplicateur_de_kappa dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> mult_kappa;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="couplage_NS_CH")
    {
      Cerr<<"On attendait couplage_NS_CH dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    {
      Motcle temp_couplage_;
      is >> temp_couplage_;
      if(temp_couplage_=="mutilde(n)")
        {
          couplage_=0;
        }
      else
        {
          couplage_=1;
        }
    }

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="implicitation_CH")
    {
      Cerr<<"On attendait implicitation_CH dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    {
      Motcle temp_implicitation_;
      is >> temp_implicitation_;
      if(temp_implicitation_=="oui")
        {
          implicitation_=1;
        }
      else
        {
          implicitation_=0;
        }
    }

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="gmres_non_lineaire")
    {
      Cerr<<"On attendait gmres_non_lineaire dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    {
      Motcle temp_gmres_;
      is >> temp_gmres_;
      if(temp_gmres_=="oui")
        {
          gmres_=1;
        }
      else
        {
          gmres_=0;
        }
    }

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="seuil_cv_iterations_ptfixe")
    {
      Cerr<<"On attendait seuil_cv_iterations_ptfixe dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> epsilon_;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="seuil_residu_ptfixe")
    {
      Cerr<<"On attendait seuil_residu_ptfixe dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> eps_;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="seuil_residu_gmresnl")
    {
      Cerr<<"On attendait seuil_residu_gmresnl dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> epsGMRES;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="dimension_espace_de_krylov")
    {
      Cerr<<"On attendait dimension_espace_de_krylov dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> nkr;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="nb_iterations_gmresnl")
    {
      Cerr<<"On attendait nb_iterations_gmresnl dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> nit;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="residu_min_gmresnl")
    {
      Cerr<<"On attendait residu_min_gmresnl dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> rec_min;

  is >> motlu;
  Cerr << motlu << finl;
  if (motlu!="residu_max_gmresnl")
    {
      Cerr<<"On attendait residu_max_gmresnl dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }
  else
    is >> rec_max;

  is>>motlu;
  Cerr << motlu << finl;
  if (motlu!="}")
    {
      Cerr<<"On attendait }  dans Source_Con_Phase_field::readOn"<<finl;
      exit();
    }

  if(kappa>0)
    {
      Cerr <<" pas de temps theorique = ?     dx4*"<< 1./alpha/kappa/2.<<" dx2*"<<1./2./kappa/beta<<finl;
    }
  return is ;

}

void Source_Con_Phase_field::associer_pb(const Probleme_base& pb)
{
  le_probleme2=pb;

  Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));
  Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& rho=eq_ns.milieu().masse_volumique().valeurs();
  int dim=rho.nb_dim();
  switch(dim)
    {
    case 1:
      rho0=rho(0);
      break;
    case 2:
      rho0=rho(0,0);
      break;
    case 3:
      rho0=rho(0,0,0);
      break;
    default:
      Cerr <<"Sou_Con_Phase_field : Pb avec la dimension de rho :"<<dim<<finl;
      exit();
      break;
    }
  rho1 = eq_c.rho1();
  rho2 = eq_c.rho2();
  mu1 = eq_c.mu1();
  mu2 = eq_c.mu2();
  drhodc();
  // Dans le modele cette derivee est constante - On la calcule au debut.

  boussi_=eq_ns.get_boussi_();
  if (boussi_!=1 && boussi_!=0)
    {
      Cerr << "Erreur dans le choix du parametre boussi_" << finl;
      exit();
    }

  diff_boussi_=eq_ns.get_diff_boussi_();
  if (diff_boussi_!=1 && diff_boussi_!=0)
    {
      Cerr << "Erreur dans le choix du parametre diff_boussi_" << finl;
      exit();
    }

  g_=eq_ns.get_g_();

  mutype_=eq_c.get_mutype_();
  if (mutype_!=0 && mutype_!=1)
    {
      Cerr << "Erreur dans le choix du parametre mutype_" << finl;
      exit();
    }

  if(implicitation_==1)
    {
      if(gmres_!=0 && gmres_!=1)
        {
          Cerr << "** Erreur du choix de la methode implicite GMRES - L'execution doit stopper **" << finl;
          Cerr << "Choix possibles : 0, 1" << finl;
          exit();
        }
    }
  else if(implicitation_!=0)
    {
      Cerr<<"================================================="<<finl;
      Cerr<<"Choix de l'implicitation incorrect ! (ni 0, ni 1)"<<finl;
      Cerr<<"================================================="<<finl;
      exit();
    }

  if(type_kappa_==1)
    {
      kappa_ind=1;
    }
  else if(type_kappa_==0)
    {
      kappa_ind=0;
    }
  else
    {
      Cerr << "Erreur dans le choix de kappa !" << finl;
      exit();
    }

  if(implicitation_==1)
    {
      if(kappa_moy_!=0 && kappa_moy_!=1 && kappa_moy_!=2)
        {
          Cerr << "Erreur dans le choix de la moyenne de kappa !" << finl;
          exit();
        }
    }

  if(mult_kappa<=0)
    {
      Cerr << "Erreur dans le choix du multiplicateur de kappa !" << finl;
    }

  // Recapitulatif des parametres

  Nom choix_boussi;
  if(boussi_==1)
    {
      choix_boussi="oui";
    }
  else
    {
      choix_boussi="non";
    }

  Nom choix_diff_boussi;
  if(diff_boussi_==1)
    {
      choix_diff_boussi="oui";
    }
  else
    {
      choix_diff_boussi="non";
    }
  if(boussi_==1 && diff_boussi_==1)
    {
      Cerr << "ATTENTION : les options 'approximation de Boussinesq' et" << "'approximation de Boussinesq dans le terme de diffusion' ne devrait pas etre utilisees simultanement'." << "(voir Source_Con_Phase_field.cpp)" << finl;
    }

  const int terme_source=eq_ns.getset_terme_source_();

  Nom choix_source;
  if(terme_source==1)
    {
      choix_source="c grad(mutilde)";
    }
  else if(terme_source==2)
    {
      choix_source="c grad(laplacien(c))";
    }
  else if(terme_source==3)
    {
      choix_source="c grad(laplacien(c))-div((grad(c))^2)/2";
    }
  else if(terme_source==4)
    {
      choix_source="-laplacien(c) grad(c)";
    }

  Nom choix_implicite;
  if(implicitation_==1)
    {
      choix_implicite="oui";
    }
  else
    {
      choix_implicite="non";
    }

  Nom type_implicite;
  if(gmres_==1)
    {
      type_implicite="GMRES non lineaire";
    }
  else
    {
      type_implicite="point fixe";
    }

  Nom mutilde_d;
  if(mutype_==1)
    {
      mutilde_d="avec le terme d'Ec";
    }
  else
    {
      mutilde_d="sans le terme d'Ec";
    }

  Nom type_couplage;
  if(couplage_==1)
    {
      type_couplage="potentiel au temps n+1/2";
    }
  else
    {
      type_couplage="potentiel au temps n";
    }

  Nom mobilite_variable;
  if(type_kappa_==1)
    {
      mobilite_variable="oui";
    }
  else
    {
      mobilite_variable="non";
    }

  Nom moyenne_kappa;
  if(kappa_moy_==0)
    {
      moyenne_kappa="arithmetique (+)";
    }
  else if(kappa_moy_==1)
    {
      moyenne_kappa="harmonique (/)";
    }
  else if(kappa_moy_==2)
    {
      moyenne_kappa="geometrique (*)";
    }

  Cerr << "" << finl;
  Cerr << "" << finl;
  Cerr << "******************************************************************************" << finl;
  Cerr << "*********************** RECAPITULATIF DU PARAMETRAGE *************************" << finl;
  Cerr << "" << finl;
  Cerr << "1) Equation de Navier-Stokes" << finl;
  Cerr << "  - approximation de Boussinesq                     : " << choix_boussi << finl;
  if(boussi_==0)
    {
      Cerr << "  - masse volumique de la phase 1                   : " << rho1 << " kg/m3" << finl;
      Cerr << "  - masse volumique de la phase 2                   : " << rho2 << " kg/m3" << finl;
      Cerr << "  - approximation de Boussinesq dans la diffusion   : " << choix_diff_boussi << finl;
    }
  else
    {
      Cerr << "  - masse volumique de reference                    : " << rho0 << " kg/m3" << finl;
    }
  Cerr << "  Dans le cas ou la viscosite dynamique est variable : " << finl;
  Cerr << "  - viscosite dynamique de la phase 1               : " << mu1 << " kg/m/s" << finl;
  Cerr << "  - viscosite dynamique de la phase 2               : " << mu2 << " kg/m/s" << finl;
  Cerr << "  - terme source                                    : " << choix_source << finl;
  Cerr << "  - intensite du champ de gravite                   : " << norme_array(g_)  << " m/s^2" << finl;
  Cerr << "" << finl;
  Cerr << "2) Equation de Cahn-Hilliard" << finl;
  Cerr << "  - discretisation implicite                        : " << choix_implicite << finl;
  if(implicitation_==1)
    {
      Cerr << "  - algorithme implicite                            : " << type_implicite << finl;
    }
  Cerr << "  - potentiel chimique generalise                   : " << mutilde_d << finl;
  Cerr << "  - couplage NS / CH via le potentiel chimique      : " << type_couplage << finl;
  Cerr << "  - mobilite variable                               : " << mobilite_variable << finl;
  if(type_kappa_==1)
    {
      Cerr << "  - type de moyenne pour la mobilite                : " << moyenne_kappa << finl;
      Cerr << "  - multiplicateur de la mobilite                   : " << mult_kappa << finl;
    }
  if(implicitation_==1 && gmres_==0)
    {
      Cerr << "  - dans le cas du point fixe ... " << finl;
      Cerr << "        seuil de convergence entre deux iterations  : " << epsilon_ << finl;
      Cerr << "        seuil du residu a convergence               : " << eps_ << finl;
    }
  else if(implicitation_==1 && gmres_==1)
    {
      Cerr << "  - dans le cas du GMRES non lineaire ... " << finl;
      Cerr << "        seuil du residu a convergence               : " << epsGMRES << finl;
      Cerr << "        dimension de l'espace de Krylov             : " << nkr << finl;
      Cerr << "        nombre d'iterations maximal de l'algorithme : " << nit << finl;
      Cerr << "        residu minimum de convergence               : " << rec_min << finl;
      Cerr << "        residu maximum de convergence               : " << rec_max << finl;
    }
  Cerr << "" << finl;
  Cerr << "********************** Ce message persistera " << tpsaff << " secondes *********************" << finl;
  Cerr << "*************** Pour continuer avant les " << tpsaff << " secondes : Ctrl + C **************" << finl;
  Cerr << "******************************************************************************" << finl;
  Cerr << "" << finl;
  Cerr << "" << finl;

  Nom tps_sleep="sleep ";
  tps_sleep+=Nom(tpsaff);
  Cerr << system(tps_sleep) << finl;
}

void Source_Con_Phase_field::associer_zones(const Zone_dis& zone_dis,
                                            const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF, zone_Cl_dis.valeur());
}


inline double mobilite(const double& c)
{
  return (1.);
  //   return (c*c);
  //   return (0.5+2*c*c);
  //  return (min(1.,4.*c*c));
  //   return(0.);
  //   const double clim1 = -0.4;
  //   const double clim2 =  0.4;
  //   if (c<=clim1)
  //     return(min(1.,-(c-clim1)/(0.5+clim1)));
  //   else if (c>=clim2)
  //     return(min(1.,(c-clim2)/(0.5-clim2)));
  //   else return(0.);
}

DoubleTab& Source_Con_Phase_field::laplacien(const DoubleTab& F, DoubleTab& resu) const
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv=eq_ns.operateur_divergence();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const DoubleVect& volumes = zone_VDF.volumes();

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;

  // Grad(F)
  resu=0.;
  opgrad.calculer(F,prov_face);
  // M*Grad(F)
  int ndeb=zone_VDF.premiere_face_int();
  int nbfaces=zone_VDF.nb_faces();
  int el0,el1;
  double cface,vol0,vol1;
  for (int fac=ndeb; fac<nbfaces; fac++)
    {
      el0=face_voisins(fac,0);
      el1=face_voisins(fac,1);
      vol0=volumes(el0);
      vol1=volumes(el1);
      cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);
      prov_face(fac)=mobilite(cface)*prov_face(fac);
    }
  //   int taille=prov_face.size();
  //   Cerr << "taille : " << taille << finl;
  //   for (int i=0;i<taille;i++)
  //     {
  //       prov_face(i)*=mobilite(c(i));
  //     }
  //   Cerr << "Fin multiplication mobilite" << finl;
  // Application solveur masse
  eq_ns.solv_masse().appliquer(prov_face);

  // Div(M*Grad(F))
  opdiv.calculer(prov_face,resu);

  return resu;
}


DoubleTab& Source_Con_Phase_field::div_kappa_grad(const DoubleTab& F, const DoubleTab& kappa_var, DoubleTab& resu) const
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv=eq_ns.operateur_divergence();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const DoubleVect& volumes = zone_VDF.volumes();

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;

  resu=0.;

  // Grad(F)
  opgrad.calculer(F,prov_face);
  // M*Kappa*Grad(F)
  int ndeb=zone_VDF.premiere_face_int();
  int nbfaces=zone_VDF.nb_faces();
  int el0,el1;
  double cface,kappa_face,vol0,vol1;
  for (int fac=ndeb; fac<nbfaces; fac++)
    {
      el0=face_voisins(fac,0);
      el1=face_voisins(fac,1);
      vol0=volumes(el0);
      vol1=volumes(el1);
      kappa_face=(vol0*kappa_var(el0)+vol1*kappa_var(el1))/(vol0+vol1);
      cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);
      prov_face(fac)=kappa_face*mobilite(cface)*prov_face(fac);
    }
  // Application solveur masse (/M)
  eq_ns.solv_masse().appliquer(prov_face);

  // Div(Kappa*Grad(F))
  opdiv.calculer(prov_face,resu);

  return resu;
}


DoubleTab& Source_Con_Phase_field::ajouter(DoubleTab& resu) const
{
  //const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));

  resu+=accr;
  return resu;
}

DoubleTab& Source_Con_Phase_field::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Source_Con_Phase_field::mettre_a_jour(double temps)
{
  Cerr << "Temps : " << temps << " s" << finl;
  Cerr << "" << finl;

  Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  Champ_Fonc& ch=eq_c.set_mutilde_();
  ch.mettre_a_jour(temps);
}


// Description:
//     Calcule le premier demi pas de temps dans le cas implicite
//     Calcule le pas de temps dans le cas explicite
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::premier_demi_dt()
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const int nb_elem = zone_VDF.nb_elem_tot();

  Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  DoubleTab& c=eq_c.inconnue().valeurs();

  DoubleTab& mutilde=eq_c.set_mutilde_().valeurs();
  DoubleTab& mutilde_demi=eq_c.set_mutilde_demi();
  DoubleTab& c_demi=eq_c.set_c_demi();

  // Utilise dans l'ancien modele de Didier Jamet
  DoubleTab& alpha_gradC_carre=eq_c.set_alpha_gradC_carre();
  calculer_alpha_gradC_carre(alpha_gradC_carre);
  DoubleTab& div_alpha_rho_gradC=eq_c.set_div_alpha_rho_gradC();
  calculer_div_alpha_rho_gradC(div_alpha_rho_gradC);
  // ---

  DoubleTab& div_alpha_gradC=eq_c.set_div_alpha_gradC();
  calculer_div_alpha_gradC(div_alpha_gradC);
  DoubleTab& pression_thermo=eq_c.set_pression_thermo();
  calculer_pression_thermo(pression_thermo);


  accr = c;
  accr = 0.;
  c_demi = c;
  c_demi = 0.;

  if(implicitation_==1)
    {
      // Cerr<<"Schema implicite"<<finl;
      Matrice_Morse matrice_diffusion_CH;
      mutilde_demi = c;
      mutilde_demi = 0.;

      // Pour le GMRES NL ----------------------------

      // Commente par DJ
      //----------------
      //       DoubleTab x1(2*nb_elem);
      //       x1=0.;
      //       if(gmres_==1)
      //         {

      //           // Ecriture de x1(0)
      //           for(int n_elem=0;n_elem<nb_elem;n_elem++)
      //             {
      //               x1(n_elem)=c(n_elem);
      //               x1(n_elem+nb_elem)=mutilde(n_elem);
      //             }
      //         }
      //----------------
      // ---------------------------------------------

      calculer_u2_elem(u_carre_);
      assembler_matrice_point_fixe(matrice_diffusion_CH);

      if(gmres_==0)
        {
          calculer_point_fixe(c, mutilde, matrice_diffusion_CH, c_demi, mutilde_demi);
        }
      else if(gmres_==1)
        {
          // Modifie par DJ
          //---------------
          //           non_lin_gmres(c, mutilde, matrice_diffusion_CH, x1);
          non_lin_gmres(c, mutilde, matrice_diffusion_CH, c_demi, mutilde_demi);
          //---------------

          const double theta=0.6;
          // On stocke les nouveaux c et mutilde
          for(int n_elem=0; n_elem<nb_elem; n_elem++)
            {
              // Commente par DJ
              //----------------
              //               c_demi(n_elem)=x1(n_elem);
              //----------------
              c_demi(n_elem)-=(1-theta)*c(n_elem);
              c_demi(n_elem)/=theta;

              // Commente par DJ
              //----------------
              //               mutilde_demi(n_elem)=x1(n_elem+nb_elem);
              //----------------
              mutilde_demi(n_elem)-=(1-theta)*mutilde(n_elem);
              mutilde_demi(n_elem)/=theta;
            }
        }

      // ATTENTION : A VERIFIER
      //=======================
      //       c_demi.echange_espace_virtuel();
      //       mutilde_demi.echange_espace_virtuel();
      //       c.echange_espace_virtuel();
      //       mutilde.echange_espace_virtuel();

      // Mise a jour

      if (couplage_==0)
        {
          // Si traitement explicite de mutilde dans le cas implicite
          calculer_mutilde(mutilde);
        }
      else if (couplage_==1)
        {
          // Si traitement implicite de mutilde dans le cas implicite
          mutilde=mutilde_demi;
        }
      mutilde.echange_espace_virtuel();

      // Calcul de l'accroissement entre n et n+1/2
      // Utile pour pouvoir utiliser n'importe quel dt
      accr=0.;
    }
  else if(implicitation_==0)
    {
      // Cerr<<"Schema explicite"<<finl;
      calculer_u2_elem(u_carre_);
      calculer_mutilde(mutilde);
      mutilde.echange_espace_virtuel();

      DoubleTab& prov_elem=ref_cast(DoubleTab, prov_elem_);
      if (prov_elem.size()==0)
        prov_elem=mutilde;

      DoubleTab kappa_var(prov_elem);
      kappa_var=0.;

      if(kappa_ind==1)
        {
          for(int ikappa=0; ikappa<nb_elem; ikappa++)
            kappa_var(ikappa)=dmax(-mult_kappa*kappa*(c(ikappa)-0.5)*(c(ikappa)+0.5),0);
          // Div(Kappa*Grad(mutilde))
          div_kappa_grad(mutilde, kappa_var, prov_elem);
        }
      else
        {
          // Kappa*Laplacien(mutilde)
          laplacien(mutilde,prov_elem);
          prov_elem*=kappa;
        }

      // Pour equation Allen-Cahn (kappa constant)
      //------------------------------------------
      //   prov_elem=F;
      //   prov_elem*=-kappa;

      accr+=prov_elem;
      //eq_c.solv_masse().appliquer(prov_elem);
      //Cerr<<"masse de la force sur c "<<prov_elem.max_abs()<<finl;

      // Utile pour pouvoir utiliser n'importe quel dt - voir Schema_Phase_Field
      c_demi=c;
    }
  else
    {
      Cerr<<"Type de schema errone !!"<<finl;
    }
}

// Description:
//     Calcul de Div(alpha*rho*Grad((C)) au centre des elements
// Precondition:
// Parametre: DoubleTab& div_alpha_gradC
//    Signification: Div(alpha*Grad((C)) au centre des elements
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::calculer_div_alpha_gradC(DoubleTab& div_alpha_gradC) const
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv= eq_ns.operateur_divergence();

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;

  // Grad(C)
  opgrad.calculer(c,prov_face);
  eq_ns.solv_masse().appliquer(prov_face);

  prov_face *= alpha;

  // Div(alpha*rho*Grad(c))
  opdiv.calculer(prov_face,div_alpha_gradC);
  eq_c.solv_masse().appliquer(div_alpha_gradC);
}

// Description:
//     Calcul de Div(alpha*rho*Grad((C)) au centre des elements
// Precondition:
// Parametre: DoubleTab& div_alpha_rho_gradC
//    Signification: Div(alpha*rho*Grad((C)) au centre des elements
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::calculer_div_alpha_rho_gradC(DoubleTab& div_alpha_rho_gradC) const
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv= eq_ns.operateur_divergence();

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const DoubleVect& volumes = zone_VDF.volumes();
  const int ndeb=zone_VDF.premiere_face_int();
  const int nbfaces=zone_VDF.nb_faces();

  int el0,el1;
  double vol0,vol1;

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;

  // Grad(C)
  opgrad.calculer(c,prov_face);
  eq_ns.solv_masse().appliquer(prov_face);

  const DoubleTab rhoPF=eq_ns.rho().valeurs();
  double rho_face;

  if (boussi_==1)
    {
      prov_face *= alpha*rho0; // Cas approximation de Boussinesq
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
          prov_face(fac) *= alpha*rho_face;
        }
    }

  // Div(alpha*rho*Grad(c))
  opdiv.calculer(prov_face,div_alpha_rho_gradC);
  eq_c.solv_masse().appliquer(div_alpha_rho_gradC);
}


// Description:
//    Assemble la matrice pour le calcul du point fixe
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::assembler_matrice_point_fixe(Matrice_Morse& matrice_diffusion_CH)
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const DoubleTab& positions = zone_VDF.xp();
  const IntVect& ori = zone_VDF.orientation();

  int compt=0;
  int compt1=0;
  int compt2=0;

  const int nb_elem = zone_VDF.nb_elem_tot();
  int nb_compo_;
  int f0;
  int min_tri;
  int old_tri;
  int voisin;
  int nb_faces_au_bord;
  int dimensionnement ;
  double dvar2=0.;
  double dvarkeep=0.;
  double dt;
  double valeur_diag=0.;
  const double theta=0.6;
  double kappa_interpolee=0.;

  //   Cerr<<"======================================"<<finl;
  //   Cerr<<"Assemblage de la matrice du point fixe"<<finl;

  // Dimension du pb
  nb_compo_=dimension;

  // Forme de la matrice : ( I A )
  //                       ( B I )

  // Assemblage de la moitie superieure de la matrice : I et A

  // On dimensionne la matrice A et donc en meme temps B (meme nombre de termes non nuls)
  // Pour cela : dimension du vecteur coeff = (2 * dim + 1) * nb_elem - nombre_de_bords

  DoubleTab kappa_var(c);

  // Cerr << "kappa_ind = " << kappa_ind << finl;
  if(kappa_ind==0)
    {
      // Cerr << "---" << finl;
      // Cerr << "kappa constant pour le point fixe" << finl;
      // Cerr << "---" << finl;
      kappa_var=kappa;
    }
  else
    {
      for(int ikappa=0; ikappa<nb_elem; ikappa++)
        {
          // Cerr << "---" << finl;
          // Cerr << "Calcul de kappa variable pour le point fixe" << finl;
          // Cerr << "---" << finl;
          kappa_var(ikappa)=dmax(-mult_kappa*kappa*(c(ikappa)-0.5)*(c(ikappa)+0.5),0);
        }
    }

  dimensionnement=(2*nb_compo_+1)*nb_elem;

  for(int elem=0; elem<nb_elem; elem++)
    {
      for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
        {
          f0 = elem_faces(elem,ncomp);
          voisin=face_voisins(f0,0);
          if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

          if (voisin==-1) dimensionnement--;
        }
    }
  // On ajoute le nombre de 1 venant de I, ceci valant pour I et A. Le dimensionnement total est donc le double
  dimensionnement+=nb_elem;
  dimensionnement*=2;

  // Allocation des tableaux specifiques a la matrice morse
  matrice_diffusion_CH.dimensionner(2*nb_elem,dimensionnement);
  DoubleVect& coeff=matrice_diffusion_CH.coeff_;
  IntVect& tab2=matrice_diffusion_CH.tab2_;
  IntVect& tab1=matrice_diffusion_CH.tab1_;

  // Boucle sur le nombre d'elements
  for(int elem=0; elem<nb_elem; elem++)
    {
      valeur_diag=0.;
      old_tri=-1;
      nb_faces_au_bord=0;

      // On ajoute le 1 de la diagonale de la sous-matrice I du haut a gauche
      coeff(compt)=1;
      tab2(compt2)=elem;
      compt++;
      compt2++;
      tab1(compt1)=compt2;
      compt1++;

      // Calcul du nombre de bords
      for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
        {
          f0 = elem_faces(elem,ncomp);
          voisin=face_voisins(f0,0);
          if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

          // On connait le voisin associe a la face en cours. On regarde s'il est au bord
          if (voisin==-1) nb_faces_au_bord++;
        }

      // Boucle sur les faces
      for(int ncomp=0; ncomp<2*nb_compo_-nb_faces_au_bord; ncomp++)
        {

          min_tri=nb_elem+1;
          dt=eq_c.schema_temps().pas_de_temps();    // on fait un calcul sur un demi pas de temps

          for(int ncomp_tri=0; ncomp_tri<2*nb_compo_; ncomp_tri++)
            {
              f0 = elem_faces(elem,ncomp_tri);

              // Voisin associe a la face - On s'assure de traiter un voisin et pas l'element resident
              voisin=face_voisins(f0,0);
              if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

              // On va compter le nombre de voisins existants (faces non aux bords de l'element elem)
              // Ceci sert a calculer la contribution au terme relatif a l'element elem
              // Cas kappa variable : on utilise la moyenne geometrique des valeurs des mobilites
              // des elements concernes par la face de calcul ("kappa_interpolee")
              if (voisin!=-1 && old_tri==-1)
                {
                  dvar2=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                  if(kappa_moy_==2)
                    {
                      kappa_interpolee=pow(dabs(kappa_var(elem)*kappa_var(voisin)),0.5);    // Moyenne geometrique
                    }
                  else if(kappa_moy_==1)
                    {
                      if (kappa_var(elem)==0 || kappa_var(voisin)==0)
                        kappa_interpolee=0;
                      else
                        kappa_interpolee=2./(1./kappa_var(elem)+1./kappa_var(voisin));  // Moyenne harmonique
                    }
                  else if(kappa_moy_==0)
                    {
                      kappa_interpolee=(kappa_var(elem)+kappa_var(voisin))/2.;        // Moyenne arithmetique
                    }
                  valeur_diag+=theta*kappa_interpolee*(dt/dvar2);
                }

              // Si le voisin a un numero plus petit que l'ancien minimum, on ne le prend pas en compte
              // En effet, le nouveau minimum est un numero d'element qui lui sera superieur
              // Rq : voisin>old_tri => voisin != -1 donc pour une face au bord le if est false
              if (voisin>old_tri)
                {
                  min_tri=min(min_tri,voisin);
                  if(min_tri==voisin)
                    {
                      dvarkeep=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                    }
                }
            }

          // Si l'ancien minimum est inferieur a elem et le nouveau superieur, il faut traiter elem juste avant le nouveau voisin
          if (old_tri<elem && min_tri>elem)
            {
              coeff(compt)=valeur_diag;
              tab2(compt2)=elem+nb_elem;
              compt++;
              compt2++;
            }

          // On complete les tableaux avec les valeurs dans les voisins non aux bords du domaine
          if (kappa_moy_==2)
            {
              kappa_interpolee=pow(dabs(kappa_var(elem)*kappa_var(min_tri)),0.5);    // Moyenne geometrique
            }
          else if(kappa_moy_==1)
            {
              if (kappa_var(elem)==0 || kappa_var(min_tri)==0)
                kappa_interpolee=0;
              else
                kappa_interpolee=2./(1./kappa_var(elem)+1./kappa_var(min_tri));  // Moyenne harmonique
            }
          else if(kappa_moy_==0)
            {
              kappa_interpolee=(kappa_var(elem)+kappa_var(min_tri))/2.;        // Moyenne arithmetique
            }
          coeff(compt)=-theta*kappa_interpolee*(dt/dvarkeep);
          tab2(compt2)=min_tri+nb_elem;
          compt++;
          compt2++;

          // Si on est a la fin et que elem n'a pas encore ete mis dans coeff, il faut le faire avant la fin de la boucle sur les faces de l'element
          // Ceci arrive si elem > tous les numeros d'elements de ses voisins
          if (ncomp==2*nb_compo_-nb_faces_au_bord-1 && min_tri<elem)
            {
              coeff(compt)=valeur_diag;
              tab2(compt2)=elem+nb_elem;
              compt++;
              compt2++;
            }

          // On sauve le numero d'element "minimum"
          old_tri=min_tri;
        }

    }

  // Assemblage de la moitie superieure de la matrice : B et I

  // Boucle sur le nombre d'elements
  for(int elem=0; elem<nb_elem; elem++)
    {
      valeur_diag=0;
      old_tri=-1;
      nb_faces_au_bord=0;

      // Calcul du nombre de bords
      for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
        {
          f0 = elem_faces(elem,ncomp);
          voisin=face_voisins(f0,0);
          if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

          // On connait le voisin associe a la face en cours. On regarde s'il est au bord
          if (voisin==-1) nb_faces_au_bord++;
        }

      // Boucle sur les faces
      for(int ncomp=0; ncomp<2*nb_compo_-nb_faces_au_bord; ncomp++)
        {

          min_tri=nb_elem+1;

          for(int ncomp_tri=0; ncomp_tri<2*nb_compo_; ncomp_tri++)
            {
              f0 = elem_faces(elem,ncomp_tri);

              // Voisin associe a la face - On s'assure de traiter un voisin et pas l'element resident
              voisin=face_voisins(f0,0);
              if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

              // On va compter le nombre de voisins existants (faces non aux bords de l'element elem)
              // Ceci sert a calculer la contribution au terme relatif a l'element elem
              if (voisin!=-1 && old_tri==-1)
                {
                  dvar2=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                  valeur_diag+=-alpha/dvar2;
                }

              // Si le voisin a un numero plus petit que l'ancien minimum, on ne le prend pas en compte
              // En effet, le nouveau minimum est un numero d'element qui lui sera superieur
              // Rq : voisin>old_tri => voisin != -1 donc pour une face au bord le if est false
              if (voisin>old_tri)
                {
                  min_tri=min(min_tri,voisin);
                  if(min_tri==voisin)
                    {
                      dvarkeep=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                    }
                }

            }

          // Si l'ancien minimum est inferieur a elem et le nouveau superieur, il faut traiter elem juste avant le nouveau voisin
          if (old_tri<elem && min_tri>elem)
            {
              coeff(compt)=valeur_diag;
              tab2(compt2)=elem;
              compt++;
              compt2++;
              if (old_tri==-1)
                {
                  // On ajoute un terme a tab1 dans le cas ou l'element diagonal est le premier terme
                  tab1(compt1)=compt2;
                  compt1++;
                  old_tri=elem;
                }
            }

          // On complete les tableaux avec les valeurs dans les voisins non aux bords du domaine
          coeff(compt)=alpha/dvarkeep;
          tab2(compt2)=min_tri;
          compt++;
          compt2++;
          if (old_tri==-1)
            {
              // On ajoute un terme a tab1 dans le cas ou le premier terme est extradiagonal
              tab1(compt1)=compt2;
              compt1++;
            }

          // Si on est a la fin et que elem n'a pas encore ete mis dans coeff, il faut le faire avant la fin de la boucle sur les faces de l'element
          // Ceci arrive si elem > tous les numeros d'elements de ses voisins
          // Bien sur dans ce cas la, l'element diagonal est forcement le dernier mais pas le premier - il n'existe pas de maillage avec des mailles seules sans voisins
          if (ncomp==2*nb_compo_-nb_faces_au_bord-1 && min_tri<elem)
            {
              coeff(compt)=valeur_diag;
              tab2(compt2)=elem;
              compt++;
              compt2++;
            }

          // On sauve le numero d'element "minimum"
          old_tri=min_tri;

        }

      // On ajoute le 1 de la diagonale de la sous-matrice I du bas a droite
      coeff(compt)=1;
      tab2(compt2)=elem+nb_elem;
      compt++;
      compt2++;


    }

  // cf Mat_Morse.h/.cpp, tab1() et tab2() sont a utiliser au sens Fortran
  // EXPLICATION :
  // coeff stocke les valeurs des coefficients
  // tab2 stocke les numeros de colonnes dans la matrice de ces coefficients, au sens FORTRAN (i-e 1 <= tab2[i] <= n)
  // tab1 stocke le rang de l'element de tab2() pour lequel on change de ligne. Ainsi pour tab1[j], on est a la i-eme ligne avec :
  //      ( tab1[i] <= j < tab1[i+1] ) qui implique qu'on considere toujours la i-eme ligne,
  //      et des que cela n'est plus respecte, le passage a la ligne suivante est effectue

  // De par l'algorithme, tab1() est deja au sens FORTRAN. Reste a le faire pour tab2()
  // De plus on doit mettre le dernier terme de tab1() a dimensionnement+1
  tab2+=1;
  tab1(compt1)=dimensionnement+1;
  compt1+=1;

  if (compt!=dimensionnement)
    {
      Cerr << "Erreur lors du calcul de la matrice du point fixe : nombre d'elements non nuls calcules different du nombre d'elements non nuls prevus" << finl;
      exit();
    }
  //   Cerr<<"Nombre d'elements non nuls calcules="<<compt<<finl;
  //   Cerr<<"Nombre d'elements non nuls prevus="<<dimensionnement<<finl;
  //   Cerr<<"Assemblage de la matrice du point fixe : OK"<<finl;
  //   Cerr<<"==========================================="<<finl;

  // Test des tableaux
  // Cerr<<"coeff="<<coeff<<finl;
  // Cerr<<"tab1="<<tab1<<finl;
  // Cerr<<"tab2="<<tab2<<finl;
}


// Description:
//    Calcul du point fixe
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::calculer_point_fixe(const DoubleTab& c, const DoubleTab& mutilde, const Matrice_Morse& matrice_diffusion_CH, DoubleTab& c_demi, DoubleTab& mutilde_demi)
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const int nb_elem = zone_VDF.nb_elem_tot();
  const double theta=0.6;
  DoubleVect residu1(nb_elem);
  DoubleVect residu2(nb_elem);

  DoubleVect term_cin(nb_elem);
  // Ceci correspond a u2/2*drhodc=mutilde_dyn-mutilde

  DoubleVect second_membre(2*nb_elem);
  DoubleVect Phi(nb_elem);
  DoubleVect M(nb_elem);
  DoubleVect old_Phi(nb_elem);
  DoubleVect old_M(nb_elem);
  DoubleVect Solu_temp(2*nb_elem);

  // Initialisation des variables
  Phi=c;
  M=mutilde;
  old_Phi=0.;
  old_M=0.;

  // Test de convergence
  for(int n_elem=0; n_elem<nb_elem; n_elem++)
    {
      residu1(n_elem)=Phi(n_elem)-old_Phi(n_elem);
      residu2(n_elem)=M(n_elem)-old_M(n_elem);
    }

  Cerr<<"==========================="<<finl;
  Cerr<<"Assemblage du second membre"<<finl;
  Cerr<<"Inversion du systeme du point fixe"<<finl;
  Cerr<<"(Attention : fonctionne seulement en sequentiel)"<<finl;
  if (Process::nproc()>1)
    exit();

  while(norme_array(residu1) > epsilon_ || norme_array(residu2)>epsilon_)
    {
      // Assemblage du second membre
      for(int n_elem=0; n_elem<nb_elem; n_elem++)
        {
          second_membre(n_elem)=c(n_elem);
          term_cin(n_elem)=0.;
          if (mutype_==1)
            {
              term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc_;
            }
          second_membre(n_elem+nb_elem)=term_cin(n_elem)+beta*dpsidc(Phi(n_elem));
        }

      // Inversion de la matrice morse
      matrice_diffusion_CH.inverse(second_membre, Solu_temp, eps_);

      old_Phi=Phi;
      old_M=M;
      for(int n_elem=0; n_elem<nb_elem; n_elem++)
        {
          Phi(n_elem)=Solu_temp(n_elem);
          M(n_elem)=Solu_temp(n_elem+nb_elem);
        }
      Phi.echange_espace_virtuel();
      M.echange_espace_virtuel();

      for(int n_elem=0; n_elem<nb_elem; n_elem++)
        {
          residu1(n_elem)=Phi(n_elem)-old_Phi(n_elem);
          residu2(n_elem)=M(n_elem)-old_M(n_elem);
        }

    }

  // Ecriture de c et mutilde au temps n+1/2, cf le theta-schema:
  // Phi=theta*c_demi+(1.-theta)*c
  // M=theta*mutilde_demi+(1.-theta)*mutilde

  if(theta!=0)
    {
      c_demi=Phi;
      c_demi.ajoute(-(1.-theta),c);
      c_demi/=theta;
      mutilde_demi=M;
      mutilde_demi.ajoute(-(1.-theta),mutilde);
      mutilde_demi/=theta;
    }
  else
    {
      Cerr<<"Attention : le traitement est explicite (theta=0) - De plus pour theta < 0.5 ou theta > 1 le systeme est instable !!"<<finl;
      c_demi=c;
      mutilde_demi=mutilde;
    }

}


// Description:
//     Construire le residu du GMRES NL
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::construire_systeme(const DoubleTab& c, const Matrice_Morse& matrice_diffusion_CH, DoubleTab& v0, const DoubleTab& x1)
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const int nb_elem = zone_VDF.nb_elem_tot();

  DoubleVect term_cin(nb_elem);
  // Ceci correspond a u2/2*drhodc=mutilde_dyn-mutilde

  DoubleVect second_membre(2*nb_elem);
  DoubleTab Ax1(2*nb_elem);

  // Assemblage du second membre

  for(int n_elem=0; n_elem<nb_elem; n_elem++)
    {
      second_membre(n_elem)=c(n_elem);

      term_cin(n_elem)=0.;
      if (mutype_==1)
        {
          term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc_;
        }

      second_membre(n_elem+nb_elem)=term_cin(n_elem)+beta*dpsidc(x1(n_elem));
    }
  //   // Ajoute par DJ pour Debog
  //   //-------------------------
  //   {
  //     DoubleTab secmem_c(c);
  //     DoubleTab secmem_mutilde(c);
  //     for(int n_elem=0; n_elem<nb_elem; n_elem++)
  //       {
  //         secmem_c(n_elem) = second_membre(n_elem);
  //         secmem_mutilde(n_elem) = second_membre(n_elem+nb_elem);
  //       }
  //     Debog::verifier("Construire Systeme secmem_c : ",secmem_c);
  //     Debog::verifier("Construire Systeme secmem_mutilde : ",secmem_mutilde);
  //   }
  //   //-------------------------

  // Calcul du produit matrice / vecteur utilise
  matrice_diffusion_CH.multvect_(x1,Ax1);
  // Modifie par DJ
  //---------------
  {
    DoubleTab Ax1_c(c);
    DoubleTab Ax1_mutilde(c);
    for(int n_elem=0; n_elem<nb_elem; n_elem++)
      {
        Ax1_c(n_elem) = Ax1(n_elem);
        Ax1_mutilde(n_elem) = Ax1(n_elem+nb_elem);
      }

    Ax1_c.echange_espace_virtuel();
    Ax1_mutilde.echange_espace_virtuel();

    for(int n_elem=0; n_elem<nb_elem; n_elem++)
      {
        Ax1(n_elem) = Ax1_c(n_elem);
        Ax1(n_elem+nb_elem) = Ax1_mutilde(n_elem);
      }
  }
  //---------------
  //   // Ajouter par DJ pour Debog
  //   //--------------------------
  //   {
  //     DoubleTab Ax1_c(c);
  //     DoubleTab Ax1_mutilde(c);
  //     for(int n_elem=0; n_elem<nb_elem; n_elem++)
  //       {
  //         Ax1_c(n_elem) = Ax1(n_elem);
  //         Ax1_mutilde(n_elem) = Ax1(n_elem+nb_elem);
  //       }
  //     Debog::verifier("Construire Systeme Ax1_c : ",Ax1_c);
  //     Debog::verifier("Construire Systemeyes Ax1_mutilde : ",Ax1_mutilde);
  //   }
  //   //-------------------------

  // Calcul de v0 = [A(xn) xn - bn]
  for(int n_elem=0; n_elem<2*nb_elem; n_elem++)
    {
      v0(n_elem)=(Ax1(n_elem)-second_membre(n_elem));
    }

  return ;

}


// Description:
//     Construire le residu du GMRES NL
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::matvect(const DoubleTab& c, const Matrice_Morse& matrice_diffusion_CH, const DoubleTab& v0, const DoubleTab& x1, DoubleTab& v1)
{
  const double delta = 1.e-5;

  DoubleTab v2(v1);
  DoubleTab x1_(v0);

  v1 = 0.;

  x1_ = v0;
  x1_ *= delta;
  x1_ += x1;

  // Construction de la differentielle (dans la direction v0)
  construire_systeme(c, matrice_diffusion_CH, v1, x1_);
  construire_systeme(c, matrice_diffusion_CH, v2, x1);

  v1 -= v2;
  v1/= delta;

  return ;

}


// Description:
//     Algorithme GMRES Non Lineaire
// Precondition:
// Parametre: mutilde
//    Signification: inverse A x = b
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: int
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
// Modifie par DJ
//---------------
int Source_Con_Phase_field::non_lin_gmres(const DoubleTab& c, const DoubleTab& mutilde, const Matrice_Morse& matrice_diffusion_CH, DoubleTab& c_demi, DoubleTab& mutilde_demi)
// int Source_Con_Phase_field::non_lin_gmres(DoubleTab& c, DoubleTab& mutilde, Matrice_Morse& matrice_diffusion_CH, DoubleTab& x1)
//---------------
/*--------------------------------------------------------------
  Solves Ax = b by GMRES with nit iterations (reinitializations)
  and nkr Krilov vectors.
*/

// Copied from Sch_Crank_Nicholson.cpp (scalar)//
{
  int i,j,nk,i0,im,it,ii;
  double tem=1.,res,ccos,ssin ;
  const int ns = 2*c.size_totale();

  // Ajoute par DJ
  //--------------
  int nb_elem_tot = c.size_totale();
  DoubleTab x1(ns);
  x1=0.;
  for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
    {
      x1(n_elem)=c(n_elem);
      x1(n_elem+nb_elem_tot)=mutilde(n_elem);
    }
  //--------------

  // A present dans le jdd
  //   double epsGMRES=1.e-10;
  //   int nkr=2;                         // dimension de l'espace de Krylov
  //   int nit=10;                        // nombre d'iterations
  //   double rec_min = 1.e-8;
  //   double rec_max = 0.1  ;

  DoubleTab v(ns,nkr);                         // Krilov vectors
  DoubleTab h(nkr+1,nkr);                // Heisenberg maatrix of coefficients
  DoubleVect r(nkr+1);
  DoubleTab v0(x1);
  DoubleTab v1(x1) ;

  // Initialisation
  v = 0. ;
  v0 = 0. ;
  v1 = 0. ;

  // v0 = -1.*construire_systeme(eqn, x1); // DJ

  // Cerr << " gmres : avant construire systeme " << finl ;
  construire_systeme(c, matrice_diffusion_CH, v0, x1) ;
  //   // Ajoute par DJ pour Debog
  //   //-------------------------
  //   {
  //     DoubleTab v0_c(c);
  //     DoubleTab v0_mutilde(c);
  //     for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
  //       {
  //         v0_c(n_elem) = v0(n_elem);
  //         v0_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
  //       }
  //     Debog::verifier("GMRES Non Lineaire v0_c apres construire_systeme initial : ",v0_c);
  //     Debog::verifier("GMRES Non Lineaire v0_mutilde apres construire_systeme initial : ",v0_mutilde);
  //   }
  //   //-------------------------
  v0 *= -1. ;
  // Cerr << " gmres : apres construire systeme " << finl ;
  //   // Ajoute par DJ pour Debog
  //   //-------------------------
  //   {
  //     DoubleTab x1_c(c);
  //     DoubleTab x1_mutilde(c);
  //     DoubleTab systeme_c(c);
  //     DoubleTab systeme_mutilde(c);
  //     for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
  //       {
  //         x1_c(n_elem) = x1(n_elem);
  //         x1_mutilde(n_elem) = x1(n_elem+nb_elem_tot);
  //         systeme_c(n_elem) = v0(n_elem);
  //         systeme_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
  //       }
  //     Debog::verifier("GMRES Non Lineaire construire_systeme x1_c : ", x1_c);
  //     Debog::verifier("GMRES Non Lineaire construire_systeme x1_mutilde : ", x1_mutilde);
  //     Debog::verifier("GMRES Non Lineaire construire_systeme systeme_c : ", systeme_c);
  //     Debog::verifier("GMRES Non Lineaire construire_systeme systeme_mutilde : ", systeme_mutilde);
  //   }
  //   //-------------------------


  res = 0. ;
  // Modifie par DJ
  //---------------
  //   for (ii=0; ii<ns;ii++) res += v0(ii) * v0(ii) ;
  for (ii=0; ii<c.size(); ii++) res += v0(ii) * v0(ii) ;
  for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) res += v0(ii) * v0(ii) ;
  //---------------
  res=mp_sum(res);
  //   Debog::verifier("GMRES Non Lineaire res debut : ",res);
  res = sqrt(res)  ;

  //Cerr<<"initial residual = "<<res<<finl;

  if(res<rec_min)
    return 0; // nothing to do

  rec_min = (rec_min<res*epsGMRES) ? res*epsGMRES : rec_min;
  rec_min = (rec_min<rec_max) ? rec_min : rec_max ;

  Cerr << "Source Concentration Phase Field - GMRES NL" << finl;
  Cerr << "Stopping rule scalar : " << rec_min << finl;

  // iterations
  for(it=0; it<nit; it++)
    {
      nk = nkr;

      //...Orthogonalisation of Arnoldi
      v0 /= res;
      r = 0. ;
      r[0] = res;
      h=0.;

      for(j=0; j<nkr; j++)
        {
          for(ii=0; ii<ns; ii++) v(ii,j) = v0(ii);
          //           // Ajoute par DJ pour Debog
          //           //-------------------------
          //           {
          //             DoubleTab v_c(c);
          //             DoubleTab v_mutilde(c);
          //             for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
          //               {
          //                 v_c(n_elem) = v(n_elem,j);
          //                 v_mutilde(n_elem) = v(n_elem+nb_elem_tot,j);
          //               }
          //             Debog::verifier("GMRES Non Lineaire v_c : ",v_c);
          //             Debog::verifier("GMRES Non Lineaire v_mutilde : ",v_mutilde);
          //           }
          //           //-------------------------
          //          v0 = a * v0; // commente par DJ
          //      Cerr << " x1 avant matvec " << x1 << finl ;
          //           MatVect(...); // Modif pour le PhF
          matvect(c, matrice_diffusion_CH, v0, x1, v1) ;
          //           // Ajoute par DJ pour Debog
          //           //-------------------------
          //           {
          //             DoubleTab v1_c(c);
          //             DoubleTab v1_mutilde(c);
          //             for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
          //               {
          //                 v1_c(n_elem) = v1(n_elem);
          //                 v1_mutilde(n_elem) = v1(n_elem+nb_elem_tot);
          //               }
          //             Debog::verifier("GMRES Non Lineaire v1_c apres matvect : ",v1_c);
          //             Debog::verifier("GMRES Non Lineaire v1_mutilde apres matvect : ",v1_mutilde);
          //           }
          //           //-------------------------
          v0 = v1 ;

          // Modifie par DJ
          //---------------
          for(i=0; i<=j; i++)
            {
              //               for (ii=0; ii<ns;ii++) h(i,j) += v0(ii) * v(ii,i);
              for (ii=0; ii<c.size(); ii++) h(i,j) += v0(ii) * v(ii,i);
              for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) h(i,j) += v0(ii) * v(ii,i);
              h(i,j)=mp_sum(h(i,j));
              //               Debog::verifier("GMRES Non Lineaire h(i,j) : ",h(i,j));
              for (ii=0; ii<ns; ii++) v0(ii) -= h(i,j) * v(ii,i);
              //               {
              //                 DoubleTab v0_c(c);
              //                 DoubleTab v0_mutilde(c);
              //                 for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
              //                   {
              //                     v0_c(n_elem) = v0(n_elem);
              //                     v0_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
              //                   }
              //                 Debog::verifier("GMRES Non Lineaire v0_c apres modif : ",v0_c);
              //                 Debog::verifier("GMRES Non Lineaire v0_mutilde apres modif : ",v0_mutilde);
              //               }
            }
          //---------------
          // tem = sqrt(v0 * v0);
          tem = 0. ;
          // Modifie par DJ
          //---------------
          //           for (ii=0; ii<ns;ii++) tem += v0(ii) * v0(ii) ;
          for (ii=0; ii<c.size(); ii++) tem += v0(ii) * v0(ii) ;
          for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) tem += v0(ii) * v0(ii) ;
          //---------------
          tem=mp_sum(tem);
          //           Debog::verifier("GMRES Non Lineaire tem apres mp_sum : ",tem);
          tem = sqrt(tem)  ;
          h(j+1,j) = tem;
          if(tem<rec_min)
            {
              nk = j+1;
              // cerr<<"tem="<<tem<<" nk="<<nk<<endl;
              goto l5;
            }
          v0 /= tem;
        }

      //...Triangularisation
l5:
      for(i=0; i<nk; i++)
        {
          im = i+1;
          tem = 1./sqrt(h(i,i)*h(i,i) + h(im,i)*h(im,i));
          ccos = h(i,i) * tem;
          ssin = - h(im,i) * tem;
          for(j=i; j<nk; j++)
            {
              tem = h(i,j);
              h(i,j) = ccos * tem - ssin * h(im,j);
              h(im,j) =  ssin * tem + ccos * h(im,j);
            }
          r[im] = ssin * r[i];
          r[i] *= ccos;
        }

      //...Solution of linear system
      for(i=nk-1; i>=0; i--)
        {
          r[i] /= h(i,i);
          for(i0=i-1; i0>=0; i0--)
            r[i0] -= h(i0,i)* r[i];
        }
      for(i=0; i<nk; i++)
        for(ii=0; ii<ns; ii++)  x1(ii) += r[i]*v(ii,i);
      //       // Ajoute par DJ pour Debog
      //       //-------------------------
      //       {
      //         DoubleTab x1_c(c);
      //         DoubleTab x1_mutilde(c);
      //         for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
      //           {
      //             x1_c(n_elem) = x1(n_elem);
      //             x1_mutilde(n_elem) = x1(n_elem+nb_elem_tot);
      //           }
      //         Debog::verifier("GMRES Non Lineaire x1_c apres modif : ", x1_c);
      //         Debog::verifier("GMRES Non Lineaire x1_mutilde apres modif : ", x1_mutilde);
      //       }
      //       //-------------------------
      //       x1.echange_espace_virtuel();

      //Cerr <<" futur in gemres : " << x1 << finl ;
      //New residual and stopping tests
      //v0 = -1.*construire_systeme(eqn, x1); //DJ
      //res = sqrt(v0 * v0);

      construire_systeme(c, matrice_diffusion_CH, v0, x1) ;
      v0 *= -1. ;

      // Cerr << "  v0 in gmres " << v0 << finl ;

      res = 0. ;
      // Modifie par DJ
      //---------------
      //   for (ii=0; ii<ns;ii++) res += v0(ii) * v0(ii) ;
      for (ii=0; ii<c.size(); ii++) res += v0(ii) * v0(ii) ;
      for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) res += v0(ii) * v0(ii) ;
      //---------------
      res=mp_sum(res);
      //       Debog::verifier("GMRES Non Lineaire res fin : ",res);
      res = sqrt(res)  ;

      Cerr<<" - At it = "<< it <<", residu scalar = "<< res << endl;

      if(res<rec_min)
        {
          // Ajoute par DJ
          //--------------
          for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
            {
              c_demi(n_elem)=x1(n_elem);
              mutilde_demi(n_elem)=x1(n_elem+nb_elem_tot);
            }
          // L'echange espace virtuel est fait dans premier_demi_dt()
          //--------------
          //         return 1;
          Cerr << "Number of iterations to reach convergence : " << it+1 << finl;
          Cerr << "" << finl;
          return it;
        }
      else if (it==nit-1)
        {
          // Ajoute par DJ
          //--------------
          // On fait le choix de mettre a jour c_demi et mutilde_demi meme s'il n'y a pas convergence...
          for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
            {
              c_demi(n_elem)=x1(n_elem);
              mutilde_demi(n_elem)=x1(n_elem+nb_elem_tot);
            }
          //--------------
          Cerr << "Stopped before convergence" << finl;
          Cerr << "" << finl;
        }
    }

  return -1;
}


// Description:
//     Calcul de mutilde au centre des elements
// Precondition:
// Parametre: mutilde
//    Signification: mutilde au centre des elements
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::calculer_mutilde(DoubleTab& mutilde) const
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const DoubleTab& div_alpha_gradC=eq_c.get_div_alpha_gradC();

  const Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));
  const DoubleTab rhoPF=eq_ns.rho().valeurs();

  // Calcul de mutilde

  mutilde = div_alpha_gradC;
  mutilde *= -1.;

  const int taille=mutilde.size();
  for (int i=0; i<taille; i++)
    {
      mutilde(i)+=beta*dpsidc(c(i));
      if(mutype_==1)
        {
          mutilde(i)+=(0.5*u_carre_(i))*drhodc_;
        }
    }
  // L'espace virtuel n'est pas a jour
  assert_invalide_items_non_calcules(mutilde);
}


// Description:
//     Calcul de u2 au centre des elements
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::calculer_u2_elem(DoubleVect& u_carre)
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& u=eq_ns.inconnue().valeurs();
  const DoubleTab& c=eq_c.inconnue().valeurs();

  DoubleVect u2;

  u_carre = c;
  u_carre = 0.;

  u2 = u;
  u2.carre(VECT_ALL_ITEMS);

  // u2 est calcule aux faces. mutildes, pression_thermo et c au centre des elements.
  // On doit donc interpoler u2 pour tout calculer au elements.

  // Interpolation de u2, sortie dans u_carre :
  // ------------------------------------------

  // Interpolation aux elements de u^2 (a partir de Discretisation_SG_VDF::calculer_tenS_capil)
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const IntTab& face_sommets = zone_VDF.face_sommets();
  const DoubleTab& positions = zone_VDF.xp();
  const Zone& zone_geom = zone_VDF.zone();
  const Domaine& dom=zone_geom.domaine();
  const int nb_elem = zone_VDF.nb_elem_tot();
  int nb_compo_;
  int f0,f1;
  int som0,som1;
  double psi,val0,val1;

  DoubleTab u2_elem(nb_elem,dimension);

  // Nombre de composantes du vecteur u2
  if(u2_elem.nb_dim()==1)
    nb_compo_ = 1;
  else
    nb_compo_ = u2_elem.dimension(1);

  // Boucle sur le nombre d'elements
  for(int elem=0; elem<nb_elem; elem++)
    {
      // Boucle sur le nombre de composantes du vecteur
      for(int ncomp=0; ncomp<nb_compo_; ncomp++)
        {
          f0 = elem_faces(elem,ncomp);
          f1 = elem_faces(elem,dimension+ncomp);

          val0 = u2(f0);
          val1 = u2(f1);

          som0 = face_sommets(f0,0);
          som1 = face_sommets(f1,0);

          psi = ( positions(elem,ncomp) - dom.coord(som0,ncomp) )
                / ( dom.coord(som1,ncomp) - dom.coord(som0,ncomp) ) ;

          if (dabs(psi) < 1.e-12)
            u2_elem(elem,ncomp) = val0 ;
          else if (dabs(1.-psi) < 1.e-12)
            u2_elem(elem,ncomp) = val1 ;
          else
            u2_elem(elem,ncomp) = val0 + psi * (val1-val0) ;
        }
    }

  // Pour chaque element, on additionne les composantes de u2_elem
  double tmp;
  for (int elem=0; elem<nb_elem; elem++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          tmp=u2_elem(elem,dim);
          u_carre(elem)+=tmp;
        }
    }

  u_carre.echange_espace_virtuel();

}


// Description:
//     Calcul de alpha*(Grad(C))^2 au centre des elements
// Precondition:
// Parametre: DoubleTab& alpha_gradC_carre
//    Signification: alpha*(Grad(C))^2 au centre des elements
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::calculer_alpha_gradC_carre(DoubleTab& alpha_gradC_carre) const
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;

  // On calcule Grad(c) que l'on met dans prov_face
  opgrad.calculer(c,prov_face);
  eq_ns.solv_masse().appliquer(prov_face);

  // On calcule (Grad(c))^2 sur les faces que l'on met dans gradc2
  DoubleVect gradc2;
  gradc2=prov_face;
  gradc2.carre();

  // Interpolation aux elements de (Grad(c))^2 (a partir de Discretisation_SG_VDF::calculer_tenS_capil)
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const IntTab& face_sommets = zone_VDF.face_sommets();
  const DoubleTab& positions = zone_VDF.xp();
  const Zone& zone_geom = zone_VDF.zone();
  const Domaine& dom=zone_geom.domaine();
  const int nb_elem = zone_VDF.nb_elem();
  int nb_compo_;
  int elem;
  int f0,f1;
  int som0,som1;
  double psi,val0,val1;

  DoubleTab gradc2_elem(nb_elem,dimension);

  // Nombre de composantes du vecteur gradc2
  if(gradc2_elem.nb_dim()==1)
    nb_compo_ = 1;
  else
    nb_compo_ = gradc2_elem.dimension(1);

  // Boucle sur le nombre d'elements
  for(elem=0; elem<nb_elem; elem++)
    {
      // Boucle sur le nombre de composantes du vecteur
      for(int ncomp=0; ncomp<nb_compo_; ncomp++)
        {
          f0 = elem_faces(elem,ncomp);
          f1 = elem_faces(elem,dimension+ncomp);

          val0 = gradc2(f0);
          val1 = gradc2(f1);

          som0 = face_sommets(f0,0);
          som1 = face_sommets(f1,0);

          psi = ( positions(elem,ncomp) - dom.coord(som0,ncomp) )
                / ( dom.coord(som1,ncomp) - dom.coord(som0,ncomp) ) ;

          if (dabs(psi) < 1.e-12)
            gradc2_elem(elem,ncomp) = val0 ;
          else if (dabs(1.-psi) < 1.e-12)
            gradc2_elem(elem,ncomp) = val1 ;
          else
            gradc2_elem(elem,ncomp) = val0 + psi * (val1-val0) ;
        }
    }

  // Pour chaque element, on additionne les composantes de gradc2_elem
  int dim;
  double tmp;
  for (elem=0; elem<nb_elem; elem++)
    {
      for (dim=0; dim<dimension; dim++)
        {
          tmp=gradc2_elem(elem,dim);
          alpha_gradC_carre(elem)+=tmp;
        }
    }

  alpha_gradC_carre *= alpha;
  // Ajout B.M suite a plantage dans assert_espace_virtuel_vect dans op_conv...
  alpha_gradC_carre.echange_espace_virtuel();
}



// Description:
//     Calcul de la pression thermodynamique aux elements
// Precondition:
// Parametre: DoubleTab& pression_thermo
//    Signification: pression thermodynamique aux elements
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: void
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::calculer_pression_thermo(DoubleTab& pression_thermo) const
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const DoubleTab& div_alpha_gradC=eq_c.get_div_alpha_gradC();

  // Recuperation de la pression en Pascal de l'etape de projection
  //---------------------------------------------------------------
  const DoubleTab& P_Pa=eq_ns.pression_pa().valeurs();

  const int taille = pression_thermo.size();
  for (int i=0; i<taille; i++)
    {
      pression_thermo(i) = P_Pa(i) - c(i)*div_alpha_gradC(i);
    }
}


// Description:
//     Derivee de rho par rapport a c = rho2-rho1 car on prend rho = rho1+(rho2-rho1)*c
// Precondition:
// Parametre: const doubleTab& rho
//    Signification: valeur de drho/dc
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Con_Phase_field::drhodc()
{
  drhodc_=rho2-rho1;
}

// Description:
//     Valeur de mu (potentiel chimique)
//     W=(1-c^2)c^2
//     Derivee de W par rapport a c = rho2-rho1 car on prend rho = rho1+(rho2-rho1)*c
// Precondition:
// Parametre: const double& c
//    Signification: concentration
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: double
//    Signification: valeur de dW/dc
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
double Source_Con_Phase_field::dpsidc(const double& c) const
{
  return (4.*c*(c+0.5)*(c-0.5));
  //   if (c<=-0.5 || c>=0.5)
  //     return(0.);
  //   else
  //     return (-1.*c);
}


