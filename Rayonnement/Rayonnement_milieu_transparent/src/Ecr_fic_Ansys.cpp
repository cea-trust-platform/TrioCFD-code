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
// File:        Ecr_fic_Ansys.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#include <Ecr_fic_Ansys.h>
#include <Zone.h>
#include <SFichier.h>

Implemente_instanciable(Ecr_fic_Ansys,"Ecrire_fichier_Ansys",Interprete);

Sortie& Ecr_fic_Ansys::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}
//
//
Entree& Ecr_fic_Ansys::readOn(Entree& is)
{
  return Interprete::readOn(is);
}
//
//
Entree& Ecr_fic_Ansys::interpreter(Entree& is)
{
  Cerr << "Debut Ecr_fic_Ansys::interpreter" << finl;
  Nom nom_dom, nom_fic;
  // Option faces non cachees, nombre de zones (precision), decoupage en teta si axi
  int non_hidden,zones,tetas;
  int i,j,k;
  is >> nom_dom >> nom_fic >> non_hidden >> zones >> tetas;
  if(! sub_type(Zone, objet(nom_dom)))
    {
      Cerr << nom_dom << " est du type " << objet(nom_dom).que_suis_je() << finl;
      Cerr << "On attendait un objet de type Zone" << finl;
      exit();
    }
  // On recupere le domaine a partir du nom de domaine
  const Zone& dom=ref_cast(Zone, objet(nom_dom));
  const Zone& zone=dom.zone(0);
  SFichier fic(nom_fic);
  // On ecrit l'en tete du fichier Ansys
  fic << "/BATCH" << finl;
  fic << "/PREP7" << finl;
  fic << "/TITLE,none" << finl;
  // Si NOPR NUMBER OF ELEMENTS=  n'est pas affiche dans le .vf
  // Faire NOPR puis GOPR ?
  // fic << "/NOPR" << finl;
  fic << "ANTYPE,STATIC" << finl;

  // Ecriture du systeme de coordonnees
  // C'est important de l'ecrire en 1er
  if (bidim_axi)
    fic << "CSYS,0" << finl;
  else if (axi)
    fic << "CSYS,1" << finl;
  else
    fic << "CSYS,0" << finl;
  fic << "/NOPR" << finl;
  // On ecrit tous les sommets
  DoubleVect xyz(dimension);
  int n=1;
  int z=1;
  if (non_hidden>1 && dimension==2)
    z=1;
  int nb_som=zone.nb_som();
  for (j=0; j<z; j++)
    {
      for (i=0; i<nb_som; i++)
        {
          for (k=0; k<dimension; k++)
            xyz(k)=dom.coord(i,k);
          // Dans Ansys, teta est en degres et non en radians comme dans TRUST
          if (axi)
            xyz(1)*=180/M_PI;
          fic<<"N,"<<n++;
          // On ecrit les coordonnees
          for (k=0; k<dimension; k++)
            {
              if (j==0)
                fic<<","<<xyz(k);
              else
                {

                }
            }
          fic << finl;
        }
    }
  Nom type=dom.zone(0).type_elem()->que_suis_je();
  if (non_hidden>1)
    {
      // Quel que soit le type et la dimension du probleme
      fic << "ET,1,SHELL57" << finl;
    }
  else
    {
      // Ecriture du type d'element pour la methode Radiation Matrix
      if (type=="Rectangle" || type=="Triangle" || type=="Hexaedre" || type=="Tetraedre" || type=="Quadrangle" || type=="Hexaedre_VEF")
        {
          if (dimension==2)
            fic << "ET,1,LINK32" << finl;
          else
            {
              fic << "ET,1,SURF152" << finl;
              fic << "KEYOPT,1,4,1" << finl;
              fic << "KEYOPT,1,5,0" << finl;
            }
        }
      else if (type=="Rectangle_axi" || type=="Hexaedre_axi")
        {
          if (dimension==2)
            fic << "ET,1,LINK32" << finl;
          else
            {
              fic << "ET,1,SURF152" << finl;
              fic << "KEYOPT,1,4,1" << finl;
              fic << "KEYOPT,1,5,0" << finl;
            }
        }
      else if (type=="Rectangle_2D_axi")
        {
          fic << "ET,1,LINK32" << finl;
        }
      else
        {
          Cerr << "Element " << type << " inexistant." << finl;
          exit();
        }
    }
  // Tableau de travail pour reordonner les sommets des faces de bord
  // pour que Ansys connaisse la face rayonnante (Voir Doc Ansys AUX12)
  IntVect elem(1);
  DoubleTab pos(1,dimension);
  double e=0.001;
  // Ecriture de la definition des faces de bord
  int num_globale=0;
  const DoubleTab& c=dom.coord_sommets();
  // Nombre de faces rayonnantes, total des bords et des raccords !
  int nombre_faces_rayonnantes=zone.nb_bords()+zone.nb_raccords();
  for (i=0; i<nombre_faces_rayonnantes; i++)
    {
      const IntTab& som=(i<zone.nb_bords()?zone.bord(i).faces().les_sommets():zone.raccord(i-zone.nb_bords()).valeur().faces().les_sommets());
      int nb_faces=som.dimension(0);
      for (j=0; j<nb_faces; j++)
        {
          if (dimension==2)
            {

              if (axi)
                {
                  double r0,t0,r1,t1;
                  r0=c(som(j,0),0);
                  t0=c(som(j,0),1);
                  r1=c(som(j,1),0);
                  t1=c(som(j,1),1);
                  if (t1<t0) t1+=2*M_PI;
                  if (est_egal(r0,r1))
                    {
                      pos(0,1)=0.5*(t0+t1);
                      if (t0<t1)
                        pos(0,0)=(1-e)*r0;
                      else
                        pos(0,0)=(1+e)*r0;
                    }
                  else if (est_egal(t0,t1))
                    {
                      pos(0,0)=0.5*(r0+r1);
                      if (r0<r1)
                        pos(0,1)=(1+e)*t0;
                      else
                        pos(0,1)=(1-e)*t0;
                    }
                }
              else
                {
                  double x0,y0,x1,y1;
                  x0=c(som(j,0),0);
                  y0=c(som(j,0),1);
                  x1=c(som(j,1),0);
                  y1=c(som(j,1),1);
                  double nx=y0-y1;
                  double ny=x1-x0;
                  double alpha=10*Objet_U::precision_geom/(sqrt(nx*nx+ny*ny));
                  pos(0,0)=0.5*(x0+x1)+alpha*nx;
                  pos(0,1)=0.5*(y0+y1)+alpha*ny;
                }
              dom.zone(0).chercher_elements(pos,elem);
              if (elem(0)!=-1)
                {
                  fic<<"EN,"<<num_globale+j+1<<","<<som(j,0)+1<<","<<som(j,1)+1;
                  if (non_hidden>1 && dimension==2)
                    fic<<",,";
                }
              else
                {
                  fic<<"EN,"<<num_globale+j+1<<","<<som(j,1)+1<<","<<som(j,0)+1;
                  if (non_hidden>1 && dimension==2)
                    fic<<",,";
                }
              fic<<finl;
            }
          else
            {
              if (axi)
                {
                  double r0,t0,z0,r3,t3,z3;
                  r0=c(som(j,0),0);
                  t0=c(som(j,0),1);
                  z0=c(som(j,0),2);
                  r3=c(som(j,3),0);
                  t3=c(som(j,3),1);
                  z3=c(som(j,3),2);
                  // Cas ou on repasse a teta=0;
                  if (t3<t0) t3+=2*M_PI;

                  if (est_egal(r0,r3))
                    {
                      pos(0,1)=0.5*(t0+t3);
                      if (t0<t3)
                        pos(0,0)=(1+e)*r0;
                      else
                        pos(0,0)=(1-e)*r0;
                      pos(0,2)=0.5*(z0+z3);
                    }
                  else if (est_egal(z0,z3))
                    {
                      pos(0,0)=0.5*(r0+r3);
                      pos(0,1)=0.5*(t0+t3);
                      if (r0<r3)
                        pos(0,2)=(1+e)*z0;
                      else
                        pos(0,2)=(1-e)*z0;
                    }
                  else if (est_egal(t0,t3))
                    {
                      pos(0,0)=0.5*(r0+r3);
                      pos(0,2)=0.5*(z0+z3);
                      if (r0<r3)
                        pos(0,1)=(1+e)*t0;
                      else
                        pos(0,1)=(1-e)*t0;
                    }
                  dom.zone(0).chercher_elements(pos,elem);
                  if (elem(0)!=-1)
                    fic<<"EN,"<<num_globale+j+1<<","<<som(j,0)+1<<","<<som(j,1)+1<<","<<som(j,3)+1<<","<<som(j,2)+1<<finl;
                  else
                    fic<<"EN,"<<num_globale+j+1<<","<<som(j,0)+1<<","<<som(j,2)+1<<","<<som(j,3)+1<<","<<som(j,1)+1<<finl;
                }
              else
                {
                  double x0,y0,z0,x1,y1,z1,x2,y2,z2;
                  x0=c(som(j,0),0);
                  y0=c(som(j,0),1);
                  z0=c(som(j,0),2);
                  x1=c(som(j,1),0);
                  y1=c(som(j,1),1);
                  z1=c(som(j,1),2);
                  x2=c(som(j,2),0);
                  y2=c(som(j,2),1);
                  z2=c(som(j,2),2);
                  double nx=(y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
                  double ny=(z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
                  double nz=(x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);
                  double alpha=100*Objet_U::precision_geom/(sqrt(nx*nx+ny*ny+nz*nz));
                  pos(0,0)=(x0+x1+x2)/3.+alpha*nx;
                  pos(0,1)=(y0+y1+y2)/3.+alpha*ny;
                  pos(0,2)=(z0+z1+z2)/3.+alpha*nz;
                  dom.zone(0).chercher_elements(pos,elem);
                  if (elem(0)!=-1)
                    {
                      if (type=="Tetraedre")
                        fic<<"EN,"<<num_globale+j+1<<","<<som(j,0)+1<<","<<som(j,1)+1<<","<<som(j,2)+1<<finl;
                      else
                        fic<<"EN,"<<num_globale+j+1<<","<<som(j,0)+1<<","<<som(j,1)+1<<","<<som(j,3)+1<<","<<som(j,2)+1<<finl;
                    }
                  else
                    {
                      if (type=="Tetraedre")
                        fic<<"EN,"<<num_globale+j+1<<","<<som(j,0)+1<<","<<som(j,2)+1<<","<<som(j,1)+1<<finl;
                      else
                        fic<<"EN,"<<num_globale+j+1<<","<<som(j,0)+1<<","<<som(j,2)+1<<","<<som(j,3)+1<<","<<som(j,1)+1<<finl;
                    }
                  if (elem(0)!=-1)
                    {
                      int el=elem(0);
                      pos(0,0)=(x0+x1+x2)/3.-alpha*nx;
                      pos(0,1)=(y0+y1+y2)/3.-alpha*ny;
                      pos(0,2)=(z0+z1+z2)/3.-alpha*nz;
                      dom.zone(0).chercher_elements(pos,elem);
                      if (elem(0)!=-1)
                        {
                          Cerr << "Cas non prevu dans l'algorithme de Ecr_fic_Ansys::interpreter ! " << finl;
                          Cerr << "La face " << j << " est entouree de deux elements " << el << " " << elem(0) << finl;
                          Cerr << alpha*nx << finl;
                          Cerr << alpha*ny << finl;
                          Cerr << alpha*nz << finl;
                          Cerr << "Contacter le support TRUST." <<finl;
                          exit();
                        }
                    }
                }
            }
        }
      num_globale+=nb_faces;
    }
  fic << "/GOPR" << finl;
  Cerr << "Nombre de faces: " << num_globale << finl;
  // Ecriture calcul des facteurs de forme
  if (non_hidden>1)
    {
      Cerr << "Methode Radiosity Solver" << finl;
      Cerr << "La RAM necessaire pour le calcul va etre de " << 16*num_globale*num_globale/1024/1024 << " Mo." << finl;
      Cerr << "26Mo <->1.8Mo" << finl;
      Cerr << "25Mo <->2.6Mo" << finl;
      fic << "SFE,ALL,2,RDSF,2,1" << finl;
      fic << "/AUX12" << finl;
    }
  else
    {
      Cerr << "Methode Radiation Matrix" << finl;
      Cerr << "La RAM necessaire pour le calcul va etre de " << 16*num_globale*num_globale/1024/1024 << " Mo." << finl;
      fic << "/AUX12" << finl;
      // A tester le nombre de zones sur la precision...
      if (non_hidden==0)
        fic << "VTYPE,0," << zones << finl;
      else if (non_hidden==1)
        fic << "VTYPE,1" << finl;
      // A tester le nombre de camenberts pour l'axi RZ (entre 6 et 90)...
      if (bidim_axi)
        fic << "GEOM,1," << tetas << finl;
      else if (dimension==2)
        fic << "GEOM,1" << finl;
      else
        fic << "GEOM,0" << finl;
      fic << "MPRINT,1" << finl;
    }
  Nom vf=nom_du_cas();
  vf+=".vf";
  if (non_hidden>1)
    {
      // On le fait en dimension 3 egalement pour eviter le message:
      // No Space Temperature or Space Node specified for open Enclosure 1.
      fic << "SPCTEMP,1,0.E+00" << finl;
      //if (dimension==2) fic << "SPCTEMP,1,0.E+00" << finl;
      fic << "HEMIOPT,"<<zones<< finl;
      fic << "TOFFST,100" << finl;
      fic << "VFOPT,NEW" << finl;
      if (non_hidden==3)
        {
          // Si vf_non_groupe rediriger vers /dev/null
          Nom vf_non_groupe(vf);
          vf_non_groupe+=".non_groupe";
          vf_non_groupe="trash";
          fic << "VFCALC," << vf_non_groupe << finl;
          fic << "/OUT,trash" << finl;
        }
      else
        {
          fic << "VFCALC," << vf << finl;
        }
      // Ecriture du regroupement des faces sur les bords
      if (non_hidden==3)
        {
          num_globale=0;
          for (i=0; i<nombre_faces_rayonnantes; i++)
            {
              const IntTab& som=(i<zone.nb_bords()?zone.bord(i).faces().les_sommets():zone.raccord(i-zone.nb_bords()).valeur().faces().les_sommets());
              int nb_faces=som.dimension(0);
              fic<<"ESEL,S,,,"<<num_globale+1<<","<<num_globale+nb_faces<<finl;
              fic<<"CM,B"<<i<<",ELEM"<<finl;
              fic<<"F"<<i<<"=0."<<finl;
              num_globale+=nb_faces;
            }
          for (i=0; i<nombre_faces_rayonnantes; i++)
            {
              for (j=0; j<nombre_faces_rayonnantes; j++)
                {
                  fic<<"VFQUERY,B"<<i<<",B"<<j<<finl;
                  fic<<"*GET,F"<<i<<j<<",RAD,,VFAVG"<<finl;
                  fic<<"F"<<i<<"=F"<<i<<"+F"<<i<<j<<finl;
                }
            }
          fic << "/OUT," << vf << finl;
          fic << "*VWRITE,1"<<finl;
          fic << "('Number of Enclosures = ',F3.1)"<<finl;
          fic << "*VWRITE,"<< nombre_faces_rayonnantes <<finl;
          fic << "('Facteurs de forme',F6.1,' Number of Surfaces = "<< nombre_faces_rayonnantes << "')";
          for (i=0; i<nombre_faces_rayonnantes; i++)
            {
              fic << finl;
              fic << "*VWRITE,F"<<i<<finl;
              fic << "('TOTAL= '1(F12.10,1X))"<<finl;
              fic << "*VWRITE";
              int colonne=0;
              for (j=0; j<nombre_faces_rayonnantes; j++)
                {
                  // On ecrit 10 colonnes par 10 colonnes car sinon Ansys tronque apres 17 colonnes
                  if (colonne == 10)
                    {
                      fic << finl;
                      fic << "("<<colonne<<"(F12.10,1X))"<<finl;
                      fic << "*VWRITE";
                      colonne=0;
                    }
                  colonne++;
                  fic << ",F" << i << j;
                }
              fic << finl;
              fic << "(" << colonne << "(F12.10,1X))";
            }
          fic << finl;
        }
    }
  else
    {
      fic << "/OUT," << vf << finl;
      fic << "WRITE,radiation_matrix" << finl;
    }
  fic << "FINISH" << finl;
  fic.close();
  Cerr << "Fin Ecr_fic_Ansys::interpreter" << finl;
  Cerr << "Lancement d'Ansys." << finl;
  Cerr << "Cela peut prendre du temps..." << finl;
  Cerr << (int) system("./.run_ansys") << finl;
  return is;
}
