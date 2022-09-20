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
// File      : Sonde_IJK.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <Sonde_IJK.h>

#include <IJK_Field.h>
#include <Ref_Probleme_FT_Disc_gen.h>
#include <Probleme_FT_Disc_gen.h>
#include <Domaine.h>
#include <Zone_VF.h>
#include <Champ_Generique_Interpolation.h>
#include <communications.h>
#include <IJK_FT.h>

// #include <string>
// #include <fstream>

static int fichier_sondes_cree = 0;

Implemente_instanciable( Sonde_IJK, "Sonde_IJK", Sonde ) ;

Sortie& Sonde_IJK::printOn(Sortie& s ) const
{
  return s << que_suis_je();
}

// Surcharge de la methode en remplacant la notion Pb.get_champ
// par
void Sonde_IJK::completer_IJK(const IJK_FT_double& ijk_ft)
{
  const Nom bidon("bidon");
  //On devrait acceder au domaine par le champ generique
  //Mais reference pas encore faite
  ref_ijk_ft_=ijk_ft;
  ref_ijk_field_ = ijk_ft.get_IJK_field(nom_champ_lu_);

  // A quoi sert ce bidon ? IJK compile meme en commentant les trois lignes suivantes
  const IJK_Splitting& splitting = ref_ijk_ft_.valeur().get_splitting_ft();
  post_bidon_.associer_nom_et_pb_base(bidon, ijk_ft.probleme(splitting));
  mon_post=post_bidon_;
  // Recherche du champ sonde
  // Remplissage de la reference au champ
  //

  // Remplissage de l'attribut ncomp (il vaut -1 par defaut)

  // const Noms nom_champ = mon_champ->get_property("nom");
  //const Noms noms_comp = mon_champ->get_property("composantes");
  ncomp =1 ;
  // Champ_Generique_base::composante(nom_champ_ref,nom_champ[0],noms_comp,mon_champ->get_property("synonyms"));

  initialiser();
}

Entree& Sonde_IJK::readOn( Entree& is )
{

  Motcle motlu;
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  int nbre_points;
  // initialisation de periode par defaut
  periode = 1.e10;
  // initialisation de periode-nodes-chsom-grav par defaut
  nodes=0;
  chsom=0;
  grav=0;
  som=0;
  dim=-1;
  numero_elem_=-1;
  is >> motlu;
  // Lecture d'un mot cle qui n'est pas le champ
  if (motlu=="nodes")
    {
      nodes=1;
      is >> motlu;
    }
  else if (motlu=="chsom")
    {
      chsom=1;
      is >> motlu;
    }
  else if (motlu=="grav")
    {
      grav=1;
      is >> motlu;
    }
  else if (motlu=="som")
    {
      som=1;
      is >> motlu;
    }
  // Affectation du nom du champ
  nom_champ_lu_ = motlu;

  //Creation des Champ_Generique_refChamp necessaire pour l initialisation de la REF a Champ_Generique_base
//Si le champ demande est un Champ_base connu du probleme on cree le Champ_Generique_refChamp correspondant



  // Lecture des caracteristiques de la sonde
  IntVect fait(2);

  Motcles les_motcles(12);
  {
    les_motcles[0] = "periode";
    les_motcles[1] = "point";
    les_motcles[2] = "points";
    les_motcles[3] = "segment";
    les_motcles[4] = "plan";
    les_motcles[5] = "volume";
    les_motcles[6] = "segmentxdx";
    les_motcles[7] = "planxdxdy";
    les_motcles[8] = "circle";
    les_motcles[9] = "position_like";
    les_motcles[10] = "numero_elem_sur_maitre";
    les_motcles[11] = "segmentpoints";
  }

  while ((fait(0) != 1) || (fait(1) != 1))
    {
      is >> motlu;
      if (motlu == accolade_fermee)
        {
          Cerr << "Error while reading the probe " << nom_ <<finl;
          Cerr << "The data of the probe were not defined" << finl;
          exit();
        }
      int rang=les_motcles.search(motlu);
      if (rang == -1)
        {
          Cerr << "Error while reading the probe " << nom_ <<finl;
          Cerr << motlu << " is not understood; the keywords understood are : " << finl;
          Cerr << les_motcles;
          exit();
        }

      switch(rang)
        {
        case 0:
          {
            is >> periode;
            break;
          }
        case 1:
        case 2:
        case 11:
          {
            type_ = les_motcles[rang];
            rang = 1;
            dim = 0;
            is >> nbre_points;
            les_positions_.resize(nbre_points,dimension);

            for (int i=0; i<nbre_points; i++)
              for (int j=0; j<dimension; j++)
                is >> les_positions_(i,j);

            break;
          }
        case 10:
          {
            type_ = les_motcles[rang];
            rang=1;
            dim=0;
            is >> numero_elem_;
            les_positions_.resize(1,dimension);
            break;
          }
        case 3:
        case 6:
          {
            std::cout<<"SONDE : CAS 3 "<<std::endl;
            std::cout<<"dimension "<<dimension<<std::endl;
            type_ = les_motcles[rang];
            int rang2=rang;
            rang = 1;
            dim = 1;
            DoubleVect origine(dimension);
            DoubleVect extremite(dimension);
            DoubleVect dx(dimension);
            int i=0,j=0;
            is >> nbre_points;
            std::cout<<"nbre_points "<<nbre_points<<std::endl;
            les_positions_.resize(nbre_points,dimension);

            for (; i<dimension; i++)
              {
                is >> origine(i);
                // std::cout<<"origine(i) "<<origine(i)<<std::endl;
                // std::cout<<"(i) "<<i<<std::endl;
                type_+=" ";
                type_+=(Nom)origine(i);
              }
            for (i=0; i<dimension; i++)
              {
                is >> extremite(i);
                // std::cout<<"extremite(i) "<<extremite(i)<<std::endl;
                type_+=" ";
                type_+=(Nom)extremite(i);
              }
            for (j=0; j<dimension; j++)
              if (rang2==6)
                dx(j)=(extremite(j))/(nbre_points-1);
              else
                dx(j)=(extremite(j)-origine(j))/(nbre_points-1);
            for (i=0; i<nbre_points; i++)
              {
                for (j=0; j<dimension; j++)
                  {
                    les_positions_(i,j)=origine(j)+i*dx(j);
                    // std::cout<<"dx(j) "<<dx(j)<<std::endl;
                  }
                // std::cout<<"les_positions_(i,0)"<<les_positions_(i,0)<<std::endl;
                // std::cout<<"les_positions_(i,1)"<<les_positions_(i,1)<<std::endl;
                // std::cout<<"les_positions_(i,2)"<<les_positions_(i,2)<<std::endl;
              }
            break;
          }
        case 4:
        case 7:
          {
            type_ = les_motcles[rang];
            int rang2=rang;
            rang = 1;
            dim = 2;
            DoubleVect origine(dimension);
            DoubleVect extremite1(dimension);
            DoubleVect extremite2(dimension);
            DoubleVect dx1(dimension);
            DoubleVect dx2(dimension);
            int i=0,j=0,k=0;
            is >> nbre_points1;
            is >> nbre_points2;
            nbre_points=nbre_points1*nbre_points2;
            les_positions_.resize(nbre_points,dimension);

            for (; i<dimension; i++)
              is >> origine(i);
            for (i=0; i<dimension; i++)
              is >> extremite1(i);
            for (i=0; i<dimension; i++)
              is >> extremite2(i);
            if (rang2==7)
              {
                for (i=0; i<dimension; i++)
                  dx1(i)=(extremite1(i))/(nbre_points1-1);
                for (i=0; i<dimension; i++)
                  dx2(i)=(extremite2(i))/(nbre_points2-1);
              }
            else
              {
                for (i=0; i<dimension; i++)
                  dx1(i)=(extremite1(i)-origine(i))/(nbre_points1-1);
                for (i=0; i<dimension; i++)
                  dx2(i)=(extremite2(i)-origine(i))/(nbre_points2-1);
              }
            for (i=0; i<nbre_points1; i++)
              for (j=0; j<nbre_points2; j++)
                for (k=0; k<dimension; k++)
                  les_positions_(i*nbre_points2+j,k)=origine(k)+i*dx1(k)+j*dx2(k);
            break;
          }

        case 5:
          {
            type_ = les_motcles[rang];
            rang = 1;
            dim = 3;
            ArrOfDouble origine(dimension);
            ArrOfDouble extremite1(dimension);
            ArrOfDouble extremite2(dimension);
            ArrOfDouble extremite3(dimension);
            ArrOfDouble dx1(dimension);
            ArrOfDouble dx2(dimension);
            ArrOfDouble dx3(dimension);
            int i=0,j=0,k=0;
            is >> nbre_points1;
            is >> nbre_points2;
            is >> nbre_points3;
            nbre_points=nbre_points1*nbre_points2*nbre_points3;
            les_positions_.resize(nbre_points,dimension);

            for (; i<dimension; i++)
              is >> origine[i];
            for (i=0; i<dimension; i++)
              is >> extremite1[i];
            for (i=0; i<dimension; i++)
              is >> extremite2[i];
            for (i=0; i<dimension; i++)
              is >> extremite3[i];
            for (i=0; i<dimension; i++)
              dx1[i]=(extremite1[i]-origine[i])/(nbre_points1-1);
            for (i=0; i<dimension; i++)
              dx2[i]=(extremite2[i]-origine[i])/(nbre_points2-1);
            for (i=0; i<dimension; i++)
              dx3[i]=(extremite3[i]-origine[i])/(nbre_points3-1);
            for (i=0; i<nbre_points1; i++)
              for (j=0; j<nbre_points2; j++)
                for (int m=0; m<nbre_points3; m++)
                  for (k=0; k<dimension; k++)
                    les_positions_(i+j*nbre_points1+m*nbre_points1*nbre_points2,k)=origine[k]+i*dx1[k]+j*dx2[k]+m*dx3[k];
            break;
          }
        case 8:
          {
            type_ = les_motcles[rang];
            // circle nbre_points x0 y0 [z0 dir] radius teta1 teta2
            rang = 1;
            dim = 1;
            int dir;
            double radius, teta1, teta2;
            DoubleVect origine(dimension);
            is >> nbre_points;
            les_positions_.resize(nbre_points,dimension);
            for (int i=0; i<dimension; i++)
              is >> origine(i);
            if (dimension==3) is >> dir;
            is >> radius >> teta1 >> teta2;
            // Ajout des informations
            for (int i=0; i<dimension; i++)
              {
                type_+=" ";
                type_+=(Nom)origine(i);
              }
            type_+=" ";
            type_+=(Nom)radius;
            type_+=" ";
            type_+=(Nom)teta1;
            type_+=" ";
            type_+=(Nom)teta2;
            // We calculate the positions
            for (int i=0; i<nbre_points; i++)
              {
                double angle=teta1+(teta2-teta1)*i/(nbre_points-1);
                angle*=M_PI/180;
                if (dimension==2)
                  {
                    les_positions_(i,0)=origine(0)+radius*cos(angle);
                    les_positions_(i,1)=origine(1)+radius*sin(angle);
                  }
                else if (dimension==3)
                  {
                    if (dir==0)
                      {
                        les_positions_(i,0)=origine(0);
                        les_positions_(i,1)=origine(1)+radius*cos(angle);
                        les_positions_(i,2)=origine(2)+radius*sin(angle);
                      }
                    else if (dir==1)
                      {
                        les_positions_(i,0)=origine(0)+radius*cos(angle);
                        les_positions_(i,1)=origine(1);
                        les_positions_(i,2)=origine(2)+radius*sin(angle);
                      }
                    else if (dir==2)
                      {
                        les_positions_(i,0)=origine(0)+radius*cos(angle);
                        les_positions_(i,1)=origine(1)+radius*sin(angle);
                        les_positions_(i,2)=origine(2);
                      }
                  }
              }
            break;
          }

        case 9:
          {
            /*
                  Motcle autre_sonde;
                  is >> autre_sonde;
                  // on cherche la sonde correspondante
                  int m=-1;
                  const Sondes& les_sondes=mon_post->les_sondes();
                  for (int i=0; i<les_sondes.size(); i++)
                    if (les_sondes(i).get_nom()==autre_sonde)
                      {
                        m=i;
                        break;
                      }
                  //else Cerr<<les_sondes(i).get_nom()<<finl;
                  if (m==-1)
                    {
                      Cerr<<" The probe name "<<autre_sonde<< " was not found"<<finl;
                      exit();
                    }
                  // on recupere  les_positions_
                  const Sonde& la_sonde_ref=les_sondes(m);
                  type_ = la_sonde_ref.get_type();
                  les_positions_=la_sonde_ref.les_positions();
                  dim =  la_sonde_ref.get_dim();
                  rang=1;
            */
            Cerr<<"pas dispo"<<finl;
            exit();
            break;
          }
        default:
          {
            Cerr << motlu <<"is not yet understood!" << finl;
            exit();
          }
        }
      fait(rang) = 1;
    }
  if ( (fait[0] == 0) || (fait[1] == 0) || (dim==-1))
    {
      Cerr << "Error while reading the probe " << nom_ << finl;
      Cerr << "The data of the probe have not been properly defined" << finl;
      exit();
    }

  // Construction du fichier associe a la sonde
  nom_fichier_=nom_du_cas();
  nom_fichier_+= "_";
  nom_fichier_+= nom_;
  if (dim==2)
    nom_fichier_+= ".plan";
  else
    nom_fichier_+= ".son";

  if (Process::je_suis_maitre())
    {
      // Le fichier nom_du_cas().son contient la liste de toutes les sondes du calcul en cours.
      // Il est reecrit pour chaque calcul au cas ou des sondes aient ete ajoutees ou supprimees.
      // Ajout du nom du fichier sonde dans le fichier listant les sondes
      // Ce fichier sera utilise par Run_sonde
      SFichier* fichier_sondes;
      Nom nom_fich=nom_du_cas();
      nom_fich+=".son";
      if (!fichier_sondes_cree)
        {
          fichier_sondes=new SFichier(nom_fich);
        }
      else
        {
          fichier_sondes=new SFichier(nom_fich,ios::app);
        }
      *fichier_sondes << nom_fichier_ << finl;
      fichier_sondes->close();
      delete fichier_sondes;
      fichier_sondes_cree = 1;
    }
  return is;
}


void Sonde_IJK::initialiser()
{
  nb_bip = 0.;
  int nbre_points = les_positions_.dimension(0);
  if(elem_.size() != nbre_points)
    elem_.resize(nbre_points);

  const IJK_Splitting& splitting = ref_ijk_ft_.valeur().get_splitting_ft();
  const Zone& zone_geom = ref_ijk_ft_.valeur().probleme(splitting).domaine().zone(0);
  if ( numero_elem_==-1)
    {
      int nb_som = zone_geom.type_elem().nb_som();
      int nb_coord = les_positions_.dimension(1);
      if (nb_coord+1>nb_som)
        {
          Cerr << "You can't specify the probe named " << nom_ << " with "<< nb_coord << " coordinates on the domain named " <<zone_geom.domaine().le_nom()<<finl;
          Cerr << "which is constituted with cells of kind " << zone_geom.type_elem().valeur().que_suis_je() << "." << finl;
          Cerr << "Change the probe coordinates or use numero_elem_sur_maitre keyword (see documentation)" << finl;
          Cerr << "to specify a cell containing the probe and not its coordinates." << finl;
          Process::exit();
        }
      zone_geom.chercher_elements(les_positions_,elem_,1);
    }
  else
    {
      Cerr<<"not implemented"<<finl;
      Process::exit();
      if (0)
        {
          elem_[0]= numero_elem_;
          const IntTab& les_elems=mon_champ->get_ref_zone_dis_base().zone().les_elems();
          const DoubleTab& coord=mon_champ->get_ref_zone_dis_base().zone().domaine().les_sommets();
          int nb_som=les_elems.dimension(1);
          if (numero_elem_<les_elems.dimension_tot(0))
            {
              for (int s=0; s<nb_som; s++)
                {
                  int soml=les_elems(numero_elem_,s);
                  for (int dir=0; dir<dimension; dir++)
                    les_positions_(0,dir)+=coord(soml,dir)/nb_som;
                }
            }
          else
            {
              if (je_suis_maitre())
                {
                  Cerr<<" The element number "<<numero_elem_<<" does not exist on the master processor "<<finl;
                  Cerr<<" we put the position to zero"<<finl;
                  les_positions_=0;
                }
              elem_[0]=-1;
            }
        }
    }
  for (int i=0; i<nbre_points; i++)
    if (mp_max(elem_[i])==-1 && je_suis_maitre())
      {
        Cerr << "WARNING: The point number " << i+1 << " of probe " << nom_ << " is outside of the computational domain." << finl;
      }
  if (je_suis_maitre()&&(nproc()>1))
    {
      if (ncomp == -1)
        {
          const Noms noms_comp = mon_champ->get_property("composantes");
          int nb_comp = noms_comp.size();
          if (nb_comp == 1)
            valeurs_sur_maitre.resize(nbre_points);
          else
            valeurs_sur_maitre.resize(nbre_points,nb_comp);
        }
      else
        valeurs_sur_maitre.resize(nbre_points);
    }
  // PQ : 07/10/04 : nodes=1 || grav=1 : relocalisation des points de sondes aux centres de gravite ou aux noeuds les plus proches

  const Zone& zone = zone_geom ; // mon_champ->get_ref_domain().zone(0);
  Noms nom_champ;
  nom_champ.dimensionner(1);
  nom_champ[0] =  nom_champ_lu_; //mon_champ->get_property("nom");
  Cerr << "Champ lu " << nom_champ_lu_ << " nodes= " << (int)nodes << " et grav= " << (int)grav << finl;
  if (grav==1)
    {
      Cerr<<"The location of probes associated to "<<nom_champ[0]<<" are modified (to centers of gravity):"<<finl;
      const Zone_VF& zoneVF = ref_cast(Zone_VF,mon_champ->get_ref_zone_dis_base());
      const DoubleTab xp = zoneVF.xp();
      for (int i=0; i<nbre_points; i++)
        {
          if(elem_[i]!=-1)
            {
              for (int dir=0; dir<dimension; dir++)
                {
                  Cerr << " x(" << dir << "): " << les_positions_(i,dir) << " -> " << xp(elem_[i],dir);
                  les_positions_(i,dir)=xp(elem_[i],dir);
                }
              Cerr << finl;
            }
        }
    }
  else if (nodes==1)
    {
      const Zone_VF& zoneVF = ref_cast(Zone_VF,mon_champ->get_ref_zone_dis_base());
      const DoubleTab xv = zoneVF.xv();
      const IntTab& elem_faces = zoneVF.elem_faces();
      if (mp_max(elem_faces.size_array())==0)
        {
          Cerr << "Error: the domain " << zoneVF.zone().domaine().le_nom() << " is not discretized." << finl;
          exit();
        }
      if (sub_type(Champ_Generique_Interpolation,mon_champ.valeur()))
        {
          Motcle dom_interp=mon_champ->get_ref_domain().le_nom();
          Cerr << finl;
          Cerr << "Error in your probe : " << nom_ << finl;
          Cerr << "You can not project to nodes, the field " << nom_champ[0] << finl;
          Cerr << "which is interpolated on the domain " << dom_interp << finl;
          exit();
        }
      Cerr<<"The location of probes associated to "<<nom_champ[0]<<" are modified (to faces):"<<finl;
      const int nfaces_par_element = zone.nb_faces_elem() ;
      for (int i=0; i<nbre_points; i++)
        {
          double dist_min=DMAXFLOAT;
          int face_min=-1;
          if(elem_[i]!=-1)
            {
              for(int fac=0; fac<nfaces_par_element; fac++)
                {
                  int face=elem_faces(elem_[i],fac);
                  double dist=0.;

                  for (int dir=0; dir<dimension; dir++)
                    dist+=(xv(face,dir)-les_positions_(i,dir))*(xv(face,dir)-les_positions_(i,dir));

                  if(dist<=dist_min)
                    {
                      dist_min=dist;
                      face_min=face;
                    }
                }

              for (int dir=0; dir<dimension; dir++)
                {
                  Cerr << " x(" << dir << "): " << les_positions_(i,dir) << " -> " << xv(face_min,dir);
                  les_positions_(i,dir)=xv(face_min,dir);
                }
              Cerr << finl;
            }
        }
    }
  else if (som==1)
    {
      if (sub_type(Champ_Generique_Interpolation,mon_champ.valeur()))
        {
          Motcle dom_interp=mon_champ->get_ref_domain().le_nom();
          Cerr << finl;
          Cerr << "Error in your probe : " << nom_ << finl;
          Cerr << "You can not project to vertexes, the field " << nom_champ[0] << finl;
          Cerr << "which is interpolated on the domain " << dom_interp << finl;
          exit();
        }
      Cerr<<"The location of probes associated to "<<nom_champ[0]<<" are modified (to vertexes):"<<finl;
      const IntTab& sommet_elem = zone.les_elems();
      const int sommets_par_element = zone.les_elems().dimension(1);
      const DoubleTab& coord = zone.domaine().les_sommets();
      for (int i=0; i<nbre_points; i++)
        {
          double dist_min=DMAXFLOAT;
          int sommet_min=-1;
          if(elem_[i]!=-1)
            {
              for(int isom=0; isom<sommets_par_element; isom++)
                {
                  int sommet=sommet_elem(elem_[i],isom);
                  double dist=0.;
                  for (int dir=0; dir<dimension; dir++)
                    dist+=(coord(sommet,dir)-les_positions_(i,dir))*(coord(sommet,dir)-les_positions_(i,dir));

                  if(dist<=dist_min)
                    {
                      dist_min=dist;
                      sommet_min=sommet;
                    }
                }

              for (int dir=0; dir<dimension; dir++)
                {
                  Cerr << " x(" << dir << "): " << les_positions_(i,dir) << " -> " << coord(sommet_min,dir);
                  les_positions_(i,dir)=coord(sommet_min,dir);
                }
              Cerr << finl;
            }
        }
    }

  // chaque processeur a regarder si il avait le point
  // le maitre construit un tableau (prop) determinant qui
  // va donner la valeur au maitre

  // GB Ajout : Corrige la position de la sonde.
  // Pour les coords demandees, on a trouver l'elem associe mais si le champ n'est pas aux elems,
  // il faut corriger sa localisation...
  // Je ne sais pas pourquoi, mais meme aux elems, la position n'est pas bonne...
  {
    IJK_Field_double& ijk_field = ref_ijk_field_.valeur();
    const IJK_Splitting::Localisation loc = ijk_field.get_localisation();

    const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
    Vecteur3 origin(0., 0., 0.), delta(0., 0., 0.);
    delta[0] = geom.get_constant_delta(DIRECTION_I);
    delta[1] = geom.get_constant_delta(DIRECTION_J);
    delta[2] = geom.get_constant_delta(DIRECTION_K);

    Int3 offset;
    offset[0] =  splitting.get_offset_local(DIRECTION_I);
    offset[1] =  splitting.get_offset_local(DIRECTION_J);
    offset[2] =  splitting.get_offset_local(DIRECTION_K);

    // L'origine est sur un noeud. Donc que la premiere face en I est sur get_origin(DIRECTION_I)
    origin[0] = geom.get_origin(DIRECTION_I)
                + ((loc==IJK_Splitting::FACES_J || loc==IJK_Splitting::FACES_K || loc==IJK_Splitting::ELEM) ? (delta[DIRECTION_I] * 0.5) : 0. ) ;
    origin[1] = geom.get_origin(DIRECTION_J)
                + ((loc==IJK_Splitting::FACES_K || loc==IJK_Splitting::FACES_I || loc==IJK_Splitting::ELEM) ? (delta[DIRECTION_J] * 0.5) : 0. ) ;
    origin[2] = geom.get_origin(DIRECTION_K)
                + ((loc==IJK_Splitting::FACES_I || loc==IJK_Splitting::FACES_J || loc==IJK_Splitting::ELEM) ? (delta[DIRECTION_K] * 0.5) : 0. ) ;

    for (int idx =0; idx < nbre_points; idx++)
      {
        const int num_elem = elem_[idx];
        if (num_elem != -1)
          {
            // L'element appartient a ce proc. On peut donc trouver son ijk...
            // Les autres positions seront mises a jour par les autres procs...
            const Int3 ijk = splitting.convert_packed_to_ijk_cell(num_elem);
            /*
              Cerr << "Avant : Sonde ?? Point " << idx
              << " x= " << les_positions_(idx, 0)
              << " y= " << les_positions_(idx, 1)
              << " z= " << les_positions_(idx, 2) << finl;
            */
            // Corrige la position de la sonde :
            for (int i = 0; i < 3; i++)
              {
                const double val = origin[i]+(ijk[i]+offset[i])*delta[i];
                /*
                // Meme aux elems, j'ai teste avec la pression, la sonde n'est pas placee a l'elem
                // Je ne sais pas pourquoi..
                if ((loc==IJK_Splitting::ELEM) && std::fabs(les_positions_(idx, i) - val)>delta[0]/10. ) {
                Cerr << "Error in Sondes_IJK::: Sonde in element with confusing position... "  << endl;
                Process::exit();
                }
                */
                les_positions_(idx, i) = val;
              }
            Cerr << "Sonde ?? Point " << idx
                 << " x= " << les_positions_(idx, 0)
                 << " y= " << les_positions_(idx, 1)
                 << " z= " << les_positions_(idx, 2) << finl;
          }
      }
  }

  // Le maitre construit aussi le vect(ArrOfInt) participant
  // lui donnant pour un proc les differents elements de celui-ci
  IntVect prop(elem_);
  if (je_suis_maitre())
    {
      IntVect elems2(elem_);

      prop=0;
      // Par defaut c'est le maitre le proprio (en particulier pour les sondes en dehors)
      IntVect elem_prov(elem_);
      DoubleTab positions_prov(les_positions_);
      int p;
      int nbproc=Process::nproc();
      for(p=1; p<nbproc; p++)
        {
          recevoir(elem_prov,p,0,2002+p);
          recevoir(positions_prov,p,0,2001+p);
          for (int el=0; el<nbre_points; el++)
            if (elem_prov[el]!=-1)
              {
                if  (elems2[el]==-1)
                  {
                    elems2[el]=elem_prov[el];
                    prop[el]=p;
                    for (int j=0; j<dimension; j++)
                      les_positions_(el,j)=positions_prov(el,j);
                  }
              }
        }
      // OK On a rempli le tableau prop;

      // le maitre dimensionne participant;
      participant.dimensionner(nbproc);
      for(p=0; p<nbproc; p++)
        {
          int size=0;
          for (int el=0; el<nbre_points; el++)
            if (prop[el]==p) size++;
          participant[p].resize_array(size);
          participant[p]=-1;
          int pos=0;
          for (int el=0; el<nbre_points; el++)
            if (prop[el]==p)
              {
                participant[p][pos]=el;
                pos++;
              }
          assert((size==0)||(min_array(participant[p])>-1));
        }
    }
  else
    {
      envoyer(elem_,Process::me(),0,2002+Process::me());
      envoyer(les_positions_,Process::me(),0,2001+Process::me());
    }

  envoyer_broadcast(prop, 0);

  reprise = ref_ijk_ft_.valeur().get_reprise();
  les_positions_sondes_=les_positions_;

  nbre_points_tot=nbre_points;

  //
  // on redimensionne les tableaux a la taille reel a partir de prop
  nbre_points=0;
  int el;
  for ( el=0; el<nbre_points_tot; el++)
    if (prop(el)==me()) nbre_points++;
  DoubleTab new_positions(nbre_points,dimension);
  IntVect new_elem(nbre_points);
  int comp=0;
  for ( el=0; el<nbre_points_tot; el++)
    if (prop(el)==me())
      {
        new_elem(comp)=elem_(el);
        for (int dir=0; dir<dimension; dir++)
          new_positions(comp,dir)=les_positions_(el,dir);
        comp++;
      }
  int test=nbre_points;
  test=mp_sum(test);
  assert(test==nbre_points_tot);
  elem_=new_elem;
  les_positions_=new_positions;
  //

  // on dimensionne le tableau valeurs_locales
  if (ncomp == -1)
    {
      const Noms noms_comp = mon_champ->get_property("composantes");
      int nb_comp = noms_comp.size();

      if (nb_comp == 1)
        valeurs_locales.resize(nbre_points);
      else
        valeurs_locales.resize(nbre_points,nb_comp);
    }
  else
    valeurs_locales.resize(nbre_points);

}
void Sonde_IJK::postraiter()
{


  if (chsom==1)
    {
      abort();
    }
  else
    {
      // GB Ajout :
      //      valeurs_locales.resize(elem_.size());
      //      valeurs_locales=ref_ijk_field_.valeur()(0,0,0);
      {
        valeurs_locales.resize(elem_.size());
        const int nb_pts = elem_.size();
        IJK_Field_double& ijk_field = ref_ijk_field_.valeur();
        const IJK_Splitting& splitting = ijk_field.get_splitting();
        // GUILLAUME, sondes. Question GAB : on inhibe juste cette condition en ajoutant le && 0
        if ((ref_ijk_ft_.valeur().get_splitting_ft() != splitting) &&
            (ref_ijk_ft_.valeur().get_splitting_extension() != 0) && 0)
          {
            Cerr << " Error in Sonde_IJK::postraiter() "
                 << " L'extension est non-nulle donc splitting_ et splitting_ft_ different."
                 << " Le champ demande n'est pas discretise sur le domaine etendu "
                 << "alors que la Zone_VF est construite sur ce domaine..."  << finl;
            Process::exit();
          }
        for (int idx =0; idx < nb_pts; idx++)
          {
            const int num_elem = elem_[idx];
            const Int3 ijk = splitting.convert_packed_to_ijk_cell(num_elem);
            valeurs_locales[idx] = ijk_field(ijk[0],ijk[1],ijk[2]);
          }
      }
    }

  //int i;
  int nb_compo=valeurs_locales.nb_dim();
  //  int nbre_points = les_positions_.dimension(0);

  // le maitre reconstruit le tableau valeurs
  // a partir des differents contributeurs

  if(je_suis_maitre())
    {
      ouvrir_fichier();
      fichier().ouvrir(nom_fichier_,ios::app);
      fichier().setf(ios::scientific);
      fichier().precision(8);
      int p;
      int nbproc = Process::nproc();
      DoubleTab valeurs_pe;


      DoubleTab& valeurs=(nbproc==1?valeurs_locales:valeurs_sur_maitre);
      if (nbproc==1)
        ;//valeurs=valeurs_locales;
      else
        {
          int nb_val=valeurs_locales.dimension(0);
          if(nb_compo==1)
            {
              for(int i=0; i<nb_val; i++)
                {
                  valeurs[participant[0][i]]=valeurs_locales(i);
                }
            }
          else
            {
              int nb_val2=valeurs_locales.dimension(1);
              int k;
              for(int i=0; i<nb_val; i++)
                for(k=0; k<nb_val2; k++)
                  {
                    valeurs(participant[0][i],k)=valeurs_locales(i,k);
                  }
            }
        }

      for(p=1; p<nbproc; p++)
        {
          // le message n'est envoye que si le proc participe
          if (participant[p].size_array()!=0)
            {
              recevoir(valeurs_pe,p,0,2002+p);
              if(nb_compo==1)
                {
                  int nb_val=valeurs_pe.dimension(0);
                  for(int i=0; i<nb_val; i++)
                    {
                      valeurs[participant[p][i]]=valeurs_pe(i);
                      //           val_max = std::max(std::fabs(valeurs(i)),std::fabs(valeurs_pe(i)));
                      //                   if(val_max==(std::fabs(valeurs_pe(i))))
                      //                     valeurs(i)=valeurs_pe(i);
                    }
                }
              else
                {
                  int nb_val1=valeurs_pe.dimension(0);
                  int nb_val2=valeurs_pe.dimension(1);
                  int k;
                  for(int i=0; i<nb_val1; i++)
                    for(k=0; k<nb_val2; k++)
                      {
                        valeurs(participant[p][i],k)=valeurs_pe(i,k);
                        //  val_max = std::max(std::fabs(valeurs(i,k)),std::fabs(valeurs_pe(i,k)));
                        //                     if(val_max==(std::fabs(valeurs_pe(i,k))))
                        //                       valeurs(i,k)=valeurs_pe(i,k);
                      }
                }
            }
        }

      double temps_courant = ref_ijk_ft_.valeur().get_current_time();

      if (dim==0 || dim==1)
        {
          fichier() << temps_courant;
          for(int i=0; i<valeurs.dimension(0); i++)
            if (nb_compo==2)
              for(int k=0; k<valeurs.dimension(1); k++)
                fichier() << " " << valeurs(i,k);
            else
              fichier() << " " << valeurs(i);
          fichier() << finl;
        }
      // Pour les sondes type plan, impression au format lml :
      // num_sommet comp1 [comp2] [comp3]
      // et dans la troisieme direction :
      else if (dim==2 || dim==3)
        {
          Nom nom_post;
          int nbre_points = les_positions_.dimension(0);
          nbre_points =nbre_points_tot;
          const Noms noms_comp = mon_champ->get_property("composantes");
          int nb_comp = noms_comp.size();
          const Noms unites = mon_champ->get_property("unites");
          if (ncomp==-1)
            {
              const Noms noms_champ = mon_champ->get_property("nom");
              nom_post = noms_champ[0];
            }
          else
            nom_post = noms_comp[ncomp];

          Nom nom_topologie("Topologie");
          nom_topologie += "_";
          nom_topologie += nom_;

          fichier() << "TEMPS " << temps_courant << "\n";
          fichier() << "CHAMPPOINT " << nom_post << " " << nom_topologie
                    << " " << temps_courant << "\n";
          fichier() << nom_post << " " << nb_comp << " " << unites[0] << "\n";

          int nbp=nbre_points;
          if (dim==2) nbp*=2;
          if (nb_comp>1)
            fichier() << "type1 " << nbp << "\n";
          else fichier() << "type0 " << nbp << "\n";
          //int i;
          for(int i=0; i<nbre_points; i++)
            {
              fichier() << i+1;
              if (nb_compo==2)
                for(int j=0; j<valeurs.dimension(1); j++)
                  fichier() << " " << valeurs(i,j);
              else
                fichier() << " " << valeurs(i);
              // Pour ne pas flusher :
              fichier() << "\n";
            }
          // Pour le 2D, on rajoute une direction
          if (dim==2)
            {
              for(int i=0; i<nbre_points; i++)
                {
                  fichier() << nbre_points+i+1;
                  if (nb_compo==2)
                    for(int j=0; j<valeurs.dimension(1); j++)
                      fichier() << " " << valeurs(i,j);
                  else
                    fichier() << " " << valeurs(i);
                  // Pour ne pas flusher :
                  fichier() << "\n";
                }
            }
        }
      fichier().flush();
      fichier().close();
    }
  else
    {
      // le processeur envoye un message que si il participe
      if (valeurs_locales.dimension(0)!=0)
        envoyer(valeurs_locales,Process::me(),0,2002+Process::me());
    }
}

void Sonde_IJK::mettre_a_jour(const double un_temps, const double dt)
{
  double nb;
  // Le *(1+Objet_U::precision_geom) est pour eviter des erreurs d'arrondi selon les machines
  if (periode<=dt)
    nb=nb_bip+1.;
  else
    modf(un_temps*(1+Objet_U::precision_geom)/periode, &nb);

  // On doit ecrire les sondes
  if (nb>nb_bip)
    {
      nb_bip=nb; // On met a jour l'attribut de classe.
      reprise=1;
      postraiter();
    }
}
