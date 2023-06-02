/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Convert_ICoCoTrioField.h>
#include <ICoCoTrioField.h>
#include <ICoCoMEDDoubleField.hxx>
#include <Domaine.h>
#include <Champ_Generique_base.h>
#include <Domaine_VF.h>
#include <PE_Groups.h>
#include <Comm_Group.h>
#include <Probleme_base.h>
#include <fstream>

void affecte_double_avec_doubletab(double** p, const ArrOfDouble& trio)
{
  *p=new double[trio.size_array()];
  memcpy(*p,trio.addr(),trio.size_array()*sizeof(double));
}

void affecte_int_avec_inttab(int** p, const ArrOfInt& trio)
{
  int sz=trio.size_array();
  *p=new int[sz];
  memcpy(*p,trio.addr(),sz*sizeof(int));
}

void build_triofield(const Champ_Generique_base& ch, ICoCo::TrioField& afield)
{
  //DEV ANTONIN LEPREVOST FOR COUPLING TRIO WITH EUROPLEXUX

  if(ch.le_nom().getString() == "neant" and ch.get_dimension() == 3) //Pourquoi neant pr le champ extraction ? -> PBM TRUST
    {
      Cerr << "Interpolation P0 -> P1 TRUST" << endl;
      build_triomesh(ch.get_ref_domaine_dis_base(), afield, 1, ch.get_localisation() == FACE);
    }
  else
    build_triomesh(ch.get_ref_domaine_dis_base(), afield, ch.get_localisation() == NODE, ch.get_localisation() == FACE);

  afield.setName(ch.le_nom().getString());
  afield._time1 = afield._time2 = ch.get_time(), afield._itnumber = 0;

  /* copie des valeurs du champ */
  afield._has_field_ownership = true;
  Champ espace_stockage;
  const Champ_base& champ_ecriture = ch.get_champ(espace_stockage);
  const DoubleTab& vals = champ_ecriture.valeurs();

  const Domaine_VF& dom_vf = ref_cast(Domaine_VF, ch.get_ref_domaine_dis_base()); //On aimerai recuperer la surface des faces et d autres info geo par Domaine_VF mais RIEN n est bien defini quand on tire le domaine dis d un champ par extraction -> PBM TRUST
  const Probleme_base& pb_base = ch.get_ref_pb_base();
  const Domaine& dom_main = pb_base.domaine();
  const Frontiere& frontiere_IFS = dom_main.frontiere(3);
  const Faces& faces_IFS = frontiere_IFS.faces();

  int nb_som = afield._nbnodes;
  int space_dim = afield._space_dim;
  int nb_face = faces_IFS.nb_faces();
  int nb_som_face = afield._nodes_per_elem;

  if(ch.le_nom().getString() == "neant" and ch.get_dimension() == 3) //Pourquoi neant pr le champ extraction ? -> PBM TRUST
    {
      afield._mesh_dim = 2; //Pourquoi mesh dim vaut 1 en 3D ? Encore pbm avec le champ extraction, si on modifie pas erreur dans le couplage -> PBM TRUST

      DoubleTab vals_som(nb_som, space_dim);
      vals_som = 0.;

      DoubleVect surface_nodes(nb_som);
      DoubleVect surface_faces(nb_face);
      faces_IFS.calculer_surfaces(surface_faces);

      DoubleVect force_tot_som(space_dim);
      DoubleVect force_tot_elem(space_dim);

			double surface_tot_elem = 0;
			double surface_tot_som = 0;
      for(int face = 0; face < nb_face; face++)
        {
					surface_tot_elem += surface_faces(face);
          for(int som_loc = 0; som_loc < nb_som_face; som_loc++)
            {
              int som = faces_IFS.sommet(face, som_loc);
              surface_nodes[som] += (1.0/3.0) * surface_faces(face);
              
							for(int k = 0; k < space_dim; k++)
                vals_som(som, k) +=  (1.0/3.0) * vals(face, k); //On ne multiplie pas par la surface_faces(face) car vals est deja multiplie par la surface -> vals(face) = int_E vals(E) dx = |E| vals(E)
            }
          for(int k = 0; k < space_dim; k++)
            force_tot_elem(k) += vals(face, k);
        }

      for(int som = 0; som < nb_som; som++)//On n a pas besoin de diviser par la surface de chaque noeud (surface_nodes) car on doit renvoyer l integral du champ et non le champ lui meme pour EPX 
        {
					surface_tot_som += surface_nodes(som);
          for(int k = 0; k < space_dim; k++)
            {
							force_tot_som(k) += vals_som(som, k);
            }
        }

      std::ofstream file_force_tot_som("force_fluid_interpolate_som.txt", ios::app);
      std::ofstream file_force_tot_elem("force_fluid_elem.txt", ios::app);

      if(file_force_tot_som)
        file_force_tot_som << ch.get_time() << " " << force_tot_som(0) << " " << force_tot_som(1) << " " << force_tot_som(2) << endl;

      if(file_force_tot_elem)
        file_force_tot_elem << ch.get_time() << " " << force_tot_elem(0) << " " << force_tot_elem(1) << " " << force_tot_elem(2) << endl;
			
			Cerr << "Surface total elem : " << surface_tot_elem << endl;
			Cerr << "Surface total som : " << surface_tot_som << endl;

      afield._nb_field_components = vals_som.nb_dim() > 1 ? vals_som.dimension(1) : 1;
      affecte_double_avec_doubletab(&afield._field, vals_som);

    }
  else
    {
      afield._nb_field_components = vals.nb_dim() > 1 ? vals.dimension(1) : 1;
      affecte_double_avec_doubletab(&afield._field, vals);
    }
}

void build_triofield(const Champ_base& ch, const Domaine_dis_base& dom_dis, ICoCo::TrioField& afield)
{
  build_triomesh(dom_dis, afield, 0, 0); // XXX seulement elem...

  afield.setName(ch.le_nom().getString());
  afield._time1 = afield._time2 = 0.0, afield._itnumber = 0;

  /* copie des valeurs du champ */
  afield._has_field_ownership = true;
  const DoubleTab& vals = ch.valeurs();
  afield._nb_field_components = vals.nb_dim() > 1 ? vals.dimension(1) : 1;
  affecte_double_avec_doubletab(&afield._field, vals);
}

void build_triomesh(const Domaine_dis_base& dom_dis, ICoCo::TrioField& afield, int type, int loc_faces)
{
  const Domaine_VF& zvf = ref_cast(Domaine_VF, dom_dis);
  const Domaine& dom = dom_dis.domaine();
  afield.clear();
  afield._type = type;

  /* tableau des sommets : copie de celui du domaine */
  const DoubleTab& coord = dom.les_sommets();
  afield._space_dim = dom.dimension;
  afield._nbnodes = coord.dimension(0);
  affecte_double_avec_doubletab(&afield._coords, coord);

  /* dimension des elements du domaine */
  Motcle type_elem_ = dom_dis.domaine().type_elem()->que_suis_je();
  Motcle type_elem(type_elem_);
  type_elem.prefix("_AXI");
  if (type_elem != Motcle(type_elem_))
    {
      if (type_elem == "QUADRILATERE_2D")
        type_elem = "SEGMENT_2D";
      if (type_elem == "RECTANGLE_2D")
        type_elem = "RECTANGLE";
    }
  if ((type_elem == "RECTANGLE") || (type_elem == "QUADRANGLE") || (type_elem == "TRIANGLE") || (type_elem == "TRIANGLE_3D") || (type_elem == "QUADRANGLE_3D") || (type_elem == "POLYGONE") || (type_elem == "POLYGONE_3D"))
    afield._mesh_dim=2;
  else if ((type_elem == "HEXAEDRE") || (type_elem == "HEXAEDRE_VEF") || (type_elem == "POLYEDRE") || (type_elem == "PRISME") || (type_elem == "TETRAEDRE"))
    afield._mesh_dim=3;
  else if ((type_elem == "SEGMENT_2D") || (type_elem == "SEGMENT"))
    afield._mesh_dim=1;
  else
    {
      Cerr << "build_triofield: " << type_elem<< " not coded" <<finl;
      Process::exit();
    }

  /* elements : ceux du domaine si le champ est aux sommets/elements, les faces si le champ est aux faces */
  if (loc_faces) afield._mesh_dim--;
  afield._nb_elems = loc_faces ? zvf.nb_faces() : dom_dis.domaine().nb_elem();
  if (loc_faces || type_elem != "POLYEDRE") //maillage de faces -> connectivity = face_sommets
    {
      const IntTab& conn = loc_faces ? dom_dis.face_sommets() : dom_dis.domaine().les_elems();
      //le seul moyen qu'on a d'eviter que des polygones soient pris pour des quadrilateres est d'avoir un tableau de connectivite de largeur > 4...
      afield._nodes_per_elem = std::max(conn.dimension(1), type_elem == "POLYGONE" || type_elem == "POLYGONE_3D"  || type_elem == "POLYEDRE"  ? (int) 5 : 0);
      afield._connectivity = new int[afield._nb_elems * afield._nodes_per_elem];
      for (int i = 0; i < afield._nb_elems; i++)
        for (int j = 0; j < afield._nodes_per_elem; j++)
          afield._connectivity[afield._nodes_per_elem * i + j] = j < conn.dimension(1) ? conn(i, j) : -1;
    }
  else //maillage de polyedres -> connectivite au format MEDCoupling, a faire a la main
    {
      afield._nodes_per_elem = std::max(zvf.elem_faces().dimension(1) * (zvf.face_sommets().dimension(1) + 1), (int) 9); //un -1 apres chaque face : au moins 9 pour eviter un papillonage
      int *p = afield._connectivity = new int[afield._nb_elems * afield._nodes_per_elem];
      for (int e = 0, f, s, i, j; e < afield._nb_elems; e++)
        {
          /* insertion de la connectivite de chaque face, suivie d'un -1 */
          for (i = 0; i < zvf.elem_faces().dimension(1) && (f = zvf.elem_faces(e, i)) >= 0; i++, *p = -1, p++)
            for (j = 0; j < zvf.face_sommets().dimension(1) && (s = zvf.face_sommets(f, j)) >= 0; j++) *p = s, p++;
          /* des -1 jusqu'a la ligne suivante */
          for ( ; p < afield._connectivity + (e + 1) * afield._nodes_per_elem; p++) *p = -1;
        }
    }
}

#ifndef NO_MEDFIELD
#include <MEDCouplingUMesh.hxx>
#include <MEDCouplingFieldDouble.hxx>
#include <MCAuto.hxx>
#include <string.h>



using ICoCo::TrioField;
using ICoCo::MEDDoubleField;
using std::vector;


/*!
 * This method is non const only due to this->_field that can be modified (to point to the same domaine than returned object).
 * So \b warning, to access to \a this->_field only when the returned object is alive.
 */
MEDDoubleField build_medfield(TrioField& triofield)
{
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> mesh(MEDCoupling::MEDCouplingUMesh::New("",triofield._mesh_dim));
  MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> coo(MEDCoupling::DataArrayDouble::New());
  coo->alloc(triofield._nbnodes,triofield._space_dim);
  mesh->setCoords(coo);
  double *ptr(coo->getPointer());
  std::copy(triofield._coords,triofield._coords+triofield._space_dim*triofield._nbnodes,ptr);
  mesh->allocateCells(triofield._nb_elems);
  INTERP_KERNEL::NormalizedCellType elemtype;
  switch(triofield._mesh_dim)
    {
    case 0 :
      {
        switch (triofield._nodes_per_elem)
          {
          case 0: // cas field vide
            elemtype=INTERP_KERNEL::NORM_SEG2; // pour eviter warning
            break;
          default:
            throw INTERP_KERNEL::Exception("incompatible Trio field - wrong nb of nodes per elem");
          }
        break;
      }
    case 1:
      {
        switch (triofield._nodes_per_elem)
          {
          case 2:
            elemtype=INTERP_KERNEL::NORM_SEG2;
            break;
          default:
            throw INTERP_KERNEL::Exception("incompatible Trio field - wrong nb of nodes per elem");
          }
        break;
      }
    case 2:
      {
        switch (triofield._nodes_per_elem)
          {
          case 3:
            elemtype=INTERP_KERNEL::NORM_TRI3;
            break;
          case 4 :
            elemtype=INTERP_KERNEL::NORM_QUAD4;
            break;
          default:
            elemtype=INTERP_KERNEL::NORM_POLYGON;
          }
        break;
      }
    case 3:
      {
        switch (triofield._nodes_per_elem)
          {
          case 4:
            elemtype=INTERP_KERNEL::NORM_TETRA4;
            break;
          case 8 :
            elemtype=INTERP_KERNEL::NORM_HEXA8;
            break;
          default:
            elemtype=INTERP_KERNEL::NORM_POLYHED;
          }
        break;
      default:
        throw INTERP_KERNEL::Exception("incompatible Trio field - wrong mesh dimension");
      }
    }
  //creating a connectivity table that complies to MED (1 indexing) <- en fait non
  //and passing it to _mesh
  MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> field;
  int *conn(new int[triofield._nodes_per_elem]);
  for (int i=0; i<triofield._nb_elems; i++)
    {

      for(int j=0; j<triofield._nodes_per_elem; j++)
        {
          conn[j]=triofield._connectivity[i*triofield._nodes_per_elem+j];
        }
      if (elemtype==INTERP_KERNEL::NORM_QUAD4)
        {
          // dans trio pas la meme numerotation
          int tmp=conn[3];
          conn[3]=conn[2];
          conn[2]=tmp;
        }
      if (elemtype==INTERP_KERNEL::NORM_HEXA8)
        {
          // dans trio pas la meme numerotation
          int tmp=conn[3];
          conn[3]=conn[2];
          conn[2]=tmp;
          tmp=conn[7];
          conn[7]=conn[6];
          conn[6]=tmp;
        }
      int size = triofield._nodes_per_elem;
      while (conn[size - 1] == -1) size--; //on enleve les -1 a la fin de la connectivite
      mesh->insertNextCell(elemtype,size,conn);
    }
  delete [] conn;
  mesh->finishInsertingCells();
  //


  std::vector<int> cells;
  // if ((mesh->getSpaceDimension() == 2 || mesh->getSpaceDimension() == 3) && mesh->getMeshDimension() == 2)
  //   mesh->checkButterflyCells(cells);
  if (!cells.empty())
    {
      Cerr<<" cells are butterflyed "<<cells[0]<<finl;
      PE_Groups::groupe_TRUST().abort();
      Process::abort();
      Process::exit();
    }
  //field on the sending end
  int nb_case=triofield.nb_values();
  if (triofield._type==0)
    {
      field =  MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
    }
  else
    {
      field =  MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES,MEDCoupling::ONE_TIME );
    }
  field->setMesh(mesh);
  field->setNature(MEDCoupling::IntensiveMaximum);
  MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> fieldArr(MEDCoupling::DataArrayDouble::New());
  fieldArr->alloc(field->getNumberOfTuplesExpected(),triofield._nb_field_components);
  field->setName(triofield.getName());
  std::string meshName("SupportOf_");
  meshName+=triofield.getName();
  mesh->setName(meshName);
  field->setTime(triofield._time1,0,triofield._itnumber);
  if (triofield._field!=0)
    {
      for (int i =0; i<nb_case; i++)
        for (int j=0; j<triofield._nb_field_components; j++)
          {
            fieldArr->setIJ(i,j,triofield._field[i*triofield._nb_field_components+j]);
          }
    }
  //field on the receiving end
  else
    {
      // the trio field points to the pointer inside the MED field
      triofield._field=fieldArr->getPointer();
      for (int i=0; i<triofield._nb_field_components*nb_case; i++)
        triofield._field[i]=0.0;
    }
  field->setArray(fieldArr);
  return MEDDoubleField(field);
}
MEDDoubleField build_medfield(const Champ_Generique_base& ch)
{
  TrioField fl;
  build_triofield(ch, fl);
  return build_medfield(fl);
}


#else
namespace ICoCo
{
class MEDDoubleField
{
};
}
ICoCo::MEDDoubleField build_medfield(ICoCo::TrioField& ch)
{
  Cerr<<"Version compiled without MEDCoupling"<<finl;
  Process::exit();
  throw;
}
ICoCo::MEDDoubleField build_medfield(const Champ_Generique_base& ch)
{
  Cerr<<"Version compiled without MEDCoupling"<<finl;
  Process::exit();
  throw;
}
#endif
