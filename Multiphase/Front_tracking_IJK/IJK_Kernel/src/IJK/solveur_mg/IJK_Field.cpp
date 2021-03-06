//TRUST_NO_INDENT
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
// File      : IJK_Field.cpp
// Directory : $IJK_ROOT/src/IJK/solveur_mg
//
/////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file IJK_Field.cpp.P
//
#include <IJK_Field.h>
#include <Comm_Group.h>
#include <PE_Groups.h>
#include <Statistiques.h>
#include <SFichier.h>
#include <TRUSTArray.h>
#include <communications.h>
#include <Comm_Group_MPI.h>
#include <simd_tools.h>
#include <IJK_Splitting.h>
#include <Noms.h>
#include <stat_counters.h>

void envoyer_recevoir(const void * send_buf, int send_buf_size, int send_proc,
		      void * recv_buf, int recv_buf_size, int recv_proc)
{
  const Comm_Group & grp = PE_Groups::current_group();
  if (!sub_type(Comm_Group_MPI, grp)) {
    if (send_proc == -1 && recv_proc == -1)
      return;
    Cerr << "Error in envoyer_recevoir: non empty message and not Comm_Group_MPI" << finl;
    Process::exit();
  }
  const Comm_Group_MPI & grpmpi = ref_cast(Comm_Group_MPI, grp);
  grpmpi.ptop_send_recv(send_buf, send_buf_size, send_proc,
			recv_buf, recv_buf_size, recv_proc);
}

Implemente_instanciable(K_Layer_Cells_Lists, "K_Layer_Cells_Lists", Objet_U);
Implemente_vect(K_Layer_Cells_Lists);

Sortie & K_Layer_Cells_Lists::printOn(Sortie & os) const { return os; }
Entree & K_Layer_Cells_Lists::readOn(Entree & is) { return is; }


Implemente_instanciable_sans_constructeur(ArrOfFloat_with_ghost, "ArrOfFloat_with_ghost", Objet_U);
Implemente_vect(ArrOfFloat_with_ghost);

Sortie & ArrOfFloat_with_ghost::printOn(Sortie & os) const { return os; }
Entree & ArrOfFloat_with_ghost::readOn(Entree & is) { return is; }

void ArrOfFloat_with_ghost::echange_espace_virtuel(int pe_min, int pe_max)
{
  statistiques().begin_count(echange_vect_counter_);
  // envoi a pe_max et reception de pe_min
  const int n = tab_.size_array() - ghost_ * 2;
  float *pdata = tab_.addr();
  envoyer_recevoir(pdata + n, ghost_ * sizeof(float), pe_max,
						 pdata, ghost_ * sizeof(float), pe_min);
  // l'autre
  envoyer_recevoir(pdata + ghost_, ghost_ * sizeof(float), pe_min,
						 pdata + n + ghost_, ghost_ * sizeof(float), pe_max);
  statistiques().end_count(echange_vect_counter_);
}

Implemente_instanciable_sans_constructeur(ArrOfDouble_with_ghost, "ArrOfDouble_with_ghost", Objet_U);
Implemente_vect(ArrOfDouble_with_ghost);

Sortie & ArrOfDouble_with_ghost::printOn(Sortie & os) const { return os; }
Entree & ArrOfDouble_with_ghost::readOn(Entree & is) { return is; }

void ArrOfDouble_with_ghost::echange_espace_virtuel(int pe_min, int pe_max)
{
  statistiques().begin_count(echange_vect_counter_);
  // envoi a pe_max et reception de pe_min
  const int n = tab_.size_array() - ghost_ * 2;
  double *pdata = tab_.addr();
  envoyer_recevoir(pdata + n, ghost_ * sizeof(double), pe_max,
						 pdata, ghost_ * sizeof(double), pe_min);
  // l'autre
  envoyer_recevoir(pdata + ghost_, ghost_ * sizeof(double), pe_min,
						 pdata + n + ghost_, ghost_ * sizeof(double), pe_max);
  statistiques().end_count(echange_vect_counter_);
}

void set_imin_nitot(int rank_dir, int nproc_dir, 
		    int pe_min, int pe_max, int n_local,
		    int & i_min, int & n_tot)
{
  int canal = 53;

  if (rank_dir == 0)
    i_min = 0;
  else
    recevoir(i_min, pe_min, canal);

  if (rank_dir < nproc_dir-1)
    envoyer(i_min + n_local, pe_max, canal);

  if (rank_dir == nproc_dir-1)
    n_tot = i_min + n_local;
  else
    recevoir(n_tot, pe_max, canal);
  if (rank_dir > 0)
    envoyer(n_tot, pe_min, canal);
}



Implemente_instanciable_sans_constructeur(IJK_Field_local_float, "IJK_Field_local_float", Objet_U);
Implemente_instanciable(IJK_Field_float, "IJK_Field_float", IJK_Field_local_float);
Implemente_vect(IJK_Field_float);

Sortie & IJK_Field_local_float::printOn(Sortie & os) const { return os; }
Entree & IJK_Field_local_float::readOn(Entree & is) { return is; }
Sortie & IJK_Field_float::printOn(Sortie & os) const { return os; }
Entree & IJK_Field_float::readOn(Entree & is) { return is; }

IJK_Field_local_float::IJK_Field_local_float()
{
  ni_ = nj_ = nk_ = ghost_size_ = 0;
  nb_compo_ = 1;
  j_stride_ = compo_stride_ = 0;
  offset_ = k_layer_shift_ = additional_k_layers_ = 0;
}


void IJK_Field_local_float::allocate(int Ni, int Nj, int Nk, int ghosts, int additional_k_layers, int nb_compos)
{
  if (ghosts > Ni || ghosts > Nj || ghosts > Nk) {
    Cerr << "Error in IJK_Field_local_float::allocate: ghostsize=" << ghosts
	 << " is larger than the local mesh size " << Ni << " " << Nj << " " << Nk << finl;
    Process::exit();
  }
  ni_ = Ni;
  nj_ = Nj;
  nk_ = Nk;
  ghost_size_ = ghosts;
  nb_compo_ = nb_compos;
  j_stride_ = ni_ + 2 * ghosts;
  while (j_stride_ % Simd_float::size() != 0)
    j_stride_++;
  compo_stride_ = j_stride_ * (nj_ + 2 * ghosts);
  offset_ = ghost_size_ * (1 + j_stride_ + compo_stride_ * nb_compo_);
  k_layer_shift_ = 0;
  additional_k_layers_ = additional_k_layers;
  int sz = (nk_ + 2 * ghosts + additional_k_layers) * compo_stride_ * nb_compo_ + 1;

  // Ensure that there is some padding data at beginning:
  offset_ += 8;

  // Align origin on cache line boundary
  const int CacheLineSizeBytes = 64;
  sz += CacheLineSizeBytes/sizeof(float);
  // Add supplemental data at end for the last SIMD instruction of the loop that might go
  // beyond the end of the usefull data
  sz += CacheLineSizeBytes/sizeof(float);
  data_.resize_array(sz);
  float *ptr = data_.addr() + offset_;
  char *cptr = (char *) ptr;
  long long iptr = (long long) cptr;
  long long decal = iptr & (CacheLineSizeBytes - 1);
  offset_ += (CacheLineSizeBytes - decal) / sizeof(float);
  
}

// Adds n * compo_stride_ * nb_compo_ to the offset
// (shifts the data by n layers in the k direction without moving memory)
// Used by the jacobi "in place" algorithm: the result replaces the source data but is shifted in k.
void IJK_Field_local_float::shift_k_origin(int n)
{
  k_layer_shift_ += n;
  // Check that the data still fits into the allocated memory block
  assert(k_layer_shift_ >= 0 && k_layer_shift_ <= additional_k_layers_);

  offset_ += compo_stride_ * nb_compo_ * n;
}

// Creates a field with nk=1, that points to the data than src(...,...,k_layer)
// (uses ArrOffloat::ref_array() to create the reference).
// ghostsize is identical to source array.
void IJK_Field_local_float::ref_ij(IJK_Field_local_float & src, int k_lay)
{
  Cerr << "Error: must implement ArrOfFloat::ref_array() and reset()" << endl;
  Process::exit();
}


void IJK_Field_float::exchange_data(int pe_send_, /* processor to send to */
			      int is, int js, int ks, /* ijk coordinates of first data to send */
			      int pe_recv_, /* processor to recv from */
			      int ir, int jr, int kr, /* ijk coordinates of first data to recv */
			      int isz, int jsz, int ksz) /* size of block data to send/recv */
{
  if (pe_send_ == Process::me() && pe_recv_ == Process::me()) {
    
    // Self (periodicity on same processor)
    float *dest = data().addr();
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++)
	  dest[linear_index(ir + i, jr + j, kr + k)] = operator()(is + i, js + j, ks + k);
    return;
  }
#if 0
  const int data_size = isz * jsz * ksz;
  const int type_size = sizeof(float);
  float *send_buffer = 0;
  float **send_buffers = &send_buffer;
  float *recv_buffer = 0;
  float **recv_buffers = &recv_buffer;
  
  if (pe_send_ >= 0) {
    send_list.resize_array(1);
    send_list[0] = pe_send_;
    send_size.resize_array(1);
    send_size[0] = data_size * type_size;
    send_buffer = new float[data_size];
  } else {
    send_list.resize_array(0);
    send_size.resize_array(0);
  }

  if (pe_recv_ >= 0) {
    recv_list.resize_array(1);
    recv_list[0] = pe_recv_;
    recv_size.resize_array(1);
    recv_size[0] = data_size * type_size;
    recv_buffer = new float[data_size];
  } else {
    recv_list.resize_array(0);
    recv_size.resize_array(0);
  }
  
  if (pe_send_ >= 0) {
    // Pack send data
    float *buf = send_buffer;
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++, buf++)
	  *buf = operator()(is + i, js + j, ks + k);
  }
  const Comm_Group & grp = PE_Groups::current_group();
  grp.send_recv_start(send_list, send_size, (char**)send_buffers,
		      recv_list, recv_size, (char**)recv_buffers,
		      );
  grp.send_recv_finish();
#else 
  const int data_size = isz * jsz * ksz;
  const int type_size = sizeof(float);
  float *send_buffer = 0;
  float *recv_buffer = 0;
  
  if (pe_send_ >= 0) {
    send_buffer = new float[data_size];
    // Pack send data
    float *buf = send_buffer;
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++, buf++)
	  *buf = operator()(is + i, js + j, ks + k);
  }
  if (pe_recv_ >= 0) {
    recv_buffer = new float[data_size];
  }

  envoyer_recevoir(send_buffer, data_size * type_size, pe_send_,
		   recv_buffer, data_size * type_size, pe_recv_);
#endif
  if (pe_recv_ >= 0) {
    // Unpack recv data
    float *buf = recv_buffer;
    float *dest = data().addr();
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++, buf++) {
	  dest[linear_index(ir + i, jr + j, kr + k)] = *buf;
	}
  }

  delete[] send_buffer;
  delete[] recv_buffer;
}

// Description: Exchange data over "ghost" number of cells.
// Precondition: ghost <= ghost_size_
void IJK_Field_float::echange_espace_virtuel(int le_ghost)
{
  statistiques().begin_count(echange_vect_counter_);
  assert(le_ghost <= ghost_size_);
  // Exchange in i direction real cells, 
  // then in j direction with ghost cells in i,
  // then in k direction, with ghost cells in i and j
  const IJK_Splitting & splitting = splitting_ref_.valeur();
  int pe_imin_ = splitting.get_neighbour_processor(0, 0);
  int pe_imax_ = splitting.get_neighbour_processor(1, 0);
  int pe_jmin_ = splitting.get_neighbour_processor(0, 1);
  int pe_jmax_ = splitting.get_neighbour_processor(1, 1);
  int pe_kmin_ = splitting.get_neighbour_processor(0, 2);
  int pe_kmax_ = splitting.get_neighbour_processor(1, 2);

  // send left layer of real cells to right layer of virtual cells
  exchange_data(pe_imin_, 0,   0, 0, 
		pe_imax_, ni_, 0, 0, 
		le_ghost, nj_, nk_); /* size of block data to send */

  // send right real cells to left virtual cells
  exchange_data(pe_imax_, ni_ - le_ghost, 0, 0,
					 pe_imin_,     - le_ghost, 0, 0,
					 le_ghost, nj_, nk_);

  exchange_data(pe_jmin_, - le_ghost, 0, 0,
					 pe_jmax_, - le_ghost, nj_, 0,
					 ni_ + 2*le_ghost, le_ghost, nk_);

  exchange_data(pe_jmax_, - le_ghost, nj_ - le_ghost, 0,
					 pe_jmin_, - le_ghost,     - le_ghost, 0,
					 ni_ + 2*le_ghost, le_ghost, nk_);

  exchange_data(pe_kmin_, - le_ghost, - le_ghost, 0,
					 pe_kmax_, - le_ghost, - le_ghost, nk_,
					 ni_ + 2*le_ghost, nj_ + 2*le_ghost, le_ghost);

  exchange_data(pe_kmax_, - le_ghost, - le_ghost, nk_ - le_ghost,
					 pe_kmin_, - le_ghost, - le_ghost,     - le_ghost,
					 ni_ + 2*le_ghost, nj_ + 2*le_ghost, le_ghost);
  statistiques().end_count(echange_vect_counter_);
}

// Initializes the field and allocates memory
// splitting: reference to the geometry of the IJK mesh and how the mesh is split on processors.
//   The field stores a reference to this IJK_Splitting object so do not delete it.
// loc: localisation of the field (elements, nodes, faces in direction i, j, or k)
//   The number of "real" items in each direction (returned by the ni(), nj() or nk() method) is obtained from
//   the IJK_Splitting object. Warning: on a processor that is in the middle of the mesh, the nodes on the
//   right of the rightmost real element are not real, they are virtual, values are copied from the neigbour
//   processor.
// ghost_size: number of ghost layers to allocate. When an exchange of ghost cells data is requested, a smaller
//   number of layers can be requested if all layers are not needed by the following operations
// additional_k_layers: allocates more layers of cells in the k direction for use with the
//   shift_k_origin() method (optimization trick used by the temporally blocked algorithms in the multigrid
//   solver
// nb_compo: number of components of the field. Warning: you cannot allocate the velocity field at faces
//   with nb_compo=3: numbers of faces in each direction differ.
//   Also, components are not grouped by node but stored by layers in k. nb_compo>1 is essentially used
//   in the multigrid solver to optimize memory accesses to the components of the matrix.
void IJK_Field_float::allocate(const IJK_Splitting & splitting, IJK_Splitting::Localisation loc, 
				int ghost_size, int additional_k_layers, int ncompo)
{
  const int ni_local = splitting.get_nb_items_local(loc, 0);
  const int nj_local = splitting.get_nb_items_local(loc, 1);
  const int nk_local = splitting.get_nb_items_local(loc, 2);
  IJK_Field_local_float::allocate(ni_local, nj_local, nk_local, ghost_size, additional_k_layers, ncompo);
  splitting_ref_ = splitting;
  localisation_ = loc;
}

double norme_ijk(const IJK_Field_float & residu)
{
  const int ni = residu.ni();
  const int nj = residu.nj();
  const int nk = residu.nk();
  double somme = 0.;
  for (int k = 0; k < nk; k++) {
    double partial1 = 0.;
    for (int j = 0; j < nj; j++) {
      double partial2 = 0.;
      for (int i = 0; i < ni; i++) {
	double x = residu(i,j,k);
	partial2 += x*x;
      }
      partial1 += partial2;
    }
    somme += partial1;
  }
  somme = Process::mp_sum(somme);
  return sqrt(somme);
}

float max_ijk(const IJK_Field_float& residu)
{
  const int ni = residu.ni();
  const int nj = residu.nj();
  const int nk = residu.nk();
  float res = -1.e30;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          float x = residu(i,j,k);
	res = std::max(x,res);
        }
  res = Process::mp_max(res);
  return res;
}

float prod_scal_ijk(const IJK_Field_float & x, const IJK_Field_float & y)
{
  const int ni = x.ni();
  const int nj = x.nj();
  const int nk = x.nk();
  float somme = 0.;
  const ArrOfFloat & xd = x.data();
  const ArrOfFloat & yd = y.data();
  for (int k = 0; k < nk; k++) {
    for (int j = 0; j < nj; j++) {
      for (int i = 0; i < ni; i++) {
	int index = x.linear_index(i,j,k);
	assert(index == y.linear_index(i,j,k));
	somme += xd[index] * yd[index];
      }
    }
  }
  somme = Process::mp_sum(somme);
  return somme;
}

double somme_ijk(const IJK_Field_float & residu)
{
  const int ni = residu.ni();
  const int nj = residu.nj();
  const int nk = residu.nk();
  double somme = 0.;
  for (int k = 0; k < nk; k++) {
    double partial1 = 0.;
    for (int j = 0; j < nj; j++) {
      double partial2 = 0.;
      for (int i = 0; i < ni; i++) {
	double x = residu(i,j,k);
	partial2 += x;
      }
      partial1 += partial2;
    }
    somme += partial1;
  }
  somme = Process::mp_sum(somme);
  return somme;
}


Implemente_instanciable_sans_constructeur(IJK_Field_local_double, "IJK_Field_local_double", Objet_U);
Implemente_instanciable(IJK_Field_double, "IJK_Field_double", IJK_Field_local_double);
Implemente_vect(IJK_Field_double);

Sortie & IJK_Field_local_double::printOn(Sortie & os) const { return os; }
Entree & IJK_Field_local_double::readOn(Entree & is) { return is; }
Sortie & IJK_Field_double::printOn(Sortie & os) const { return os; }
Entree & IJK_Field_double::readOn(Entree & is) { return is; }

IJK_Field_local_double::IJK_Field_local_double()
{
  ni_ = nj_ = nk_ = ghost_size_ = 0;
  nb_compo_ = 1;
  j_stride_ = compo_stride_ = 0;
  offset_ = k_layer_shift_ = additional_k_layers_ = 0;
}


void IJK_Field_local_double::allocate(int Ni, int Nj, int Nk, int ghosts, int additional_k_layers, int nb_compos)
{
  if (ghosts > Ni || ghosts > Nj || ghosts > Nk) {
    Cerr << "Error in IJK_Field_local_double::allocate: ghostsize=" << ghosts
	 << " is larger than the local mesh size " << Ni << " " << Nj << " " << Nk << finl;
    Process::exit();
  }
  ni_ = Ni;
  nj_ = Nj;
  nk_ = Nk;
  ghost_size_ = ghosts;
  nb_compo_ = nb_compos;
  j_stride_ = ni_ + 2 * ghosts;
  while (j_stride_ % Simd_double::size() != 0)
    j_stride_++;
  compo_stride_ = j_stride_ * (nj_ + 2 * ghosts);
  offset_ = ghost_size_ * (1 + j_stride_ + compo_stride_ * nb_compo_);
  k_layer_shift_ = 0;
  additional_k_layers_ = additional_k_layers;
  int sz = (nk_ + 2 * ghosts + additional_k_layers) * compo_stride_ * nb_compo_ + 1;

  // Ensure that there is some padding data at beginning:
  offset_ += 8;

  // Align origin on cache line boundary
  const int CacheLineSizeBytes = 64;
  sz += CacheLineSizeBytes/sizeof(double);
  // Add supplemental data at end for the last SIMD instruction of the loop that might go
  // beyond the end of the usefull data
  sz += CacheLineSizeBytes/sizeof(double);
  data_.resize_array(sz);
  double *ptr = data_.addr() + offset_;
  char *cptr = (char *) ptr;
  long long iptr = (long long) cptr;
  long long decal = iptr & (CacheLineSizeBytes - 1);
  offset_ += (CacheLineSizeBytes - decal) / sizeof(double);
  
}

// Adds n * compo_stride_ * nb_compo_ to the offset
// (shifts the data by n layers in the k direction without moving memory)
// Used by the jacobi "in place" algorithm: the result replaces the source data but is shifted in k.
void IJK_Field_local_double::shift_k_origin(int n)
{
  k_layer_shift_ += n;
  // Check that the data still fits into the allocated memory block
  assert(k_layer_shift_ >= 0 && k_layer_shift_ <= additional_k_layers_);

  offset_ += compo_stride_ * nb_compo_ * n;
}

// Creates a field with nk=1, that points to the data than src(...,...,k_layer)
// (uses ArrOfdouble::ref_array() to create the reference).
// ghostsize is identical to source array.
void IJK_Field_local_double::ref_ij(IJK_Field_local_double & src, int k_lay)
{
  data_.reset();
  ni_ = src.ni_;
  nj_ = src.nj_;
  nk_ = 1;
  ghost_size_ = src.ghost_size_;
  j_stride_ = src.j_stride_;
  nb_compo_ = 1;
  compo_stride_ = 0;
  k_layer_shift_ = 0;
  additional_k_layers_ = 0;

  const int i_first = src.linear_index(-ghost_size_, -ghost_size_, k_lay);
  const int i_last = src.linear_index(ni_ + ghost_size_ - 1, nj_ + ghost_size_ - 1, k_lay);
  offset_ = src.linear_index(0,0,k_lay) - i_first;
  data_.ref_array(src.data_, i_first, i_last - i_first + 1);
}


void IJK_Field_double::exchange_data(int pe_send_, /* processor to send to */
			      int is, int js, int ks, /* ijk coordinates of first data to send */
			      int pe_recv_, /* processor to recv from */
			      int ir, int jr, int kr, /* ijk coordinates of first data to recv */
			      int isz, int jsz, int ksz) /* size of block data to send/recv */
{
  if (pe_send_ == Process::me() && pe_recv_ == Process::me()) {
    
    // Self (periodicity on same processor)
    double *dest = data().addr();
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++)
	  dest[linear_index(ir + i, jr + j, kr + k)] = operator()(is + i, js + j, ks + k);
    return;
  }
#if 0
  const int data_size = isz * jsz * ksz;
  const int type_size = sizeof(double);
  double *send_buffer = 0;
  double **send_buffers = &send_buffer;
  double *recv_buffer = 0;
  double **recv_buffers = &recv_buffer;
  
  if (pe_send_ >= 0) {
    send_list.resize_array(1);
    send_list[0] = pe_send_;
    send_size.resize_array(1);
    send_size[0] = data_size * type_size;
    send_buffer = new double[data_size];
  } else {
    send_list.resize_array(0);
    send_size.resize_array(0);
  }

  if (pe_recv_ >= 0) {
    recv_list.resize_array(1);
    recv_list[0] = pe_recv_;
    recv_size.resize_array(1);
    recv_size[0] = data_size * type_size;
    recv_buffer = new double[data_size];
  } else {
    recv_list.resize_array(0);
    recv_size.resize_array(0);
  }
  
  if (pe_send_ >= 0) {
    // Pack send data
    double *buf = send_buffer;
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++, buf++)
	  *buf = operator()(is + i, js + j, ks + k);
  }
  const Comm_Group & grp = PE_Groups::current_group();
  grp.send_recv_start(send_list, send_size, (char**)send_buffers,
		      recv_list, recv_size, (char**)recv_buffers,
		      Comm_Group::DOUBLE
		      );
  grp.send_recv_finish();
#else 
  const int data_size = isz * jsz * ksz;
  const int type_size = sizeof(double);
  double *send_buffer = 0;
  double *recv_buffer = 0;
  
  if (pe_send_ >= 0) {
    send_buffer = new double[data_size];
    // Pack send data
    double *buf = send_buffer;
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++, buf++)
	  *buf = operator()(is + i, js + j, ks + k);
  }
  if (pe_recv_ >= 0) {
    recv_buffer = new double[data_size];
  }

  envoyer_recevoir(send_buffer, data_size * type_size, pe_send_,
		   recv_buffer, data_size * type_size, pe_recv_);
#endif
  if (pe_recv_ >= 0) {
    // Unpack recv data
    double *buf = recv_buffer;
    double *dest = data().addr();
    for (int k = 0; k < ksz; k++)
      for (int j = 0; j < jsz; j++)
	for (int i = 0; i < isz; i++, buf++) {
	  dest[linear_index(ir + i, jr + j, kr + k)] = *buf;
	}
  }

  delete[] send_buffer;
  delete[] recv_buffer;
}

// Description: Exchange data over "ghost" number of cells.
// Precondition: ghost <= ghost_size_
void IJK_Field_double::echange_espace_virtuel(int le_ghost)
{
  statistiques().begin_count(echange_vect_counter_);
  assert(le_ghost <= ghost_size_);
  // Exchange in i direction real cells, 
  // then in j direction with ghost cells in i,
  // then in k direction, with ghost cells in i and j
  const IJK_Splitting & splitting = splitting_ref_.valeur();
  int pe_imin_ = splitting.get_neighbour_processor(0, 0);
  int pe_imax_ = splitting.get_neighbour_processor(1, 0);
  int pe_jmin_ = splitting.get_neighbour_processor(0, 1);
  int pe_jmax_ = splitting.get_neighbour_processor(1, 1);
  int pe_kmin_ = splitting.get_neighbour_processor(0, 2);
  int pe_kmax_ = splitting.get_neighbour_processor(1, 2);

  // send left layer of real cells to right layer of virtual cells
  exchange_data(pe_imin_, 0,   0, 0, 
		pe_imax_, ni_, 0, 0, 
		le_ghost, nj_, nk_); /* size of block data to send */

  // send right real cells to left virtual cells
  exchange_data(pe_imax_, ni_ - le_ghost, 0, 0,
					 pe_imin_,     - le_ghost, 0, 0,
					 le_ghost, nj_, nk_);

  exchange_data(pe_jmin_, - le_ghost, 0, 0,
					 pe_jmax_, - le_ghost, nj_, 0,
					 ni_ + 2*le_ghost, le_ghost, nk_);

  exchange_data(pe_jmax_, - le_ghost, nj_ - le_ghost, 0,
					 pe_jmin_, - le_ghost,     - le_ghost, 0,
					 ni_ + 2*le_ghost, le_ghost, nk_);

  exchange_data(pe_kmin_, - le_ghost, - le_ghost, 0,
					 pe_kmax_, - le_ghost, - le_ghost, nk_,
					 ni_ + 2*le_ghost, nj_ + 2*le_ghost, le_ghost);

  exchange_data(pe_kmax_, - le_ghost, - le_ghost, nk_ - le_ghost,
					 pe_kmin_, - le_ghost, - le_ghost,     - le_ghost,
					 ni_ + 2*le_ghost, nj_ + 2*le_ghost, le_ghost);
  statistiques().end_count(echange_vect_counter_);
}

// Initializes the field and allocates memory
// splitting: reference to the geometry of the IJK mesh and how the mesh is split on processors.
//   The field stores a reference to this IJK_Splitting object so do not delete it.
// loc: localisation of the field (elements, nodes, faces in direction i, j, or k)
//   The number of "real" items in each direction (returned by the ni(), nj() or nk() method) is obtained from
//   the IJK_Splitting object. Warning: on a processor that is in the middle of the mesh, the nodes on the
//   right of the rightmost real element are not real, they are virtual, values are copied from the neigbour
//   processor.
// ghost_size: number of ghost layers to allocate. When an exchange of ghost cells data is requested, a smaller
//   number of layers can be requested if all layers are not needed by the following operations
// additional_k_layers: allocates more layers of cells in the k direction for use with the
//   shift_k_origin() method (optimization trick used by the temporally blocked algorithms in the multigrid
//   solver
// nb_compo: number of components of the field. Warning: you cannot allocate the velocity field at faces
//   with nb_compo=3: numbers of faces in each direction differ.
//   Also, components are not grouped by node but stored by layers in k. nb_compo>1 is essentially used
//   in the multigrid solver to optimize memory accesses to the components of the matrix.
void IJK_Field_double::allocate(const IJK_Splitting & splitting, IJK_Splitting::Localisation loc, 
				int ghost_size, int additional_k_layers, int ncompo)
{
  const int ni_local = splitting.get_nb_items_local(loc, 0);
  const int nj_local = splitting.get_nb_items_local(loc, 1);
  const int nk_local = splitting.get_nb_items_local(loc, 2);
  IJK_Field_local_double::allocate(ni_local, nj_local, nk_local, ghost_size, additional_k_layers, ncompo);
  splitting_ref_ = splitting;
  localisation_ = loc;
}

double norme_ijk(const IJK_Field_double & residu)
{
  const int ni = residu.ni();
  const int nj = residu.nj();
  const int nk = residu.nk();
  double somme = 0.;
  for (int k = 0; k < nk; k++) {
    double partial1 = 0.;
    for (int j = 0; j < nj; j++) {
      double partial2 = 0.;
      for (int i = 0; i < ni; i++) {
	double x = residu(i,j,k);
	partial2 += x*x;
      }
      partial1 += partial2;
    }
    somme += partial1;
  }
  somme = Process::mp_sum(somme);
  return sqrt(somme);
}

double max_ijk(const IJK_Field_double& residu)
{
  const int ni = residu.ni();
  const int nj = residu.nj();
  const int nk = residu.nk();
  double res = -1.e30;
  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        {
          double x = residu(i,j,k);
	res = std::max(x,res);
        }
  res = Process::mp_max(res);
  return res;
}

double prod_scal_ijk(const IJK_Field_double & x, const IJK_Field_double & y)
{
  const int ni = x.ni();
  const int nj = x.nj();
  const int nk = x.nk();
  double somme = 0.;
  const ArrOfDouble & xd = x.data();
  const ArrOfDouble & yd = y.data();
  for (int k = 0; k < nk; k++) {
    for (int j = 0; j < nj; j++) {
      for (int i = 0; i < ni; i++) {
	int index = x.linear_index(i,j,k);
	assert(index == y.linear_index(i,j,k));
	somme += xd[index] * yd[index];
      }
    }
  }
  somme = Process::mp_sum(somme);
  return somme;
}

double somme_ijk(const IJK_Field_double & residu)
{
  const int ni = residu.ni();
  const int nj = residu.nj();
  const int nk = residu.nk();
  double somme = 0.;
  for (int k = 0; k < nk; k++) {
    double partial1 = 0.;
    for (int j = 0; j < nj; j++) {
      double partial2 = 0.;
      for (int i = 0; i < ni; i++) {
	double x = residu(i,j,k);
	partial2 += x;
      }
      partial1 += partial2;
    }
    somme += partial1;
  }
  somme = Process::mp_sum(somme);
  return somme;
}


