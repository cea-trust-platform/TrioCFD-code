// #ifndef DEF_FORCE_PH
// #define DEF_FORCE_PH
//
// #include <iostream>
// #include <string>
// #include <vector>
// #include <Force_sp.h>
// #include <FixedVector.h>
// #include <IJK_Field.h>
// #include <fftw3.h>
// #include <communications.h>
//
// class Force_ph
// {
// 	public:
//
// 	Force_ph();
// 	~Force_ph();
// 	void initialise(int nproc_tot, int ni, int nj, int nk, int nl,int nm,int nn,
// 			double Lx, double Ly, double Lz, double Ox,double Oy,double Oz, int momin, int momax, double kmin, double kmax,
// 		  std::string nom_fichier, const IJK_Splitting& splitting);
// 	// void initialise(int ni, int nj, int nk, int nl,int nm,int nn,
// 	// 		double Lx, double Ly, double Lz, double kmin, double kmax);
// 	void from_spect_to_phys(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
// 	void from_spect_to_phys2(const std:: vector <double >& coeff_force);
// 	void from_spect_to_phys_opti(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
// 	void from_spect_to_phys_opti2(const std:: vector <double >& coeff_force);
// 	void from_spect_to_phys_bis(const std::vector< std::vector< std:: vector <double >>>& coeff_force);
// 	void cheat_function();
// 	// void from_spect_to_phys_fftw(fftw_complex* in);
// 	void write(std::string nom_fichier_sortie, double t);
// 	void write_separate(std::string nom_fichier_sortie, double t);
// 	void write_offset_index_position( const IJK_Splitting& my_splitting);
// 	void compute_energie();
// 	double get_energie();
// 	FixedVector<IJK_Field_double, 3> get_force_attribute();
//
// 	void gbz_gather(FixedVector<IJK_Field_double, 3> force_);
//
// 	private:
// 		int nproc_tot;
//
// 	int ni,nj,nk,n_ijk;
// 	int nl,nm,nn,n_lmn;
// 	double Lx, Ly, Lz;
// 	double Ox, Oy, Oz;
// 	double kmin,kmax;
// 	int momin,momax;
// 	FixedVector<IJK_Field_double, 3> force_;
// 	std::vector<std::vector< std:: vector < double > > > force;
// 	double energie;
//
// 	int nproc_i,nproc_j,nproc_k;
//
//
// };
//
// std::vector< std::vector< std:: vector <double >>> set_dimensions(std::vector< std::vector< std:: vector <double >>> the_vector, int dim_one, int dim_two, int dim_three);
//
// // FixedVector<IJK_Field_double, 3> set_to_zero(FixedVector<IJK_Field_double, 3> vector)
// // {
// // 		for (int dir=0; dir<3; dir++)
// // 			for (int i=0; i<vector[0].ni(); i++)
// // 				for (int j=0; j<vector[1].nj(); j++)
// // 					for (int k=0; k<vector[2].nk(); k++)
// // 						vector[dir](i,j,k) = 0.;
// // 		return vector;
// // }
//
// #endif
