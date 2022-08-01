#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#include "classes_paires.h"
using namespace std;

// a very simple c++ file to show how to read the snapshot



int main () {
	
  double L(0),db(0),repulsion_range(0);
  int Nb(0), Nt(0);
  int id_b_min(0),id_b_max(0);
  int t_min(0),t_max(0);
  int Ncat;
  string chemin_out;

  tableau_positions bulles;
  tableau_distances distances;
  statistiques dist_pdf;

  bulles.initialise(Nb,id_b_min,id_b_max,t_min,t_max,Nt,L,db,repulsion_range);
  bulles.get_positions(chemin_out);
  
  distances.initialise(&bulles);
  distances.compute_cartesian_distance();
  distance.compute_cyl_distance();
  distance.compensate_cyl_distance();
  distance.compute_sphe_distance();
  distance.compensate_sphe_distance();
  
  dist_pdf.initialise(&distance,Ncat);
  dist_pdf.compute_pdf();
  dist_pdf.compute_stats();
  dist_pdf.write();
    
  cout << "Done" << endl;
  return 0;
}
