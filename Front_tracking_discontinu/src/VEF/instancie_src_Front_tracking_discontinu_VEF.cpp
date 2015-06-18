//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <SourceFiltree_FT_disc_VEF_P1NC.h>
#include <Source_Reaction_Particules_VEF.h>
void instancie_src_Front_tracking_discontinu_VEF() {
Cerr << "src_Front_tracking_discontinu_VEF" << finl;
SourceFiltree_FT_disc_VEF_P1NC inst1;verifie_pere(inst1);
Source_Reaction_Particules_VEF inst2;verifie_pere(inst2);
}
