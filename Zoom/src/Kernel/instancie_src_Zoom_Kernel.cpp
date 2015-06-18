//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Champ_front_zoom.h>
#include <Connectivites.h>
#include <Connectivites_base.h>
#include <Pb_2G.h>
#include <Pb_MG.h>
#include <Prolongement.h>
#include <Restriction.h>
void instancie_src_Zoom_Kernel() {
Cerr << "src_Zoom_Kernel" << finl;
Champ_front_zoom inst1;verifie_pere(inst1);
Connectivites inst2;verifie_pere(inst2);
Connectivites_base inst3;verifie_pere(inst3);
Pb_2G inst4;verifie_pere(inst4);
Pb_MG inst5;verifie_pere(inst5);
Prolongement inst6;verifie_pere(inst6);
Restriction inst7;verifie_pere(inst7);
}
