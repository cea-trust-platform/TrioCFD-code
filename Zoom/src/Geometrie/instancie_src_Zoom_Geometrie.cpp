//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Connectivites_DWF.h>
#include <Connectivites_IndGros.h>
#include <Connectivites_faces_couple.h>
#include <ConstruireDomaine.h>
#include <PaveCoincidant.h>
void instancie_src_Zoom_Geometrie() {
Cerr << "src_Zoom_Geometrie" << finl;
Connectivites_DWF inst1;verifie_pere(inst1);
Connectivites_IndGros inst2;verifie_pere(inst2);
Connectivites_faces_couple inst3;verifie_pere(inst3);
ConstruireDomaine inst4;verifie_pere(inst4);
PaveCoincidant inst5;verifie_pere(inst5);
}
