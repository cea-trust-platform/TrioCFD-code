//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Assemblage.h>
#include <Coeur.h>
#include <Modif_Bord_Assemblage.h>
void instancie_src_UtilitairesAssemblages() {
Cerr << "src_UtilitairesAssemblages" << finl;
Assemblage inst1;verifie_pere(inst1);
Coeur inst2;verifie_pere(inst2);
Modif_Bord_Assemblage inst3;verifie_pere(inst3);
}
