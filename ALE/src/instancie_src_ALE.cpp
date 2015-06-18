//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Champ_front_ALE.h>
#include <Domaine_ALE.h>
#include <Imposer_vit_bords_ALE.h>
#include <Op_Conv_ALE_VDF.h>
#include <Op_Conv_ALE_VEF.h>
void instancie_src_ALE() {
Cerr << "src_ALE" << finl;
Champ_front_ALE inst1;verifie_pere(inst1);
Domaine_ALE inst2;verifie_pere(inst2);
Imposer_vit_bords_ALE inst3;verifie_pere(inst3);
Op_Conv_ALE_VDF inst4;verifie_pere(inst4);
Op_Conv_ALE_VEF inst5;verifie_pere(inst5);
}
