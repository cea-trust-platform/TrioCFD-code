#!/bin/bash
cas="cst_growth"
err=$cas".err"

trust $cas 1>$cas.out 2> $err
[ $? != 0 ] && echo "Abnormal termination of calculation!" && exit

\rm -f sum_dIdt_or_dV.txt ai.txt vol_vap.txt sum_dIdt_after-PCH.txt vi.txt
grep "^Volume_phase_0" $err | awk '{print $4, $2}' > vol_vap.txt
grep "^Surface_Totale_Interface" $err | awk '{print $4, $6}' > ai.txt
grep "sum_dI_dt" $err | awk '{print $2, $4, $6, $8, $10, $12}' > sum_dIdt_or_dV.txt
grep "^Interfacial_velocity" $err | awk '{print $6, $8, $9, $10, $12, $13, $14}' > vi.txt

grep "R-PCH].* The sum" $err | awk '{print $2, $7}' >  sum_dIdt_after-PCH.txt
python plot.py && display plot.png

echo """
Many new options in development : 
If ns.newi_mass_source, secmem2 is computed from ai*mp in the cell instead of div(delta_u0). 
                (and also ajouter_contribution_saut_vitesse(deplacement) is modified to use a 1D-interp 
                 per component instead of the standard MultiD).
Local extension of gas velocity : vitesse_0_
Available types of methods for calculer_dI_dt { INTERP_STANDARD, INTERP_MODIFIEE, INTERP_AI_BASED }
                (Changes are NOT effective for (phase_pilote != 1))
Available types of methods in calculer_indicatrices_faces  { STANDARD, MODIFIEE, AI_BASED }
Additional keyword to prescribe an ad-hoc phase-change rate from the datafile (prescribed_mpoint). 
Correction of post-pro vitesse interfaciale : ajouter_contribution_saut_vitesse(vit) was missing.

dI_dt already has volume in it. calculer_dI_dt has been corrected to add mpoint[elem]*interfacial_area[elem]*un_sur_rho_0 to the result

Temporary : Flag DEBUG_CONSERV_VOLUME to 1.

"""
