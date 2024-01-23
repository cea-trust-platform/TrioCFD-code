#!/bin/bash
grep "Volume_phase_1" FTD_remaillage_vdf_PAR/PAR_FTD_remaillage_vdf.err | awk '{print $4,$2}' > par_vdf_v.txt
grep "Volume_phase_1" FTD_remaillage_vdf_SEQ/FTD_remaillage_vdf.err | awk '{print $4,$2}' > vdf_v.txt
grep "Surface_Totale_Interface" FTD_remaillage_vdf_SEQ/FTD_remaillage_vdf.err | awk '{print $4,$6}' > vdf_s.txt
grep "Surface_Totale_Interface" FTD_remaillage_vdf_PAR/PAR_FTD_remaillage_vdf.err | awk '{print $4,$6}' > par_vdf_s.txt

grep "Volume_phase_1" FTD_remaillage_vef_PAR/PAR_FTD_remaillage_vef.err | awk '{print $4,$2}' > par_vef_v.txt
grep "Volume_phase_1" FTD_remaillage_vef_SEQ/FTD_remaillage_vef.err | awk '{print $4,$2}' > vef_v.txt
grep "Surface_Totale_Interface" FTD_remaillage_vef_SEQ/FTD_remaillage_vef.err | awk '{print $4,$6}' > vef_s.txt
grep "Surface_Totale_Interface" FTD_remaillage_vef_PAR/PAR_FTD_remaillage_vef.err | awk '{print $4,$6}' > par_vef_s.txt
