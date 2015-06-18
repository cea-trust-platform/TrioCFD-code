# FTD_2D_Axi
echo Oscilllation d\'une goutte rho_l=1000 rho_v=100 mu_l=0.002 mu_v=0.0002 sigma=0.05 r=7.5mm
echo Periode theorique d\`apres eq6 de \"Test-case number 5: Oscillation of an inclusion immersed in a quiescent  uid \(PA\)\"
echo 'omega^2=24*sigma/((3*rho_l+2*rho_v)*r^3)'
echo t = 2*pi/omega = 0.209 s
cat <<END | gnuplot
set output "compare_sonde.eps"
set terminal postscript color eps
plot "FTD_2D_Axi_VITESSE_reference.son" title "reference" w l,\
"FTD_2D_Axi_VITESSE.son" title "version_actuelle" w l,\
0.01*sin(29.8*x) title "periode theorique pour mu=0" w l
END
echo Pour voir le fichier resultat: gv compare_sonde.eps
gv compare_sonde.eps
