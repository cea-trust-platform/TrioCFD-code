sed 'n; d' standard/composantes_connexes.txt | awk '{print $2}' > temps
sed '1d; n; d' standard/composantes_connexes.txt | awk '{print $9}' > vitesse_standard
sed '1d; n; d' arithmetic/composantes_connexes.txt | awk '{print $9}' > vitesse_aritmetic
sed '1d; n; d' harmonic/composantes_connexes.txt | awk '{print $9}' > vitesse_harmonic
paste -d ' ' temps vitesse_standard vitesse_aritmetic vitesse_harmonic > resultat.data
sed -i '1i # temps standard arithmetic harmonic' resultat.data
rm temps vitesse_standard vitesse_aritmetic vitesse_harmonic
