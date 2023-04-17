#!/bin/bash

Calcul_Uplus()
{
   jdd=$ESSAI/$SCHEME

   NAMEs="branch main0.23" && [ "$ESSAI" == "Wall_Jet" ] && NAMEs="branch main1.46"

   for NAME in $NAMEs
   do
      mu=`grep "mu" Watlon.data | head -1 | $TRUST_Awk '{print $NF}'`
      rho=`grep "rho" Watlon.data | tail -1 | $TRUST_Awk '{print $NF}'`
      Dh=`grep "diametre_hydraulique_branch" Watlon.data | head -1 | $TRUST_Awk '{print $(NF-1)}'` && [ "$NAME" != "branch" ] && Dh=`grep "diametre_hydraulique_main" Watlon.data | head -1 | $TRUST_Awk '{print $(NF-1)}'`

      cd $jdd/$NAME\_box
      
      # Je cree le profil de vitesse
      extrait_coupe $NAME\_box SONDE_V_BOX
      
      # Correlation de Blasius sur Utau
      [ "$NAME" == "branch" ] && U=1
      [ "$NAME" == "main0.23" ] && U=0.23
      [ "$NAME" == "main1.46" ] && U=1.46
      
      Re=`echo "" | $TRUST_Awk -v U=$U -v D=$Dh -v mu=$mu -v rho=$rho '{print D*U/(mu/rho)}'`
      Cf=`echo "" | $TRUST_Awk -v Re=$Re '{print 0.079*Re**-0.25}'`
      Utau_correl=`echo "" | $TRUST_Awk -v Cf=$Cf -v U=$U '{print (0.5*Cf)^0.5*U}'`

      # Je recupere le u+, y+ et u* calcule par Trio_U
      Utau=`grep "<u\*>" *_pb_Ustar.face | tail -1 | $TRUST_Awk '{print $6}'`
      uplus=`grep "<u\*>" *_pb_Ustar.face | tail -1 | $TRUST_Awk '{print $2}'`
      yplus=`grep "<u\*>" *_pb_Ustar.face | tail -1 | $TRUST_Awk '{print $4}'`
      echo "Boite		Utau correlation	Utau Trio_U" > Utau.txt
      echo "$NAME   	$Utau_correl		$Utau" >> Utau.txt
      
      # Je calcul le y+ = Y.u*/nu et u+ = U/u*
      $TRUST_Awk -v Utau=$Utau -v D=$Dh -v mu=$mu -v rho=$rho '{print ($1+D/2)*Utau/(mu/rho) "\t" sqrt($2*$2+$3*$3+$4*$4)/Utau }' $NAME\_box_SONDE_V_BOX.coupe > Uplus.out
      echo "# Y_plus-trio	U_plus-trio" > Uplus_trio.out
      echo "$yplus	$uplus" >> Uplus_trio.out
      
      # Je supprime les 0 pour l'affichage log
      sed 's?0$?supprimer?g' Uplus.out > Uplus.tmp
      sed '/supprimer/d' Uplus.tmp > Uplus.out

      # Comme mon profil va d'une paroi a l'autre, je reduit mon profil au trajet paroi -> centre du canal
      nbl=`wc -l Uplus.out |awk '{print $1/2}'`
      head -$nbl Uplus.out > Uplus.tmp
      mv -f Uplus.tmp Uplus.out
      rm -f Uplus.tmp

      cd - 1>/dev/null 2>&1
   done
}


ESSAIs="Impinging_Jet Wall_Jet"
SCHEMEs="Muscl Ef_Stab1.0 Ef_Stab0.2"

# Je boucle pour entrer dans tout les cas test Trio_U
for ESSAI in $ESSAIs
do
   for SCHEME in $SCHEMEs
   do
      Calcul_Uplus
   done
done
