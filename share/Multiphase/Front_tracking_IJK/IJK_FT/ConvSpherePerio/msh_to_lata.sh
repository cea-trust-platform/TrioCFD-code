input=$1
lata=$2
export connex_compo=$3
nodefile=$2.nodes
elemfile=$2.elem
compofile=$2.connex_compo
nnodes=`awk 'flag==1{print $1;exit}$1~"Nodes"&&(!($1~"EndNodes")){flag=1}' $input`
awk '$1~"EndNodes"{flag=0}flag==2{print $2,$3,$4}flag==1{n=$1;flag=2}$1~"Nodes"&&(!($1~"EndNodes")){flag=1}' $input>$nodefile
awk '$1~"EndElements"{flag=0}flag==2&&NF==8{print $6,$7,$8}flag==1{n=0;flag=2}$1~"Elements"&&(!($1~"EndElements")){flag=1}' $input>$elemfile
nelem=`awk '$1~"EndElements"{flag=0}flag==2&&NF==8{n=n+1}flag==1{n=0;flag=2}$1~"Elements"&&(!($1~"EndElements")){flag=1}END{print n}' $input`
seq $nelem | awk 'BEGIN{x=ENVIRON["connex_compo"]}{print x}'>$compofile
echo nbnodes=$nnodes
echo nbelem=$nelem
echo LATA_V2.1 >$lata
echo titi >>$lata
echo Trio_U >>$lata
echo Format ASCII,F_INDEXING,C_ORDERING,F_MARKERS_NO,INT32,REAL32 >>$lata
echo TEMPS 0 >>$lata
echo Geom FTMESH type_elem=TRIANGLE_3D >>$lata
echo Champ SOMMETS $nodefile geometrie=FTMESH size=$nnodes composantes=3 >>$lata
echo Champ ELEMENTS $elemfile geometrie=FTMESH size=$nelem composantes=3 FORMAT=INT32 >>$lata
echo Champ COMPO_CONNEXE $compofile geometrie=FTMESH size=$nelem composantes=1 FORMAT=INT32,NO_INDEXING >>$lata
echo FIN >>$lata
