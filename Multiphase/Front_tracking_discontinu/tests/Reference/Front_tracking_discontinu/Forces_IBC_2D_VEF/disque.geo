Include "donnees.geo" ;
Include "disque_interne.geo" ;
// Surface du cercle
Plane Surface(1) = {1};
// Definition du maillage lineique (reunion des lignes)
Disque=1000;
Physical Line(Disque)={1,2,3,4};
