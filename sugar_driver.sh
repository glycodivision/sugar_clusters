#!/bin/sh
#/usr/bin/env vmd -dispdev text -e 
#/usr/bin/env vmd -e `dirname $0`/sugar_1_1.tcl -args $@
/usr/bin/env vmd -dispdev text -e `dirname $0`/sugar_1.5.tcl -args $@

if [ "$1" = "-h" ] ; then
echo " "
echo " SUGARCLUSTERS es un script para buscar sitios de hidratacion (Watersites)"
echo " "
echo " Utiliza VMD en modo texto."
echo " "
echo " - Requiere que los archivos de topologia, trayectoria y referencia esten ubicados en el mismo directorio donde se ejecuta."
echo " - Lee los parametros especificados en el archivo 'parameters.in'"
echo " - La ultima linea del archivo parameters.in, es un texto plano con el mask de residuos a clusterizar en nomenclatura VMD"
echo "   (Ej: resid 1 47 78 80 121 122 123 125)"
echo " "
echo " Comando: sugarclustdebug -a parameters.in"
echo " "
echo " "
echo "   trayectoria.nc    # TRAJECTORY FILE( *.nc or *.binpos)" >> parameters.in
echo "   topologia.parm7   # TOPOLOGY FILE (*.parm7 or *.prmtop)" >> parameters.in
echo "   referencia.pdb    # REFERENCE FILE (*.pdb)" >> parameters.in
echo "   100               # PERCENT (% of trajectory frames to be used)" >> parameters.in
echo "   0.6   		   # CLUSTER RADIUS (float)" >> parameters.in
echo "   sitio_union.dat"  >> parameters.in
exit
fi


