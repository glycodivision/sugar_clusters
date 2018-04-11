package provide clustering 1.0

namespace eval ::clustering {
    namespace export asignar_corte ldelete clusterizar
}


proc ::clustering::asignar_corte { atom } {

	if { [string first "O" $atom] > -1 } { 
		return 1.4 
	} elseif { [string first "H" $atom] > -1 } {
		return 1
	} elseif { [string first "C" $atom] > -1 } {
		return 1.96
	} else {
		puts "CORTE : atomo $atom no encontrado, se asigna 0.6 como radio de corte"
		return 0.6
	}

}

# borra el elemento $indice de la lista
proc ::clustering::ldelete { lista indice } {
	return [ concat [ lrange $lista 0 [expr { $indice - 1  } ] ] [ lrange $lista [ expr { $indice + 1  } ] end  ] ] 
}

proc ::clustering::clusterizar { distcut ncut pdb_overlap } {
	
	
	
	set overlap [mol new $pdb_overlap filebonds off autobonds off waitfor all]
	
	

	puts "cargue archivo $pdb_overlap"

	set todo [ atomselect $overlap "all" ]
	
	set lista_precluster {}	
	
	############ AVANCE #######################################################
	puts "Cantidad mínima de elementos por precluster: [format %0.0f $ncut] "

	puts "Configurando preclusters\n0  10  20  30  40  50  60  70  80  90  100\n|---|---|---|---|---|---|---|---|---|---|"
	
	set flag_total [$todo num ]
	set flag_marca [expr { $flag_total / 41 }]
	set flag_avance 1 
	###########################################################################

	foreach atomo [ $todo list ] {

		set punto [ atomselect $overlap "within $distcut of index $atomo" ]
		
		set cantidad_puntos [ $punto num ]

		if { $cantidad_puntos > $ncut } {
			lappend lista_precluster [ list $atomo $cantidad_puntos ]
		}
		
		############ AVANCE ##############	
		incr flag_avance
		if { [expr { $flag_avance % $flag_marca  } == 0] } {
			puts -nonewline "#"
			flush stdout
		}
		##################################

		$punto delete 

	}
		

	set lista_precluster [lsort -decreasing -integer -index 1 $lista_precluster ]
	
	set lista_clusters {}

	set i 0
	while { $i < [ llength $lista_precluster ] } {

		puts "$iº  cluster incorporado "
		
		lappend lista_clusters [ lindex $lista_precluster $i 0 ]
	    
	    set punto [ atomselect $overlap " within [expr { $distcut * 2 } ] of index [lindex $lista_precluster $i 0 ] "]
		
		set puntos_en_cluster [ $punto list ] 

		set j [expr { $i + 1 }] 
		
		while { $j  < [ llength $lista_precluster ] } {
			

			set indice_atom [ lindex $lista_precluster $j 0 ]
			
			if { [lsearch -sorted -increasing -integer -exact $puntos_en_cluster $indice_atom ] > -1 } {
					
				set lista_precluster [ ldelete $lista_precluster $j ]
			    
				continue 
			}
			incr j
		
		}
				 
		$punto delete
		incr i		
	}

	#####AVANCE#####
	#puts "\n\n"
	################

	
	mol delete $overlap

	return $lista_clusters
}
