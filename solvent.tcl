
package provide solvent 

namespace eval ::solvent {
    namespace export WPF calcular_parametros_SS escribir_pdb r90
}

proc ::clustering::r90 { $centro_index $radio $overlap_pdb } {



	set centro [atomselect ]

	set gorfi [measure gofr $centro $cluster delta .005 rmax $radio first 0 last 0]

	set posicion [lindex gorfi 0]

	set integral_ac [lindex gorfi 2]

	set max_ac [lindex integral_ac end]

	set i [llength $posicion]
	
	set flag_max False
	
	while { $i > 0  && !flag_max } {

		set proporcion [expr { [ lindex $integral_ac $i ] / $max_ac } ]

		if { !flag_max && proporcion <  0.9 } {

			set r90_min [ lindex $integral_ac $i ]
			set r90_max [lindex $integral_ac [expr {$i + 1}]]
			set flag_max True
		}		
	}	

	set i [expr { $i - 1 } ]

	return [expr {( r90_min + r90_max ) / 2 } ]

}


proc ::solvent::WFP { n_w_cluster N_fotos_total WFRr } { ; return [expr { $n_w_cluster / ( $N_fotos_total * (( $WFRr ) ** 3) * 4/3 * 3.1416 * 0.0334 )}]  }

# entran los indices y devuelve una lista ws = {  { nWS cantidad_aguas WFP R index }  }
proc ::solvent::calcular_parametros_SS { indices  num_frames R90 WFRr pdb_overlap mol_probe atom } {
	
	puts "calculando parametros Solv Sites"

	puts "archivo : $pdb_overlap"
	set id_overlap [mol new $pdb_overlap filebonds off autobonds off waitfor all]
	puts "mol mem $id_overlap"	


	set f [open "MolSites/solvent_site_$mol_probe.$atom.csv" w ]
	
	puts $f "SS;x;y;z;#atomos;SPF;R90;indice_overlap"
	
	set lista_SS {}

	set i 0
	foreach indice_atom $indices {
		
		set cluster_atoms 	[atomselect $id_overlap "within 0.6 of index $indice_atom" ]
		set cantidad_aguas  [$cluster_atoms num]
		$cluster_atoms delete
		
		set punto 			[atomselect $id_overlap "index $indice_atom"]

		set x_punto 		[$punto get x]
		set y_punto 		[$punto get y]
		set z_punto 		[$punto get z] 

		set tipo			[$punto get name]

		set wfp 			[WFP  $cantidad_aguas $num_frames $WFRr ]
		
		
		lappend lista_SS [list  "SS_$i"					  \
								[format "%.3f" $x_punto ] \
								[format "%.3f" $y_punto ] \
								[format "%.3f" $z_punto ] \
								$cantidad_aguas			  \
								[format "%.2f" $wfp ] 	  \
								$R90 					  \
								$indice_atom ]
		
		

		$punto delete
		incr i
	}
	puts $f [join $lista_SS ";"]
	mol delete $id_overlap
	close $f
	return $lista_SS
}


proc ::solvent::escribir_pdb { lista_SS overlap_pdb mol_probe atom } {
	
	puts "escribiendo pdb"
	
	set mol_overlap 	[mol new $overlap_pdb filebonds off autobonds off waitfor all]
	#set mol_referencia	[mol new $referencia  filebonds off autobonds off waitfor all]
	
	foreach SS $lista_SS {

		set atomo [atomselect $mol_overlap "index [lindex $SS 7 ]"]
		
		#asigno numero de SS al centro del cluster
		$atomo set resid [string range [lindex $SS 0] 3 end]
		
		#asigno el valor wfp al beta
		$atomo set beta  [lindex $SS 5]	

		$atomo delete
		
	}

	set indices [list ]

	foreach SS $lista_SS { 
		puts "$SS"
		lappend indices [lindex $SS 7] } 
	
	puts $indices
	set sel [atomselect $mol_overlap "index $indices" ]
	
	puts "atomo: $atom , mol : $mol_probe "

	$sel writepdb "MolSites/SS_${mol_probe}_$atom.pdb"
	
	$sel delete 

	return 1
}