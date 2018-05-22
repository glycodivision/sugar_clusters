
package provide solvent 

namespace eval ::solvent {
    namespace export WPF calcular_parametros_SS escribir_pdb r90
}

proc ::solvent::r90 { punto radio_max cluster } {
	
 	#set id_overlap [mol new $pdb_overlap filebonds off autobonds off waitfor all]

	#set centro [atomselect $id_overlap "index $centro_index" ]

	set gorfi [measure gofr $punto $cluster delta .005 rmax $radio_max first 0 last 0]

	set posicion [lindex $gorfi 0]

	set integral_ac [lindex $gorfi 2]

	set max_ac [lindex $integral_ac end]

	set i [expr { [llength $posicion] - 1 } ]
	
	set flag_max False
	
	while { $i > 0  && !$flag_max } {

		puts "$i [ lindex $integral_ac $i ]"
		set proporcion [expr { [ lindex $integral_ac $i ] / ( $max_ac * 1.0) } ]


		if { !$flag_max && $proporcion <  0.9} {

			set r90_min [lindex $posicion $i ]
			set r90_max [lindex $posicion [expr {$i + 1}] ]
			set flag_max True
		}
		set i [expr { $i - 1 } ]
	}	

	

	return [expr { ( $r90_min + $r90_max ) / 2.0 } ]

}
# 0.00224 es la densidad
# rescatar "r 60"
# en 1 A cubico hay 0.030219679 aguas en cualquier momento dado
# para una MD de 500 fotos por ns, en 20ns deberian haber 302 aguas en el bulk
proc ::solvent::WFP { n_w_cluster N_fotos_total WFRr } { ; return [expr { $n_w_cluster / ( $N_fotos_total * (( $WFRr ) ** 3) * 4/3 * 3.1416 * 0.0334 )}]  }

proc ::solvent::number_element_cluster_06 { punto radio_max cluster_atoms } {
	
 	#set id_overlap [mol new $pdb_overlap filebonds off autobonds off waitfor all]

	#set centro [atomselect $id_overlap "index $centro_index" ]

	set gorfi [measure gofr $punto $cluster_atoms delta .005 rmax $radio_max first 0 last 0]

	set posicion [lindex $gorfi 0]

	set flag_max False
	
	set integral_ac [lindex $gorfi 2]
	
	set j 0
	foreach i $posicion {
		if { $i == 0.6025} {

			return [ lindex $integral_ac $j ]
		}
		incr j	
	}
}


# entran los indices y devuelve una lista ws = {  { nWS cantidad_aguas WFP R index }  }
proc ::solvent::calcular_parametros_SS { indices  num_frames WFRr pdb_overlap mol_probe atom radio } {
	
	puts "calculando parametros Solv Sites"

	puts "archivo : $pdb_overlap"
	set id_overlap [mol new $pdb_overlap filebonds off autobonds off waitfor all]
	puts "mol mem $id_overlap"	


	set f [open "MolSites/solvent_site_$mol_probe.$atom.csv" w ]
	
	puts $f "SS;x;y;z;#atomos;SPF;R90;indice_overlap"
	
	set lista_SS {}

	set i 1
	foreach indice_atom $indices {
		
		set cluster_atoms 	[atomselect $id_overlap "within $radio of index $indice_atom" ]		
		
		set punto 			[atomselect $id_overlap "index $indice_atom"]

		set x_punto 		[$punto get x]
		set y_punto 		[$punto get y]
		set z_punto 		[$punto get z] 

		set tipo		[$punto get name]
                set number_cluster_elemt	[llength [ $cluster_atoms get index]] 
		set n_s_cluster 	[number_element_cluster_06 $punto $radio $cluster_atoms]
		puts "calcule n_s_cluster da = $n_s_cluster"

		puts "  $n_s_cluster  $num_frames $WFRr  "
		set wfp 			[WFP $n_s_cluster $num_frames $WFRr ]
		#punto radio cluster
		set R90 			[r90 $punto $radio $cluster_atoms]
		$cluster_atoms delete
		
		lappend lista_SS    [ list  "SS_$i" 				  \
							[format "%.3f" $x_punto ]  \
							[format "%.3f" $y_punto ]  \
							[format "%.3f" $z_punto ]  \
							[format "%.0f" $number_cluster_elemt ]  \
							[format "%.2f" $wfp ]  	  \
							[format "%.2f" $R90 ]  	  \
							$indice_atom]

		lappend lista_SS_csv [join [ list  "SS_$i" 				  \
								[format "%.3f" $x_punto ]  \
								[format "%.3f" $y_punto ]  \
								[format "%.3f" $z_punto ]  \
								[format "%.0f" $number_cluster_elemt ]  \
								[format "%.2f" $wfp ]  	  \
								[format "%.2f" $R90 ]  	  \
								$indice_atom]  ";"]
								

		$punto delete
		incr i
	}

	puts $f [join $lista_SS_csv "\n"]
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
