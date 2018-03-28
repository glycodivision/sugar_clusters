
package provide overlap

namespace eval ::overlap {
    namespace export solapamiento_dinamica procesar_temporal
}


proc ::overlap::solapamiento_dinamica { dinamica referencia salto binding_site mols_probe } {

	set rmsd_out [open "ws/rmsd_fit.csv" w]

	#rmsd_fit mask
	puts $rmsd_out "Frame;RMSD"

	
	set reference	[atomselect $referencia "$binding_site" frame 0]
	set compare 	[atomselect $dinamica "$binding_site" ] 
	set num_steps	[molinfo $dinamica get numframes]
	
	
	#inicializo los archivos 
	foreach mol_probe [dict key $mols_probe] { 			
		foreach atom [dict get $mols_probe $mol_probe] {
			set current_overlap [open "ws/overlap_$mol_probe.$atom.pdb" w]
			close $current_overlap
		}
	}

	for {set frame 0} {$frame < $num_steps} {incr frame $salto } {

		

        # compute the transformation
		set trans_mat [measure fit $compare $reference]
		# do the alignment
		set compare_bs [atomselect $dinamica "all" frame $frame ]
		
		$compare_bs move $trans_mat
		
		set rmsd [measure rmsd $compare $reference]

		puts $rmsd_out "$frame;$rmsd"
		puts "RMSD of $frame is $rmsd"

		
        foreach mol_probe [dict key $mols_probe] { 
					
			foreach atom [dict get $mols_probe $mol_probe] {
       			

				

       			
       			# aca va el reparto a distintos mapas....
       			set compare_probe [atomselect $dinamica "(resname $mol_probe and name $atom) and within 5 of ( $binding_site )" frame $frame ]
        		
        		#std_salida "./stdout"
        		$compare_probe writepdb ws/temp_$mol_probe.$atom.pdb
				#std_salida "/dev/tty"
        		
        		set current_overlap [open "ws/overlap_$mol_probe.$atom.pdb" a]
        		
        		foreach linea [procesar_temporal "ws/temp_$mol_probe.$atom.pdb" ] { 
        			puts $current_overlap $linea 
        		}       			
       			close $current_overlap	
       		}       		
		}

	}

	foreach mol_probe [dict key $mols_probe] {			
		foreach atom [dict get $mols_probe $mol_probe] {
			set current_overlap [open "ws/overlap_$mol_probe.$atom.pdb" a]
			puts $current_overlap "END"
			close $current_overlap
		}
	}

    close $rmsd_out

	return 1
}

proc ::overlap::procesar_temporal { temp_pdb } {	
	set temporal [open $temp_pdb r]
	gets $temporal line
	set salida {}	
	while { [gets $temporal line ] >= 0 } {
		if { $line != "END" } {
			lappend salida $line			
		}
	}
	close $temporal
	return $salida
}
