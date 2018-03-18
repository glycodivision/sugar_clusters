#######################	DEFINICION DE PARSERS ############################

proc parsear_parametros { path_archivo } {

	# funcion para parsear el archivo de parametros
	set dic_argumentos [dict create]

	set f [open $path_archivo]

	while { [gets $f line] >= 0 } {
    	puts "$line"
    	set dic_argumentos [parser $line $dic_argumentos]
    	puts dic_argumentos
	}

	close $f
	return $dic_argumentos 
}

proc string_op { palabra } {



	# separo por '=', queda todo el lado derecho
	set der [lindex [split $palabra "=" ] 1 ]

	# separo linea de comentario '#' y quedo con argumento
	set d_arg [lindex [split $der "#" ] 0 ]

	# saco espacios en blanco
	set t_arg [string trim $d_arg]

	return $t_arg
}


proc parser { linea dic_argumentos} {

	# toma una linea del archivo de parametros y parsea la linea, devuelve una actualizacion de diccionario
	puts [string trim [lindex [split $linea "=" ] 0 ] ]
	

	if { [string match "TRAJECTORY" [string trim [lindex [split $linea "=" ] 0 ] ] ] == 1 } {

		set argumento [ string_op $linea ]
		dict set dic_argumentos "trayectoria" $argumento

	}	elseif { [string match "TOPOLOGY" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {

			set argumento [ string_op $linea ]
			dict set dic_argumentos "topologia" $argumento
		

	}	elseif { [string match "REFERENCE" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {
		
			set argumento [ string_op $linea ]
			dict set dic_argumentos "referencia" $argumento

	}	elseif { [string match "STEP" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {

			set argumento [ string_op $linea ]
			dict set dic_argumentos "step" $argumento

	}	elseif { [string match "CLUSTER_RADIUS" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {
			set argumento [ string_op $linea ]
			dict set dic_argumentos "radio" $argumento

	}	elseif { [string match "BINDING_SITE" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {
			
			set argumento [ string_op $linea ]
			dict set dic_argumentos "binding_site" $argumento

	}	elseif { [string match "SOLVENT" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {
			
			set argumento [ string_op $linea ]
			dict set dic_argumentos "solvent_probe" $argumento

	}	elseif { [string match "COSOLVENT" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {
		
			set argumento [ string_op $linea ]
			dict set dic_argumentos "cosolvent_mol_probe" $argumento

	}	elseif { [string match "CS_ATOM_PROBE" [string trim [lindex [split $linea "=" ] 0 ] ] ] == 1 } {

			set argumento [ string_op $linea ]

			
			set argumento [string map {" " ""} $argumento]

			set argumento [split $argumento ","]

			dict set dic_argumentos "cosolvent_atoms_probe" $argumento
		

	}	elseif { [string match "SOLVENT_THRESHOLD" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {

			set argumento [ string_op $linea ]
			dict set dic_argumentos "solvent_threshold" $argumento

		

	}	elseif { [string match "COSOLVENT_THRESHOLD" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {

			set argumento [ string_op $linea ]
			dict set dic_argumentos "cosolvent_threshold" $argumento

		

	}	elseif { [string match "N_CUT_RATIO" [string trim [lindex [split $linea "=" ] 0 ] ] ]  == 1 } {

			set argumento [ string_op $linea ]
			dict set dic_argumentos "n_cut_ratio" $argumento

		

	}

	return $dic_argumentos
}




############################ OVERLAPING #################################

proc solapamiento_dinamica { dinamica referencia salto binding_site mols_probe } {

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

proc procesar_temporal { temp_pdb } {	
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

############################ FIN OVERLAPING #############################

############################# CLUSTERING ################################


# borra el elemento $indice de la lista
proc ldelete { lista indice } {
	return [ concat [ lrange $lista 0 [expr { $indice - 1  } ] ] [ lrange $lista [ expr { $indice + 1  } ] end  ] ] 
}

proc clusterizar { distcut ncut pdb_overlap } {
	
	
	
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


proc std_salida { archivo } {
	flush stdout
	close stdout
   	open $archivo w
   	return 1
}

#################################### FIN CLUSTERING ####################################

#################################### CALCULO DE WS #####################################
proc WFP { n_w_cluster N_fotos_total WFRr } { ; return [expr { $n_w_cluster / ( $N_fotos_total * (( $WFRr ) ** 3) * 4/3 * 3.1416 * 0.0334 )}]  }

# entran los indices y devuelve una lista ws = {  { nWS cantidad_aguas WFP R index }  }
proc calcular_parametros_SS { indices distcut num_frames R90 WFRr pdb_overlap mol_probe atom } {
	
	puts "calculando parametros Solv Sites"

	puts "archivo : $pdb_overlap"
	set id_overlap [mol new $pdb_overlap filebonds off autobonds off waitfor all]
	puts "mol mem $id_overlap"	


	set f [open "MolSites/solvent_site_$mol_probe.$atom.csv" w ]
	
	puts $f "SS;x;y;z;#atomos;SPF;R90;indice_overlap"
	
	set lista_SS {}

	set i 0
	foreach indice_atom $indices {
		
		set cluster_atoms 	[atomselect $id_overlap "within $distcut of index $indice_atom" ]
		set cantidad_aguas  [$cluster_atoms num]
		$cluster_atoms delete
		
		set punto 			[atomselect $id_overlap "index $indice_atom"]

		set x_punto 		[$punto get x]
		set y_punto 		[$punto get y]
		set z_punto 		[$punto get z] 

		set tipo			[$punto get name]

		set wfp 			[WFP  $cantidad_aguas $num_frames $WFRr ]
		
		
		lappend lista_SS [list  "SS_$i"					\
								$x_punto				\
								$y_punto				\
								$z_punto				\
								$cantidad_aguas			\
								[format "%.2f" $wfp ] 	\
								$R90 					\
								$indice_atom ]
		
		puts $f [join $lista_SS ";"]

		$punto delete
		incr i
	}
	
	mol delete $id_overlap
	close $f
	return $lista_SS
}

############################# FIN CALCULO SS ############################


############################# PDBS SOLVENT S ############################

proc escribir_pdb { lista_SS overlap_pdb mol_probe atom } {
	
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

############################# ARGUMENTOS ################################

set dict_parametros [parsear_parametros "parameters.in"]

set binding_site	[dict get $dict_parametros binding_site]

set trayectoria		[dict get $dict_parametros trayectoria ]
set topologia		[dict get $dict_parametros topologia ]
set referencia		[dict get $dict_parametros referencia ]

set n_cut_ratio		[dict get $dict_parametros n_cut_ratio ]
set radio			[dict get $dict_parametros radio ]

set binding_site			[dict get $dict_parametros binding_site ]
set cosolvent_mol_prueba	[dict get $dict_parametros cosolvent_mol_probe ]
set cosolvent_atoms_probe	[dict get $dict_parametros cosolvent_atoms_probe ]	
set salto 					[dict get $dict_parametros step ]

set aguas 			[dict get $dict_parametros solvent_probe]

puts "\n\n\n"
puts "binding_site: $binding_site"
puts "trayectoria : $trayectoria"
puts "referencia  : $referencia"
puts "salto  	  : $salto"
puts "radio       : $radio "
puts "binding_site: $binding_site"	
puts "mol_prueba  : $cosolvent_mol_prueba "
puts "atomos prueb: $cosolvent_atoms_probe"
puts "n_cut_ratio : $n_cut_ratio"
puts "aguas		  : $aguas"
puts "\n\n\n"
################### APERTURA DINAMICA Y REFERENCIA #####################

set id_referencia	[mol new $referencia autobonds off ]


set id_dinamica		[mol new $topologia filebonds off autobonds off]
mol addfile $trayectoria step $salto waitfor all molid $id_dinamica

set num_frames [molinfo $id_dinamica get numframes]
set ncut [expr { $num_frames * $n_cut_ratio}]

file mkdir "MolSites"
file mkdir "ws"

################### ARMADO DE OVERLAPING ################################

if { $aguas == "TRUE"} {

	set lista_pruebas [list $cosolvent_mol_prueba $cosolvent_atoms_probe "WAT" {O} ]
	
} else {

	set lista_pruebas [list $cosolvent_mol_prueba $cosolvent_atoms_probe ]
}
puts "llegue a solapamiento"
solapamiento_dinamica $id_dinamica $id_referencia 1 $binding_site $lista_pruebas

mol delete $id_dinamica
mol delete $id_referencia

################### CLUSTERING Y CALCULO DE WS ##########################

set WFRr 0.6
set R90 1

foreach mol_probe [dict key $lista_pruebas] {
	foreach atom [dict get $lista_pruebas $mol_probe] {
		puts "mol : $mol_probe, atom: $atom "
		puts "overlap pdb : ws/overlap_$mol_probe.$atom.pdb"
		set overlap_pdb "ws/overlap_$mol_probe.$atom.pdb"
		
		set lista_indices_cluster [clusterizar $radio $ncut $overlap_pdb]

		#calcular_parametros_WS { 					indices 			distcut num_frames R90  WFRr id_overlap }

		set solvents_sites [calcular_parametros_SS $lista_indices_cluster $radio $num_frames $R90 $WFRr $overlap_pdb $mol_probe $atom]

		puts "$solvents_sites"

		
		if {[llength $solvents_sites] > 0} {
			
			escribir_pdb $solvents_sites $overlap_pdb $mol_probe $atom
			puts "CLUSTERING '$mol_probe $atom' terminado. "
	
		}
	}
}	
