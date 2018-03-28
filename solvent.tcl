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