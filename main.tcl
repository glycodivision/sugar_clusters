#set dir [lindex $argv 2]
set dir "/home/glyco/Dropbox/git-projects/sugar_clusters/sugar_clusters"

foreach src [list solvent.tcl qt_clustering.tcl overlap.tcl parser.tcl residence_time.tcl ] {
	source "$dir/$src"
}

############################# ARGUMENTOS ################################

set dict_parametros [parser::parsear_parametros "parameters.in"]

set foldname		[dict get $dict_parametros foldname]
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
puts "foldname	  : $foldname"
puts "binding_site: $binding_site"
puts "trayectoria : $trayectoria"
puts "referencia  : $referencia"
puts "salto  	  : $salto"
puts "radio       : $radio"
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

file delete -force "overlaps"
file delete -force "molsites"

file mkdir "overlaps"
file mkdir "molsites"

################### ARMADO DE OVERLAPING ################################

if { $aguas == "TRUE" 	} {

	set lista_pruebas [list $cosolvent_mol_prueba $cosolvent_atoms_probe "WAT" {O} ]
	
} else {

	set lista_pruebas [list $cosolvent_mol_prueba $cosolvent_atoms_probe ]
}
puts "llegue a solapamiento"

overlap::solapamiento_dinamica $id_dinamica $id_referencia 1 $binding_site $lista_pruebas

mol delete $id_dinamica
mol delete $id_referencia

################### CLUSTERING Y CALCULO DE WS ##########################

set WFRr 0.62035


foreach mol_probe [dict key $lista_pruebas] {
	foreach atom [dict get $lista_pruebas $mol_probe] {
		puts "mol : $mol_probe, atom: $atom "
		puts "overlap pdb : overlaps/overlap_$mol_probe.$atom.pdb"
		set overlap_pdb "overlaps/overlap_$mol_probe.$atom.pdb"
	    
	    set radio [ ::qt_clustering::asignar_corte $atom]	
		set lista_indices_cluster [qt_clustering::clusterizar $radio $ncut $overlap_pdb]
				# $centro_index $radio $pdb_overlap 
						  # solvent::calcular_parametros_SS indices                 num_frames WFRr pdb_overlap mol_probe atom radio
		set solvents_sites [solvent::calcular_parametros_SS $lista_indices_cluster $num_frames $WFRr $overlap_pdb $mol_probe $atom $radio]

		puts "$solvents_sites"

		
		if {[llength $solvents_sites] > 0} {
			
			::solvent::escribir_pdb $solvents_sites $overlap_pdb $mol_probe $atom
			puts "CLUSTERING '$mol_probe $atom' terminado. "
	
		}
	}
}


################### TIEMPO DE RESIDENCIA #################################


# archivos en dir con los csvs
set solvent_sites_files [glob molsites/*.csv]

############################################################################################## ELIJO STEP !!!
set step 10

# abro referencia
set id_referencia	[mol new $referencia autobonds off ]

# abro topologia y agrego dinamica
set id_dinamica [mol new $topologia filebonds off autobonds off]
mol addfile $trayectoria step $salto waitfor all molid $id_dinamica

set num_frames [molinfo $id_dinamica get numframes]

# 
::residence_time::tiempo_residencia $solvent_sites_files $id_dinamica $binding_site


file delete -force $foldname
file mkdir $foldname
exec mv overlaps molsites $foldname 

quit	
