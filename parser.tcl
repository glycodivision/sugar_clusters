package provide parser

namespace eval ::parser:: {
    namespace export parsear_parametros string_op parser 
}


proc ::parser::parsear_parametros { path_archivo } {

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

proc ::parser::string_op { palabra } {



	# separo por '=', queda todo el lado derecho
	set der [lindex [split $palabra "=" ] 1 ]

	# separo linea de comentario '#' y quedo con argumento
	set d_arg [lindex [split $der "#" ] 0 ]

	# saco espacios en blanco
	set t_arg [string trim $d_arg]

	return $t_arg
}


proc ::parser::parser { linea dic_argumentos} {

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
