
# alineamiento, tiempo acumulado por cada index, 



set topologia "./../1SL4-A_apo_fen_HMR.prmtop"
set trayectoria "./../MD-1-1SL4-A_apo_fen_HMR-200-ns-WAT,FEN.nc"
set step 10
# en ps
set time_interval 2

proc asignar_corte { atom } {

  if { [string first "O" $atom] > -1 } { 
    return 1.4 
  } elseif { [string first "H" $atom] > -1 } {
    return 1
  } elseif { [string first "C" $atom] > -1 } {
    return 1.96
       } elseif { [string first "DU" $atom] > -1 } {
                return 1.96
  } else {
    puts "CORTE : atomo $atom no encontrado, se asigna 0.6 como radio de corte"
    return 0.6
  }

}


# recibe las coordenadas xyz del solvent site
# arma un atomselect que define un cubo de volumen ocupancia^3
# la idea es:lo que queda dentro del cubo podria estar contenido en un solvent site
# despues de acotar de esta forma, se computa la distancia al centro del solvent site
# si esa distancia es menor o igual que la ocupancia, el atomo se considera dentro del solvent site
proc armar_seleccion { coordenadas_xyz mol atom } {

  set corte_atom [asignar_corte $atom]

  set x_min [expr {[lindex $coordenadas_xyz 0] - $corte_atom } ]
  set x_max [expr {[lindex $coordenadas_xyz 0] + $corte_atom } ]

  set y_min [expr {[lindex $coordenadas_xyz 1] - $corte_atom } ]
  set y_max [expr {[lindex $coordenadas_xyz 1] + $corte_atom } ]

  set z_min [expr {[lindex $coordenadas_xyz 2] - $corte_atom } ]
  set z_max [expr {[lindex $coordenadas_xyz 2] + $corte_atom } ]

  return " resname ${mol} and name ${atom} and x > ${x_min} and x < ${x_max} and y > ${y_min} and y < ${y_max} and z > ${z_min} and z < ${z_max} "
}

proc corroborar_distancia { coordenadas_xyz_atom coordenadas_ss index_atom atom_type} {

  set corte [asignar_corte $atom_type]

  set indices_atom_in_ss {}
  if { [llength coordenadas_xyz_atom] > 0 } {
    set i 0
    foreach coord $coordenadas_xyz_atom {
      
      set distancia [expr { ((([lindex ${coord} 0] - [lindex ${coordenadas_ss} 0]) ** 2 ) + (([lindex ${coord} 1] - [lindex ${coordenadas_ss} 1]) ** 2 ) + (([lindex ${coord} 2] - [lindex ${coordenadas_ss} 2]) ** 2 ) ) ** 1 / 2 } ]
      
      if { $distancia <= $corte} {
        lappend indices_atom_in_ss [lindex $index_atom $i]
      }
      incr i
    }
  }
  
  return $indices_atom_in_ss
}

set archivos [glob *.csv]

set general_sites {}
foreach path_archivo $archivos {

    set nombre_archivo [lindex [split $path_archivo "/"] end]

  set atom [lindex [split $nombre_archivo "."] 1 ]

  set mol [lindex [split [lindex [split $nombre_archivo "."] 0] "_"] end]

  set f [open $path_archivo]

  set linea_0 1
  
  while { [gets $f linea] >= 0 } {
      
      if {$linea_0} {
        set linea_0 0

      } else {

          #puts "$linea"

          set s_linea [split $linea ";"]
          
          set SS [lindex $s_linea 0]

          set xyz [lrange $s_linea 1 3]

          lappend general_sites [list $SS $xyz $mol $atom ]

          
        }
  }
  puts $general_sites
  close $f
  
}


set id_dinamica [mol new $topologia filebonds off autobonds off]
mol addfile $trayectoria step $step waitfor all molid $id_dinamica

set num_steps [molinfo $id_dinamica get numframes]


# se recorre por foto y una vez en la foto
# se observa el atom index por cada SS 
for {set frame 0} {$frame < $num_steps} {incr frame 1 } {

  foreach ss $general_sites {

    set SS   [lindex $ss 0]
    set xyz  [lindex $ss 1]
    set mol  [lindex $ss 2]
    set atom [lindex $ss 3]

    set seleccion [armar_seleccion $xyz $mol $atom]

    set atoms_in_ss [atomselect top $seleccion frame $frame]

    set f [open "${SS}_${mol}_${atom}.site" a]

    #sumo el step al frame para obtener el numero correcto de frame en la dinamica
    #puts $f "[expr {${frame} * ${step}}] ; [expr {${frame}*${step}*${time_interval}}] ; [$atoms_in_ss get index]"
    
    set atom_index_in_ss [corroborar_distancia [$atoms_in_ss get { x y z }] $xyz [$atoms_in_ss get index] $atom ] 
    
    puts $f "[expr {${frame} * ${step}}] ; [expr {${frame}*${step}*${time_interval}}] ; ${atom_index_in_ss}"

    close $f
    #$seleccion delete

  }
}

exit
