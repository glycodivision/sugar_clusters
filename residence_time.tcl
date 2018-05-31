        # toma como argumentos las vector coordenadas del sitio y el molid sobre 
        # el que tiene que buscar los solventes. y el tipo de solvente y el nombre del grupo.
        # retorna una lista de frames que contienen al menos un solvente en un dado volumen.
proc solvent_finder {coordenadas_sitio molid solvent name R90} {
        #  seteo el numero de frames totales.
        set num_steps [molinfo $molid get numframes]
        #defino las coordenadas x y z  de cada sitio:
        set x [lindex $coordenadas_sitio 0]
        set y [lindex $coordenadas_sitio 1]
        set z [lindex $coordenadas_sitio 2]
        # defino una lista para ser usada luego.
        set frame_list {}
        set serial_sel {}
        for {set f 0} {$f<$num_steps} {incr f} {
                #       defino un volumen de R90A de radio sobre la que busco el atomo de interes de la moleculas de solvente
                set sel [[atomselect $molid "(sqrt( sqr(x-$x) + sqr(y-$y) + sqr(z-$z) ) < $R90) and (resname $solvent and name $name)" frame $f] get serial]
                #       si encuentro en un frame en particular una molecula de solvente en ese volumen entonces agrego el frame en una lista. 
                if {[llength $sel] != 0} {
                        append frame_list "$f "
                        append serial_sel "$sel "
                }
        }
        # retorno un vector que tiene tiene los frames para los cuales tengo un site.
        lappend outlist $frame_list $serial_sel
	package require multiplot
	# This plot will be immediately created because we specified -plot
	set plothandle [multiplot -x $frame_list -y $serial_sel -title "Residence time plot" -ymin 0 -ymax auto -lines -linewidth 3 -marker point -plot]
        return $outlist
}


set c  {  37.298     30.449     52.684}
# Modo de uso.
set vector [solvent_finder $c 0 WAT O 1.00]
# La primer lista son los frames que contienen un solvente en el volumen sugerido.
set frame_list [lindex $vector 0]
# La segunda lista retornada son los serials del solvente que encontro.
set serial_list [lindex $vector 1]

