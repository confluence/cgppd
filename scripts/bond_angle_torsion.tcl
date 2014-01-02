set bb [atomselect top "name CA"]
set bbi [$bb list]
set len [llength $bbi]

puts "\nPseudo-bonds:"

for {set i 0} {$i < [expr $len - 1]} {incr i} {
    puts [measure bond [lrange $bbi $i [expr $i + 1]]]
}

puts "\nPseudo-angles:"

for {set i 0} {$i < [expr $len - 2]} {incr i} {
    puts [measure angle [lrange $bbi $i [expr $i + 2]]]
}

puts "\nPseudo-torsions:"

for {set i 0} {$i < [expr $len - 3]} {incr i} {
    puts [measure dihed [lrange $bbi $i [expr $i + 3]]]
}