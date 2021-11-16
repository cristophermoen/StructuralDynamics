 #########################################################################################################################
 # Johns Hopkins University
 # Thin-walled Structures Group
 # Model of CFS Shear Walls
 # developed by Mohammed Eladly and Ben Schafer									      															 #
 #########################################################################################################################

wipe
wipe all
model BasicBuilder -ndm 1 -ndf 1
set EQ_ID EQ45
set dataDir_EQ EQ45
set model_id Model-01
set dt 0.02
set Analysis_Duration 37
set story_Height 2.4384
set bay_Width 0.6096
set mass 3648.58361367032
set weight 35780.382495
set Tn 0.5
set dampRatio 0.02
set Fyield 389.256
set eigensolver_type -genBandArpack

#Units

set N 1.
set sec 1.
set meter 1.
set mm [expr $meter/1000.]
set kN [expr $N*1000.]
set U 1.e2
set E [expr 2.03e+11*$N/pow($meter,2)]
set v 0.3
set G [expr $E/(2*(1+$v))]
set A 10000.0; #rigid 
set Abeam 10000.0; #rigid
set Atruss 10000.0; #rigid
set g 9.80665

set finish_massage "      Ground Motion Successfully Applied!"
set Time_in_finish_massage "End Time: "
set three_spaces "   "

set half_V_force [expr $weight/2]
 
file mkdir Results-HW5/$dataDir_EQ/; 			# create data directory
set GMdir "../GMfiles";			# ground-motion file directory
 
    node 1 0.0
	node 2 1
	node 3 1


# nodal masses:
mass 3 $mass;



# we need to set up parameters that are particular to the model.
set IDctrlNode 3;			# node where displacement is read for displacement control
set IDctrlDOF 1;			# degree of freedom of displacement read for displacement control

 


 # add the material to domain through the use of a procedure

 
 #omega calculation
set Pi [expr 2*asin(1.0)]; 		# define constants 
set wn [expr 2*$Pi/$Tn]; 		# natural circular frequency
set wn2 [expr pow($wn,2)];       # natural circular frequency squared
set Ki [expr $wn2*$mass];            # intial stiffness

 
#uniaxialMaterial CFSWSWP 1 $hight $width $fuf $tf $Ife $Ifi $ts $np $ds $Vs $screw_Spacing $nc $type $opening_Area $opening_Length 
uniaxialMaterial Elastic 1 $Ki
#uniaxialMaterial Steel01 1 $Fyield $Ki 0.001
uniaxialMaterial Elastic 2 2.03e+11

#element truss 1 1 2 1 1 -doRayleigh 1

 element truss 1 1 2 $Atruss 2
 element zeroLength 2 2 3 -mat 1 -dir 1 -doRayleigh 1


 # set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?

 fix 1 1

#recorder

 recorder Node -file Results-HW5/$dataDir_EQ/Xdisp-$EQ_ID-$model_id.out -time -node 3 -dof 1 disp
 recorder Element -file Results-HW5/$dataDir_EQ/XForce-$EQ_ID-$model_id.out -time -ele 2 localForce

#Additional_recorder

recorder Node -file Cd-CFS-SW/Results/by-models/X-for-auto-Excel/Xdisp-Cd-001.out -time -node 3 -dof 1 disp
recorder Element -file Cd-CFS-SW/Results/by-models/X-for-auto-Excel/XForce-Cd-001.out -time -ele 2 localForce

 

source Dynamic.EQ.Uniform.tcl

puts $model_id$three_spaces$EQ_ID$three_spaces$Time_in_finish_massage[getTime]$finish_massage