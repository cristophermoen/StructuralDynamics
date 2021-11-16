# --------------------------------------------------------------------------------------------------
# Verified



# DYNAMIC ground-motion analysis -------------------------------------------------------------
# create load pattern
 


#SFDE_here
set SF_DE_for_EQ [expr 1*1]
#tohere
#puts $SF_DE_for_EQ
#set GMfatt [expr $NF*$SF_DE_for_EQ*$g]
#set GMfatt [expr $NF*$SF_DE*$g]

set GMfatt [expr 1*$g]

#set GMfatt [expr $NF*$SF_DE*$g];
set GMfile "ElCentro.txt" ;

timeSeries Path 2 -dt $dt -filePath $GMfile -factor $GMfatt; # define acceleration vector from file (dt=0.005 is associated with the input file gm)
pattern UniformExcitation 2 1 -accel 2;		         # define where and how (pattern tag, dof) acceleration is applied

# set damping based on first eigen mode
set freq [expr [eigen -fullGenLapack 1]**0.5]
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq]


set DtAnalysis	[expr 0.01*$sec];	# time-step Dt for lateral analysis
set TmaxAnalysis	[expr $Analysis_Duration*$sec];	# maximum duration of ground-motion analysis -- should be 50*$sec
set Nsteps [expr int($TmaxAnalysis/$DtAnalysis)];


# create the analysis
wipeAnalysis;					     # clear previously-define analysis parameters
constraints Plain;     				 # how it handles boundary conditions
numberer Plain;					     # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;					 # how to store and solve the system of equations in the analysis
algorithm Linear					 # use Linear algorithm for linear analysis
integrator Newmark 0.5 0.25 ;	     # determine the next time step for an analysis
analysis Transient;					 # define type of analysis: time-dependent
analyze $Nsteps $DtAnalysis;					 # apply 3995 0.01-sec time steps in analysis



