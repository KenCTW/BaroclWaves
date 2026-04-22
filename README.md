# BaroclWaves

Coastal trapped waves with stratification, topography and a mean flow


Detailed theory, examples and instructions can be found in the manual: ManjHWaves.pdf .

Each .jl file is a self-contained set of functions that can be used to calculate modal structures, dispersion curves and other properties. Their use also requires installation of various Julia packages which are listed in the manual.

Both of these files can be used to compute straight-coast coastal-trapped waves for subinertial frequencies.

The two ,jl files in this repository are:

jHwaves: Limited to real frequencies, and it can compute frictional and coupling coefficients when the coastal long wave approximation is imposed. An asymmetric Gaussian mean alongshore flow can be applied if desired.

jHWavesC: Allows a complex wave frequency, and so is applicable to problems involving lowest-order bottom friction and/or an unstable mean alongshore flow. The mean flow has an asymmetric Gaussian form.



Notation:
	j for Julia,
	H for hybrid (both stratification and topography are important),
	C for complex,
	M for Main

These programs are very similar to Matlab codes (bigr* and bigc*) that can be found at https://www2.whoi.edu/staff/kbrink/software/matlab-links/ .


-	K.H. Brink
