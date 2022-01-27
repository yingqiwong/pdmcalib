This directory stores calibration functions for the three-phase 
(solid-liquid-gas) mixture.

Each function takes in a [NPHS x Npts] matrix of phase fractions 
f, and any other pure-phase material properties. 

Each constraint has a range of validity, that must be defined
in the function. The output variable should be nan in the 
regions where the constraint is invalid.   