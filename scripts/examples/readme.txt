
GENSPECT EXAMPLES:

This folder contains script files or user interfaces that demonstrate how the 
GENSPECT tools are used in line-by-line calculations. Before running any of these 
scripts, be sure that the '\\file_path\genspect' folder is also on your MATLAB path, as the 
scripts here have function calls to numerous GENSPECT functions.

-BENCHTST: 				script file showing line-by-line calculation of absorption coefficients
-EXPT_CELL_EXAMPLE: 	script illustrating calculation of radiance through a simple cell with a Sun-like source
-REFLECTED_PATH_EXAMPLE:script file showing a radiance calculation for a reflected path [note that this
						function does not include a "diffusivity" factor and therefore will not accurately
						account for emission reflection from a lambertian scatterer.
-DOWNWARD_PATH_EXAMPLE: script for simple downward radiance computation.
-LIMB_VIEW_EXAMPLE: 	Script with limb viewing geometry.
-PRESS_VS_ALT_EXAMPLE: 	Script that compares a pressure division atmoshere calculation with 
						an altitude division.

(C) Ben Quine, 2001.