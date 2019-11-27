# Spectrometer-GenSpect
A matlab code that uses the GenSpect line-by-line code to get the mixing ratios of different atmospheric gases out of the count output of the Argus-2000 spectrometer. This is part of the MeznSat mission (a green house gas monitoring satellite equipped with an Argus-2000 spectrometer and an RGB camera)
To run it, add all the files to the Matlab path, Change the paths in looper.m to a path of your choosing, then run the script Main.m through Matlab.
The results will be saved in batches of 100 combination each. You can then run errorCalculation.m to get the best combination out of the results.
