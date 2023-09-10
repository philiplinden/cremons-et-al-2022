This document describes the MATLAB scripts and spectra files associated with the simulations described in "Simulated Lunar Surface Hydration Measurements using Multispectral Lidar at 3 µm" by Cremons and Honniball (2022).

These scripts were created using MATLAB 2021a and were not tested on other versions.

Our goal was to simulate lunar spectra in the 1 to 4 micron wavelength range with a controllable hydration signature to determine if a lidar measuring at limited laser wavlelengths could constrain the total water abundance accurately and precisely. We developed this code to both simulate lunar spectra using laboratory endmember spectra and also to perform the lidar retrieval using a non-negative least squares algorithm. More details on the motivation and background for this project can be found in the open access article mentioned above.

How to Run:
1. Download all files to your local machine and place them together in a directory, noting the directory location.
2. Run MATLAB and open the main script: "Hydrated_Regolith_Spectra_Generator_and_Retrieval.m".
3. Modify line 7 of the script using the "cd" command to set the active directory to the location from Step 1 that includes the .txt files, .csv files, and other .m files.
4. Key variables to change are "NumSpectra" which controls how many simulations to run, abundance ranges for the endmember spectra, "WLS" which sets the wavelength range and spectral resolution of the simulations.
5. Run script
6. Four figures should be generated:
	First is a figure comparing the interpolated MORB spectra to the experimental MORB spectra. This figure is also present in the manuscript. 
	Second is a subset of the noisy lidar reflectance spectra that have been generated. By default this plot will include all spectra or the first 200 spectra, whichever is lower.
	Third is a scatter plot (input total water vs. retrieved total water) and histogram (total water error) of the mare results.
	Fourth is a scatter plot and histogram of the highlands results.

Modifying mixtures with new endmembers:
New endmembers can be included in mixtures by downloading data from the RELAB database (or other). Modifications to the source data may be required to change the data to wavelegth vs. reflectance. If the endmember spectra were measured in a geometry other than incidence angle = 30 degrees, esmisison angle = 0 degrees, phase angle = 30 degrees, the script "Hapke_Inverse_Function_Passive" will need to be modified. This includes the variables "mu0" (incidence angle), "mu" (emission angle), and "g" (phase angle).




