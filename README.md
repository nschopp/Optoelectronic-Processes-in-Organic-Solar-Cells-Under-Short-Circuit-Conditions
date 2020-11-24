# Optoelectronic-Processes-in-Organic-Solar-Cells-Under-Short-Circuit-Conditions

Please reference the following work:
N. Schopp, V. Brus*, J. Lee, G. C. Bazan, T.-Q. Nguyen*, Advanced Energy Materials 2020, DOI 10.1002/aenm.202002760.


The python-based software created for this work operates in conjunction with the transfer matrix method software that was developed by George F. Burkhard, Eric T. Hoke at Stanford University.[2]
Our part of the code includes:
the visualization of the 3D Generation rate
the calculation of the EQE in dependence on the position x in the device, the wavelengths and the μτ
fitting of the calculated EQEs to the experimental EQE and the output of the μτ
the calculation of the Jsc from measured EQE for AM1.5 illumination
a visualization of the extraction efficiency η(x) of electrons and holes at 0 V in dependence of the position in the device based on the determined μτ
a visualization of the extraction efficiency η(x) of electrons and holes at -3 V in dependence of the position in the device based on the determined μτ
the average extraction efficiency η ̅ for different active layer thicknesses (required for device optimization)

To use the code for device optimization, follow the step-by-step guide below.

- Recommendations for Python Beginners:

Install Anaconda and Spyder, import the .py (code) file into Spyder.
Download the required nk data in one of your folders, as well as the AM1.5 file and the EQE file. 
You can replace the EQE file with your own EQE file and use your own optical properties as desired. 
Please keep in mind that your files need to be exactly in the same format for the code to work without modification (eg. csv or tab delimiter, number of lines in header, no empty line at the bottom of the file). 

Make sure you plug in the correct file names (without the nk_ prefix) in the header of the trasnfer matrix code.
Plug in your file directories in the code (pathcopycopy is a useful tool to easily obtain the path of any folder you are in on the windows file explorer)
l 87 savepath= "..." (where the program saves output data and graphics)
l 435 path_EQE="" (your experimental EQE data)
  path_AM15 = path_AM15=' ' (spectrum file)


- Step-by Step Guide for the Application of our Approach to the Device Optimization:

After obtaining the optical constants of the active layer and of other layers that are not provided in this work or other literature, the following steps lead to the prediction of the thickness-dependent Jsc.

Experimental steps:
Fabrication of at least one solar cell  (thickness measurement should be done after completing the other experimental steps) 
Measurement of the JV characteristics at 1 sun and in the dark (to obtain Vbi) 
Measurement of the EQE
Measuring JV-characteristics under monochromatic illumination

Simulation and calculation steps:

Simulating the Jsc,theo under the monochromatic illumination used in step 4. (It is important to set other wavelengths to zero in the spectral distribution file.) 
Calculation of Pg 
Using the provided software (see below) to obtain μτ and the extraction efficiency η by fitting to measured EQE.
Calculation of the extraction probability η(x) for each thickness L of interest, then averaging η(x) to get extraction efficiency η ̅  for each L.
Simulation of the Jsc,ideal for each L (no losses) with TMM under 1 sun AM 1.5 or the spectrum of your interest. (blue curve)
Product of Pg, η(L) and simulated Jsc,ideal give predicted Jsc. (red curve)

The combined software can be accessed publicly and used under the GNU license agreement after the date of publication of this work. It can be found together with the optical constants provided in the manuscript under https://github.com/nschopp/Optoelectronic-Processes-in-Organic-Solar-Cells-Under-Short-Circuit-Conditions. 

