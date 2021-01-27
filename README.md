# Optoelectronic-Processes-in-Organic-Solar-Cells-Under-Short-Circuit-Conditions

Please reference the following work:
N. Schopp, V. Brus*, J. Lee, G. C. Bazan, T.-Q. Nguyen*, Advanced Energy Materials 2020, DOI 10.1002/aenm.202002760.


The python-based software created for this work operates in conjunction with the transfer matrix method software that was developed by George F. Burkhard, Eric T. Hoke at Stanford University.

Our part of the code includes:
- visualization of the 3D Generation rate
- calculation of the EQE in dependence on the position x in the device, the wavelengths, and the <img src="https://latex.codecogs.com/svg.latex?\mu\tau"> 
- fitting of the calculated EQEs to the experimental EQE and the output of the <img src="https://latex.codecogs.com/svg.latex?\mu\tau"> 
- calculation of the J<sub>sc</sub> from measured EQE for AM1.5 illumination
- visualization of the extraction efficiency <img src="https://latex.codecogs.com/svg.latex?\eta(x)">  of electrons and holes at 0 V in dependence of the position in the device based on the determined μτ
- visualization of the extraction efficiency <img src="https://latex.codecogs.com/svg.latex?\eta(x)"> of electrons and holes at -3 V in dependence of the position in the device based on the determined μτ
- the average extraction efficiency <img src="https://latex.codecogs.com/svg.latex?\bar{\eta}"> for different active layer thicknesses (required for device optimization)

To use the code for device optimization, follow the step-by-step guide below.

## Recommendations for Python Beginners

Install Anaconda and Spyder, import the .py (code) file into Spyder.
Download the required nk data in one of your folders, as well as the AM1.5 file and the EQE file. 
You can replace the EQE file with your own EQE file and use your own optical properties as desired. 
Please keep in mind that your files need to be exactly in the same format for the code to work without modification (eg. csv or tab delimiter, number of lines in header, no empty line at the bottom of the file).

By default, running the code will output the graphs into the root directory of the repository.

### Alternative using `pipenv`

To install `pipenv` (a Python package and virtual environment manager), follow the instructions for your operating system [here](https://pipenv.pypa.io/en/latest/install/).

Download the repository from Github, then navigate to the root directory in a terminal prompt. Run `pipenv install` to automatically create the virtual environment and install the packages from `Pipfile`. After installing, run `pipenv shell` to activate the environment.


## Step-by Step Guide for the Application of our Approach to the Device Optimization

After obtaining the optical constants of the active layer and of other layers that are not provided in this work or other literature, the following steps lead to the prediction of the thickness-dependent J<sub>sc</sub>.

**Experimental steps:**

Fabrication of at least one solar cell (thickness measurement should be done after completing the other experimental steps) 

Measurement of the JV characteristics at 1 sun and in the dark (to obtain V<sub>bi</sub>) 

Measurement of the EQE

Measuring JV-characteristics under monochromatic illumination

Simulation and calculation steps:

**Simulating the J<sub>sc</sub>**

theo under the monochromatic illumination used in step 4. (It is important to set other wavelengths to zero in the spectral distribution file.)

**Calculation of Pg**

Using the provided software (see below) to obtain μτ and the extraction efficiency η by fitting to measured EQE.

Calculation of the extraction probability η(x) for each thickness L of interest, then averaging η(x) to get extraction efficiency <img src="https://latex.codecogs.com/svg.latex?\bar{\eta}"> for each L.
Simulation of the Jsc, ideal for each L (no losses) with TMM under 1 sun AM 1.5 or the spectrum of your interest. (blue curve)

Product of Pg, η(L) and simulated J<sub>sc</sub>,ideal give predicted J<sub>sc</sub>. (red curve)
