# Optoelectronic-Processes-in-Organic-Solar-Cells-Under-Short-Circuit-Conditions

Please reference the following work:
N. Schopp, V. Brus*, J. Lee, G. C. Bazan, T.-Q. Nguyen*, Advanced Energy Materials 2020, DOI 10.1002/aenm.202002760.

And reference for the Transfer matrix part 1:
G. F.[1]G. F. Burkhard, E. T. Hoke, M. Group, 2011, 6.
Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater. 2010, 22, 3293.


To use the code for device optimization, follow the step-by-step guide in the SI.

Recommendations for Python Beginners:

Install Anaconda and Spyder, import the .py (code) file into Spyder.
Download the required nk data in one of your folders, as well as the AM1.5 file and the EQE file. 
You can replace the EQE file with your own EQE file and use your own optical properties as desired. 
Please keep in mind that your files need to be exactly in the same format for the code to work without modification (eg. csv or tab delimiter, number of lines in header, no empty line at the bottom of the file). 
Make sure you plug in the correct file names (without the nk_ prefix) in the header of the trasnfer matrix code.
Plug in your file directories in the code (pathcopycopy is a useful tool to easily obtain the path of any folder you are in on the windows file explorer)
l 87 savepath= "..." (where the program saves output data and graphics)
l 435 path_EQE="" (your experimental EQE data)
  path_AM15 = path_AM15=' ' (spectrum file)

