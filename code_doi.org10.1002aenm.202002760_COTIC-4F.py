    # 2020, Nora Schopp, Univeristy of California, Santa Barbara

    # The following code consists of two parts, as indicated by 'Part 1' and 'Part 2' in the script.
    # Part 1 is the Optical Modeling software based on the original code from
    # 2010 George F. Burkhard, Eric T. Hoke, Stanford University. Minor modification were made. 
    # Copyright details of this part of the code can be found under "Part 1").

    #Part 2:
    #Copyright 2020, Nora Schopp, Viktor Brus, Univeristy of California, Santa Barbara
    
    #     This program is free software: you can redistribute it and/or modify
    #     it under the terms of the GNU General Public License as published by
    #     the Free Software Foundation, either version 3 of the License, or
    #    any later version. 
    #
    #     This program is distributed without any warranty and
    #     without even the implied warranty of MERCHANTABILITY or
    #     FITNESS FOR A PARTICULAR PURPOSE. See the
    #     GNU General Public License for more details.
    #     You should have received a copy of the GNU General Public License
    #     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    #     To refer to the approach demonstrated with this software,
    #     please cite the original pulication:
    #  
    
  #%% ______Part1__________:
  
# Modified by Nora Schopp, 12/2019, University of California, Santa Barbara
    
    #	'Copyright' 2012 Kamil Mielczarek (kamil.m@utdallas.edu), University of Texas at Dallas
    #	Modifications:
    #		10/2012 - Matlab code was ported to python, which is free and readily accessible to all.
    #                   3/2015 Improved accuracy of fundamental constants, fix plotting bugs on some platforms
    #
    #	Installation Notes:
    #		Requires use of the 2.x Branch of python and
    #		Matplotlib 	: http://matplotlib.org/
    #		Scipy		: http://www.scipy.org/
    #		Numpy		: http://numpy.scipy.org/
    #       Pylab             : wiki.scipy.org/PyLab
    #  
  
    # Original MATLAB CODE WAS :
    # Copyright 2010 George F. Burkhard, Eric T. Hoke, Stanford University
    
    #     This program is free software: you can redistribute it and/or modify
    #     it under the terms of the GNU General Public License as published by
    #     the Free Software Foundation, either version 3 of the License, or
    #     (at your option) any later version.
    #
    #     This program is distributed in the hope that it will be useful,
    #     but WITHOUT ANY WARRANTY; without even the implied warranty of
    #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #     GNU General Public License for more details.
    #
    #     You should have received a copy of the GNU General Public License
    #     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    # This program calculates the calculated short circuit current (assuming an
    # internal quantum efficiency of 100%) for a device where the thickness
    # of one of the layers is varied. The program uses the transfer matrix
    # method.
    
    # The procedure was adapted from J. Appl. Phys Vol 86, No. 1 (1999) p.487
    # and JAP 93 No. 7 p. 3693.
    # George Burkhard and Eric Hoke February 5, 2010
    
    # Build up the device, light enters from left to right, the layer names
    # correspond to the csv filenames in the 'matDir' ignoring the material prefix.
    # example : you have a material file named 'nk_P3HT.csv' the layer name would be 'P3HT'
    # Layer thicknesses are in nm.
 
import pylab as pl
from scipy.interpolate import interp1d
from math import ceil
from os.path import join,isfile
import numpy as np
  
savepath= "C:/Users/schop/Google Drive/Hecht_/COTIC-4F N101_/gRec-20190603-prefactor spectral baster 2.0-and simulation- current best fit/N235 20190912/"
    
layers 			= ['NS_Glassubstrate-E' , 'NS_James_ITO_Genosc', 'NS_ZnO-GenOsc-E-James', 'NS_NE9_COTIC4F-PCE10-James-ordinary' , 'NS_MoOx-James' , 'NS_20190603-Ag-E-with-T-glass-substrate-fitted-as-cauchy-a']
thicknesses		= [700000 , 130 , 42, 140 , 7 , 100 ] #for N183 d=93?

plotGeneration    = True   
plotWavelengths    = [400 , 500 , 600] 
  
activeLayer		= 3		# indexing starts from 0 
layerVaried		= 3		# indexing starts from 0


varyStart		= 0		# starting thickness of ›the layer to be varied (nm)
varyFinish		= 	thicknesses[activeLayer] # final thickness of the layer to be varied (nm)
varyStep		= 1		# thickness increments (nm)


lambda_start	= 301	# build a wavelength range to calculate over, starting wavelength (nm)
lambda_stop		= 1091	# final wavelength (nm)
lambda_step		= 5	# wavelength step size
lambda_diff        = ((lambda_stop - lambda_start)/lambda_step)+1 

x_step=1.0 # grid spacing of device cross section simulation in nm 
# this is where the n/k optical constants for each material are stored.
# file format is : CSV (comma separated values)
# each materials file is prefixed with 'nk_' and the first 'matHeader' number of lines are skipped.
matDir			= 'C:/Users/schop/Google Drive/OpticalModeling/TransferMatrix/matdata'	# materials data
matPrefix		= 'nk_'		# materials data prefix
matHeader		= 1				# number of lines in header

### PROGRAM BEGINS ###

# Constants
h 	= 	6.62606957e-34 	# Js Planck's constant
c	= 	2.99792458e8	# m/s speed of light
q	=	1.60217657e-19	#C electric charge

lambdas			= np.arange(lambda_start,lambda_stop+lambda_step,lambda_step,float)
varyThickness	= np.arange(varyStart,varyFinish+varyStep,varyStep,float)

## HELPER FUNCTIONS
def openFile(fname):
	"""
	opens files and returns a list split at each new line
	"""
	fd = []
	if isfile(fname):
		fn = open(fname, 'r')
		fdtmp = fn.read()
		fdtmp = fdtmp.split('\n')
		# clean up line endings
		for f in fdtmp:
			f = f.strip('\n')
			f = f.strip('\r')
			fd.append(f)
		# make so doesn't return empty line at the end
		if len(fd[-1]) == 0:
			fd.pop(-1)
	else:
		print("%s Target is not a readable file" % fname)
	return fd


def get_ntotal(matName,lambdas):
	fname = join(matDir,'%s%s.csv' % (matPrefix,matName))
	fdata = openFile(fname)[matHeader:]
	# get data from the file
	lambList	= []
	nList		= []
	kList		= []
	for l in fdata:
		wl , n , k = l.split(',')
		wl , n , k = float(wl) , float(n) , float(k)
		lambList.append(wl)
		nList.append(n)
		kList.append(k)
	# make interpolation functions
	int_n	= interp1d(lambList,nList)
	int_k	= interp1d(lambList,kList)
	# interpolate data
	kintList	= int_k(lambdas)
	nintList	= int_n(lambdas)
	# make ntotal
	ntotal = []
	for i,n in enumerate(nintList):
		nt = complex(n,kintList[i])
		ntotal.append(nt)
	return ntotal

def I_mat(n1,n2):
	# transfer matrix at an interface
	r = (n1-n2)/(n1+n2)
	t = (2*n1)/(n1+n2)
	ret = np.array([[1,r],[r,1]],dtype=complex)
	ret = ret / t
	return ret

def L_mat(n,d,l):
	# propagation matrix
	# n = complex dielectric constant
	# d = thickness
	# l = wavelength
	xi = (2*np.pi*d*n)/l
	L = np.array( [ [ np.exp(complex(0,-1.0*xi)),0] , [0,np.exp(complex(0,xi))] ] )
	return L

## / HELPER FUNCTIONS

# load AM1.5G Spectrum
am15_file = join(matDir,'AM15G.csv')
am15_data = openFile(am15_file)[1:]
am15_xData = []
am15_yData = []
for l in am15_data:
	x,y = l.split(',')
	x,y = float(x),float(y)
	am15_xData.append(x)
	am15_yData.append(y)
am15_interp = interp1d(am15_xData,am15_yData,'linear')
am15_int_y  = am15_interp(lambdas)

# initialize an array
n = np.zeros((len(layers),len(lambdas)),dtype=complex)

# load index of refraction for each material in the stack
for i,l in enumerate(layers):
	ni = np.array(get_ntotal(l,lambdas))
	n[i,:] = ni

# calculate incoherent power transmission through substrate

T_glass = abs((4.0*1.0*n[0,:])/((1+n[0,:])**2))
R_glass = abs((1-n[0,:])/(1+n[0,:]))**2

t 		= thicknesses
t[0] 	= 0 #glass does not hange much hen inlcuded but slows own the calculation. include fr final calculation but not for testing.
Jsc		= 0*varyThickness

for thickInd,thick in enumerate(varyThickness):
	t[layerVaried] = varyThickness[thickInd]
	# calculate transfer marices, and field at each wavelength and position
	t_cumsum	= np.cumsum(t)
	x_pos		= np.arange((x_step/2.0),sum(t),x_step)
	# get x_mat
	comp1	= np.kron(np.ones( (len(t),1) ),x_pos)
	comp2	= np.transpose(np.kron(np.ones( (len(x_pos),1) ),t_cumsum))
	x_mat 	= sum(comp1>comp2,0) 	# might need to get changed to better match python indices
	#
	R		= lambdas*0.0
	T		= lambdas*0.0
	E		= np.zeros( (len(x_pos),len(lambdas)),dtype=complex )
	# start looping
	for ind,l in enumerate(lambdas):
		# calculate the transfer matrices for incoherent reflection/transmission at the first interface
		S = I_mat(n[0,ind],n[1,ind])
		for matind in np.arange(1,len(t)-1):
			mL = L_mat( n[matind,ind] , t[matind] , lambdas[ind] )
			mI = I_mat( n[matind,ind] , n[matind+1,ind])
			S  = np.asarray(np.mat(S)*np.mat(mL)*np.mat(mI))
		R[ind] = abs(S[1,0]/S[0,0])**2
		T[ind] = abs((2/(1+n[0,ind])))/np.sqrt(1-R_glass[ind]*R[ind])
		# good up to here
		# calculate all other transfer matrices
		for material in np.arange(1,len(t)):
			xi = 2*np.pi*n[material,ind]/lambdas[ind]
			dj = t[material]
			x_indices	= np.nonzero(x_mat == material)
			x			= x_pos[x_indices]-t_cumsum[material-1]
			# Calculate S_Prime
			S_prime		= I_mat(n[0,ind],n[1,ind])
			for matind in np.arange(2,material+1):
				mL = L_mat( n[matind-1,ind],t[matind-1],lambdas[ind] )
				mI = I_mat( n[matind-1,ind],n[matind,ind] )
				S_prime  = np.asarray( np.mat(S_prime)*np.mat(mL)*np.mat(mI) )
			# Calculate S_dprime (double prime)
			S_dprime	= np.eye(2)
			for matind in np.arange(material,len(t)-1):
				mI	= I_mat(n[matind,ind],n[matind+1,ind])
				mL	= L_mat(n[matind+1,ind],t[matind+1],lambdas[ind])
				S_dprime = np.asarray( np.mat(S_dprime) * np.mat(mI) * np.mat(mL) )
			# Normalized Electric Field Profile
			num = T[ind] * (S_dprime[0,0] * np.exp( complex(0,-1.0)*xi*(dj-x) ) + S_dprime[1,0]*np.exp(complex(0,1)*xi*(dj-x)))
			den = S_prime[0,0]*S_dprime[0,0]*np.exp(complex(0,-1.0)*xi*dj) + S_prime[0,1]*S_dprime[1,0]*np.exp(complex(0,1)*xi*dj)
			E[x_indices,ind] = num / den
	# Absorption coefficient in 1/cm
	a = np.zeros( (len(t),len(lambdas)) )
	for matind in np.arange(1,len(t)):
		a[matind,:] = ( 4 * np.pi * np.imag(n[matind,:]) ) / ( lambdas * 1.0e-7 )
	#
	ActivePos = np.nonzero(x_mat == activeLayer)
	tmp1	= (a[activeLayer,:]*np.real(n[activeLayer,:])*am15_int_y)
	Q	 	= np.tile(tmp1,(np.size(ActivePos),1))*(abs(E[ActivePos,:])**2)
	# Exciton generatoion are
	Gxl		= (Q*1.0e-3)*np.tile( (lambdas*1.0e-9) , (np.size(ActivePos),1))/(h*c)
	if len(lambdas) == 1:
		lambda_step = 1
	else:
		lambda_step = (sorted(lambdas)[-1] - sorted(lambdas)[0])/(len(lambdas) - 1)
	Gx		= np.sum(Gxl,2)*lambda_step
	# calculate Jsc
	Jsc[thickInd] = np.sum(Gx)*x_step*1.0e-7*q*1.0e3

# Plotting Portion
fig0 = pl.figure(0)
fig0.clf()
ax0 = fig0.add_subplot(111)
ax0.set_ylim(ymin=0)
ax0.set_title('Current Density obtained from 100% IQE')
ax0.set_ylabel('Current Density , mA/cm$^2$')
ax0.set_xlabel('Layer Thickness , nm')
ax0.set_ylim(top=ceil(sorted(Jsc)[-1]))
ax0.plot(varyThickness,Jsc)
pl.savefig(savepath+'Current Density obtained from 100% IQE.tif')
fig0.show()


# overall Reflection from device with incoherent reflections at first interface
Reflection = R_glass+T_glass**2*R/(1-R_glass*R)
    
# Absorption coefficient in 1/cm
a = np.zeros( (len(t),len(lambdas)) )
for matind in np.arange(1,len(t)):
    a[matind,:] = ( 4 * np.pi * np.imag(n[matind,:]) ) / ( lambdas * 1.0e-7 )



#---------------

#---------------

if plotGeneration:
    
    fig2 = pl.figure(2)
    fig2.clf()
    ax2 = fig2.add_subplot(111)
    ax2.set_title('Fraction of Light Absorbed')
    ax2.set_xlabel('Wavelength (nm)')
    ax2.set_ylabel('Light Intensity Fraction')
    Absorption = np.zeros( (len(t),len(lambdas)) )
    for matind in np.arange(1,len(t)):
        Pos         = np.nonzero(x_mat == matind)
        AbsRate     = np.tile( (a[matind,:] * np.real(n[matind,:])),(len(Pos),1)) * (abs(E[Pos,:])**2)
        Absorption[matind,:] = np.sum(AbsRate,1)*x_step*1.0e-7
        ax2.plot(lambdas,Absorption[matind,:],label=layers[matind])
    ax2.plot(lambdas,Reflection,label='Reflection')
    art = []
    lgd = pl.legend(loc=9,  bbox_to_anchor=(0.5, -0.15), ncol=3)
    art.append(lgd)
    pl.savefig(
    savepath+'Fraction of Light Absorbed or Reflected.tif', dpi=600,
    additional_artists=art,
    bbox_inches="tight"
    )
    

    Absorption_T=Absorption.T

   
    
    # load and interpolate AM1.5G Data
    am15_file = join(matDir,'AM15G.csv')
    am15_data = openFile(am15_file)[1:]
    am15_xData = []
    am15_yData = []
    for l in am15_data:
        x,y = l.split(',')
        x,y = float(x),float(y)
        am15_xData.append(x)
        am15_yData.append(y)
    am15_interp = interp1d(am15_xData,am15_yData,'linear')
    am15_int_y  = am15_interp(lambdas)
    #
    ActivePos = np.nonzero(x_mat == activeLayer)
    tmp1    = (a[activeLayer,:]*np.real(n[activeLayer,:])*am15_int_y)
    Q         = np.tile(tmp1,(np.size(ActivePos),1))*(abs(E[ActivePos,:])**2)
    # Exciton generatoion are
    Gxl        = (Q*1.0e-3)*np.tile( (lambdas*1.0e-9) , (np.size(ActivePos),1))/(h*c)
    if len(lambdas) == 1:
        lambda_step = 1
    else:
        lambda_step = (sorted(lambdas)[-1] - sorted(lambdas)[0])/(len(lambdas) - 1)
    Gx        = np.sum(Gxl,2)*lambda_step

    

    # generation rate in whole device
    fig3 = pl.figure(3)
    fig3.clf()
    ax3 = fig3.add_subplot(111)
    ax3.set_title('Generation Rate in Device')
    ax3.set_xlabel('Position in Device (nm)')
    ax3.set_ylabel('Generation Rate/sec-cm$^3$')
    ax3.set_xlim(0,t_cumsum[-1])
    ax3.plot(x_pos[ActivePos[0]] , Gx[0])
    # Layer Bars
    for matind in np.arange(2,len(t)+1):
        xpos = (t_cumsum[matind-2]+t_cumsum[matind-1])/2.0
        ax3.axvline(np.sum(t[0:matind]),linestyle=':')
    ax3.text(xpos,sorted(Gx[0])[0],layers[matind-1],va='bottom',ha='center')
    pl.savefig(savepath+'Generation rate in device.tif')
    fig3.show()
    
   #Generation rate in active layer
    varyThickness_regrid = np.linspace(np.min(varyThickness), np.max(varyThickness), len(varyThickness)-1)
    # plot
    fig4 = pl.figure(4)
    fig4.clf()
    ax4 = fig4.add_subplot(111)
    ax4.set_title('Generation Rate in Active Layer')
    ax4.set_xlabel('Position in Active Layer (nm)')
    ax4.set_ylabel('Generation Rate/sec-cm$^3$')
    ax4.plot(varyThickness_regrid , Gx[0])
    pl.savefig(savepath+'Generation rate in active layer.tif')
    fig4.show()
    #End of Part 1
#%% ______Part2__________:
    
import pandas as pd
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #(necessary for 3D plot, do not delete even if it shows that it is unused)

#3D plot of generation rate
Gxl_df=pd.DataFrame(Gxl[0])
Gxl_2D= Gxl_df.values
x = varyThickness_regrid
y = lambdas
X, Y = np.meshgrid(x, y)
Z = np.transpose(Gxl_2D)
#X, Y = np.meshgrid(x, y)
fig = plt.figure(figsize=(10, 7))
ax = plt.axes(projection='3d')
ax.dist = 10
plot=ax.plot_surface(X,Y,Z, cmap=plt.cm.jet, linewidth=0, edgecolor=None, rstride=1, cstride=1)
fig.colorbar(plot, shrink=0.4, aspect=15, pad=0)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

ax.set_xlabel('x (nm)', size=12)
ax.set_ylabel(' $\lambda$(nm) ',size=12, )
ax.set_zlabel('G (cm$^{-3}$ s$^{-1}$)',size=12)
fig.savefig(savepath+'PCE10_COTIC-4F_3Dgenerationrate.tif', transparent=True, dpi=400)
plt.show()

#%%____________________INPUT:_______________________________________________
# measured EQE spectrum and AM1.5 spectrum in your directory:
path_EQE="C:/Users/schop/Google Drive/Data/PCE10-COTIC-4F, i/20190605-N230-N235-PCE10-COTIC-4F-for prefctor and impedance/EQE/EQE_N235_AExport_EQE.csv"
path_AM15 = path_AM15='C:/Users/schop/Google Drive/OpticalModeling/TransferMatrix/matdata/AM15G.csv'

    
gRec= 0.76 #<------- geminate recombination prefactor, needs to be determined for each blend system, see main text
V_bi=0.623 #V <-----can be obtained from JV curve fittitng, see SI
samples=100 # <---- define how many mu and tau values in the given range should be tested. 
#recommended: set mutau  range to be tested intitially wide by setting a wide mutau_exponent range (below), then based on the 'best fit' mutau that is given out by this program narrow your range in few iterations,
# this reduces computation time and can make the result mire accuruate
wl_stop=1100 #wavelength at which the fitting stops, could be smaller than what has been measured if high wl region is hard to fit

#e.g. test mutau values from 1E-5 to 1E-9: (plug in -9,-5, samples)
mutau_exponent= np.linspace(-9.4,-9, samples) #<--- set range in which to test value, start with wide range then narrow range down for results of higher accuracy
mutau= 10**(mutau_exponent) # the mutau values that will be tested

#___________________________________________________________________
L_cm=(thicknesses[activeLayer])*1e-7 #<----- thickness of device in cm
x=np.arange(0, thicknesses[activeLayer])
x_cm=x*1e-7

wavelength_up_to_stop= lambdas[np.where(lambdas<=wl_stop)]
stop_fit_index=len(wavelength_up_to_stop)
wavelength = lambdas[:stop_fit_index]

EQE_array_raw=pd.read_csv((path_EQE), delimiter=',',header=0).values #importing experimental EQE
all_wavelength_EQE= np.around(EQE_array_raw[:,0])
all_EQE=EQE_array_raw[:stop_fit_index,1]/100

EQE_measured_rough=np.zeros([len(wavelength), 1]) # bring EQE measured into the same shape as calculated EQE (determined by optical modeling wavelength)
wavelength_measured=np.zeros([len(wavelength), 1])
for j in range (0,len(wavelength)):
                EQE_measured_rough[j:,0]=all_EQE[np.where(all_wavelength_EQE==wavelength[j])]

from scipy.signal import savgol_filter
EQE_measured = savgol_filter(np.ravel(EQE_measured_rough), 21, 7 ) #smoothing measured EQE data, not always necessary

   
df_AM15=pd.read_csv((path_AM15), delimiter=',',  header=None) #importing AM1.5 spectrum
AM15=df_AM15.values#convert pandas dataframe into np array
I_fullspectrum= np.asarray(AM15[1:,1], dtype=float)
wl_fullspectrum= np.asarray(AM15[1:,0], dtype=float)
wl_and_I_needed=np.zeros([len(wavelength), 2]) #Filter: get the intensities from AM 15 file for only those wavelength used in Optical Modeling
for j in range (0,len(wavelength)):
                wl_and_I_needed[j:,0]=wl_fullspectrum[np.where(wl_fullspectrum==wavelength[j])]
                wl_and_I_needed[j:,1]=I_fullspectrum[np.where(wl_fullspectrum==wavelength[j])]
        
I=wl_and_I_needed[:,1] #intensities of AM1.5 for the wavelength used in OM
photonflux=(I*1e-3)/((1240/wavelength)*sc.e)
L_nm= L_cm*1e7 

#%%________________Theoretical EQE Calculation Loop or 'Hecht Calulation loop', EQE_hecht = EQE that has been calulated with Hecht equation, Pg for all mutau___________________________________________________
Jsc_ideal=np.zeros([len(wavelength)]) # ideal = for 100% extraction
EQE_ideal=np.zeros([len(wavelength)])
Jsc_ideal_reduced_by_gem_rec=np.zeros([len(wavelength)]) # ideal = for 100% extraction
EQE_ideal_reduced_by_gem_rec=np.zeros([len(wavelength)])

Every_f = np.zeros([len(wavelength),len(x)])
Every_eta= np.zeros([len(mutau) ,len(x)])
Every_eta_n= np.zeros([len(mutau) ,len(x)])
Every_eta_p= np.zeros([len(mutau) ,len(x)])
Every_w_nm= np.zeros([len(mutau)])
Every_EQE_hecht= np.zeros([len(mutau) ,len(wavelength)])

Every_deviation= np.zeros([len(mutau), len(wavelength)]) #1
Every_average_deviation= np.zeros([len(mutau), 2]) 

best_mutau_for_each_wavelength=np.zeros([2, len(wavelength)]) #2
Every_EQE_hecht_smooth= np.zeros([len(mutau), len(wavelength)]) #3

for l in np.arange(0, len(wavelength)):
     for mt in np.arange(0, len(mutau)):
         for p in np.arange(0, len(x)):
             
             Jsc_ideal[l]= np.sum(Gxl.T[l,:])*1e-7*sc.e #100% extaction and no geminate recmbination
             EQE_ideal[l]= (Jsc_ideal[l]/sc.e)/photonflux[l]
             
             Jsc_ideal_reduced_by_gem_rec[l]= gRec*np.sum(Gxl.T[l,:])*1e-7*sc.e #includes geminate recombinatin prefactor g, determinded separetly
             EQE_ideal_reduced_by_gem_rec[l]= (Jsc_ideal_reduced_by_gem_rec[l]/sc.e)/photonflux[l]
             
             w_cm = mutau[mt]*(V_bi)/(L_cm) 
             w_nm=w_cm*1e7 # path 
             Every_w_nm[mt]=w_nm
             
             eta_n= w_nm/L_nm*(1-np.exp(-(L_nm-x[p])/w_nm))
             eta_p= w_nm/L_nm*(1-np.exp(-x[p]/w_nm))                           
             eta= eta_n+eta_p #extraction efficiency
                 
             Every_eta[mt, p]=   eta
             Every_eta_n[mt, p]=   eta_n
             Every_eta_p[mt, p]=   eta_p
             Every_f[l,p]=Gxl.T[l,p]*eta #extracted charge carriers
             Jsc_hecht= gRec * np.sum(Every_f[l,:])*1e-7*sc.e #now integrated over x and only wavelength, mu, tau dependant
             EQE_hecht=(Jsc_hecht/sc.e)/photonflux[l]
             Every_EQE_hecht[mt, l]= EQE_hecht
             
             deviation = np.absolute(EQE_hecht - EQE_measured[l]) 
             Every_deviation[mt, l] = deviation #mutau values * len(wavelength)
             
             #get best EQE based on minimizing the overall deviation at each wavelength from measured:
             average_deviation= np.sum(Every_deviation[mt,:])/len(wavelength)
             Every_average_deviation[mt,0] = mutau[mt]
             Every_average_deviation[mt,1] = average_deviation #array to correlate average deviation with corresponding mu and tau values
            
#%%________ find smallest average deviation for all mu-tau-products__________________________________________________________            
#_____least square fit____
Every_r_squared= np.zeros([len(mutau)])
for mt in np.arange(0, len(mutau)):
    residuals= EQE_measured- Every_EQE_hecht[mt,:]
    ss_res=np.sum(residuals**2)
    ss_tot=  np.sum((Every_EQE_hecht[mt,:]-np.mean(EQE_measured))**2)
    r_squared_loop = 1 - (ss_res / ss_tot)
    Every_r_squared[mt]=r_squared_loop
    
index_best_fit_1=np.argmax(Every_r_squared)
mutau_best_1= mutau[index_best_fit_1]
w_nm_best_1= mutau_best_1*(V_bi)/L_cm
EQE_hecht_best_1=Every_EQE_hecht[index_best_fit_1,:] #the EQE that gives the best fit
r_squared = Every_r_squared[index_best_fit_1]


#____output and__plots______________________________________________________________
fig=plt.figure()
plt.title('EQE calculated and measured, best fits', {'color': 'black', 'fontsize': 10})
plt.plot(wavelength, EQE_measured, color='black',linewidth= 2, label='measured')
plt.plot(wavelength,  EQE_ideal,'navy',  linestyle='-.' ,linewidth= 2, label=r'$\eta$=1, $P_g$=1')
plt.plot(wavelength,  EQE_ideal_reduced_by_gem_rec,'yellowgreen',linewidth= 2,  linestyle='-.', label=(r'$\eta$=1, $P_g$='+str(gRec)))
plt.xlim(290, 1100)
plt.xlabel('wavelength (nm)', {'color': 'black', 'fontsize': 10})
plt.ylabel('EQE', {'color': 'black', 'fontsize': 10})
plt.ylim(0, 0.9)#plt.plot(wavelength,  EQE_hecht_best_2, color= 'darkorange',  label='Hecht, Approach 2')
plt.plot(wavelength,  EQE_hecht_best_1, color= 'red', linestyle='-.', linewidth= 2, label='Fit')# approahx 1b #from A to mA
art = []
lgd = pl.legend(loc=9,  bbox_to_anchor=(0.5, -0.15), ncol=3)
art.append(lgd)
pl.savefig(
savepath+str(thicknesses)+'EQE_calc.tif', dpi=400,
additional_artists=art,
bbox_inches="tight"
)

output_EQEs= pd.DataFrame({'wavelength in nm ':wavelength, 'EQE measured':EQE_measured, 'EQE eta=1, Pg=1':EQE_ideal, ("EQE eta=1, Pg="+str(gRec)):EQE_ideal_reduced_by_gem_rec, 'EQE Simulation Fit':EQE_hecht_best_1})
output_EQEs.to_csv(savepath+'calculated_EQEs.csv') 

#%%
Jsc_Am15= np.sum(EQE_measured*photonflux*5*sc.e *1000)

print("Jsc from EQE Am1.5=")
print(Jsc_Am15)

#%%#extraction efficiency at 0 V
Every_eta_0V= np.zeros([len(x)])
Every_eta_n_0V= np.zeros([len(x)])
Every_eta_p_0V= np.zeros([len(x)])

for p in np.arange(0, len(x)):
    w_cm_0V = mutau_best_1*(V_bi)/(L_cm) #Vbi plus applies voltage
    w_nm_0V=w_cm_0V*1e7 #
    eta_n_0V= w_nm_0V/L_nm*(1-np.exp(-(L_nm-x[p])/w_nm_0V))
    eta_p_0V= w_nm_0V/L_nm*(1-np.exp(-x[p]/w_nm_0V))                           
    eta_0V= eta_n_0V+eta_p_0V
    Every_eta_0V[p]= eta_0V
    Every_eta_n_0V[p]= eta_n_0V
    Every_eta_p_0V[p]= eta_p_0V

average_extraction_eff_at_0_V= np.sum(Every_eta_0V)/len(Every_eta_0V)
 
plt.title('Extraction efficiencies at 0 V, inerted structure', {'color': 'black', 'fontsize': 10})
plt.plot(x, Every_eta_n_0V, label='electron' ) #choose any index of the samples just as exmple visualiztion
plt.plot(x, Every_eta_p_0V, label='hole')
plt.plot(x, Every_eta_0V, label='sum')
plt.ylim(0,1.05)
plt.xlabel('Postion in active layer (distance from Cathode) in nm' )
plt.ylabel('extarction efficency')
plt.legend()
plt.savefig(savepath+'Extraction Efficency at 0 V.tif', dpi=300)
plt.show()
#%%
#extraction efficiency at -3V
Every_eta_m3V= np.zeros([len(x)])
Every_eta_n_m3V= np.zeros([len(x)])
Every_eta_p_m3V= np.zeros([len(x)])

for p in np.arange(0, len(x)):
    w_cm_m3V = mutau_best_1*(3+V_bi)/(L_cm) #Vbi plus applies voltage
    w_nm_m3V=w_cm_m3V*1e7 #
    eta_n_m3V= w_nm_m3V/L_nm*(1-np.exp(-(L_nm-x[p])/w_nm_m3V))
    eta_p_m3V= w_nm_m3V/L_nm*(1-np.exp(-x[p]/w_nm_m3V))                           
    eta_m3V= eta_n_m3V+eta_p_m3V
    Every_eta_m3V[p]= eta_m3V
    Every_eta_n_m3V[p]= eta_n_m3V
    Every_eta_p_m3V[p]= eta_p_m3V

average_extraction_eff_at_minus_3_V= np.sum(Every_eta_m3V)/len(Every_eta_m3V)
 
plt.title('Extraction efficiencies that lead to best fit, Appr. 1, at -3 V', {'color': 'black', 'fontsize': 10})
plt.plot(x, Every_eta_n_m3V, label='electron' ) #choose any index of the samples just as exmple visualiztion
plt.plot(x, Every_eta_p_m3V, label='hole')
plt.plot(x, Every_eta_m3V, label='sum')
plt.ylim(0,1.05)
plt.xlabel('Postion in active layer (distance from Cathode) in nm' )
plt.ylabel('extarction efficency')
plt.legend()
plt.savefig(savepath+'Extraction Efficency at -3 V.tif', dpi=300)
plt.show()

#%%#%%
#(spatially averaged) extraction efficiency at 0V with optimized mutau

L_nm=np.arange(0,200)
Every_average_extraction_eff_vary_L= np.zeros([len(L_nm)])

for l in np.arange(0, len(L_nm)):
    
    x= np.arange(0, l) #these are the possible positions for each L_nm
    
    Every_eta_0V= np.zeros([len(x)])
    Every_eta_n_0V= np.zeros([len(x)])
    Every_eta_p_0V= np.zeros([len(x)])
    
    for p in np.arange(0, len(x)):
        w_cm_0V = mutau_best_1*(V_bi)/(L_nm[l]*1e-7) #Vbi plus applies voltage
        w_nm_0V=w_cm_0V*1e7 #
        eta_n_0V= w_nm_0V/L_nm[l]*(1-np.exp(-(L_nm[l]-x[p])/w_nm_0V))
        eta_p_0V= w_nm_0V/L_nm[l]*(1-np.exp(-x[p]/w_nm_0V))                           
        eta_0V= eta_n_0V+eta_p_0V
        Every_eta_0V[p]= eta_0V
        Every_eta_n_0V[p]= eta_n_0V
        Every_eta_p_0V[p]= eta_p_0V
    
        average_extraction_eff_at_0_V= np.sum(Every_eta_0V)/len(Every_eta_0V)
        
    Every_average_extraction_eff_vary_L[l]=average_extraction_eff_at_0_V
    
 
plt.title('Extraction efficiencies at different thickness', {'color': 'black', 'fontsize': 10})
plt.plot(L_nm[1:], Every_average_extraction_eff_vary_L[1:], label='eta')
plt.xlabel('Thickness of active layer (nm)' )
plt.ylabel('extraction efficency')
plt.legend()
#plt.savefig(savepath+'Ex_varyLtraction Efficency at -3 V.tif', dpi=300)
plt.show()
        

#%%
output= ["Simulation Fit N.S.:", "mobility-lifeftime-product (cmˆ2/Vs): " , "{:.2e}".format(mutau_best_1), 
         'average_extraction_eff_at_minus_3_V approach 1 =',
         str("{:.3f}".format(average_extraction_eff_at_minus_3_V)), "r2=", str("{:.2f}".format(r_squared)) ]

np.savetxt(savepath+'Output-tau_and_mu_values.txt', output , fmt='%s')