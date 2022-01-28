#python 3
#2019-2020
#CC-By-SA
#Old code from Carron Leopold
#Modified by Damien Legros, Cédric Cornede, Arnaud Quelin, Rouquaya Mouss, Hamid Hachemi

#####
#Imports

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats
from scipy.spatial import distance
from tqdm import tqdm
from hmmlearn import hmm
import seaborn as sns

np.seterr(divide = 'ignore') 

#####
#ADDITIONAL LOADER

def EpiGbyres(EpiGfilename,res,achr,sizeatres,NbEpi):
	"""
	Generate a matrix of repartition of epiG in a chr at the given resolution
	Do it the hard way : parse the file to the desired resolution directly
	
	Entries : EpiGfilename -> String
		  res -> int
		  achr -> int
		  sizeatres -> int
		  NbEpi -> int
	Return : Epimat -> list
	"""
	Epimat=np.zeros((NbEpi,sizeatres))
	o=open(EpiGfilename,'r')
	l=o.readline()
	while l:
		ls=l.split()
		if ls[0]==achr:
			#force to scale bp cover in bin cover here
			EpiV=int(ls[3])-1
			b=np.float(ls[1])
			e=np.float(ls[2])
			begin=int(np.floor(b/res)) #force to cast for indexing error
			end=int(np.floor(e/res)) #force to cast for indexing error
			eC=end*res
			bC=begin*res
			val=1
			if begin==end: #1 bin
				Epimat[EpiV,begin]+=val*(e-b)
			else:
				if (end-begin)==1: #2 bin
					Epimat[EpiV,begin]+=val*(eC-b)
					Epimat[EpiV,end]+=val*(e-eC)
				else: #more than 2 bin
					Epimat[EpiV,begin]+=val*(res-(b-bC))
					while begin<end:
						Epimat[EpiV,begin]+=val*res
						begin+=1
					if (e-eC)>0:
						Epimat[EpiV,begin]+=val*(e-eC) #here begin=end
		l=o.readline()
	o.close()
	return Epimat

#####
#Building tools

def buildMatrix(HiCfilename, printer=True):
    """
    Build the matrix from the HiC data
    
    Entries : HiCfilename -> String
	      printer -> boolean
    Return : A -> list
    """
    if printer:
        bar = tqdm(range(1), desc="Loading Matrix ")
    else:
        bar = range(1)

    for i in bar:
        pass
    try:
        A=np.loadtxt(HiCfilename)
        A=np.int_(A)
    except:
        print("Data was not correctly downloaded, please delete the folder and relaunch the main.py")
    
    #Build array at pb resolution
    A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)
    A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
    return A

def buildColors(EpiGfilename, chromosome_number, LENTEST, printer=True):
	"""
	Build the colors for the epigenetic marks
	
	Entries : HiCfilename -> String
		  chromosome_number -> int
		  LENTEST -> int
	      	  printer -> boolean
    	Return : color_vec -> list
	"""
	if printer:
		bar = tqdm(range(1), desc="Loading Colors ")
	else:
		bar = range(1)

	for i in bar:
		pass
	color=pd.read_csv(EpiGfilename,delimiter='\t',header=None,names=[1,2,3,4])
	#Take only chr of interest
	color=color[color[1]==chromosome_number]
	#Number of color in the file
	number=color[4].max()
	#Build array at pb resolution LENchr * number of color
	color_vec=np.zeros((LENTEST,number+1), dtype='uint8') 
	i=0
	while i<np.shape(color)[0]:
		color_vec[color[2].iloc[i]:color[3].iloc[i],color[4].iloc[i]]=1
		i+=1
	return color_vec

#####
#Plotting tools

def plotter(data2D, nameFig="Plot", plotType="Matrix", cmap="hot_r", vmin=None, vmax=None, jupyter=False, nameFile="Plot.pdf", centro_start=None, centro_end=None, Data_type='Contact'):
    """
    Plotter used in the analysis
    
    Entries : data2D -> list
	      nameFig -> String
	      plotType -> String
	      cmap -> String
	      vmin -> float
	      vmax -> float
	      jupyter -> boolean
	      nameFile -> String
	      centro_start -> int
	      centro_end -> int
	      Data_type -> String
    
    Make the various following plots :
    - Matrix plot
    - AB Compartments plot
    - Density plot
    - HMM Score plot
    - All Results plot
    """
    if plotType == "Matrix":
        custom_params = {"axes.spines.right": False, 
                         "axes.spines.top": False}
        sns.set_theme(style="ticks", rc=custom_params)
        fig = plt.figure()
        fig.set_dpi(300)
        im = plt.imshow(data2D, cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, drawedges=False)
        plt.xlabel("Position")
        plt.ylabel("Position")
        plt.title(nameFig)

    elif plotType == "AB":
        custom_params = {"axes.spines.right": False, 
                         "axes.spines.top": False,}
        sns.set_theme(style="ticks", rc=custom_params)
        x = np.arange(0.0, np.shape(data2D)[0], 1)
        y = data2D
        fig, ax = plt.subplots()
        fig.set_dpi(300)
        ax.fill_between(x, y, where=(data2D>0), color='red')
        ax.axhline(0, color="black", linewidth=0.75)
        ax.fill_between(x, y, where=(data2D<0), color='blue')
        plt.axvspan(centro_start, centro_end, facecolor='red', alpha=0.2)
        plt.xlabel("Position")
        plt.ylabel("Activation")
        plt.title(nameFig)

    elif plotType == "Density": #Density plot
        custom_params = {"axes.spines.right": False, 
                         "axes.spines.top": False}
        sns.set_theme(style="ticks", rc=custom_params)
        x = np.arange(0,len(data2D))
        y = data2D
        fig = plt.figure()
        fig.set_dpi(300)
        plt.plot(x, y, color='blue', alpha=1, linewidth=0.1)
        plt.axvspan(centro_start, centro_end, facecolor='red', alpha=0.2)
        plt.xlabel("Position")
        plt.ylabel("Density")
        plt.title(nameFig)

    elif plotType == "HMMScore":
        custom_params = {"axes.spines.right": False, "axes.spines.top": False}
        sns.set_theme(style="ticks", rc=custom_params)
        x = np.arange(2,16)
        y = data2D
        fig = plt.figure()
        fig.set_dpi(300)
        plt.plot(x, y, color='blue', alpha=0.4)
        best = 2
        if (Data_type=='Contact'):
            threshold = 1
        elif Data_type=='Epigenetic':
            threshold = 0.5
        for i in range(1,len(data2D)):
            temp = ((data2D[i]-data2D[i-1])/np.abs(data2D[i]))*100
            if temp >= threshold:
                best += 1
            else: 
                break
        plt.scatter(best, data2D[best-2], alpha=1, color='red', label="predicted number of compartments : " + str(best))
        plt.xlabel("Number of compartments")
        plt.ylabel("HMM Score")
        plt.legend(loc='upper left')
        plt.title(nameFig)

    elif plotType == "Barcode":
        custom_params = {"axes.spines.left": False, 
                         "axes.spines.right": False, 
                         "axes.spines.top": False}
        sns.set_theme(style="ticks", rc=custom_params)
        pixel_per_bar = 1
        dpi = 100
        fig = plt.figure(figsize=(len(data2D) * pixel_per_bar / dpi, 2))
        fig.set_dpi(300)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(data2D.reshape(1, -1), cmap='bwr', aspect='auto', interpolation='nearest')
        ax.axes.get_yaxis().set_visible(False)
        plt.xlabel("Position")
        plt.title(nameFig)

    elif plotType == "Visualization":
        [barcode, corr, density] = data2D
        fig = plt.figure(constrained_layout=True)
        fig.set_dpi(300)
        gs = fig.add_gridspec(ncols=2, 
                              nrows=2,
                              width_ratios=[10, 6.25], 
                              wspace=0.,
                              hspace=0., 
                              height_ratios=[1, 13])
        custom_params = {"axes.spines.right": False, 
                         "axes.spines.top": False}
        sns.set_theme(style="ticks", rc=custom_params)

        #Barcode
        f_ax1 = fig.add_subplot(gs[0, 0])
        f_ax1.set_xlim([0,np.shape(barcode)[0]])
        f_ax1.grid(False)
        f_ax1.set(yticklabels=[], xticklabels=[])
        f_ax1.tick_params(left=False, right=False)
        x = np.arange(0, np.shape(barcode)[0], 1)
        f_ax1.set_title(nameFig, x=0.5, y=1.2)
        f_ax1.imshow(barcode.reshape(1, -1), cmap='bwr', aspect='auto', interpolation='nearest')

        #Correlation Matrix
        f_ax2 = fig.add_subplot(gs[1, 0])
        im = f_ax2.imshow(corr, cmap='seismic', vmin=-np.amax(corr), vmax=np.amax(corr))
        f_ax2.grid(False)
        plt.colorbar(im, location='left')

        #Gene density
        f_ax3 = fig.add_subplot(gs[1, 1])
        f_ax3.axis('off')
        x = density
        y = np.flip(np.arange(0,np.shape(density)[0]))
        f_ax3.axhspan(np.shape(density)[0]-centro_start, np.shape(density)[0]-centro_end, facecolor='red', alpha=0.2)
        f_ax3.plot(x, y, color='blue', alpha=1, linewidth=0.1)

    if jupyter:
        plt.show()
    else:
        plt.savefig(nameFile)
    plt.clf()

#####
#Filtering tools

def filtering(Hicmat,factor=1.5, printer=True):
	"""
	Filter the matrix and return the matrix filtered and the list of indexs removed
	
	Entries : Hicmat -> list
		  factor -> float
		  printer -> boolean
	Return : Hicmatreduce -> list
		 segmenter1 -> list
	"""
	if printer:
		bar = tqdm(range(1), desc="Filtering ")
	else:
		bar = range(1)

	for i in bar:
		pass

	Filterextremum=True
	Hicmatreduce=Hicmat
	#first step : filter empty bin
	sumHicmat=Hicmat.sum(axis = 0)
	segmenter1=sumHicmat>0
	A=np.where(segmenter1)
	Hicmatreduce=Hicmatreduce[A[1],:]
	Hicmatreduce=Hicmatreduce[:,A[1]]
	if Filterextremum:
		#second step : filter lower bin
		sumHicmat=np.sum(Hicmatreduce,0)
		msum=np.mean(sumHicmat)
		mstd=np.std(sumHicmat)
		mini = msum-mstd*factor
		maxi = msum+mstd*factor
		#Make the bolean condition
		newcond=mini < sumHicmat
		newcond2=sumHicmat < maxi
		newcond=np.logical_and(newcond,newcond2)
		B=np.where(newcond)
		#Filter
		Hicmatreduce=Hicmatreduce[B[1],:]
		Hicmatreduce=Hicmatreduce[:,B[1]]
		segmenter1=A[1][B[1]] #Create the binsaved index
	return Hicmatreduce,segmenter1
	
def unfiltering(binsaved, contact_map, shape, printer=True):
	"""
	Unfilter the matrix and return the matrix unfiltered
	
	Entries : binsaved -> list
		  contact_map -> list
		  shape -> list
		  printer -> boolean
	Return : unfiltered_map -> list
	"""
	if printer:
		bar = tqdm(range(1), desc="Unfiltering ")
	else:
		bar = range(1)

	for i in bar:
		pass

	unfiltered_map = contact_map
	binunsaved = []
	for i in range(shape[0]):
		if ((i in binsaved) == False):
			binunsaved.append(i)
	for i in binunsaved:
		unfiltered_map = np.insert(unfiltered_map, i, 0, axis= 0)
		unfiltered_map = np.insert(unfiltered_map, i, 0, axis= 1)
	return unfiltered_map


#####
#Binning tools

def bin2d(Data,p,q, printer=True):    
	"""   
	Data = input matrix
	p,q rescaling factors
	Written for sparse 
	"""
	if printer:
		bar = tqdm(range(1), desc="Changing Resolution ")
	else:
		bar = range(1)

	for i in bar:
		pass

	n,m=np.shape(Data);
	s=(int(np.ceil(n/p)),int(np.ceil(m/q)))
	i,j,d = sparse.find(Data);
	i=np.int_(np.ceil(i/p))
	j=np.int_(np.ceil(j/q))
	M=sparse.csr_matrix((d,(i,j)))
	return M

def bin1D(anumpyarray, resolutionfrom, resolutionto):
	"""
	in : A numpy array , number of bin in raw and in col
	out : the matrix binned
	"""
	print(resolutionto,resolutionfrom)
	if resolutionto>resolutionfrom:
		convertionfactor=np.ceil(resolutionto/resolutionfrom)
		s=anumpyarray.shape
		print("dimension du vecteur:",s)
		#has to be identical as result in other function like chrsizedict)
		newsizei=np.ceil(s[0]*resolutionfrom/resolutionto)
		newarray=np.zeros(int(newsizei))
		print("taille du vecteur appres rescale :",newarray.shape)
		i=0
		while i<newsizei:
			ifrom=int(i*convertionfactor)
			ito=int((i+1)*convertionfactor)
			if i==newsizei-1:
				asum=np.sum(anumpyarray[ifrom:])
			else:
				asum=np.sum(anumpyarray[ifrom:ito])
			newarray[i]=asum
			i+=1
		return newarray
	elif resolutionto==resolutionfrom:
		print("no binning")
		return anumpyarray
	else:
		print("wrong resolution parameter in bin1D")

def bin2dfullmat(anumpyarray, resolutionfrom, resolutionto):
	"""
	in : A numpy array , number of bin in raw and in col
	out : the matrix binned
	Written for full
	"""
	print('change of resolution from ',resolutionfrom,' to ',resolutionto)
	if resolutionto>resolutionfrom:
		convertionfactor=np.ceil(resolutionto/resolutionfrom)
		s=anumpyarray.shape
		print("Initial HiC size before binning:",s)
		#has to be identical as result in other function like chrsizedict)
		newsizei=np.ceil(s[0]*resolutionfrom/resolutionto)
		newsizej=np.ceil(s[1]*resolutionfrom/resolutionto)
		newarray=np.zeros((int(newsizei),int(newsizej)))
		print("HiC size after binning :",newarray.shape)
		i=0
		j=0
		while i<newsizei:
			while j<newsizej:
				ifrom=int(i*convertionfactor)
				ito=int((i+1)*convertionfactor)
				jfrom=int(j*convertionfactor)
				jto=int((j+1)*convertionfactor)
				if i==newsizei-1:
					asum=np.sum(anumpyarray[ifrom:,jfrom:jto])
				elif j==newsizej-1:
					asum=np.sum(anumpyarray[ifrom:ito,jfrom:])
				elif i==newsizei-1 and j==newsizej-1:
					asum=np.sum(anumpyarray[ifrom:,jfrom:])
				else:
					asum=np.sum(anumpyarray[ifrom:ito,jfrom:jto])
				newarray[i,j]=asum
				newarray[j,i]=asum
				j+=1
			i+=1
			j=0
		return newarray
	elif resolutionto==resolutionfrom:
		print("No binning")
		return anumpyarray
	else:
		print("Wrong resolution parameter")


#####
#Operations tools

def SVD(D):
	"""
	Calculate the SVD and return the first eigenvector
	
	out : SVD(D)
	""" 
	eigens_values,eigens_vectors = np.linalg.eig(D)
	return eigens_vectors[:,0]

def Corr(D, printer=True):
	"""
	Calculate the pearson correlation and return the matrix
	
	out : Corr(D)
	""" 
	lines, columns = np.shape(D)
	C = np.zeros((lines, columns))

	if printer:
		bar = tqdm(range(lines), desc="Pearson Correlation ")
	else:
		bar = range(lines)

	for i in bar:
		for j in range(columns):
			C[i, j] = stats.pearsonr(D[i, :], D[:, j])[0]
	return C
	
def OE(D, printer=True):
	"""
	Calculate the OE matrix and return it
	out : OE(D)
	"""  
	i=0
	j=0
	L=len(D)

	if printer:
		bar = tqdm(range(1), desc="O/E ")
	else:
		bar = range(1)
	
	for i in bar:
		pass
		
	while j<L:
		thediag=np.diag(D,k=j)
		mtg=np.mean(thediag)
		if mtg == 0:
			mtg = 1e-8
		while i<(L-j):
			v=D[i,i+j]/mtg
			D[i,i+j]=v
			D[i+j,i]=v
			i+=1
		i=0
		j+=1
	return D

def SCN(D, max_iter = 10, printer=True):
	"""
	Calculate the SCN and return it
	Out  : SCN(D)
	Code version from Vincent Matthys
	"""    
	# Iteration over max_iter

	if printer:
		bar = tqdm(range(max_iter), desc="SCN ")
	else:
		bar = range(max_iter)

	for i in bar:        
		D /= np.maximum(1, D.sum(axis = 0))       
		D /= np.maximum(1, D.sum(axis = 1))    
		# To make matrix symetric again   
	return (D + D.T)/2 

def fastFloyd(contact, printer=True):
	"""
	out : FF(contact)
	Code version from Vincent Matthys
	"""      
	n = contact.shape[0]    
	shortest = contact    

	if printer:
		bar = tqdm(range(n), desc="Fast Floyd ")
	else:
		bar = range(n)

	for k in bar:        
		i2k = np.tile(shortest[k,:], (n, 1))        
		k2j = np.tile(shortest[:, k], (n, 1)).T        
		shortest = np.minimum(shortest, i2k + k2j)    
	return shortest


def filteramat(Hicmat, Filterextremum=True, factor=1.5, printer=True):
	"""
	in : a HiCmat without any transformation, factor of reduction
	out : the HiCmatreduce,thevector of his transformation
	THE filter part from the main in one function
	"""
	if printer:
		bar = tqdm(range(1), desc="Matrix Filtering ")
	else:
		bar = range(1)

	for i in bar:
		pass
	Hicmatreduce=Hicmat
	#first step : filter empty bin
	sumHicmat=Hicmat.sum(axis = 0)
	segmenter1=sumHicmat>0
	A=np.where(segmenter1)
	Hicmatreduce=Hicmatreduce[A[1],:]
	Hicmatreduce=Hicmatreduce[:,A[1]]
	if Filterextremum:
		#second step : filter lower bin
		sumHicmat=np.sum(Hicmatreduce,0)
		msum=np.mean(sumHicmat)
		mstd=np.std(sumHicmat)
		mini = msum-mstd*factor
		maxi = msum+mstd*factor
		#Make the bolean condition
		newcond=mini < sumHicmat
		newcond2=sumHicmat < maxi
		newcond=np.logical_and(newcond,newcond2)
		B=np.where(newcond)
		#Filter
		Hicmatreduce=Hicmatreduce[B[1],:]
		Hicmatreduce=Hicmatreduce[:,B[1]]
		segmenter1=A[1][B[1]] #Create the binsaved index
	return Hicmatreduce,segmenter1


def sammon(x, n, display = 2, inputdist = 'raw', maxhalves = 20, maxiter = 500, tolfun = 1e-9, init = 'default'):


    """Perform Sammon mapping on dataset x

    y = sammon(x) applies the Sammon nonlinear mapping procedure on
    multivariate data x, where each row represents a pattern and each column
    represents a feature.  On completion, y contains the corresponding
    co-ordinates of each point on the map.  By default, a two-dimensional
    map is created.  Note if x contains any duplicated rows, SAMMON will
    fail (ungracefully). 

    [y,E] = sammon(x) also returns the value of the cost function in E (i.e.
    the stress of the mapping).

    An N-dimensional output map is generated by y = sammon(x,n) .

    A set of optimisation options can be specified using optional
    arguments, y = sammon(x,n,[OPTS]):

       maxiter        - maximum number of iterations
       tolfun         - relative tolerance on objective function
       maxhalves      - maximum number of step halvings
       input          - {'raw','distance'} if set to 'distance', X is 
                        interpreted as a matrix of pairwise distances.
       display        - 0 to 2. 0 least verbose, 2 max verbose.
       init           - {'pca', 'cmdscale', random', 'default'}
                        default is 'pca' if input is 'raw', 
                        'msdcale' if input is 'distance'

    The default options are retrieved by calling sammon(x) with no
    parameters.

    File        : sammon.py
    Date        : 18 April 2014
    Authors     : Tom J. Pollard (tom.pollard.11@ucl.ac.uk)
                : Ported from MATLAB implementation by 
                  Gavin C. Cawley and Nicola L. C. Talbot

    Description : Simple python implementation of Sammon's non-linear
                  mapping algorithm [1].

    References  : [1] Sammon, John W. Jr., "A Nonlinear Mapping for Data
                  Structure Analysis", IEEE Transactions on Computers,
                  vol. C-18, no. 5, pp 401-409, May 1969.

    Copyright   : (c) Dr Gavin C. Cawley, November 2007.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

    """

    # Create distance matrix unless given by parameters
    if inputdist == 'distance':
        D = x
        if init == 'default':
            init = 'cmdscale'
    else:
        D = distance.cdist(x, x)
        if init == 'default':
            init = 'pca'

    if inputdist == 'distance' and init == 'pca':
        raise ValueError("Cannot use init == 'pca' when inputdist == 'distance'")

    if np.count_nonzero(np.diagonal(D)) > 0:
        raise ValueError("The diagonal of the dissimilarity matrix must be zero")

    # Remaining initialisation
    N = x.shape[0]
    scale = 0.5 / D.sum()
    D = D + np.eye(N)     

    if np.count_nonzero(D<=0) > 0:
        raise ValueError("Off-diagonal dissimilarities must be strictly positive")   

    Dinv = 1 / D
    if init == 'pca':
        [UU,DD,_] = np.linalg.svd(x)
        y = UU[:,:n]*DD[:n] 
    elif init == 'cmdscale':
        from cmdscale import cmdscale
        y,e = cmdscale(D)
        y = y[:,:n]
    else:
        y = np.random.normal(0.0,1.0,[N,n])
    one = np.ones([N,n])
    d = distance.cdist(y,y) + np.eye(N)
    dinv = 1. / d
    delta = D-d 
    E = ((delta**2)*Dinv).sum()

    # Get on with it
    for i in range(maxiter):

        # Compute gradient, Hessian and search direction (note it is actually
        # 1/4 of the gradient and Hessian, but the step size is just the ratio
        # of the gradient and the diagonal of the Hessian so it doesn't
        # matter).
        delta = dinv - Dinv
        deltaone = np.dot(delta,one)
        g = np.dot(delta,y) - (y * deltaone)
        dinv3 = dinv ** 3
        y2 = y ** 2
        H = np.dot(dinv3,y2) - deltaone - np.dot(2,y) * np.dot(dinv3,y) + y2 * np.dot(dinv3,one)
        s = -g.flatten(order='F') / np.abs(H.flatten(order='F'))
        y_old    = y

        # Use step-halving procedure to ensure progress is made
        for j in range(maxhalves):
            s_reshape = np.reshape(s, (-1,n),order='F')
            y = y_old + s_reshape
            d = distance.cdist(y, y) + np.eye(N)
            dinv = 1 / d
            delta = D - d
            E_new = ((delta**2)*Dinv).sum()
            if E_new < E:
                break
            else:
                s = 0.5*s

        # Bomb out if too many halving steps are required
        if j == maxhalves-1:
            print('Warning: maxhalves exceeded. Sammon mapping may not converge...')

        # Evaluate termination criterion
        if abs((E - E_new) / E) < tolfun:
            if display:
                print('TolFun exceeded: Optimisation terminated')
            break

        # Report progress
        E = E_new
        if display > 1 and printer:
            print('epoch = %d : E = %12.10f'% (i+1, E * scale))

    if i == maxiter-1:
        print('Warning: maxiter exceeded. Sammon mapping may not have converged...')

    # Fiddle stress to match the original Sammon paper
    E = E * scale
    
    return [y,E]

#####
#Predicting sub-compartments tools

def gaussianHMM(vector, nb_comp=2, n_iter=100):
	"""
	generate a simple hmm to have compartment on correlation map
	Return the model of the hmm and the states of the prediction at given N
	
	Entries : vector -> list
		  nb_comp -> int
		  n_iter -> int
	Return : labels_list -> list
		 scores_list -> list
	"""
	# Run Gaussian HMM
	model = hmm.GaussianHMM(nb_comp, "tied",n_iter=n_iter)
	model.fit(vector)
	labels = model.predict(vector) #etats/compartiments
	score = model.score(vector)

	return labels, score

def multiplegaussianHMM(vector, nb_comp_max=16, n_iter=100):
	"""
	generate the hmms of each compartment on correlation map
	Return the model of all hmm and the states of the prediction at given N
	
	Entries : vector -> list
		  nb_comp_max -> int
		  n_iter -> int
	Return : labels_list -> list
		 scores_list -> list
	"""
	labels_list = []
	scores_list = []
	for i in tqdm (range(2, nb_comp_max), desc="Finding Compartments with HMM "):
		labels, score = gaussianHMM(vector, i, n_iter=n_iter)
		labels_list.append(labels)
		scores_list.append(score)

	return labels_list, scores_list

def expr_repr_scoring(color_bins, marks, scores):
	"""
	Calculate the expression repression score and return it
	
	Entries : color_bins -> list
		  marks -> list
		  scores -> list
	Return : float
	"""
	output = color_bins[:,0].toarray() #only zeros
	for i in range(len(marks)):
		mark = marks[i]
		score = scores[i]
		X = color_bins[:,mark].toarray()
		output+=X*score

	return output

def similarity_score(val_data, data):
	"""
	Calculate the similarity score and return it
	
	Entries : val_data -> list
		  data -> list
	Return : float
	"""
	cnt = 0
	for i in range(np.shape(data)[0]):
		if val_data[i] == data[i]:
			cnt += 1

	similarity = cnt/np.shape(data)[0] * 100
	
	return round(similarity, 3)

#####
#PDB maker tool


def writePDB(fnameout, value, EpiValue, printer=True):
	"""
	Write a PDB file from values
	
	Entries : fnameout -> String
		  value -> int
		  EpiValue -> list  
	"""
	if printer:
		bar = tqdm(range(1), desc="Writing PDB ")
	else:
		bar = range(1)

	for i in bar:
		pass
	Sh=np.shape(value)
	#out
	fout=open(fnameout,'w')
	i=1
	while i<=Sh[0]:
		S="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}" #PDB format in python
		#print("here",value[i-1,0],value[i-1,1],EpiValue[i-1])
		S=S.format('ATOM',i,'CA','','ALA','A',i,'',float(value[i-1,0]),float(value[i-1,1]),float(value[i-1,2]),1,float(EpiValue[i-1]),'C','')
		#print(S)
		fout.write(S+"\n")
		i+=1
	fout.close()

