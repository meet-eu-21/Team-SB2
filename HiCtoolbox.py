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
from sklearn.manifold import TSNE
from sklearn.decomposition import KernelPCA
import umap
import seaborn as sns

np.seterr(divide = 'ignore') 

#####
#ADDITIONAL LOADER

def EpiGbyres(EpiGfilename,res,achr,sizeatres,NbEpi):
	"""
	Generate a matrix of repartition of epiG in a chr at the given resolution
	Do it the hard way : parse the file to the desired resolution directly
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
	in :
	out : 
	"""
	if printer:
		bar = tqdm(range(1), desc="Loading Matrix ")
	else:
		bar = range(1)

	for i in bar:
		pass
	A=np.loadtxt(HiCfilename)
	A=np.int_(A)
	#Build array at pb resolution
	A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)
	A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
	return A

def buildColors(EpiGfilename, chromosome_number, LENTEST, printer=True):
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

def plotter(data2D, nameFig="Plot", plotType="Matrix", cmap="hot_r", vmin=None, vmax=None, jupyter=False, nameFile="Plot.pdf", centro_start=None, centro_end=None):
	"""
	in :
	out : 
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
		custom_params = {"axes.spines.right": False, 
						 "axes.spines.top": False}
		sns.set_theme(style="ticks", rc=custom_params)
		x = np.arange(2,16)
		y = data2D
		fig = plt.figure()
		fig.set_dpi(300)
		plt.plot(x, y, color='blue', alpha=0.4)
		threshold = 0.5
		best = 2
		for i in range(1,len(data2D)):
			temp = np.abs(((data2D[i]-data2D[i-1])/data2D[i])*100)
			if temp >= threshold:
				best += 1
			else: 
				break
		plt.scatter(best, data2D[best-2], alpha=1, color='red', label="predicted number of compartments")
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

def filtering(binned_map, filter_ratio, printer=True):
	"""
	in :
	out : 
	"""
	if printer:
		bar = tqdm(range(1), desc="Filtering ")
	else:
		bar = range(1)

	for i in bar:
		pass

	sumHicmat=np.sum(binned_map,0) 
	mini = np.mean(sumHicmat)-np.std(sumHicmat)*filter_ratio #min value of filtering
	maxi = np.mean(sumHicmat)+np.std(sumHicmat)*filter_ratio #max value of filtering
	binsaved=np.where(np.logical_and(mini < sumHicmat,sumHicmat < maxi)) #coord of bin to save
	filtered_map=binned_map[binsaved[1],:] #save on raw
	filtered_map=filtered_map[:,binsaved[1]] #save on col
	return filtered_map, binsaved
	
def unfiltering(binsaved, contact_map, shape, printer=True):
	"""
	in :
	out : 
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
		if ((i in binsaved[1]) == False):
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
	out : SVD(D)
	Code version from Damien LEGROS
	""" 
	eigens_values,eigens_vectors = np.linalg.eig(D)
	return eigens_vectors[:,0]

def Corr(D, printer=True):
	"""
	out : Corr(D)
	Code version from Damien LEGROS
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
	out : OE(D)
	Code version from Damien LEGROS
	"""  
	lines, columns = np.shape(D)
	means = np.zeros(lines*2-1)
	
	for diag in range(-lines+1, lines):
		means[diag] = np.mean(np.diagonal(D, diag))

	if printer:
		bar = tqdm(range(lines), desc="O/E ")
	else:
		bar = range(lines)
	
	for i in bar:
		for j in range(columns):
			D[i, j] = D[i, j] / means[j-i]
			
	return D

def SCN(D, max_iter = 10, printer=True):
	"""
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
	generate a simple hmm to have compartiment on correlation map
	out : the model of the hmm, AND the states of the prediction at given N
	"""
	# Run Gaussian HMM
	model = hmm.GaussianHMM(nb_comp, "tied",n_iter=n_iter)
	model.fit(vector)
	labels = model.predict(vector) #etats/compartiments
	score = model.score(vector)

	return labels, score

def multiplegaussianHMM(vector, nb_comp_max=16, n_iter=100):
	"""
	generate a simple hmm to have compartiment on correlation map
	out : the model of all hmm, AND the states of the prediction at given N
	"""
	labels_list = []
	scores_list = []
	for i in tqdm (range(2, nb_comp_max), desc="Finding Compartments with HMM "):
		labels, score = gaussianHMM(vector, i, n_iter=n_iter)
		labels_list.append(labels)
		scores_list.append(score)

	return labels_list, scores_list

def autoencoder(): #TODO
	"""
	in:
	out:
	"""
	pass

def expr_repr_scoring(color_bins, marks, scores):
	"""
	in:
	out:
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
	in:
	out:
	"""
	cnt = 0
	for i in range(np.shape(data)[0]):
		if val_data[i] == data[i]:
			cnt += 1

	similarity = cnt/np.shape(data)[0] * 100
	
	return round(similarity, 3)
		
#####
#Visualisation tools 

def tsne(chrMatrix, labels): #NOT USED
	"""
	generate data for tsne visualisation from labels
	"""
	data_map = chrMatrix
	states = labels

	transformer = KernelPCA(n_components=30, kernel='precomputed')
	transformed_data = transformer.fit_transform(data_map)

	X_embedded = TSNE(n_components=2, learning_rate=200, init='pca', perplexity=1000, n_iter=500).fit_transform(transformed_data)

	print(X_embedded[:,0])

	df = pd.DataFrame({'Score1':list(X_embedded[:,0]), 'Score2':list(X_embedded[:,1]), 'Color':list(states)})

	return df

def k_means(chrMatrix, labels): #NOT USED
	"""
	generate data for k_means visualisation from labels
	"""
	data_map = chrMatrix
	states = labels

	transformer = KernelPCA(n_components=2, kernel='precomputed')
	X_embedded = transformer.fit_transform(data_map)

	c, label, i = k_means(X_embedded, n_clusters=10)

	df = pd.DataFrame({'Score1':list(X_embedded[:,0]), 'Score2':list(X_embedded[:,1]), 'Color':list(label)})

	return df

def umap(chrMatrix, labels): #NOT USED
	"""
	generate data for umap visualisation from labels
	"""
	data_map = chrMatrix
	states = labels

	reducer = umap.UMAP(n_neighbors=15,
        				min_dist=0.1,
        				n_components=2,
        				metric='euclidean')

	X_embedded = reducer.fit_transform(data_map)

	df = pd.DataFrame({'Score1':list(X_embedded[:,0]), 'Score2':list(X_embedded[:,1]), 'Color':list(states)})

	return df

#####
#Similarity tree compartments tool

def calculateSimilarity(scenario, All_Preds, nb_cpt, nb_methods): #TO TEST
    dico_count = {}
    for cpt in range(len(All_Preds[0])):
        current_scenario = []
        for i in range(len(All_Preds)):
            current_scenario.append(All_Preds[i][cpt])
        if tuple(current_scenario) in dico_count:
            dico_count[tuple(current_scenario)]+=1
        else:
            dico_count[tuple(current_scenario)] = 0
    
    good_cla = 0
    bad_cla = 0
    for element in scenario:
        good_cla+=dico_count[tuple(element)]
    
    total = np.sum(list(dico_count.values()))
    bad_cla = total-good_cla
    similarity = good_cla/total
    return similarity

def find_bestScenario(All_tests, nb_cpt): #TO TEST
    best = All_tests[0]
    similarity = 0
    for scenario in All_tests:
        good_association = True
        current_similarity = scenario[1]
        for i in range(len(scenario[0][0])):
            test = set(np.array(scenario[0])[:,i])
            if (len(test)<nb_cpt):
                good_association = False
        if (good_association==True) and (current_similarity>similarity):
            best = scenario
            similarity = current_similarity

    return best

#####
#Files finder tools

def findfiles(filesToFind): #NOT USED 
	"""
	return paths to chr
	TODO : possibility to make it shorter and cleaner
	"""
	pathToData = "../HiCdata"
	pathToResults ="../HiCresults"
	cellNameList = os.listdir(pathToData)
	pathToChrList = []
	pathToPlotList = []
	os.makedirs(pathToResults, exist_ok = True)

	for cellName in tqdm(cellNameList, desc="Getting paths to work with "):
		os.makedirs(pathToResults+"/"+cellName, exist_ok = True)

		#"intra_100kb": work with every interchromosomic 100kb files
		if filesToFind == 'inter_100kb':
			typeNresList = os.listdir(pathToData+"/"+cellName)
			typeNresList = [element for element in typeNresList 
							if (filesToFind.split("_")[0] and 
								filesToFind.split("_")[1]) in element]
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[0], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)

		#"intra_100kb": work with every intrachromosomic 100kb files
		elif filesToFind == 'intra_100kb':
			typeNresList = os.listdir(pathToData+"/"+cellName)
			typeNresList = [element for element in typeNresList 
							if (filesToFind.split("_")[0] and 
								filesToFind.split("_")[1]) in element]
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[0], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)

		#"intra_25kb": work with every intrachromosomic 25kb files
		elif filesToFind == 'intra_25kb':
			typeNresList = os.listdir(pathToData+"/"+cellName)
			typeNresList = [element for element in typeNresList 
							if (filesToFind.split("_")[0] and 
								filesToFind.split("_")[1]) in element]
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[0], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)

		#"100kb": work with every 100kb files
		elif filesToFind == '100kb':
			typeNresList = os.listdir(pathToData+"/"+cellName)
			typeNresList = [element for element in typeNresList 
							if filesToFind in element]
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[0], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[1], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[1]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[1]+"/"+chromosome)

		#"25kb": work with every 25kb files
		elif filesToFind == '25kb':
			typeNresList = os.listdir(pathToData+"/"+cellName)
			typeNresList = [element for element in typeNresList 
							if filesToFind in element]
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[0], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[1], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[1]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[1]+"/"+chromosome)
		
		#"all": work with every file
		elif filesToFind == 'all':
			typeNresList = os.listdir(pathToData+"/"+cellName)
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[0], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[0]+"/"+chromosome)
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[1], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[1]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[1]+"/"+chromosome)
			os.makedirs(pathToResults+"/"+cellName+"/"+typeNresList[2], exist_ok = True)
			chrList = os.listdir(pathToData+"/"+cellName+"/"+typeNresList[0])
			for chromosome in chrList:
				pathToChrList.append(pathToData+"/"+cellName+
									 "/"+typeNresList[2]+"/"+chromosome)
				pathToPlotList.append(pathToResults+"/"+cellName+
									 "/"+typeNresList[2]+"/"+chromosome)

	return pathToChrList


#####
#PDB maker tool


def writePDB(fnameout, value, EpiValue, printer=True):
	"""
	write a PDB from value contain in value
	just bored to read 200 pages of biopython for a simple writing script
	I use tips from pierre poulain blog
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

