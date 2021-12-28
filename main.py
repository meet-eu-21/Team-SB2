#python 3
#2020
#CC-By-SA
#Old code from Carron Leopold, Julien Mozziconacci
#Modified by Damien Legros, Cédric Cornede, Arnaud Quelin, Rouquaya Mouss, Hamid Hachemi
#####
#Imports

import h5py
import sys
import numpy as np
from scipy.spatial import ConvexHull
from scipy import sparse
import pandas as pd
from tqdm import tqdm
from hmmlearn.hmm import GaussianHMM
from sklearn.manifold import MDS #If you want the scikit learn mds
import HiCtoolbox
from All_possi import *
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def main(HiCfilename, nb_comp_max, chromosome_name="", printer=True):
	#####
	#Variables
	R=100000 #Resolution
	NbmaxEpi=15 #Epi states go from 0 to 15
	alpha=0.227 #Alpha
	selectedmark=1 #Index of the selected mark
	HiCfilename=HiCfilename #HiC File

	chromosome_number="chr"+HiCfilename.split("chr")[1].split("_")[0]
	EpiGfilename='E116_15_coreMarks_dense' #Epigenic Marks File
	Densityfilename = chromosome_number+'.hdf5'
	pathToDataLeopold = "../HiCdataLeopold/GM12878/100kb/"+chromosome_number+"_compartiment.txt"

	#####
	#Plot results Leopold

	chr2 = pd.read_table(pathToDataLeopold)
	chr2F = chr2[chr2.iloc[:,0]!= -1] # -1 = régions filtrées
	chr2F=chr2F.to_numpy()

	pixel_per_bar = 0.5
	dpi = 100

	fig = plt.figure(figsize=(len(chr2F) * pixel_per_bar / dpi, 2), dpi=dpi)
	ax = fig.add_axes([0, 0, 1, 1]) 
	ax.imshow(chr2F.reshape(1, -1), cmap='bwr', aspect='auto',interpolation='nearest')
	ax.axes.get_yaxis().set_visible(False)
	plt.savefig("BarcodeLeopold_"+chromosome_name+'.pdf')
	plt.figure().clear()

	del chr2
	del chr2F
	del fig
	del ax

	#####
	#Plot Density

	data = h5py.File(Densityfilename, 'r')
	data2D = data['data'][:]
	density = data2D
	HiCtoolbox.plotter(data2D, 'Density_plot_'+chromosome_name+'.pdf', 'Density')
	
	del data

	#####
	#Build matrix
	
	print('#####')
	A=HiCtoolbox.buildMatrix(HiCfilename, printer)
	print('#####')
	print('Original resolution : ', np.shape(A))
	print('#####')
	#!become csr sparse array
	binned_map=HiCtoolbox.bin2d(A,R,R, printer)
	print('#####')
	LENTEST=np.shape(A)[0]
	goodShape = np.shape(binned_map)
	print('Input at the good resolution : ', goodShape)
	print('#####')
	
	#####
	#Build color annotation at desired resolution
	
	color_vec = HiCtoolbox.buildColors(EpiGfilename, chromosome_number, LENTEST, printer)
	print('#####')
	color_bins=HiCtoolbox.bin2d(color_vec,R,1, printer)
	print('#####')
	color_bins=color_bins/np.amax(color_bins)
	print('Bp cover by this mark, has to be >0 :',np.sum(color_bins[:,selectedmark]) )
	print('#####')
	
	#####
	#Filtering
	
	print("Before filtering : ",np.shape(binned_map))
	print('#####')
	data2D = np.log2(1 + sparse.csr_matrix.toarray(binned_map))
	HiCtoolbox.plotter(data2D, 'BeforeFiltering_plot_'+chromosome_name+'.pdf', 'Matrix', "hot_r")
	#Filtering
	filtered_map, binsaved = HiCtoolbox.filtering(binned_map, 1.5, printer)
	print('#####')
	print("After filtering : ", np.shape(filtered_map)) #,np.shape(color_vecseg))
	print('#####')
	data2D = np.log2(1 + sparse.csr_matrix.toarray(filtered_map))
	HiCtoolbox.plotter(data2D, 'AfterFiltering_plot_'+chromosome_name+'.pdf', 'Matrix', "hot_r")
	
	#####
	#Update color annotation at desired resolution
	
	#Filter the epi by removed bin in HiC
	color2=color_bins[binsaved[1]]
	#Now color2 is 1D
	color2=color2[:,selectedmark]
	#Type issue
	color2=np.float64(color2.todense())
	
	######
	#SCN

	scn_map = HiCtoolbox.SCN(filtered_map.copy(), printer)
	print('#####')
	#Now we are not sparse at all
	scn_map = np.asarray(scn_map)**alpha
	data2D = np.log2(scn_map)
	HiCtoolbox.plotter(data2D, 'SCN_plot_'+chromosome_name+'.pdf', 'Matrix' , "hot_r")
	
	#####
	#Distance Matrix
	
	#Shortest path on the matrix
	dist_matrix = HiCtoolbox.fastFloyd(1/scn_map, printer)
	print('#####')
	#Remove the diagonal
	dist_matrix = dist_matrix-np.diag(np.diag(dist_matrix))
	#Just to be sure that the matrix is symetric, not really usefull in theory
	dist_matrix = (dist_matrix+np.transpose(dist_matrix))/2
	data2D = np.log(dist_matrix)
	HiCtoolbox.plotter(data2D, 'DistMatrix_plot_'+chromosome_name+'.pdf', 'Matrix', "hot_r")
	
	#####
	#O/E
	
	contact_map = HiCtoolbox.OE(scn_map, printer)
	print('#####')
	data2D = np.log(contact_map)
	HiCtoolbox.plotter(data2D, 'OE_plot_'+chromosome_name+'.pdf', 'Matrix', "seismic", vmin=-np.amax(data2D), vmax=np.amax(data2D))
	
	#####
	#Correlation Matrix
	
	contact_map = HiCtoolbox.Corr(contact_map, printer)
	print('#####')
	data2D = contact_map
	HiCtoolbox.plotter(data2D, 'Corr_plot_'+chromosome_name+'.pdf', 'Matrix', "seismic", vmin=-np.amax(data2D), vmax=np.amax(data2D))
	
	#####
	#Unfiltering
	
	unfiltered_map = HiCtoolbox.unfiltering(binsaved, scn_map, goodShape, printer)
	print('#####')
	data2D = unfiltered_map
	HiCtoolbox.plotter(data2D, 'UnfilteredSCN_plot_'+chromosome_name+'.pdf', 'Matrix', "hot_r")
	
	unfiltered_map = HiCtoolbox.unfiltering(binsaved, contact_map, goodShape, printer)
	print('#####')
	data2D = unfiltered_map
	Corrr = data2D
	HiCtoolbox.plotter(data2D, 'UnfilteredCorr_plot_'+chromosome_name+'.pdf', 'Matrix', "seismic", vmin=-np.amax(data2D), vmax=np.amax(data2D))
	
	#####
	#AB Compartments
	
	data2D = HiCtoolbox.SVD(unfiltered_map)
	EigenV = data2D.reshape(-1, 1)
	HiCtoolbox.plotter(data2D, 'SVD_plot_'+chromosome_name+'.pdf', 'AB')
	
	#####
	#MDS
	
	#embedding = MDS(n_components=3)#LOAD the MDS #With scikit-learn mds
	#XYZ = embedding.fit_transform(dist_matrix) #Make the transform
	#XYZ = np.float64(XYZ)
	XYZ,E = HiCtoolbox.sammon(dist_matrix, 3, printer)#with the one from tom j pollard
	print("Output shape : ",np.shape(XYZ),np.shape(color2))
	print('#####')
	#Point rescale	
	hull=ConvexHull(XYZ)
	scale=100/hull.area**(1/3)
	XYZ=XYZ*scale

	#####
	#Gaussian HMM #USE THE HiCtoolbox VERSION

	nb_etats = np.arange(2,16,1)
	Scores_Contact = []
	Preds_Contact = []
	X = EigenV.reshape(-1, 1)
	for n in nb_etats:
		model = GaussianHMM(n, "diag",n_iter=100)
		model.fit(X)
		hidden_states = model.predict(X)#etats/compartiments
		Preds_Contact.append(hidden_states)
		score = model.score(X)
		Scores_Contact.append(score)

	best_nb = 2
	threshold = 1
	for i in range(1,len(Scores_Contact)):
		temp = ((Scores_Contact[i]-Scores_Contact[i-1])/Scores_Contact[i])*100
		print('Nombre de compartiments :', i+2,'\nVariation par rapport au score precedent : ', temp,'%')
		if temp >= threshold: #we have found a best number of compartment -> score increase significantly
			best_nb+=1
		else:
			print('The Score no longer increases enough')
			break
	print('\nBest compartment number : ', best_nb)

	#####
	#Plot results Barcode HMM #USE THE HiCtoolbox VERSION

	chr2 = pd.DataFrame(Preds_Contact[0])
	chr2F = chr2[chr2.iloc[:,0]!= -1] # -1 = régions filtrées
	chr2F=chr2F.to_numpy()

	pixel_per_bar = 0.5
	dpi = 100

	fig = plt.figure(figsize=(len(chr2F) * pixel_per_bar / dpi, 2), dpi=dpi)
	ax = fig.add_axes([0, 0, 1, 1]) 
	ax.imshow(chr2F.reshape(1, -1), cmap='bwr', aspect='auto',interpolation='nearest')
	ax.axes.get_yaxis().set_visible(False)
	plt.savefig("BarcodeHMM_"+chromosome_name+'.pdf')
	plt.figure().clear()

	del chr2
	del chr2F
	del fig
	del ax

	Pred_HMM_Contact = Preds_Contact[best_nb-2]

	#Writing PDB File
	HiCtoolbox.writePDB('HMMcolors_Contact_'+chromosome_name+'_'+str(best_nb)+'_compartments.pdb',XYZ,Pred_HMM_Contact)

	HiCtoolbox.plotter(Scores_Contact, 'HMMScore_plot_'+chromosome_name+'.pdf', 'HMMScore')

	#####
	#Gaussian HMM with Epigenetic marks #USE THE HiCtoolbox VERSION

	marks = np.arange(1,16)
	split = 21.43/15
	scores = np.arange(-10,11,split)
	EpiMark = score_EpiMark(color_bins, marks, scores)

	nb_etats = np.arange(2,16,1)
	Scores_Epi = []
	Preds_Epi = []
	X = EpiMark
	for n in nb_etats:
		model = GaussianHMM(n, "diag",n_iter=100)
		model.fit(X)
		hidden_states = model.predict(X)#etats/compartiments
		Preds_Epi.append(hidden_states)
		score = model.score(X)
		Scores_Epi.append(score)

	best_nb = 2
	threshold = 1
	for i in range(1,len(Scores_Epi)):
		temp = (np.abs((Scores_Epi[i]-Scores_Epi[i-1]))/np.abs(Scores_Epi[i]))*100
		print('Nombre de compartiments :', i+2,'\nVariation par rapport au score precedent : ', temp,'%')
		if temp >= threshold: #we have found a best number of compartment -> score increase significantly
			best_nb+=1
		else:
			print('The Score no longer increases enough')
			break
	print('\nBest compartment number : ', best_nb)

	#####
	#Plot results Barcode HMM Epi #USE THE HiCtoolbox VERSION

	chr2 = pd.DataFrame(Preds_Epi[0])
	chr2F = chr2[chr2.iloc[:,0]!= -1] # -1 = régions filtrées
	chr2F=chr2F.to_numpy()

	pixel_per_bar = 0.5
	dpi = 100

	fig = plt.figure(figsize=(len(chr2F) * pixel_per_bar / dpi, 2), dpi=dpi)
	ax = fig.add_axes([0, 0, 1, 1]) 
	ax.imshow(chr2F.reshape(1, -1), cmap='bwr', aspect='auto',interpolation='nearest')
	ax.axes.get_yaxis().set_visible(False)
	plt.savefig("BarcodeEpi_"+chromosome_name+'.pdf')
	plt.figure().clear()

	del chr2
	del chr2F
	del fig
	del ax

	Pred_HMM_Epi = Preds_Epi[best_nb-2]

	#Writing PDB File
	HiCtoolbox.writePDB('HMMcolors_EpiMarks_'+chromosome_name+'_'+str(best_nb)+'_compartments.pdb',XYZ,Pred_HMM_Epi)

	HiCtoolbox.plotter(Scores_Epi, 'HMMEpiScore_plot_'+chromosome_name+'.pdf', 'HMMScore')

	
	#####
	#Writing PDB File
	
	HiCtoolbox.writePDB('3Dcolors_'+str(alpha)+'_'+chromosome_name+'.pdb', XYZ, color2, printer)
	print('#####')
	print("Done :D")
	print('#####')


	#####
	#Similarity tools #TO PUT IN HiCtoolbox AND MAKE IT WORK

	#/!\ if nb_methods=1 -> we compare 2 methods
	#nb_scenarios = (nb_cpt^(nb_cpt+nb_methods))^2
	best_scenarios = []
	similarities = []
	for n in range(2,16):
		nb_methods = 1
		nb_cpt = n
		All_possi = get_All_possiblesPaths(nb_cpt, nb_methods)

		All_Preds = [Preds_Contact[nb_cpt-2],Preds_Epi[nb_cpt-2]]
		All_tests = []
		for scenario in All_possi:
			similarity = calculateSimilarity(scenario, All_Preds, nb_cpt, nb_methods)
			All_tests.append([scenario, similarity])

		best_scenario, similarity = find_bestScenario(All_tests, nb_cpt)
		best_scenarios.append(best_scenario)
		similarities.append(similarity)

	print('Best Scenario : ', best_scenarios)
	#scenario : ([labelcptPred1, labelcptPred2, ...], ...)
	print('Similarity : ', similarities)
	print('#####')

	HiCtoolbox.plotter(similarities, 'Similarity_Score_plot_'+chromosome_name+'.pdf', 'HMMScore')

#####
#Main Function

if __name__ == '__main__': #TO MODIFY -> USE ArgumentParser
	
	# Initialize the Parser
	parser = argparse.ArgumentParser(description ='Compartments Script : ')

	
	# Adding Argument
	files = ['inter_100kb', 'intra_100kb', 'intra_25kb', '100kb', 'intra', 'all']
	#print(sys.argv[1][-12:-1])
	#Usage for each file of HiC Folder
	if sys.argv[1] in files:
		if len(sys.argv)==3:
			nb_comp_max = sys.argv[2]
		elif len(sys.argv)==2:
			nb_comp_max = 15
		pathToChrList = HiCtoolbox.findfiles(sys.argv[1])
		print(pathToChrList)
		for pathToChr in pathToChrList:
			chrName = pathToChr.replace("/","_")[:-12]
			print(chrName)
			for i in tqdm (range(1), desc="Analazing "+str(chrName)+" "):
				main(pathToChr, nb_comp_max, chromosome_name=chrName, printer=False)

	#Usage for one file
	elif sys.argv[1][-12:-1] == ".RAWobserve":
		if len(sys.argv)==3:
			nb_comp_max = sys.argv[2]
		elif len(sys.argv)==2:
			nb_comp_max = 15
		main(sys.argv[1], nb_comp_max, printer=True)

	#Help guide
	else:
		print("====Help Guide====\n")
		print("USAGE FOR ONE FILE : python ", sys.argv[0], "<HiC filename> <OPTIONAL: nb compartments max (default : 15)>\n")
		print("USAGE FOR EACH FILE OF HiCdata FOLDER : python ", sys.argv[0], "<files to work with ", str(files), "> <OPTIONAL: nb compartments max (default : 15)>")
		print("\n==================")
		sys.exit(1)

