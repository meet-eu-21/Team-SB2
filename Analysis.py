#Imports used
import h5py
import sys
import os
import numpy as np
from scipy.spatial import ConvexHull
from scipy import sparse
from sklearn.manifold import MDS
from hmmlearn import hmm

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import random as rd
import seaborn as sns

import HiCtoolbox #functions in our file HiCtoolbox.py

from All_possi import * #Functions to create association between labels of different predictions

import urllib.request #package to download data

#Function to create folders
def create_Folders(cell_line, chromo, resol, FOLDER):
    """
    Function to make the folders
    
    Enters : cell_line -> []String
    	     chromo > []String
    	     resol -> []String
    	     FOLDER -> String
    Returns : folders_exist -> Boolean
    
    Create the folders for the data (if FOLDER=Data) or for the results (if FOLDER=Results)
    make the subfolders for the cell line, then the chromosomes in the end the resolution
    """
    chromo = str(chromo)
    if chromo=='23':
        chromo = 'X'
    subfolders = [ f.path for f in os.scandir(FOLDER) if f.is_dir() ]
    flag_cell = False #to check if cell line folder exist
    flag_chr = False #to check if chromosome folder exist
    flag_resol = False #to check if resolution folder exist
    folders_exist = False # to check if all folders already exist
    #Check which folder exist
    for folder in subfolders:
        cell_folder = folder[5:]
        if cell_line in cell_folder:
            flag_cell = True
    if flag_cell:
        subfolders = [ f.path for f in os.scandir(FOLDER+cell_line+'/') if f.is_dir() ]
        for folder in subfolders:
            chromo_folder = folder[-2:]
            for chromo_fold in chromo_folder.split('\n'):
                try:
                    test = int(chromo_fold)
                except:
                    chromo_fold = chromo_fold[1]
                if chromo == chromo_fold:
                    flag_chr = True
    if flag_chr:
        subfolders = [ f.path for f in os.scandir(FOLDER+cell_line+'/'+chromo+'/') if f.is_dir() ]
        for folder in subfolders:
            resol_folder = folder[-3:]
            if resol in resol_folder:
                flag_resol = True
           
    #Create necessary folders
    if (flag_cell==False): #Create all folders
        #Cell folder
        cell_path = os.path.join(FOLDER, cell_line)
        os.mkdir(cell_path)
        #Chromosome folder
        parent_dir = FOLDER+cell_line+'/'
        chr_path = os.path.join(parent_dir, chromo)
        os.mkdir(chr_path)
        #Resolution folder
        parent_dir = FOLDER+cell_line+'/'+chromo+'/'
        resol_path = os.path.join(parent_dir, resol)
        os.mkdir(resol_path)
        
    elif (flag_cell==True) and (flag_chr==False):
        #Chromosome folder
        parent_dir = FOLDER+cell_line+'/'
        chr_path = os.path.join(parent_dir, chromo)
        os.mkdir(chr_path)
        #Resolution folder
        parent_dir = FOLDER+cell_line+'/'+chromo+'/'
        resol_path = os.path.join(parent_dir, resol)
        os.mkdir(resol_path)
        
    elif (flag_chr==True) and (flag_resol==False):
        #Resolution folder
        parent_dir = FOLDER+cell_line+'/'+chromo+'/'
        resol_path = os.path.join(parent_dir, resol)
        os.mkdir(resol_path)
        
    if (flag_cell==True) and (flag_chr==True) and (flag_resol==True): #Analysis already done
        folders_exist = True
    return folders_exist

#Function to download data and put it in the corresponding folder
def Download_data(cell_line, chromo, resol):
    """
    Function to download the data
    
    Enters : cell_line -> []String
    	     chromo > []String
    	     resol -> []String
    Returns : folders_exist -> Boolean
    
    Download the chromosome data from the cell line given with the given resolution
    Download the following files from lcqb.upmc.fr :
    - the validation data of Leopold Carron
    - the gene density data
    - the HiC data in .RAWobserved
    """
    print("###################################Downloading Data Chr",chromo," Cell line ",cell_line, " Resolution ",resol,"###################################")
    resolution = resol+'kb'
    contact_type = 'intra'
    chromo = str(chromo)
    if chromo=='23':
        chromo = 'X'
    
    #download Leopold compartment file
    path = 'Data/'+cell_line+'/'+chromo+'/'+resol+'/chr'+chromo+'_compartiment.txt'
    if (os.path.isfile(path)==False): #check if file exist
        Path = 'http://www.lcqb.upmc.fr/meetu/dataforstudent/comp/'+cell_line+'/'+resolution+'/chr'+chromo+'_compartiment.txt'
        urllib.request.urlretrieve(Path,path)
    
    #download gene density file
    path = 'Data/'+'chr'+chromo+'.hdf5'    
    if (os.path.isfile(path)==False): #check if file exist
        Path = 'http://www.lcqb.upmc.fr/meetu/dataforstudent/genedensity/'+resolution+'/chr'+chromo+'.hdf5'
        urllib.request.urlretrieve(Path,path)

    #download HiC file
    if cell_line=='GM12878':
        Path = 'http://www.lcqb.upmc.fr/meetu/dataforstudent/HiC/'+cell_line+'/'+resolution+'_resolution_'+contact_type +'chromosomal/chr'+chromo+'_'+resolution+'.RAWobserved'
        Sequences = 'Data/'+cell_line+'/'+chromo+'/'+resol+'/chr'+chromo+'_'+resolution+'.RAWobserved'
        urllib.request.urlretrieve(Path,Sequences)
    else:
        Path = 'http://www.lcqb.upmc.fr/meetu/dataforstudent/HiC/'+cell_line+'/'+resolution+'_resolution_'+contact_type +'chromosomal/chr'+chromo+'/MAPQGE30/chr'+chromo+'_'+resolution+'.RAWobserved'
        Sequences = 'Data/'+cell_line+'/'+chromo+'/'+resol+'/chr'+chromo+'_'+resolution+'.RAWobserved'
        urllib.request.urlretrieve(Path,Sequences)


def Analysis(gene, nb, resol, filter_ratio=0.5, nb_max_epi=15, alpha=0.227, selected_mark=1):
    """
    Function to make the analysis
    
    Enters : gene -> String (cell line)
    	     nb -> int (chromosome number)
    	     filter_ratio -> float (filtering ratio)
    	     nb_max_epi -> int (number of marks epigenetics studied
    	     alpha -> float
    	     selected_mark -> int
    Returns : folders_exist -> Boolean
    
    Do the analysis : 
    - Build the matrix and filter the data
    - Build the SCN, distance and O/E matrix
    - Build the correlation matrix and unfilter the data
    - Make the A/B compartments and the MDS
    - Calculate the HMM from correlation matrix (First method)
    - Calculate the HMM form epigenetic marks (Second method)
    - Calculate the marks with WB1 code
    - Make the consenus labels and calculate the similarity
    """

    ###################################COMPUTING BASIC ANALYSIS###################################
    print("###################################Analysing Chr",nb," Cell line ",gene, " Resolution ",resol,"kb###################################")
    
    print("###################################COMPUTING BASIC ANALYSIS###################################")
     
    #####Global variables
    print("#####Global variables")
    
    nb = str(nb)
    if nb=='23':
        nb = 'X'
    res = int(resol+'000') #Resolution
    resolution = resol+'kb'
    
    result_path = 'Results/'+gene+'/'+nb+'/'+resol+'/'
    hic_filename = 'Data/'+gene+'/'+nb+'/'+resol+'/chr'+nb+'_'+resolution+'.RAWobserved' #HiC File
    chr_nb = 'chr'+nb
    epig_filename='Data/E116_15_coreMarks_dense' #Epigenic Marks File
    density_filename = 'Data/'+'chr'+nb+'.hdf5'
    path_to_leo_data = 'Data/'+gene+'/'+nb+'/'+resol+'/chr'+nb+'_compartiment.txt'

    centro_pos = pd.read_csv("Data/centromericpos_hg19.txt", sep="\t", header=None)
    centro_start = int(centro_pos.loc[centro_pos[0] == chr_nb][1])/100000
    centro_end = int(centro_pos.loc[centro_pos[0] == chr_nb][2])/100000

    expr_repr_scores_file = pd.read_csv("Data/expr_repr_scores.csv", sep="\t", header=None)
    expr_scores = expr_repr_scores_file[2].to_numpy()
    repr_scores = expr_repr_scores_file[4].to_numpy()
    expr_repr_scores = expr_scores - repr_scores

    validation_data = pd.read_table(path_to_leo_data)
    data = validation_data.to_numpy()
    data[data == 0.] = 2.
    data[data == -1.] = 0.
    data[data == 2.] = -1.
    val_data = data
    
    #####Density of the chromosome
    print("#####Density of the chromosome")
    density_data = h5py.File(density_filename, 'r')
    data = density_data['data'][:]
    density = data
    title = 'Density of chromosome '+nb+' from '+gene
    path = result_path+'Density'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Density', nameFile=path, centro_start=centro_start, centro_end=centro_end)

    #####Building the matrix
    print("#####Building the matrix")
    matrix = HiCtoolbox.buildMatrix(hic_filename, printer=False)
    binned_map = HiCtoolbox.bin2d(matrix, res, res, printer=False)
    old_shape = np.shape(matrix)[0]
    new_shape = np.shape(binned_map)
    
    #####Building the colors
    print("#####Building the colors")
    color_vec = HiCtoolbox.buildColors(epig_filename, chr_nb, old_shape, printer=False)
    color_bins = HiCtoolbox.bin2d(color_vec, res, 1, printer=False)
    color_bins = color_bins/np.amax(color_bins)
    
    #####Plot the matrix before filtering
    print("#####Plot the matrix before filtering")
    data = np.log2(1 + sparse.csr_matrix.toarray(binned_map))
    title = 'Matrix of chromosome '+nb+' from '+gene+' before filtering'
    path = result_path+'Matrix before filtering'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Matrix', nameFile=path)
    
    #####Plot the matrix after filtering
    print("#####Plot the matrix after filtering")
    filtered_map, bin_saved = HiCtoolbox.filtering(binned_map, filter_ratio, printer=False)    
    data = np.log2(1+sparse.csr_matrix.toarray(filtered_map))
    title = 'Matrix of chromosome '+nb+' from '+gene+' after filtering'
    path = result_path+'Matrix after filtering'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Matrix', nameFile=path)
    
    #####Plot the SCN matrix
    print("#####Plot the SCN matrix")
    scn_map = HiCtoolbox.SCN(filtered_map.copy(), printer=False)
    scn_map = np.asarray(scn_map)**alpha
    data = np.log2(scn_map)
    title = 'SCN matrix of chromosome '+nb+' from '+gene
    path = result_path+'SCN matrix'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Matrix', nameFile=path)
    
    #####Plot the distance matrix
    print("#####Plot the distance matrix")
    dist_matrix = HiCtoolbox.fastFloyd(1/scn_map, printer=False)
    dist_matrix = dist_matrix-np.diag(np.diag(dist_matrix))
    dist_matrix = (dist_matrix+np.transpose(dist_matrix))/2
    data = np.log(dist_matrix)
    title = 'Distance matrix of chromosome '+nb+' from '+gene
    path = result_path+'Distance matrix'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Matrix', nameFile=path)
    
    #####Plot the O/E matrix
    print("#####Plot the O/E matrix")
    contact_map = HiCtoolbox.OE(scn_map, printer=False)
    data = np.log(contact_map)
    vmin = -np.amax(data)
    vmax = np.amax(data)
    title = 'O/E matrix of chromosome '+nb+' from '+gene
    path = result_path+'OE matrix'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Matrix', cmap='seismic', vmin=vmin, vmax=vmax, nameFile=path)
    
    #####Plot the correlation matrix
    print("#####Plot the correlation matrix")
    corr_map = np.corrcoef(contact_map)
    data = corr_map
    vmin = -np.amax(data)
    vmax = np.amax(data)
    title = 'Correlation filtered matrix of chromosome '+nb+' from '+gene
    path = result_path+'Correlation filtered matrix matrix'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Matrix', cmap='seismic', vmin=vmin, vmax=vmax, nameFile=path)
    
    #####Plot the correlation unfiltered matrix
    print("#####Plot the correlation unfiltered matrix")
    unfiltered_corr_map = HiCtoolbox.unfiltering(bin_saved, corr_map, new_shape, printer=False)
    data = unfiltered_corr_map
    corr = data
    title = 'Correlation unfiltered matrix of chromosome '+nb+' from '+gene
    path = result_path+'Correlation unfiltered matrix matrix'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Matrix', cmap='seismic', vmin=vmin, vmax=vmax, nameFile=path)
    
    #####Plot the A/B compartments
    print("#####Plot the A/B compartments")
    data = HiCtoolbox.SVD(unfiltered_corr_map)
    title = 'AB compartments of chromosome '+nb+' from '+gene
    path = result_path+'AB compartments'
    HiCtoolbox.plotter(data, nameFig=title, plotType='AB', nameFile=path, centro_start=centro_start, centro_end=centro_end)
    #Barcode
    title = 'AB compartments barcode of chromosome '+nb+' from '+gene
    path = result_path+'AB compartments barcode'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Barcode', nameFile=path)
    
    #####Load the MDS
    print("#####Load the MDS")
    XYZ, E = HiCtoolbox.sammon(dist_matrix, 3, display=0) #with the one from tom j pollard
    hull = ConvexHull(XYZ)
    scale = 100/hull.area**(1/3)
    XYZ = XYZ*scale
    
    ###################################PREDICTIONS###################################
    print("###################################PREDICTIONS###################################")
    
    path = result_path+"Predictions.txt"
    f = open(path,"w") #prediction txt file to write data
    
    ##########Gaussian HMM with eigenvector from correlation matrix##########
    print("##########Gaussian HMM with eigenvector from correlation matrix##########")
    f.write('###################################HMM Contact###################################\n')
    vector = HiCtoolbox.SVD(corr_map).reshape(-1, 1)
    labels_Contact, scores_Contact = HiCtoolbox.multiplegaussianHMM(vector)
    Preds_Contact = labels_Contact
    data = scores_Contact
    title = 'HMM Score depending of chromosome '+nb+' from '+gene
    path = result_path+'HMM Contact compartments'
    HiCtoolbox.plotter(data, nameFig=title, plotType='HMMScore', nameFile=path, Data_type='Contact')
    
    #####HMM Barcode with eigenvector from correlation matrix
    print("#####HMM Barcode with eigenvector from correlation matrix")
    data = labels_Contact[0]
    vector = HiCtoolbox.SVD(unfiltered_corr_map).reshape(-1, 1).flatten()
    indexs = np.argwhere(vector == 0).flatten()
    data = list(data)
    for i in indexs:
        data.insert(i, -1.)
    data = np.array(data)
    with open(result_path+'Labels_contacts_2compartments.txt', "w") as txt_file:
    	for line in data:
            txt_file.write(str(line) + "\n")
    data = labels_Contact[0]
    data[data == 0.] = -1.
    vector = HiCtoolbox.SVD(unfiltered_corr_map).reshape(-1, 1).flatten()
    indexs = np.argwhere(vector == 0).flatten()
    data = list(data)
    for i in indexs:
        data.insert(i, 0.)
    data = np.array(data)
    if HiCtoolbox.similarity_score(val_data, data) > HiCtoolbox.similarity_score(val_data, -data):
        sim_score = HiCtoolbox.similarity_score(val_data, data)
    else:
        data = -data
        sim_score = HiCtoolbox.similarity_score(val_data, data)
    title = 'Eigenvector HMM barcode of chromosome '+nb+' from '+gene+' - Similarity score : '+str(sim_score)+' %'
    path = result_path+'HMM Contact Barcode'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Barcode', nameFile=path)
    
    #####Validation Barcode from Leopold Carron
    print("#####Validation Barcode from Leopold Carron")
    title = 'Validation barcode from Leopold Carron of chromosome '+nb+' from '+gene
    path = result_path+'HMM Contact Validation Barcode'
    HiCtoolbox.plotter(val_data, nameFig=title, plotType='Barcode', nameFile=path)
    
    #####Similarity Score with eigenvector from correlation matrix
    print("#####Similarity Score with eigenvector from correlation matrix")
    text = "Similarity score with Leopold: "+str(sim_score)+' %\n'
    f.write(text)
    
    #####Visualization of the results
    print("#####Visualization of the results")
    data = [data, corr, density]
    title = 'Results of chromosome '+nb+' from '+gene
    path = result_path+'HMM Contact All Results'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Visualization', nameFile=path, centro_start=centro_start, centro_end=centro_end)
    
    #####Save PDB File
    print("#####Save PDB File")
    best_nb_Contact = 2
    threshold = 1
    for i in range(1,len(scores_Contact)):
        temp = ((scores_Contact[i]-scores_Contact[i-1])/np.abs(scores_Contact[i]))*100
        if temp >= threshold: #we have found a best number of compartment -> score increase significantly
            best_nb_Contact+=1
        else:
            break #The Score no longer increases enough
    data = labels_Contact[best_nb_Contact-2]
    data = list(data)
    for i in indexs:
        data.insert(i, -1.)
    data = np.array(data)
    with open(result_path+'Labels_contacts_'+str(best_nb_Contact)+'compartments.txt', "w") as txt_file:
    	for line in data:
            txt_file.write(str(line) + "\n")
    Pred_HMM_Contact = labels_Contact[best_nb_Contact-2]
    HiCtoolbox.writePDB(result_path+'HMM_Contact.pdb',XYZ,Pred_HMM_Contact)
    f.write('Best compartment number : '+ str(best_nb_Contact)+'\n')
    
    #####Warsaw Code -> to validate our results
    x = np.arange(1,16)
    nameMark = ['']
    #Epigenetic marks color and labels
    Colors = []
    for i in range(0,len(x)):
        r = rd.random()  ;  b = rd.random()  ;  g = rd.random()  ;  color = (r, g, b)
        Colors.append(color)
        
    nameMark = ['TssA', 'TssAFlnk', 'TxFlnk', 'Tx', 'TxWk', 'EnhG', 'Enh', 'ZNF/Rpts', 
              'Het', 'TssBiv', 'BivFlnk', 'EnhBiv', 'ReprPC', 'ReprPCWk', 'Quies']
    
    Pred = Pred_HMM_Contact
    nb_cpt = best_nb_Contact
    freqs_Epi = [[0]*15]*nb_cpt
    for i in range(len(Pred)):
        Bin = Pred[i]
        if Bin==-1:
            Bin = 0
        freq = color_bins[:,][i].toarray()[0]
        freqs_Epi[Bin] = [x + y for x, y in zip(freqs_Epi[Bin], freq[1:])]

    if nb_cpt<=3:
        needed_lines = 1
    else:
        if nb_cpt%3==0:
            needed_lines = int(nb_cpt/3)
        else:
            needed_lines = 1+int((nb_cpt-nb_cpt%3)/3)

    plt.figure(figsize=(16,9))
    for i in range(len(freqs_Epi)):
        plt.subplot(needed_lines, 3, i+1)
        freq = freqs_Epi[i]
        norm_freq = freq/np.linalg.norm(freq)
        N, bins, patches = plt.hist(x, np.arange(1,len(x)+2) - 0.5, weights=norm_freq, edgecolor = 'red', align = 'mid', rwidth = 0.4)
        for j in range(0,len(x)):
            patches[j].set_facecolor(Colors[j])

        plt.xlabel('Marks', size=10)  ;  plt.xticks(np.arange(1,16))#, nameMark)
        plt.ylabel('Frequency', size=10)
        ax = plt.gca()  ;  ax.set_ylim([0, 1])
        plt.title("Compartment "+str(i+1), size=10)
        plt.tight_layout()

    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in Colors]
    plt.legend(handles, nameMark, bbox_to_anchor=(needed_lines+0.5,1))#bbox = (x,y)

    path = result_path+'Warsaw Validation_Contact'
    plt.savefig(path)
    
    ##########Gaussian HMM with expr/repr scores and epigenetic marks##########
    print("##########Gaussian HMM with expr/repr scores and epigenetic marks##########")
    f.write('\n###################################HMM Epigenetic###################################\n')
    marks = np.arange(1,16)
    vector = HiCtoolbox.expr_repr_scoring(color_bins, marks, expr_repr_scores)
    labels_Epi, scores_Epi = HiCtoolbox.multiplegaussianHMM(vector)
    Preds_Epi = labels_Epi
    data = scores_Epi
    title = 'HMM Score depending of chromosome '+nb+' from '+gene
    path = result_path+'HMM Epigenetic compartments'
    HiCtoolbox.plotter(data, nameFig=title, plotType='HMMScore', nameFile=path, Data_type='Epigenetic')
    
    #####HMM Barcode with expr/repr scores and epigenetic marks
    print("#####HMM Barcode with expr/repr scores and epigenetic marks")
    data = labels_Epi[0]
    vector = HiCtoolbox.SVD(unfiltered_corr_map).reshape(-1, 1).flatten()
    indexs = np.argwhere(vector == 0).flatten()
    data = list(data)
    for i in indexs:
        del data[i]
        data.insert(i, -1.)
    data = np.array(data)
    with open(result_path+'Labels_epi_2compartments.txt', "w") as txt_file:
    	for line in data:
            txt_file.write(str(line) + "\n")
    data = labels_Epi[0]
    data[data == 0.] = -1.
    vector = HiCtoolbox.SVD(unfiltered_corr_map).reshape(-1, 1).flatten()
    indexs = np.argwhere(vector == 0).flatten()
    data = list(data)
    for i in indexs:
        del data[i]
        data.insert(i, 0.)
    data = np.array(data)
    if HiCtoolbox.similarity_score(val_data, data) > HiCtoolbox.similarity_score(val_data, -data):
        sim_score = HiCtoolbox.similarity_score(val_data, data)
    else:
        data = -data
        sim_score = HiCtoolbox.similarity_score(val_data, data)
    title = 'Expr/Repr scores HMM barcode of chromosome '+nb+' from '+gene+' - Similarity score : '+str(sim_score)+' %'
    path = result_path+'HMM Epigenetic Barcode'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Barcode', nameFile=path)
    
    #####Similarity Score with expr/repr scores and epigenetic marks
    print("#####Similarity Score with expr/repr scores and epigenetic marks")
    text = "Similarity score with Leopold: "+str(sim_score)+' %\n'
    f.write(text)
    
    #####Visualization of the results
    print("#####Visualization of the results")
    data = [data, corr, density]
    title = 'Results of chromosome '+nb+' from '+gene
    path = result_path+'HMM Epigenetic All Results'
    HiCtoolbox.plotter(data, nameFig=title, plotType='Visualization', nameFile=path, centro_start=centro_start, centro_end=centro_end)
    
    #####Save PDB File
    print("#####Save PDB File")
    best_nb_Epi = 2
    threshold = 0.5
    for i in range(1,len(scores_Epi)):
        temp = ((scores_Epi[i]-scores_Epi[i-1])/np.abs(scores_Epi[i]))*100
        if temp >= threshold: #we have found a best number of compartment -> score increase significantly
            best_nb_Epi+=1
        else:
            break #The Score no longer increases enough
    Pred_HMM_Epigenetic = labels_Epi[best_nb_Epi-2]
    data = labels_Epi[best_nb_Epi-2]
    vector = HiCtoolbox.SVD(unfiltered_corr_map).reshape(-1, 1).flatten()
    indexs = np.argwhere(vector == 0).flatten()
    data = list(data)
    for i in indexs:
        del data[i]
        data.insert(i, -1.)
    data = np.array(data)
    if best_nb_Epi != 2:
        with open(result_path+'Labels_epi_'+str(best_nb_Epi)+'compartments.txt', "w") as txt_file:
            for line in data:
                txt_file.write(str(line) + "\n")
    HiCtoolbox.writePDB(result_path+'HMM_Epigenetic.pdb',XYZ,Pred_HMM_Epigenetic)
    f.write('Best compartment number : '+ str(best_nb_Epi)+'\n')
    
    #####Warsaw Code -> to validate our results
    Pred = Pred_HMM_Epigenetic
    nb_cpt = best_nb_Epi
    freqs_Epi = [[0]*15]*nb_cpt
    for i in range(len(Pred)):
        Bin = Pred[i]
        if Bin==-1:
            Bin = 0
        freq = color_bins[:,][i].toarray()[0]
        freqs_Epi[Bin] = [x + y for x, y in zip(freqs_Epi[Bin], freq[1:])]

    if nb_cpt<=3:
        needed_lines = 1
    else:
        if nb_cpt%3==0:
            needed_lines = int(nb_cpt/3)
        else:
            needed_lines = 1+int((nb_cpt-nb_cpt%3)/3)

    plt.figure(figsize=(16,9))
    for i in range(len(freqs_Epi)):
        plt.subplot(needed_lines, 3, i+1)
        freq = freqs_Epi[i]
        norm_freq = freq/np.linalg.norm(freq)
        N, bins, patches = plt.hist(x, np.arange(1,len(x)+2) - 0.5, weights=norm_freq, edgecolor = 'red', align = 'mid', rwidth = 0.4)
        for j in range(0,len(x)):
            patches[j].set_facecolor(Colors[j])

        plt.xlabel('Marks', size=10)  ;  plt.xticks(np.arange(1,16))#, nameMark)
        plt.ylabel('Frequency', size=10)
        ax = plt.gca()  ;  ax.set_ylim([0, 1])
        plt.title("Compartment "+str(i+1), size=10)
        plt.tight_layout()

    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in Colors]
    plt.legend(handles, nameMark, bbox_to_anchor=(needed_lines+0.5,1))#bbox = (x,y)

    path = result_path+'Warsaw Validation_Epi'
    plt.savefig(path)
    
    ##########Similarity##########
    print("##########Similarity##########")
    f.write('\n###################################Similarity###################################\n')
    #/!\ if nb_methods=1 -> we compare 2 methods
    #nb_scenarios = (nb_cpt^(nb_cpt+nb_methods))^2
    #scenario : ([labelcptPred1, labelcptPred2, ...], ...)
    nb_methods = 1
    nb_cpt = min(best_nb_Contact,best_nb_Epi) #In case the nb cpt predicted is not the same, we put more confidence in less compartment number
    Pred_HMM_Contact = labels_Contact[nb_cpt-2]
    Pred_HMM_Epigenetic = labels_Epi[nb_cpt-2]
    
    #####Compute all scenario (labels association)
    print("#####Compute all scenario (labels association)")
    All_possi = get_All_possiblesPaths(nb_cpt, nb_methods)
    All_Preds = [Pred_HMM_Contact,Pred_HMM_Epigenetic]
    All_tests = []
    for scenario in All_possi:
        similarity = calculateSimilarity(scenario, All_Preds, nb_cpt, nb_methods)
        All_tests.append([scenario, similarity])
    
    #####Find best scenario (labels association)
    print("#####Find best scenario (labels association)")
    best_scenario = find_bestScenario(All_tests, nb_cpt)
    f.write('Best Scenario : '+str(best_scenario[0])+'\n')
    f.write('Similarity : '+str(best_scenario[1])+'\n')
    
    #####Get Consensus labels
    print("#####Get Consensus labels")
    Consensus = []
    Consensus_labels = np.arange(len(best_scenario[0]))
    for i in range(len(Pred_HMM_Contact)):
        label1 = Pred_HMM_Contact[i]
        label2 = Pred_HMM_Epigenetic[i]
        flag = False #to check if we have found the association
        for j in range(len(best_scenario[0])):
            asso = best_scenario[0][j]
            if [label1,label2]==asso:
                Consensus.append(Consensus_labels[j])
                flag = True
                break
        if flag==False:
            Consensus.append(-100) #-100=arbitrary value to show that there is no consensus
        
    #####Save PDB
    print("#####Save PDB")
    HiCtoolbox.writePDB(result_path+'Consensus_PDB.pdb',XYZ,Consensus)
        
    print("###################################DONE###################################")
    
    #Delete temporary file created
    try:
        filename = "save_variable.txt"
        os.remove(filename) 
    except:
        pass
