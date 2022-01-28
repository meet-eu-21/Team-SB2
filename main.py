from GUI_CPT import *
from Analysis import *

#Main file to launch the graphic interface and the analysis:
# Be careful -> Selecting to many chromosomes at a time can
# make runtime memory error /!\
#
# Be careful -> Sometimes you need to delete the files in
# /Data and /Results to switch between resolution and make
# new analysis of the same chromosome

#Launch the GUI
chosen_chroms, chosen_cell_line, chosen_resol = GUI_Start()

Cell_Lines = []
filter_ratio = 1.5 #Filtering ratio used

#Select the input cell lines
cell_lines = ['GM12878', 'HMEC', 'HUVEC', 'IMR90', 'NHEK' ]
for i in range(len(chosen_cell_line)):
    if chosen_cell_line[i]==True:
        Cell_Lines.append(cell_lines[i])
   
#Select the input resolution     
all_resol = ['100', '25']
Resols = []
for i in range(len(chosen_resol)):
    if chosen_resol[i]==True:
        Resols.append(all_resol[i])

#Select the input chromosomes
chosen_chroms = chosen_chroms.split(',')
Chromos = []
for i in range(len(chosen_chroms)):
    try:
        Chromos.append(int(chosen_chroms[i]))
    except:
        continue

#Main loop :
for cell_line in Cell_Lines: #Select the cell line
    for chromo in Chromos: #Select the chromosome
        for resol in Resols: #Select the resolution
            #Create the Data folders
            FOLDER = 'Data/'
            Data_available = create_Folders(cell_line, chromo, resol, FOLDER)
            #Create the Results folders
            FOLDER = 'Results/'
            Results_folders_created = create_Folders(cell_line, chromo, resol, FOLDER)
            #Download the data if you don't have it already
            if (Data_available==False): 
                Download_data(cell_line, chromo, resol)
            #Launch the analysis and save the results
            Analysis(cell_line, chromo, resol, filter_ratio)
