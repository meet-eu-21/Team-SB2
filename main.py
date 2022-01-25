from GUI_CPT import *
from Analysis import *


chosen_chroms, chosen_cell_line, chosen_resol = GUI_Start()

Cell_Lines = []
filter_ratio = 1.5
cell_lines = ['GM12878', 'HMEC', 'HUVEC', 'IMR90', 'NHEK' ]
for i in range(len(chosen_cell_line)):
    if chosen_cell_line[i]==True:
        Cell_Lines.append(cell_lines[i])
        
all_resol = ['100', '25']
Resols = []
for i in range(len(chosen_resol)):
    if chosen_resol[i]==True:
        Resols.append(all_resol[i])
        
chosen_chroms = chosen_chroms.split(',')
Chromos = []
for i in range(len(chosen_chroms)):
    try:
        Chromos.append(int(chosen_chroms[i]))
    except:
        continue

for cell_line in Cell_Lines:
    for chromo in Chromos:
        for resol in Resols:
            FOLDER = 'Data/'
            Data_available = create_Folders(cell_line, chromo, resol, FOLDER)
            FOLDER = 'Results/'
            Results_folders_created = create_Folders(cell_line, chromo, resol, FOLDER)
            if (Data_available==False): #download it
                Download_data(cell_line, chromo, resol)
            Analysis(cell_line, chromo, resol, filter_ratio)
