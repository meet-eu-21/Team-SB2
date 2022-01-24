from GUI_CPT import *
from Analysis import *


chosen_chroms, chosen_cell_line, chosen_resol = GUI_Start()

filters = 
cell_lines = ['GM12878', 'HMEC', 'HUVEC', 'IMR90', 'NHEK' ]
Cell_Lines = [1.8,2.6,2.8,2.6,2.4,2.6,2.4,2.6,1.4,2.4,2.4,2.7,1.45,1.6,1.35,1.72,2.1,1.96,1.95,1.85,1.2,1.21,2.4]
for i in range(len(chosen_cell_line)):
    if chosen_cell_line[i]==True:
        Cell_Lines.append(cell_lines[i])
        
all_resol = ['100', '25']
Resols = []
for i in range(len(chosen_resol)):
    if chosen_resol[i]==True:
        Resols.append(all_resol[i])
        
if 'all' in chosen_chroms:
    chosen_chroms = range(1, 23)
else:       
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
            nb = str(nb)
    	    if chromo=='X':
        	filter_ratio = filters[23-1]
            else:
            	filter_ratio = filters[chromo-1]
            Analysis(cell_line, chromo, resol, filter_ratio)
