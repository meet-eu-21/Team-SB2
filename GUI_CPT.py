from tkinter import *
import sys
import time
import os
from PIL import ImageTk, Image

def GUI_Start():
    """
    GUI starter to launch the graphic interface
    Buttons for :
    - Choosing the chromosomes, (23 for X) use comma to separate the chromosomes (ex: 23, 2, 15)
    - Choosing the cell lines, tick the cell line you want to study
    - Choosing the resolution, iick the resolution you want to study
    - Launching the analysis
    - Quitting the interface
    """
    #Windows options
    root = Tk(className="CPT Predict")
    #root.state('zoomed')
    root.attributes("-topmost", True)
    
    #Chromosomes to be chosen
    l = Label( root, text = "Enter all chromosomes you want to analyse\n\n separate each chromosome with a comma (23 for X)" ).pack()
    Chromo = Entry( root, bd = 5 )  ;  Chromo.pack()

    #Cells lines to be chosen
    l = Label( root, text = "\nSelect the cell lines you want to analyse" ).pack()
    all_cells_line = ['GM12878', 'HMEC', 'HUVEC', 'IMR90', 'NHEK' ]
    check_cell = []
    for i in range(len(all_cells_line)):
        var = BooleanVar()  ;  var.set(False)    
        check_button = Checkbutton(root, variable=var, text=all_cells_line[i], width=10)
        check_button.pack()
        check_cell.append(var)
    
    #Resolution to be chosen
    l = Label( root, text = "\nSelect the Resolution you want to analyse" ).pack()
    all_resol = ['100kb', '25kb' ]
    check_resol = []
    for i in range(len(all_resol)):
        var = BooleanVar()  ;  var.set(False)    
        check_button = Checkbutton(root, variable=var, text=all_resol[i], width=10)
        check_button.pack()
        check_resol.append(var)
    
    Send = Button(root, text = "Process this/these chromosome(s)", command=lambda:choose(Chromo, check_cell, check_resol, root)).pack()
    
    #Exit
    Exit = Button(root, text="Quitter", fg="red", command=root.destroy).pack()
    
    root.mainloop()
    return chosen_chrom, cell_line, resol
    
def choose(Chromo, check_cell, check_resol, root):
    """
    Function to update the inputs made in the graphic interface
    Return :
    - Chromosomes (List of strings)
    - Cell lines (List of strings)
    - Resolution (List of strings)
    """
    global chosen_chrom  ;  global cell_line  ;  global resol
    chosen_chrom = Chromo.get()  ;  cell_line = []  ;  resol = []
    for i in check_cell:
        cell_line.append(i.get())
    for i in check_resol:
        resol.append(i.get())
    root.destroy()
    return chosen_chrom, cell_line, resol

def update(ind, frames, frameCnt, label, root):
    """
    Function to update the graphic interface when an input is made
    Update in realtime the selected buttons and written text
    """
    frame = frames[ind]
    ind += 1
    if ind == frameCnt:
        ind = 0
    label.configure(image=frame)
    root.after(100, update, ind, frames, frameCnt, label, root)
