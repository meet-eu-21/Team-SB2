# Team-SB2

Meet-EU Team SB2 - 
**[Link to MEET-EU website](https://hdsu-bioquant.github.io/meet-eu-2021/)**

Topic B : Chromosome compartments

# Chromosome compartments project

Project of the MEET-U course for ou university __Sorbonne Université__



> **Meetings :**
> One meeting each week every wednesday afternoon from 4pm to 6pm

**Responsability list :**
Role | Name
------------ | -------------
Manager Expert | Rouquaya Mouss
Technician Expert | Damien Legros & Cédric Cornede
Scientific Expert | Hamid Hachemi
Delivrable Expert | Arnaud Quelin

# TODO List

- Add descriptions of each function in *.pynb* notebook
- Add cleaned versions of the autoencoder and similarity tree to the *.pynb* notebook
- Clean *main.py* file and use HiCtoolbox functions
- Add comments to every functions in *HiCtoolbox.py* file
- Make the similarity tree inside the *HiCtoolbox.py* file
- Add the autoencoder to *HiCtoolbox.py* file
- Make the autoencoder with expr/repr scores
- Make the *main.py* file work with **argparse**
- Make the *main.py* work automatically with the all *.RAWobserved* files
- Generate the documentation of the HiCtoolbox
- Separate the functions in the HiCtoolbox to make it cleaner

# Readme

For now, the best is to launch our notebook (*Compartment SB2 - Jupyter Script.ipynb*)

To work on a chromosome from a gene, you need to add 3 files to the folder. For example if you want to work on chromosome 18 of gene GM12878 :
- *chr18.hdf5*
- *chr18_100kb.RAWobserved*
- *chr18_compartiment.txt*

Respectively, the density data, the HiC data and the validation data

You then need to change the variables in **Variable to change** in the notebook
- nb = '18'
- gene = 'GM12878'

The filtering ratio (filter_ratio) is also define in **Variable to change**, if you have some errors it is probably because the filtering ratio is not small enough and you should reduce it to resolve the errors.

# Example of results

Example of data obtained with the corrélation matrix of chromosome 18 from GM12878 :
![Results]()
