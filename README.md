# Team-SB2

Meet-EU Team SB2 - 
**[Link to MEET-EU website](https://hdsu-bioquant.github.io/meet-eu-2021/)**

Topic B : Chromosome compartments

[Digest](https://github.com/meet-eu-21/Team-SB2/blob/main/Team_SB2_-_Digest.pdf)

[Report](https://github.com/meet-eu-21/Team-SB2/blob/main/Team_SB2_-_Report.pdf)


# Chromosome compartments project

Project of the MEET-U course for our university __Sorbonne Université__

> **Meetings :**
> One meeting each week every wednesday afternoon from 4pm to 6pm

**Responsability list :**
Role | Name
------------ | -------------
Manager Expert | Rouquaya Mouss
Technician Expert | Damien Legros & Cédric Cornede
Scientific Expert | Hamid Hachemi
Delivrable Expert | Arnaud Quelin

# Install and launch

First you need to clone our repository :
```
git clone https://github.com/meet-eu-21/Team-SB2
```
After that you have to install the requirements :
```
pip install -r requirements.txt
```
Then to launch the code :
```
python main.py
```
The following panel will come up, you can enter the numbers of the chromosomes you want to analyse and select the cell lines and the resolution **(be careful ! : launching too many chromosomes can make a runtime memory error)**.

![panel](https://github.com/meet-eu-21/Team-SB2/blob/main/images/panel.png?raw=true)

Here you have an example of the selection to run an analysis for **GM12878**, **HUVEC** and **IMR90** with the chromosomes **11**, **X** and **4** in **100kb**.

The notebook **Plot_histo.pynb** can then be used to make the histograms of the results to compare the cell lines.

# Our work

Class diagram of our code :

![diagram](https://github.com/meet-eu-21/Team-SB2/blob/main/images/class_diagram.png?raw=true)

Our results for 100kb are available in the folder **/Results**

All **.py** files except **color_b.py** are used for our analysis. Two methods are used :
- With the first method, we search the number of compartments with an HMM using the data from the correlation matrix
- With the second method, we search the number of compartments with an HMM using the data from epigenetic marks.

The file **color_b.py** can be used with **Pymol** to see the colors of the compartments of the **.pdb** files of our results. (**color_b.py** can be put in **Pymol** folder, no need of the path before color_b.py in the following command if that's the case). You need to run the following commands in **Pymol** after opening the **.pdb** result you want to see : 
```
run path_to_our_folder\Team-SB2\color_b.py
color_b()
```
Here is an example of the colors that can be seen on the **.pdb** file of chromosome **22** of **HMEC** with the compartments of the HMM method with correlation matrix : 
![pdb](https://github.com/meet-eu-21/Team-SB2/blob/main/images/chr22_HMEC_HMMContact.png?raw=true)



# Example of results

## Example of data obtained with the correlation matrix method of chromosome 16 from GM12878 :

![corrcontact](https://raw.githubusercontent.com/meet-eu-21/Team-SB2/main/Results/GM12878/16/100/HMM%20Contact%20All%20Results.png)

## Example of data obtained with the epigenetic marks method of chromosome 16 from GM12878 :

![correpigenetic](https://raw.githubusercontent.com/meet-eu-21/Team-SB2/main/Results/GM12878/16/100/HMM%20Epigenetic%20All%20Results.png)

### Our predictions :
```
###################################HMM Contact###################################
Similarity score with Leopold: 74.972 %
Best compartment number : 7

###################################HMM Epigenetic###################################
Similarity score with Leopold: 64.673 %
Best compartment number : 5
```

### HMM score with the correlation matrix method :

![hmmcontact](https://github.com/meet-eu-21/Team-SB2/blob/main/Results/GM12878/16/100/HMM%20Contact%20compartments.png?raw=true)

### HMM score with epigenetic marks method :

![hmmepigenetic](https://github.com/meet-eu-21/Team-SB2/blob/main/Results/GM12878/16/100/HMM%20Epigenetic%20compartments.png?raw=true)

## Example of the histogram of results for **HUVEC** cell line :

### For the correlation matrix method :

![hmmcontact](https://github.com/meet-eu-21/Team-SB2/blob/main/Results/Results%20-%20Histograms/HUVEC/100/Best_cpt_Contact.jpg?raw=true)


![hmmcontact2](https://github.com/meet-eu-21/Team-SB2/blob/main/Results/Results%20-%20Histograms/HUVEC/100/similarity_Leo_Contact.jpg?raw=true)

### For the epigenetic method :

![hmmepigenetic](https://github.com/meet-eu-21/Team-SB2/blob/main/Results/Results%20-%20Histograms/HUVEC/100/Best_cpt_Expr-Repr.jpg?raw=true)

![hmmepigenetic2](https://github.com/meet-eu-21/Team-SB2/blob/main/Results/Results%20-%20Histograms/HUVEC/100/similarity_Leo_Expr-Repr.jpg?raw=true)
