import os
import numpy as np
from itertools import chain, combinations

#Building tree of possibilities
class TreeNode:
    def __init__(self,root):
        self.root = root
        self.leaves = []
        
def count_leaves(tree):
    nb_leaves = 0
    for i in range(len(tree.leaves)):
        if tree.leaves[i].leaves==[]:
            nb_leaves+=1
        else:
            nb_leaves+=count_leaves(tree.leaves[i])
    return nb_leaves
    
def add_leaves(tree, nb_cpt, nb_methods, depth=1):
    if depth==nb_methods:
        for i in range(nb_cpt):
            tree.leaves.append(TreeNode(i))
        return tree
    else:
        for i in range(nb_cpt):
            tree.leaves.append(TreeNode(i))
        for j in range(nb_cpt):
            tree.leaves[j] = add_leaves(tree.leaves[j], nb_cpt, nb_methods, depth=depth+1)
    return tree

def create_tree(nb_cpt, nb_methods):
    All_trees = []
    for i in range(nb_cpt):
        tree = TreeNode(i)
        tree = add_leaves(tree, nb_cpt, nb_methods)
        All_trees.append(tree)
    return All_trees
    
# function to get all path from root to leaf
def get_Paths(tree):
    # list to store path
    path = []
    get_PathsRec(tree, path, 0)
    return None

#Helper function to get path from root to leaf
def get_PathsRec(tree, path, pathLen):
    #print(All_paths)
    # if length of list is gre
    if(len(path) > pathLen):
        path[pathLen] = tree.root
    else:
        path.append(tree.root)
 
    # increment pathLen by 1
    pathLen = pathLen + 1
 
    if tree.leaves==[]:
        file_save = open("save_variable.txt", "a")
        str_save = repr(path)
        file_save.write(str_save + "\n")
        file_save.close()
        return path
    else:
        # try for each subtree
        for subtree in tree.leaves:
            path2 = get_PathsRec(subtree, path, pathLen)
            
    return path2
    
def loadPaths(filename):
    Paths = []
    lecture = np.loadtxt(filename, dtype=object)
    for line in lecture:
        path = []
        for cpt in line:
            node = ''
            for letter in cpt:
                try:
                    node+=str(int(letter))
                except:
                    pass
            path.append(int(node))
        Paths.append(path)
    return Paths

def all_subsets(ss):
        return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))
    
def get_All_possiblesPaths(nb_cpt, np_methods):
    tree = create_tree(nb_cpt,np_methods)
    All_Paths = []
    filename = "save_variable.txt"
    for cpt_subtree in tree:
        get_Paths(cpt_subtree)
        Paths = loadPaths(filename)
        os.remove(filename)
        All_Paths.append(Paths)
    
    All_possi = []

    joined_list = []
    for path in All_Paths:
        joined_list += path
    
    for subset in all_subsets(joined_list):
        if len(subset)==nb_cpt:
            if len(set(np.array(subset)[:,0]))==nb_cpt:
                All_possi.append(subset)
    
    return All_possi
    
def calculateSimilarity(scenario, All_Preds, nb_cpt, nb_methods):
    #calculate the similarity of 2 methods for a scenario given (labels association)
    dico_count = {}
    for cpt in range(len(All_Preds[0])):
        current_scenario = []
        for i in range(len(All_Preds)):
            current_scenario.append(All_Preds[i][cpt])
        for i in range(len(current_scenario)):
            if current_scenario[i]==-1:
                current_scenario[i] = 0
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
    
def find_bestScenario(All_tests, nb_cpt):
    #find best scenario of labels association
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