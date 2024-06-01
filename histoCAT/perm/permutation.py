import scipy.io
import numpy as np 
import math
import itertools
import numpy.matlib
import time
import multiprocess as mp
import glob
import os 

'''
function handling phenograph computations
'''
def pheno(Neighbor_Matrix, allLabels, allLabelsInd, nOfCells, cellTypes):
    Neighbor_Matrix_index = Neighbor_Matrix + 1
    combos_oneside = np.array(list(itertools.combinations(allLabelsInd, 2)))
    combos_oneside = combos_oneside[:,:,0]

    flipped = np.flip(combos_oneside, axis=1)
    toSelf = np.concatenate((allLabelsInd, allLabelsInd), axis = 1)
    combos_all = np.concatenate((combos_oneside, flipped, toSelf), axis = 0)

    Phenograph_Vector = np.zeros((nOfCells[0][0],1))

    for i in range(cellTypes.shape[0]):
        # print(i)
        if cellTypes[i][0].size <= 0:
            index = -1
            if 'NONE' in np.unique(allLabels.tolist()):
                index = np.argwhere(allLabels=='NONE')
        else:
            curType = cellTypes[i][0][0]
            ## HANDLE error in celltype names when needed
            # if curType == "Cl MDM":
            #     index = np.argwhere(allLabels=="CL MDM")
            # elif curType == "Cl Mono":
            #     index = np.argwhere(allLabels=="Cl Mo")
            # else:
            index = np.argwhere(allLabels==curType)
        # check if allLabels is horizontal or vertical
        if allLabels.shape[1] > allLabels.shape[0]:
            Phenograph_Vector[i] = allLabelsInd[index[0][1]]
        else:
            Phenograph_Vector[i] = allLabelsInd[index[0][0]] 

    Phenograph_Vector_index = np.concatenate(([[0]], Phenograph_Vector), axis = 0)

    d = dict(enumerate(Phenograph_Vector_index[:,0],1))
    Phenograph_Neighor_Matrix = np.vectorize(d.get)(Neighbor_Matrix_index)
    return Neighbor_Matrix_index, combos_all, Phenograph_Vector, Phenograph_Neighor_Matrix


'''
function computing all histcounts with inputted celltype combos 
'''
def Calculate_STDandMean(combos_all, Phenograph_Neighor_Matrix, Phenograph_Vector):
    
    combos_all_histcount = np.zeros((len(combos_all), 3))
    
    for i in range(len(combos_all)):
        idx = np.where(Phenograph_Neighor_Matrix == combos_all[i,0])
        l1 = idx[0]
        idx = np.where(Phenograph_Vector == combos_all[i,1])
        l2 = idx[0]
        
        Get_unique_rows = set(l1) & set(l2)
        combos_all_histcount[i, :2] = combos_all[i, :]
        
        if len(Get_unique_rows) <= 0:
            combos_all_histcount[i, 2] = 0
        else:
            intersectNeighbours = Phenograph_Neighor_Matrix[list(Get_unique_rows), :]
            eachCount = np.sum(intersectNeighbours == combos_all[i,0], axis=1, dtype=np.int64)
            combos_all_histcount[i, 2] = np.sum(eachCount) / len(Get_unique_rows)
            if not combos_all_histcount[i, 2]:
                combos_all_histcount = 0
                
    return combos_all_histcount

'''
function performing permutation with specified number of iterations 
'''
def permutation(permutations):
    combos_all_histcount_Perm_multiple = np.zeros((len(combos_all),permutations))
    for p in range(permutations):
        Phenograph_Vector_perm = Phenograph_Vector[np.random.permutation(len(Phenograph_Vector))]
        Phenograph_Vector_index_perm = np.concatenate(([[0]], Phenograph_Vector_perm), axis = 0)
        Phenograph_Neighor_Matrix_perm = Phenograph_Vector_index_perm[Neighbor_Matrix_index-1]
        Phenograph_Neighor_Matrix_perm = Phenograph_Neighor_Matrix_perm[:,:,0]
        combos_all_histcount_Perm_
        single = Calculate_STDandMean(combos_all,Phenograph_Neighor_Matrix_perm,\
            Phenograph_Vector_perm)
        combos_all_histcount_Perm_multiple[:,p] = combos_all_histcount_Perm_single[:,2]
    return combos_all_histcount_Perm_multiple


if __name__ == '__main__':

    ## config parameters
    patch_det = 0 
   
    permutations_lst = [1_000, 5_000, 10_000, 20_000, 30_000, 50_000]
    # permutations_lst = [50_000, 30_000, 20_000, 10_000, 5_000, 1_000]
    # permutations_lst = [30_000]

    # num_processes = mp.cpu_count()
    num_processes = 16

    for permutations in permutations_lst:
        print("-----Running "+str(permutations)+" permutations-----")
        perm_perProcess = int(permutations/num_processes)+1

        f_src = '../inter/inter_nada/*.mat'
        f_save_root = '../save/nada/nada-'+ str(int(permutations/1000)) + 'k/'
        if not os.path.exists(f_save_root):
            os.makedirs(f_save_root)

        fileList = []
        for f in glob.glob(f_src):
            fileList.append(f)

        for idx, files in enumerate(fileList):
            
            ## load the computed matrices from matlab
            scanName = files.split('/')[-1][:-4]
            f_mat = f_src.replace('*', scanName)
            inputs = scipy.io.loadmat(f_mat)

            for k in inputs.keys():
                if not k.startswith("_"):
                    vars()[k] = inputs[k]
                    
            ## phenoGraph computations 
            [Neighbor_Matrix_index, combos_all, Phenograph_Vector, Phenograph_Neighor_Matrix] = pheno(Neighbor_Matrix, allLabels, allLabelsInd, nOfCells, cellTypes)
            combos_all_histcount_real = Calculate_STDandMean(combos_all, Phenograph_Neighor_Matrix, Phenograph_Vector)

            ## running permutations with multiprocessing
            start_time = time.time()
            with mp.Pool(num_processes) as pool:
                t = pool.map(permutation,[perm_perProcess]*num_processes)
            print("file %d takes %s seconds ---" % (idx,(time.time() - start_time)))
            combos_all_histcount_Perm = np.hstack(t)
            combos_all_histcount_Perm = combos_all_histcount_Perm[:,:permutations]
            print(combos_all_histcount_Perm.shape)
            
            ## process perm results and save to .mat files 
            real_data_mean = combos_all_histcount_real[:,2]
            perm_data_mean = combos_all_histcount_Perm
            # real_data_mean.shape

            Higher_perm_test = np.transpose(np.matlib.repmat(real_data_mean,perm_data_mean.shape[1],1)) <= perm_data_mean
            Lower_perm_test  = np.transpose(np.matlib.repmat(real_data_mean,perm_data_mean.shape[1],1)) >= perm_data_mean
            Amount_higher = np.sum(Higher_perm_test, axis = 1)
            Amount_lower = np.sum(Lower_perm_test, axis = 1)

            pValue_higher = (Amount_higher+1)/(permutations+1)
            pValue_lower = (Amount_lower+1)/(permutations+1)

            f_save = f_save_root + scanName + '_perm.mat'

            mdic = {"pValue_higher": pValue_higher, 
                    "pValue_lower": pValue_lower, 
                    "real_data_mean" : real_data_mean, 
                    "combos_all" :  combos_all}
            scipy.io.savemat(f_save, mdic)
        

