import pandas as pd
import numpy as np
import h5py

def pipeline(file_name):
    d, gene_dict, gene_to_idx, ph_to_idx = preprocessing(file_name)
    phenotypes = list(d.keys())
    gslib_genes = list(gene_dict.keys())
    bm = binary_matrix(phenotypes, gene_dict)
    gslib = gene_set_library(bm, gslib_genes, phenotypes, d, gene_to_idx)
    return new_gslib(phenotypes, gslib_genes, gslib)

def preprocessing(file_name):
    d = {} # convert text file into dict, with (phenotype, [genes])
    with open(file_name) as file:
        for line in file:
            lst = line.strip().split("\t")
            ph = lst[0]
            lst = lst[2:]
            d[ph] = lst
    gene_dict = {} # reverse of d 
    for k,v in d.items(): 
        for gene in v: 
            if gene not in gene_dict: 
                gene_dict[gene] = []
            gene_dict[gene].append(k) 
    functions = list(d.keys()) 
    genes = list(gene_dict.keys())
    gene_to_idx = {}
    f_to_idx = {}
    # to get index associated with a particular gene
    for i in range(len(genes)): 
        gene_to_idx[genes[i]] = i
    # to get index associated with particular function
    for j in range(len(functions)):
        f_to_idx[functions[j]] = j
    return d, gene_dict, gene_to_idx, f_to_idx
    

def binary_matrix(phenotypes, gene_dict):
    genes = list(gene_dict.keys())
    gene_set = np.zeros((len(genes), len(phenotypes)))
    for row in range(len(gene_set)): 
        gene = genes[row] 
        for col in range(len(phenotypes)):
            f = phenotypes[col]
            if f in gene_dict[gene]: 
                gene_set[row][col] = 1
    check = h5py.File("rna_level_predictions.hdf5", "r")
    mgi_gslib = check['data']['mgi_gslib']
    print(pd.DataFrame(mgi_gslib) == pd.DataFrame(gene_set))
    return gene_set


def gene_set_library(gene_set, genes, functions, d, gene_to_idx):
    cor = np.corrcoef(gene_set)
    # Initialize mouse gene set library 
    sim = np.zeros((len(genes), len(functions)))
    # Convert to Pandas DataFrame to easily use 
    # .iloc function, which allows row selection
    cor = pd.DataFrame(cor)
    # dict with (phenotype, [gene indices]) 
    gene_indices = {}
    for func in d:
        gene_indices[func] = [] 
        for gene in d[func]:
            gene_indices[func].append(gene_to_idx[gene])
    for j in range(len(functions)):
        f = functions[j]
        indices = gene_indices[f] 
        n = len(indices)-1
        temp = cor.iloc[indices]
        for i in range(len(indices)):
            gene_idx = indices[i]
            gene_cor = temp.iloc[:,i]
            sim[gene_idx][j] = (sum(gene_cor)-1)/n
    for row in sim:
        if sum(row) == 0:
            return "gslib broken"
    return sim 
    
    
def new_gslib(phenotypes, gslib_genes, mgsl):
    # Initialize new gene set library to contain TCGA genes as rows and mouse phenotypes as columns.
    fil = h5py.File("tcga.hdf5", "r")
    genes = fil['meta']['genes']
    corr = fil['correlation_matrix']
    gene_names = [ str(g[0])[2:-1] for g in genes ]
    ex_mgsl = np.zeros((len(gene_names), len(phenotypes)))
    # TCGA gene to index dictionary to help fill in expanded mouse gene set library 
    tcga_to_idx = {} 
    for g_idx in range(len(gene_names)): 
        g = gene_names[g_idx]
        tcga_to_idx[g] = g_idx
    for m in range(len(gslib_list)): 
        mouse_gene = gslib_list[m]
        if mouse_gene in tcga_to_idx:
            idx = tcga_to_idx[mouse_gene]
            ex_mgsl[idx] = mgsl[m]
    sums = sum(ex_mgsl)
    for s in sums:
        if s == 0: 
            return "new_gslib broken"
    new_gslib = np.matmul(corr, ex_mgsl)
    pheno_sums = []
    for col in np.transpose(new_gslib):
        pheno_sums.append(sum(col))
    for i in trange(len(new_gslib)):
        for j in range(len(phenotypes)):
            sub = ex_mgsl[i][j]
            denom = pheno_sums[j]
            new_gslib[i][j] = new_gslib[i][j]/(denom-sub)
    return new_gslib 
    
    
    
    
    
    
    
    
    
    
    
    