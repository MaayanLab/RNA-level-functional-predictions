import h5py
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score
from matplotlib import pyplot


def workflow(lib, tcga_genes, gene):
    f = h5py.File("auc_data.hdf5", "r")
    gslib, curr_genes, curr_pheno, binary_matrix = get_variables(lib, f)
    tcga_idx = get_idx(tcga_genes, gene)
    lib_idx = get_idx(curr_genes, gene)
    top_pheno = get_top_phenotypes(gslib[tcga_idx], curr_pheno)
    if lib_idx: auc(binary_matrix[lib_idx], gslib[tcga_idx], gene)
    else: print("Not enough gene annotations available.")
    f.close()
    return top_pheno

def get_variables(curr, f):
    data = f['data']
    meta = f['meta']
    gslib = data['tcga_' + curr + 'gslib']
    curr_genes = [str(g[0])[2:-1] for g in meta[curr + 'genes']]
    curr_pheno = [str(g[0])[2:-1] for g in meta[curr + 'pheno']] 
    binary_matrix = data[curr + "bin_mat"]
    return gslib, curr_genes, curr_pheno, binary_matrix


def get_idx(gene_list, gene):
    return np.where(np.array(gene_list) == gene)[0][0]
    

def get_top_phenotypes(row, curr_pheno):
    tups = list(dict(zip(range(len(row)), row)).items())
    tups.sort(key = lambda t: t[1], reverse = True)
    top_rank = [ t[0] for t in tups[: 50] ]
    top_pheno = [ curr_pheno[i] for i in top_rank ]
    return top_pheno 


def auc(y_true, probs, gene):
    ns_probs = [0 for _ in range(len(y_true))]
    fpr, tpr, _ = roc_curve(y_true, probs)
    ns_fpr, ns_tpr, _ = roc_curve(y_true, ns_probs)
    auc = roc_auc_score(y_true, probs)
    pyplot.plot(ns_fpr, ns_tpr, linestyle='--')
    pyplot.plot(fpr, tpr, marker='.')
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    pyplot.title(gene)
    pyplot.text(0.75, 0.05, 'AUC: %.3f' % auc, fontsize=12)
    pyplot.show()
    
