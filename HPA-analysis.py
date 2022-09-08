'''
Ran code on madame-curie using tmux.
Input files link: https://www.proteinatlas.org/about/download
''' 

## CONVERT CSV TO ANNDATA OBJECT
import numpy as np
import anndata as ad
import array
import csv
from scipy.sparse import csr_matrix

input_file_name = 'rna_single_cell_read_count.tsv'

tissues=[]
cells=[]
clusters=[]
# array of type unsigned short
data = array.array("H")
# array of type unsigned int
indices = array.array("I")
indptr = array.array("I", [0])
for i, row in enumerate(csv.reader(open('rna_single_cell_read_count.tsv'),delimiter='\t')):
	if i==0:
		columns=row
	else:
		tissues.append(row[0])
        	cells.append(row[1])
        	clusters.append(row[2])
		row = np.array(row[3:],dtype=np.ushort)
		shape1 = len(row)
		nonzero = np.where(row)[0]
		data.extend(row[nonzero])
		indices.extend(nonzero)
		indptr.append(indptr[-1]+len(nonzero))
		if i%10000==0:
			print("Completed 10k rows")
matrix = csr_matrix((data, indices, indptr),
                      dtype=float, shape=(i, shape1))
df = {columns[0]:tissues, columns[1]:cells, columns[2]:clusters}
adata = ad.AnnData(X=matrix, obs=df, var=columns[3:])
