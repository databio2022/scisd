#!/user/bin/python3
# from counts matrix removing doublet cells
#self.py countsmatrixfile 
import pandas as pd
import scrublet as scr
import sys

counts_matrix=pd.read_csv(sys.argv[1],sep="\t",index_col=0)
counts_matrix=counts_matrix.T
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
X=counts_matrix
X['doublet']=predicted_doublets
X=X.loc[X['doublet'] == False]
X=X.drop(columns = ['doublet'])
X.T.to_csv(sys.argv[1].replace(".txt","_doubletfilter.txt"),sep="\t")

