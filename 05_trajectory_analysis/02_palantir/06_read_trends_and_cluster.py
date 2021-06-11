from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import palantir
import scanpy as sc

adata = sc.read(filename="/home/mainciburu/scRNA/palantir/processed/mds5.h5ad")

# Read trends 
file_pr = open('/home/mainciburu/scRNA/palantir/results/mds5/pr_res.obj', 'rb') 
pr_res = pickle.load(file_pr)

results = OrderedDict()
res_path = "/home/mainciburu/scRNA/palantir/gene_trends/mds5/"

lineages = ["Erythroid_late"]
for branch in lineages:
	results[branch] = OrderedDict()
	trend_file = res_path + branch + "_" + "trends.csv"
	std_file = res_path + branch + "_" + "std.csv"
	results[branch]['trends'] = pd.read_csv(trend_file, index_col = 0)
	results[branch]['std'] = pd.read_csv(std_file, index_col = 0)

# Cluster trends 

for branch in lineages:
	trends = results[branch]['trends']   
	gene_clusters = palantir.presults.cluster_gene_trends(trends)
	csvname="/home/mainciburu/scRNA/palantir/results/mds5/clusters_"+branch+".csv"
	gene_clusters.to_csv(csvname)
	palantir.plot.plot_gene_trend_clusters(trends, gene_clusters)
	figname="/home/mainciburu/scRNA/palantir/pics/clusters_"+branch+"_mds5.pdf"
	plt.savefig(figname)