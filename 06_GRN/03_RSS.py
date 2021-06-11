########### RSS ########################
import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.cli.utils import save_enriched_motifs

from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize


auc_mtx = pd.read_csv("/home/mainciburu/scRNA/scenic/regulons/AML_regulons_AUCMat.csv",index_col=0).T
CellTypes = pd.read_csv("/home/mainciburu/scRNA/scenic/data/aml_CellType.txt" ,index_col=0, sep = " ")
CellTypes.columns=["prediction"]
rss_cellType = regulon_specificity_scores(auc_mtx, CellTypes['prediction'] )
rss_cellType = rss_cellType.T
rss_cellType.to_csv('/home/mainciburu/scRNA/scenic/regulons/AML_rssCelltype.csv',index=True,header=True)

