#####################################
#########   Run Pyscenic   ##########
#####################################

#!/bin/sh

conda activate scenic

############ Variables ##############
## Repeat for every sample
ExpMatrix="/home/mainciburu/scRNA/scenic/data/mds1_norm_mat_full.csv"
TFs="/home/mainciburu/scRNA/scenic/resources/hs_hgnc_tfs.txt"
OutputName="MDS1"
#####################################

################ Run pyscenic ########################
# Step 1: Inference of co-expression modules
# Run GRNboost from arboreto to infer co-expression modules
# Output -> adjacencies. df with 3 columns TF / target / importance of the link
# **problems with dask versions => arboreto_with_multiprocessing.py


#pyscenic grn ${ExpMatrix} \    
#      ${TFs} \
#      -o "/home/mainciburu/scRNA/scenic/adjacencies/${OutputName}_adjacency.tsv" \
#      --num_workers 8 \
#      --transpose


echo -e "GRNboost ------------------------------- \n"
/home/mainciburu/scRNA/scenic/arboreto_with_multiprocessing.py \
    ${ExpMatrix} \
    ${TFs} \
    --method grnboost2 \
    --output "/home/mainciburu/scRNA/scenic/adjacencies/${OutputName}_adjacency.tsv" \
    --num_workers 6 \
    --transpose
    #--seed 777

# STEP 2-3: Regulon prediction aka cisTarget 
# Calculate a list of enriched motifs and the corresponding target genes for all modules
# Create regulons from this table of enriched motifs
# input
	# adjacency table
	# databases
	# motif annotations
	# expression matrix
# Output - regulon list
echo -e "Regulon prediction ------------------------------- \n"
pyscenic ctx "/home/mainciburu/scRNA/scenic/adjacencies/${OutputName}_adjacency.tsv" \
       "/home/mainciburu/scRNA/scenic/resources/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather" \
       "/home/mainciburu/scRNA/scenic/resources/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather" \
       --annotations_fname "/home/mainciburu/scRNA/scenic/resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
       --expression_mtx_fname ${ExpMatrix} \
       --mode "dask_multiprocessing" \
       --chunk_size 1 \
       --output "/home/mainciburu/scRNA/scenic/regulons/${OutputName}_regulons.gmt" \
       --num_workers 8 \
       -t 

# STEP 4: Cellular enrichment (aka AUCell)
echo -e "AUC ------------------------------- \n"
pyscenic aucell ${ExpMatrix} "/home/mainciburu/scRNA/scenic/regulons/${OutputName}_regulons.gmt" \
        --transpose \
        --num_workers 8 \
        -o "/home/mainciburu/scRNA/scenic/regulons/${OutputName}_regulons_AUCMat.csv"


conda deactivate
