#### STEAM Trajectory Analysis


# Load Required Packages
import stream as st
import numpy as np

data_path = "/home/mainciburu/scRNA/stream/data/"
res_path = "/home/mainciburu/scRNA/stream/results/young/"

Young=st.read(file_name=data_path + 'ExpData.csv',
              file_feature=data_path + 'Genes.tsv',
              file_sample=data_path + 'barcodes.tsv',
              delimiter=',',workdir=res_path)    

# Load the Metadata
st.add_metadata(Young,file_name=data_path + 'metadata.tsv')

# select variable genes (we use 2000) and PCA analysis
st.select_variable_genes(Young,loess_frac=0.009,n_genes=2000,save_fig=True)
st.select_top_principal_components(Young,feature='var_genes',first_pc=True,n_pc=20,save_fig=True)

# Dimensional reduction - mlle
st.dimension_reduction(Young,method='mlle',feature='top_pcs',n_components=3,n_neighbors=30,n_jobs=12)
st.plot_dimension_reduction(Young,color=['label'],show_graph=False,show_text=False,save_fig=True)

# Plot the UMAP
st.plot_visualization_2D(Young,method='umap',n_neighbors=75,color=['label'],save_fig=True,use_precomputed=False)

# Create the EPG 
st.seed_elastic_principal_graph(Young,n_clusters=20)
st.elastic_principal_graph(Young,epg_alpha=0.02,epg_mu=0.05,epg_lambda=0.01)

# Plot dimenasional reduction with graph
st.plot_dimension_reduction(Young,color=['label'],n_components=3,show_graph=True,show_text=False,save_fig=True)
st.plot_branches(Young,show_text=True,save_fig=True)

# Plot flat tree to select root branch
st.plot_flat_tree(Young,color=['label','branch_id_alias'], dist_scale=0.5,show_graph=True,show_text=True,save_fig=True)

st.plot_stream(Young,root='S1',color=['label'],save_fig=True,log_scale=True,dist_scale=0.5,log_scal=True)
st.plot_stream(Young,root='S1',color=['CD79A','HBB','CSTA','GATA1','EBF1','IRF8'],save_fig=True,log_scale=True,dist_scale=0.5)

st.detect_transition_markers(Young,marker_list=Young.uns['var_genes'],cutoff_spearman=0.4,cutoff_logfc=0.25,root='S1',n_jobs=8)

st.write(Young,file_name=res_path + 'young.pkl')


########## Senior mapping ###########
Young = st.read('/home/mainciburu/scRNA/stream/results/young/stream_result_Young.pkl')

res_path = '/home/mainciburu/scRNA/stream/results/senior_mapping/'

Senior=st.read(file_name=data_path + 'ExpDataSenior.csv',
              file_feature=data_path + 'GeneSenior.tsv',
              file_sample=data_path + 'barcodesSenior.tsv',
              delimiter=',', workdir=res_path)    

st.add_metadata(Senior,file_name=data_path + 'metadataSenior.tsv')

Integrated = st.map_new_data(Young,Senior)

Integrated.obs.head()

Integrated.obs['Condition'] = ''
Integrated.obs.loc[Young.obs_names+'-ref','Condition'] = 'Young'
Integrated.obs.loc[Senior.obs_names+'-new','Condition'] = 'Senior'
Integrated.uns['Condition_color'] = {'Young':'#b8342b','Senior':'#5775b0'}

Integrated.obs['label_new'] = ''
Integrated.obs.loc[Young.obs_names+'-ref','label_new'] = 'ref'
Integrated.obs.loc[Senior.obs_names+'-new','label_new'] = Integrated.obs.loc[Senior.obs_names+'-new','label']
Integrated.uns['label_new_color'] = {'ref':'gray',**Senior.uns['label_color']}

st.plot_dimension_reduction(Integrated,color=['Condition','label_new'],show_graph=True,show_text=False,save_fig=True)
st.plot_visualization_2D(Integrated,method='umap',n_neighbors=50,color=['Condition','label_new'],save_fig=True,use_precomputed=False)
st.plot_flat_tree(Integrated,color=['Condition','label_new'], dist_scale=0.5,show_graph=True,show_text=True,save_fig=True)
st.plot_stream(Integrated,root='S1',color=['label_new','Condition'],dist_scale=0.5,save_fig=True,log_scale=True)

st.plot_stream(Integrated[Integrated.obs['Condition']=='Senior'],root='S1',color=['label_new'],dist_scale=0.5,save_fig=True,log_scale=True, fig_size = (7,3))

st.write(Integrated,file_name=res_path + 'stream_result_integrated_senior.pkl')


##### MDS mapping #########

Young = st.read('/home/mainciburu/scRNA/stream/results/stream_result_Young_jp.pkl')

res_path = '/home/mainciburu/scRNA/stream/results/MDS_mapping/'

MDS=st.read(file_name=data_path + 'ExpDataMDS5.csv',
              file_feature=data_path + 'GeneMDS5.tsv',
              file_sample=data_path + 'barcodesMDS5.tsv',
              delimiter=',',workdir=res_path)    

st.add_metadata(MDS,file_name=data_path + 'metadataMDS5.tsv')


Integrated = st.map_new_data(Young,MDS)

Integrated.obs.head()

Integrated.obs['Condition'] = ''
Integrated.obs.loc[Young.obs_names+'-ref','Condition'] = 'Young'
Integrated.obs.loc[MDS.obs_names+'-new','Condition'] = 'MDS'
Integrated.uns['Condition_color'] = {'Young':'#FC766AFF','MDS':'#4e0191'}

Integrated.obs['label_new'] = ''
Integrated.obs.loc[Young.obs_names+'-ref','label_new'] = 'ref'
Integrated.obs.loc[MDS.obs_names+'-new','label_new'] = Integrated.obs.loc[MDS.obs_names+'-new','label']
Integrated.uns['label_new_color'] = {'ref':'gray',**MDS.uns['label_color']}

st.plot_dimension_reduction(Integrated,color=['Condition','label_new'],show_graph=True,show_text=False,save_fig=True)
st.plot_visualization_2D(Integrated,method='umap',n_neighbors=50,color=['Condition','label_new'],save_fig=True,use_precomputed=True)
st.plot_flat_tree(Integrated,color=['Condition','label_new'], dist_scale=0.5,show_graph=True,show_text=True,save_fig=True)
st.plot_stream(Integrated,root='S1',color=['label_new','Condition'],dist_scale=0.5,save_fig=True,log_scale=True)

st.plot_stream(Integrated[Integrated.obs['Condition']=='MDS'],root='S1',color=['label_new'],dist_scale=0.5,save_fig=True,log_scale=True, fig_size = (7,3))

st.write(Integrated,file_name=res_path + 'stream_result_integrated_MDS5.pkl')
