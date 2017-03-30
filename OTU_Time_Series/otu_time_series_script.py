#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:10:12 2016

@author: Keith 

First, we will load in the data:

Then we can add some alpha diversity metrics, introduce the CLR transform,
and drop a bunch of columns that are either bad replicates or not informative.

The three main efforts of the analysis below include: 
    1. What are the correlates of the PCA components? Why does Chao1 so tightly
       correlate to the main PCA component? Why is it significant that depth
       tightly correlates with the second principal component? 
    2. What are the combination of priors and components maximize the log-
       likelihood of a Dirichlet Mixture Model for OTUs ? How do the topics
       and their distributions 
    3. What do the beta diversity comparisons tell us about how stable the 
       community within the lake is over time ? 

       
http://nbviewer.jupyter.org/github/dsquareindia/gensim/blob/
a4b2629c0fdb0a7932db24dfcf06699c928d112f/docs/notebooks/topic_coherence_tutorial.ipynb


https://jonlefcheck.net/2015/03/31/how-much-is-enough
-a-new-technique-for-quantifying-precision-of-community-surveys/

http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

"""
import pickle
import ecopy as ep
from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys, os
from ChemDataVectors import LoadChemDFs
from importCoverage import loadFullMetadata
from otu_ts_support import listRepGroups, ReplicateReport, varStabTransform
from otu_ts_support import numericEncodings# , beta_wrapper
from otu_ts_support import importratesandconcentrations_obs, importratesandconcentrations_mod
from otu_ts_support import scalePCAcorrelate, parseBiosample#, plotHeatmap
from otu_ts_support import add_quadrants, plotCountTotalsByMetadata
from otu_ts_support import analyze_alpha_diversity, alpha_diversity
from otu_ts_support import plotCommonTaxa, dropBadReps, collapseBdiversity
from otu_ts_support import centeredLogRatio, JensenShannonDiv_Sqrt
from otu_ts_support import time_scale_modeled_chem_data, append_depths
from sklearn.decomposition import TruncatedSVD #, LatentDirichletAllocation

# Relevant File Names OTU Table
tree_file='unique.dbOTU.tree'
seq_file='unique.dbOTU.nonchimera.fasta'
otu_matrix_file = 'unique.dbOTU.nonchimera.mat.rdp'

print "\nReading in OTU Table"
otu_table = pd.read_csv(otu_matrix_file, sep="\t", index_col=0)
taxa_series = otu_table.ix[:, -1]

print "\nImporting metadata"
metadata_df = loadFullMetadata(verbose=False)

print "\nDropping -DO NOT USE- samples and those not associated with depths"
# Drop Recommended Columns
otu_time_table = otu_table.drop([otu_table.columns[-1]], axis=1)
do_not_use = ['SB100912TAWMD14VV4TMR1', 'SB061713TAWMD22VV4TMR1',
              'SB011413TAWMD22VV4TMR1', 'SB011413TAWMDSBVV4TMR1',
              'SB011413TAWMDEBVV4TMR1', 'SB011413TAWMD22VV4TMR2',
              'SB011413TAWMDSBVV4TMR2']
              
otu_time_table.drop(do_not_use, axis=1, inplace=True)

print "\nParsing biosample strings into metadata columns"
# Add Biosample metadata & convert date format
samps_by_feats = parseBiosample(otu_time_table.T)
undated = samps_by_feats[samps_by_feats.date == "NA"].index
samps_by_feats.drop(undated, inplace=True)
samps_by_feats['date'] = pd.to_datetime(samps_by_feats['date'])
squishy_depths = ['bottom', 'SB', 'control1', 'control3', 'control6', 'mid1', 
                  'mid2', 'mid3', 'river', 'upper', 'River', 'Neg', 'NA', 'EB', 
                  'FB', 'CR', 'Dock', 'Bridge']

# Drop samples of unknown provenance & sort by depth and date
off_target = samps_by_feats[samps_by_feats.depth.isin(squishy_depths)].index
only_depths = samps_by_feats.drop(off_target, inplace=False)
only_depths.sort_values(['date', 'depth'], ascending=[True, True],
                        inplace=True)

print "\nIdenitfying replicate groups"
# Identify replicate groupings & create OTU only matrix
rep_groups = listRepGroups(only_depths)
metadata = ['date', 'primers', 'kit', 'replicates', 'depth']
only_depth_otus = only_depths.drop(metadata, axis=1)

print "\nAdding quadrant metadata"
# Add lake quadrant metadata variable
only_depths_q = add_quadrants(only_depths)
metadata.append("Quadrants")

print "\n Adding remaining metadata by sample ID"
# set the sample id as the index of the metadata_df 
metadata_df.set_index('#SampleID', inplace=True)
# drop all rows that don't correspond to the existing index
lost_idxs = set(list(metadata_df.index)).difference(set(list(only_depths_q.index)))
metadata_df.drop(list(lost_idxs), inplace=True, axis=0)
# concatenate along columnar axis 
only_depths_all = only_depths_q.join(metadata_df)
# add categorical columns to the list for encoding
metadata += ['Sequencing Date', 'Sequencing platform', 'Forward read length',
             'Index read length']

print "\nNumerically encoding categorical metadata"
# Numerically encode categoricals
only_depths_n, decoder = numericEncodings(only_depths_all, metadata, True)
metadata += ['Coverage', 'TotalSeqs', 'BadBarcodes', 'GoodBarcodes']

print "\nCalculating Alpha Diversity"
# Add metadata columns pertaining to alpha diversity
alphaList = ['shannon', 'chao1', 'enspie']
only_depth_otus, only_depths_n = alpha_diversity(only_depth_otus, 
                                                 only_depths_n,
                                                 alphaList)
metadata+=alphaList

print "\nTransforming counts into centered log ratios and model based methods"
# Perform centered log-ratio transform 
clr_T_otus, clr_T_m = centeredLogRatio(only_depth_otus, only_depths_n)

to_transform = "/Users/login/Desktop/shared_otus.csv"

tmm_vst_otus, tmm_vst_m = varStabTransform(to_transform, only_depth_otus, 
                                           only_depths_n, 'TMM')
rle_vst_otus, rle_vst_m = varStabTransform(to_transform, only_depth_otus, 
                                           only_depths_n, 'RLE')

print "\nChecking for vetted replicates file"

shitty_reps_path = "shitty_reps.p"
do_rep_rep = False
if not os.path.exists(shitty_reps_path) or do_rep_rep:
    print "\nMeasuring B-diversity (sqrt(JSD)) distances between replicates"
    shitty_reps, broken_groups = ReplicateReport(only_depths_n, 
                                                       only_depth_otus, 
                                                       rep_groups, False, "JSD")
        
    kept_pct = broken_groups/float(len(rep_groups))
    print "\nGroups with replicates removed: {0:.3f}%".format(kept_pct*100)
    pickle.dump( shitty_reps, open( shitty_reps_path, "wb" ) )
else:
    print "\tLoading pre-approved replicates"
    shitty_reps = pickle.load( open( shitty_reps_path , "rb" ) )
    

print "\nDropping most distant replicates"
dr_tmm_vst_otus = tmm_vst_otus.drop(shitty_reps)
dr_tmm_vst_m = tmm_vst_m.drop(shitty_reps)
dr_rle_vst_otus = rle_vst_otus.drop(shitty_reps)
dr_rle_vst_m = rle_vst_m.drop(shitty_reps)
derepped_clr_m = clr_T_m.drop(shitty_reps)
dereplicated_clr = clr_T_otus.drop(shitty_reps)
derepped_otu_m = only_depths_n.drop(shitty_reps)
dereplicated_otu = only_depth_otus.drop(shitty_reps)

new_rep_groups = dropBadReps(shitty_reps, rep_groups)

print "\nDescriptive Statistics of Alpha Diversity for var. Sample Groupings"

valPairs = [(1,1),(1,2),(2,1),(2,2)]
            
alpha_df, primer_outgroup = analyze_alpha_diversity(decoder, derepped_otu_m, 
                                                    valPairs)

print "\nDropping V4V5 primer samples"
ingroup_tmm_vst_otus = dr_tmm_vst_otus.drop(primer_outgroup)
ingroup_tmm_vst_m = dr_tmm_vst_m.drop(primer_outgroup)
ingroup_rle_vst_otus = dr_rle_vst_otus.drop(primer_outgroup)
ingroup_rle_vst_m = dr_rle_vst_m.drop(primer_outgroup)
ingroup_otu_m = derepped_otu_m.drop(primer_outgroup)
ingroup_otu = dereplicated_otu.drop(primer_outgroup)
ingroup_clr_m = derepped_clr_m.drop(primer_outgroup)
ingroup_clr = dereplicated_clr.drop(primer_outgroup)

print "\Lets see how the signal holds up through PCA:"
_ = scalePCAcorrelate(ingroup_otu, ingroup_otu_m, metadata, True)
_ = scalePCAcorrelate(ingroup_tmm_vst_otus, ingroup_tmm_vst_m, metadata, True)
_ = scalePCAcorrelate(ingroup_rle_vst_otus, ingroup_rle_vst_m, metadata, True)
_ = scalePCAcorrelate(ingroup_clr, ingroup_clr_m, metadata, True)

X_std2 = ingroup_clr.values
svd = TruncatedSVD(n_components=2, n_iter=20, random_state=42)
y_std2 = svd.fit_transform(X_std2)

with plt.style.context('seaborn-whitegrid'):
    plt.figure(4, figsize=(9, 9))
    for lab, lab_en, col in zip(('Q1', 'Q2', 'Q3', 'Q4'),range(1,5),
                        ('blue', 'red', 'green', 'magenta')):
        plt.scatter(y_std2[(ingroup_clr_m.Quadrants==lab_en).values, 0],
                    y_std2[(ingroup_clr_m.Quadrants==lab_en).values, 1],
                    label=lab,c=col, s=40)
    ax = plt.gca()
    for axis in [ax.xaxis, ax.yaxis]:
        for tick in axis.get_major_ticks():
            tick.label.set_fontsize(16)
    plt.xlabel('Principal Component 1', fontsize=24)
    plt.ylabel('Principal Component 2', fontsize=24)
    plt.legend(loc='lower left', fontsize=16)
    plt.tight_layout()
    plt.show()

print "\nPlot Alpha Diversity Heat Map"
a_metrics = ['shannon', 'chao1', 'enspie']
for a_metric in a_metrics:
    depths = np.unique(ingroup_otu_m.depth)
    dates = np.unique(ingroup_otu_m.date)
    alpha_plot_df = pd.DataFrame(index=depths, columns=dates)
    for i in ingroup_otu_m.index:
        this_depth = ingroup_otu_m.ix[i, 'depth']
        this_date = ingroup_otu_m.ix[i, 'date']
        alpha_plot_df.ix[this_depth, this_date] = ingroup_otu_m.ix[i, a_metric]
    
    decoded_depths = [str(decoder['depth'][i]) for i in depths]
    decoded_dates = [str(decoder['date'][i]).split("T")[0] for i in dates]
    alpha_plot_df.index = decoded_depths
    alpha_plot_df.columns = decoded_dates
    
    alpha_plot_df.to_csv(a_metric+'.csv', index_label='depths')
    # make an autoplotting function for 
    
#plotHeatmap(alpha_plot_df, 2)
#al_ax = plt.gca()
#al_ax.set(xlabel='date', ylabel='depth')

print "\nCleaning out any OTUs that only appeared in dropped samples"
var_thresh = (.9 * (1 - .9))
otu_var = ingroup_otu.var(axis=0)
var_otu_var = otu_var[otu_var < var_thresh].index
infrequent_otus = list(var_otu_var)

print "\nVariance Thresholding OTUS"

shared_otus_m = ingroup_otu_m.drop(infrequent_otus, axis=1)
shared_otus = ingroup_otu.drop(infrequent_otus, axis=1)
shared_clr = ingroup_clr.drop(infrequent_otus, axis=1)
shared_clr_m = ingroup_clr_m.drop(infrequent_otus, axis=1)
notzero = (shared_otus != 0).sum().sum()
sparsity = float(notzero) / shared_otus.values.flatten().shape[0]
print "Probability that a given OTU val is nonzero: {}".format(sparsity)

rate_conc_f = os.path.join(os.getcwd(), 'ChemData', 'mystic_model_data')
obs_conc_f = os.path.join(os.getcwd(), 'ChemData', 'mystic_measured_data')

depthGrad_df, surfaceM_df = LoadChemDFs()
surfaceM_df['date'] = surfaceM_df.index
depthGrad_df.set_index(['date', 'depth'], inplace=True)
    
surfaceM_df['depth'] = np.zeros((surfaceM_df.shape[0], ))
depths_vector = np.arange(22)
surfaceM_df = append_depths(surfaceM_df, depths_vector)
surfaceM_df.set_index(['date', 'depth'], inplace=True)

obs_conc_df = importratesandconcentrations_obs(obs_conc_f)

rate_dict, conc_dict = importratesandconcentrations_mod(rate_conc_f)
model_proc_df = time_scale_modeled_chem_data(rate_dict, conc_dict, 146, 
                                             '03/23/2013', '08/15/2013')

shared_clr_m = ingroup_clr_m.drop(infrequent_otus, axis=1)
def transform_date_to_periodic(date_time):
    return (pd.to_datetime(decoder['date'][date_time]).month - 6.) / float(6)
def transform_date(date_num):
    return pd.to_datetime(decoder['date'][date_num])
def transform_depth(depth_num):
    return float(decoder['depth'][depth_num])
    
shared_clr_m['seasonality'] = shared_clr_m.date.apply(transform_date_to_periodic)
shared_clr_m['Date'] = shared_clr_m.date.apply(transform_date)
shared_clr_m['depth'] = shared_clr_m.depth.apply(transform_depth)

clr_idx = shared_clr_m.set_index(['Date', 'depth'])
remaining_metadata = clr_idx.columns[-16:]
otu_columns = clr_idx.columns[:-16]
clr_y = clr_idx.drop(otu_columns, 1)
clr_x = clr_idx.drop(remaining_metadata, 1)

parade_of_responses = [model_proc_df, obs_conc_df, surfaceM_df, depthGrad_df, clr_y]

#for baton_twirler in parade_of_responses:

# TODO: -Ask how to aggregate model data (monthly? daily? intermediate locality? )
# TODO: set up properly formatted files for LSA analysis

from sklearn.model_selection import train_test_split as tts
from sklearn.metrics import mean_squared_error, r2_score, explained_variance_score
from sklearn.linear_model import RandomizedLasso as rlasso
from sklearn.ensemble import RandomForestRegressor as rForestReg

baton_twirler = clr_y
score_fxns = ['r2', 'evs', 'mse']
score_df = pd.DataFrame(data = np.zeros((len(baton_twirler.columns), len(score_fxns))),
                        index = baton_twirler.columns,
                        columns = ['r2', 'evs', 'mse'])

feat_dict = {}
for baton in baton_twirler.columns:
    y = baton_twirler.ix[:, baton].values
    X = clr_x.values
    assert list(baton_twirler.index) == list(clr_x.index)
    
    X_tr, X_te, y_tr, y_te = tts(X, y, test_size=0.33)
    
    lookatTheTrees = rForestReg(n_estimators=X_tr.shape[1],
                                min_samples_split = 2,
                                n_jobs=2, verbose=10)
    
    lookatTheTrees.fit(X_tr, y_tr)
    y_pred = lookatTheTrees.predict(X_te)
    this_mse = mean_squared_error(y_te, y_pred)
    this_r2 = r2_score(y_te, y_pred)
    this_evs = explained_variance_score(y_te, y_pred)
    
    feat_dict[baton] = (lookatTheTrees.estimators_, 
                        lookatTheTrees.feature_importances_)
    score_df.ix[baton, :] = np.array([this_r2, this_evs, this_mse])
    
sys.exit()

rlasso_params = {}
rlasso_params['alpha'] = np.exp(np.arange(-2,2, 0.1))
rlasso_params['l1_ratio'] = np.arange(0,1.01,0.01)
rlasso_params['normalize'] = [True, False]
rlasso_params['selection'] = ['random']
rlasso_cols = ['alpha', 'l1_ratio', 'normalize', 'EVR', 'R2', 'MSE']
rlasso_data = {}
trial_counter = 0
for alpha in rlasso_params['alpha']:
    for l1_ratio in rlasso_params['l1_ratio']:
        for normalize in rlasso_params['normalize']:
            trial_counter+=1
            mses_, evrs_, r2s_ = [],[],[]
            for folds in range(5):                
                X_tr, X_te, y_tr, y_te = tts(X, y, test_size=0.33)            
                an_rlasso = rlasso(alpha=alpha, l1_ratio=l1_ratio, normalize=normalize, 
                               selection=rlasso_params['selection'])
                y_pred_rlasso = an_rlasso.fit(X_tr).predict(X_te)
                r2s_.append(r2_score(y_te, y_pred_rlasso))
                evrs_.append(explained_variance_score(y_te, y_pred_rlasso))
                mses_.append(mean_squared_error(y_te, y_pred_rlasso))
                
            mse_arr, evr_arr, r2_arr = np.array(mses_), np.array(evrs_), np.array(r2s_)
            
            rlasso_data[trial_counter] = (alpha, l1_ratio, normalize, 
                                          evr_arr.mean(), r2_arr.mean(), 
                                          mse_arr.mean())


            
rlasso
            













print "Loading Stored OTU matrices if available"
shared_otu_path = 'shared_otu.csv'
jsd_dist_path = 'jsd_dist_mat.csv'
bc_dist_path = 'bray_curtis_dist_mat.csv'


if os.path.exists(shared_otu_path):
    print "Loading OTU Table used for Distance Calculations"
    loaded_otu = pd.read_csv(shared_otu_path, index_col=0)
else:
    print "No previous distance calculation attempted"
    shared_otus.to_csv(shared_otu_path)
    loaded_otu = pd.DataFrame()
    
if loaded_otu.equals(shared_otus):
    print "\n Input Equivalence detected"
    if os.path.exists(jsd_dist_path):
        print "Loading Jensen Shannon Distance matrix" 
        jsd_mat = pd.read_csv(jsd_dist_path, index_col=0)
    else:
        print "Calculating sqrt of Jensen Shannon Divergence of all rows"
        jsd_mat = JensenShannonDiv_Sqrt(shared_otus)
        jsd_mat.to_csv(jsd_dist_path)
        print "Saving distance matrix and shared otu matrix"
        
    if os.path.exists(bc_dist_path):
        print "Loading Bray-Curtis Distances"
        brayDF = pd.read_csv(bc_dist_path, index_col=0)
    else:
        print "Calculating sqrt of Jensen Shannon Divergence of all rows"
        brayDist = ep.distance(shared_otus, method='bray')
        brayDF = pd.DataFrame(data=brayDist,
                      index=shared_otus.index,
                      columns=shared_otus.index)
        brayDF.to_csv(bc_dist_path)
        print "Saving distance matrix and shared otu matrix"
    del loaded_otu
else:
    print "\n Previous calculations derived from distinct OTU table"
    print "Calculating sqrt of Jensen Shannon Divergence of all rows"
    jsd_mat = JensenShannonDiv_Sqrt(shared_otus)
    print "Caclulating Bray Curtis Distance of all rows"
    brayDist = ep.distance(shared_otus, method='bray')
    brayDF = pd.DataFrame(data=brayDist,
                  index=shared_otus.index,
                  columns=shared_otus.index)
    print "Saving distance matrices"
    brayDF.to_csv(bc_dist_path)
    jsd_mat.to_csv(jsd_dist_path)

mD_date, mDev_date = collapseBdiversity(jsd_mat, shared_otus_m, 'date')
mD_depth, mDev_depth = collapseBdiversity(jsd_mat, shared_otus_m, 'depth')

plt.figure(5, figsize=(12,9))
mD_date_df = pd.DataFrame(data=mD_date, index=decoded_dates, columns=decoded_dates)
ax = sns.heatmap(mD_date_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)
plt.title("Mean sqrt(JSD) distance by date")
    
plt.figure(6, figsize=(12,9))
mDev_date_df = pd.DataFrame(data=mDev_date, index=decoded_dates, columns=decoded_dates)
ax = sns.heatmap(mDev_date_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)   
plt.title("Std Dev of sqrt(JSD) distance by date")

plt.figure(7, figsize=(12,9))
mD_depth_df = pd.DataFrame(data=mD_depth, index=decoded_depths, 
                           columns=decoded_depths)
ax = sns.heatmap(mD_depth_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)   
plt.title("Mean sqrt(JSD) distance by depth (meters)")

plt.figure(8, figsize=(12,9))
mDev_depth_df = pd.DataFrame(data=mDev_depth, index=decoded_depths, 
                             columns=decoded_depths)
ax = sns.heatmap(mDev_depth_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)   
plt.title("Std Dev sqrt(JSD) distance by depth (meters)")

bcD_date, bcDev_date = collapseBdiversity(brayDF, shared_otus_m, 'date')
bcD_depth, bcDev_depth = collapseBdiversity(brayDF, shared_otus_m, 'depth')

bcD_date_df = pd.DataFrame(data=bcD_date, index=decoded_dates, 
                          columns=decoded_dates)
bcDev_date_df = pd.DataFrame(data=bcDev_date, index=decoded_dates,
                            columns=decoded_dates)
bcD_depth_df = pd.DataFrame(data=bcD_depth, index=decoded_depths, 
                           columns=decoded_depths)
bcDev_depth_df = pd.DataFrame(data=bcDev_depth, index=decoded_depths, 
                             columns=decoded_depths)

plt.figure(9, figsize=(12,9))
ax = sns.heatmap(bcD_date_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)
plt.title("Mean Bray Curtis distance by date")
    
plt.figure(10, figsize=(12,9))
ax = sns.heatmap(bcDev_date_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)   
plt.title("Std Dev of Bray Curtis distance by date")

plt.figure(11, figsize=(12,9))
ax = sns.heatmap(bcD_depth_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)   
plt.title("Mean Bray Curtis distance by depth (meters)")

plt.figure(12, figsize=(12,9))
ax = sns.heatmap(bcDev_depth_df)
ylabs, xlabs = ax.get_yticklabels(), ax.get_xticklabels()
for ytem, xtem in zip(ylabs, xlabs):
    ytem.set_rotation(0)
    xtem.set_rotation(90)   
plt.title("Std Dev Bray Curtis distance by depth (meters)")


mD_depth_df.to_csv('Mean_root_JSD_distance_by_depth.csv')
mD_date_df.to_csv('Mean_root_JSD_distance_by_date.csv')
sys.exit()



print "\nCreating All Plots"

plotCommonTaxa(taxa_series)

plt.figure(3, figsize=(12,9))
colors_ = list(colors.cnames)
for i in [7, 10, 11]:
    del colors_[i]
for date_, col in zip(np.unique(ingroup_otu_m.date), colors_, ):
    date_bool = ingroup_otu_m.date == date_
    date_subdf = ingroup_otu_m[date_bool]
    if date_subdf.shape[0] > 5:
        plt.scatter(date_subdf['depth'].values, 
                 date_subdf['enspie'].values, 
                 c=col, label=str(decoder['date'][date_]).split("T")[0])
plt.legend(loc='best')



print "\nCreate sample summary plots"
plotCountTotalsByMetadata(shared_otus_m, decoder, metadata, 'depth', 6)
plotCountTotalsByMetadata(shared_otus_m, decoder, metadata, 'date', 7)

#print "Plotting UPGMA clustering of Unifrac Distance between samples (skipped)"
#from cStringIO import StringIO
#from Bio import Phylo
#handle = StringIO(unifracDendogram)
#tree = Phylo.read(handle, "newick")
#tree.ladderize()   # Flip branches so deeper clades are displayed at top
#Phylo.draw(tree)


#print "\n Performing LDA decomposition"
#
#lda = LatentDirichletAllocation(n_topics=3, doc_topic_prior=0.2, 
#                                topic_word_prior=0.1, learning_method='batch', 
#                                random_state=42, max_iter=20)
#X = shared_otus.values
#lda.fit(X)
#topic_word_dist = lda.components_
#twd_df = pd.DataFrame(data=topic_word_dist, 
#                      index = ['t1', 't2', 't3'], 
#                      columns = shared_otus.columns)

#distances = JensenShannonDiv_Sqrt(twd_df)
#twd_df_t = twd_df.T
#abundance = shared_otus.mean(axis=0)
#twd_df_t['importance'] = abundance
#twd_df_t.sort_values(['importance'], ascending=False, inplace=True)
# nDMS of topic seperation
# evaluate mixture model redundancy and distribution across time points


sys.exit()



