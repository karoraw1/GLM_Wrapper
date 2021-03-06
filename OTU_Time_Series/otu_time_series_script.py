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
from otu_ts_support import score_clusters, matchXandYbyIndex, extract_linkages
def date_to_periodic(date_time):
    return (pd.to_datetime(decoder['date'][date_time]).month - 6.) / float(6)
def transform_date(date_num):
    return pd.to_datetime(decoder['date'][date_num])
def transform_depth(depth_num):
    return float(decoder['depth'][depth_num])
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
from otu_ts_support import listRepGroups, ReplicateReport, edgeRtmm, DESeqRlog, standardScaled
from otu_ts_support import numericEncodings, originate_rep_groupings# , beta_wrapper
from otu_ts_support import importratesandconcs_obs, importratesandconcs_mod
from otu_ts_support import scalePCAcorrelate, parseBiosample, plotHeatmap
from otu_ts_support import add_quadrants, plotCountTotalsByMetadata
from otu_ts_support import analyze_alpha_diversity, alpha_diversity
from otu_ts_support import plotCommonTaxa, dropBadReps, collapseBdiversity
from otu_ts_support import centeredLogRatio, JensenShannonDiv_Sqrt
from otu_ts_support import time_scale_modeled_chem_data, append_depths
from sklearn.decomposition import TruncatedSVD #, LatentDirichletAllocation
from otu_ts_support import bz2wrapper, simpleProportions, rarefyOTUtables
from otu_ts_support import plot_interreplicate_distances, distant_replicate_covariates


# Relevant File Names OTU Table
tree_file='unique.dbOTU.tree'
seq_file='unique.dbOTU.nonchimera.fasta'
otu_matrix_file = 'unique.dbOTU.nonchimera.mat.rdp'

chk_pt_test_f = 'shared_test_x.csv'
chk_pt_train_f = 'shared_test_x.csv'
skip_ahead = False
checkpoint = os.path.exists(chk_pt_test_f) and os.path.exists(chk_pt_train_f)
    
print "\nReading in OTU Table"
if not os.path.exists(otu_matrix_file):
    # if unzipped not present add prefix and unzip
    otu_to_unzip = otu_matrix_file + ".bz2"
    zOut = bz2wrapper(otu_to_unzip)
else:
    # if present, make sure to check for archive later to zip
    otu_to_unzip = otu_matrix_file + ".bz2"
    
otu_table = pd.read_csv(otu_matrix_file, sep="\t", index_col=0)

#if not os.path.exists(otu_to_unzip):
#    zOut = bz2wrapper(otu_matrix_file)

taxa_series = otu_table.ix[:, -1]

if (skip_ahead == False):
    print "\nCheckpoint requirements not met"
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
    plot_interreplicate_distances(only_depth_otus, rep_groups, 66)

    print "\nMeasuring B-diversity (sqrt(JSD)) distances between replicates"
    shitty_reps, broken_groups_n = ReplicateReport(only_depths, 
                                                 only_depth_otus, 
                                                 rep_groups, False, "JSD")
        
    kept_pct = broken_groups_n/float(len(rep_groups))
    flat_replicates = [item for sublist in rep_groups for item in sublist]
        
    print "\nAdding quadrant column and full metadata df"
    # Add lake quadrant metadata variable
    only_depths_q = add_quadrants(only_depths)
    metadata.append("Quadrants")
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
    
    print "Adding solo replicates to shitty replicates & binning distances"
    new_rep_groups, broken_groups = dropBadReps(shitty_reps, rep_groups)   
    meta_dict_u = distant_replicate_covariates(broken_groups, only_depth_otus, 
                                               only_depths_all, metadata)
    covars = [(i, j.mean(), len(j)) for i, j in meta_dict_u.items()]
    covars.sort(key=lambda x: (x[1]), reverse=True)
    
    print "\nDropping most distant replicates"
    all_shitty = [i for sl in broken_groups for i in sl]
    alt_kept_pct = len(all_shitty) / float(len(flat_replicates))
    print "\nGroups with replicates removed: {0:.3f}%".format(kept_pct*100)
    print "% of all replicates dropped : {0:.3f}%".format(alt_kept_pct*100)
    dereplicated_otu = only_depth_otus.drop(all_shitty)
    dereplicated_m = only_depths_all.drop(all_shitty)
    
    jsd_dist_mat_fn = "jsd_dist_full.p"
    if os.path.exists(jsd_dist_mat_fn):
        dist_mat = pickle.load( open(jsd_dist_mat_fn, "rb" ) )
    else:
        dist_mat = JensenShannonDiv_Sqrt(dereplicated_otu)
        pickle.dump( dist_mat, open(jsd_dist_mat_fn, "wb") )
        
    plt.figure(67)
    plt.clf()
    plt.hist(dist_mat.values.flatten(), facecolor='g', bins=100)
    plt.tick_params(labelsize=14)
    plt.xlabel("root JS distance (all)", fontsize=14)
    plt.ylabel("N", fontsize=14)
    
    print "\nNumerically encoding categorical metadata"
    # Numerically encode categoricals
    metadata += ['Coverage', 'TotalSeqs', 'BadBarcodes', 'GoodBarcodes']
    dereplicated_m2, decoder = numericEncodings(dereplicated_m, metadata, False)
    metadata += ['seasonality']
    dereplicated_m2['seasonality'] = dereplicated_m2.date.apply(date_to_periodic)
    
    print "\nCalculating Alpha Diversity"
    # Add metadata columns pertaining to alpha diversity
    alphaList = ['shannon', 'chao1', 'enspie']
    derepped_otus, derepped_m = alpha_diversity(dereplicated_otu, 
                                                dereplicated_m2,
                                                alphaList)
    metadata+=alphaList

    print "\nDescriptive Statistics of Alpha Diversity for var. Sample Groupings"
    
    valPairs = [(1,1),(1,2),(2,1),(2,2)]
    chao1_df, primer_outgroup = analyze_alpha_diversity(decoder, derepped_m, 
                                                        valPairs, 'chao1')
    enspie_df, primer_outgroup = analyze_alpha_diversity(decoder, derepped_m, 
                                                        valPairs, 'enspie')
    
    ingroup_otu_m = derepped_m.drop(primer_outgroup)
    ingroup_otu = derepped_otus.drop(primer_outgroup)
    
    print "\nTransforming counts into centered log ratios and model based methods"
    # Perform centered log-ratio transform 
    ss_T_otus, ss_T_m = standardScaled(ingroup_otu, ingroup_otu_m)
    clr_T_otus, clr_T_m = centeredLogRatio(ingroup_otu, ingroup_otu_m)
    rare_t_otus, rare_t_m = rarefyOTUtables(ingroup_otu, ingroup_otu_m)
    prop_T_otus, prop_T_m = simpleProportions(ingroup_otu, ingroup_otu_m)
    tmm_vst_otus, tmm_vst_m = edgeRtmm(ingroup_otu, ingroup_otu_m)
    rlog_vst_otus, rlog_vst_m = DESeqRlog(ingroup_otu, ingroup_otu_m,
                                          saved_file=True)
    
    print "\n Generating superficial comparisons of transforms"
    
    transforms = {"ss": (ss_T_otus, ss_T_m, 0, 0),
                  "clr": (clr_T_otus, clr_T_m, 0, 1),
                  "rare": (rare_t_otus, rare_t_m, 1, 0),
                  "props": (prop_T_otus, prop_T_m, 1, 1),
                  "tmm": (tmm_vst_otus, tmm_vst_m, 2, 0),
                  "rlog": (rlog_vst_otus, rlog_vst_m, 2, 1)
                  }

    # of dropped OTUS
    f, axarr = plt.subplots(3, 2)
    for t_name, t_dats in transforms.items():
        print "Transformation Name: {}".format(t_name)
        t_otus = t_dats[0]
        t_otu_set = set(t_otus.columns)
        base_otu_set = set(ingroup_otu.columns) 
        otus_removed = list(t_otu_set.symmetric_difference(base_otu_set))
        print "\t# OTUS dropped: {}".format(otus_removed)
    # average change per OTU
        med_otu_change = np.median((t_otu - ingroup_otu).mean().values)
        total_otu_change = (t_otu - ingroup_otu).sum().sum()
        print "\tTotal change in OTUS counts: {}".format(total_otu_change)
        print "\tMedian change in OTUS counts: {}".format(med_otu_change)
    # histogram of sample sizes 
        this_ax = axarr[t_dats[2], t_dats[3]]
        this_ax.set_title(t_name, fontsize=14)
        if (sr != 0) or (sr != 1):
            this_ax.set_xlabel('sample size (bins)', fontsize=14)
        this_ax.set_ylabel('n', fontsize=14)
        this_ax.hist(t_otus.sum(axis=1).values, bins=20)
        for axis in [axarr[sr, sc].xaxis, axarr[sr, sc].yaxis]:
            for tick in axis.get_major_ticks():
                tick.label.set_fontsize(14)
        
    
    print "\Lets see how the signal holds up through PCA:"
    corrComp0, fcomp0 = scalePCAcorrelate(ingroup_otu, ingroup_otu_m, metadata, False)
    corrComp1, fcomp1 = scalePCAcorrelate(ingroup_otu, ingroup_otu_m, metadata, True)
    corrComp2, fcomp2 = scalePCAcorrelate(clr_T_otus, clr_T_m, metadata, True)
    corrComp3, fcomp3 = scalePCAcorrelate(prop_T_otus, prop_T_m, metadata, True)
    corrComp4, fcomp4 = scalePCAcorrelate(tmm_vst_otus, tmm_vst_m, metadata, True)
    corrComp5, fcomp5 = scalePCAcorrelate(rlog_vst_otus, rlog_vst_m, metadata, True)
    corrComp6, fcomp6 = scalePCAcorrelate(rare_t_otus, rare_t_m, metadata, True)
    
    t_types = [prop_T_otus, clr_T_otus, rare_t_otus, tmm_vst_otus]
    mdata_dfs = [prop_T_m, clr_T_m, rare_t_m, tmm_vst_m]
    t_titles = ["Relative Abundances", "Center Log Ratios", "Rarefaction", "TMM (edgeR)"]
    qtup, coltup = ('Q1', 'Q2', 'Q3', 'Q4'), ('blue', 'red', 'green', 'magenta')
    sprow, spcol = [0,0,1,1], [0, 1, 0, 1]
    with plt.style.context('seaborn-whitegrid'):
        f, axarr = plt.subplots(2, 2)
        for tt, title_, ttm, sr, sc in zip(t_types, t_titles, mdata_dfs, sprow, spcol):
            X_std2 = tt.values
            svd = TruncatedSVD(n_components=2, n_iter=20, random_state=42)
            y_std2 = svd.fit_transform(X_std2)
            axarr[sr, sc].set_title(title_, fontsize=14)
            if sr != 0:
                axarr[sr, sc].set_xlabel('PC 1', fontsize=14)
            axarr[sr, sc].set_ylabel('PC 2', fontsize=14)
            for lab, lab_en, col in zip(qtup,range(1,5),coltup):
                axarr[sr, sc].scatter(y_std2[(ttm.Quadrants==lab_en).values, 0],
                                      y_std2[(ttm.Quadrants==lab_en).values, 1],
                                      label=lab,c=col, s=40)
    
            for axis in [axarr[sr, sc].xaxis, axarr[sr, sc].yaxis]:
                for tick in axis.get_major_ticks():
                    tick.label.set_fontsize(14)
        plt.legend(loc='upper right', fontsize=16)


    print "\nPlot Alpha Diversity Heat Map"
    a_metrics = ['chao1', 'enspie']
    fig_nums = [203,333,444]
    for a_metric, fig_n in zip(a_metrics, fig_nums):
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
        plotHeatmap(alpha_plot_df, fig_n)
        al_ax = plt.gca()
        al_ax.set(xlabel='date', ylabel='depth')
    
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
    
    print "\nReassigning replicate numberings to all start at 1"
    final_rep_groups = dropBadReps((primer_outgroup+infrequent_otus), new_rep_groups)
    final_rep_dict = originate_rep_groupings(final_rep_groups)
    for group in final_rep_dict:
        for idx, num in group.items():
            shared_clr_m.ix[idx, 'replicates'] = num
    
    print "\nSplitting design matrix into test and train according to replicate numberings"
    
    idx_of_replicates = list(shared_clr_m[shared_clr_m.replicates != 1].index)
    idx_of_test = list(shared_clr_m[shared_clr_m.replicates == 2].index)
    idx_of_train = list(set(list(shared_clr_m.index)).symmetric_difference(set(idx_of_test)))
    shared_test_x = shared_clr_m.drop(idx_of_train)
    shared_train_x = shared_clr_m.drop(idx_of_replicates)
    shared_test_x.to_csv('shared_test_x.csv')
    shared_train_x.to_csv('shared_train_x.csv')
else:
    print "Skipping ahead"
    shared_test_x = pd.read_csv('shared_test_x.csv', index_col=0)
    shared_train_x = pd.read_csv('shared_train_x.csv', index_col=0)
    

print "\n Pull in all metadata columns"

rate_conc_f = os.path.join(os.getcwd(), 'ChemData', 'mystic_model_data')
obs_conc_f = os.path.join(os.getcwd(), 'ChemData', 'mystic_measured_data')

depthGrad_df, surfaceM_df = LoadChemDFs()
surfaceM_df['date'] = surfaceM_df.index
depthGrad_df.set_index(['date', 'depth'], inplace=True)
    
surfaceM_df['depth'] = np.zeros((surfaceM_df.shape[0], ))
depths_vector = np.arange(22)
surfaceM_df = append_depths(surfaceM_df, depths_vector)
surfaceM_df.set_index(['date', 'depth'], inplace=True)

obs_conc_df = importratesandconcs_obs(obs_conc_f)
rate_dict, conc_dict = importratesandconcs_mod(rate_conc_f)
model_proc_df = time_scale_modeled_chem_data(rate_dict, conc_dict, 146, 
                                             '03/23/2013', '08/15/2013')

print "\n Force pulling any metadata columns out into lists"
remaining_metadata, otu_columns = [],[]
for col_ in shared_test_x.columns:
    if not col_.startswith("seq"):
        remaining_metadata.append(col_)
    else:    
        otu_columns.append(col_)

print "\n Using date + depth multi-index for subsecting x & y vectors"

for design_matrix in [shared_test_x, shared_train_x]:
    design_matrix['Date'] = design_matrix.date.apply(transform_date)
    design_matrix['Depth'] = design_matrix.depth.apply(transform_depth)
    design_matrix.set_index(['Date', 'Depth'], inplace=True)
    design_matrix.drop(remaining_metadata, 1, inplace=True)


## Heavy Hitters
print "\n Lets demo the covariance approach and see how traces follow in tandem"
x_to_print = shared_train_x.copy()
totals = x_to_print.sum(axis=0)
heavy_hitters = totals[totals > 550].index
heavy_hitter_pct = (totals[heavy_hitters].sum() / totals.sum())*100.
heavy_hitter_pct2 = only_depth_otus.T.ix[heavy_hitters, :].sum().sum() 
heavy_hitter_pct2 = heavy_hitter_pct2/only_depth_otus.T.sum().sum()*100.
print "\tThese represent {:.2f}% of the processed data & {:.2f}% of the raw counts".format(heavy_hitter_pct,
                                                                                 heavy_hitter_pct2)
mode1 = ["seq79", "seq6", "seq172"]
mode2 = ["seq24", "seq50", "seq47"]
x_to_print.T.ix[mode1, :].T.plot()
x_to_print.T.ix[mode2, :].T.plot()

from scipy.spatial.distance import pdist, squareform

heavy_x = x_to_print.T.ix[heavy_hitters, :]
test_corr_mat = squareform(pdist(heavy_x, metric='correlation'))
a_df = pd.DataFrame(data=test_corr_mat, 
                    columns = heavy_hitters, 
                    index = heavy_hitters)
plotHeatmap(a_df, 22)
test_corr_mat = np.corrcoef(heavy_x)
a_df = pd.DataFrame(data=test_corr_mat, 
                    columns = heavy_hitters, 
                    index = heavy_hitters)
plotHeatmap(a_df, 23)
test_corr_mat = np.cov(heavy_x)
a_df = pd.DataFrame(data=test_corr_mat, 
                    columns = heavy_hitters, 
                    index = heavy_hitters)
plotHeatmap(a_df, 24)



## Hierarchical agglomorative clustering
# http://sebastianraschka.com/Articles/heatmaps_in_r.html
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram

# perform clustering on rows
test_labels = heavy_x.index
test_corr_df = pd.DataFrame(test_corr_mat, index=test_labels, columns=test_labels)
test_clusters = linkage(test_corr_df.values, method='ward')
test_dendro = dendrogram(test_clusters, labels=test_labels.values)
test_cluster_dict = extract_linkages(test_clusters, test_labels)

# create a sorted index of all the levels of taxanomy
all_taxa = set()
for column_ in x_to_print.columns:
    taxa_str = taxa_series[column_]
    all_taxa.update(taxa_str.split(";"))

print "lets test clustering"
test_clust_df = score_clusters(test_cluster_dict, all_taxa, 
                               taxa_series, test_labels)

test_clust_df.to_csv('test_clust_df.tsv', sep="\t", index_label='Iteration')

# cluster everything
print "Lets cluster all the data"
corr_mat = np.cov(x_to_print.T)
corr_labels = x_to_print.T.index
corr_mat_2A = pdist(x_to_print.values, metric='correlation')
corr_mat_2 = squareform(corr_mat_2A)
row_clusters = linkage(corr_mat, method='ward')
cluster_dict = extract_linkages(row_clusters, corr_labels)
row_dendr = dendrogram(row_clusters, labels=corr_labels, orientation='top')

full_cluster_df = score_clusters(cluster_dict, all_taxa, taxa_series, corr_labels)
full_cluster_df.to_csv('full_cluster_df.tsv', sep="\t", index_label='Iteration')

sys.exit("awaiting your command")

print "\n Lets fit some models"

from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LassoCV as lassoCV
from sklearn.ensemble import RandomForestRegressor as rForestReg

baton_twirler = model_proc_df.copy()
score_fxns = ['r2', 'mse', 'cv-oob']
models = ['lasso', 'rf']
dual_score_fxns = []
for sf in score_fxns:
    for m in models:
        dual_score_fxns.append(m+"-"+sf)
        
score_matrix = np.zeros((len(baton_twirler.columns), 
                         len(dual_score_fxns)))

score_df = pd.DataFrame(data = score_matrix, 
                        index = baton_twirler.columns, 
                        columns = dual_score_fxns)

feat_dict = {}
for baton in baton_twirler.columns:
    y_vector = baton_twirler.ix[:, baton]
    X_train, y_train = matchXandYbyIndex(shared_train_x, y_vector)
    X_test, y_test = matchXandYbyIndex(shared_test_x, y_vector)
    
    lookatTheTrees = rForestReg(n_estimators=500, min_samples_split = 2,
                                n_jobs=1, oob_score=True)
    rodeo_clown = lassoCV(cv=10, max_iter=10000)
    
    lookatTheTrees.fit(X_train, y_train)
    rodeo_clown.fit(X_train, y_train)
    
    y_predicted_rf = lookatTheTrees.predict(X_test)
    y_predicted_lasso = rodeo_clown.predict(X_test)
    
    this_oob = lookatTheTrees.oob_score_
    this_cv = rodeo_clown.mse_path_.min()
    
    rf_mse = mean_squared_error(y_test, y_predicted_rf)
    lasso_mse = mean_squared_error(y_test, y_predicted_lasso)
    
    rf_r2 = r2_score(y_test, y_predicted_rf)
    lasso_r2 = r2_score(y_test, y_predicted_lasso)
    
    feat_dict[baton] = {'rf_estimators': lookatTheTrees.estimators_, 
                        'rf_feat_import':lookatTheTrees.feature_importances_,
                        'lasso_coeffs': rodeo_clown.coef_}
    
    scores = [lasso_r2, rf_r2, lasso_mse, rf_mse, this_cv, this_oob]
    score_df.ix[baton, :] = np.array(scores)
    
score_df.to_csv(os.path.join(os.getcwd(),"model_scores.tsv"), sep="\t")

with open(os.path.join(os.getcwd(), "model_feature_import.p"), "wb") as rf_feat_f:
     pickle.dump( feat_dict, rf_feat_f )

sys.exit()

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




