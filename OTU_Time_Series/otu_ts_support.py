#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:15:28 2016

@author: login
"""
import pandas as pd
import numpy as np
import sys, os
from sklearn.decomposition import PCA, TruncatedSVD, LatentDirichletAllocation
from sklearn.preprocessing import StandardScaler

def removeZeroCols(df):
    return (df.T[(df != 0).any()]).T

def parseBiosample(df):
    """
    This function accepts a df with samples in rows and OTUs in columns
    and parses the biosample key 
    """
    biosamples = list(df.index)
    dates_, primers_, kits_, replicates_, depths_ = [], [], [], [], []
    for bs in biosamples:
        if bs[:2] == "SB":
            clipped = bs[2:]
        else:
            sys.exit("Non-SB start to biosample")
        
        if 'TAWMD' in clipped:
            date, rest = clipped.split("TAWMD")
        elif "TAWWD" in clipped:
            date, rest = clipped.split("TAWWD")
        else:
            sys.exit("Bridge Sequence not Detected")
            
        if "VV4T" in rest:
            primer = "VV4"
            depth, rest2 = rest.split("VV4T")
        elif "VV4V5T" in rest:
            primer = "V4V5"
            depth, rest2 = rest.split("VV4V5T")
        else:
            sys.exit("Primer Sequence not Detected")
        
        if rest2[0] == "M":
            kit = "MolBio"
        elif rest2[0] == "Q":
            kit = "Qiagen"
        elif rest2[:4] == "filt" and rest2[4] == "M":
            kit = "MolBio"
        elif rest2[:2] == "NA":
            kit = "NA"
        else:
            print clipped
            print rest2[0]
            sys.exit("kit type not detected")
        
        if rest2[-2] == "R":
            replicate = rest2[-1]
        else:
            sys.exit("replicate signifier not detected")
        
        if depth == '015':
            depth = '01.5'
        
        dates_.append(date)
        primers_.append(primer)
        kits_.append(kit)
        replicates_.append(replicate)
        depths_.append(depth)
    
    df['date'] = dates_
    df['primers'] = primers_
    df['kit'] = kits_
    df['replicates'] = replicates_
    df['depth'] = depths_
    return df

def add_quadrants(dF):
    depth_list = list(dF.depth)
    depth_ns = np.array([float(n) for n in depth_list])
    quad_ns = np.array(["  "]*len(depth_ns))
    quad_ns[depth_ns < 5] = "Q1"
    quad_ns[(depth_ns > 4) & (depth_ns < 11)]= "Q2"
    quad_ns[(depth_ns > 10) & (depth_ns < 17)]= "Q3"
    quad_ns[depth_ns > 16]= "Q4"
    dF["Quadrants"] = quad_ns
    return dF.copy()
    
def numericEncodings(df, metadata_cols, verbose=True):
    print "Changing {} metadata columns".format(len(metadata_cols))
    unq_metadata = {i:{} for i in metadata_cols}
    for md, unqs in unq_metadata.items():
        print "\n", md
        for num, unq in enumerate(np.unique(df[md])):
            num+=1
            unqs[num] = unq
            if verbose == True:
                print "Encoding", unq, "as", num
            bool_ = df[md] == unq 
            df.ix[bool_, md] = num

    for i in metadata_cols:
        df[[i]] = df[[i]].apply(pd.to_numeric)
        
    df2 = df.copy()
    return df2, unq_metadata
    
from itertools import groupby

def listRepGroups(df):
    # make a list of the index 
    dfindex = list(df.index)
    # create a list for each grouping and a list of already matched samples
    rep_groups, consumed = [], []
    # Start with one index
    for first in dfindex:
        this_group = []
        # Check the entire list for members that match 
        if first not in consumed:
            for second in dfindex:
                # If a sample wasn't already consumed, and isn't the exact same sample
                if second not in consumed and first != second:
                    # check for a match
                    if first.split("VV")[0] == second.split("VV")[0]:
                        # if detected, add to the already consumed list 
                        consumed.append(second)
                        this_group.append(second)
        if len(this_group) != 0:
            this_group.append(first)
            this_group.sort()
            rep_groups.append(this_group)
    rep_groups.sort()
    unq_rep_groups = list(rg for rg,_ in groupby(rep_groups))
    print len(unq_rep_groups), "groups of replicates detected"
    return unq_rep_groups

def JensenShannonDiv_Sqrt(df_otu):
    ps_df = df_otu.copy().replace(0.0, 1.0)
    ps_n_df = ps_df.divide(ps_df.sum(axis=1), axis=0)
    shape_sq = len(ps_n_df.index)    
    dist_dat = np.zeros((shape_sq, shape_sq))
    
    for r_idx, r in enumerate(ps_n_df.index):
        for c_idx, c in enumerate(ps_n_df.index):
            x_ = ps_n_df.ix[r, :]
            y_ = ps_n_df.ix[c, :]
            m_ = (x_+y_) / 2
            jsd_sqrt = (0.5*(x_*np.log2(x_/m_) + y_*np.log2(y_/m_)).sum())**0.5
            dist_dat[r_idx, c_idx] = jsd_sqrt

    dist_mat = pd.DataFrame(index=ps_n_df.index, columns=ps_n_df.index,
                            data=dist_dat)
    return dist_mat
    
def ReplicateReport(df, df_otus, rep_groups, verbose=True, metric="JSD"):
    print "REPLICATE REPORT"
    
    in_rep_distances, all_dists, worst_reps = [], [], []
    broken_groups = 0
    
    for idx, group in enumerate(rep_groups):
        print "Group {}".format(idx)
        
        this_grps = df_otus.ix[group, :]
        dist_mat = JensenShannonDiv_Sqrt(this_grps)
        for a_d in np.unique(dist_mat.values):
            if a_d != 0:
                in_rep_distances.append(a_d)
        
        if verbose == True:
            print dist_mat.max()
            
        most_distant = dist_mat.max().max()
        
        print "Most Distant: {}".format(most_distant)
        if most_distant > 0.3:
            broken_groups+=1
            while most_distant > 0.3:
                # find one of the pair that are most divergent
                bad_reps_bool = (dist_mat.max() == dist_mat.max().max())
                bad_means = dist_mat[bad_reps_bool].mean(axis=1)
                worst_rep = bad_means.argmax()
                worst_reps.append(worst_rep)
                print "\tdropping {}".format(worst_rep)
                group.remove(worst_rep)
                this_grps = df_otus.ix[group, :]
                dist_mat = JensenShannonDiv_Sqrt(this_grps)          
                most_distant = dist_mat.max().max()
                print "\tmost distant now: {}".format(most_distant)
    
    in_rep_distances = np.array(in_rep_distances)
    if len(in_rep_distances) > (len(df_otus.index)-2):
        random_samples = len(df_otus.index)-2
    else:
        random_samples = len(in_rep_distances)
        
    rand_samp_idx = list(np.random.choice(df_otus.index, random_samples, 
                                          replace=False))
    
    print "{} random vectors selected".format(len(rand_samp_idx))
    df_otus_rand = df_otus.ix[rand_samp_idx, :]
    rand_dist_mat = JensenShannonDiv_Sqrt(df_otus_rand)
    for all_d in np.unique(rand_dist_mat.values):
        if all_d != 0:
            all_dists.append(all_d)
        
    all_dists = np.array(all_dists)
    for lab, dist_arr in zip(['In Rep','Rand Dists'], [in_rep_distances, all_dists]):
        print lab
        print "Mean: {}".format(dist_arr.mean())
        print "Variance: {}".format(np.var(dist_arr))
        print "CV: {}".format(np.std(dist_arr)/dist_arr.mean())
        print "n: {}".format(len(dist_arr))
     
    
    return worst_reps, broken_groups
    
"""
        this_sum = this_grps        
        best_rep = list(this_sum[ this_sum == this_sum.max()].index)
        bad_reps = list(this_sum[this_sum != this_sum.max()].index)        
        farthest_rep = dist_mat.max().max()

        if farthest_rep > 0.25:
            while farthest_rep > 0.25:
                distances = euclidean_distances(df.ix[best_rep, :], 
                                                df.ix[bad_reps, :])
                worst_rep = bad_reps[np.argmax(distances)]
                bad_reps.remove(worst_rep)
                worst_reps.append(worst_rep)
                rem_reps = bad_reps + best_rep
                rem_grps = df.ix[rem_reps, 'enspie']
                rem_sum = rem_grps
                this_cv = (rem_sum.std() / rem_sum.mean())
        else:
            rem_reps = group
            
        this_grps = df_otus.ix[rem_reps, :]
        this_bool = this_grps != 0
        this_sum = this_bool.T.sum()
        
        if verbose == True:
            print group[0][:-1]
            print "Mean Count Total:", this_sum.mean()
            print "CV:", (this_sum.std() / this_sum.mean())
            reps, _ = this_grps.shape
            n_ = range(1,reps+1)
            totalTaxa = (this_bool.sum() > 0).sum()
            print "Total Taxa Discovered in {} reps: {}".format(reps ,totalTaxa)
            for n in n_:
                overlap = ((this_bool.sum() == n).sum()/float(totalTaxa))*100.0
                print "{0:.3f}% among {1} replicates".format(overlap, n)
            
    return (list(set(worst_reps)), cvs)
"""

def dropBadReps(less_diverse_reps, rep_groups):
    new_rep_groups = []
    for g in rep_groups:
        for l in less_diverse_reps:
            if l in g:
                g.remove(l)
            else:
                pass
        if len(g) > 1:
            new_rep_groups.append(g)
        else:
            pass
        
    return new_rep_groups
    
import matplotlib.pyplot as plt
import seaborn as sns
    
def plotHeatmap(df, fignum):
    plt.figure(fignum, figsize=(12,9))
    ax = sns.heatmap(df)
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    
    for item in ax.get_xticklabels():
        item.set_rotation(90)
    #plt.savefig('seabornPandas.png', dpi=100)
    plt.show()
    
    
import ecopy as ep

def beta_wrapper(df, var_key):
    print "\n", var_key
    brayDist = ep.distance(df, method='bray')
    groups = list(df[var_key])
    rand_groups = list(np.random.choice(np.unique(groups), 
                                    size=np.array(groups).shape))

    ep.beta_dispersion(brayDist, groups, test='anova', 
                   center='median', scores=False)
    ep.beta_dispersion(brayDist, rand_groups, test='anova', 
                   center='median', scores=False)
    return brayDist
    
def scalePCAcorrelate(df_numerical, df_w_mdata, metadata_cols, transformed):
    if transformed:
        X_std2 = df_numerical.values.T
    else:
        X_std2 = StandardScaler().fit_transform(df_numerical.values.T)
    
    rows_n, cols_n = X_std2.shape
    print "\nPerforming PCA"
    print "Matrix has {} samples and {} features".format(rows_n, cols_n)
    pca2 = TruncatedSVD(n_components=100, n_iter=20, random_state=42)
    pca2.fit(X_std2)
    print "Top two components explain {} and {} of variance.".format(pca2.explained_variance_ratio_[0],
                                                                     pca2.explained_variance_ratio_[1])
    for mdata in metadata_cols:
        md_arr = np.array(df_w_mdata[mdata])
        
        corrs = [np.corrcoef(pca2.components_[i, :], md_arr)[0, 1] for i in range(100)]
        print "Metadata Variable:", mdata
        print "Max Correlation:", np.array(corrs).max()
        pca_comp_no = np.argmax(np.array(corrs))
        print "PCA component No:", pca_comp_no
        print "Explained Variance:", pca2.explained_variance_ratio_[pca_comp_no]
        print ""
    return pca2
        
def readChemData(chem_path, units, ftype, plotbool=False):
    print os.path.basename(chem_path)
    print units
    
    if ftype == 'depth_profile':
        site_chem_spec_df = pd.read_csv(chem_path, index_col=0, 
                                        parse_dates=True, 
                                        infer_datetime_format=True)
        new_idx = []
        for i in site_chem_spec_df.index:
            new_idx.append(pd.Period(i, 'M'))
        site_chem_spec_df.index = new_idx
        
    elif ftype == 'surface_measure':
        chem_spec_csv = pd.read_csv(chem_path, sep=",", index_col=0)
        print "Null Values: {}".format(chem_spec_csv.isnull().sum().sum())
        print "Database Shape: {}".format(chem_spec_csv.shape)
        new_cols = [pd.Period(i) for i in chem_spec_csv.columns]
        chem_spec_csv.columns = new_cols
        site_chem_spec_df = chem_spec_csv.T.interpolate()
    else:
        sys.exit("invalid ftype")
        
    if plotbool == True:
        site_chem_spec_df.plot()
    return site_chem_spec_df
    
def plotCommonTaxa(taxa_series):
    
    taxa_list = list(taxa_series.values)
    taxa_decomposition = []
    for tax_str in taxa_list:
        this_tl = tax_str.split(";")
        clean_list = [i for i in this_tl if i[-2:] != "__" ]
        taxa_decomposition += clean_list
    unq_taxa = np.unique(taxa_decomposition, return_counts=True)
    dfindex, dfdata = unq_taxa[0], unq_taxa[1]
    taxa_df = pd.DataFrame(data=dfdata, index=dfindex, columns=['count'])
    taxa_df.sort_values('count', ascending=False, inplace=True)
    rtaxa_df = taxa_df.divide(taxa_df.sum(), axis=1)
    total_called = []
    t_levels_s = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    t_levels = ['kingdom', 'pylum', 'class', 'order', 'family', 'genus', 'species']
    rtaxa_df = taxa_df.divide(taxa_df.sum(), axis=1)
    for tL in t_levels_s:
        posHits = [i for i in rtaxa_df.index if i[0] == tL]
        subdf = rtaxa_df.ix[posHits, :]
        print tL, subdf.sum()
        total_called.append(subdf.sum())
    t_level_rel = np.array(total_called)
    width = 0.35
    ind = np.arange(len(t_levels))
    print "total pct%", t_level_rel.sum()
    fig, ax = plt.subplots(1, 1)
    ax.bar(ind + width, t_level_rel,width)
    ax.set_xticks(ind + width)
    ax.set_xticklabels(t_levels)
    
def inferTaxaLevel(taxa_series):
    addSemicolon = lambda x: x+";"
    taxa_edit = taxa_series.copy().apply(addSemicolon)
    taxa_dict = taxa_edit.to_dict()
    taxa_frame = taxa_edit.copy().to_frame()
    taxa_depth = np.zeros(taxa_series.shape)
    taxa_frame["Taxa Depth"] = taxa_depth
    for seq, ts in taxa_dict.items():
        taxa_frame.ix[seq, "Taxa Depth"] = 7 - ts.count("_;")
    return taxa_frame

    
def analyze_alpha_diversity(decoder, derepped_otus_m, valPairs):
    ad_cats =  []
    ad_cols = {"Mean":[], "Median":[], "Std":[], "N":[]}
                
    for var in decoder.keys():
        codes = decoder[var].keys()
        for code in codes:
            if var == 'date':
                ad_cats.append(str(decoder[var][code]).split("T")[0]+" ("+var+")")
            elif var == 'replicates':
                pass
            else:
                ad_cats.append(str(decoder[var][code])+" ("+var+")")
            
            if var == 'replicates':
                pass
            else:
                sub_bool = derepped_otus_m[var] == code
                subdf = derepped_otus_m[sub_bool]
                ad_cols["Median"].append(np.median(subdf.enspie))
                ad_cols["Std"].append(subdf.enspie.std())
                ad_cols["N"].append(subdf.shape[0])
                ad_cols["Mean"].append(subdf.enspie.mean())
    
    for idx, vp in enumerate(valPairs):
        kitT, primT = vp[0], vp[1]
        bool1 = derepped_otus_m.primers == primT
        bool2 = derepped_otus_m.kit == kitT    
        subdf2 = derepped_otus_m[bool1 & bool2]
        if idx == 2:
            primer_outgroup = list(subdf2.index)
            
        ad_cats.append(str(decoder['primers'][primT])+" & "+str(decoder['kit'][kitT]))
        ad_cols["Median"].append(np.median(subdf2.enspie))
        ad_cols["Std"].append(subdf2.enspie.std())
        ad_cols["N"].append(subdf2.shape[0])
        ad_cols["Mean"].append(subdf2.enspie.mean())
                
    alpha_df = pd.DataFrame(data=ad_cols, index=ad_cats)
    return alpha_df, primer_outgroup
    
def alpha_diversity(dereplicated_otus, derepped_otus_m, metrics):
    row_sum = dereplicated_otus.copy().sum(axis=1)
    row_rel = dereplicated_otus.copy().divide(row_sum, axis=0).astype('float64')
    if 'enspie' in metrics:
        enspie_sq = row_rel.apply(np.square)
        enspie_dom = enspie_sq.sum(axis=1)
        enspie_ = enspie_dom**-1
        derepped_otus_m['enspie'] = enspie_
    if 'shannon' in metrics:
        entrop = lambda n: n*np.log(n)
        shannon_ = row_rel.replace({ 0 : np.nan }).applymap(entrop).T.sum()*-1.0
        derepped_otus_m['shannon'] = shannon_.apply(np.exp)
    if 'chao1' in metrics:
        total_s =  (dereplicated_otus > 0).T.sum()
        singletons = (dereplicated_otus == 1).T.sum()
        doubletons = (dereplicated_otus == 2).T.sum()
        numerator = singletons.multiply(singletons-1)
        denominator = 2*(doubletons+1)
        chao1_ = total_s + numerator.divide(denominator, axis=0)
        derepped_otus_m['chao1'] = chao1_
    
    return dereplicated_otus, derepped_otus_m
    
    
def plotCountTotalsByMetadata(df_m, decoder, mList, mVar, fignum): 
    ## Get Normalized Count Totals
    # drop everything but the key grouping variable and sum
    mList.remove(mVar)
    counts_Depth = df_m.drop(mList, axis=1)
    depthGroup = counts_Depth.groupby([mVar]).sum()
    # find the number of samples per grouping
    (_, n_per_depth) = np.unique(df_m[mVar].values, return_counts=True)
    # average the total counts per group by the number of samples
    mean_counts_per_group = depthGroup.T.sum().divide(n_per_depth)
    ## Get Standard Deviation of Count Totals
    # Drop depths & sum each sample before grouping
    just_counts = counts_Depth.drop([mVar], axis=1)
    depthSum = just_counts.T.sum().to_frame()
    # Convert Series into DataFrame, add col names, and modify dtype
    depthSum[mVar] = df_m[mVar].values
    depthSum.columns = ['counts', mVar]
    depthSum = depthSum.applymap(pd.to_numeric)
    # group each sum by depth and flatten by std
    depthStd = depthSum.groupby([mVar]).std()
    # convert labels for display
    if mVar == 'date':        
        decoded_labs = [str(decoder[mVar][i]).split("T")[0] for i in list(np.unique(df_m[mVar].values))]
    else:
        decoded_labs = [str(decoder[mVar][i]) for i in list(np.unique(df_m[mVar].values))]
    # Setup Plot
    width = 0.35; ind = np.arange(len(n_per_depth));
    plt.figure(fignum, figsize=(14, 8))
    ax = plt.gca()
    ax.bar(ind + width, mean_counts_per_group, width,
            yerr=depthStd.values.flatten())
    ax.set_xticks(ind + width)
    ax.set_xticklabels(decoded_labs)
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45)
    ax.set_xlim(0.0, float(len(list(np.unique(df_m[mVar].values)))))
    ax.set_xlabel(mVar.capitalize()+" (m)")
    ax.set_ylabel("Average Total OTUs (n)")
    # Add back metadata variable
    mList.append(mVar)
    
def replicateAlphaDiv(df, metric, rep_groups):
    enspie_cv = []
    enspie_1 = []
    enspie_2 = []
    for g in rep_groups:
        this_grp = df.ix[g, metric]
        enspie_cv.append(this_grp.std() / this_grp.mean())
        justTwo = list(np.random.choice(g, size=2, replace=False))
        enspie_1.append(df.ix[justTwo[0], metric])
        enspie_2.append(df.ix[justTwo[1], metric])

    enspie_1 = np.array(enspie_1)
    enspie_2 = np.array(enspie_2)

    return (enspie_cv, enspie_1, enspie_2)

def centeredLogRatio(otu_table, otu_table_m=None):
    from scipy.stats.mstats import gmean
    noZeros = otu_table.copy().replace(0, np.nan)
    geomeans = np.repeat(np.nan, repeats = noZeros.shape[0])
    for i in range(0, noZeros.shape[0]):
        geomeans[i] = gmean(noZeros.ix[i, :].dropna())
    clr_table = np.log(noZeros.divide(geomeans, axis=0))
    clr_table.replace(np.nan, 0, inplace=True)
    if otu_table_m is not None:
        clr_table_m = otu_table_m.copy()
        clr_table_m.ix[:, otu_table.columns] = clr_table
        return clr_table, clr_table_m
    else:
        return clr_table

from sklearn.model_selection import train_test_split

def lda_tuner(ingroup_otu, best_models):

    best_score = -1*np.inf
    dtp_series = [0.0001, 0.001, 0.01, 0.1, 0.2]
    twp_series = [0.0001, 0.001, 0.01, 0.1, 0.2]
    topic_series = [3]
    X = ingroup_otu.values
    eval_counter = 0

    for topics in topic_series: 
        for dtp in dtp_series:
            for twp in twp_series:
                eval_counter +=1
                X_train, X_test = train_test_split(X, test_size=0.5)
                lda = LatentDirichletAllocation(n_topics=topics, 
                                                doc_topic_prior=dtp, 
                                                topic_word_prior=twp, 
                                                learning_method='batch',
                                                random_state=42,
                                                max_iter=20)
                lda.fit(X_train)
                this_score = lda.score(X_test)
                this_perplexity = lda.perplexity(X_test)
                if this_score > best_score:
                    best_score = this_score
                    print "New Max Likelihood: {}".format(best_score)

                print "#{}: n:{}, dtp:{}, twp:{}, score:{}, perp:{}".format(eval_counter, 
                                                                 topics, dtp, twp,
                                                                 this_score, this_perplexity)

                best_models.append({'n': topics, 'dtp': dtp, 'twp': twp,
                                    'score': this_score, 'perp': this_perplexity})
                if (dtp == dtp_series[-1]) and (twp == twp_series[-1]):
                    eval_counter +=1
                    X_train, X_test = train_test_split(X, test_size=0.5)
                    lda = LatentDirichletAllocation(n_topics=topics, 
                                                    doc_topic_prior=1./topics, 
                                                    topic_word_prior=1./topics, 
                                                    learning_method='batch',
                                                    random_state=42,
                                                    max_iter=20)
                    lda.fit(X_train)
                    this_score = lda.score(X_test)
                    this_perplexity = lda.perplexity(X_test)
                    if this_score > best_score:
                        best_score = this_score
                        print "New Max Likelihood: {}".format(best_score)

                    print "#{}: n:{}, dtp:{}, twp:{}, score:{} perp: {}".format(eval_counter, 
                                                                                topics, 
                                                                                (1./topics), 
                                                                                (1./topics),
                                                                                this_score,
                                                                                this_perplexity)

                    best_models.append({'n': topics, 'dtp': (1./topics), 
                                        'twp': (1./topics), 'score': this_score,
                                        'perp': this_perplexity})
    return best_models
    

    
def collapseBdiversity(dist_mat, raw_data_m, metaData_var):
    metaOptions = np.unique(raw_data_m.ix[:, metaData_var])
    n_ = len(metaOptions)
    metaDistance = np.full((n_, n_), np.nan)
    metaDeviation = np.full((n_, n_), np.nan)
    for r_i, r in enumerate(metaOptions):
        for c_i, c in enumerate(metaOptions):
            dist_copy = dist_mat.copy()
            dist_copy[metaData_var] = raw_data_m.ix[:, metaData_var]
            dist_filt_1 = dist_copy[dist_copy[metaData_var] == r]
            dist_filt_1.drop([metaData_var], axis=1, inplace=True)
            dist_filt_t = dist_filt_1.T
            dist_filt_t[metaData_var] = raw_data_m.ix[:, metaData_var]
            dist_filt_2 = dist_filt_t[dist_filt_t[metaData_var] == c]
            dist_filt = dist_filt_2.drop([metaData_var], axis=1)
            dist_filt_flat = dist_filt.values.flatten()
            dist_filt_nz = dist_filt_flat[dist_filt_flat != 0]
            mD = np.median(dist_filt_nz)
            mDev = dist_filt_nz.std()
            print "{} versus {} metadistance".format(r, c)
            print "\t {} ({})".format(mD, mDev)
            metaDistance[r_i,c_i] = mD
            metaDeviation[r_i,c_i] = mDev
    return metaDistance, metaDeviation

"""    
http://qiime.org/scripts/make_otu_network.html
http://qiime.org/scripts/differential_abundance.html

# alternative to random forest classifier 

from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

clf = GaussianProcessClassifier(1.0 * RBF(1.0), warm_start=True,
                                n_jobs=-1)
X = StandardScaler().fit_transform(X)
X_train, X_test, y_train, y_test = \
        train_test_split(X, y, test_size=.4, random_state=42)
clf.fit(X_train, y_train)
score = clf.score(X_test, y_test)

"""

import subprocess as sp

def varStabTransform(path, df_otus):
    """
    Uses edgeR's variance stabilizing transformation to transform sequence counts
    into OTU abundances. Involves adding and subtracting psuedocounts. 
    """
    wrapper_file = os.path.join(os.getcwd(), "edgeRwrapper.R")
    if not os.path.exists(wrapper_file):
        sys.exit("Accessory script missing")
        
    # direct export path
    to_transform = path
    shared_otus = df_otus.copy()
    arg_1 = os.path.dirname(to_transform)
    arg_2 = to_transform.split("/")[-1]
    # export data to disk
    shared_otus.to_csv(to_transform, sep=",")
    # run R script 
    base_cmd = "Rscript edgeRwrapper.R"
    cmd = " ".join([base_cmd, arg_1, arg_2])
    p = sp.Popen(cmd, cwd=os.getcwd(), shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
    stdout, stderr = p.communicate()
    if "Execution halted" in stderr:
        sys.exit("R wrapper failed")
    
    to_read_back = os.path.join(arg_1, arg_2.split(".")[0]+"_vst.csv")
    shared_vst = pd.read_csv(to_read_back, index_col = 0)
    for i, j in zip(shared_vst.index, shared_otus.index):
        assert i == j
    for i, j in zip(shared_vst.columns, shared_otus.columns):
        assert i == j 
    
    for temps in [to_read_back, to_transform]:
        os.remove(temps)
    return shared_vst