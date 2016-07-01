# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 18:28:28 2016

@author: login
"""



if plotting == True:
    old_cols = Hobbs.df.columns.values
    old_cols[7] = 'Water Temperature'
    old_cols[15] = 'Air Temperature'
    old_cols[9] = 'Specific conductance'
    old_cols[17] = 'Precipitation'
    
    Hobbs.df.columns = old_cols
    to_drop = list(old_cols)
    to_drop.remove(old_cols[7])
    to_drop.remove(old_cols[15])
    to_drop.remove(old_cols[9])
    to_drop.remove(old_cols[17])
    temp_df = Hobbs.df.drop(to_drop, axis=1)
    temp_df=temp_df.fillna(method="pad")
    temp_ser3 = temp_df['Precipitation'].resample('D').sum()
    temp_df = temp_df.resample('D').mean()
    temp_df2 = LakeModel.subsectbydate_1(FreshPond.df, temp_df)
    temp_df2=temp_df2.fillna(method="pad")
    
    air = LakeModel.error_metrics(temp_df2['Temperature, air'].values,
                        temp_df['Air Temperature'].values)
    cond = LakeModel.error_metrics(temp_df2['Specific conductance'].values,
                         temp_df['Specific conductance'].values)
    water = LakeModel.error_metrics(temp_df2['Temperature, water'].values,
                          temp_df['Water Temperature'].values)
    precip = LakeModel.error_metrics(temp_df2['Precipitation, total'].values,
                           temp_ser3.values)


    plt.style.use('ggplot')
    f, ((ax3, ax2), (ax1, ax4)) = plt.subplots(2, 2, sharex='col', figsize=(24,12))
    ax4.plot(temp_ser3)
    ax4.plot(temp_df2['Precipitation, total'])
    ax4.set_title('Precipitation', fontsize=20)
    ax3.plot(temp_df['Air Temperature'])
    ax3.plot(temp_df2['Temperature, air'])
    ax2.plot(temp_df['Water Temperature'])
    ax2.plot(temp_df2['Temperature, water'])
    ax1.plot(temp_df['Specific conductance'])
    ax1.plot(temp_df2['Specific conductance'])
    ax1.set_title('Specific conductance', fontsize=20)
    ax2.set_title('Water Temperature', fontsize=20)
    ax3.set_title('Air Temperature', fontsize=20)
    ax4.set_alpha(0.5)
    
    
    ax1.text(733238, 1700, "KSE: "+"{:04.2f}".format(cond['KGE']), fontsize=18)
    ax1.text(733230.5, 1600, "NSE: "+"{:04.2f}".format(cond['NSE']), fontsize=18)
    ax1.set_ylabel('microSiemens/cm @ 25C', fontsize=20)
    
    ax2.text(732953, 2.5, "KSE: "+"{:04.2f}".format(water['KGE']), fontsize=18)
    ax2.text(732953, 0.5, "NSE: "+"{:04.2f}".format(water['NSE']), fontsize=18)
    ax2.set_ylabel('degrees C', fontsize=20)
    
    ax3.text(733238, -12, "KSE: "+"{:04.2f}".format(air['KGE']), fontsize=18)
    ax3.text(733238, -14.5, "NSE: "+"{:04.2f}".format(air['NSE']), fontsize=18)
    ax3.set_ylabel('degrees C', fontsize=20)
    
    ax4.text(732953, 2.8, "KSE: "+"{:04.2f}".format(precip['KGE']), fontsize=18)
    ax4.text(732953, 2.62, "NSE: "+"{:04.2f}".format(precip['NSE']), fontsize=18)
    ax4.set_ylabel('inches (sum)', fontsize=20)
    f.legend( ax1.get_lines(), ['Hobbs Brook', 'Fresh Pond'], loc = (0.405, .465), 
             fontsize=18, ncol=2)
    
    ## Saving is broken for some reason. Copying and pasting the above and below
    ## code into ipython sequentially works though
    #labs = ax4.get_xticklabels()
    #ax1.set_xticklabels(labs, fontsize=12, color='k')
    #ax4.set_xticklabels(labs, fontsize=12, color='k')
    #f.savefig("HobbsandFresh.png")


#test.mergedata()

if plotting == True:
    plt.style.use('fivethirtyeight')
    adaptive_n = pd.read_csv(adapt_data, index_col = 0, 
                             parse_dates = ['start', 'end'], 
                             infer_datetime_format = True)
    adaptive_3n = pd.read_csv(adapt_data3n, index_col = 0, 
                              parse_dates = ['start', 'end'], 
                              infer_datetime_format = True)
    
    adaptive = pd.concat([adaptive_3n, adaptive_n], ignore_index=True)
    flat_n = test.geodisc_clouds.sample(n=160)
    xlabs = ["Std. Deviation per Swath(n)", "Measurements per Swath (n)" ]
    col_names = ['cc_stds', 'cc_n']
    fignames = ['adaptiveVflat_stds.png', 'adaptiveVflat_n.png']
    colors1 = ['r', 'g']
    
    for xlab, col_name, fig_n, cs in zip(xlabs, col_names, fignames, colors1):
        s1 = flat_n[col_name].values
        s2 = adaptive[col_name].values
        plot_vals = pd.Series(np.concatenate((s1,s2)))
        labs1 = np.tile(np.array('0.3 Degree Threshold'), (160,))
        labs2 = np.tile(np.array('Non-Zero Minimum Threshold'), (160,))
        plot_grps = np.concatenate((labs1, labs2))
        ax1, ax2 = plot_vals.hist(by=plot_grps, figsize=(12,6), color=cs)
        ax1.set_xlabel(xlab)
        ax2.set_xlabel(xlab)
        ax1.set_ylabel("Swath Count")
        plt.savefig(fig_n, dpi=300)


def getStartTimes(df, col1_i, i):
    a_time  = df.iloc[i, col1_i]
    total_s = a_time.hour * 3600 + a_time.minute * 60 + a_time.second
    return total_s/3600.

def getPeriodicity(df, col1_i, time_d, i):
    date_format = '%Y-%m-%d %H:%M:%S'
    if i < time_d:
        return np.nan
    else:
        end_s = str(test.geodisc_clouds.iloc[i, col1_i])
        end_t = pd.datetime.strptime(end_s, date_format)
        start_s = str(test.geodisc_clouds.iloc[i-time_d, col1_i])
        start_t = pd.datetime.strptime(start_s, date_format)
    return (end_t - start_t).seconds/60./60.


if plotting == True:
    nondataCols = ['start', 'end', 'hdf5', 'i_out_2', 'checksum']
    cc_plt_df = test.geodisc_clouds.drop(nondataCols, axis=1)
    cc_plt_df.columns = ['Mean Cloud Cover', 'Cloud Cover STD',
                         'Cloud Cover Sample Size (n)', 'Latitude Std (deg)', 
                         'Longitude Std (deg)' ]
    
    rows, cols = cc_plt_df.shape
    cc_plt_df.insert(0, 'Time of Day (hrs)', np.zeros(rows))
    
    for n in range(rows):
        cc_plt_df.iloc[n, 0] = getStartTimes(test.geodisc_clouds, 1, n)
    
    cc_plt_df.hist(bins=30)

sorted_df = test.geodisc_clouds.sort_values('start', axis = 0, ascending = True)
Geodisc_clouds = sorted_df.copy()
Geodisc_clouds.set_index('start', drop=True, inplace=True, verify_integrity=True)


lon2df = CERES_data.df2.cldarea_total_daily.fillna(method='pad')
lon_dfs = [lon1df,lon2df ]
lon_keys = list(CERES_data.lons_n)
Geo_filled=Geodisc_clouds.interpolate()
for srs, key in zip(lon_dfs, lon_keys):
    cld_C, cld_G = LakeModel.subsectbydate_2(srs,Geo_filled.cc_means)
    cld_1 = cld_G.resample('D').mean()
    cld_1 = cld_1.interpolate()*100.
    errors[key] = LakeModel.error_metrics(cld_1, cld_C)

lon1df_t, _ = LakeModel.subsectbydate_2(lon1df,Geo_filled.cc_means)
lon2df_t, _ = LakeModel.subsectbydate_2(lon2df,Geo_filled.cc_means)
errors['other'] = LakeModel.error_metrics(lon1df_t, lon2df_t)


SSF_D = Geodisc_clouds.resample('D').mean()
cld_GE, cld_SSF = LakeModel.subsectbydate_2(Geo_filled,SSF_Filled.fillna(method='bfill'))
cld_GE, cld_SSF = LakeModel.subsectbydate_2(Geo_filled,SSF_Filled.fillna(method='bfill'))
errors['SSFvGEO'] = LakeModel.error_metrics(cld_SSF.cloud_frac.values, 
                                  (cld_GE.cc_means*cld_SSF.cloud_frac.max()).values)
lon_keys2 = ['df1', 'df2']
for srs, key in zip(lon_dfs, lon_keys2):
    cld_C, cld_G = LakeModel.subsectbydate_2(srs,filledagain.cloud_frac)
    cld_1 = cld_G.resample('D').mean()
    cld_1 = cld_1.interpolate()*cld_C.max()
    errors[key] = LakeModel.error_metrics(cld_1, cld_C)

if plotting == True:    
    f , ax5 = plt.subplots(1, 1, figsize=(12,12))
    ax5.plot(cld_SSF.cloud_frac, label="CERES SSF Data")
    ax5.plot((cld_GE.cc_means*cld_SSF.cloud_frac.max()), label="Goddard OMI Data")
    ax5.legend()
    ax5.set_ylim((-10, 100))
    ax5.set_xlim(734503+90,734503+185)
    f.savefig("CERES_SSFvOMMYCLD0.png")

    f, axarr = plt.subplots(5, 1, sharex='all', figsize=(24,12))
    axarr[0].plot(inflow.df.index, np.ones(inflow.df.index.shape[0]), linewidth=35)
    axarr[0].set_title('Inflow Params', fontsize=22)
    axarr[1].plot(SSF_Filled.index, np.ones(SSF_Filled.shape[0]), linewidth=35)
    axarr[1].plot(SSF_Filled.index, np.ones(SSF_Filled.shape[0]), linewidth=35)
    axarr[1].set_title('Humidity', fontsize=22)
    axarr[2].plot(BOS_weather.df.index, np.ones(BOS_weather.df.index.shape[0]), linewidth=35)
    axarr[2].plot(BOS_weather.df.index, np.ones(BOS_weather.df.index.shape[0]), linewidth=35)
    axarr[2].plot(BOS_weather.df.index, np.ones(BOS_weather.df.index.shape[0]), linewidth=35)
    axarr[2].set_title('Wind, Precipitation, Snow, Temperature', fontsize=22)
    axarr[3].plot(lon1df.index, np.ones(lon1df.index.shape[0]), linewidth=35)
    axarr[3].plot(lon1df.index, np.ones(lon1df.index.shape[0]), linewidth=35)
    axarr[3].plot(lon1df.index, np.ones(lon1df.index.shape[0]), linewidth=35)
    axarr[3].plot(lon1df.index, np.ones(lon1df.index.shape[0]), linewidth=35)
    axarr[3].set_title('Cloud Area Fraction, SW & LW Radiation', fontsize=22)
    axarr[4].plot(mystic.df.index, np.ones(mystic.df.index.shape[0]), linewidth=35)
    axarr[4].plot(mystic.df.index, np.ones(mystic.df.index.shape[0]), linewidth=35)
    axarr[4].plot(mystic.df.index, np.ones(mystic.df.index.shape[0]), linewidth=35)
    axarr[4].plot(mystic.df.index, np.ones(mystic.df.index.shape[0]), linewidth=35)
    axarr[4].plot(mystic.df.index, np.ones(mystic.df.index.shape[0]), linewidth=35)
    axarr[4].set_title('Outflow Params', fontsize=22)
    f.savefig("TimeCoveredByDataset.jpg")


if plotting == True:    
    sorted_df.insert(0, 'periodicity (hours)', np.zeros(rows))
    
    for m in range(rows):
        sorted_df.iloc[m, 0] = getPeriodicity(test.geodisc_clouds, 1, 1, m)
    
    sorted_df.set_index('start', drop=True, inplace=True, verify_integrity=True)
    sorted_df_pre_agg, new_idx_cols = LakeModel.TimeIdx(sorted_df)
    null_col = sorted_df_pre_agg.notnull().cc_means
    sorted_df_pre_agg['Null Values'] = null_col
    cc_aggs_mu = LakeModel.makeAggregations(sorted_df_pre_agg, new_idx_cols, np.mean)
    cc_aggs_std = LakeModel.makeAggregations(sorted_df_pre_agg, new_idx_cols, np.std )
    cc_aggs_sum = LakeModel.makeAggregations(sorted_df_pre_agg, new_idx_cols, np.sum)
    
    plt.figure(3)
    plt.clf()
    plt.scatter(flat_n.cc_stds , flat_n.cc_n, c = 'r')
    plt.scatter(adaptive.cc_stds , adaptive.cc_n, c = 'b')
    plt.ylabel('Measurements (n)')
    plt.xlabel('Standard Deviation (n)')