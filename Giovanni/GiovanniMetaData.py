# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:41:42 2016

This downloads xml metadata files from a text list 
It reads the start, end, and corresponding data file name.
It writes this data out to a CSV
If passed a command line argument, a test will be performed on the first
16 files

@author: Keith
"""
import subprocess
import pandas as pd
import numpy as np
import os, sys
import time
import mmap
from joblib import Parallel, delayed
from Giovanni_SpatialFilter import dl_parse_hdf


def textDict(path, fname):
    f_path = os.path.join(path, fname)
    f = open(f_path, 'r')
    newDict = {}
    for idx, line in enumerate(f):
        newDict[idx] = line.strip()
    f.close()
    return newDict

def getmetadata(xmlfile):
    # this parses the metadata xml file in one pass and return a list
    try:
        f = open(xmlfile)
    except IOError:
        xmlfile = xmlfile + ".1"
        f = open(xmlfile)
        
    s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    parsed = [xmlfile]
    head1 = s.find('RangeDateTime')
    head2 = s.find('GranuleID')
    
    # Grab time & date info
    s.seek(head1)
    s.readline()
    start_date = s.readline()[24:34]
    start_time = s.readline()[24:32]
    end_date = s.readline()[21:31]
    end_time = s.readline()[21:29]
    start_string = start_date+" "+start_time
    end_string = end_date+" "+end_time
    parsed.append((start_string, end_string))
    
    # Grab file name 
    s.seek(head2+10)
    parsed.append((s.readline().split('<')[0]))
    s.close()
    f.close()
    return parsed

def dlANDparseXML(meta_urls, out_df, i, test):
    m_u = meta_urls[i]
    """pass this a metadata url and it downloads and pulls out fname & dates"""
    if (m_u[0:3] == 'ftp') and (m_u[:-4:-1] == 'lmx'):

        if test == True:
            start_dl = time.time()
            print "Start Download: \n\t%s" % (m_u)

        wget = subprocess.Popen("wget -nv " + m_u, shell = True, executable = '/bin/bash')
        wget.communicate()
        #open .xml & get [(filename, (start,_end)]
        xmlfile = os.path.basename(m_u)
        mdata = getmetadata(xmlfile)
        xmlfile = mdata.pop(0)            
        mdata.insert(0, i)
        subprocess.call('rm '+xmlfile, shell = True, executable = '/bin/bash')

        if test == True:
            print "Deleting temp file. This file took %r" % (time.time()-start_dl)
            
        return mdata
    else:
        print "dlANDparseXML error: not an xml file"
        mdata = []
        mdata.append(i)
        mdata.append((np.nan, np.nan))
        mdata.append(np.nan)
        return mdata
        
if len(sys.argv) > 1:
    test = True
else:
    test = False

print "Is this a test run? ", test

# start with data folder & point to needed files
giovanni_path = os.getcwd()

# both .he5 & xml files
all_text = 'ListOfCloudFiles_Giovanni.txt'

# just metadata (xml) files
meta_text = 'GiovanniMetadataOnlyList.txt'
meta_test = 'test_xml_list.txt'

# create handles

if test == False:
    meta_urls = textDict(giovanni_path, meta_text)
else:
    meta_urls = textDict(giovanni_path, meta_test)

# Print the number of files to be downloaded
samples = meta_urls.keys()
print "%r metadata files located" % len(samples)

# loop through list of all urls
all_urls = textDict(giovanni_path, all_text)

# pre-allocate space for parallel jobs
cloud_recs = []
time_m_dl = time.time()

# Try to find previous evidence of download
if test == False:
    m_csv_p = os.path.join(giovanni_path, "cloud_meta.csv")
else:
    m_csv_p = os.path.join(giovanni_path, "test_meta.csv")

make_db = False

try:
    cloud_df = pd.read_csv(m_csv_p, index_col=0)
    print "Skipping metadata download. Database parse detected"
except IOError:
    make_db = True
    print "No previous metadata detected. Downloading metadata"

# If no previous build detected, download proceeds
if make_db == True:
    print "Downloading began at ", time.ctime()
    if __name__ == '__main__':
        accumulator = []
        n_iter = 0
        with Parallel(n_jobs=-1) as parallel:
            cloud_recs = parallel(delayed(dlANDparseXML)(meta_urls, cloud_recs, i, test) for i in samples)
            [accumulator.append(l) for l in cloud_recs]
            n_iter += 1
            
    print "Downloading ended after %r s." % (time.time()-time_m_dl)
    print "Total of %r iterations" % n_iter
    print "Length of Accumulator: %r" % len(accumulator)
    
    # pull out individual lists from list of lists
    
    if len(cloud_recs) != len(meta_urls.keys()):
        print "Not enough records returned by joblib"
        sys.exit()        
        
    idx_out, times, files = zip(*cloud_recs)
            
    starts, ends = pd.datetools.parse_time_string(np.array(zip(*times)))
    
    # Put each list into numpy matrix and then dataframe
    row_n = starts.shape[0]
    cloud_data = np.hstack((starts, ends, np.asarray(files), 
                            np.asarray(idx_out))).reshape((4, row_n))
                            
    cloud_df = pd.DataFrame(index = samples, data = cloud_data.T, 
                              columns=['start', 'end', 'hdf5', 'idx_out'])
    
    # print parsed tags to a csv file
    if test == False:
        out_csv = os.path.join(giovanni_path, "cloud_meta.csv")
    else:
        out_csv = os.path.join(giovanni_path, "test_meta.csv")
    
    cloud_df.to_csv(out_csv)

# Else we move ahead to next set of downloads
else:
    clean_clouds = cloud_df[cloud_df.hdf5.notnull()]
    data_urls = []
    all_files = {os.path.basename(k):k for k in all_urls.values()}
    bad_urls = []
    
    for f_n in list(clean_clouds.hdf5):
        if f_n in all_files.keys():
            data_urls.append(all_files[f_n])
        else:
            bad_urls.append(all_files[f_n])
            
    archives = len(data_urls)
    try:
        assert archives == clean_clouds.shape[0]
    except AssertionError:
        print archives, " archives"
        print clean_clouds.shape[0], "metadata recs"
        assert archives == clean_clouds.shape[0]
    
    data_vals = []
    time_d_dl = time.time()
    print "Data download began at ", time.ctime()
    if __name__ == '__main__':
        accumulator2 = []
        n_iter2 = 0
        with Parallel(n_jobs=-1) as parallel:
            data_vals = parallel(delayed(dl_parse_hdf)(data_urls, data_vals, j) for j in range(archives))
            [accumulator2.append(l) for l in data_vals]
            n_iter2 += 1
        
    print "Data download ended after %r s." % (time.time()-time_d_dl)
    print "Total of %r iterations" % n_iter2
    print "Length of Accumulator: %r" % len(accumulator2)
    
    i_out, cc_means, cc_stds, cc_n, lat_stds, lons_stds, checksum, min_lim = np.array(zip(*data_vals))
    new_cols = [i_out, cc_means, cc_stds, cc_n, lat_stds, lons_stds, checksum, min_lim]
    new_labels = ["i_out_2", "cc_means", "cc_stds", "cc_n", "lat_stds", 
                  "lons_stds", "checksum", "min_lim"]
    for colu, labe in zip(new_cols, new_labels):
        clean_clouds[labe] = colu

    cleanup = clean_clouds.checksum != '.he5'
    
    if cleanup.sum() != 0:
        cleanuplist = list(clean_clouds.hdf5[cleanup].values)
        for leftover in cleanuplist:
            os.remove(leftover)
    
    if clean_clouds.columns[0][:7] == 'Unnamed':
        clean_clouds.drop(clean_clouds.columns[0], axis=1)
        
    
    if test == False:
        d_csv_p = os.path.join(giovanni_path, "cloud_data_adaptive_filter.csv")
        print "cloud fetch complete"
    else:
        d_csv_p = os.path.join(giovanni_path, "test_data.csv")
        print "testing complete"
    
    clean_clouds.to_csv(d_csv_p)
    
    
    
