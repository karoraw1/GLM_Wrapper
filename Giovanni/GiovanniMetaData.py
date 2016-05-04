# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:41:42 2016

This downloads xml metadata files from a text list 
It reads the start, end, and corresponding data file name.
It writes this data out to a CSV

@author: Keith
"""
import subprocess
import pandas as pd
import numpy as np
import os
import time
import mmap
from joblib import Parallel, delayed


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
    f = open(xmlfile)
    s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    parsed = []
    head1 = s.find('RangeDateTime')
    head2 = s.find('GranuleID')
    print "head1: %r" %head1
    s.seek(head1)
    s.readline()
    start_date = s.readline()[24:34]
    start_time = s.readline()[24:32]
    end_date = s.readline()[21:31]
    end_time = s.readline()[21:29]
    start_string = start_date+" "+start_time
    end_string = end_date+" "+end_time
    parsed.append((start_string, end_string))
    print "head2: %r" %head2
    s.seek(head2+10)
    parsed.append((s.readline().split('<')[0]))
    s.close()
    f.close()
    return parsed

def dlANDparseXML(meta_urls, out_df, i):
    m_u = meta_urls[i]
    """pass this a metadata url and it downloads and pulls out fname & dates"""
    if (m_u[0:3] == 'ftp') and (m_u[:-4:-1] == 'lmx'):
        start_dl = time.time()
        print "Start Download: \n\t%s" % (m_u)
        #xmlfile, headers = urllib.urlretrieve(m_u)
        wget = subprocess.Popen("wget -nv " + m_u, shell = True, executable = '/bin/bash')
        wget.communicate()
        #open .xml & get [(filename, (start,_end)]
        xmlfile = os.path.basename(m_u)
        mdata = getmetadata(xmlfile)        
        #print "\t Nulls in `cloud_recs` before", i,": ", out_df.isnull().sum()
        #out_df.loc[i]['start'] = mdata[0][0]
        #out_df.loc[i]['end'] = mdata[0][1]
        #out_df.loc[i]['hdf5'] = mdata[1]
        #print out_df.loc[i]
        #print "\t Nulls in `cloud_recs` After", i,": ", out_df.shape
        subprocess.call('rm '+xmlfile, shell = True, executable = '/bin/bash')
        print "Deleting temp file. This file took %r" % (time.time()-start_dl)
        return mdata
    else:
        print "dlANDparseXML error: not an xml file"


test = False

# start with data folder & point to needed files
giovanni_path = '/home/login/GDrive/Documents/Landscape_Hydrology/Final Project/Giovanni'

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

if __name__ == '__main__':
    cloud_recs = Parallel(n_jobs=-1)(delayed(dlANDparseXML)(meta_urls, cloud_recs, i) for i in samples)

# pull out individual lists from list of lists
times, files = zip(*cloud_recs)
starts, ends = pd.datetools.parse_time_string(np.array(zip(*times)))

# Put each list into numpy matrix and then dataframe
row_n = starts.shape[0]
cloud_data = np.hstack((starts, ends, np.asarray(files))).reshape((3, row_n))
cloud_df = pd.DataFrame(index = samples, data = cloud_data.T, 
                          columns=['start', 'end', 'hdf5'])

# print parsed tags to a csv file
if test == True:
    out_csv = os.path.join(giovanni_path, "cloud_meta.csv")
else:
    out_csv = os.path.join(giovanni_path, "test_meta.csv")

cloud_df.to_csv(out_csv)