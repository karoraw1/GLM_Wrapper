# -*- coding: utf-8 -*-
"""
Created on Sun May 29 15:36:20 2016

@author: login
"""

from LakeModel import make_dir
import os, mmap, time
import urllib2
import pandas as pd
import numpy as np
from multiprocessing import Pool
import shutil
from contextlib import closing

def is_archive(extension):
    if extension.upper() == '.HDF':
        return True
    elif extension.upper() == '.HE5':
        return True
    else:
        return False

def is_xml(extension):
    if extension.upper() == '.XML':
        return True
    else:
        return False

def maybe_download(url):
    """Download a file if not present, and make sure it's the right size."""
    filename = os.path.basename(url)
    if filename[-3:] == 'xml' or filename[-3:] == 'he5':
        if not os.path.exists(filename):
            start_dl = time.time()

            with closing(urllib2.urlopen(url)) as r:
                with open(filename, 'wb') as f:
                    shutil.copyfileobj(r, f)
            print time.time()-start_dl
            return filename
        else:
            return "AlreadyDownloaded"
    else:
        print "unexpected file name skipped"
        pass
        return "NotAFile"

def textDict(path, newDict):
    "pass in the path to a text file and a dictionary to-be-filled with lines"
    f = open(path, 'r')
    for idx, line in enumerate(f):
        newDict[idx] = line.strip()
    f.close()
    return newDict
    
def parseXML(xmlfile):
    tags = ['RangeDateTime', 'GranuleID']
    f = open(xmlfile)    
    s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

    ## Initialize variable to return with observed file name
    head1 = s.find('RangeDateTime')
    head2 = s.find('GranuleID')    
    
    mdata = [xmlfile]
    if 'RangeDateTime' in tags:        
        ## Grab time & date info & package into tuple
        s.seek(head1)
        s.readline()
        start_date = s.readline()[24:34]
        start_time = s.readline()[24:32]
        end_date = s.readline()[21:31]
        end_time = s.readline()[21:29]
        start_string = start_date + " " + start_time
        end_string = end_date + " " + end_time
        mdata.append(start_string)
        mdata.append(end_string)
    else:
        mdata.append("")
        mdata.append("")

    if 'GranuleID' in tags:

        ## Grab file name and append to list       
        s.seek(head2+10)
        mdata.append(s.readline()[:-13])
    else:
        mdata.append("")

    s.close()
    f.close()
    return mdata


class GEO_Database(object):
        
    def __init__(self, name, TestBool=True, BuildBool=True):
        self.meta_df = None        
        print "Is this a test run? ", TestBool
        print "Is this a full build? ", BuildBool
        self.test = TestBool
        self.build = BuildBool
        self.meta_db = {}
        self.data_db = {}
        self.name = os.path.join(os.path.split(os.getcwd())[0], 'weatherData', 
                                 name)
        make_dir(self.name)        
        self.mdata_csv = os.path.join(self.name, "meta_table.csv")
        if os.path.exists(self.mdata_csv):
            print "There appears to be an existing metadata table in this location"
            print "To use this intermediate build, use read_metadata_table()`"
            
        #get the test path
        self.test_path = os.path.join(os.path.split(os.getcwd())[0], 
                                 'test_files')        
        
    def read_metadata_URL_list(self, meta_text = None):
        # `meta_text` can be modified by the user to a distinct list of meta
        # data files posted at GEODISC, or another ftp server. The two samples
        # differ considerably in size, to speed up testing. 

        if meta_text is None:
            if self.test is False:
                self.meta_text = os.path.join(self.test_path, 
                                              'GEO_MetaDataOnly.txt')
            elif self.test is True:
                self.meta_text = os.path.join(self.test_path, 
                                              'test_xml_list120.txt')
        else:
            self.meta_text = meta_text
            
        self.meta_db = {}
        self.meta_db = textDict(self.meta_text, self.meta_db)

    
    def testsplit(self):
        
        self.read_metadata_URL_list()
        self.keepers = self.meta_db.values()
        mixed_urls = self.all_urls.items()

        self.keeper_files = set([os.path.basename(i)[:-4] for i in self.keepers])
        self.mixed_files = {os.path.basename(v): (i,v) for i, v in mixed_urls}        
        for kf in self.keeper_files:
            l, u = self.mixed_files[kf]
            self.data_db[l] = u
                        
    def read_combined_URL_list(self, combo_text=None):
        if combo_text is None:
            self.combo_text = os.path.join(self.test_path, 
                                              'GEO_DataAndMetaData.txt')                                              
        else:
            self.combo_text = combo_text
            
        self.all_urls = {}
        self.all_urls = textDict(self.combo_text, self.all_urls)
        old_mdata = self.meta_db.values()

        if self.test is True:
            self.testsplit()
        else:
            for k,v in self.all_urls.items():
                if is_archive(v[-4:]):
                    self.data_db[k] = v
                elif is_xml(v[-4:]) and v[-4:] not in old_mdata:
                    self.meta_db[k] = v
                    
    
    def read_metadata_table(self, meta_csv_path):
        if self.meta_df is None:
            self.meta_df = pd.read_csv(meta_csv_path, index_col=0)
        else:
            print "This one is created already"
    
    def download_data(self, type):
        if type == 'metadata':
            samples = self.meta_db.values()
            print "%r metadata files required" % len(samples)
            fns = {os.path.basename(i):i for i in samples}
            for fn in fns.keys():
                dest = os.path.join(self.name, fn)
                if os.path.exists(dest):
                    samples.remove(fns[fn])
            procs = 10
            print "%r metadata files not located" % len(samples)
                    
        # If no previous build detected, download proceeds
        print "Downloading began at ", time.ctime()
        time_m_dl = time.time()
        
        pool = Pool(processes=procs)
        self.results = pool.map(maybe_download, samples)
        print "Downloading ended after %r s." % (time.time()-time_m_dl)
        s_right = [os.path.join(self.name, os.path.basename(i)) for i in samples]
        self.s_wrong = [os.path.join(os.getcwd(), os.path.basename(i)) for i in samples]
        for idx, s in enumerate(self.s_wrong):
            if os.path.exists(s):
                shutil.move(s, s_right[idx])
        self.meta_db_files = s_right
    
    def check_table_to_database(self, table_choice):
        ## Remove records without a file name
        clean_clouds = cloud_df[cloud_df.hdf5.notnull()]
        
        ## Verify that metadata records are matched to existing archive URLs
        for f_n in list(clean_clouds.hdf5):
            if f_n in all_files.keys():
                data_urls.append(all_files[f_n])
            else:
                bad_urls.append(all_files[f_n])
        
        ## If any metadata records are not located, this assertion fails
        archives = len(data_urls)
        try:
            assert archives == clean_clouds.shape[0]
        except AssertionError:
            print archives, " archives"
            print clean_clouds.shape[0], "metadata recs"
        assert archives == clean_clouds.shape[0]
        
    
    def write_metadata_table(self, tags = ['RangeDateTime', 'GranuleID']):    
        if len(self.meta_db_files) == 0:
            for f in self.meta_db.values():
                shouldbe = os.path.join(self.name, os.path.basename(f))    
                if os.path.exists(shouldbe):
                    self.meta_db_files.append(shouldbe)
                    
        
        pool = Pool(processes=len(self.meta_db_files))
        self.timeIDlist = pool.map(parseXML, self.meta_db_files)
        row_n = len(self.timeIDlist)
        col_n = len(self.timeIDlist[0])
        cols = ['Path', 'Start', 'End', 'file']
        cloud_data = np.array(self.timeIDlist).reshape((row_n,col_n))
        self.meta_df = pd.DataFrame(index = self.meta_db_files, columns = cols,
                                    data = cloud_data)
        self.out_csv = os.path.join(self.name, "meta_db.csv")    
        self.meta_df.to_csv(self.out_csv)