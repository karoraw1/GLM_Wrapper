# -*- coding: utf-8 -*-
"""
Created on Tue May 31 16:56:03 2016

@author: login
"""

from GEO_Database import GEO_Database

mystic_lake = GEO_Database(name="test_db", TestBool=True, BuildBool=True)
mystic_lake.read_combined_URL_list()
mystic_lake.download_data(type="metadata")
mystic_lake.write_metadata_table(tags = ['RangeDateTime', 'GranuleID'])
mystic_lake.read_metadata_table(mystic_lake.out_csv)

other_lake = GEO_Database(name="test_db", TestBool=False, BuildBool=False)
other_lake.read_metadata_table(mystic_lake.out_csv)

