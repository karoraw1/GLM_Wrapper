# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 10:58:24 2016

@author: login
"""
import matplotlib.pyplot as plt

results = {}

class MARCCRun(object):
    def __init__(self, job_num):
        self.job_num = job_num
        self.nodes = None
        self.calls = None
        self.partition = ''
        self.queueTime = None
        self.runTime = None
        self.completed = None
        self.lineR = []
        self.lineQ = []
    


with open("MARCC_Test_Results.txt", "r") as f:
    for line in f:
        
        splitLine = line.split()
        splitLine.remove('SLURM')
        
        
        job_ = splitLine[0].split('=')[1]
        splitLine[1] = splitLine[1].split('=')[1]
        splitLine[1] = splitLine[1].split('.')[0]
        
        if job_ not in results.keys():
            results[job_] = MARCCRun(job_)
            results[job_].nodes = int(splitLine[1].split('_')[1])
            results[job_].calls = int(splitLine[1].split('_')[3])
            results[job_].partition = splitLine[1].split('_')[5]
            if results[job_].partition == 'parallel':
                results[job_].nodes*=24
        
        therest = splitLine[2:]
        
        if 'Ended,' in therest:
            results[job_].lineR = therest
            rT_noF = therest[3][:-1].split(':')
            rT_F = int(rT_noF[0])*60*60 + int(rT_noF[1])*60 + int(rT_noF[2])
            results[job_].runTime = rT_F
            
        if 'TIMEOUT' in therest:
            results[job_].completed = False
        elif 'COMPLETED,' in therest:
            results[job_].completed = True
            
        if 'Began,' in therest:
            results[job_].lineQ = therest
            qT_noF = therest[3].split(':')
            qT_F = int(qT_noF[0])*60*60 + int(qT_noF[1])*60 + int(qT_noF[2])
            results[job_].queueTime = qT_F

import pandas as pd

cols = ['nodes', 'calls', 'partition', 'queueTime', 'runTime', 'completed']

mTest_df = pd.DataFrame(index=results.keys(), columns=cols)

for job in results.keys():
    run = results[job]
    mTest_df.loc[job]['nodes'] = run.nodes
    mTest_df.loc[job]['calls'] = run.calls
    mTest_df.loc[job]['partition'] = run.partition
    mTest_df.loc[job]['queueTime'] = run.queueTime
    mTest_df.loc[job]['runTime'] = run.runTime
    mTest_df.loc[job]['completed'] = run.completed
    
mTest_df['tot_time_per_executions'] = (mTest_df.runTime+mTest_df.queueTime*1.0)/mTest_df.calls
mTest_df['run_time_per_executions'] = (mTest_df.runTime*1.0)/mTest_df.calls
success=mTest_df[mTest_df.completed == True]
success.sort_values('run_time_per_executions', inplace=True)

parallel_jobs = mTest_df[mTest_df.partition == 'parallel']

#plt.scatter(parallel_jobs.nodes, parallel_jobs.queueTime/60)
#plt.ylabel('Time in job queue (min)')
#plt.xlabel('CPUs requested (24 per node)')
#plt.title('Queue Time on Parallel Partition')

shared_success = success[success.partition != 'parallel']
parallel_success = success[success.partition == 'parallel']

plt.figure(2)
plt.scatter(shared_success.nodes, shared_success.runTime/60./shared_success.calls, 
            c='r', label='shared')
plt.scatter(parallel_success.nodes, parallel_success.runTime/60./parallel_success.calls, 
            c='b', label='parallel')
plt.xlabel('CPUS')
plt.ylabel('Run Time per Call (min)')
plt.legend(loc='best', frameon=True, shadow=True)
