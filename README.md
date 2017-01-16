## GLM Wrapper

This is a pipeline from various data sources into the General Lake Model among many other things.

The following is a manifest containing explanations of the file structure above. 

1. **GLM_Executables/** : Contains all binaries & source code for GLM as well as **example configurations**
2. **GLM_documentation/** : Contains two official GLM informational documents, the writup and presentation I did for the class, and useful papers
3. **MARCCTEST/** : Contains scripts for performance testing on MARCC with diff. combos of nodes, GLM function calls, and partitions. Lessons learned include:
  * Queue times for each partition are highly variable depending on other users demands, but usually one is free enough to run immediately
  * Parallel automatically assigns 24 CPUS per node, while others assign 1
  * GLM runs >1000x faster on Mac OS compared to Linux
4. **OTU_Time_Series/** : contains OTU data for modelled year, see README within that folder for more info
5. **bin/** : contains scripts for rapid execution of multiple cases of GLM for calibration purposes, each script should/will contain explanations within.
6. **donePlots/** : contains many plots used for GLM writeup and other generated during testing
7. **test_files/** : contains many files used to test GLM wrapper functionality (including the all important "optimize_these.txt")
8. **waterData/** : contains all USGS discharge, salinity, rainfall data for many different streams. Only some are used
9. **weatherData/** : contains meterological data from CERES and USGS Hydrological Service
10. **morphometry/** : contains base image used for deriving morphometry data

The files contained within this root directory include descriptions of the work created at various stages in development.
These require revision and coallation. 



