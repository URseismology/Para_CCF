# Para_CCF
parallel the computing of CCF in Matlab

## Run the code
run get_BH_CCF.m in matlab

## Code description
### get_BH_CCF.m
1. need to use Matlab2019a, the detail of conneting to clusters can be found at: http://www.sas.rochester.edu/ees/urseismo/projects-urseismo/

2. input for get_BH_CCF.m should be a txt file with columns [net1, sta1, comp1, net2, sta2, comp2]

3. get_BH_CCF.m parallel the jobs by parallely calls the function a1_ccf_ambnoise_RTZ_NE_Para.m

### a1_ccf_ambnoise_RTZ_NE_Para.m
1. aa1_ccf_ambnoise_RTZ_NE_Para.m is a matlab function, and it takes input: sta1, comp1, net1, sta2, comp2, net2, PATHS
 
2. need to set the workingdir, datadir, and log_path (log_path is for debug purpose)
 
3. the datapath should be the path to the Data folder. The structure of the Data folder: Data/network/station/actual_data
 
4. set IsOutputFullstack = 1 to get the output of yearly-stack CCF in matlab format
  
5. set IsFigure1 = 0 to get plots of the CCF result (not recommended to do this in paralle job; it might crush matlab)

6. set IsSaveTxt = 1 to get the text files needed for the future study (zero_crossing, fitting Bessel function...) using AkiEstimate
