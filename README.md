# analysis
Main updated code in V2

getTRdata

  Extract trial by trial data with DA all together, and Nlx/ephys separated by channel
  Combined good/bad and all trial conditions into one file for each channel
  Refer to trlist in trlists for trial by trial info
  Outputs trlists and tr_XX mat files for each Nlx channel and FSCV
  *NEED to update the folders in the top of the code*
  e.g. getTRdata(83)
  
  
updateTRlists

  Same as above, but just update trlist in trlists file and resave in order to update manually selected good/bad trials list (saved in excel in original FSCV dir)
  
  
analyzeSes

  Open up trlists (as created in above functions) and create workspace for analyzing data within
  e.g. analyzeSes(83,'fscvchs',[2 4],'nlxch',{'cl4-cl6','pulse','lickx'})
  *NEED to update the folders in the top of the code*
  
Procedure for use:
Extract Session Data Scripts:
First change chronicXXchconfigsimple.m (in MATLAB analysis folder) to define ncschannels to extract and the correct path for paths{1} (can delete all current variables used to define previous paths{1} and do not need paths{2}) (this would be the path containing original FSCV data, e.g. Y:\data_MIT\patra_fscv\patra_chronic94_06252018)  Update file
Run getTRdata for target session, e.g. getTRdata(94)  This will store analysis extracted trial by trial data into the analysis folder, e.g.  C:\Users\putamen\Documents\MATLAB\analysis\chronic94
After data has been extracted can load analysis scripts
Session Analysis Scripts:
Run analyzeSes , e.g. analyzeSes(83,'fscvchs',[2 4],'nlxch',{'eyed','pulse','lickx’}) to load trlists into directory containing all targeted data variables for analysis
Raster plots, trial by trial:
Run tr_raster(….)
![image](https://user-images.githubusercontent.com/23349223/170110037-1685cfca-733c-4704-930a-ac3420d7e133.png)

