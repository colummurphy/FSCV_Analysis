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
