%record settings for fscv & artifact removal
diffChannel=16;
resort=[];
badChannels=[];

parameters.calcdata=1;          %for calculating processed.info parameters
parameters.traditionalQthres=1;     %use calculated Qa as thres, not standrad deviation like with monkey
parameters.samplesperscan=214;
parameters.glitchThres=50;       %

parameters.badoverlap=[20 25];      %-/+ # of fscv samples around artifact to devalidate for peaks
%ie
%time from badoverlap(1) - artifact is invalid (or aka maxTS+badoverlap(1))
%time from artifact to + badoverlap(2) is invalid (aka maxTS-badoverlap(2))
parameters.numpcs=2;            %1/2019 change back to 3 pcs rather than 2 /////09/06/2018 just DA & pH components for plotting concentrations to remove 

parameters.glitchThresT=0.15;       %high frequency departure along time rather than CV
parameters.RThresOut=0.85;              %if signal correlated by > this to movement template remove

parameters.MThres=0.1;             %new threshold 10/03/2017, automatically remove if movement PCR component greater than
parameters.RThres=0.8;             %smaller R for scaled K version works better
parameters.RThres=0.75;             %smaller R for scaled K version works better

parameters.peakDiff=2.5;             %nM difference must be reached to distinguish peak
parameters.peakthreshold=1; %# of std's for threshold of peak, new detectdatransients.m 05/28/2018
parameters.peakthreshold=5; %# of std's for threshold of peak, new detectdatransients.m 05/28/2018
parameters.peakthreshold=2; %# of std's for threshold of peak, 09/25/2018
parameters.peakthreshold=1; %# of std's for threshold of peak, 09/25/2018

parameters.peakDiff=5;             %nM difference must be reached to distinguish peak, alternative to std

parameters.minPeakWidth=3;     %peak must reach < max-peakdif*std above this many samples, other wise transient "glitch" 
parameters.maxDurationPeak=20;         %time to find half life, in either direction
parameters.maxDurationPeak=25;         %limit on what defines transient 05/28/2018
parameters.maxDurationPeak=50;         %limit on what defines transient 05/28/2018
parameters.maxDurationPeak=100;         %changed 10/31/2018, see longer transeitns eg chronic 67

parameters.tracepad=[-10 10];        %seconds before/after transient window for storing traces


parameters.chunksize=30;           %300 samplse per chunk
parameters.peakWindow=10;         %number of samples +/- around tscanID to check for peak peak detected to check for half 
parameters.meanBGThres=.1;          %nA change in mean BG current between samples otherwise don't use as BG
parameters.deflectionThres=1;
parameters.stableLevel=0.05;
parameters.scanWindowArtifact=10;
parameters.scanWindowMaxDetected=7;
parameters.noiselevel=0.1;
parameters.saturationLevel=8;      %peak current to accept
parameters.window_range_samples=200;       % samples +/- to check around BG sub point for correlated CV's
parameters.noiseThres=0.1;
parameters.BGWin=[100 100];               % samples +/- to check around BG sub point for correlated CV's
parameters.noiseIThres=0.5;             %noise threshold in retrieved samples
parameters.noiseDAThres=60;             %da noise threshold in retrieved samples


parameters.includepcam=0;           %05/16/2018, do not include movement PCA component in current contributions calculation