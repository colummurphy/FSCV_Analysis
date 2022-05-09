function plotsummarybeh(xinfo,trialinfo,plotparam,varargin)
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%same as plotx, but plot big/small/targ on same plots
%for plotting multiple sessions
%1/3/2018, udpates as in plotx/plotxallsess
%PLOT TRACES for DA GROUPS, mean
%OVERLAY BIG VS SAMLL and defined groups
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-2 2];     %+/-2 s from aln idx
interval=1;       %in seconds
ratelfp=1000;
argnum=1;
datype='dapos';
evtype={};
sesstypes={};
trtypes={'all'};
binfo=[];
bplot=0;
ntype=[];

%summary of behaviors big vs small, mean from specified b sess id's & conditions ie left



