function plotParam=getplotsettings2(filtLFP,ncschannels,events,samplerate,settingsfile)

    %default settings
    plotParam.plotbad=0;
   
    plotParam.blsub=1;
    plotParam.fontsize=8;
    plotParam.freqlim=[5 100];       %morlet gram
    plotParam.fftclim=[-110 -70];      %power limits color scale log sclae
    plotParam.winlength=.25;         %smooth window for enveloped lfp
    plotParam.winlengthphys=0.3;    %for pulse/lick/eye
    plotParam.dispLFP=1;                %0 do not display LFP, but save, 1 display LFP
    plotParam.displaypH=[0 0 0 0];            %[0 0 1] display pH, BG, or movement or iox
    plotParam.scaleBar=0;     % scaleBar=cmax;
    plotParam.scaleBarConc=10;
    plotParam.colorplotsize=[750 150];
    plotParam.itplotsize=plotParam.colorplotsize;
    plotParam.ncsplotsize=plotParam.colorplotsize;

    plotParam.fftsize=[500 175];
    plotParam.widen=400;
    plotParam.BGstartpoint=20;          %default subtraction sample index
    plotParam.CVavg=2;        
    plotParam.BGavg=2;               %# samples to average
    plotParam.buttonm=0;                %plot morletgram
    plotParam.mcsc='cl1';                %plot morletgram
    plotParam.cmax=1;
    plotParam.zoomTS=[20 40];           %zoom in window to look at LFP TS's
    plotParam.zoomTS=[0 60];
    plotParam.filtbetal=[12 17];        %low beta
    plotParam.filtbetah=[17 33];        %high beta
    plotParam.filtgammal=[40 58];           %low gamma
    plotParam.filtlfp=filtLFP;        %freq bands of LFP filter 1 high beta
    plotParam.envwin=round(samplerate*.5/mean(plotParam.filtlfp));       %envelop over 1/2 cycle of targetd osc 
    plotParam.filtlfp2=[10 17];         %freq band of lfp filter 2 low beta
    plotParam.envwin2=round(samplerate*.5/mean(plotParam.filtlfp2));       %envelop over 1/2 cycle of targetd osc 
    plotParam.filtlick=[0 10];        %freq bands of lick filter
    plotParam.fscvscale=[-5 10];
    plotParam.fscvscale=[];         %concentration scale
    plotParam.LFPscale=[-150 150];              %u-volts
    plotParam.powerscale=[0 6].*1e-3;    %u-volts^2
    plotParam.powerscale=[0 .5].*1e3;    %u-volts^2 p2
    plotParam.eyescale=[-1.5 2.5].*1e3;
    plotParam.physscale=[]; 
    plotParam.unitslfpfiltenv='env-filt LFP \muV^2';
    plotParam.unitsz='env-filt signal z-score';
    plotParam.unitslfp='LFP \muV';
    plotParam.unitslfpfilt='filtered LFP \muV';
    plotParam.colormaptab=[1,0,0;0.223175965665236,1,0;0,0.935622317596566,1;0.609442060085836,0,1;1,0.618025751072961,0;0,1,0.549356223175966;0,0.0343347639484977,1;1,0,0.875536480686695];
    plotParam.colormaptab=[plotParam.colormaptab; plotParam.colormaptab];
        plotParam.colormaptab=[plotParam.colormaptab; plotParam.colormaptab];
    plotParam.colormaptab=plotParam.colormaptab.*.7;
    numch=4;
    getColorCV;
    colorFSCV=colorCV; 
    colorFSCV(2,:)=[0.3 .7 0];
    colorFSCV(3,:)=[0.0 .5 1];
    colorFSCV(3,:)=[0.0 .5 1];
    colorFSCV(4,:)=[0.48  0   1];
    plotParam.colorFSCV=colorFSCV;

    [plotParam.map, plotParam.cmin]=setColorScale(plotParam.cmax);
    plotParam.cscale=[plotParam.cmin plotParam.cmax];
    plotParam.event_codes=events;
    if ~isfield(plotParam,'disp_events')
        %display events
        plotParam.disp_events=[]; idx=1;
        for ii=1:length(events)/2
            plotParam.disp_events(ii)=str2num(events{ii*2-1});
        end
    end
    
    plotParam.xsel=plotParam.BGstartpoint;
    plotParam.ysel=50;
    plotParam.clktype=0;        %mouse clikc type, 0 is left, 1 is right

plotParam.ratelfp=samplerate;
%variables that change upon calling
plotParam.cscNames=ncschannels;      %csc names (labels from map)
plotParam.eyeid=[];
plotParam.pulseid=[];
plotParam.lickid=[];
plotParam.physid=[];
plotParam.lfpid=[];
if ~isempty(ncschannels)
%get ch id's for categories of ncs data
plotParam.eyeid=find(ismember(ncschannels,{'eyed','eyex','eyey'})==1);
plotParam.pulseid=find(ismember(ncschannels,'pulse')==1);
plotParam.lickid=find(ismember(ncschannels,{'lick','lickx','licky','lickz'})==1);
plotParam.physid=[plotParam.eyeid plotParam.pulseid plotParam.lickid];
plotParam.lfpid=find(~ismember(1:length(ncschannels),plotParam.physid));
end
end

