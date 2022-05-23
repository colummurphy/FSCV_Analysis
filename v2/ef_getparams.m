function [parameters,csc_map,eventcodes]=ef_getparams(subjectid,recordid,ncschannels,varargin)
%05/02/2022 Merge param with plotparam from getplotsettings.m
%get plotParam and parameters variables as well as maps for subject &
%recording setting
rateEphys=1e3;  %assume all 1kHz downconverted ep hys
csc_map={};
eventcodes={};
parameters.ratelfp=rateEphys;  %default sampling rate
%default fscv params
anodal_lim=1.3;
cathodal_lim=-0.4;
vlim=anodal_lim-cathodal_lim;
Vox=0.6;
sampling_freq=10;   %10 hz sampling freq 
offset=0;
argnum=1;
sessionnum=[];
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'sessnum'
            argnum=argnum+1;
            sessionnum=varargin{argnum};
    end
    argnum=argnum+1;
end
if isempty(sessionnum)
    %delete variable, so does not get fed into patra_map_biplar to det
        %which csc_map to use
        clear sessionnum;
end
%get csc_map, PCR file dir, eventcodes for subject
switch subjectid
    case 'patra'
        patra_map
    case 'patrabipolar'
        patra_map_bipolar
    case 'patrabichunky'
        patra_map_bipolar_chunky
    case 'cleo'
        cleo_map_bipolar
    case 'cleo16'
        cleo_map_bipolar16
    case 'cleo17'
        cleo_map_bipolar17
    case 'cleo25'
        cleo_map_bipolar25
    case 'cleo14'
        cleo_map_bipolar14   
    case 'cleo18'
        cleo_map_bipolar18
    case 'cleo19b'
        cleo_map_bipolar19b
    case 'cleo19'
        cleo_map_bipolar19
    case 'cleo20'
        cleo_map_bipolar20
    case 'cleo22'
        cleo_map_bipolar22
    case 'cleo24'
        cleo_map_bipolar24
    case 'cleo28'
        cleo_map_bipolar28    
    case 'cleo15'
        cleo_map_bipolar15
    case 'cleo23'
        cleo_map_bipolar23    
    case 'cleo22b'
        cleo_map_bipolar22b
    case 'cleo19b'
        cleo_map_bipolar19b
    case 'cfmea'
        rodent_map
    otherwise
        %no ncs channels, maybe da recording only
end

if exist('recordid')
    if ~isempty(recordid)
%get settings for recording session, for artifact removal
%load parameters variables
switch recordid
    case 'default'
        settings_default
    case 'noisycleo'
        settings_cleo_noisy
    case 'cleo'
        settings_cleo_noisy
        if exist(['settings_' subjectid])>0
            run(['settings_' subjectid]);
        end
    case 'cleo14'
        settings_cleo14
    case 'rodent'
        settings_rodent
        parameters.differentialChannel=diffChannel;
        parameters.badChannels =badChannels;        
        parameters.resort=resort;
        parameters.IVconv=1e9/4.99e6; 
end
end
else
    settings_default
end

if exist('ncschannels')
    if ~isempty(ncschannels)
    %targeted ncs to display
    %get cscid's in the order they are given in ncschannels
    for ii=1:length(ncschannels)
        cscid=find(ismember(csc_map,ncschannels{ii}))-1;
        parameters.NCSchannels(ii)=str2num(csc_map{cscid});
    end
    parameters.CSCnames=ncschannels;    %Store CSC labels, replaces orig plotParam.CSCnames
    else
        parameters.NCSchannels=[];
    end
end
    
%PCR templates, subject-specific depending on referencing method
%MAKE SEPARATE FUNCTION HERE
%paramsFSCV=getPCAparams(PCR_srcdir,PCR_template);
parameters.templatesrcdir=PCR_srcdir;   %from subjectid loaded map file
parameters.PCR_template=PCR_template; 

%default fscv paramters based on templates
load([parameters.templatesrcdir parameters.PCR_template])      %loads Ats_mat & Cts_mat
[index1,~]=find(Cts_mat~=0); index_da=find(index1==1); index_ph=find(index1==2);
Cts_mat_DAPH=Cts_mat(1:2,1:index_ph(end)); CV_mat_DA=Ats_mat(:,index_da); CV_mat_pH=Ats_mat(:,index_ph);

parameters.DAcvs=CV_mat_DA;
parameters.Cts=Cts_mat;
parameters.CV=Ats_mat;

if strcmp(recordid,'rodent')
    %scale pca templates
    [index1,index2]=find(Cts_mat~=0); index_da=find(index1==1); index_ph=find(index1==2);
    Cts_mat_DA=Cts_mat(1,index_da); Cts_mat_pH=Cts_mat(2,index_ph); Cts_mat_DAPH=Cts_mat(1:2,1:index_ph(end));
    
    resampled=resampleCV(Ats_mat,length(Ats_mat),parameters.samplesperscan);       %resample templates to match recorded length
    Ats_mat=resampled;
    CV_mat_DA=Ats_mat(:,index_da);
    CV_mat_pH=Ats_mat(:,index_ph);
    %calculate PCR parameters from templates
    [Vc F Aproj Qa]=generatePCR(Ats_mat,Cts_mat);
    [VcDAPH FDAPH AprojDAPH QaDAPH]=generatePCR([CV_mat_DA CV_mat_pH],Cts_mat_DAPH);
    QaDAPH=Qa*5;            %just for CFMEA
    K=pinv(FDAPH*transpose(VcDAPH));        %na/uM, ideal CV's for each analyte
    parameters.DAcvs=CV_mat_DA;
    parameters.Cts=Cts_mat_DAPH;
    parameters.CV=[CV_mat_DA CV_mat_pH];
end

%calculate PCR parameters from templates
[Vc, F, Aproj, Qa]=generatePCR(Ats_mat,Cts_mat);
[VcDAPH, FDAPH, AprojDAPH, QaDAPH]=generatePCR([CV_mat_DA CV_mat_pH],Cts_mat_DAPH);
K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte

parameters.Vc=Vc;  parameters.F=F;    parameters.Aproj=Aproj;    
parameters.Qa=Qa;  parameters.QaDAPH=QaDAPH;  parameters.K=K;

%default info about fscv data based on templates
sizeData=size(Ats_mat); 
samplesperscan=sizeData(1); 
stepsizeCV=((vlim-vlim/samplesperscan)*2/samplesperscan);
Vrange=(cathodal_lim:stepsizeCV:(anodal_lim-vlim*2/samplesperscan))+offset;
Vrange_cathodal=(anodal_lim:-stepsizeCV:(cathodal_lim+vlim*2/samplesperscan))+offset;
if (length(Vrange)+length(Vrange_cathodal))<samplesperscan
    Vrange=[Vrange 1.3+offset];
    Vrange_cathodal=(anodal_lim:-stepsizeCV:(cathodal_lim+vlim*2/samplesperscan))+offset;
end
Voxid=find(Vrange<=Vox+0.01); Voxid=Voxid(end);
Voxid=length(Vrange_cathodal)+length(Vrange)-Voxid;

parameters.Voxid=Voxid;
parameters.Vrange=Vrange;
parameters.Vrange_cathodal=Vrange_cathodal;

parameters.LFPterm='.ncs';
parameters.includeMovementThres=1;          %include this threshold to remove singals
parameters.samplesperscan=samplesperscan;
parameters.rateFSCV=sampling_freq;


%Get plotparam (e.g. from getplotsettings) and merge with parameters above
    %default settings
   
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
    plotParam.filtlfp=[13 28];        %default BETA?????
    plotParam.envwin=round(rateEphys*.5/mean(plotParam.filtlfp));       %envelop over 1/2 cycle of targetd osc 
    plotParam.filtlfp2=[10 17];         %freq band of lfp filter 2 low beta
    plotParam.envwin2=round(rateEphys*.5/mean(plotParam.filtlfp2));       %envelop over 1/2 cycle of targetd osc 
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

plotParam.colormap=parula;
%plotparam.markercolors=[1 1 1; 0 1 1; 1 1 1; 1 0 1; 1 1 0]; %reward color yel
plotParam.markercolors=[1 1 1; 0 1 1; 1 1 1; 1 0 1; 1 0 0];
plotParam.cminshow=0;
%plotParam.alnevt=alnevt;


%plotParam.ratelfp=rateEphys;   %Moved to param.ratelfp
%variables that change upon calling
%plotParam.cscNames=ncschannels;      %csc names (labels from map)   %NOT NEEDED in param.cscNames now;
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
plotParam.physid=find(~ismember(1:length(ncschannels),plotParam.lfpid));
end

parameters.plotParam=plotParam;
code=reshape(eventcodes,2,length(eventcodes)/2);
map=reshape(csc_map,2,length(csc_map)/2);
parameters.eventcodes=code';
parameters.cscmap=map';
end