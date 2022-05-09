function autotrialdirlfponly(pathdir,varargin)
global plotParam processed parameters
%need patra_map_xx config file in config folder
%path dir should be pointed to dir of bigreward_pro or other dir with fscv
%split by trials aligned to mid of file (30s default)
%automatically load saved trial multifscv files 
%create ipca data automatically and find correlated cv's within time window

argnum=1;
selch=[ 1 2 3 4];            %default fscv channels to plot
configmode='patrabipolar';     %for pc, default patra
%configmode='patrabichunky';     %for chunky, needs to specify dir for pcr temp

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'fscvchs'
            %user provided fscv selected channels
            argnum=argnum+1;
            selch=varargin{argnum};
        case 'config'
            %user provided config file
            argnum=argnum+1;
            configmode=varargin{argnum};
    end
    argnum=argnum+1;
end
d=dir(pathdir);
cd(pathdir);
filenames={d.name};
targfiles=strfind(filenames,'fscv_multi');
processfiles=find(~cellfun(@isempty,targfiles));
%contains function does not work for 2013 for cells 
files=filenames(processfiles);
nfiles=length(files);

settingsfile=[];                    %load settings file (plotParam), if empty default
filtLFP=[17 34];  
csc_map={};

tParam.timeWin=[];
tParam.timeWin=[20 40];       %window to search for cv's if empty, will search all file

tParam.samplingRate=10;        %10 samples per sec
t_start=0;  t_end=60;           %in seconds
BGavg=3;                        %average BG +/- samples
event_ref_code=4;               %fix on appear cue (4) for bg averaging.
offsetbg=10;                    %use 5 s offset before event ref for bg averaging
%initialize variables that will be saved
BG=[];
Idata={};
Ipcr={};
Isub={};
Vrange=[];
Vrange_cathodal=[];
chnum=selch;
events=[];
info={};
ncsread={};
parameters={};
plotParam={};


%[parameters,csc_map,event_codes]=getparams('patra');
%[parameters,csc_map,event_codes]=getparams('patrabipolar');
[parameters,csc_map,event_codes]=getparams(configmode);

Vrange=parameters.Vrange;
Vrange_cathodal=parameters.Vrange_cathodal;

parameters.NCSchannels=[];           %load all csc data
parameters.badChannels=[];

sep=findstr(filesep,files{1});
filename=files{1};
%load single file to get default parameters
[Iread,LFPread,samplesNCS]=loadall(pathdir, filename,parameters,selch);
loadedCSCchs=[];
for ii=1:length(LFPread.LFPchNames)
    getnum=strfind(LFPread.LFPchNames{ii},'.ncs');
    loadedCSCchs(ii)=str2num(LFPread.LFPchNames{ii}(4:getnum-1));
end

pathdirNew=[pathdir 'analyzed_lfp' filesep ];
if ~isdir(pathdirNew)
    status = mkdir(pathdirNew);
end

sampleratencs=LFPread.LFPsamplingfreq;
nlx.cscNames=getcscnames(loadedCSCchs,csc_map);     
%remove licky & lickz, since lickx has magnitude (generated in
%reconvertFSCV), and remove eyeY (do not use)
removeIDs=find(ismember(nlx.cscNames,{'licky','lickz','eyey'})==1);
nlx.cscNames(removeIDs)=[];
%Need to remove above everytime file loaded in directory

%load ch names from patra_map based on loaded csc chs
%populate plotParam
getplotsettings(filtLFP,nlx.cscNames,event_codes,sampleratencs,settingsfile);       
nlx=plotParam;      %nlx.lfpid nlx.eyeid, nlx.. have the id's of categorized chs
nlx.selectCSC=[nlx.lfpid nlx.eyeid nlx.lickid nlx.pulseid];
parameters.BGavg=nlx.BGavg;

parameters.sampleratencs=sampleratencs;

rateLFP=LFPread.LFPsamplingfreq;
tsLFP=LFPread.LFPts;
nlx.filtlick=[0 10];        %freq bands of lick filter
    nlx.filtbetal=[13 18];        %low beta
    nlx.filtbetah=[18 33];        %high beta
    nlx.filtgammal=[40 58];           %low gamma
status=0;
%param.NCSchannels=[];
for ifile=1:nfiles
    filename=files{ifile};
    [processed.Iread,LFPread,samplesNCS]=loadall(pathdir, filename,parameters,selch);
    samplesNCS(removeIDs,:)=[];     %remove unused CSC channels, defined above
    processed.LFPread=LFPread;
    ncsread=processed.LFPread;
        ncsread.LFPchNames(removeIDs)=[];

    %get duration info on fscv data loaded
    sizeData=size(Iread(1).data); 
    samplesperscan=sizeData(1);
    lengthData=sizeData(2);
    t_end=lengthData;
    nlx_events=[];
    events=[];

    %for each channel of data (in selch) get detected da 
    %process nlx data
    %get nlx events into fscv domain
    if size(Iread(1).events,1)>0
        [nlx_events,events]=readEvents(Iread(1).events,event_ref_code);
    end
    %get CSC signals 
    for ich=1:length(nlx.selectCSC)
        %
        nlx.resampled(ich,:)=samplesNCS(ich,:);
        %lick signal --> filter 
        if ismember(ich,nlx.lickid)
            nlx.resampled(ich,:)=filterLFP(samplesNCS(ich,:),rateLFP,nlx.filtlick);
        %eye/pulse signal --> dont do anything already filt in reconvertFSCV.m            
        
        end
    end
    info.behav=calcBehav(processed,tParam.timeWin,event_codes);    

    saveName=['lfpdata_' filename];
    save([pathdirNew saveName],'nlx','ncsread','parameters','info','events');
    %display status
    percentcomplete=round(ifile/nfiles*1000)/10;
    if percentcomplete>=25 && status==0
        disp([num2str(percentcomplete) '%']);
        status=1;
    elseif percentcomplete>=50 && status==1
        disp([num2str(percentcomplete) '%']);
        status=2;
    elseif percentcomplete>=75 && status==2
        disp([num2str(percentcomplete) '%']);
        status=3;
    end

end

end
