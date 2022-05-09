%automatically load saved trial multifscv files 
%create ipca data automatically and find correlated cv's within time window
%around target at 30s? usually
%run R command to find corralated cv's for each file, assume highest CV
%04/08/2018 for patra behavior trials
clear all;
selch=1;
global plotParam processed parameters
selch=[ 1 2 3 4];            %fscv channels to plot
settingsfile=[];                    %load settings file (plotParam), if empty default
filtLFP=[17 34];  
csc_map={};
%[parameters,csc_map,event_codes]=getparams('patra');
[parameters,csc_map,event_codes]=getparams('patrabipolar','default');

tParam.timeWin=[20 40];       %window to search for cv's if empty, will search all file
tParam.samplingRate=10;        %10 samples per sec
t_start=0;  t_end=60;           %in seconds
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

Vrange=parameters.Vrange;
Vrange_cathodal=parameters.Vrange_cathodal;

parameters.NCSchannels=[];           %load all csc data
parameters.badChannels=[];

files = uigetdir2(pwd,'Select files'); nFiles=size(files,2);
sep=findstr(filesep,files{1});
PathName=files{1}(1:sep(end)-1);
FileName=files{1}(sep(end)+1:end);
%load single file to get default parameters
PathName=[PathName '\'];
[Iread,LFPread,samplesNCS]=loadAll(PathName, FileName,parameters,selch);
loadedCSCchs=[];
for ii=1:length(LFPread.LFPchNames)
    getnum=strfind(LFPread.LFPchNames{ii},'.ncs');
    loadedCSCchs(ii)=str2num(LFPread.LFPchNames{ii}(4:getnum-1));
end
PathNameNew={};
for ich=1:length(selch)
PathNameNew{ich}=[PathName 'analyzed\ch' num2str(selch(ich)) '\' ];
if ~isdir(PathNameNew{ich})
    status = mkdir(PathNameNew{ich});
end
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

rateLFP=LFPread.LFPsamplingfreq;
tsLFP=LFPread.LFPts;
%nlx.filtlfp=[16 34];        %standard freq band
nlx.filtlick=[0 10];        %freq bands of lick filter

%nlx.filtlfp=[22 34];            %04/16/2018 analyze high freq beta only ch3
nlx.filtlfp=[17 34];            %chronic 58
nlx.filtlfp2=[10 16];            %chronic 58

status=0;
%param.NCSchannels=[];
for idFile=1:nFiles
    FileName=files{idFile}(sep(end)+1:end);
    [processed.Iread,LFPread,samplesNCS]=loadAll(PathName, FileName,parameters,selch);
    samplesNCS(removeIDs,:)=[];     %remove unused CSC channels, defined above
    processed.LFPread=LFPread;
    ncsread=processed.LFPread;
    %remove unused CSC channel names, defined above
    ncsread.LFPchNames(removeIDs)=[];

    %get duration info on fscv data loaded
    sizeData=size(Iread(1).data); 
    samplesperscan=sizeData(1);
    lengthData=sizeData(2);
    t_end=lengthData;
    nlx_events=[];
    events=[];
    Idata.anodal=[];
    Idata.cathodal=[];
    %process nlx data
    %get nlx events into fscv domain
    if size(Iread(1).events,1)>0
        [nlx_events,events]=readEvents(Iread(1).events,event_ref_code);
    end
    filtLFP=[];
    nlx.resampled=[];
    nlx.resampled2=[];
    %process CSC signals with input filter
    for ich=1:length(nlx.selectCSC)
        %lfp signal --> filter, square, envelope
        if ismember(ich,nlx.lfpid)
            %filter at "beta-band" as defined in nlx.filtlfp
            nlx.resampled(ich,:)=filterLFP(samplesNCS(ich,:),rateLFP,nlx.filtlfp);
            %square & envelope signal
            nlx.resampled(ich,:)=nlx.resampled(ich,:).^2;   %get power V^2
            winlength=round(rateLFP*.5/mean(nlx.filtlfp));
            nlx.resampled(ich,:)=smoothwin(nlx.resampled(ich,:),winlength);   %smoothing
            if isfield(nlx,'filtlfp2')
                %if another filter band specified 05/03/2018
                %filter at "low beta-band" as defined in nlx.filtlfp2
                nlx.resampled2(ich,:)=filterLFP(samplesNCS(ich,:),rateLFP,nlx.filtlfp2);
                %square & envelope signal
                nlx.resampled2(ich,:)=nlx.resampled2(ich,:).^2;   %get power V^2
                winlength=round(rateLFP*.5/mean(nlx.filtlfp2));
                nlx.resampled2(ich,:)=smoothwin(nlx.resampled2(ich,:),winlength);   %smoothing 
            end
        %lick signal --> filter 
        elseif ismember(ich,nlx.lickid)
            nlx.resampled(ich,:)=filterLFP(samplesNCS(ich,:),rateLFP,nlx.filtlick);
        %eye/pulse signal --> dont do anything already filt in reconvertFSCV.m            
        else
            nlx.resampled(ich,:)=samplesNCS(ich,:);
        end
    end
    
    %for each channel of data (in selch) get detected da 
    samples={};         %temporary storage of samples detected
    for ich=1:length(selch)
        chnum=selch(ich);
        processed.Iread(chnum).data = filterDisplay(samplesperscan,[],processed.Iread(chnum).data);
        %determine bg subtraction point
        [aa, closeevent]=min(abs(events-mean(tParam.timeWin).*tParam.samplingRate));
        idxbgref=events(closeevent)-offsetbg;     %1 s before first appropriate event (display fix)
        [processed.Isub(chnum).data,processed.BG(chnum).data,BGrefmatrix]=...
            refreshPlots(processed.Iread(chnum).data,tParam.samplingRate,idxbgref,nlx.BGavg,t_end); 
        [processed.Isub(chnum).rawdata,processed.BG(chnum).rawdata,bb]=...
            refreshPlots(processed.Iread(chnum).rawdata,tParam.samplingRate,idxbgref,nlx.BGavg,t_end); 
        %compute pca based on fixed BG sub point
        Ipcr = getpct(processed, parameters,chnum);
            %compute pca components of data
        %get behavioral parameters based on task events & actions
        info{chnum}.behav=calcBehav(processed,tParam.timeWin,event_codes);    
        %find ids where HF glitches, skip these data points
        glitchids=findglitches(processed.Iread(chnum).rawdata,parameters.glitchThres);
        %get pct dopamine signal to remove M/Q thres & glitches
        Ipcr = getpctdirect(processed.Isub(chnum).data,...
            'ats',parameters.CV,'cts',parameters.Cts,...
            'qthres',parameters.QaDAPH,'rthresm',parameters.RThresOut,...
            'imaxbg',max(abs(processed.Iread(chnum).data(:,optimalidx))),...
            'glitchids',glitchids);
        %find local "transient' increases in DA & their parameters
        detectedda=detectdatransients(processed.Iread(chnum).data,parameters,glitchids,tParam);       %input dataFSCV, output samplesTemp
        
        
        if ~isempty(detectedda.cv)
            rr=abs(corr(detectedda.cv,parameters.K(:,1)));         %correlate retrieved CV's to da template & keep best
            maxCV=detectedda.cv(:,find(rr==max(rr)));  
            Idata.anodal=maxCV(1:round(samplesperscan/2));
            Idata.cathodal=maxCV(round(samplesperscan/2+1):round(samplesperscan));          
        end
        Iread=processed.Iread(chnum);
        Isub=processed.Isub(chnum);
        info=info{chnum};
        saveName=['CVITdata_ch' num2str(chnum) '_' FileName];
        BG=processed.BG(chnum).data;
        save([PathNameNew{ich} saveName],'Ipcr','nlx','ncsread','parameters','Idata','BG',...
            'Vrange','Vrange_cathodal','Isub','info','events','chnum');
    end
    %display status
    percentcomplete=round(idFile/nFiles*1000)/10;
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
