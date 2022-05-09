%re-convert LFP trial-split data as generated by getFSCV.m
%spectral sub and downsample to smaller sized files
function reconvertfscv(pathdir,chs,varargin)
%4/2019 add nofscv flag for no spec interpolation if no fscv recorded
%run chronicXXchconfig first to load ncschannels for that session
%get csc #'s & their categories from channel labels provided for given session
%ncschannels={'p2-p3','p2','p3','p3-p5','p5','pl1-pl2','pl1','pl2','cl1','cl1-cl4','cl4',...
% %   'cl5','cl4-cl5','s6','s5','s6-s5','s4','s5-s4','s3','s2','s3-s2',...
%    's1','s2-s1','eyed','eyex','lickx','pulse'};          %chronic 38
%1/10/19
%add option
%only reconvert select nlx chs using 'nlxchsel', no da
                %already previously reconverted everything lese
                %already ran syncsigs again (by default merges everything..
                %so have to run all chs in this step)
ncschannels=chs;
d=dir(pathdir);
cd(pathdir);
filenames={d.name};
targfiles=strfind(filenames,'fscv_multi');
processfiles=find(~cellfun(@isempty,targfiles));
%contains function does not work for 2013 for cells 
files=filenames(processfiles);
fscvchnames={};
patra_map_bipolar;            %get default csc_map
cleom=0;
argnum=1;
movecsc=0;
cflag='';
splitfiles=0;
mergeflag=0;
getda=1;
getlfps=1;
nospecfun=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'nofscv'
            %no spectralsubtraction flag
            nospecfun=1;
        case 'map'
            %user provided csc map
            argnum=argnum+1;
            csc_map=varargin{argnum};
        case 'nlxchsel'
            %user provided select nlx chs in ncschannels to reconvert
            %don't convert da in this case
            getda=0;
        case 'fscvchnames'
            %user provided select nlx chs to reconvert
            argnum=argnum+1;
            fscvchnames=varargin{argnum};
       case 'split'
            %split files for separate chs
            splitfiles=1;
        case 'merge'
            %merge files in directory to single file for each ch
            mergeflag=1;
        case 'cleomoves'
            %use scripts to remove artifacts due to movement as done with
            %Ophelia in past 
            cleom=1;
        case 'movecsc'
            %channel for movement detection
            argnum=argnum+1;
            movecsc=varargin{argnum};
        case 'cleolick37'
            %cleo lick channel on 37 (single ch)
            cflag='cleolick37';
    end
    argnum=argnum+1;
end
cscexists=1;
[lfpchs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map,cflag);

%lfpchs=[3, 4, 5, 8,9, 11, 14, 15, 27:32];    %chronic38
%bipolarchs={[3, 4], [4, 5],[8,9],[11, 14], [14, 15], [27, 28], [28, 29], [30, 31], [31 32]};
%otherchs=[33:35,39:41,42];   %eye, lick, pulse
%sumchs=[39:41];            %sum up signals from lick accel x,y,z data
selectCSC=[lfpchs otherchs];
filtPhys=[0 100];        %filter physiological signals in this range
cscNames={};
dsN=30;         %factor for downsampling

%get directory for csc data & make new directory for processed data
%files = uigetdir2(pwd,'Select multifscv trial-split recording files'); 
%sep2=findstr(filesep,files{1});
%PathName=files{1}(1:sep2(end));
%pathSave=[pathdir(1:end-1) '_pro\'];
pathSave=[pathdir(1:end-1) '_pro' filesep];
%pathDown=[PathName(1:end-1) '_downsampled\'];
if ~isdir(pathSave) 
    status = mkdir(pathSave)
  %  status2 = mkdir(pathDown)
end

samplingFreq=30e3;
load(files{1});
nlxIDs=[];
if ~exist('tsLFP')
    cscexists=0;        %if no tsLFP variable in file, means no csc signal
else
    samplingFreq=getSampleRate(tsLFP);
	smoothWindow=round(samplingFreq*0.0384);       %.0384 s = .5 cycle beta 
    digids=regexp(nlxFileNames,'\d');       %assume csc chs are same for all files
    for ii=1:length(nlxFileNames)
        ab=nlxFileNames{ii};
        nlxIDs(ii)=str2num(ab(digids{ii}));
    end
end
alreadydsts=0;
%%
for ifile=1:size(files,2)
    disp(['processing file # ' ...
        num2str(ifile) '/' num2str(size(files,2))]);
    load(files{ifile});
    filename=files{ifile};
if cscexists
    alreadydsts=0;
    ts=[];
    samples={};
    samples_dsonly={};      %unprocessed, just downsamples only like eye
    renlxFileNames={};
    sumchsflag=0;        %already summed sum ch
    bipolarFlag=0;      %already processed bipolar pair flag
    saveID=1;           %save channel id (ie. added chs will be saved because of created bipolar chs)
    bchs=[bipolarchs{:}];       %reset bipolar each file, since will remove from this list as processed
    bchsfirstids=bchs(1:2:end);     %just first ones, ie for p1-p2, get p1
    refdata=[];
    %cycle all targeted channels
    %{
    glitchdata=[];      %samples for use as glitch/movement detector
    if movecsc~=0
        mid=find(nlxIDs==movecsc);
        glitchtemp=spectralSub(samplesLFP{mid},ceil(samplingFreq));
        mdifsamp = abs(diff(glitchtemp));   %ie differential of data to look for change in slope
        thresh=mean(mdifsamp)*5;     %arbitrarily set threshold for given signal based on mean differential
        mdata=deglitchDiff(glitchtemp, 30,'thres',thresh);
        mdataf=filterLFP(mdata,ceil(samplingFreq),[60 1000]);
        mdatafs=mdataf.^2;
        glitchdata=smoothwin(mdatafs,ceil(.1*samplingFreq));      %100 ms
    end
    %}
    for ich=1:size(selectCSC,2)
        %get chid of stored nlx samples as matched to targeted ich
        chid=find(nlxIDs==selectCSC(ich));
        %get new ch name if stored in different order
        renlxFileNames{saveID}=nlxFileNames{chid};
        if ~ismember(selectCSC(ich),otherchs)
            %is LFP signal, process
            %need to spectral sub/deglitch, etc.
            samplestemp=samplesLFP{chid};
            if ~nospecfun
                %da recorded with lfp, spectral interpolation needed
                samplestemp=spectralSub(samplesLFP{chid},ceil(samplingFreq)); %spectral subtraction
            end
            %samplesdg=deglitchDiff(samplestemp, 30);       %deglitch with
            %30 pt padding, line interpolation
            
            %fillgap function in deglitchandsmooth too time consuming
            %revert
           % samplesdg=deglitchandsmooth(samplestemp, 20);       %deglitch with 20 pt padding & smooth/not just line
              samplesdg=deglitchDiff(samplestemp, 30);       %deglitch with 20 pt padding & smooth/not just line

            %downsample processed signal
            subsamples = downsampleSingleBatch(samplesdg,dsN);
            if alreadydsts==0
                %downsample timestamps once per file
                ts = reshape(tsLFP(1:dsN:end), 1, []);
                alreadydsts=1;
            end
            %figure check
            %{
            fig2=figure;
            set(fig2,'color',[1 1 1], 'units','pixels','position',[100 100 800 400]);
            hold on;
            plot(tsLFP-tsLFP(1),samplesLFP{chid})
            xlabel('time (s)')
            ylabel('amplitude (V)')
            title('LFP (time-domain)');
            ylim([-2e-4 2e-4]);
            xlim([35 36]);
           % plot(tsLFP-tsLFP(1),samplestemp);
            plot(ts-tsLFP(1),subsamples);
            samplesfilt=filterLFP(subsamples,samplingFreq/dsN,[13 36]);
            %if want envelop for this channel, envelope
            samplesfiltp=samplesfilt.^2;   %get power V^2
            samplesfiltps=smoothwin(samplesfiltp,round(samplingFreq/dsN*0.0384));   %smoothing
            yyaxis right;
            plot(ts-tsLFP(1),samplesfiltp,'linewidth',1,'linestyle','--','color',[.6 .6 .6]);
                   plot(ts-tsLFP(1),samplesfiltps,'linewidth',1,'color',[.4 .4 .4]);
             %}
            %store processed/downsampled LFP
            samples{saveID}=subsamples;
            refdata=subsamples;     %reference below in bipolar artifact remoal
            %store unprocessed/downsampled LFP
        elseif ismember(selectCSC(ich),sumchs) && sumchsflag==0
            %if summing channel for lick xyz
            %take absolute mag of each x,y,z signal component & sum
            chidsum=[];
            subsamples=abs(samplesLFP{chid});
            if length(sumchs)>1
                %patra 3 channels for lick xyz
                for ichsum=1:length(sumchs)-1
                    %summate with next 2 consecutive channels y/z lick
                    %store just downsampled signals for channels y/z
                    chidsum=find(nlxIDs==selectCSC(ich+ichsum));
                    subsamples=subsamples+abs(samplesLFP{chidsum});
                    samples_dsonly{saveID+ichsum}=downsampleSingleBatch(samplesLFP{chidsum},dsN);
                    %save raw data for other y/z channels 
                    %skip repeating this in next cycles because of sumchsflag
                    samples{saveID+ichsum}=samples_dsonly{saveID+ichsum};
                end
            else
                %single lick channel (strain sensor cleo)
                %keep subsamples
            end
            %process summed samples
            subsamples=filterLFP(subsamples,samplingFreq,filtPhys);
            %save summed & processed samples
            samples{saveID}=downsampleSingleBatch(subsamples,dsN);
            sumchsflag=1;        %turn on flag so do not process other sum channels
                        if alreadydsts==0
                %downsample timestamps once per file
                ts = reshape(tsLFP(1:dsN:end), 1, []);
                alreadydsts=1;
            end
        elseif ismember(selectCSC(ich),otherchs) && ~ismember(selectCSC(ich),sumchs)
            %just downsample if not LFP signal (ie. pulse/eye diam)
            samples{saveID}=downsampleSingleBatch(samplesLFP{chid},dsN);
                        if alreadydsts==0
                %downsample timestamps once per file
                ts = reshape(tsLFP(1:dsN:end), 1, []);
                alreadydsts=1;
            end
        end
        %if bipolar selected channel, subtract channels & process LFP
        %csc naming convention here is to use 0 between ch #'s
        %eg 2-3 is csc203.ncs
        if ismember(selectCSC(ich),bchs)
         %   bid=find(selectCSC(ich)==bchs);
           % chid=find(nlxIDs==bchs(bid));
            %07/04/2018 just find first ids for each csc selected
            bid=find(selectCSC(ich)==bchsfirstids);
                                    numb=length(bid);

            %assume could be more than one bipolar sub for this csc
            countb=0;
            bid=find(selectCSC(ich)==bchs); %ids from master list
            bid=bid(find(mod(bid,2)==1));   %only odd #'s (ie first ids);
            idstodelete=[];
            while countb<numb
                saveID=saveID+1; %increment idx to save created samples
                countb=countb+1;
                %chid=find(nlxIDs==bchsfirstids(bid(countb)));                
                chid=find(nlxIDs==bchs(bid(countb)));              
                chid2=find(nlxIDs==bchs(bid(countb)+1));         
                %chid=find(nlxIDs==bchsfirstids(bid));
                %assume next bipolar listed ch is the one to subtract, skip next cycle
                %chid2=find(nlxIDs==bchs(bid+1));        
                %create bipolar save name
                renlxFileNames{saveID}=['csc' num2str(bchs(bid(countb))) '0' ...
                     num2str(bchs(bid(countb)+1)) '.ncs'];
                samplesbipolar=samplesLFP{chid}-samplesLFP{chid2};
                samplestemp=samplesbipolar;
                if ~nospecfun
                    %spec sub for fscv coupled lfp data
                %same processing for LFP
                if cleom==0
                samplestemp=spectralSub(samplesbipolar,samplingFreq); %spectral subtraction
                else
                    samplestemp=spectralSub(samplesbipolar,samplingFreq,'cleo');
                end
                end
               % samplesdg=deglitchandsmooth(samplestemp, 30);       %deglitch with 30 pt padding
                badids=[];
                badidsh=[];
               if cleom==1
                    difsamp = abs(diff(samplestemp));   %ie differential of data to look for change in slope
                    thresh=mean(difsamp)*4;     %arbitrarily set threshold for given signal based on mean differential
                    samplesdg=deglitchDiff(samplestemp, 30,'thres',thresh);       %deglitch with 20 pt padding & smooth/not just line
                   % samplesdg=deglitchDiff(samplesdg, 100,'thres',thresh*2); 
                   % glitchdataspt=filterLFP(samplesdg,ceil(samplingFreq),[0 5]);
                    %glitchdatasp=smoothwin(glitchdataspt.^2,ceil(.1*samplingFreq));      %100 ms
                    
                    glitchdataspt2=filterLFP(samplesdg,ceil(samplingFreq),[60 10000]);
                    glitchdatasp2=smoothwin(glitchdataspt2.^2,ceil(.1*samplingFreq));
                    %threshml=mean(glitchdatasp)*4;   
                    threshmh=mean(glitchdatasp2)*4;                    
                    %badids=find(glitchdatasp>threshml);
                    badidsh=find(glitchdatasp2>threshmh);
                   % badids=sort([badids; badidsh]);
                    samplesdg=filterLFP(samplesdg,ceil(samplingFreq),[2 10000]);        %remove low freq fluctuations , change low cut to 2 Hz from 5 Hz, 04/19/2019
                 %   deglitchm=deglitchIds(samplesdg,badids,ceil(.01*samplingFreq));
                else
                    samplesdg=deglitchDiff(samplestemp, 30);       %deglitch with 20 pt padding & smooth/not just line
                end
                
                subsamples = downsampleSingleBatch(samplesdg,dsN);
                
                if cleom==1
                    badisd=badidsh./dsN;
                    badidsd=round(badisd);
                    satudata=smoothwin(refdata'.^2,ceil(.1*samplingFreq/dsN));
                    thres=2.5e-7;           %for fscv coupled data, remove too high dc fluctuations
                    if nospecfun
                        thres=2.5e-6;
                    end
                    satuids=find(satudata>thres); %too high signal original
                    if size(badidsd,1)<size(badidsd,2)
                        badidsd=badidsd';
                    end
                    badids=sort([badidsd; satuids]);
                    badids=unique(badids);
                    deglitchm=deglitchIds(subsamples,badids,ceil(.03*samplingFreq/dsN));
                    subsamples=deglitchm;
                end
                if alreadydsts==0
                    %downsample timestamps once per file
                    ts = reshape(tsLFP(1:dsN:end), 1, []);
                    alreadydsts=1;
                end              
                samples{saveID}=subsamples;
                idstodelete=[idstodelete bid(countb) bid(countb)+1];
            end            
            bchs(idstodelete)=[];
           % bchs=bchs(numb*2+1:end);       %take out processed pairs from list, since processed
        end
        
        saveID=saveID+1;

    end
    
    %save processed data in new paths with original variables
    samplesLFP=samples;
    tsLFP=ts;           %store signals in original variable names for comptability with other funcs.
    nlxFileNames=renlxFileNames;
    if splitfiles==0
        %default save merged file
        save([pathSave, filename],'fscv','samplesLFP','tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','fscvchnames')
        disp(['saved to: ' pathSave filename]);
    end
    if splitfiles==1
        %save nlx signals as separate files and fscv as separate file
        %save in nlx format compatible with dg / lfp lib
        dg_Nlx2Mat_Timestamps=tsLFP.*1e6;     %convert to Microseconds
        dg_Nlx2Mat_Samples={};
        dg_Nlx2Mat_SamplesUnits='V';
        filenum=filename(regexp(filename,'\d'));
        for nlxch=1:size(samplesLFP,2)
            %save each nlx channel as separate file
            dg_Nlx2Mat_Samples=samplesLFP{nlxch};   %already in volts
            sepidx=strfind(renlxFileNames{nlxch},'.ncs');
            cscname=renlxFileNames{nlxch}(1:sepidx-1);
            savecsc=[cscname '_' filenum];
            save([pathSave, savecsc],'dg_Nlx2Mat_Samples','dg_Nlx2Mat_SamplesUnits','dg_Nlx2Mat_Timestamps','-v7.3');
        end
        %save fscv file without nlx signals, but with events
        if getda
            if ~nospecfun
                %da & lfp signals
                save([pathSave, filename],'fscv','tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','fscvchnames')
            else
                %only lfp signals
                 save([pathSave, filename],'tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','fscvchnames')
            end
        end
        prompt = ['saved as, split ' num2str(size(samplesLFP,2)) ' chs: ' pathSave filename];
    end       
    
    else
        %no lfp signals, just switching folders....
      if splitfiles==1        
        filenum=filename(regexp(filename,'\d'));        
        %save fscv file without nlx signals, but with events
        if getda
        save([pathSave, filename],'fscv','NlxEventTTL','NlxEventTS','fscvchnames')
        end
        prompt = ['saved as, ' pathSave filename];
    end  
        
        
end


end


end
