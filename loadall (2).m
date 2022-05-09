function [Iread,LFPread,samplesLFP]=loadall(PathName,FileName,param,selch,varargin)
%load all ephys and fscv data
%selch= array of channels selected for multichannel recording system to
%limit to 4
Iread=[];
LFPread=[];
samplesLFP=[];
numch=4;
isTXT=0;
argnum=1;
D={};
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dir'
            argnum=argnum+1;
            %provide directory instead of loading each time
            D=varargin{argnum};
    end
    argnum=argnum+1;
end
if isfield(param,'badChannels')
for ii=1:numch
        if ~any(ii==param.badChannels)
        Iread(ii).LPFdata=nan(1000,3000);
        Iread(ii).data=Iread(ii).LPFdata;
        Iread(ii).events=0;
        end
end
end
actualNumCh=1;
if isempty(D)
D=dir(PathName);
end
FileNames{1}=FileName;       %Store original filename, open mat file (fscvData)
PathNameorig=PathName;
%determine if original tarheel file format selected
%if there is a .txt file in same folder with same name
istar=find(ismember({D.name},[FileName '.txt'])==1);
if ~isempty(istar)
    %go to \cvtotxt folder
    PathName=[PathName 'cvtotxt\'];
    D=dir(PathName);
    dirid=find(ismember({D.name},['0_' FileName '_cvtotxt'])==1);
    if isempty(dirid)
        %try with caps
        dirid=find(ismember({D.name},['0_' FileName '_CVtoTXT'])==1);
        FileName=['0_' FileName '_CVtoTXT'];
            %PathName=[PathNameorig 'CVtoTXT\'];
    else
        FileName=['0_' FileName '_cvtotxt'];
        PathName=[PathName 'cvtotxt\'];
    end
end
    
fscvdata=load([PathName,FileName]);
FileNum=strsplit(FileName,'_');
FileNum=FileNum(2); FileNum=char(FileNum);
%contains does not work in 2013 matlab in chunky, switched to cellfun
%isTXT=contains(FileName,'TXT','IgnoreCase',true) ;  
FileNamel=lower(FileName);  %make lowercase
isTXTa=strfind(FileNamel,'txt');
if ~isempty(isTXTa)
    isTXT=1;
end

%converted exported data usually starting with CVITdata...
isColor=strfind(FileName,'Colordata');
if isempty(isColor)
    isColor=strfind(FileName,'CVITdata');
end
if ~isempty(isColor)
        load([PathName,FileName]);
    for ii=1:numch
         if exist('Isubdata')
           Iread(ii).LPFdata=Isubdata;
           Iread(ii).data=Isubdata;
           Iread(ii).events=0;
         end
         if exist('Isubavg')
          Iread(ii).LPFdata=Isubavg;
           Iread(ii).data=Isubavg;
           Iread(ii).events=0;
         end
         if exist('Isub')       %%check here 06/20 added
          Iread(ii).LPFdata=Isub.avgIsub;
           Iread(ii).data=Isub.avgIsub;
           Iread(ii).events=0;
         end
    end
end

%Synchronized converted fscv/behavior data (from getFSCV)
%if isempty(isTXT) &&  isempty(isColor)
splitfiles=0;
if isfield(fscvdata,'fscv') 
    %check if csc files in directory, need to open these separately later..
    %D=dir(PathName);
    if ~isempty(param.NCSchannels)
    matchcsc=regexpi({D.name},['csc' num2str(param.NCSchannels(1))]);
    idscsc=~cellfun('isempty',matchcsc);
    if any(idscsc)==1
        splitfiles=1;       %csc trial split files exist separately        
    end
    end
   % load([PathName,FileName]);
   if length(selch)==1
        if intersect(selch,1:4)
            switch selch
                case 1
                    Iread(1).data=flipud((fscvdata.fscv.ch0)');
                case 2
                    Iread(2).data=flipud((fscvdata.fscv.ch1)');
                case 3
                    Iread(3).data=flipud((fscvdata.fscv.ch2)');
                case 4
                    Iread(4).data=flipud((fscvdata.fscv.ch3)');  
             end
            Iread(selch).LPFdata=Iread(selch).data;
            Iread(selch).rawdata=Iread(selch).data;
            Iread(selch).events=(fscvdata.fscv.events)';
                        Iread(1).LPFdata=Iread(1).data;
            Iread(1).rawdata=Iread(1).data;
            Iread(1).events=(fscvdata.fscv.events)';
        end
   else
        for ii=1:length(selch)
            switch selch(ii)
                case 1
                    Iread(1).data=flipud((fscvdata.fscv.ch0)');
                    Iread(1).LPFdata=Iread(1).data;
                    Iread(1).rawdata=Iread(1).data;
                case 2
                    Iread(2).data=flipud((fscvdata.fscv.ch1)');
                    Iread(2).LPFdata=Iread(2).data;
                    Iread(2).rawdata=Iread(2).data;

                case 3
                    Iread(3).data=flipud((fscvdata.fscv.ch2)');
                    Iread(3).LPFdata=Iread(3).data;
                    Iread(3).rawdata=Iread(3).data;

                case 4
                    Iread(4).data=flipud((fscvdata.fscv.ch3)');
                    Iread(4).LPFdata=Iread(4).data;
                    Iread(4).rawdata=Iread(4).data;

                    
            end
            Iread(selch(ii)).events=(fscvdata.fscv.events)';
        end
   end
end

%Synchronized ephys data with converted fscv/behavior data
if isfield(fscvdata,'samplesLFP') 
    count=1;
    samplesLoaded=[];
    if ~isempty(param.NCSchannels)
        for idx=1:length(param.NCSchannels)
            %load in order of appearing in NCSchannels list
            matchID=regexpi(fscvdata.nlxFileNames,['csc' num2str(param.NCSchannels(idx)) param.LFPterm]);
            empties = cellfun('isempty',matchID);   %get those not matching in logical array           
            targetID=find(empties~=1);      %get index for targeted channel
            samplesLoaded(:,idx)=fscvdata.samplesLFP{targetID};  %store sample for a selected index
           % LFPread(idx).samplesLFP=samplesLFP{targetID};
           LFPread.LFPchNames{idx}=fscvdata.nlxFileNames{targetID};
        end
    else
       % for idx=1:size(samplesLFP,2)
           % LFPread(idx).samplesLFP=samplesLFP{idx};
         %   count=count+1;
         LFPread.LFPchNames=fscvdata.nlxFileNames;
            samplesLoaded=vertcat(fscvdata.samplesLFP{:});   %store all samples
        %end
    end
       % a={LFPread(:).samplesLFP};
      %  a={LFPread(:).NlxEventTTL};    b=vertcat(a{:});    eventTTL=b(1,:);
        LFPread.LFPeventTTL=fscvdata.NlxEventTTL;        %NLX TTL encoded w/ high nlx sampling rate 
       % a={LFPread(:).NlxEventTS};    b=vertcat(a{:});    eventTS=b(1,:);
        LFPread.LFPeventTS=fscvdata.NlxEventTS;
        %a={LFPread(:).tsLFP};    b=vertcat(a{:});    tsLFP=b(1,:);
        LFPread.LFPts=fscvdata.tsLFP;
        rateLFP=getSampleRate(fscvdata.tsLFP);
        LFPread.LFPsamplingfreq=rateLFP;
     %   a={LFPread(:).nlxFileNames};
       % LFPread.LFPchNames=fscvdata.nlxFileNames;
        samplesLFP=samplesLoaded;
else
    %set defaults if no ephys
    LFPread.LFPsamplingfreq=[];
    LFPread.LFPeventTS=[];
    LFPread.LFPeventTTL=[];
    LFPread.LFPchNames=[];
    LFPread.LFPts=0;
end

%if split fscv/ephys data in directory
if isfield(fscvdata,'fscv') && splitfiles==1
    %get csc files separately in direcotry
    samplesLoaded=[];
    namx=strfind(FileName,'.mat');
    trialnum=FileName(namx-3:namx-1);      %get trialnum to find corresponding files    
    %D=dir(PathName);
    %get csc signals defined in NCSchannels
    if ~isempty(param.NCSchannels)
        for idx=1:length(param.NCSchannels)
            %load in order of appearing in NCSchannels list
            matchcsc=regexpi({D.name},['csc' num2str(param.NCSchannels(idx)) '_']);
            idscsc=~cellfun('isempty',matchcsc);       %make array of zeros/1 for matching
            matchtrial=regexpi({D.name},['_' trialnum]);
            idstrial=~cellfun('isempty',matchtrial);  
            fileid=find((and(idscsc,idstrial))==1);     %get file id that matches both target trial/cscnum
            tempload=load([PathName,D(fileid).name]);
            samplesLoaded(:,idx)=tempload.dg_Nlx2Mat_Samples;  %store samples, 
            %timestamps in tempload are in microseconds and should match tsLFP already loaded 
            %get/store csc name
            cscnameb=strfind(D(fileid).name,'_');
            cscname=[D(fileid).name(1:cscnameb-1) param.LFPterm];                
            LFPread.LFPchNames{idx}=cscname;
        end    
    end
        LFPread.LFPeventTTL=fscvdata.NlxEventTTL;        %NLX TTL encoded w/ high nlx sampling rate 
        LFPread.LFPeventTS=fscvdata.NlxEventTS;
        LFPread.LFPts=fscvdata.tsLFP;
        rateLFP=getSampleRate(fscvdata.tsLFP);
        LFPread.LFPsamplingfreq=rateLFP;
        samplesLFP=samplesLoaded;
end


%Tarheel system recording
if isTXT &&  isempty(isColor)
    numch=4;
    FileNum=strsplit(FileName,'_');
    FileNum=FileNum(2); FileNum=char(FileNum);
    dirIndices=findOtherChannelFiles(D,FileName);
    numch=length(dirIndices);
    if numch<length(selch)
        warning('channels selected less than available')
        selch=1:numch;
    end
    %actualNumCh=length(dirIndices);
    Iread(1).data=flipud(dlmread([PathName,D(dirIndices(1)).name],'\t',0,0)');
    sizeMeas=size(Iread(1).data);
    Iread(1).events=[];
    if length(dirIndices)<numch
        dirIndices=repmat(dirIndices,1,4);
    end
    if isfield(param,'differentialChannel')
        %channel to subtract from
    if any(param.differentialChannel) %if a bad channel exists set this array to zero
        Iread(param.differentialChannel).data=flipud(dlmread([PathName,D(dirIndices(param.differentialChannel)).name],'\t',0,0)');
       for ii=1:length(selch)
        Iread(selch(ii)).data=Iread(selch(ii)).data-Iread(param.differentialChannel).data;
       end
    end
    end

    for ii=1:length(selch)
        Iread(selch(ii)).data=flipud(dlmread([PathName,D(dirIndices(selch(ii))).name],'\t',0,0)');
        aa=find(~isnan(Iread(selch(ii)).data(:,1)));
        aaa=find(~isnan(Iread(selch(ii)).data(1,:)));
        aaaa=(Iread(selch(ii)).data(aa,aaa));
        Iread(selch(ii)).data=aaaa;
        Iread(selch(ii)).rawdata=aaaa;
    end
end
%custom multichannel fscv recordings
if isfield(fscvdata,'recordedData')
    fscvData=fscvdata.recordedData;
    sizefscvData=size(fscvData);
    numcols=sizefscvData(2);
    %numcols=numcols-1;              %discard last col (clk in)
    ts1=fscvData(1,1);      %1st timestamp
    %if first timestamp is > 1, then corresponds to new unpaused recording
    %system where timestampsa re accumulated, when loading want everything
    %relative to 0, so subtract htis
    if ts1>0.5
        fscvData(:,1)=fscvData(:,1)-round(fscvData(1,1));   
    end
    res=find(fscvData(:,2)<=-0.39);
    res=find(fscvData(:,1)<0.199);        %10 Hz cycle frequency so 200 ms corresponds to second cycle
    samplesperscan=length(res);
    firstCH=3;          %1st col = timestamp, 2nd col = filt ramp in, 3rd col = ch 1 (ai1), 19th col = clk in
    for ii=1:length(selch)
        %acquire data from selected channels as indicated in selch arg
        %reshape to map current across time (x-axis) w/ voltage distribution on y-axis
        arrayedData=reshape(fscvData(:,selch(ii)+firstCH-1),samplesperscan,length(fscvData)/samplesperscan);
        Iread(ii).data=flipud(arrayedData);
        Iread(ii).data=Iread(ii).data.*param.IVconv;      %convert to current based on gain setting IVconv
        Iread(ii).LPFdata=Iread(ii).data;
        if selch(ii)==param.badChannels 
            %if a bad channel exists set this array to zero
            Iread(ii).data=ones(size(Iread(ii).data));
            Iread(ii).LPFdata=ones(size(Iread(ii).data));
        end

    end
    if param.differentialChannel>0 
        %if differential ch exists
        for ii=1:length(selch)
            if selch(ii)~=param.differentialChannel
                Iread(ii).data=Iread(ii).data-Iread(param.differentialChannel).data;
                Iread(ii).LPFdata=Iread(ii).LPFdata-Iread(param.differentialChannel).LPFdata;
            end
        end
    end
end


end
    