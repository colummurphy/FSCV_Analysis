function menuselect(source,event)
global hgui plotParam parameters processed csc_map settings
val = source.Value;     %selected # setting on drop down list
listoptions = source.String;   %all strings on drop down list    
% For R2014a and earlier: 
% val = get(source,'Value');
% maps = get(source,'String'); 
disp(listoptions{val});
selch=plotParam.selch;
switch val
    case 2
        %select file
        [hgui.FileName,hgui.PathName] = uigetfile('*.*','Select file');
        cd(hgui.PathName)
        %get file info, global hgui to get info on file/path, update hgui
        getfileinfo;    
        cursess=hgui.cname;
        hgui.cname=[hgui.subject num2str(hgui.sessionnum) ];
        if ~isequal(hgui.cname,cursess)
            if ~isempty(hgui.ephysid)
                %ephys session
                hgui.cname=[hgui.subject num2str(hgui.sessionnum) hgui.ephysid];
            end
            ncschannels={};
            if strcmp(hgui.subject,'cleo')
                [parameters,csc_map,eventcodes]=getparams(hgui.cname,'cleo',ncschannels);
            else
                [parameters,csc_map,eventcodes]=getparams('patrabipolar','default',ncschannels,'sessnum',hgui.sessionnum);
            end
        end
        plotParam.sites=getfscvsitenames(hgui.PathName);
        processed={}; parameters.badzone=[];
        %load file
        plotParam.dir=dir(hgui.PathName);
        [processed.Iread,processed.LFPread,processed.samplesNCS]=...
            loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch,'dir',plotParam.dir);

        parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;
        %update plotParam (incase channels changed update)
        getplotsettings(plotParam.filtlfp,plotParam.cscNames,...
            plotParam.event_codes,parameters.sampleratencs,settings);       
        disp(['loading: ' hgui.FileName]);


        delete(hgui.title)
        %compile/organize data for plotting as loaded above in hgui
                if isempty(hgui.ephysid)
        compileloaded(hgui);    %update processed, parameters, plotParam globals
                end
        %plot ephys signals if exist
        if ~isempty(processed.samplesNCS) && ~isempty(plotParam.cscNames)
            if isfield(plotParam,'event_codes')
                processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
            end
            guiprolfp('getbursts');  %process ncs signals for plotting
            if ~isempty(plotParam.lfpid)
            setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
                        'sqenv','useprocessed','bursts');  
                    setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
                'scale',plotParam.LFPscale,'bursts');
            end
            setguicloseup(hgui.closeup{3},'ncsids',[plotParam.eyeid...
                plotParam.pulseid],...
                'norm',1,'xlab'); 
            setguicloseup(hgui.closeup{3},'ncsids',plotParam.lickid,...
                'norm',1,'xlab','env','winlength',plotParam.winlengthphys,'hold');  
            
            if plotParam.buttonm==1
                setspectrum(processed,plotParam,hgui.fftplot,plotParam.mcsc); 
            end
        end
        hgui.title=text(hgui.titletext,0, 0,...
            [hgui.trialtype ' | trial #: ' hgui.trialnum ' | session #: '...
            num2str(hgui.sessionnum)],...
            'interpreter','none');     
        
        
    case 3
        %load settings
        loadfile = uigetfile('*.*', 'Select settings file');
        load(loadfile)
        disp(['loaded ' loadfile ' settings']);  
        
        
    case 4
        %save settings
        savedir = uigetdir('*.*', 'Select directory to save plotParam settings');
        saveName=[savedir '\settings_' hgui.trialtype];
        save(saveName,'plotParam','parameters');
        disp(['saved ' saveName]);
        
        
    case 5
        %select ephys chs
        %update parameters.NCSchannels based on string labels in 
        %plotParam.cscNames (taken from ncschannels)
        %get ch selection from chronicXXchconfig file or cleo_chronicXX
        chmap=hgui.chmap;
        if ~isempty(chmap)
           setfscvguichannelselection(chmap,'selectchs');    %updates plotParam.cscNames
        else
            setfscvguichannelselection(csc_map);    %updates plotParam.cscNames
        end
        

        %next time files are loaded, will update ncs plotted channels
    case 6
        %change color scale of morletgram
        prompt = 'morletgram scale?';
        plotParam.fftclim=  input(prompt) ;
        cla(hgui.closeup{1});
        setspectrum(processed,plotParam,hgui.closeup{1},plotParam.mcsc);  
        
        
    case 7
        %detect da in selected time window only, whether correlated or not
        t1=plotParam.t_start;
        t2=plotParam.t_end;
        plotParam.timeWin=([t1 t2]-1)...
            ./parameters.samplerate;    %window in seconds
        cla(hgui.itplot)
        
        %global processed, parameters, plotParam
        %run findsigartifacts.m, findglitches.m, & detectdatransients.m
        %updates processed.glitchids also
        getdasignals([t1 t2],'uncorr');
            
        for ii=1:length(selch)
            %refresh pca computed concentrations w/ detection marks
            setguipct(hgui.itplot,processed.Ipcr{selch(ii)},plotParam, ...
                parameters,plotParam.colorFSCV(ii-4*(ceil(ii/4)-1),:),...
                'detected',round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -t1+2),'plotnum',ii);       
        end
        
    case 8
        %detect da for entire recording loaded, whether correlated or not
        t1=1;
        t2=size(processed.Iread(selch(1)).rawdata,2);   %entire length
        plotParam.timeWin=([t1 t2]-1)...
            ./parameters.samplerate;    %window in seconds
        
        %global processed, parameters, plotParam
        %run findsigartifacts.m, findglitches.m, & detectdatransients.m
        getdasignals([t1 t2],'uncorr');
        cla(hgui.itplot)

        for ii=1:length(selch)
            plotax=hgui.itplot;  
                if selch(ii)<=8 && selch(ii)>4  
                    plotax=hgui.itplotx{1};
                elseif selch(ii)<=12 && selch(ii)>8
                    plotax=hgui.itplotx{2};
                elseif selch(ii)<=16 && selch(ii)>12
                    plotax=hgui.itplotx{3};
                end
            %refresh pca computed concentrations w/ detection marks
                setguipct(plotax,processed.Ipcr{selch(ii)},plotParam, ...
                parameters,plotParam.colorFSCV((selch(ii)-4*(ceil(selch(ii)/4)-1)),:),...
                'detected',round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -plotParam.t_start+2),'plotnum',ii);                  
        end
        
        disp('finished detection')
        
    case 9
        %get cross covariance of signals within window
        %for each da signal, get xcov with all ephys channels selected
        t1=1;
        t2=size(processed.Iread(selch(1)).rawdata,2);   %entire length
        plotParam.timeWin=([t1 t2]-1)...
            ./parameters.samplerate;    %window in seconds
        
        %global processed, parameters, plotParam
        %run findsigartifacts.m, findglitches.m, & detectdatransients.m
        getdasignals([t1 t2],'xcov');
        
        disp('finished detection')     
        
    case 10
        %get xcov for entire recording loaded, whether correlated or not
        t1=1;
        t2=size(processed.Iread(selch(1)).rawdata,2);   %entire length
        plotParam.timeWin=([t1 t2]-1)...
            ./parameters.samplerate;    %window in seconds
        
        %global processed, parameters, plotParam
        %run findsigartifacts.m, findglitches.m, & detectdatransients.m
        getdasignals([t1 t2],'xcov');
        
        disp('finished detection')
   
    case 11
        %merge files, multi-channel system only working
        [filenames, pathname] = uigetfile( ...
       {'*.mat','MAT-files (*.mat)'; ...
        '*.mdl','Models (*.mdl)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Select files to merge, must be selected in order of the serial recording, will NOT re-sort', ...
        'MultiSelect', 'on');
        firstFileName=filenames{1};
        pathdir=pathname;
        filenums=[];
        nfiles=length(filenames);
        tempdata=[];
        for ifile=1:nfiles
            load([pathdir filenames{ifile}]);       %loads 'recordedData' variable
            disp(['merging file ' filenames{ifile}]);
            tempdata=[tempdata; recordedData];                  
        end
        savedrive=[pathname 'merged' filesep];
        if ~isdir(savedrive)
            mkdir(savedrive);
        end
        recordedData=tempdata;
         numid=strfind(filenames{1},'_');
         numid2=strfind(filenames{1},'.mat');
         numid3a=strfind(filenames{end},'_');
         numid3=strfind(filenames{end},'.mat');
        save([savedrive filenames{1}(1:numid2-1) '-' filenames{end}(numid3a(end)+1:numid3-1) '_merged'],'recordedData','filenames','pathname','-v7.3');   
        
   case 12
        %detect da for entire directory    %window in seconds        
        %run scrips file by file
        %user selects files and these are sorted
        %assumes parameters are the ones already loaded before loading new
        %files        
        [filenames, pathname, filterindex] = uigetfile( ...
       {'*.mat','MAT-files (*.mat)'; ...
        '*.mdl','Models (*.mdl)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Select files to detect da', ...
        'MultiSelect', 'on');
        files=dir(pathname);
        ispresent=contains(filenames,'_','IgnoreCase',true);       
        firstFileName=filenames{1};
        pathdir=pathname;
        hgui.PathName=pathdir;
        cd(hgui.PathName);
        filenums=[];
        nfiles=length(filenames);
        if nfiles>1
        for ifile=1:nfiles
            numid=strfind(filenames{ifile},'_');
            numid2=strfind(filenames{ifile},'.mat');
            filenums(ifile)=str2num(filenames{ifile}(numid(end)+1:numid2-1));
        end
        if ~issorted(filenums)
            %sort by date
            [xx,sortids]=sort(filenums);
            filenamesorted=filenames(sortids);
            filenames=filenamesorted;
        end       
        end
        for ifile=1:nfiles
            hgui.FileName=filenames{ifile};
            processed={};
            [processed.Iread,processed.LFPread,processed.samplesNCS]=...
                loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch);
            disp(['loading: ' hgui.FileName]);
            compileloaded(hgui);    %update processed, parameters, plotParam globals

            t1=1;
            t2=size(processed.Iread(selch(1)).rawdata,2);   %entire length
            plotParam.timeWin=([t1 t2]-1)...
                ./parameters.samplerate;        
            getdasignals([t1 t2],'uncorr');
            
            detected=processed.detected;
            for ii=1:length(selch)
                saveName=['da_ch' num2str(selch(ii)) '_' hgui.FileName];
                tlab=['_tracewin_' num2str(parameters.tracepad(2)) 's'];
                savedir=[hgui.PathName 'auto' tlab filesep 'ch' num2str(selch(ii)) '\'];
                if ~isdir(savedir)
                    mkdir(savedir);   
                end
                da=detected{selch(ii)};
                if ~isempty(da.datrace)
                disp(['ch ' num2str(selch(ii)) ' detected ' num2str(length(da.maxDA)) ...
                    ' da signals ' ]);
                else
                    disp(['ch ' num2str(selch(ii)) ' nothing detected'])
                end
                       
                %save process settings / parameters/plotParam & detected
                if isfield(hgui,'session')
                    parameters.session=hgui.session;
                    parameters.date=hgui.date;
                end
                events.LFPeventTS=processed.LFPread.LFPeventTS;
                events.LFPeventTTL=processed.LFPread.LFPeventTTL;
                save([savedir saveName],'parameters','plotParam','events','da','-v7.3');   
                disp(['saved ' saveName]);
            end
        end
        disp('finished detection')
        
      case 13
        %detect da get xcov save for entire directory
    %window in seconds        
        %run scrips file by file
        %user selects files and these are sorted
        %assumes parameters are the ones already loaded before loading new
        %files
        files = uigetdir2(pwd,'Select successive files for det');
        nfiles=size(files,2);
        sep=findstr(filesep,files{1});
        ispresent=contains(files,'fscv_multi_','IgnoreCase',true);
        gfiles=find(ispresent>0);         %find first rec files
        if isempty(gfiles)
            error('fscv_multi_ format files required');
        end
        files=files(gfiles);
        firstFileName=files{1}(sep(end)+1:end);
        pathdir=files{1}(1:sep(end));
        hgui.PathName=pathdir;
        currentdir=dir(hgui.PathName);
        arrayedFiles =vertcat(files{:});        %only works when length of strings same
        arrayedNames=arrayedFiles(:, sep(end)+1:end);
        endNameIdx=regexpi(arrayedNames(1,:),'.mat');
        filenums=arrayedNames(:,endNameIdx-3:endNameIdx-1);
        %file #'s of continuous fscv recordings, should start from 100 typically
        fileNums2=str2num(filenums);         
        %sort filenames
        [sorted,sortID]=sort(fileNums2,'ascend');
        sortedfiles=files(sortID);
        sortednums=filenums(sortID,:);
        sortednames=arrayedNames(sortID,:);
        %check if already sorted (not really needed..)
        xx=sorted-fileNums2;
        if sum(xx)==0
            display('names already sorted as loaded');
        end
        for ifile=1:nfiles
            hgui.FileName=sortednames(ifile,:);
            processed={};
            [processed.Iread,processed.LFPread,processed.samplesNCS]=...
                loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch);
            disp(['loading: ' hgui.FileName]);
            compileloaded(hgui);    %update processed, parameters, plotParam globals

            t1=1;
            t2=size(processed.Iread(selch(1)).rawdata,2);   %entire length
            plotParam.timeWin=([t1 t2]-1)...
                ./parameters.samplerate;        
            getdasignals([t1 t2],'uncorr');
            getdasignals([t1 t2],'xcov');
            
            detected=processed.detected;
            for ii=1:length(selch)
                saveName=['da_ch' num2str(selch(ii)) '_' hgui.FileName];
                tlab=['_tracewin_' num2str(parameters.tracepad(2)) 's'];
                savedir=[hgui.PathName 'auto' tlab filesep 'ch' num2str(selch(ii)) '\'];
                if ~isdir(savedir)
                    mkdir(savedir);   
                end
                da=detected{selch(ii)};
                if ~isempty(da.datrace)
                disp(['ch ' num2str(selch(ii)) ' detected ' num2str(length(da.maxDA)) ...
                    ' da signals ' ]);
                else
                    disp(['ch ' num2str(selch(ii)) ' nothing detected'])
                end
                       
                %save process settings / parameters/plotParam & detected
                if isfield(hgui,'session')
                    parameters.session=hgui.session;
                    parameters.date=hgui.date;
                end
                events.LFPeventTS=processed.LFPread.LFPeventTS;
                events.LFPeventTTL=processed.LFPread.LFPeventTTL;
                save([savedir saveName],'parameters','plotParam','events','da','-v7.3');   
                disp(['saved ' saveName]);
            end
        end
        disp('finished detection')
        
    case 14
        %load detected da files in auto_ dir for corresponding loaded file
        %ie _100 or _101.. etc. to modify or check detected ts's with color
        [dapath] = uigetdir('*.*','Select auto dir with detected signals');
        %find other ch's in back dir        
        dirpath=dir(dapath);
        for ii=1:length(selch)
            chf=find(ismember({dirpath.name},['ch' num2str(selch(ii))])==1);
            if ~isempty(chf)
                %load selected da auto detect file
                %in format da_chXX_fscv_multi_1XX
                fileload=[dirpath(chf).name...
                    filesep 'da_ch' num2str(selch(ii)) '_' ...
                    hgui.FileName];
                pathload=[dirpath(chf).folder filesep ];
                load([pathload fileload],'da'); 
                %insert into processed.detected variable
                processed.detected{selch(ii)}=da;        
                disp(['loading: ' fileload])
                if ii==1
                    %load parameters and badzone ids..
                    paramtemp=load([pathload fileload],'parameters');
                    processed.badids=paramtemp.parameters.badzone;
                    parameters.badzone=paramtemp.parameters.badzone;
                    disp(['loaded parameters bad zone, stored badzone in processed badids'])
                end
                hgui.loaded{selch(ii)}=[pathload fileload];
            end
        end

    case 15
        %save detected da transients in processed winwo
        %ie detected parameters (processed.detected) &
        %for defined window around each transient, save:
        %pca da traces, & if included:
        %lfp signals selected, filtered signals selected,
        %events in period
        detected=processed.detected;
        tracepad=parameters.tracepad;
        bf=[];
        tempload={};
        da={};
        for ii=1:length(selch)
            saveName=['da_ch' num2str(selch(ii)) '_' hgui.FileName];
            tlab=['_tracewin_' num2str(parameters.tracepad(2)) 's'];
            savedir=[hgui.PathName 'auto' tlab filesep 'ch' num2str(selch(ii)) filesep];
            da=detected{selch(ii)};
            if ~isdir(savedir)
                mkdir(savedir);   
            else
                if exist([savedir saveName])>0
                    %file already exists
                    if isempty(bf)
                        %query to replace or not parameters
                        bf=input(['1. replace existing \n'...
                            '2. replace detected da and bad ids only (keep nlx chs/etc same) \n'...
                            '3. keep existing bad ids/parameters & replace detected da only \n']...
                        ,'s') ;
                    end
                    %user already provided, use uniform for all chs
                    switch bf
                        case '1'
                            %don't do anything, ie don't load anything
                            %existing
                        case '2'
                            %only replace detected da & bad ids 
                            %bad ids in parameters
                            %da in da
                            tempload=load(hgui.loaded{selch(ii)},'plotParam','parameters','da');
                            %keep existing plot param
                            plotParam=tempload.plotParam;
                            %keep existing parameters/bad ids
                            newbad=parameters.badzone;
                            %load existing parameters
                            parameters=tempload.parameters;
                            %replace only bad ids
                            parameters.badzone=newbad;
                            newTS=da.maxTS;
                            %keep existing da file
                            da=tempload.da;
                            %check that existing length of detected da matches new
                            if length(newTS)~=length(da.maxTS)
                                error('# of detected da does not match existing file')
                            end
                            %just replace maxTS
                            da.maxTS=newTS;                            
                        case '3'
                            %only replace existing da.maxTS
                            %tempload=load(hgui.loaded{selch(ii)},'plotParam','parameters','da');
                            
                            %plotParam=tempload.plotParam;
                            tempload=load(hgui.loaded{selch(ii)},'parameters','da');
                            %keep existing parameters/bad ids
                            parameters=tempload.parameters;
                            newTS=da.maxTS;
                            %load existing da data
                            da=tempload.da;
                            %check that existing length of detected da matches new
                            if length(newTS)~=length(da.maxTS)
                                error('# of detected da does not match existing file')
                            end
                            %replace existing detected timestamps only
                            da.maxTS=newTS;
                    end
                end
            end            
            
            %save process settings / parameters/plotParam & detected
            if isfield(hgui,'session')
                parameters.session=hgui.session;
                if isfield(hgui,'date')
                parameters.date=hgui.date;
                end
            end
            %get nlx events info
            events.LFPeventTS=processed.LFPread.LFPeventTS;
            events.LFPeventTTL=processed.LFPread.LFPeventTTL;
            save([savedir saveName],'parameters','plotParam','events','da','-v7.3');   
            disp(['saved ' saveName]);
        end
        
    case 16
        %load bad timestamps
        loadfile = uigetfile('*.*', 'Select badzone file');
        load(loadfile)
        disp(['loaded ' loadfile ' bad ids in badts']); 
        tempbadids=[];
        %merge with existing bad ids
        if isfield(processed,'badids')
            tempbadids=processed.badids;
        end    
        disp(['merge with current processed and parameters']);
        processed.badids=sort(unique([tempbadids badts]));
        parameters.badzone=processed.badids;
        disp(['nanning detected peaks around bad zone for sel fscv chs']);
        %NAN any detected peaks around badzones, just maxTS variable
        if isfield(processed,'detected')
            badoverlap=[10 25];      %after these samples of movement artifact do not search
            if isfield(parameters,'badoverlap')
                badoverlap=parameters.badoverlap;
            end
            for ii=1:length(plotParam.selch)
                ich=plotParam.selch(ii);
                maxts=round((processed.detected{ich}.maxTS-...
                    processed.LFPread.LFPts(1)).*parameters.samplerate...
                    +1);
                for ip=1:length(maxts)
                    if any(ismember(maxts(ip)-badoverlap(2):maxts(ip)+badoverlap(1),...
                            parameters.badzone))
                        processed.detected{ich}.maxTS(ip)=nan;
                    end                    
                end
            end
        end
        
    case 17   
        %save bad timestamps
        badts=parameters.badzone;
        if size(badts,2)<size(badts,1)
            %if not row vector make row vector first
            badts=badts';
        end
        filenum=hgui.trialnum;
        saveName=['badzone_' filenum];
        do='1';
        if exist([hgui.PathName saveName '.mat'])>0
            do=input(['file exists \n' ...
            '1. replace \n'...
            '2. merge \n'...
            '3. new file \n'],'s');
        end
        switch do
            case '1'
                %replace
                save([hgui.PathName saveName],'badts');  
                        disp(['saved ' saveName]);       
            case '2'
                %merge
                badtstemp=badts;
                %load stored
                load([hgui.PathName saveName],'badts');
                badts=[badts badtstemp];   %merge
                badts=unique(sort(badts));  %sort
                save([hgui.PathName saveName],'badts');     
                        disp(['saved ' saveName]);       
            case '3'
                %new file
                saveName=input(['new file name: \n'],'s');
                save([hgui.PathName saveName],'badts');   
                        disp(['saved ' saveName]);       
            otherwise
                disp('bad input');
        end
            
    case 18
        %save bad trials log text files to excel sheet in triallogs folder
        savebadtrials(hgui.subject,hgui.sessionnum);

            
end

        set(source,'Value',1);

end