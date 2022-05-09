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
        plotParam.sites=getfscvsitenames(hgui.PathName);
        processed={}; parameters.badzone=[];
        %load file
        [processed.Iread,processed.LFPread,processed.samplesNCS]=...
            loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch);
        parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;
        %update plotParam (incase channels changed update)
        getplotsettings(plotParam.filtlfp,plotParam.cscNames,...
            plotParam.event_codes,parameters.sampleratencs,settings);       
        disp(['loading: ' hgui.FileName]);
        delete(hgui.title)
        %compile/organize data for plotting as loaded above in hgui
        compileloaded(hgui);    %update processed, parameters, plotParam globals
        %plot ephys signals if exist
        if ~isempty(processed.samplesNCS)
            if isfield(plotParam,'event_codes')
                processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
            end
            %setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
           %     'filt',[16 33],'sqenv','units','env-filt z score',...
            %    'scale',plotParam.powerscale,'norm',1,'winlength',plotParam.winlength);  
            setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
                'filt',plotParam.filtlfp,'sqenv',...
                'scale',plotParam.powerscale,'winlength',plotParam.winlength);  
            setguicloseup(hgui.closeup{3},'ncsids',[plotParam.eyeid...
                plotParam.pulseid],'norm',1,'xlab'); 
            setguicloseup(hgui.closeup{3},'ncsids',plotParam.lickid,...
                'norm',1,'xlab','env','winlength',plotParam.winlengthphys,'hold');  
            if plotParam.buttonm==0
                setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
                    'scale',plotParam.LFPscale);
            else
                setspectrum(processed,plotParam,hgui.closeup{1},plotParam.mcsc); 
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
        setfscvguichannelselection(csc_map);    %updates plotParam.cscNames
        

        %next time files are loaded, will update ncs plotted channels
    case 6
        %change color scale of morletgram
        prompt = 'morletgram scale?';
        plotParam.fftclim=  input(prompt) ;
        cla(hgui.closeup{1});
        setspectrum(processed,plotParam,hgui.closeup{1},plotParam.mcsc);  
        
        
    case 7
        %detect da in selected time window only
        t1=plotParam.t_start;
        t2=plotParam.t_end;
        plotParam.timeWin=([t1 t2]-1)...
            ./parameters.samplerate;    %window in seconds
        cla(hgui.itplot)
        
        %global processed, parameters, plotParam
        %run findsigartifacts.m, findglitches.m, & detectdatransients.m
        %updates processed.glitchids also
        getdasignals([t1 t2],'corr');
            
        for ii=1:length(selch)
            %refresh pca computed concentrations w/ detection marks
            setguipct(hgui.itplot,processed.Ipcr{selch(ii)},plotParam, ...
                parameters,plotParam.colorFSCV(selch(ii),:),...
                'detected',round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -t1+2),'plotnum',ii);       
        end
        
        
  
    case 8
        %detect da for entire recording loaded
        t1=1;
        t2=size(processed.Iread(selch(1)).rawdata,2);   %entire length
        plotParam.timeWin=([t1 t2]-1)...
            ./parameters.samplerate;    %window in seconds
        
        %global processed, parameters, plotParam
        %run findsigartifacts.m, findglitches.m, & detectdatransients.m
        getdasignals([t1 t2],'corr');
        cla(hgui.itplot)

        for ii=1:length(selch)
                        
            %refresh pca computed concentrations w/ detection marks
            setguipct(hgui.itplot,processed.Ipcr{selch(ii)},plotParam, ...
                parameters,plotParam.colorFSCV(selch(ii),:),...
                'detected',round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -plotParam.t_start+2),'plotnum',ii);   

        end
        
        disp('finished detection')
        
    case 9
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
                parameters,plotParam.colorFSCV(selch(ii),:),...
                'detected',round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -t1+2),'plotnum',ii);       
        end
        
    case 10
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
                        
            %refresh pca computed concentrations w/ detection marks
            setguipct(hgui.itplot,processed.Ipcr{selch(ii)},plotParam, ...
                parameters,plotParam.colorFSCV(selch(ii),:),...
                'detected',round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -plotParam.t_start+2),'plotnum',ii);   

        end
        
        disp('finished detection')
        
    case 11
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
        
    case 12
        %get xcov for entire recording loaded, whether correlated or not
        t1=1;
        t2=size(processed.Iread(selch(1)).rawdata,2);   %entire length
        plotParam.timeWin=([t1 t2]-1)...
            ./parameters.samplerate;    %window in seconds
        
        %global processed, parameters, plotParam
        %run findsigartifacts.m, findglitches.m, & detectdatransients.m
        getdasignals([t1 t2],'xcov');
        
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
                    ' da signals & ' num2str(length(da.xcov{1}.sel)) ...
                    ' valid xcov traces per site pair' ]);
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
                    load([pathload fileload],'parameters');
                    processed.badids=parameters.badzone;
                    disp(['loaded parameters, stored badzone in processed badids'])
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
        %{
        if ~isempty(processed.samplesNCS)
            %get events/lfp signals & other ncs signals filtered/enveloped
            if isfield(processed.detected{selch(ii)},'maxTS')
                if ~isempty(processed.detected{selch(ii)}.maxTS)
                for ii=1:length(selch)
                    disp(['processing ch ' num2str(selch(ii))]);
                    %global processed plotParam parameters
                    detected{selch(ii)}.nlx=reconvertncs([],...
                        processed.detected{selch(ii)}.maxTS,tracepad);
                end
                end
            end
        end
        %}
        for ii=1:length(selch)
            saveName=['da_ch' num2str(selch(ii)) '_' hgui.FileName];
            tlab=['_tracewin_' num2str(parameters.tracepad(2)) 's'];
            savedir=[hgui.PathName 'auto' tlab filesep 'ch' num2str(selch(ii)) '\'];
            if ~isdir(savedir)
            mkdir(savedir);   
            else
                if exist([savedir saveName])>0
                    %file already exists
                    bf=input(['ch ' num2str(selch(ii)) ...
                    ': replace stored parameters with current settings?']...
                    ,'s') ;
                    if strcmpi(bf,'n') || strcmpi(bf,'no')
                        %load previously saved parameters so that ncs channels
                        %are same and just replace bad ids/detected
                        load(hgui.loaded{selch(ii)},'parameters','plotParam');
                    end
                end
            end

            da=detected{selch(ii)};
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
        
    case 16 
        %load bad timestamps
        loadfile = uigetfile('*.*', 'Select badzone file');
        load(loadfile)
        disp(['loaded ' loadfile ' bad ids']);   
        
    case 17   
        %save bad timestamps
        badts=parameters.badzone;
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
                badts=[badts; badtstemp];   %merge
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
            

            
end

        set(source,'Value',1);

end