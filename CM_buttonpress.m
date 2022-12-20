function []= CM_buttonpress(key, appDesColorFig, ...
        appDesDopaminePlot, appDesTitleText)

global plotParam processed hgui parameters settings
%functions to run when button pressed

 
sessid=hgui.sessionnum;
if isfield(hgui,'ephysid')
    if ~isempty(hgui.ephysid)
        sessid=[num2str(hgui.sessionnum) hgui.ephysid];
    end
end

switch key
    case 'return'
        %same function as refresh button
        refreshbutton
        
    case {'leftarrow','rightarrow'}
    %load previous/next file
        
        % clear processed data
        processed={};
        
        % clear title text
        appDesTitleText.Value = '';

        if strcmpi(key,'leftarrow')
            % left - decrement index
            hgui.currentIndex=hgui.currentIndex-1;
        else
            % right - increment index
            hgui.currentIndex=hgui.currentIndex+1;
        end

        % if currentIndex is outside the index of files
        if hgui.currentIndex>=length(hgui.dirIndex) || hgui.currentIndex<1
            % set it to the first file
            hgui.currentIndex=1; 
        end
        
        % get the FileName at the new index
        hgui.FileName=hgui.dirIndex(hgui.currentIndex).name; 
        

        % get the index of '.mat'
        namx=strfind(hgui.FileName,'.mat');
        

        % if '.mat' found get the trial number eg. 101 
        if ~isempty(namx)
            hgui.trialnum=hgui.FileName(namx-3:namx-1);
        else
            hgui.trialnum=[];
        end

        % if '.mat' not found, trial type = File Name
        if isempty(namx)
            %not behavior file, maybe raw data recorded so no name
            hgui.trialtype=hgui.FileName;
        end

        % Returns vector of indexes, each index - location of '/' 
        namy=strfind(hgui.PathName, ("/"|"\"));
        
        % Return index of chronic
        s1=strfind(hgui.PathName,'chronic');
        
        % Returns vector of indexes, each index - location of '_' 
        ss1=strfind(hgui.PathName,'_');
        
        % if pathname contains 'chronic', s1 not empty
        if ~isempty(s1)
            % Vector of index/indexes, of '_', which occur after 'chronic'
            s2=ss1(ss1>s1);
            % session string set from 'chronic' to '_' eg. chronic125 
            hgui.session=hgui.PathName(s1:s2(1)-1);
            % get the session number from the session string
            hgui.sessionnum=str2num(hgui.session(regexpi(hgui.session,'\d')));
            
            % Vector of index/indexes of '/', which occur after 
            % the '_' after 'chronic'
            s3=namy(namy>s2(1));
            
            % get the date
            hgui.date=hgui.PathName(s2(1)+1:s3(1)-1);
        end

        % load all ephys and fscv data
        [processed.Iread,processed.LFPread,processed.samplesNCS]=...
            loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch,'dir',plotParam.dir);
        
        % parameters.samplingratencs = sampling rate from fscv file
        parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;
        
        %update plotParam (incase channels changed update)
        getplotsettings(plotParam.filtlfp,plotParam.cscNames,...
            plotParam.event_codes,parameters.sampleratencs,settings);   
        
        disp(['loading: ' hgui.FileName]);
        
        % compile/organize data for plotting as loaded above in hgui
        if isempty(hgui.ephysid)
            CM_compileloaded(hgui, appDesColorFig, appDesDopaminePlot); 
        end
        
        %plot ephys signals if exist
        % CM - Comment out ephys code 
        %{
        if ~isempty(processed.samplesNCS) && ~isempty(plotParam.cscNames)
            if isfield(plotParam,'event_codes')
                processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
            end
            guiprolfp('getbursts');  %process ncs signals for plotting
            if ~isempty(plotParam.lfpid)
            setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
                        'sqenv','useprocessed','bursts','scale',plotParam.powerscale);   
            setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
                'scale',plotParam.LFPscale, 'bursts');
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
        %}

        % update title text
        appDesTitleText.Value = [ hgui.trialtype ' | trial #: '...
            hgui.trialnum ' | session #: ' num2str(hgui.sessionnum) ];


% CM - Commented Out        
%{
    case 's'
        %save checked color plot data
        savechs=[];        
        %get common data across all channels
        nlx=plotParam;
        events=plotParam.events;
        ncsread=processed.LFPread;
        %filter ncs channels based on eye/lick/lfp categories
        groups{1}=plotParam.lfpid;      %chs for lfps
        groups{2}=[plotParam.eyeid plotParam.pulseid];  %do not filt, already filted
        groups{3}=plotParam.lickid;
        fgroups{1}=plotParam.filtlfp;
        fgroups{3}=plotParam.filtlick;

        %get lfp signals
        for ich=1:length(processed.LFPread.LFPchNames)
            %lfp signal --> filter, square, envelope
            samplesncs=processed.samplesNCS(:,ich)';
            if ismember(ich,plotParam.lfpid)
                %filter at "beta-band" as defined in nlx.filtlfp
                %default filter is filtbetah now                
                nlx.resampled(ich,:)=samplesncs;   %store unfiltered data
                    %if another filter band specified 07/07/18
                    nlx.resampledb(ich,:)=filterLFP(samplesncs,plotParam.ratelfp,plotParam.filtlfp);
                    %square & envelope signal
                    winlength=round(plotParam.ratelfp*.5/mean(plotParam.filtlfp));
                    if isfield(plotParam,'winlength')
                    %if window already defined in parameters load value
                    winlength=round(plotParam.ratelfp*plotParam.winlength);
                    end
                    nlx.resampledb(ich,:)=nlx.resampledb(ich,:).^2;   %get power V^2
                    nlx.resampledb(ich,:)=smoothwin(nlx.resampledb(ich,:),winlength);   %smoothing 
            %lick signal --> filter 
            elseif ismember(ich,plotParam.lickid)
                nlx.resampled(ich,:)=filterLFP(samplesncs,plotParam.ratelfp,plotParam.filtlick);
                winlength=round(plotParam.ratelfp*.5/mean(plotParam.filtlick));
                if isfield(plotParam,'winlengthphys')
                    %if window already defined in parameters load value
                    winlength=round(plotParam.ratelfp*plotParam.winlengthphys);
                end
                nlx.resampledb(ich,:)=smoothwin(nlx.resampledb(ich,:),winlength);   %smoothing

            %eye/pulse signal --> dont do anything already filt in reconvertFSCV.m            
            else
                nlx.resampled(ich,:)=samplesncs;
            end
        end
       % nlx.resampled=filterncssignals(processed.samplesNCS,...
        %    parameters.sampleratencs,groups,fgroups,'envgroups',[1 0 0]);
        
        %create output path
        for iii=1:length(plotParam.selch)
            %hgui.PathNameOut{iii}=[hgui.PathName 'out\ch' num2str(plotParam.selch(iii)) '\'];
                        hgui.PathNameOut{iii}=[hgui.PathName 'out\'];
            if ~isdir(hgui.PathNameOut{iii})
                status = mkdir(hgui.PathNameOut{iii});
            end
            hgui.PathNameOutSingles{iii}=[hgui.PathName 'out\singles\'];
            if ~isdir(hgui.PathNameOutSingles{iii})
                status = mkdir(hgui.PathNameOutSingles{iii});
            end
        end
        tParam.timeWin=[20 40]; 
        infobehav={};
        if ~contains(hgui.subject,'cfmea')
        infobehav=calcBehav(processed,tParam.timeWin,plotParam.event_codes);   
        end
        %save selected channels
        cvdata=[];
        for iii=1:length(plotParam.selch)
            chnum=plotParam.selch(iii);
            saveName=['CVITdata_ch' num2str(chnum) '_' hgui.FileName];
            processed.info(iii).behav=infobehav;    

            Ipcr=processed.Ipcr{chnum};        
            Idata.anodal=processed.cv{chnum}(1:length(parameters.Vrange));
            Idata.cathodal=processed.cv{chnum}(length(parameters.Vrange)+1:end);
            cvdata(chnum,:)=[Idata.anodal; Idata.cathodal];
            BG=processed.BG(chnum).data;
            Vrange=(parameters.Vrange); Vrange_cathodal=parameters.Vrange_cathodal;
            Isub=processed.Isub(chnum).data;
            info=processed.info(iii);
            save([hgui.PathNameOutSingles{iii} saveName],'Ipcr','nlx','ncsread','parameters','Idata','BG',...
                'Vrange','Vrange_cathodal','Isub','info','events','chnum');
        end
            Ipcr=processed.Ipcr;
            BG=processed.BG;
            Isub=processed.Isub;
            info=processed.info;
            saveName=['out_' hgui.FileName];
            save([hgui.PathNameOut{iii} saveName],'Ipcr','cvdata','parameters','BG',...
                'Vrange','Vrange_cathodal','Isub','info','events');     %merged save
            
    case 'f'
        %plot morlet
        prompt = 'morletgram channel?';
        plotParam.mcsc = hgui.mcsc;
        cla(hgui.fftplot);
        cla(hgui.data);
        cla(hgui.cv);
        set(hgui.cv, 'visible','off', 'box','off','xtick',[],'ytick',[]);  %make invisble until prompted
        set(hgui.data, 'visible','off');  %make invisble until prompted
        setspectrum(processed,plotParam,hgui.fftplot,plotParam.mcsc); 
        plotParam.buttonm=1;
        set(hgui.checkfft,'Value',1);
        
     case 'g'
        %no morlet
      % cla(hgui.closeup{1},'reset')
        set(hgui.cv, 'visible','on');  %make invisble until prompted
        set(hgui.data, 'visible','on');  %make invisble until prompted
        cla(hgui.fftplot,'reset')
        set(hgui.fftplot, 'visible','off', 'box','off','xtick',[],'ytick',[]);  %make invisble until prompted
        plotParam.buttonm=0; 
        set(hgui.checkfft,'Value',0);
        
     case 'm'
        %identify movement artifact periods interactively
        %on currently clicked plot
        axes(hgui.itplot);
        axsel=hgui.itplot;
        if length(plotParam.selch)>4
            axsel=gca;
            axes(axsel);            
            itplotpos=hgui.itplotpos;
            cp = get(gcf,'CurrentPoint');
            if cp(1)>=itplotpos{2}(1) && cp(1)<=itplotpos{2}(1)+itplotpos{2}(3)
                axsel=hgui.itplotx{1};
            elseif cp(1)>=itplotpos{3}(1) && cp(1)<=itplotpos{3}(1)+itplotpos{3}(3)
                axsel=hgui.itplotx{2};
            elseif cp(1)>=itplotpos{4}(1) && cp(1)<=itplotpos{4}(1)+itplotpos{4}(3)
                axsel=hgui.itplotx{3};
            else
                axsel=hgui.itplot;
            end
            axes(axsel);
        end
        bf='';
        x=[];y=[];
        bids=[];
        tempbadids=processed.badids;
        movsel=getrect(axsel); %[xmin ymin width height];    
        

        badts(1)=round(movsel(1))+plotParam.t_start;
        badts(2)=round(movsel(1)+movsel(3))+plotParam.t_start;
        disp(['bad ts sel: ' num2str(badts(1)/10) ' - ' num2str(badts(2)/10) ' s'])
        bids=[badts(1):badts(2)];
        if ~isempty(tempbadids) 
            %already stored bad ids
            ovids=intersect(tempbadids,bids);
            if ~isempty(ovids)
                %if already identified, remove these, inverse operation
                nib=find(ismember(tempbadids,ovids)==1); 
                %for ib=1:length(nib)
                %processed.badids(nib)=[];
                tempbadids(nib)=[];
                    %delete(hgui.bad{nib(ib)});  
                %end
               % parameters.badzone=processed.badids;
               processed.badids=tempbadids;
               parameters.badzone=tempbadids;
               togglem(hgui.badplot,[]);        %display bad points
            else
                %new ids not in tempbadids
                if size(tempbadids,2)>size(tempbadids,1)
                    tempbadids=[tempbadids bids];
                else
                    %make row vector
                    %flip matrix..somehow got flipped somewhere..need fix
                    tempbadids=tempbadids';
                    tempbadids=[tempbadids bids];
                end
                processed.badids=tempbadids;
                       %parameters.badzone=processed.badids;
                parameters.badzone=tempbadids;
                togglem(hgui.badplot,[]);   
               %{
                togglem( hgui.badplot,[]);
                if hgui.badplot.Value==1
                    nib=1;
                    if isfield(hgui,'bad')
                        nib=length(hgui.bad)+1;
                    end
                    yp=get(hgui.itplot,'ylim');
                    hgui.bad{nib}=text(hgui.itplot,xsel-plotParam.t_start,yp(2),'x','fontsize',12,'horizontalalignment',...
                        'center','color',[1 0 0]);
                end
                %}
            end                
        else
            tempbadids=bids;
            parameters.badzone=tempbadids;
               %     parameters.badzone=processed.badids;
            togglem( hgui.badplot,[]);
        end                    
        %store bad ids into parameters also
        processed.badids=sort(unique(tempbadids));
        parameters.badzone=processed.badids;
        %NAN any detected peaks around badzones, just maxTS variable
        if isfield(processed,'detected')
            badoverlap=[1 1];      %after these samples of movement artifact do not search
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
%}
% CM - End of Temp Comment Out        

   case 'b'
       %bad trial all channels, log       
       logbadtrials(hgui.subject,sessid,hgui.trialtype,str2num(hgui.trialnum),1:4)
    case 'l'
       %bad trial all channels, log       ephys
       logbadtrials(hgui.subject,[num2str(sessid) 'l'],hgui.trialtype,str2num(hgui.trialnum),1)
   case 'u'
       %bad trial for channel 1
       logbadtrials(hgui.subject,hgui.sessionnum,hgui.trialtype,str2num(hgui.trialnum),1)
   case 'i'
       %bad trial for channel 2
       logbadtrials(hgui.subject,hgui.sessionnum,hgui.trialtype,str2num(hgui.trialnum),2)
   case 'o'
       %bad trial for channel 3
       logbadtrials(hgui.subject,hgui.sessionnum,hgui.trialtype,str2num(hgui.trialnum),3)
   case 'p'
       %bad trial for channel 4
       logbadtrials(hgui.subject,hgui.sessionnum,hgui.trialtype,str2num(hgui.trialnum),4)


    otherwise
        warning(['unknown key: ' key]);
    
end

key=[];

end