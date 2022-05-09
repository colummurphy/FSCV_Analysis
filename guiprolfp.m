function guiprolfp(varargin)
%process ncs data for gui plots
global processed plotParam 
bandfilt=plotParam.filtlfp;
getburst=0;
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'getbursts'
            %get beta bursts based on threshold
            getburst=1;
    end
    argnum=argnum+1;
end

for ilfp=1:length(processed.LFPread.LFPchNames)
    sitelfp=[];
    samplesncs=processed.samplesNCS(:,ilfp)';
    if ismember(ilfp,plotParam.lfpid)
        sitelfp=plotParam.cscNames{ilfp};
        %{
        bandid=find(strcmp(sitelfp,plotParam.fbands)==1);
        if ~isempty(bandid)            
        %use first listed filter band in the list
        bandfilt=plotParam.fbands{bandid(1)+1};
        end
        %}
        %filter at "beta-band" as defined in fbands     
        processed.resampled(ilfp,:)=filterLFP(samplesncs,plotParam.ratelfp,bandfilt);
        %square & envelope signal
        winlength=round(plotParam.ratelfp*.5/mean(bandfilt));
        if isfield(plotParam,'winlength')
            %if window already defined in parameters load value
            winlength=round(plotParam.ratelfp*plotParam.winlength);
        end
        squaredata=processed.resampled(ilfp,:).^2;   %get power V^2
                       envdata=envwave(squaredata);     % 1/2019

        sqenvdata=smoothwin(envdata,winlength);
        processed.resampledb(ilfp,:)=sqenvdata;   %smoothing
        processed.sitelfp{ilfp}=sitelfp;
        processed.sitefilts{ilfp}=bandfilt;
        if getburst==1
            %find "bursts" of beta
            %thres=100;
            thres=std(sqenvdata)*3.5;
            idsbeta=find(sqenvdata>=thres);
            breakpts=find(abs(diff(idsbeta))>1);
            idsplot=[];
            if ~isempty(breakpts)
                for ib=1:length(breakpts)+1
                    if ib==1
                        idsplot{ib}=idsbeta(1:breakpts(ib));
                    elseif ib>1 && ib<=length(breakpts)
                        idsplot{ib}=idsbeta(breakpts(ib-1)+1:breakpts(ib));
                    elseif ib==length(breakpts)+1
                        idsplot{ib}=idsbeta(breakpts(ib-1)+1:end);
                    end
                   % plot(a1,ts(idsplot),plotdata(idsplot),'color',[1 0 0],'linestyle','-','marker','none');
                end
                processed.betaids{ilfp}=idsplot;
            else
               % betaids=sort(unique([idsplots idsbeta]));
                %plot(a1,ts(idsbeta),plotdata(idsbeta),'color',[1 0 0],'linestyle','-','marker','none');
                processed.betaids{ilfp}=idsbeta;
            end
        end
    elseif ismember(ilfp,plotParam.lickid)
        sitelfp=plotParam.cscNames{plotParam.lickid};
        %lick filter
        processed.resampled(ilfp,:)=filterLFP(samplesncs,plotParam.ratelfp,plotParam.filtlick);
        %square & envelope signal
        winlength=round(plotParam.ratelfp*.5/mean(plotParam.filtlick));
        if isfield(plotParam,'winlengthphys')
            %if window already defined in parameters load value
            winlength=round(plotParam.ratelfp*plotParam.winlengthphys);
        end
        processed.resampledb(ilfp,:)=smoothwin(processed.resampled(ilfp,:),winlength);   %smoothing
        processed.sitelfp{ilfp}=sitelfp;
        processed.sitefilts{ilfp}=plotParam.filtlick;
    %eye/pulse signal --> dont do anything already filt in reconvertFSCV.m            
    else
        sitelfp=plotParam.cscNames{ilfp};
        processed.resampled(ilfp,:)=samplesncs;
        processed.resampledb(ilfp,:)=samplesncs;
        processed.sitelfp{ilfp}=sitelfp;
        processed.sitefilts{ilfp}=[0 0];
    end
end

