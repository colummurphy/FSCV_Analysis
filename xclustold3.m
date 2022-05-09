function xinfo=xclust(xcovdata,plotparam,varargin)
%updated 12/29/2018 to look at timing relationships of raw & lags
xinfo=[];
fwidth=.5;      %fractional width of peak
thres=7;
sessiontype='';
if ~isempty(plotparam)
sessiontype=plotparam.trialtype;
end
bdata={};
alllag=0;
argnum=1;
eventtype=[];
wins=[-20 30];      %search win for peaks (in samples, 10hz rate assum)
pade=2;
basesamples=3;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'behav'
            %get behavioral data and group by lfp pairs/events/groups
            argnum=argnum+1;
            bdata=varargin{argnum};
        case 'event'
            %get event name
            argnum=argnum+1;
            eventtype=varargin{argnum};
    end
    argnum=argnum+1;
end

ratelfp=plotparam.ratelfp;

for ilfp=1:length(xcovdata)
    %scroll all non-empty lfp chs for given da ch
    if isempty(xcovdata{ilfp})
        continue
    end
rate=xcovdata{ilfp}.rate;
tsx=xcovdata{ilfp}.tsx;
xdata=xcovdata{ilfp}.xcovda;
seltrials=xcovdata{ilfp}.seltrials;
winids=xcovdata{ilfp}.winids;
xcovlag=xcovdata{ilfp}.xcovlag;
xcovlaganti=xcovdata{ilfp}.xcovlaganti;
%note da/lfps are indexed by original trial nums (seltrials),
%other data already only contain seltrials only
da=xcovdata{ilfp}.dadata;       
lfps=xcovdata{ilfp}.lfpdata;

%get trials only when lags not nan & da not nan
notnantrials=find(~isnan(xcovlaganti(:,1))  & ~isnan(xcovlag(:,1)));    %ids from seltrials
%get trials that don't have peaks above  cutoff exclusion threshold
gpos=find(xcovlag(:,2)<=nanmean(xcovlag(:,2))*thres);
gneg=find(xcovlaganti(:,2)>=nanmean(xcovlaganti(:,2))*thres);
gsel=intersect(notnantrials,unique([gpos; gneg]));
goodtrials=seltrials(gsel);     %trial ids for orig data

%just get good trials without nan lags
da=xcovdata{ilfp}.dadata(goodtrials,:);
dabase=xcovdata{1}.daavgstds;
lfp=xcovdata{ilfp}.lfpdata(goodtrials,:);
lfpbase=xcovdata{ilfp}.betabaseline;
lag=xcovlag(gsel,:);
lagneg=xcovlaganti(gsel,:);
winids=winids(gsel,:);
xcovda=xdata(gsel,:);
alignidx=median(1:size(winids,2));
meanalnwin=round(mean(winids(:,alignidx)));

%shifted da/lfp signals to mean align idx (middle of winids);
daaln=[];
lfpaln=[];
for ii=1:size(da,1)
    %numshift to align to meanalnwin
    shiftpos=round((meanalnwin-winids(ii,alignidx)));
    daaln(ii,:)=circshift(da(ii,:),shiftpos,2);
    lfpaln(ii,:)=circshift(lfp(ii,:),shiftpos,2);
end
%baseline immediately prior to aln idx w/ pad
baseids=meanalnwin-pade-basesamples:meanalnwin-pade;   
baseline=nanmean(daaln(:,baseids),2);
%baseline subtracted DA
daaln=daaln-baseline;
%get da trials where positive da increase
%previously determined in plotxvar not most accurate way
maxdas=max(daaln(:,meanalnwin:meanalnwin+10),[],2,'omitnan');
abovethres=find(maxdas>dabase*2);
mindas=min(daaln(:,meanalnwin:meanalnwin+10),[],2,'omitnan');
belowthres=find(mindas<-dabase*2);
rebounds=intersect(abovethres,belowthres);  %trials that show min/max are rebound usu
decsonly=setdiff(belowthres,rebounds);  %dec da trials without rebound only
dapos=daaln(abovethres,:);
lfppos=lfpaln(abovethres,:);
daneg=daaln(decsonly,:);
lfpneg=lfpaln(decsonly,:);
dareb=daaln(rebounds,:);
lfpreb=lfpaln(rebounds,:);

%targda=xcovdata{ilfp}.posttargda(gsel);
xinfo(ilfp).dapos.trials=goodtrials(abovethres);
xinfo(ilfp).dapos.datracesaln=dapos;
xinfo(ilfp).dapos.lfptracesaln=lfppos;
xinfo(ilfp).dapos.mididx=meanalnwin;
xinfo(ilfp).daneg.trials=goodtrials(decsonly);
xinfo(ilfp).daneg.datracesaln=daneg;
xinfo(ilfp).daneg.lfptracesaln=lfpneg;
xinfo(ilfp).daneg.mididx=meanalnwin;
xinfo(ilfp).dareb.trials=goodtrials(rebounds);
xinfo(ilfp).dareb.datracesaln=dareb;
xinfo(ilfp).dareb.lfptracesaln=lfpreb;
xinfo(ilfp).dareb.mididx=meanalnwin;
xinfo(ilfp).daall.trials=goodtrials;
xinfo(ilfp).daall.datracesaln=daaln;
xinfo(ilfp).daall.lfptracesaln=lfpaln;
xinfo(ilfp).daall.mididx=meanalnwin;
%get characteristics of da/lfp peaks for pos/neg groups
countpos=0;
countneg=0;
countreb=0;
params={'damax','damaxts','damin','damints','dafallts','darisets',...
    'lfpmax','lfpmaxts','lfpmin','lfpmints','lfprisets','lfpfallts',...
    'lfppostmax','lfppostmaxts','lfppostmin','lfppostmints',...
    'lfppostrisets','lfppostfallts','zlagcoef','maxprelagts',...
    'minprelagts','maxpostlagts','minpostlagts','maxprecoef',...
    'minprecoef','maxpostcoef','minpostcoef'};
    
%initialize parameters storage w/nans default
for ii=1:length(params)
    xinfo(ilfp).dapos=setfield(xinfo(ilfp).dapos,params{ii},nan(1,size(dapos,1)));
    xinfo(ilfp).daneg=setfield(xinfo(ilfp).daneg,params{ii},nan(1,size(daneg,1)));
    xinfo(ilfp).dareb=setfield(xinfo(ilfp).dareb,params{ii},nan(1,size(dareb,1)));
    xinfo(ilfp).daall=setfield(xinfo(ilfp).daall,params{ii},nan(1,size(daaln,1)));
end
%dawin=dapos(:,meanalnwin+wins(1):meanalnwin+wins(2));
%lfpwin=lfppos(:,meanalnwin+wins(1):meanalnwin+wins(2));
dawin=daaln(:,meanalnwin+wins(1):meanalnwin+wins(2));
lfpwin=lfpaln(:,meanalnwin+wins(1):meanalnwin+wins(2));
dadiff=diff(dawin,1,2);
lfpdiff=diff(lfpwin,1,2);
curmid=-wins(1)+1; %alignment idx for wins windowed signal, 
%also reference point for all ts's curr group, relative wins(1)
winoffset=meanalnwin+wins(1)-1;     %add to wins(1) relative ts's found to derive abs ts

for ii=1:size(dawin,1)
    %da dynamics around provided cue idx
    damax=nan;
    damaxts=nan;
    damin=nan;
    damints=nan;
    dafallts=nan;
    darisets=nan;
    %beta-lfp dynamics preceding da rise bounded by cue idx
    lfpmax=nan;     %default if no da rise found
    lfpmaxts=nan;
    lfpmin=nan;
    lfpmints=nan;
    lfprisets=nan;
    lfpfallts=nan;
    %lfp dynamics proceeding after da rise
    lfppostmax=nan;     %default if no da rise found
    lfppostmaxts=nan;
    lfppostmin=nan;
    lfppostmints=nan;
    lfppostrisets=nan;
    lfppostfallts=nan;
    
    %da dynamics
    %max after within 1.5 s of aln idx
        %for intertargeye, could occur sligthly before eye on cue if responsive
    %NENED TO UPDATE
    %MAYBE SHOULD NOT CALCULATE ANY RELATIVE TO DSRISE OTHERWISE PRODUCES
    %FALSE DELT"S DEPENDING ON EVENT MARKER ALIGNMENT
    %RATHER MAKE ALL DIFERNCES RELATIVE TO CUE SO UNIFORM REFERNCING IN
    %TIME
    damax=max(dawin(ii,curmid:curmid+15),[],'omitnan');
    %damaxts=find(dawin(ii,curmid:end)==damax)+meanalnwin-1;
    damaxts=find(dawin(ii,curmid:curmid+15)==damax)+curmid-1;      %relative to wins(1)
    %min within 1s 
    damin=min(dawin(ii,curmid:curmid+10),[],'omitnan');
    damints=find(dawin(ii,curmid:curmid+10)==damin)+curmid-1;
    %slope peaks, extend win since rise/falls occur earlier than peaks
    dafallmaxts=find(dadiff(ii,damints-10:damints+10)==...
        min(dadiff(ii,damints-10:damints+10)))+damints-10-1;
    if length(dafallmaxts)>1
        dafallmaxts=dafallmaxts(1);
    end
    %if min at either boundary, nan undefined
    if dafallmaxts<=curmid-10 || dafallmaxts>=curmid+10
        dafallmaxts=nan;
    else
        %find peak in diff prior to fall as window boundary
        daprepeakdiffts=find(dadiff(ii,damints-10:dafallmaxts)>=0)...
            +damints-10-1;
        if ~isempty(daprepeakdiffts)
            daprepeakdiffts=daprepeakdiffts(end);
            %fall time actually earlier than diff peak, get this initial point
            %based on window of prior diff pos peak to fall peak
            if isequal(daprepeakdiffts,dafallmaxts)
                dafallts=daprepeakdiffts;
            else
                dafallts=find(dadiff(ii,daprepeakdiffts:dafallmaxts)<0)+daprepeakdiffts-1;
                if ~isempty(dafallts)
                    dafallts=dafallts(1);
                end
            end
            if dafallts<=curmid-10 || dafallts>=curmid+10
                    %if min at either boundary, nan undefined
                dafallts=nan;
            end
        end
    end
    
    %find peak rise slope of da, end point is previously found max
    %for intertargeye, could occur sligthly before eye on cue if responsive
    %NENED TO UPDATE
    darisemaxts=find(dadiff(ii,curmid:damaxts)==...
        max(dadiff(ii,curmid:damaxts),[],'omitnan'))+curmid-1;
    if length(darisemaxts)>1
        darisemaxts=darisemaxts(1);
    end
    %if max at either boundary, nan undefined
    if darisemaxts<=curmid || darisemaxts>=curmid+15
        darisemaxts=nan;
    else
    %find initial da increase point by first finding prior min in diff da
    dapremindiffts=find(dadiff(ii,curmid-10:darisemaxts)<=0)...
        +curmid-10-1;
    if ~isempty(dapremindiffts)
        dapremindiffts=dapremindiffts(end);
        %rise time initial point based on window of prior diff min peak 
        if isequal(dapremindiffts,darisemaxts)
            darisets=dapremindiffts;
        else
            darisets=find(dadiff(ii,dapremindiffts:darisemaxts)>0)+dapremindiffts-1;
            darisets=darisets(1);
        end
        if darisets<=curmid-10 || darisets>=curmid+14
                %if max at either boundary, nan undefined
            darisets=nan;
        end
    end 
    end
    
    %beta-lfp dynamics preceding da rise
    if ~isnan(darisets)
        lfpwinend=darisets;     %set da local peak as ts end to search beta burst preceding 
        %get max beta 'preceding' da max diff within 0.75s of cue
        lfpmax=max(lfpwin(ii,curmid-10:lfpwinend),[],'omitnan');     %default maxpreceding
        lfpmaxts=find(lfpwin(ii,curmid-10:lfpwinend)==lfpmax)+curmid-10-1;
        %get min beta suppresion 
        lfpmin=min(lfpwin(ii,curmid-10:lfpwinend),[],'omitnan');     %default maxpreceding
        lfpmints=find(lfpwin(ii,curmid-10:lfpwinend)==lfpmin)+curmid-10-1; 
        %lfp min could be preceding beta peak (c's) or after (p's)
        if lfpmaxts==curmid-10
            %at starting edge, see if another peak exists after min
            lfpmaxtemp=max(lfpwin(ii,lfpmints:lfpwinend),[],'omitnan');     %default maxpreceding
            lfpmaxtstemp=find(lfpwin(ii,lfpmints:lfpwinend)==lfpmaxtemp)+lfpmints-1;
            lfpmax=lfpmaxtemp;
            lfpmaxts=lfpmaxtstemp;
        end
        %get max beta rise bounded by lfpmin preceding lfpmaxts
        if lfpmints<lfpmaxts
            %min precedes peak
            startwin=lfpmints-10;
            if startwin<curmid-17
                startwin=curmid-17;
            end
        else
            startwin=curmid-15;
        end
        lfprisemaxts=find(lfpdiff(ii,startwin:lfpmaxts)==...
            max(lfpdiff(ii,startwin:lfpmaxts),[],'omitnan'))+startwin-1;
        lfprisemaxts=lfprisemaxts(1);
        %find initial lfp increase point by first finding prior min in diff
        lfppremindiffts=find(lfpdiff(ii,curmid-17:lfprisemaxts)<=0)...
            +curmid-17-1;
        if ~isempty(lfppremindiffts)
            lfppremindiffts=lfppremindiffts(end);
            %rise time initial point based on window of prior diff min peak 
            if isequal(lfppremindiffts,lfprisemaxts)
                lfprisets=lfppremindiffts;
            else
                lfprisets=find(lfpdiff(ii,lfppremindiffts:lfprisemaxts)>0)+lfppremindiffts-1;
                lfprisets=lfprisets(1);
            end
            if lfprisets<=curmid-17
                %on boundary
                lfprisets=nan;
            end
        end
        
        %get max beta fall after lfpmaxts (after peak, peak fall)
        startwin=lfpmaxts-5;
        if startwin>=curmid-15 || startwin<=curmid+5
            %only if within boundaries, can be defined
            lfpfallmaxts=find(lfpdiff(ii,startwin:startwin+10)==...
                min(lfpdiff(ii,startwin:startwin+10),[],'omitnan'))+startwin-1;
            lfpfallmaxts=lfpfallmaxts(1);
            %find initial lfp increase point by first finding prior min in diff
            lfppremaxdiffts=find(lfpdiff(ii,curmid-15:lfpfallmaxts)>=0)...
                +curmid-15-1;
            if ~isempty(lfppremaxdiffts)
                lfppremaxdiffts=lfppremaxdiffts(end);
                if isequal(lfppremaxdiffts,lfpfallmaxts)
                    lfpfallts=lfppremaxdiffts;
                else
                    %fall time initial point based on window of prior diff min peak 
                    lfpfallts=find(lfpdiff(ii,lfppremaxdiffts:lfpfallmaxts)<0)+lfppremaxdiffts-1;
                    lfpfallts=lfpfallts(1);
                end
                if lfpfallts<=curmid-15 || lfpfallts>=curmid+5
                    %if on boundaries, undefined
                    lfpfallts=nan;
                end
            end 
        end
    end
    %end lfp dynamics preceding da rise
    
    %lfp dynamics proceeding after da rise
    if ~isnan(darisets)
        lfpwinstart=darisets;     %set da local peak as ts start to search beta bu
        lfpwinend=curmid+17;
        curwin=lfpwinstart:lfpwinend;
        curwin=curwin(curwin<=size(lfpwin,2) & curwin>=1);
        %get max beta after da max diff within 1.5s after cue
        lfppostmax=max(lfpwin(ii,curwin),[],'omitnan');     
        lfppostmaxts=find(lfpwin(ii,curwin)==lfppostmax)+curwin(1)-1;
        %get min beta suppresion preceding burst but after cue
        curwin=curwin(curwin>=curmid);
        lfppostmin=min(lfpwin(ii,curwin),[],'omitnan');     %default maxpreceding
        lfppostmints=find(lfpwin(ii,curwin)==lfppostmin)+curwin(1)-1;        
        %get max beta rise before lfppostmax
        lfprisemaxts=find(lfpdiff(ii,curmid:lfppostmaxts)==...
            max(lfpdiff(ii,curmid:lfppostmaxts),[],'omitnan'))+curmid-1;
        if ~isempty(lfprisemaxts)
            lfprisemaxts=lfprisemaxts(1);
            %find initial lfp increase point by first finding prior min in diff
            lfppremindiffts=find(lfpdiff(ii,curmid-5:lfprisemaxts)<=0)...
            +curmid-5-1;
            if ~isempty(lfppremindiffts)
                lfppremindiffts=lfppremindiffts(end);
                %rise time initial point based on window of prior diff min peak 
                if isequal(lfppremindiffts,lfprisemaxts)
                    lfppostrisets=lfppremindiffts;
                else
                    lfppostrisets=find(lfpdiff(ii,lfppremindiffts:lfprisemaxts)>0)+lfppremindiffts-1;
                    lfppostrisets=lfppostrisets(1);
                end
                if lfppostrisets>=curmid+17 || lfppostrisets<=curmid
                    %on boundaries, undefined
                    lfppostrisets=nan;
                end
            end
        end
        
        %get max beta fall before beta min (could be before/after lfp max)
        lfpwinstart=darisets;     %set da local peak as ts start
        %COULD BE BEFORE< EXTEND WINDOW
        lfpwinstart=darisets-5;
        lfpwinend=lfppostmints+5;     %beta min as end
        curwin=lfpwinstart:lfpwinend;
        %get max beta fall slope before beta min
        lfpfallmaxts=find(lfpdiff(ii,curwin)==...
            min(lfpdiff(ii,curwin),[],'omitnan'))+curwin(1)-1;
        if ~isempty(lfpfallmaxts)
            lfpfallmaxts=lfpfallmaxts(1);
            %find initial lfp decrease point by first finding prior pos diff
            lfppremaxdiffts=find(lfpdiff(ii,curwin(1):lfpfallmaxts)>=0)...
                +curwin(1)-1;
            if ~isempty(lfppremaxdiffts)
                lfppremaxdiffts=lfppremaxdiffts(end);
                if isequal(lfppremaxdiffts,lfpfallmaxts)
                    lfppostfallts=lfppremaxdiffts;
                else
                    %fall time initial point based on window of prior diff min peak 
                    lfppostfallts=find(lfpdiff(ii,lfppremaxdiffts:lfpfallmaxts)<0)+lfppremaxdiffts-1;
                    lfppostfallts=lfppostfallts(1);
                end
                if lfppostfallts<=curmid-5 || lfppostfallts>=curmid+15
                    %on edges, make undefined
                    lfppostfallts=nan;
                end
            end 
        end
    end
    
    %get xcov parameters
    %get mean xcoeff near zero, (ie mean xcoeff at 0+/-100 ms
    midxcov=median(1:size(xcovda,2));
    zlagcoef=nanmean(xcovda(ii,midxcov-1:midxcov+1));
    %get leading lag peaks (ie lags < 0 s) both cor & anti
    maxprecoef=max(xcovda(ii,1:midxcov-1),[],'omitnan');
    maxprelagts=nan;
    if ~isempty(maxprecoef)
        maxprelagts=find(xcovda(ii,1:midxcov-1)==maxprecoef);
        maxprelagts=tsx(maxprelagts);
    end
    %anti-correlated peak leading
    minprecoef=min(xcovda(ii,1:midxcov-1),[],'omitnan');
    minprelagts=nan;
    if ~isempty(minprecoef)
        minprelagts=find(xcovda(ii,1:midxcov-1)==minprecoef);
        minprelagts=tsx(minprelagts);
    end
    %get lagging lag peaks (ie lags < 0 s) both cor & anti
    maxpostcoef=max(xcovda(ii,midxcov+1:end),[],'omitnan');
    maxpostlagts=nan;
    if ~isempty(maxpostcoef)
        maxpostlagts=find(xcovda(ii,midxcov+1:end)==maxpostcoef)+midxcov;
        maxpostlagts=tsx(maxpostlagts);
    end
    %anti-correlated peak leading
    minpostcoef=min(xcovda(ii,midxcov+1:end),[],'omitnan');
    minpostlagts=nan;
    if ~isempty(minpostcoef)
        minpostlagts=find(xcovda(ii,midxcov+1:end)==minpostcoef)+midxcov;
        minpostlagts=tsx(minpostlagts);
    end
    
    %store in dapos or daneg cell according to da value for trial
    group='';
    curridx=nan;
    %exclude nans & 0s
    if ismember(ii,abovethres)
        group='dapos';
        countpos=countpos+1;
        curridx=countpos;
    elseif ismember(ii,decsonly)
        group='daneg';
        countneg=countneg+1;
        curridx=countneg;
    elseif ismember(ii,rebounds)
        group='dareb';
        countreb=countreb+1;
        curridx=countreb;
    end
    groups{1}=group;
    groups{2}='daall';      %include all signals as well
    for igroup=1:length(groups)
      if ~isempty(groups{igroup})
        if strcmp(groups{igroup},'daall')
            curridx=ii;
        end
        for iparam=1:length(params)
            storevalue=eval(params{iparam});
            if contains(params{iparam},'ts') && ~contains(params{iparam},'lag')
                storevalue=eval(params{iparam})+winoffset;
            end
            xinfo(ilfp)=setfield(xinfo(ilfp),groups{igroup},params{iparam},{curridx},storevalue);
        end
      end
    end
end

%get ts differences between da peaks and beta bursts/suppresions
dagroups={'dapos','daneg','dareb','daall'};
for ii=1:length(dagroups)
    %lfp max/min/rise/fall, preceding da max, rel damax/damin, darise/dafall
    lfpmax=getfield(xinfo(ilfp),dagroups{ii},'lfpmaxts');
    lfpmin=getfield(xinfo(ilfp),dagroups{ii},'lfpmints');
    lfprise=getfield(xinfo(ilfp),dagroups{ii},'lfprisets');
    lfpfall=getfield(xinfo(ilfp),dagroups{ii},'lfpfallts');
    damax=getfield(xinfo(ilfp),dagroups{ii},'damaxts');
    delt_preda_lfpmax_damax=damax-lfpmax;
    delt_preda_lfpmin_damax=damax-lfpmin;
    delt_preda_lfprise_damax=damax-lfprise;
    delt_preda_lfpfall_damax=damax-lfpfall;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpmax_damax',delt_preda_lfpmax_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpmin_damax',delt_preda_lfpmin_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfprise_damax',delt_preda_lfprise_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpfall_damax',delt_preda_lfpfall_damax);
    
    darise=getfield(xinfo(ilfp),dagroups{ii},'darisets');
    delt_preda_lfpmax_darise=darise-lfpmax;
    delt_preda_lfpmin_darise=darise-lfpmin;
    delt_preda_lfprise_darise=darise-lfprise;
    delt_preda_lfpfall_darise=darise-lfpfall;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpmax_darise',delt_preda_lfpmax_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpmin_darise',delt_preda_lfpmin_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfprise_darise',delt_preda_lfprise_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpfall_darise',delt_preda_lfpfall_darise);
%skip damin for now..
    dafall=getfield(xinfo(ilfp),dagroups{ii},'dafallts');
    delt_preda_lfpmax_dafall=dafall-lfpmax;
    delt_preda_lfpmin_dafall=dafall-lfpmin;
    delt_preda_lfprise_dafall=dafall-lfprise;
    delt_preda_lfpfall_dafall=dafall-lfpfall;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpmax_dafall',delt_preda_lfpmax_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpmin_dafall',delt_preda_lfpmin_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfprise_dafall',delt_preda_lfprise_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_preda_lfpfall_dafall',delt_preda_lfpfall_dafall);

    %lfp max/min/rise/fall, after da rise, rel damax/damin, darise/dafall
    lfppostmaxts=getfield(xinfo(ilfp),dagroups{ii},'lfppostmaxts');
    lfppostmints=getfield(xinfo(ilfp),dagroups{ii},'lfppostmints');
    lfppostrisets=getfield(xinfo(ilfp),dagroups{ii},'lfppostrisets');
    lfppostfallts=getfield(xinfo(ilfp),dagroups{ii},'lfppostfallts');
    delt_postda_lfpmax_damax=lfppostmaxts-damax;
    delt_postda_lfpmin_damax=lfppostmints-damax;
    delt_postda_lfprise_damax=lfppostrisets-damax;
    delt_postda_lfpfall_damax=lfppostfallts-damax;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpmax_damax',delt_postda_lfpmax_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpmin_damax',delt_postda_lfpmin_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfprise_damax',delt_postda_lfprise_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpfall_damax',delt_postda_lfpfall_damax);
    
    delt_postda_lfpmax_darise=lfppostmaxts-darise;
    delt_postda_lfpmin_darise=lfppostmints-darise;
    delt_postda_lfprise_darise=lfppostrisets-darise;
    delt_postda_lfpfall_darise=lfppostfallts-darise;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpmax_darise',delt_postda_lfpmax_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpmin_darise',delt_postda_lfpmin_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfprise_darise',delt_postda_lfprise_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpfall_darise',delt_postda_lfpfall_darise);
%skip damin for now..
    delt_postda_lfpmax_dafall=lfppostmaxts-dafall;
    delt_postda_lfpmin_dafall=lfppostmints-dafall;
    delt_postda_lfprise_dafall=lfppostrisets-dafall;
    delt_postda_lfpfall_dafall=lfppostfallts-dafall;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpmax_dafall',delt_postda_lfpmax_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpmin_dafall',delt_postda_lfpmin_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfprise_dafall',delt_postda_lfprise_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_postda_lfpfall_dafall',delt_postda_lfpfall_dafall);
    
end    
  
%save variables
xinfo(ilfp).sessiontype=sessiontype;
xinfo(ilfp).sessionperiod=plotparam.triallabel;
xinfo(ilfp).event=xcovdata{ilfp}.eventnames;
if isfield(xcovdata{ilfp},'sitenameda')
    xinfo(ilfp).siteda=xcovdata{ilfp}.sitenameda;
    xinfo(ilfp).sitelfp=xcovdata{ilfp}.sitename;
else
    %new format
    xinfo(ilfp).siteda=xcovdata{ilfp}.sitename{2};
    xinfo(ilfp).sitelfp=xcovdata{ilfp}.sitename{1};
end
xinfo(ilfp).goodtrials=goodtrials;        %based on original trial id's
xinfo(ilfp).freq=xcovdata{ilfp}.freqband;
xinfo(ilfp).tsx=tsx;
xinfo(ilfp).xcovda=xcovda;


end

end
