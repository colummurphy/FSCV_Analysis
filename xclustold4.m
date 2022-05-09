function xinfo=xclust(xcovdata,plotparam,varargin)
%PRE 2/4/2019 after making major changes to analysis
%
%updated 12/29/2018 to look at timing relationships of raw & lags
%revised 1/11/2019 (old one saved as xclustold3)
%xclustold3 looked at lfp dynamics relative to da rise
%now look at relative to cue provided for all signals
%realized that when relative to da rise get increases in delt post for
%intertargeye, but decreases for intertarg
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
pbound= 15;     %max pos boundary to search for peaks relative to cue idx 1/13/2018
pade=-2;        %Changed 1/26/2019 SO THAT STRADDLING AROUND ALN EVENT NOT PRIOR
basesamples=3;
basepad=1;          %+/- .1 s for baseline 1/26/2019
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
eyetype=0;
noff=0;     %negative offset
if contains(eventtype,'eye')
   %need to offset some search win more negative if relative to eye
    eyetype=1;
    noff=-2;
end
if contains(eventtype,'fix')
    pbound=12;
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
%baseline immediately prior to aln idx w/ pad PREVIOUSLY 4 SAMPLES WIDE
%FIX
%baseids=meanalnwin-pade-basesamples:meanalnwin-pade;   
%baseids=meanalnwin-pade-basesamples:meanalnwin-pade-1;   
baseids=meanalnwin-basepad:meanalnwin+basepad;   

baseline=nanmean(daaln(:,baseids),2);
%baseline subtracted DA
daaln=daaln-baseline;
%get da trials where positive da increase
%previously determined in plotxvar not most accurate way
%IF DEFINE MAX DA WITHIN 1S BUT SEARCH WIN BELOWS ARE > 1S???????????
maxdas=max(daaln(:,meanalnwin:meanalnwin+pbound),[],2,'omitnan');
abovethres=find(maxdas>dabase*2);
mindas=min(daaln(:,meanalnwin:meanalnwin+pbound),[],2,'omitnan');
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
    'lfppostmax','lfppostmaxts','zlagcoef','maxprelagts',...
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
maxbound=size(dawin,2);
dadiff=diff(dawin,1,2);
lfpdiff=diff(lfpwin,1,2);
curmid=-wins(1)+1; %alignment idx for wins windowed signal, 
%also reference point for all ts's curr group, relative wins(1)
winoffset=meanalnwin+wins(1)-1;     %add to wins(1) relative ts's found to derive abs ts
if contains(eventtype,'fix')
   %need to bound to 1.2 s not size of win defined above
    maxbound=curmid+pbound;
end

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
    lfpmin=nan; %proceeds or precedes cue (in eye case)
    lfpmints=nan;       %proceeds or precedes cue (in eye case)
    lfprisets=nan;
    lfpfallts=nan;      %precedes cue
    lfppostmax=nan;     %after min
    lfppostmaxts=nan;
    
    if ~ismember(ii,decsonly)
        %if da pos, look for max rise around  following cue
        %or if any other type of signal besides decreases
        risets=[];
        swin=curmid-5+noff:curmid+10;
        [darisemax,risets]=max(dadiff(ii,swin),[],'omitnan');
        if ~isnan(darisemax) 
            darisets=risets(1)+swin(1)-1;
            if darisets<swin(end) && darisets>swin(1)
                %make sure not at edge
                %look for min prior to da rise
                minwin=swin(1):darisets;
                if length(minwin)>3
                    damin=min(dawin(ii,minwin),[],'omitnan');
                    %look for max fall prior to min da
                    if ~isnan(damin)
                        damints=find(dawin(ii,minwin)==damin)+minwin(1)-1;
                        if damints>swin(1)
                            %make sure not at edge of boundaries, false min
                            fwin=curmid-10+noff:damints;
                            if length(fwin)>3
                                dafall=min(dadiff(ii,fwin),[],'omitnan');
                                if ~isnan(dafall)                    
                                    dafallts=find(dadiff(ii,fwin)==dafall)+fwin(1)-1;
                                end
                            end
                        else
                            %at edge boundary, make undefined
                            damints=nan;
                            damin=nan;
                        end
                    end
                end
            %look for max da after da rise up to 1.5s after or max win
            maxwin=darisets:darisets+15;
            maxwin=maxwin(maxwin<=maxbound);
            damax=max(dawin(ii,maxwin),[],'omitnan');
            if ~isnan(damax)
                damaxts=find(dawin(ii,maxwin)==damax)+maxwin(1)-1;
                if damaxts>=maxbound
                    %if at boundary, may continue increasing, false max
                    damaxts=nan;
                    damax=nan;
                end
            end            
        else
             %if darise near bounds, all nan
            darisets=nan;
        end
        end
    end
    
    if ismember(ii,decsonly) 
        %if da neg, look for max fall around  following cue first
        fwin=curmid-10+noff:curmid+10;
        [dafall,fallts]=min(dadiff(ii,fwin),[],'omitnan');
        if ~isnan(dafall)
         dafallts=fallts(1)+fwin(1)-1;
        if dafallts>fwin(1) && dafallts<fwin(end)
            %make sure not at edges
            minwin=dafallts:fwin(end);
            %find min da after max da fall
            if length(minwin)>3
                damin=min(dawin(ii,minwin),[],'omitnan');
                if ~isnan(damin)
                    damints=find(dawin(ii,minwin)==damin)+minwin(1)-1;
                    if damints<fwin(end)
                        %if not at edges
                    %look for subsequent rise da max
                    swin=damints:damints+15;
                    swin=swin(swin<=maxbound);
                    if length(swin)>3
                        darise=max(dadiff(ii,swin),[],'omitnan');
                        if ~isnan(darise)
                            darisets=find(dadiff(ii,swin)==darise)+swin(1)-1;
                            if darisets>=maxbound-1
                                darisets=nan;
                            end
                        end
                    end
                    else
                        %if damints  at boundaires, all nan
                        damints=nan;
                        damin=nan;
                        dafallts=nan;
                    end
                end
            end
            %look for subsequent max da, only if prev var found ie non-nan
            if ~isnan(darisets) && ~isempty(darisets)
                maxwin=darisets:darisets+15;
                maxwin=maxwin(maxwin<=maxbound);
                if length(maxwin)>3
                    damax=max(dawin(ii,maxwin),[],'omitnan');
                    if ~isnan(damax)
                        damaxts=find(dawin(ii,maxwin)==damax)+...
                            maxwin(1)-1;
                    end
                    if damaxts>=maxbound
                        %at boundary, all nan
                        damaxts=nan;
                        damax=nan;
                    end
                end
            end
        else
            %if at edges, all nan
            dafallts=nan;
        end 
        end            
    end
    
    %get lfp dynamics, (now irregardless of da dynamic)
    %first search for min around cue
    minwin=curmid-1+noff:curmid+12;
    lfpmin=min(lfpwin(ii,minwin),[],'omitnan');
    fwin=[];
    if ~isnan(lfpmin)
        lfpmints=find(lfpwin(ii,minwin)==lfpmin)+minwin(1)-1;
        fwin=lfpmints-10:lfpmints;
        if lfpmints>minwin(1) && lfpmints<minwin(end) 
        %make sure not at edges
            if length(fwin)>3
                lfpfall=min(lfpdiff(ii,fwin),[],'omitnan');
                if ~isnan(lfpfall)
                    lfpfallts=find(lfpdiff(ii,fwin)==lfpfall)+fwin(1)-1;
                    if lfpfallts>fwin(1)
                        %if not at edges
                        %get max peak pre fall
                        maxwin=curmid-10:lfpfallts;
                        if length(maxwin)>3
                            lfpmax=max(lfpwin(ii,maxwin),[],'omitnan');
                            if ~isnan(lfpmax)
                                lfpmaxts=find(lfpwin(ii,maxwin)==lfpmax)+maxwin(1)-1;
                                if lfpmaxts<=maxwin(1)
                                    lfpmaxts=nan;
                                    lfpmax=nan;
                                end
                            end
                        end
                    else
                        lfpfallts=nan;
                    end
                end
            end
            %get max rise after lfp min
            %get max lfp first
            rwin=lfpmints:lfpmints+15;
            rwin=rwin(rwin<=maxbound);
            if length(rwin)>3
                lfppostmax=max(lfpwin(ii,rwin),[],'omitnan');
                if ~isnan(lfppostmax)
                    lfppostmaxts=find(lfpwin(ii,rwin)==lfppostmax)+rwin(1)-1;
                    if lfppostmaxts>=maxbound
                        %if at edge
                        lfppostmax=nan;
                        lfppostmaxts=nan;
                    end
                end                        
                lfprise=max(lfpdiff(ii,rwin),[],'omitnan');
                if ~isnan(lfprise)
                    lfprisets=find(lfpdiff(ii,rwin)==lfprise)+rwin(1)-1;
                    if lfprisets>=maxbound
                        %if at edge
                        lfprisets=nan;
                    end
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
    %lfp max/min/rise/fall, preceding cue (lfp) (da is after), rel damax/damin, darise/dafall
    lfpmax=getfield(xinfo(ilfp),dagroups{ii},'lfpmaxts');
    lfpmin=getfield(xinfo(ilfp),dagroups{ii},'lfpmints');
    lfprise=getfield(xinfo(ilfp),dagroups{ii},'lfprisets');
    lfpfall=getfield(xinfo(ilfp),dagroups{ii},'lfpfallts');
    lfppostmax=getfield(xinfo(ilfp),dagroups{ii},'lfppostmaxts');
    damax=getfield(xinfo(ilfp),dagroups{ii},'damaxts');
    delt_lfpmax_damax=damax-lfpmax;
    delt_lfpmin_damax=damax-lfpmin;
    delt_lfprise_damax=damax-lfprise;
    delt_lfpfall_damax=damax-lfpfall;
    delt_lfppostmax_damax=damax-lfppostmax;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmax_damax',delt_lfpmax_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpfall_damax',delt_lfpfall_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmin_damax',delt_lfpmin_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfprise_damax',delt_lfprise_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfppostmax_damax',delt_lfppostmax_damax);
    
    darise=getfield(xinfo(ilfp),dagroups{ii},'darisets');
    delt_lfpmax_darise=darise-lfpmax;
    delt_lfpmin_darise=darise-lfpmin;
    delt_lfprise_darise=darise-lfprise;
    delt_lfpfall_darise=darise-lfpfall;
        delt_lfppostmax_darise=darise-lfppostmax;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmax_darise',delt_lfpmax_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpfall_darise',delt_lfpfall_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmin_darise',delt_lfpmin_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfprise_darise',delt_lfprise_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfppostmax_darise',delt_lfppostmax_darise);
    
    damin=getfield(xinfo(ilfp),dagroups{ii},'damints');
    delt_lfpmax_damin=damin-lfpmax;
    delt_lfpmin_damin=damin-lfpmin;
    delt_lfprise_damin=damin-lfprise;
    delt_lfpfall_damin=damin-lfpfall;
            delt_lfppostmax_damin=damin-lfppostmax;
        xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmax_damin',delt_lfpmax_damin);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpfall_damin',delt_lfpfall_damin);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmin_damin',delt_lfpmin_damin);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfprise_damin',delt_lfprise_damin);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfppostmax_damin',delt_lfppostmax_damin);

%skip damin for now..
    dafall=getfield(xinfo(ilfp),dagroups{ii},'dafallts');
    delt_lfpmax_dafall=dafall-lfpmax;
    delt_lfpmin_dafall=dafall-lfpmin;
    delt_lfprise_dafall=dafall-lfprise;
    delt_lfpfall_dafall=dafall-lfpfall;
            delt_lfppostmax_dafall=dafall-lfppostmax;
        xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmax_dafall',delt_lfpmax_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpfall_dafall',delt_lfpfall_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmin_dafall',delt_lfpmin_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfprise_dafall',delt_lfprise_dafall);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfppostmax_dafall',delt_lfppostmax_dafall);
    
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
