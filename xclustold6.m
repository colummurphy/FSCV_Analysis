function xinfo=xclust(xcovdata,plotparam,varargin)
%2/18/2019, SESSNUM NOT DEFINED CORRECTLY JUST TAKES FIRST 2 DIGITS FROM
%PATH NAME ! SO PBOUND INCORRECT SESS 113 & 114
%2/12/2019, revised parameters for timing, now working, show sig
%delt_lfpmin_damax (& delt_lfppostmax_damax) 73,83, 92 patra
%2/11/2019, revise again, make windows longer for searching beta rebounds &
%da max's, save old one as xclustold5
%2/4/2019 after making major changes to analysis
%NEED TO MAKE WINDOWS WIDER BASED ON PROLONGED DYNAMICS SSEEN IN TARG
%now all input xcovdata is prealigned
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
wins=[-20 50];      %search win for peaks (in samples, 10hz rate assum)
pbound= 15;     %max pos boundary to search for peaks relative to cue idx 1/13/2018
pade=-2;        %Changed 1/26/2019 SO THAT STRADDLING AROUND ALN EVENT NOT PRIOR
basesamples=3;
basepad=1;          %+/- .1 s for baseline 1/26/2019
sessfind=strfind(plotparam.pathname,'chronic');
sessnum=str2num(plotparam.pathname(sessfind+7:sessfind+8)); %sess num indicates length pbound for target
if sessnum<20
    sessnum2=str2num(plotparam.pathname(sessfind+7:sessfind+9));
    if isempty(sessnum2)
        %indeed sessnum is < 20
    else
        sessnum=sessnum2;
    end
end
isbeh=0;    
nolfp=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'beh'
            %get behavioral data and group by lfp pairs/events/groups
            isbeh=1;
        case 'event'
            %get event name
            argnum=argnum+1;
            eventtype=varargin{argnum};
        case 'nolfp'
            nolfp=1;
    end
    argnum=argnum+1;
end
eyetype=0;
noff=0;     %negative offset
if contains(eventtype,'eye')
   %need to offset some search win more negative if relative to eye
    eyetype=1;
    noff=-4;        %change from 2 to 4 2/6/2019
end
if contains(eventtype,'fix')
    pbound=12;
end
if contains(eventtype,'targ')
    pbound=35;      %changed from 35 to 39 2/11/2019,back to 35 2/12/
    if sessnum<67
        pbound=18; %NEED TO EXTEND OUT FOR SMALLER TARG SINCE PEAK OCCURS LATER IN BR both for da & beta reobund
    end
end

ratelfp=plotparam.ratelfp;
if nolfp
    xcovtemp=xcovdata;
    xcovdata={};
    xcovdata{1}=xcovtemp;
end

for ilfp=1:length(xcovdata)
    %scroll all non-empty lfp chs for given da ch
    if isempty(xcovdata{ilfp})
        continue
    end
rate=xcovdata{ilfp}.rate;
tsx=xcovdata{ilfp}.tsx;
seltrials=xcovdata{ilfp}.seltrials;
winids=xcovdata{ilfp}.winids;
%note da/lfps are indexed by original trial nums (seltrials),
%other data already only contain seltrials only
da=xcovdata{ilfp}.dadata;       
goodtrials=seltrials;
gsel=1:length(seltrials);
xdata=[];
xcovlag=[];
xcovlaganti=[];
lfp=[];
xcovda=[];
if ~nolfp
xdata=xcovdata{ilfp}.xcovda;
xcovlag=xcovdata{ilfp}.xcovlag;
xcovlaganti=xcovdata{ilfp}.xcovlaganti;
lfps=xcovdata{ilfp}.lfpdata;
if ~isbeh
    %get trials only when lags not nan & da not nan
    notnantrials=find(~isnan(xcovlaganti(:,1))  & ~isnan(xcovlag(:,1)));    %ids from seltrials
    %get trials that don't have peaks above  cutoff exclusion threshold
    gpos=find(xcovlag(:,2)<=nanmean(xcovlag(:,2))*thres);
    gneg=find(xcovlaganti(:,2)>=nanmean(xcovlaganti(:,2))*thres);
    gsel=intersect(notnantrials,unique([gpos; gneg]));
    goodtrials=seltrials(gsel);     %trial ids for orig data
end
lfp=xcovdata{ilfp}.lfpdata(goodtrials,:);
xcovda=xdata(gsel,:);

end

%just get good trials without nan lags
da=xcovdata{ilfp}.dadata(goodtrials,:);
dabase=xcovdata{1}.daavgstds;
%lfpbase=xcovdata{ilfp}.betabaseline;
%lag=xcovlag(gsel,:);
%lagneg=xcovlaganti(gsel,:);
alignidx=median(1:size(winids,2));  %middle of win set as align ts during xcov -----NO NOT ANYMORE 04/28/2019
alignidx=xcovdata{ilfp}.alnidx;     %ADDED 04/28/2019
%meanalnwin=round(mean(winids(:,alignidx)));
meanalnwin=alignidx;        %added 04/28/2019
%baseline immediately prior to aln idx w/ pad PREVIOUSLY 4 SAMPLES WIDE
%FIX
%baseids=meanalnwin-pade-basesamples:meanalnwin-pade;   
%baseids=meanalnwin-pade-basesamples:meanalnwin-pade-1;   
baseids=meanalnwin-basepad:meanalnwin+basepad;   

daaln=da;
lfpaln=lfp;
baseline=nanmean(daaln(:,baseids),2);
%baseline subtracted DA
daaln=daaln-baseline;
%get da trials where positive da increase
%previously determined in plotxvar not most accurate way
%IF DEFINE MAX DA WITHIN 1S BUT SEARCH WIN BELOWS ARE > 1S???????????
maxdas=max(daaln(:,meanalnwin:meanalnwin+pbound),[],2,'omitnan');
abovethresmax=find(maxdas>dabase*2);
avgdas=nanmean(daaln(:,meanalnwin+basepad:meanalnwin+pbound),2);
avgabovethres=find(avgdas>dabase*.5);      %now also look at atverage win to get pos signals
abovethres=intersect(abovethresmax,avgabovethres); %after storing xinfo_noff4 chronic 92, 2/6/2019, saved new xinfo_noff4_ab
abovethres=abovethresmax;           %revert 2/6/2019 5:20 PM
mindas=min(daaln(:,meanalnwin:meanalnwin+pbound),[],2,'omitnan');
belowthresmin=find(mindas<-dabase*2);
rebounds=intersect(abovethresmax,belowthresmin);  %trials that show min/max are rebound usu
decsonly=setdiff(belowthresmin,rebounds);  %dec da trials without rebound only
dapos=daaln(abovethres,:);
daneg=daaln(decsonly,:);
dareb=daaln(rebounds,:);
if ~nolfp
lfppos=lfpaln(abovethres,:);
lfpneg=lfpaln(decsonly,:);
lfpreb=lfpaln(rebounds,:);
end

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
maxbound=curmid+pbound;

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
        swin=curmid-5+noff:curmid+pbound;
        %TAKTE OUT RISE CONTINGENCY 2/18/2019, max cna occur any period
        [damax, maxts]=max(dawin(ii,swin),[],'omitnan');
        if ~isnan(damax) 
           damaxts=maxts+swin(1)-1;
        end 
        rwin=swin;
        if ~isnan(damaxts)
            rwin=swin(swin<=damaxts);
        end
        [darisemax,risets]=max(dadiff(ii,rwin),[],'omitnan');
        if ~isnan(darisemax) 
            % && risets>1 take off edge constrain 2/12/2019
            darisets=risets(1)+rwin(1)-1;
            %make sure not at edge
            %look for min prior to da rise
            minwin=swin(1):darisets;
            if length(minwin)>3
                [damin,mints]=min(dawin(ii,minwin),[],'omitnan');
                %look for max fall prior to min da
                if ~isnan(damin) && mints>1 
                        %not nan and not at edge
                    damints=mints+minwin(1)-1;
                    %make sure not at edge of boundaries, false min
                    fwin=curmid-10+noff:damints;
                    if length(fwin)>3
                        [dafall, fallts]=min(dadiff(ii,fwin),[],'omitnan');
                        if ~isnan(dafall)                    
                            dafallts=fallts(1)+fwin(1)-1;
                        end
                    end                        
                else
                    %at edge boundary, make undefined
                    damints=nan;
                    damin=nan;
               end
            end
        end
        %{
        TAKE OUT 2/18/2019
            %look for max da after da rise up to 1.5s after or max win
            maxwin=darisets:darisets+25;        %change from 15 to 25 2/12/2019
            maxwin=maxwin(maxwin<=maxbound);
            [damax, maxts]=max(dawin(ii,maxwin),[],'omitnan');
            if ~isnan(damax) 
                %&& maxts<length(maxwin) take out edge constratin 2/12/2019
                %not nan and not at edge
                damaxts=maxts+maxwin(1)-1;
            end        
        %}   
    end
    
    if ismember(ii,decsonly) 
        %if da neg, look for max fall around  following cue first
        fwin=curmid-10+noff:curmid+pbound;
        [dafall,fallts]=min(dadiff(ii,fwin),[],'omitnan');
        if ~isnan(dafall) && fallts<length(fwin) && fallts>1
         dafallts=fallts(1)+fwin(1)-1;
            %make sure not at edges
            minwin=dafallts:fwin(end);
            %find min da after max da fall
            if length(minwin)>3
                [damin,mints]=min(dawin(ii,minwin),[],'omitnan');
                if ~isnan(damin) && mints>1
                    damints=mints+minwin(1)-1;
                    %look for subsequent rise da max
                    swin=damints:damints+15;
                    swin=swin(swin<=maxbound);
                    if length(swin)>3
                        [darise,risets]=max(dadiff(ii,swin),[],'omitnan');
                        if ~isnan(darise) && risets<length(swin)
                            darisets=risets+swin(1)-1;
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
                    [damax,maxts]=max(dawin(ii,maxwin),[],'omitnan');
                    if ~isnan(damax) && maxts<length(maxwin)
                        damaxts=maxts+maxwin(1)-1;
                    end
                end
            end
        end 
    
    %get lfp dynamics, (now irregardless of da dynamic)
    %first search for min around cue
    %minwin=curmid-1+noff:curmid+13;     %RIGHT NEXT TO CUE, EXPECT TRANSIENT, changed from 12 to 10 2/11/2019, changed to 13 2/12/2019 (pl1 longer time to suppress)
    %change to 20 maybe 2/17/2019 >>>>>>>>>>>>>>>>>>>>>>>>>
    minwin=curmid-1+noff:curmid+20;
    %2/6/2019, noff needs to be wider for targeye, but then might increase
    %to larger values for targ eye..
    [lfpmin,mints]=min(lfpwin(ii,minwin),[],'omitnan');
    fwin=[];
    if ~isbeh
    if ~isnan(lfpmin) 
        %&& mints>1 && mints<length(minwin) take off edge constrain
        %2/12/2019
                %make sure not at edges
        lfpmints=mints+minwin(1)-1;
        fwin=lfpmints-10:lfpmints;
        if length(fwin)>3
            [lfpfall, fallts]=min(lfpdiff(ii,fwin),[],'omitnan');
            if ~isnan(lfpfall) && fallts>1
                lfpfallts=fallts+fwin(1)-1;
                %if not at edges get max peak pre fall
                maxwin=curmid-10:lfpfallts;
                if length(maxwin)>3
                    [lfpmax, maxts]=max(lfpwin(ii,maxwin),[],'omitnan');
                    if ~isnan(lfpmax) && maxts<length(maxwin)
                        lfpmaxts=maxts+maxwin(1)-1;
                    end
                end
            end
        end
        %get max rise after lfp min
        %get max lfp first
        %rwin=lfpmints:lfpmints+20;      %changed from 15 to 20, 02/6/2019
        %rwin=lfpmints:lfpmints+15;      %revert 2/6/2019 5:20 pm
        %rwin=lfpmints:lfpmints+pbound;      %change 2/11/2019
        %rwin=lfpmints:lfpmints+17;      %change 2/12/2019 given avg characteristics of min and max (want to lok for max rebound related to targ not max before rew)
        rwin=lfpmints:lfpmints+20;      %change again 2/12/2019 longer ts for rebound 73, also remove max lim ts constratin that has to be less than window, can be on edge
        %MAYGE SHOULD CHANGE TO +18 or PBOUND GIVEN DELAYS REB W SMALL TRS
        rwin=rwin(rwin<=maxbound);
        if length(rwin)>3
            [lfppostmax,maxts]=max(lfpwin(ii,rwin),[],'omitnan');
            if ~isnan(lfppostmax) 
                lfppostmaxts=maxts+rwin(1)-1;
                risewin=rwin(rwin<=lfppostmaxts);
                [lfprise,risets]=max(lfpdiff(ii,risewin),[],'omitnan');
                if ~isnan(lfprise) && risets<length(risewin)
                    lfprisets=risets+risewin(1)-1;  %LOOK FOR INITIAL RISE PAST MIN NOT MAX RISE
                end
            end
        end
    end      
    else
        %beh signal
        minwin=curmid-1+noff:curmid+pbound;
       [lfpmin,mints]=min(lfpwin(ii,minwin),[],'omitnan');
        if ~isnan(lfpmin) 
            lfpmints=mints+minwin(1)-1;
            fwin=lfpmints-10:lfpmints;
            if length(fwin)>3
                [lfpfall, fallts]=min(lfpdiff(ii,fwin),[],'omitnan');
                if ~isnan(lfpfall) && fallts>1
                    lfpfallts=fallts+fwin(1)-1;
                end
            end
            %max & post-max same for beh
            rwin=minwin;      %revert 2/6/2019 5:20 pm
            rwin=rwin(rwin<=maxbound);
            if length(rwin)>3
                [lfppostmax,maxts]=max(lfpwin(ii,rwin),[],'omitnan');
                lfpmax=lfppostmax;
                if ~isnan(lfppostmax) 
                    lfppostmaxts=maxts+rwin(1)-1;
                    lfpmaxts=lfppostmaxts;
                end
                mwin=minwin(1):lfppostmaxts;
                [lfprise,risets]=max(lfpdiff(ii,mwin),[],'omitnan');
                if ~isnan(lfprise) && risets<length(mwin)
                    lfprisets=risets+mwin(1)-1;
                end
            end
        end
    end
    
    %get xcov parameters
    %get mean xcoeff near zero, (ie mean xcoeff at 0+/-100 ms
        maxprelagts=nan;
    zlagcoef=nan;
    minprelagts=nan;
    maxpostlagts=nan;
    minpostlagts=nan;
    maxprecoef=nan;
    minprecoef=nan;
    maxpostcoef=nan;
    minpostcoef=nan;
    midxcov=median(1:size(xcovda,2));
    if nanmean(xcovda(ii,:))~=0 && ~isnan(nanmean(xcovda(ii,:)))
        zlagcoef=nanmean(xcovda(ii,midxcov-1:midxcov+1));
        %get leading lag peaks (ie lags < 0 s) both cor & anti
        [maxprecoef, maxprelagid]=max(xcovda(ii,1:midxcov-1),[],'omitnan');
        if ~isempty(maxprecoef) && ~isnan(maxprelagid)
            maxprelagts=tsx(maxprelagid(1));
        end
        %anti-correlated peak leading
        [minprecoef,minprelagid]=min(xcovda(ii,1:midxcov-1),[],'omitnan');
        if ~isempty(minprecoef) && ~isnan(minprelagid)
            minprelagts=tsx(minprelagid(1));
        end
        %get lagging lag peaks (ie lags < 0 s) both cor & anti
        [maxpostcoef,maxpostlagid]=max(xcovda(ii,midxcov+1:end),[],'omitnan');
        if ~isempty(maxpostcoef) && ~isnan(maxpostlagid)
            maxpostlagts=tsx(maxpostlagid(1)+midxcov);
        end
        %anti-correlated peak leading
        [minpostcoef,minpostlagid]=min(xcovda(ii,midxcov+1:end),[],'omitnan');
        if ~isempty(minpostcoef) && ~isnan(minpostlagid)
            minpostlagts=tsx(minpostlagid(1)+midxcov);
        end
    end
    
    %store in dapos or daneg cell according to da value for trial
    group='';
    curridx=nan;
    %exclude nans & 0s
    %OVERLY REDUNDANT< JUST SAVE POS/NEG TRIALS IDS AND GET DATA FROM DAALL
    %2/17/2019 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
xinfo(ilfp).postrials=goodtrials(abovethres);
xinfo(ilfp).negtrials=goodtrials(decsonly);
xinfo(ilfp).freq=xcovdata{ilfp}.freqband;
xinfo(ilfp).tsx=tsx;
xinfo(ilfp).xcovda=xcovda;


end

end
