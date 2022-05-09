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

%get da trials where positive da increase
%previously determined in plotxvar not most accurate way

targda=xcovdata{ilfp}.posttargda(gsel);
posdastrials=find(targda>0);
negdastrials=find(targda<0);
dapos=daaln(posdastrials,:);
lfppos=lfpaln(posdastrials,:);
daneg=daaln(negdastrials,:);
lfpneg=lfpaln(negdastrials,:);
xinfo(ilfp).dapos.trials=goodtrials(posdastrials);
xinfo(ilfp).dapos.datracesaln=dapos;
xinfo(ilfp).dapos.lfptracesaln=lfppos;
xinfo(ilfp).dapos.mididx=meanalnwin;
xinfo(ilfp).daneg.trials=goodtrials(negdastrials);
xinfo(ilfp).daneg.datracesaln=daneg;
xinfo(ilfp).daneg.lfptracesaln=lfpneg;
xinfo(ilfp).daneg.mididx=meanalnwin;

%get characteristics of da/lfp peaks for pos/neg groups
countpos=0;
countneg=0;
params={'damaxts','damints','damax','damin','darisets','dafallts',...
    'dafallts','damaxrisets','damaxrise','lfpmaxts','lfpmints',...
    'lfpmax','lfpmin','lfpmaxpostts','lfpmaxpost',...
    'zlagcoef','maxprecoef','maxprelagts','minprecoef','minprelagts',...
    'maxpostcoef','maxpostlagts','minpostcoef','minpostlagts'};
%initialize parameters storage w/nans default
for ii=1:length(params)
    xinfo(ilfp).dapos=setfield(xinfo(ilfp).dapos,params{ii},nan(1,size(dapos,1)));
    xinfo(ilfp).daneg=setfield(xinfo(ilfp).daneg,params{ii},nan(1,size(daneg,1)));
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

for ii=1:size(targda,2)
    %da dynamics
    damax=max(dawin(ii,curmid:end),[],'omitnan');
    %damaxts=find(dawin(ii,curmid:end)==damax)+meanalnwin-1;
    damaxts=find(dawin(ii,curmid:end)==damax)+curmid-1;      %relative to wins(1)
    %1s window around aln idx for minima
    damin=min(dawin(ii,curmid-10:curmid+10),[],'omitnan');
    damints=find(dawin(ii,curmid-10:curmid+10)==damin)+curmid-10-1;
    %slope peaks
    dafallts=find(dadiff(ii,curmid-10:curmid+10)==...
        min(dadiff(ii,curmid-10:curmid+10)))+curmid-10-1;
    if length(dafallts)>1
        dafallts=dafallts(1);
    end
    darisets=find(dadiff(ii,curmid:end)==...
        max(dadiff(ii,curmid:end),[],'omitnan'))+curmid-1;
    if length(darisets)>1
        darisets=darisets(1);
    end
    %max da immediately after max slope & before fall
    falls=find(dadiff(ii,darisets:end)<0)+darisets-1;
    damaxrise=nan;
    damaxchange=nan;
    damaxrisets=nan;
    if ~isempty(falls)
        damaxrise=max(dawin(ii,darisets:falls(1)+1));   %max da after max slope/plateau
        damaxchange=damaxrise-damin;        %max minus minimum previouslyf ound
        damaxrisets=find(dawin(ii,darisets:falls(1)+1)==damaxrise)+darisets-1;
    end
    
    %beta-lfp dynamics
    lfpwinend=darisets;     %set da local peak as ts end to search beta burst preceding 
    if damaxts<darisets
        lfpwinend=damaxts;
    end
    %get max beta 'preceding' da max diff or closest beta may not be
    %preceding
    %first get point closest preceding that exceeds baseline level
    abovethres=find(lfpwin(ii,1:lfpwinend)>lfpbase);
    lfpmax=max(lfpwin(ii,1:lfpwinend),[],'omitnan');     %default max, any preceding
    if ~isempty(abovethres)
        closets=abovethres(end);        %last ts, closest to darisets
        localpeak=find(lfpdiff(ii,1:closets)>0);   %find plateau, when diff nearly above 0
        if ~isempty(localpeak) && localpeak(end)<=curmid+5
            %still must be within 0.5 s of mididx, otherwise arbitrary max 
            %maybe background stronger in trial
            peaksat=localpeak(end)+1;        %when actually dips, is when peak occurs usually
            lfpmax=lfpwin(ii,peaksat);
        end
    end
    swin=1:lfpwinend;
    if curmid+6>lfpwinend
        swin=1:curmid+6;
    end
    swin=swin(swin<=length(lfpwin(ii,:)));
    lfpmaxts=find(lfpwin(ii,swin)==lfpmax);

    %get max beta suppression (min) (1.5s after mid win since sometimes delayed after da max)
    %, but after lfp max (ie. suppression after burst)
    swin=lfpmaxts:curmid+15;
    lfpmints=nan;
    lfpmin=nan;
    if ~isempty(swin)
        swin=swin(swin<=length(lfpwin(ii,:)));
        lfpmin=min(lfpwin(ii,swin),[],'omitnan');
        lfpmints=find(lfpwin(ii,swin)==lfpmin)+lfpmaxts-1;  
        %get beta burst post da max diff post suppression
        if lfpmints>lfpwinend
            %if beta suppression ts occurs later than da rise
            %new starting search ts
            lfpwinend=lfpmints;
        end
    end
    lfpmaxpost=max(lfpwin(ii,lfpwinend:end),[],'omitnan');     %default max, any subsequent
    lfpmaxpostts=find(lfpwin(ii,lfpwinend:end)==lfpmaxpost)+lfpwinend-1;
    
    %get xcov parameters
    %get mean xcoeff near zero, (ie mean xcoeff at 0+/-200 ms
    midxcov=median(1:size(xcovda,2));
    zlagcoef=nanmean(xcovda(ii,midxcov-2:midxcov+2));
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
    group=[];
    curridx=nan;
    %exclude nans & 0s
    if targda(ii)>0
        group='dapos';
        countpos=countpos+1;
        curridx=countpos;
    elseif targda(ii)<0
        group='daneg';
        countneg=countneg+1;
        curridx=countneg;
    end
    if ~isempty(group)
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'damax',{curridx},damax);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'damaxts',{curridx},damaxts+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'damin',{curridx},damin);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'damints',{curridx},damints+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'dafallts',{curridx},dafallts+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'darisets',{curridx},darisets+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'damaxrise',{curridx},damaxchange);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'damaxrisets',{curridx},damaxrisets+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'lfpmaxts',{curridx},lfpmaxts+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'lfpmax',{curridx},lfpmax);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'lfpmints',{curridx},lfpmints+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'lfpmin',{curridx},lfpmin);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'lfpmaxpostts',{curridx},lfpmaxpostts+winoffset);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'lfpmaxpost',{curridx},lfpmaxpost);

        xinfo(ilfp)=setfield(xinfo(ilfp),group,'zlagcoef',{curridx},zlagcoef);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'maxprecoef',{curridx},maxprecoef);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'maxprelagts',{curridx},maxprelagts);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'minprecoef',{curridx},minprecoef);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'minprelagts',{curridx},minprelagts);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'maxpostcoef',{curridx},maxpostcoef);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'maxpostlagts',{curridx},maxpostlagts);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'minpostcoef',{curridx},minpostcoef);
        xinfo(ilfp)=setfield(xinfo(ilfp),group,'minpostlagts',{curridx},minpostlagts);
    end
end

%get ts differences between da peaks and beta bursts/suppresions
dagroups={'dapos','daneg'};
for ii=1:length(dagroups)
    lfpmax=getfield(xinfo(ilfp),dagroups{ii},'lfpmaxts');
    lfpmin=getfield(xinfo(ilfp),dagroups{ii},'lfpmints');
    damax=getfield(xinfo(ilfp),dagroups{ii},'damaxts');
    delt_lfpmax_damax=damax-lfpmax;
    delt_lfpmin_damax=damax-lfpmin;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmax_damax',delt_lfpmax_damax);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmin_damax',delt_lfpmin_damax);
    darise=getfield(xinfo(ilfp),dagroups{ii},'darisets');
    delt_lfpmax_darise=darise-lfpmax;
    delt_lfpmin_darise=darise-lfpmin;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmax_darise',delt_lfpmax_darise);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmin_darise',delt_lfpmin_darise);
    damin=getfield(xinfo(ilfp),dagroups{ii},'damints');
    delt_lfpmax_damin=damin-lfpmax;
    delt_lfpmin_damin=damin-lfpmin;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmax_damin',delt_lfpmax_damin);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_lfpmin_damin',delt_lfpmin_damin);
    lfpmaxpost=getfield(xinfo(ilfp),dagroups{ii},'lfpmaxpostts');    %beta rebound
    delt_damin_lfpmaxpost=lfpmaxpost-damin;
    delt_darise_lfpmaxpost=lfpmaxpost-darise;
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_damin_lfpmaxpost',delt_damin_lfpmaxpost);
    xinfo(ilfp)=setfield(xinfo(ilfp),dagroups{ii},'delt_darise_lfpmaxpost',delt_darise_lfpmaxpost);

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

%behavior data if needed
if ~isempty(bdata) && ilfp==1
    %get behaviors for both burst/no burst groups based on lfp pair
    %reaction times (left/right split)
    rsidetrials=find(bdata.rewardside==1);
    lsidetrials=find(bdata.rewardside==0);
    selsidel=intersect(goodtrials,lsidetrials);
    selsider=intersect(goodtrials,rsidetrials);
    target_lrt=bdata.target_rt(selsidel);
    if any(target_lrt<=0)
        target_lrt(find(target_lrt<=0))=nan;
    end

    target_rrt=bdata.target_rt(selsider);
    if any(target_rrt<=0)
        target_rrt(find(target_rrt<=0))=nan;
    end
    fix_rt=bdata.fix_rt(goodtrials);
    if any(fix_rt<=0)
        fix_rt(find(fix_rt<=0))=nan;
    end
    %pupil d (left/right split) when eye at indicated cue
    leyed=[];
    reyed=[];
    eyed=[];
    if strcmp(eventtype,'intertarg')
        winidsl=find(ismember(goodtrials,lsidetrials)==1);
        for ileft=1:length(selsidel)
            tsid=round(mean(winids(winidsl(ileft),:))*ratelfp/rate+target_lrt(ileft)*ratelfp);
            if tsid+1000>length(bdata.eyedata(selsidel(ileft),:)) || isnan(target_lrt(ileft))
                eyed(ileft)=nan;
            else
            leyed(ileft)=nanmean(bdata.eyedata(selsidel(ileft),tsid:tsid+1000));  %1000 ms average post
            end
        end
        winidsr=find(ismember(goodtrials,rsidetrials)==1);
        for iright=1:length(selsider)
            tsid=round(mean(winids(winidsr(iright),:))*ratelfp/rate+target_rrt(iright)*ratelfp);
             if tsid+1000>length(bdata.eyedata(selsider(iright),:)) || isnan(target_rrt(iright))
                eyed(iright)=nan;
            else
            reyed(iright)=nanmean(bdata.eyedata(selsider(iright),tsid:tsid+1000));  %200 ms average
             end
        end
    end
    if strcmp(eventtype,'interfix')
        for ifix=1:length(goodtrials)
            tsid=round(mean(winids(ifix,:))*ratelfp/rate+fix_rt(ifix)*ratelfp);
            if tsid+1000>length(bdata.eyedata(goodtrials(ifix),:)) || isnan(fix_rt(ifix))
                eyed(ifix)=nan;
            else
                eyed(ifix)=nanmean(bdata.eyedata(goodtrials(ifix),tsid:tsid+1000));  %200 ms average
            end
        end
    end
    %pulse/lick arnd event window
    pulse=[];
    lickpre=[];
    lickpost=[];
    for iev=1:length(goodtrials)
        tsid=round(mean(winids(iev,:)));
        pulse(iev)=nanmean(bdata.pulsedata(goodtrials(iev),tsid:tsid+20));  %2000 ms + win
        tsid=round(mean(winids(iev,:))*ratelfp/rate);
        lickpost(iev)=nanmean(bdata.lickdata(goodtrials(iev),tsid:tsid+1000));  %1000 ms post
        lickpre(iev)=nanmean(bdata.lickdata(goodtrials(iev),tsid-1000:tsid));  %1000 ms pre
    end
    xinfo(ilfp).target_lrt=target_lrt;
    xinfo(ilfp).target_rrt=target_rrt;
    xinfo(ilfp).fix_rt=fix_rt;
    xinfo(ilfp).leyed=leyed;
    xinfo(ilfp).reyed=reyed;
    xinfo(ilfp).eyed=eyed;
    xinfo(ilfp).pulse=pulse;
    xinfo(ilfp).lickpre=lickpre;
    xinfo(ilfp).lickpost=lickpost;
end


end

end
