function xinfo=xclust(xcovdata,plotparam,varargin)
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
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'alllags'
            %do not group lags/pos/neg clusters, since already defined in
            %burst groups
            alllag=1;
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

%get trials only when lags not nan
notnantrials=find(~isnan(xcovlaganti(:,1))  & ~isnan(xcovlag(:,1)));    %ids from seltrials
%get trials that don't have peaks above  cutoff exclusion threshold
gpos=find(xcovlag(:,2)<=nanmean(xcovlag(:,2))*thres);
gneg=find(xcovlaganti(:,2)>=nanmean(xcovlaganti(:,2))*thres);
gsel=intersect(notnantrials,unique([gpos; gneg]));
goodtrials=seltrials(gsel);     %trial ids for orig data

%just get good trials without nan lags
da=xcovdata{ilfp}.dadata(goodtrials,:);
lfp=xcovdata{ilfp}.lfpdata(goodtrials,:);
lag=xcovlag(gsel,:);
lagneg=xcovlaganti(gsel,:);
winids=winids(gsel,:);
xcovda=xdata(gsel,:);

%get mean lag
meanlag=median(lag(:,1));
meanlagneg=median(lagneg(:,1));

%numshift to align to peak xcov
shiftpos=round((meanlag-lag(:,1))*rate);
shiftneg=round((meanlagneg-lagneg(:,1))*rate);

%shifted signals to mean lag
xcovdashiftpos=[];
xcovdashiftneg=[];
for ii=1:size(xcovda,1)
    xcovdashiftpos(ii,:)=circshift(xcovda(ii,:),shiftpos(ii),2);
    xcovdashiftneg(ii,:)=circshift(xcovda(ii,:),shiftneg(ii),2);
end

%mean aligned signal
alignedpos=mean(xcovdashiftpos);
alignedneg=mean(xcovdashiftneg);

%full width at half max
pospeak=max(alignedpos);
negpeak=min(alignedneg);
pospeakid=find(alignedpos==pospeak);
negpeakid=find(alignedneg==negpeak);
hwposidleft=find(alignedpos(1:pospeakid)<=(1-fwidth)*pospeak);
%hwposidleft=find(alignedpos(1:pospeakid)<=0);
hwposidleft=hwposidleft(end);
hwposidright=find(alignedpos(pospeakid:end)<=(1-fwidth)*pospeak);
%hwposidright=find(alignedpos(pospeakid:end)<=0);
hwposidright=hwposidright(1)+pospeakid-1;
hwnegidleft=find(alignedneg(1:negpeakid)>=(1-fwidth)*negpeak);
%hwnegidleft=find(alignedneg(1:negpeakid)>=0);
hwnegidleft=hwnegidleft(end);
hwnegidright=find(alignedneg(negpeakid:end)>=(1-fwidth)*negpeak);
%hwnegidright=find(alignedneg(negpeakid:end)>=0);
hwnegidright=hwnegidright(1)+negpeakid-1;
poshwl=pospeakid-hwposidleft;       %mean rise width in samples
poshwr=hwposidright-pospeakid;      %fall width in samples
neghwl=negpeakid-hwnegidleft;       %rise width in samples
neghwr=hwnegidright-negpeakid;      %fall width in samples

%find xcov trials with peaks in mean peak zone +/- half widths
%ie. find "shifts" that are less than hw
poswnmean=find((shiftpos<=0 & shiftpos>=-poshwl) | (shiftpos>=0 & shiftpos<=poshwr));
negwnmean=find((shiftneg<=0 & shiftneg>=-neghwl) | (shiftneg>=0 & shiftneg<=neghwr));

%determine if peak of each trial is negative or positive (abs)
abss=max(abs(xcovda),[],2);
peakids=[];
peakispos=[];
for itrial=1:size(xcovda,1)
    peakids(itrial)=find(abs(xcovda(itrial,:))==abss(itrial));
    peakispos(itrial)=xcovda(itrial,peakids(itrial))>0;
end
peakpos=find(peakispos==1);
peakneg=find(peakispos==0);
posids=intersect(peakpos,poswnmean);
negids=intersect(peakneg,negwnmean);

if alllag
    %do not cluster
    posids=peakpos;
    negids=peakneg;
end
    
%lags/xcov's of defined xcov's
%previously used poswnmean to find pos peaks for all trials
%now use posids for only trials where peak of abs is indeed pos
lagsdefpos=lag(posids,:);
lagsdefneg=lagneg(negids,:);
%lagsdefpos=lag(poswnmean,:);
%lagsdefneg=lagneg(negwnmean,:);
ratiopos=length(lagsdefpos)/length(lag);
rationeg=length(lagsdefneg)/length(lag);
xpos=xcovda(posids,:);
xneg=xcovda(negids,:);
%xpos=xcovda(poswnmean,:);
%xneg=xcovda(negwnmean,:);

%get waveform variances
varwf=var(xcovda,0,1);
varwfposaln=var(xcovdashiftpos,0,1);        %pos peak aligned
varwfnegaln=var(xcovdashiftneg,0,1);        %pos peak aligned

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
xinfo(ilfp).ratiopos=ratiopos;    %ratio pos peak xcovs within defined lag range
xinfo(ilfp).rationeg=rationeg;    %ratio neg peak xcovs within defined lag range
xinfo(ilfp).lagsdefpos=lagsdefpos;    %positive peak xcov lags within defined lag range
xinfo(ilfp).lagsdefneg=lagsdefneg;    %negativ peak xcov lags within defined lag range
xinfo(ilfp).xpos=xpos;            %x cov traces pos peak within defined lag range
xinfo(ilfp).xneg=xneg;    %x cov traces neg peak within defined lag range
xinfo(ilfp).poshwl=poshwl./rate;      %half rise width (seconds);
xinfo(ilfp).poshwr=poshwr./rate;      %half fall  width (seconds);
xinfo(ilfp).neghwl=neghwl./rate;      %half rise width (seconds);
xinfo(ilfp).neghwr=neghwr./rate;      %half fall width (seconds);
xinfo(ilfp).goodtrials=goodtrials;        %based on original trial id's
xinfo(ilfp).trialspos=seltrials(poswnmean);   %based on original trial id's
xinfo(ilfp).trialsneg=seltrials(negwnmean);   %based on original trial id's
xinfo(ilfp).freq=xcovdata{ilfp}.freqband;
xinfo(ilfp).tsx=tsx;
xinfo(ilfp).lagspos=lag;
xinfo(ilfp).lagsneg=lagneg;
xinfo(ilfp).xcovda=xcovda;
xinfo(ilfp).varwf=varwf;
xinfo(ilfp).varwfposaln=varwfposaln;
xinfo(ilfp).varwfnegaln=varwfnegaln;
xinfo(ilfp).xcovdashiftpos=xcovdashiftpos;
xinfo(ilfp).xcovdashiftneg=xcovdashiftneg;
xinfo(ilfp).posids=posids;
xinfo(ilfp).negids=negids;


%behavior data if needed
if ~isempty(bdata)
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
