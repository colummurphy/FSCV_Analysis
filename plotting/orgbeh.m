function binfo=orgbeh(bdata,plotparam,eventtype,seltrials,winids,varargin)
ratelfp=ceil(plotparam.ratelfp);      %original bdata rate
rate=ceil(plotparam.samplespersec);            %rate of winids
argnum=1;
win=[];
avgdur=1*ratelfp;        %samples to average
while argnum<=1
    switch varargin{argnum}
        case 'win'
            %user specifies new win to extract data rath tan provided by 
            %xcov
            argnum=argnum+1;
            win=varargin{argnum};
    end
    argnum=argnum+1;
end
midids=median(winids,2);
if ~isempty(win)
    winids=[];
    %expand or shrink provided winids by value provided in win
    for ii=1:length(midids)
        winids(ii,:)=midids(ii)+win(1)*rate:midids(ii)+win(2)*rate;
    end
end
ts=bdata.relts;
tsr=round(ts*ratelfp)./ratelfp;       %ts aligned to outcome for each trial synchronized to ts listed for ttarg/etc event marks
            
%get behaviors reaction times (left/right split)
%seltrials is ids for good trials for original data with all trials
%winids is win ids for targeted event, just seltrials
rsidetrials=find(bdata.rewardside==1);
lsidetrials=find(bdata.rewardside==0);
selsidel=intersect(seltrials,lsidetrials);
selsider=intersect(seltrials,rsidetrials);
target_lrt=bdata.target_rt(selsidel);
if any(target_lrt<=0)
    target_lrt(find(target_lrt<=0))=nan;
end

target_rrt=bdata.target_rt(selsider);
if any(target_rrt<=0)
    target_rrt(find(target_rrt<=0))=nan;
end
fix_rt=bdata.fix_rt(seltrials);
if any(fix_rt<=0)
    fix_rt(find(fix_rt<=0))=nan;
end
%all rts' in order of trials irregardless of side
rts=bdata.target_rt(seltrials);

%pupil d (left/right split) when eye at indicated cue
leyed=[];
reyed=[];
eyed=[];
eyedist=[];
eyedisttrace=[];
leyedisttrace=[];
reyedisttrace=[];
eyextrace=[];
leyextrace=[];
reyextrace=[];
leyedtrace=[];
reyedtrace=[];
eyedtrace=[];
alignidx=median(1:size(winids,2));  %middle of win set as align ts during xcov
meanalnwin=round(mean(winids(:,alignidx)));
winlen=(meanalnwin-winids(1))./rate*ratelfp;  %in samples
ttids=[];
if ~isempty(bdata.eyedata) 
if contains(eventtype,'targ')
    ttids=bdata.ttarg;
    if strcmp(eventtype,'targeye')
        ttids=bdata.ttargeye;
    end
    lids=find(ismember(seltrials,lsidetrials)==1);
    for ileft=1:length(selsidel)
        tsid=ttids(selsidel(ileft));            %time stamp of event
        tsidr=round(tsid*ratelfp)./ratelfp;        %rounded to nearest 1/rate
        %tsid=round(mean(winids(winidsl(ileft),:))*ratelfp/rate+target_lrt(ileft)*ratelfp);
        if tsidr*ratelfp+avgdur>length(bdata.eyedata(selsidel(ileft),:)) || isnan(target_lrt(ileft))
            eyed(lids(ileft))=nan;
            leyed(ileft)=nan;
        else
            leyed(ileft)=nanmean(bdata.eyedata(selsidel(ileft),tsidr*ratelfp:tsidr*ratelfp+avgdur));  %1000 ms average post
            eyed(lids(ileft))=leyed(ileft);
        end
        leyextrace(ileft,:)=bdata.eyexdata(selsidel(ileft),tsidr*ratelfp-winlen:tsidr*ratelfp+winlen);
        leyedtrace(ileft,:)=bdata.eyedata(selsidel(ileft),tsidr*ratelfp-winlen:tsidr*ratelfp+winlen);
        eyedtrace(lids(ileft),:)=leyedtrace(ileft,:);
        eyextrace(lids(ileft),:)=leyextrace(ileft,:);
    end
    rids=find(ismember(seltrials,rsidetrials)==1);
    for iright=1:length(selsider)
        tsid=ttids(selsider(iright));            %time stamp of event
        tsidr=round(tsid*ratelfp)./ratelfp;        %rounded to nearest 1/rate
         if tsidr*ratelfp+avgdur>length(bdata.eyedata(selsider(iright),:)) || isnan(target_rrt(iright))
            eyed(rids(iright))=nan;
            reyed(iright)=nan;
         else
            reyed(iright)=nanmean(bdata.eyedata(selsider(iright),tsidr*ratelfp:tsidr*ratelfp+avgdur));  %1000 ms average post
         end
         reyextrace(iright,:)=bdata.eyexdata(selsider(iright),tsidr*ratelfp-winlen:tsidr*ratelfp+winlen);
         reyedtrace(iright,:)=bdata.eyedata(selsider(iright),tsidr*ratelfp-winlen:tsidr*ratelfp+winlen);
         eyedtrace(rids(iright),:)=reyedtrace(iright,:);
        eyextrace(rids(iright),:)=reyextrace(iright,:);
    end
end

if contains(eventtype,'fix')    
    ttids=bdata.tfix;
    if strcmp(eventtype,'fixeye')
        ttids=bdata.tfixeye;
    end
    for ifix=1:length(seltrials)
        tsid=ttids(seltrials(ifix));
        tsidr=round(tsid*ratelfp)./ratelfp;        %rounded to nearest 1/rate
        if tsidr*ratelfp+avgdur>length(bdata.eyedata(seltrials(ifix),:)) || isnan(fix_rt(ifix))
            eyed(ifix)=nan;
        else
            eyed(ifix)=nanmean(bdata.eyedata(seltrials(ifix),tsidr*ratelfp:tsidr*ratelfp+avgdur));
        end
        eyextrace(ifix,:)=bdata.eyexdata(seltrials(ifix),tsidr*ratelfp-winlen:tsidr*ratelfp+winlen);
        eyedtrace(ifix,:)=bdata.eyedata(seltrials(ifix),tsidr*ratelfp-winlen:tsidr*ratelfp+winlen);
    end
end
if strcmp(eventtype,'outcome')    
    ttids=meanalnwin./rate;         %since all trials aligned to outcome at 30s
    for it=1:length(seltrials)
        if ttids*ratelfp+avgdur>length(bdata.eyedata(seltrials(it),:))
            eyed(it)=nan;
        else
            eyed(it)=nanmean(bdata.eyedata(seltrials(it),ttids*ratelfp:ttids*ratelfp+avgdur));
        end
        eyextrace(it,:)=bdata.eyexdata(seltrials(it),ttids*ratelfp-winlen:ttids*ratelfp+winlen);
        eyedtrace(it,:)=bdata.eyedata(seltrials(it),ttids*ratelfp-winlen:ttids*ratelfp+winlen);
    end
end
end
%pulse/lick arnd event window
pulse=[];
lickpre=[];
lickpost=[];
licktrace=[];
pulsetrace=[];
tsid=[];
if ~isempty(bdata.pulsedata) 
for iev=1:length(seltrials)
    if strcmp(eventtype,'outcome')
        tsid=ttids;
    else
        tsid=ttids(seltrials(iev));
    end
    tsidr=round(tsid*rate)./rate; 
    pulse(iev)=nanmean(bdata.pulsedata(seltrials(iev),tsidr*rate:tsidr*rate+avgdur/ratelfp*rate*2));
    pulsetrace(iev,:)=bdata.pulsedata(seltrials(iev),tsidr*rate-winlen/ratelfp*rate:tsidr*rate+winlen/ratelfp*rate);
    samplength=length(tsidr*rate-winlen/ratelfp*rate:tsidr*rate+winlen/ratelfp*rate);
    if any(pulsetrace(iev,:)>230 | pulsetrace(iev,:)<60)
        %non physiological pulse remove
        pulsetrace(iev,:)=nan(1,samplength);
    end
end
end
if  ~isempty(bdata.lickdata)
for iev=1:length(seltrials)
    if strcmp(eventtype,'outcome')
        tsid=ttids;
    else
        tsid=ttids(seltrials(iev));
    end
    tsidr=round(tsid*ratelfp)./ratelfp; 
    lickpost(iev)=nanmean(bdata.lickdata(seltrials(iev),tsidr*ratelfp:tsidr*ratelfp+avgdur));
    lickpre(iev)=nanmean(bdata.lickdata(seltrials(iev),tsidr*ratelfp-avgdur:tsidr*ratelfp));
    licktrace(iev,:)=bdata.lickdata(seltrials(iev),tsidr*ratelfp-winlen:tsidr*ratelfp+winlen);
end
end


binfo.evt=eventtype;
binfo.ttids=ttids;
binfo.seltrials=seltrials;
binfo.seltrialsl=selsidel;
binfo.seltrialsr=selsider;
binfo.target_lrt=target_lrt;
binfo.target_rrt=target_rrt;
binfo.target_rts=rts;
binfo.fix_rt=fix_rt;
binfo.leyed=leyed;
binfo.reyed=reyed;
binfo.eyed=eyed;
binfo.pulse=pulse;
binfo.lickpre=lickpre;
binfo.lickpost=lickpost;
%raw data vs time around cue event
binfo.eyedtrace=eyedtrace;
binfo.eyextrace=eyextrace;
binfo.leyextrace=leyextrace;
binfo.reyextrace=reyextrace;
binfo.leyedtrace=leyedtrace;
binfo.reyedtrace=reyedtrace;
binfo.licktrace=licktrace;
binfo.pulsetrace=pulsetrace;

end
