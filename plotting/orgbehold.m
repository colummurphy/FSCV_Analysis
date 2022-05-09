function binfo=orgbeh(bdata,plotparam,eventtype,seltrials,winids,varargin)
ratelfp=ceil(plotparam.ratelfp);      %original bdata rate
rate=ceil(plotparam.samplespersec);            %rate of winids
argnum=1;
win=[];
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
if contains(eventtype,'targ')
    winidsl=find(ismember(seltrials,lsidetrials)==1);
    for ileft=1:length(selsidel)
        tsid=round(mean(winids(winidsl(ileft),:))*ratelfp/rate+target_lrt(ileft)*ratelfp);
        if tsid+1000>length(bdata.eyedata(selsidel(ileft),:)) || isnan(target_lrt(ileft))
            eyed(ileft)=nan;
        else
            leyed(ileft)=nanmean(bdata.eyedata(selsidel(ileft),tsid:tsid+1000));  %1000 ms average post
        end
        leyextrace(ileft,:)=bdata.eyexdata(selsidel(ileft),winids(winidsl(ileft),1)*ratelfp/rate:winids(winidsl(ileft),end)*ratelfp/rate);
        leyedtrace(ileft,:)=bdata.eyedata(selsidel(ileft),winids(winidsl(ileft),1)*ratelfp/rate:winids(winidsl(ileft),end)*ratelfp/rate);
    end
    winidsr=find(ismember(seltrials,rsidetrials)==1);
    for iright=1:length(selsider)
        tsid=round(mean(winids(winidsr(iright),:))*ratelfp/rate+target_rrt(iright)*ratelfp);
         if tsid+1000>length(bdata.eyedata(selsider(iright),:)) || isnan(target_rrt(iright))
            eyed(iright)=nan;
        else
        reyed(iright)=nanmean(bdata.eyedata(selsider(iright),tsid:tsid+1000));  %200 ms average
         end
         reyextrace(iright,:)=bdata.eyexdata(selsider(iright),winids(winidsr(iright),1)*ratelfp/rate:winids(winidsr(iright),end)*ratelfp/rate);
         reyedtrace(iright,:)=bdata.eyedata(selsider(iright),winids(winidsr(iright),1)*ratelfp/rate:winids(winidsr(iright),end)*ratelfp/rate);
    end
end

if strcmp(eventtype,'fix')
    for ifix=1:length(seltrials)
        tsid=round(mean(winids(ifix,:))*ratelfp/rate+fix_rt(ifix)*ratelfp);
        if tsid+1000>length(bdata.eyedata(seltrials(ifix),:)) || isnan(fix_rt(ifix))
            eyed(ifix)=nan;
        else
            eyed(ifix)=nanmean(bdata.eyedata(seltrials(ifix),tsid:tsid+1000));  %200 ms average
        end
        eyextrace(ifix,:)=bdata.eyexdata(seltrials(ifix),winids(ifix,1)*ratelfp/rate:winids(ifix,end)*ratelfp/rate);
        eyedtrace(ifix,:)=bdata.eyedata(seltrials(ifix),winids(ifix,1)*ratelfp/rate:winids(ifix,end)*ratelfp/rate);
    end
end
%pulse/lick arnd event window
pulse=[];
lickpre=[];
lickpost=[];
licktrace=[];
pulsetrace=[];
for iev=1:length(seltrials)
    tsid=round(mean(winids(iev,:)));
    pulse(iev)=nanmean(bdata.pulsedata(seltrials(iev),tsid:tsid+20));  %2000 ms + win
    pulsetrace(iev,:)=bdata.pulsedata(seltrials(iev),winids(iev,1):winids(iev,end));
    tsid=round(mean(winids(iev,:))*ratelfp/rate);
    lickpost(iev)=nanmean(bdata.lickdata(seltrials(iev),tsid:tsid+1000));  %1000 ms post
    lickpre(iev)=nanmean(bdata.lickdata(seltrials(iev),tsid-1000:tsid));  %1000 ms pre
    licktrace(iev,:)=bdata.lickdata(seltrials(iev),winids(iev,1)*ratelfp/rate:winids(iev,end)*ratelfp/rate);
end

binfo.seltrials=seltrials;
binfo.seltrialsl=selsidel;
binfo.seltrialsr=selsider;
binfo.target_lrt=target_lrt;
binfo.target_rrt=target_rrt;
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
