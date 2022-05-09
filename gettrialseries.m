function [listout,sitemap]=gettrialseries(sessnum,trlists,xinfos,xbinfos,binfos,datms,lfptms,plotparam,varargin)

sessiontypes={'bigreward','smallreward'};
fs=10;          %rate for all xinfos & xbinfos variables
sessid=num2str(sessnum);
event='targ';
trlistid=find(strcmp({trlists.sessid},sessid));
trlist=trlists(trlistid).list;
numtrials=length(trlist);
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,['trialseries_' sessid '_']);
sessnums=plotparam.sessnums;
targdasites=plotparam.dasites;
sites=getsites(sessnum,targdasites);
targdasites={sites.site};
targdanames={sites.probeid};
targdachs={sites.ch};
[dapair,lfppair]=getsitepairs(targdanames,'numpairs',2);     %just get first 2 prioritized lfp pairs
targlfpnames=unique(lfppair);
sitemap=sites;
lfpsites=getlfpsites(sessnum,lfppair);
for ix=1:length(sitemap)
    curda=sitemap(ix).probeid;
    curids=find(strcmp(dapair,curda));
    sitemap(ix).lfp=lfppair(curids);
    curlfps=find(contains({lfpsites.probeid},lfppair(curids)));
    sitemap(ix).lfpsite={lfpsites(curlfps).site};
    sitemap(ix).lfpids=find(contains(targlfpnames,lfppair(curids)));
end
targses=find(strcmp({trialgrps.sessid},sessid));
trialinfo=trialgrps(targses).trialinfo;

%initialize trlist
for ix=1:length(trlist)
    trlist(ix).side=nan;
    trlist(ix).rts=nan;
     trlist(ix).lick=nan;
     trlist(ix).eye=nan;
          trlist(ix).eyelate=nan;

    trlist(ix).hr=nan;
    trlist(ix).rrstd=nan;
        trlist(ix).da=[];
    trlist(ix).lfp=[];  
        trlist(ix).lfplate=[];  

end

argnum=1;
siteda='';
sitelfp='';
bflag=0;
linfos={};
metric='targpeak';
metriclfp='targwin';

hrwin=[0 4];
rrwin=[0 4];
eyewin=[0 1];       %3 to 4 has inverse characteristics now using z-score
eyewinlate=[3 4];
while argnum<=length(varargin)
    switch varargin{argnum}        
        case 'metric'
            %for datm values
            argnum=argnum+1;
            metric=varargin{argnum};
        case 'metriclfp'
            %for betatm values
            argnum=argnum+1;
            metriclfp=varargin{argnum};        
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};     %alignment event/window for values
        case 'tlim'
            argnum=argnum+1;
            triallim=varargin{argnum};
    end
    argnum=argnum+1;
end     
%compile data
%first get good trials list for trial type big/small/targ
maplist={};     %listing of rows mapped by sel trials; ie. seltrials{1}(1) --> maplist{1}(1) = targ row for specific type & trial
seltrials={};
dtm={};
lfptm=[];
xtargda=[];
xtarglfp=[];
for itype=1:length(sessiontypes)
    goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
    incswgrp=find(strcmp(trialinfo(itype).names,'switch')==1);
    goodtrials=sort(unique([trialinfo(itype).nums{goodtrialgrp} trialinfo(itype).nums{incswgrp}]));  
    seltrials{itype}=goodtrials;
    tlabelsgood=goodtrials+99;              %convert to trial labels ie 100,101,etc.
    typelistids=find(strcmp({trlist.type},sessiontypes{itype}));
    maplist(itype).trlist=find(ismember([trlist.id],tlabelsgood) & ...
        strcmp({trlist.type},sessiontypes{itype}));
    maplist(itype).trialnums=goodtrials;
    maplist(itype).type=sessiontypes{itype};
    
%get side of target
    leftids=find(contains(trialinfo(itype).names,'left')==1);
    lefttrials=sort(unique(trialinfo(itype).nums{leftids}));
    rightids=find(contains(trialinfo(itype).names,'right')==1);
    righttrials=sort(unique(trialinfo(itype).nums{rightids}));
    for ida=1:length(maplist(itype).trlist)
        if ismember(maplist(itype).trialnums(ida),righttrials)
            trlist(maplist(itype).trlist(ida)).side=1;
        else
            trlist(maplist(itype).trlist(ida)).side=0;
        end
    end

%get reaction times from binfos
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},'fix'));            %event shouldn't matter since getting rts defined at target, fix may have more "good" trials where calculated
    rts=binfos(btarg).target_rts;
    rttrials=binfos(btarg).seltrials;
    for ida=1:length(rttrials)
        tid=find(ismember(maplist(itype).trialnums,rttrials(ida)));
        if ~isempty(tid)
            trlist(maplist(itype).trlist(tid)).rts=rts(ida);
        end
    end
    
%get hr  from xbinfos
targrows=find(strcmp({xbinfos.sitelfp},'pulse') & ...
    strcmp({xbinfos.event},event) & ...
    contains({xbinfos.sessionid},sessid) &...
    contains({xbinfos.sessiontype},sessiontypes{itype}));
    btarg=targrows(1);            %event matters since where calculated
    curdata=getfield(xbinfos(btarg),'daall');
    hr=curdata.lfptracesaln;
    alnidx=curdata.mididx;
    hrtrials=curdata.trials;
    tids=find(ismember(maplist(itype).trialnums,hrtrials));
    hrdata=hr(:,hrwin(1)*fs+alnidx:hrwin(2)*fs+alnidx);
    hrtemp=nanmean(hrdata,2);
    rrtemp=1./hr.*60.*1e3; 
    rrdata=rrtemp(:,rrwin(1)*fs+alnidx:rrwin(2)*fs+alnidx);
    rrstdtemp=nanstd(rrdata,[],2);
   for ida=1:length(hrtrials)
        tid=find(ismember(maplist(itype).trialnums,hrtrials(ida)));
        if ~isempty(tid)
            trlist(maplist(itype).trlist(tid)).hr=hrtemp(ida);
            trlist(maplist(itype).trlist(tid)).rrstd=rrstdtemp(ida);            
        end
   end   
   
%get eye info from from xbinfos
targrows=find(strcmp({xbinfos.sitelfp},'eye') & ...
    strcmp({xbinfos.event},'targeye') & ...
    contains({xbinfos.sessionid},sessid) &...
    contains({xbinfos.sessiontype},sessiontypes{itype}));
glitchwidth=2;
    btarg=targrows(1);            %event matters since where calculated
    curdata=getfield(xbinfos(btarg),'daall');
    eye=curdata.lfptracesaln;
    alnidx=curdata.mididx;
    etrials=curdata.trials;
    tids=find(ismember(maplist(itype).trialnums,etrials));
    edata=-eye;
    basewin=[alnidx-7*fs:alnidx+4*fs];          %fix baseline window for eye below
    %zscore eye wrt baseline as done in olivier et al 2014 frontier
    %remove blinks
    edata=deglitchinterpsimple(edata,3.5,glitchwidth);        
    basedata=-eye(:,basewin);   %normalizing to pre-fix and targ window allows normalization trial to trial to look at faster changes within 1 sec 06/04 without phase trend eg. 67/68/69
    basedata=basedata;
    for ix=1:size(basedata,1)
        [xx, xb]=rmoutliers(basedata(ix,:));
        basedata(ix,xb)=nan;
    end
    baseline=nanstd(basedata,[],2);          %normalized to each trial so independent of phase changes
    dataz=(edata-nanmean(basedata,2))./baseline;    
    etemp=nanmean(edata(:,eyewin(1)*fs+alnidx:eyewin(2)*fs+alnidx),2);   
    etemplate=nanmean(edata(:,eyewinlate(1)*fs+alnidx:eyewinlate(2)*fs+alnidx),2); 
   for ida=1:length(etrials)
        tid=find(ismember(maplist(itype).trialnums,etrials(ida)));
        if ~isempty(tid)
            trlist(maplist(itype).trlist(tid)).eye=etemp(ida);
            trlist(maplist(itype).trlist(tid)).eyelate=etemplate(ida);            
        end
   end             

%get signal values 
    if ~isempty(targdanames)
        %targeted da values from datm
        dtm={};
        for ich=1:length(targdanames)   
            dtm{itype}=datms{targses}{itype}{targdachs{ich}};           %get all da targeted values for current trial type big/small
            datrials=dtm{itype}.trialnums;
            datatemp=getfield(dtm{itype},metric);
            for ida=1:length(datrials)
                tid=find(ismember(maplist(itype).trialnums,datrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).da(ich)=datatemp(ida);
                end
            end
        end
    end

    if ~isempty(targlfpnames) 
        %targeted beta values from betatm
        %mapped to lfpids in sitemap
        ltm={};
        for ich=1:length(targlfpnames)
            lfptm=lfptms{targses}{itype};
            for ii=1:length(lfptm)
                if strcmp(lfptm{ii}.site,targlfpnames{ich})
                    lch=ii;
                end
            end
            ltm{itype}=lfptm{lch};           %get all da targeted values for current trial type big/small           
            lfptrials=ltm{itype}.trialnums;
            ltempearly=getfield(ltm{itype},'targimwin');
            ltemplate=getfield(ltm{itype},'targwin');

            for ida=1:length(lfptrials)
                tid=find(ismember(maplist(itype).trialnums,lfptrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).lfp(ich)=ltempearly(ida);
                    trlist(maplist(itype).trlist(tid)).lfplate(ich)=ltemplate(ida);
                end
            end
        end
    end
end

%fill any nans for sides that need to be filled.  %%%FIX LATER MAYBE MORE
%ACCURATE TO USE ORIGINAL SIDE LIST OR GET THIS FROM LIST OF EVENTS RATHER
%THAN GET FROM GOOD TRIALS ONLY
filledrows=find(~isnan([trlist.side]));
for ix=1:length(trlist)
    if isnan(trlist(ix).side)
        nearestnextside=min(filledrows(filledrows>ix));
        if ~isempty(nearestnextside)
            getside=trlist(nearestnextside).side;
            trlist(ix).side=getside;
        end
    end
end

listout=trlist;

save([plotparam.savepath 'trialseries_sess' sessid],'sitemap','trlist');
%print(hf,savename,'-painters','-depsc');
end


