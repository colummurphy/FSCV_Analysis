function trialgrps=getmulttrials(plotparam,xinfos,binfo,sessnum,siteda,type,event,varargin)
fs=10;
trialgrps=[];
sessnums=plotparam.sessnums;
trgrps=plotparam.trialgrps;
if length(sessnum)>1
    %multiple arguments in cell
    argtemp=sessnum;
    sessnum=argtemp{1};
    siteda=argtemp{2};
    type=argtemp{3};
    event=argtemp{4};
    varargin=argtemp(5:end);
end
sessid=find(sessnums==sessnum);
argnum=1;
label='goodtrials';
pertype='';       %dopamine  from xinfos(x).daall.datracesaln
pertile=15;     %10 percentile
pwin=[];
lfpflag=0;
if contains(event,'fix')
    pwin=1.2;
end
if contains(event,'targ')
    pwin=3.9;      %changed from 35 to 39 2/11/2019,back to 35 2/12/
    if sessnum<67
        pwin=1.8; %NEED TO EXTEND OUT FOR SMALLER TARG SINCE PEAK OCCURS LATER IN BR both for da & beta reobund
    end
end
bwin=.1;
ttypes={};
trtable={};
trialinfo=plotparam.trialgrps(contains({plotparam.trialgrps.sessid},num2str(sessnum))).trialinfo;
infonames={'bigreward','smallreward','targetbreak','fixbreak'};
condtypes={};
trtypes={};
win=[];
fixedcond=0;
fixedgroup={};
bstatic=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dapos'
            label='postrials';
        case 'daneg'
            label='negtrials';
        case 'percentile'
            %percentile of data, peak for da, win for lfp
            argnum=argnum+1;
            pertype=varargin{argnum}{1};
            pertile=varargin{argnum}{2}; 
        case 'trtable'
            %supply trial list table with labels of which groups to plot
            %(as generaetd by gettrtable & maketrorg)
            %each {} group separated by session type
            argnum=argnum+1;
            trtable=varargin{argnum};
            argnum=argnum+1;
            ttypes=varargin{argnum}; %eg. {{'-3break','-2break','-1break'}}; for big / smalll
        case 'fixedcond'
            %fix trial type conditions, ie right only
            fixedcond=1;
            argnum=argnum+1;
            fixedgroup=varargin{argnum};
        case 'conditions'
            argnum=argnum+1;
            condtypes=varargin{argnum};
        case 'trialtypeconds'
            %grps for individual session types (ie. big vs small) & the
            %supplied conditions
            argnum=argnum+1;
            trtypes=varargin{argnum};
        case 'targimwin'
            pwin=1;
        case 'win'
            argnum=argnum+1;
            win=varargin{argnum};
            
    end
    argnum=argnum+1;
end
if contains(pertype,'lfp') || contains(pertype,'beh')
    lfpflag=1;
end
if contains(siteda,'-')
    lfpflag=1;
end
if strcmp(siteda,'trt')
    %reaction times, not time data
    targsig='target_rts';
    bstatic=1;
end
if strcmp(siteda,'rt')
    targsig='fix_rt';
    bstatic=1;
end

trialgrptypes={'big','small'};
targvals={};
allvals=[];
%multiple types
if isempty(trtable) && isempty(trtypes) && isempty(condtypes)
for itype=1:length(type)
    trialnums=[];
    targrow=0;
    trialinfo=trgrps(sessid).trialinfo;    
    if ~bstatic
    targrow=getrowxinfos(xinfos,sessnum,siteda,type{itype},event,lfpflag);
    trialnums=unique(sort([getfield(xinfos(targrow),label)]));
    end
    if fixedcond
        %fixed trial type categories only
        grptrialnums=[];
        trialgrpid=find(strcmp(trialgrptypes,type));
        for condtype=1:length(fixedgroup)
            idtype=find(contains(trialinfo(trialgrpid).names,fixedgroup(condtype)));
            if condtype==1
                grptrialnums=trialinfo(trialgrpid).nums{idtype};
            else
                grptrialnums=[grptrialnums trialinfo(trialgrpid).nums{idtype}];
            end
        end
        trialnums=grptrialnums;
    end
    %get signal values for targeted trials
    switch pertype
        case 'beh'
            dawin=[];
            %supply xbinfos instead of xinfos
            if ~bstatic
            trialnums=intersect(xinfos(targrow(1)).daall.trials,trialnums);
            trialids=find(ismember(xinfos(targrow(1)).daall.trials,trialnums));
            temp=getfield(xinfos(targrow(1)).daall,'lfptracesaln');
            alnidx=getfield(xinfos(targrow(1)).daall,'mididx');
            datraces=temp(trialids,:);
            if ~isempty(win)                
                win=[alnidx+win(1)*fs:alnidx+win(2)*fs];
            else
                win=[alnidx:alnidx+pwin*fs];
            end
            if strcmp(siteda,'eye')
                basewin=[alnidx-3*fs:alnidx+2*fs];          %fix baseline window for eye below
                %zscore eye wrt baseline as done in olivier et al 2014 frontier
                %remove blinks
                glitchwidth=2;
                datraces=-datraces;       %flip so original sign
                datraces=deglitchinterpsimple(datraces,3.5,glitchwidth);
                targrows=find(strcmp({xinfos.sitelfp},siteda) & ...
                    strcmp({xinfos.event},'fix') & ...
                    contains({xinfos.sessionid},num2str(sessnum)) &...
                    contains({xinfos.sessiontype},type{itype}));
                targrow=targrows(1);
                fixdata=getfield(xinfos(targrow),'daall');
               % seltrials=find(ismember(fixdata.trials,trialtypes.nums{ttype}) & ismember(fixdata.trials,datatrials) & ismember(fixdata.trials,trialnums(igrp).trialgrps(itype).trials));
                fixdata=getfield(fixdata,'lfptracesaln');
                basedata=fixdata(:,basewin);
                basedata=-basedata;
                 basedata=deglitchinterpsimple(basedata,3.5,4);
                for ix=1:size(basedata,1)
                        glitchpts=find(abs(basedata(ix,:))>nanstd(basedata(ix,:))*2.5+nanmean(abs(basedata(ix,:))));
                        if ~isempty(glitchpts)
                            glitchall=sort(unique([glitchpts glitchpts-1 glitchpts-2 glitchpts-3 glitchpts+1 glitchpts+2 glitchpts+3]));
                            glitchall=glitchall(glitchall>0 & glitchall<size(basedata,2));
                            basedata(ix,glitchpts)=nan;
                        end
                end
                baseline=nanmean(nanstd(basedata(:,1:5*fs),[],2));
                dataz=(datraces-nanmean(datraces,2))./baseline;
                datraces=dataz;
            end    
            dawin=nanmean(datraces(:,win),2);
            else
                    %reaction time static data
                targrows=find(strcmp({binfo.event},event) & ...
                contains({binfo.sessionid},num2str(sessnum)) &...
                contains({binfo.sessiontype},type{itype}));
                targrow=targrows(1);
                seltrials=find(ismember(binfo(targrow).seltrials,trialnums));
                targdata=getfield(binfo(targrow),targsig);
                dawin=targdata(seltrials).*1e3;       %in ms
            end
            targvals{itype}=dawin;
            allvals=[allvals; dawin];            
        case 'da-abs'
            trialnums=intersect(xinfos(targrow(1)).daall.trials,trialnums);
            trialids=find(ismember(xinfos(targrow(1)).daall.trials,trialnums));
            temp=getfield(xinfos(targrow(1)).daall,'datracesaln');
            alnidx=getfield(xinfos(targrow(1)).daall,'mididx');
            datraces=temp(trialids,:);
            win=[alnidx:alnidx+pwin*fs];
            baseids=alnidx-bwin*fs:alnidx+bwin*fs;   
            baseline=nanmean(datraces(:,baseids),2);
            dasub=datraces-baseline;
            [damax,maxid]=max(abs(dasub(:,win)),[],2);
            damax=dasub(:,maxid);
            targvals{itype}=damax;
            allvals=[allvals; damax];
        case 'da'
            trialnums=intersect(xinfos(targrow(1)).daall.trials,trialnums);
            trialids=find(ismember(xinfos(targrow(1)).daall.trials,trialnums));
            temp=getfield(xinfos(targrow(1)).daall,'datracesaln');
            alnidx=getfield(xinfos(targrow(1)).daall,'mididx');
            datraces=temp(trialids,:);
            win=[alnidx:alnidx+pwin*fs];
            baseids=alnidx-bwin*fs:alnidx+bwin*fs;   
            baseline=nanmean(datraces(:,baseids),2);
            dasub=datraces-baseline;
            [damax,maxid]=max(dasub(:,win),[],2);
            targvals{itype}=damax;
            allvals=[allvals; damax];
        case 'lfp'
            trialnums=intersect(xinfos(targrow(1)).daall.trials,trialnums);
            trialids=find(ismember(xinfos(targrow(1)).daall.trials,trialnums));
            temp=getfield(xinfos(targrow(1)).daall,'lfptracesaln');
            alnidx=getfield(xinfos(targrow(1)).daall,'mididx');
            datraces=temp(trialids,:);
            win=[alnidx:alnidx+pwin*fs];
            %baseids=alnidx-bwin*fs:alnidx+bwin*fs;   
            %baseline=nanmean(datraces(:,baseids),2);
            dasub=datraces;
            dawin=nanmean(dasub(:,win),2);
            targvals{itype}=dawin;
            allvals=[allvals; dawin];
    end

    trialgrps(1).site=siteda;
    trialgrps(1).trialgrps(itype).trials=trialnums;
    trialgrps(1).trialgrps(itype).type=xinfos(targrow(1)).sessiontype;
    if bstatic
            trialgrps(1).trialgrps(itype).type=binfo(targrow).sessiontype;
    end
    if isempty(targvals)
        trialgrps(1).cat=[label '_' event];
    end
 %   trialgrps(itype).trials=trialnums;
   % trialgrps(itype).type=xinfos(targrow(1)).sessiontype;
end
end
if strcmp(pertype,'beh')
    pertype=siteda;
end
if ~isempty(trtable)
for itype=1:length(type)
    trialnums=[];
    trialinfo=trgrps(sessid).trialinfo;
    %trialgrps(1).trialgrps(itype).trials=[];
    
    %first get good trials for current trial type big/small/targ
    goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
    goodtrials=trialinfo(itype).nums{goodtrialgrp};  
    %get trial list
    listses=find(strcmp({trtable.sessid},num2str(sessnum)));
    listtyp=find(contains({trtable(listses).trorg.type},type{itype}));
    for tt=1:length(ttypes)       
            trialgrps(tt).trialgrps(itype).trials=[];       %initialize empty array
        trlisttarg=find(contains({trtable(listses).trorg(listtyp).grp.label},ttypes{tt}));
        trialnumtemp=trtable(listses).trorg(listtyp).grp(trlisttarg).trials;
        trialnumtemp=intersect(trialnumtemp,goodtrials);
        %trialnums=[trialnums trialnumtemp];
        trialgrps(tt).site=siteda;
        trialgrps(tt).trialgrps(itype).trials=sort(unique([trialgrps(tt).trialgrps(itype).trials trialnumtemp]));
        trialgrps(tt).trialgrps(itype).type=type{itype};
            trialgrps(tt).cat=[ttypes{tt}];
    end
end
end
if ~isempty(condtypes)
for itype=1:length(type)
    trialnums=[];
    trialinfo=trgrps(sessid).trialinfo;
    %trialgrps(1).trialgrps(itype).trials=[];
    
    %first get good trials for current trial type big/small/targ
    goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
    goodtrials=trialinfo(itype).nums{goodtrialgrp};  
    %get trial list    
    trtype=find(contains(infonames,type{itype}));
    trialtypes=trialinfo(trtype); 
    for tt=1:length(condtypes)       
        trialgrps(tt).trialgrps(itype).trials=[];       %initialize empty array        
        ttype=find(contains(trialtypes.names,condtypes{tt})==1);
        typename=trialtypes.names{ttype};
        seltrials=trialtypes.nums{ttype};
        trialnumtemp=seltrials;
        trialgrps(tt).site=siteda;
        trialgrps(tt).trialgrps(itype).trials=sort(unique([trialgrps(tt).trialgrps(itype).trials trialnumtemp]));
        trialgrps(tt).trialgrps(itype).type=type{itype};
        trialgrps(tt).cat=[condtypes{tt}];
    end
end
end
if ~isempty(trtypes)
for itype=1:length(type)
    trialnums=[];
    trialinfo=trgrps(sessid).trialinfo;
    %trialgrps(1).trialgrps(itype).trials=[];
    
    %first get good trials for current trial type big/small/targ
    goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
    goodtrials=trialinfo(itype).nums{goodtrialgrp};  
    %get trial list    
    trtype=find(contains(infonames,type{itype}));
    trialtypes=trialinfo(trtype); 
    trialgrps(itype).trialgrps(1).trials=[];       %initialize empty array 
     trialgrps(itype).cat=[type{itype} ];
    for tt=1:length(trtypes)       
        ttype=find(contains(trialtypes.names,trtypes{tt})==1);
        typename=trialtypes.names{ttype};
        seltrials=trialtypes.nums{ttype};
        trialnumtemp=seltrials;
        trialgrps(itype).site=siteda;
        trialgrps(itype).trialgrps(tt).trials=sort(unique([trialgrps(itype).trialgrps(tt).trials trialnumtemp]));
        trialgrps(itype).trialgrps(tt).type=type{itype};
             trialgrps(itype).cat=[trialgrps(itype).cat '_'  typename];

    end
           

end
end
if ~isempty(targvals)
    %get trialnums according to category of values for each type
    %'big'/sm/targ/break of trials, retaining origianl trial nums
    grpstemp=trialgrps;
    totalnums=length(allvals);
    [maxdasort,sortid]=sort(allvals);
    targtypes=[];
    trialids=[];
    for itype=1:length(type)
        if itype==1
        targtypes(1:length(grpstemp(1).trialgrps(itype).trials))=itype;
        trialids=1:length(grpstemp(1).trialgrps(itype).trials);
        else
        targtypes(length(targtypes)+1:length(targtypes)+length(grpstemp(1).trialgrps(itype).trials))=itype;
        trialids=[trialids 1:length(grpstemp(1).trialgrps(itype).trials)];
        end
    end
    sorttrtypes=targtypes(sortid);      %get trial type of sorted nums to differentiate trial group nums
    nums=round(pertile/100*totalnums);
    %top percentile
    pertrls{1}=sortid(end-nums+1:end);
    pertypes{1}=sorttrtypes(end-nums+1:end);
    pertrls{2}=sortid(1:nums);
    pertypes{2}=sorttrtypes(1:nums);
    trialgrps=[];
    cats={'top','bot'};
    for ix=1:length(pertrls)
        for itype=1:length(type)
            trialidstype=pertrls{ix}(pertypes{ix}==itype);
            trialgrps(ix).site=siteda;
            trialgrps(ix).cat=[cats{ix}  num2str(pertile)  pertype '_' event];
            trialgrps(ix).trialgrps(itype).trials=grpstemp(1).trialgrps(itype).trials(trialids(trialidstype));
            trialgrps(ix).trialgrps(itype).type=grpstemp(1).trialgrps(itype).type;
           % trialgrps{ix}(itype).trials=grpstemp(itype).trials(trialids(trialidstype));
            %trialgrps{ix}(itype).type=grpstemp(itype).type;
        end
    end       
end