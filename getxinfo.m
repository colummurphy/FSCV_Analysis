function xdata=getxinfo(sessnum,targda,xinfo,plotparam,varargin)
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'summary_da_');
sessnums=plotparam.sessnums;
daregions={'c','p'};
dalabels={'cn','put'};
conditions={'left','right'};
condlabels={'L','R'};
ratefscv=10;
plotvar={};
event='targ';
alldasites=plotparam.dasites;

normda=0;
meanda=0;
grpreg=0;
scattypes=0;
plotlfp=0;
datype='goodtrials';        %all good trials for xinfo data
xflag=0;            %get xinfo data
sessiontypes={};
xcflag=0;
subj='patra';
gain=0;
offset=[0 -0.2];
xwin=[0 4];
ratefscv=10;
plotvar={};
event='targ';
ttypes={};
argnum=1;
trtable={};
cellarray=0;
smoothdata=0;
smoothlength=2;         %200 ms
while argnum<=length(varargin)
    switch varargin{argnum}       
        case 'xwin'
            argnum=argnum+1;
            xwin=varargin{argnum};       
        case 'event'
            argnum=argnum+1;    
            event=varargin{argnum};         %event period, eg. fix, targ, outcome
        case 'metric'
            argnum=argnum+1;
            targval=varargin{argnum};       %datm variable, eg. targpeak, targwin, etc.
        case 'ttypes'
            argnum=argnum+1;
            ttypes=varargin{argnum};        %condition types grouped eg. {{'big','left'},{'small','left','aftersm'}}
        case 'trtable'
            %supply trial list table with labels of which groups to plot
            %(as generaetd by gettrtable & maketrorg)
            argnum=argnum+1;
            trtable=varargin{argnum};
            argnum=argnum+1;
            ttypes=varargin{argnum}; %eg. {{'big','-3break'},{'big','-2break'},{'big','-1break'}};
        case 'cleo'
            subj='cleo';            
        case 'cellarray'
            %cell array format, each trial separate cell
            cellarray=1;
        case 'smooth'
            smoothdata=1;
    end
    argnum=argnum+1;
end

sites=getsites(sessnum,{targda},subj);

[dapair,lfppair]=getsitepairs(targda,subj);
lfpsite=getlfpsites(sessnum,lfppair(1));

yvals={};
lfpvals={};
davals={};

targses=find(strcmp({trialgrps.sessid},num2str(sessnum)));
trialinfo=trialgrps(targses).trialinfo;
dasite=sites.site;
daregion=contains(targda,'c');     %1=='c'   

trialnumsall={};            %for each ttype
for icond=1:length(ttypes)
%get trials corresponding to groups of trial conditions in current group ttypes{icond}
%get trial types (big/small/break) for current group
targtrialtypes=[any(strcmp('big',ttypes{icond}))  ...
    any(strcmp('small',ttypes{icond})) ...
    any(strcmp('targetbreak',ttypes{icond})) ...
    any(strcmp('fixbreak',ttypes{icond})) ]; %get logical array for if want big/small/targ/fix
targtrialtypes=find(targtrialtypes==1);
nottypes=find(~strcmp(ttypes{icond},'big') & ...
    ~strcmp(ttypes{icond},'small') &...
    ~strcmp(ttypes{icond},'targetbreak') &...
   ~strcmp(ttypes{icond},'fixbreak'));
if isempty(sessiontypes)
    sessiontypes={'big','small'};
end
%get trial #'s for current condition, ie. not trial type identifier
%get da values for current group condition
%initialize each sess, as different avg  even if same
lfpvals{icond}=[];
davals{icond}=[];
cc=1;
counttype=1;
for itype=targtrialtypes
    trialnums=[];
    %get xinfo data
    targrow=find((contains({xinfo.siteda},targda) & ...
        contains({xinfo.sitelfp},lfppair(1)) & ...
        strcmp({xinfo.event},event)) & ...
        strcmp({xinfo.sessionid},num2str(sessnum)) &...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);        
    xtarg=getfield(xinfo(targrow),'daall');    %'dapos' or 'daneg' types get below   
for tt=nottypes
    if isempty(trtable)
        %get list from trialinfo that has coarse groups defined
        curcond=ttypes{icond}{tt};
        targt=find(contains(trialinfo(itype).names,curcond)==1);
        trialnums=[trialnums trialinfo(itype).nums{targt}];  
    else
        %first get good trials for current trial type big/small/targ
        goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
        goodtrials=trialinfo(itype).nums{goodtrialgrp};  
        %get trial list
        listses=find(strcmp({trtable.sessid},num2str(sessnum)));
        listtyp=find(contains({trtable(listses).trorg.type},sessiontypes(itype)));
        trlisttarg=find(contains({trtable(listses).trorg(listtyp).grp.label},ttypes{icond}{tt}));
        trialnumtemp=trtable(listses).trorg(listtyp).grp(trlisttarg).trials;
        trialnumtemp=intersect(trialnumtemp,goodtrials);
        trialnums=[trialnums trialnumtemp];
    end
end
trialnums=unique(trialnums);        %trial nums for current trial type big/small
if ~isempty(xtarg)
    trialids=[];
    trialnums=intersect(trialnums,getfield(xinfo(targrow),datype));     %get dapos/daall/neg trials
    trialids=find(ismember(xtarg.trials,trialnums)==1); %trial ids for targeted metric
    alnidx=getfield(xtarg,'mididx');  
    xwinids=xwin(1)*ratefscv+alnidx+offset(1)*ratefscv:xwin(2)*ratefscv+alnidx+offset(2)*ratefscv;        
    tempdata=[];
    lfpdata=getfield(xtarg,'lfptracesaln');
    dadata=getfield(xtarg,'datracesaln');     
    if ~cellarray
        if itype==1
            lfpvals{icond}=lfpdata(trialids,xwinids);
            davals{icond}=dadata(trialids,xwinids);
        else
            lfpvals{icond}=[lfpvals{icond}; lfpdata(trialids,xwinids)];         %da values for specific group of conditions (ie. trial nums) 
            davals{icond}=[davals{icond}; dadata(trialids,xwinids)];         %da values for specific group of conditions (ie. trial nums) 
        end
    else
        ltemp=lfpdata(trialids,xwinids);
        dtemp=dadata(trialids,xwinids);
        ltemp=ltemp';
        dtemp=dtemp';
        for ix=1:size(ltemp,2)
            if smoothdata
                ltempsmooth=smoothwin(ltemp(:,ix)',smoothlength);   %smoothing 
                ltemp(:,ix)=ltempsmooth';
            end
            ltdata=ltemp(:,ix);
            dtdata=dtemp(:,ix);
            if any(isnan(dtdata))
                dtdata2=dtdata(~isnan(dtdata));
                ltdata2=ltdata(~isnan(dtdata));
                ltdata=ltdata2;
                dtdata=dtdata2;
            end
            if any(isnan(ltdata))
                dtdata2=dtdata(~isnan(ltdata));
                ltdata2=ltdata(~isnan(ltdata));
                ltdata=ltdata2;
                dtdata=dtdata2;
            end
            if any(isnan(ltdata)) || any(isnan(dtdata)) || any(isinf(dtdata)) || any(isinf(ltdata))
                warning('');
                pause(10);
            end
            
            lfpvals{icond}{cc}=ltdata;
            davals{icond}{cc}=dtdata;
            cc=cc+1;
        end
    end
            
        
end
trialnumsall{counttype}=trialnums;
counttype=counttype+1;
end
end
%finish getting condition types


xdata.lfp=lfpvals;
xdata.da=davals;
xdata.offset=offset;
xdata.win=xwin;
xdata.event=event;
xdata.ttypes=ttypes;
xdata.trialnums=trialnumsall;

end
