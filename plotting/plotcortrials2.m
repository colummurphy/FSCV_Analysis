function dataout=plotcortrials2(sessnum,trlists,datms,lfptms,binfos,plotparam,varargin)
%call from corsmultiple
%4/2020, output relevant parameters for tabulating
%supply xinfo for manual lfp win's and trial #'s comparison
%noticed differences in R, between lfptm targimwin & 0-1 win xinfo maybe
%because more trials in lfptm than xinfo

sessiontypes={'bigreward','smallreward'};
condtype='all';
sessid=num2str(sessnum);
event='targ';
trlistid=find(strcmp({trlists.sessid},sessid));
trlist=trlists(trlistid).list;
numtrials=length(trlist);
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'cor_scatter',filesep,['trialcorr_' sessid]);
if ~isdir(savepath)
mkdir(savepath);
end
sessnums=plotparam.sessnums;
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites);
uniquesites=unique({sites(1:end).site});
[dapair,lfppair]=getsitepairs(targdasites);
targses=find(strcmp({trialgrps.sessid},sessid));
trialinfo=trialgrps(targses).trialinfo;

%initialize trlist
for ix=1:length(trlist)
    trlist(ix).da=nan;
    trlist(ix).lfp=nan;   
     trlist(ix).hr=nan;   
          trlist(ix).eye=nan; 
          trlist(ix).eyev=nan;   
           trlist(ix).eyedist=nan;  
     trlist(ix).lick=nan;   
     trlist(ix).rr=nan;   
     trlist(ix).rrstd=nan;   
     trlist(ix).rmssd=nan;   
end

figpos=[50,50,1400,900];
xoff=100;
yoff=100;
mar=25;
axsize=[400,250];
fontsize=15;
axpos={};
axa={};
triallim=[];
xdata=1:numtrials;
bigtype=[];
smalltype=[];
targtype=[];
fixtype=[];
marksize=[50,50,20,30];
mark={'.','.','o','x'};
plotvar={};
argnum=1;
siteda='';
sitelfp='';
bflag=0;
plotrt=0;
plothr=0;
ploteye=0;
linfos={};
xinfo={};
xbinfos={};
hfsametype=0;
metric='targpeak';
metriclfp='targwin';
win = [0 4];
offset=[0 0];
winb= [.1 4];
rate=10;            %sample rate for xinfo data
%metriclfp='targpeakabs';
plotnum=0;
plothist=[];
plottypes=[];
wine=[];
winhr=[];
winl=[];
dataout={};
relfp=0;
noplot=0;
nodata=0;
xvarflag=0;
daonly=0;
lfponly=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dax'
            %values from xinfo
            argnum=argnum+1;
            siteda=varargin{argnum};
            argnum=argnum+1;
            xinfo=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'datm'
            %values from datm, supplied instead of xinfo
            argnum=argnum+1;
            siteda=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'lfpx'
            argnum=argnum+1;
            sitelfp=varargin{argnum};            
            argnum=argnum+1;
            xinfo=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='lfp';
            relfp=1;
        case 'lfptm'
            %values from betatm
            argnum=argnum+1;
            sitelfp=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='lfp';
        case 'behx'
            %values from xbinfos
            argnum=argnum+1;
            xbinfos=varargin{argnum};
        case 'metric'
            %for datm values
            argnum=argnum+1;
            metric=varargin{argnum};
        case 'metriclfp'
            %for betatm values
            argnum=argnum+1;
            metriclfp=varargin{argnum};
        case 'trialseq'
            plotnum=plotnum+1;
            plotvars{plotnum}='trialnum';
        case 'rt'
            plotrt=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='trt';
        case 'frt'
            plotrt=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='frt';
        case 'hr'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='hr';
        case 'rr'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rr';
        case 'rrstd'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rrstd';
        case 'rmssd'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rmssd';
        case 'lick'
            plotnum=plotnum+1;
            plotvars{plotnum}='lick';
        case 'eye'
            plotnum=plotnum+1;
            plotvars{plotnum}='eye';
        case 'eyev'
            %max eye velocity within ewin (1 s)            
            plotnum=plotnum+1;
            plotvars{plotnum}='eyev';
        case 'eyedist'
            %max of integral of abs eye velocity within ewin (1 s)            
            plotnum=plotnum+1;
            plotvars{plotnum}='eyedist';
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};    
        case 'hf'
            argnum=argnum+1;
            plothist=varargin{argnum};
        case 'hfsametype'
            hfsametype=1;       %constrict types of futre/past trials
        case 'types'
            argnum=argnum+1;
            plottypes=varargin{argnum};
        case 'condition'
            argnum=argnum+1;
            condtype=varargin{argnum};
        case 'bwin'
            argnum=argnum+1;
            winb=varargin{argnum};
        case 'win'
            argnum=argnum+1;
            win=varargin{argnum};
            metriclfp=['win' num2str(win(1)) '_' num2str(win(2))];
        case 'offset'
            argnum=argnum+1;
            offset=varargin{argnum};
            metriclfp=['win' num2str(win(1)+offset(1)) '_' num2str(win(2)+offset(2))];
        case 'lwin'
            argnum=argnum+1;
            winl=varargin{argnum};
        case 'hrwin'
            argnum=argnum+1;
            winhr=varargin{argnum};
        case 'ewin'
            argnum=argnum+1;
            wine=varargin{argnum};
        case 'noplot'
            noplot=1;
        case 'xinfo'
            %supply xinfo for trial #'s comparison
            argnum=argnum+1;
            xinfo=varargin{argnum};
        case 'nodata'
            %don't save data in struct just # trials
            nodata=1;
        case 'xvar'
            %use metriclfp variable to load variable in xinfos, ignore win
            %use with lfpx argument
            xvarflag=1;
            win=[0 0];
            offset=[0 0];
        case 'daonly'
            %only da and beh
            daonly=1;
        case 'lfponly'
            %only lfp and beh
            lfponly=1;
    end
    argnum=argnum+1;
end   
hf={};
if ispc && ~noplot
hf=figure('visible','off');     %figure for each channe
set(hf,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',hf);    %set figure handle to current figure
set(hf,'visible','on');     %figure for each channel
else
   % set(hf,'visible','off');
end

if any(win~=0) && relfp && ~xvarflag
    metriclfp=['win' num2str(win(1)+offset(1)) '_' num2str(win(2)+offset(2))];
end
%set up plots
if isempty(plothist)
plothist=zeros(1,length(plotvars));       %1 means 1 trial in future value
end
if isempty(plottypes)
plottypes=repmat({'all'},1,length(plotvars));
end

numplots=sum([length(plotvars)-1:-1:1]);        % # possible unique combinations correlations, for 3 variables, this is da vs beta, da vs trial, beta vs trials, num vars being correlated to single var + remaining vars after first iteration
labels=[];
if ~noplot
clf(hf,'reset');
set(hf,'color',[1 1 1]);
axpos={};
for ip=1:numplots
    axa{ip}=subplot(ceil(numplots/2),2,ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
        axpos{ip}=get(axa{ip},'position');    

    %{
    set(axa{ip},'position',[xoff,figpos(4)-axsize(2)*ip-yoff-mar*(ip-1),axsize(1),axsize(2)])
    if ip>1
        set(axa{ip},'position',[xoff,axpos{ip-1}(2)-axsize(2)-mar,axsize(1),axsize(2)])
    end
    %}
    %ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
end
end
titletext=['sess ' sessid ' | ' event ' | '  plottypes{1} ' | ' condtype ' | ' siteda ' ' metric ' | ' sitelfp ' ' metriclfp ];  
for ix=1:length(plotvars)
    if ix<length(plotvars)
    labels=[labels plotvars{ix} '_'];
    else
    labels=[labels plotvars{ix}];
    end        
end
labels=[labels '_' event '_' plottypes{1} '_' condtype ];
if ~isempty(siteda)
    labels=[labels '_' metric '_' siteda];
end
if ~isempty(sitelfp)
    if ~relfp
    labels=[labels '_' metriclfp '_' sitelfp];
    elseif relfp && any(win~=0) && ~xvarflag
        labels=[labels '_win_' num2str((win(1)+offset(1))*rate) '_' num2str((win(2)+offset(2))*rate) '_' sitelfp];
    elseif relfp && xvarflag
        labels=[labels '_' metriclfp '_' sitelfp];
    end  
end
if ~isempty(winhr)
    labels=[labels '_winhr_' num2str(winhr(1)*rate) '_' num2str(winhr(2)*rate)];
        titletext=[titletext ' | ' 'winhr' num2str(winhr(1)) '-' num2str(winhr(2)) 's' ];  
end
if ~isempty(wine)
    labels=[labels '_wine_' num2str(wine(1)*rate) '_' num2str(wine(2)*rate)];
    titletext=[titletext ' | ' 'wine' num2str(wine(1)) '-' num2str(wine(2)) 's' ];  
end
if ~isempty(winl)
    labels=[labels '_winl_' num2str(winl(1)*rate) '_' num2str(winl(2)*rate)];
    titletext=[titletext ' | ' 'winl' num2str(winl(1)) '-' num2str(winl(2)) 's' ];  
end
if isempty(wine)
    wine=winb;
end
if isempty(winhr)
    winhr=winb;
end
if isempty(winl)
    winl=winb;
end
%compile data
%first get good trials list for trial type big/small/targ
maplist={};     %listing of rows mapped by sel trials; ie. seltrials{1}(1) --> maplist{1}(1) = targ row for specific type & trial
seltrials={};
dtm={};
lfptm=[];
xtargda=[];
xtarglfp=[];
goodtrials=[];
labx='';
laby='';
for itype=1:length(sessiontypes)
    if strcmp(condtype,'all')
    goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
    incswgrp=find(strcmp(trialinfo(itype).names,'switch')==1);                                                  %includes switch trials.. maybe need to remove to match plotdasummary2
    goodtrials=sort(unique([trialinfo(itype).nums{goodtrialgrp} trialinfo(itype).nums{incswgrp}]));  
    else
    goodtrialgrp=find(contains(trialinfo(itype).names,condtype)==1);
    goodtrials=sort(unique([trialinfo(itype).nums{goodtrialgrp}]));          
    end
    seltrials{itype}=goodtrials;
    tlabelsgood=goodtrials+99;              %convert to trial labels ie 100,101,etc.
    typelistids=find(strcmp({trlist.type},sessiontypes{itype}));
    maplist(itype).trlist=find(ismember([trlist.id],tlabelsgood) & ...
        strcmp({trlist.type},sessiontypes{itype}));
    maplist(itype).trialnums=goodtrials;
    maplist(itype).type=sessiontypes{itype};    
    
%get hr  from binfos
hr=[];
rr=[];
rrstd=[];
rmssd=[];
hrtrials=[];
if isempty(xbinfos)
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    hr=binfos(btarg).pulse;
    hrtrials=binfos(btarg).seltrials;
else 
    %get hr from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'pulse') &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));            %event matters since where calculated
    if ~isempty(targrows)
    temp=xbinfos(targrows(1)).daall.lfptracesaln;
    hrtrials=xbinfos(targrows(1)).daall.trials;
        alnidx=xbinfos(targrows(1)).daall.mididx;
            hhwin=alnidx+winhr(1)*rate:alnidx+winhr(2)*rate;
    hrwin=temp(:,hhwin); 
    hr=nanmean(hrwin,2);
   rrtemp=1./temp.*60; 
   rrwin=rrtemp(:,hhwin);
   rr=nanmean(rrwin,2);
    rrstd=nanstd(rrwin,[],2);
    %outofrange=find(rrstd<.01);
   % rrstd(outofrange)=nan;
    rmssd=rms(diff(rrtemp(:,hhwin),1,2),2);    
    end
    
end
for ida=1:length(hrtrials)
    tid=find(ismember(maplist(itype).trialnums,hrtrials(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).hr=hr(ida);
        if ~isempty(xbinfos)
        trlist(maplist(itype).trlist(tid)).rr=rr(ida);
        trlist(maplist(itype).trlist(tid)).rrstd=rrstd(ida);
        trlist(maplist(itype).trlist(tid)).rmssd=rmssd(ida);
        end
    end
end

%get lick  from binfos
lick=[];
licktrials=[];
if isempty(xbinfos)
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    lick=binfos(btarg).lickpost;
    licktrials=binfos(btarg).seltrials;  
else
    %get  from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'lick') &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));            %event matters since where calculated
        if ~isempty(targrows)

    temp=xbinfos(targrows(1)).daall.lfptracesaln;
    licktrials=xbinfos(targrows(1)).daall.trials;
            alnidx=xbinfos(targrows(1)).daall.mididx;

    licks=temp(:,alnidx+winl(1)*rate:alnidx+winl(2)*rate); 
    lick=nanmean(licks,2);
        end
end
for ida=1:length(hrtrials)
    tid=find(ismember(maplist(itype).trialnums,licktrials(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).lick=lick(ida);
    end
end 
    
%eye diameter, only right side valid
eye=[];
etrials=[];
eyevmax=[];
eyedist=[];
etrialsall=[];
btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
    strcmp({binfos.sessionid},sessid) &...
    strcmp({binfos.evt},event)); 
if isempty(xbinfos)    
           %event matters since where calculated
    eye=-binfos(btarg).reyed;
    etrials=binfos(btarg).seltrialsr;    
else
    %get  from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'eye') &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));   
    %pupil diameter, only right side, since opened door on left
            %many times   
        if ~isempty(targrows)
            etrialsall=xbinfos(targrows(1)).daall.trials;                     
            etrialsr=binfos(btarg).seltrialsr;      %right side only
            etrialsel=find(ismember(etrialsall,etrialsr));
            etrials=intersect(etrialsall,etrialsr);
            temp=xbinfos(targrows(1)).daall.lfptracesaln(etrialsel,:);
            alnidx=xbinfos(targrows(1)).daall.mididx;
            eewin=[alnidx+wine(1)*rate:alnidx+wine(2)*rate];
            eyewin=-temp(:,eewin); 
            eye=nanmean(eyewin,2);            
            %Z SCORE EYE 072019 from plottracesel
            glitchwidth=2;
            basewin=[alnidx-7*rate:alnidx+4*rate];          %fix baseline window for eye below
            %zscore eye wrt baseline as done in olivier et al 2014 frontier
            %remove blinks
            extwin=[alnidx-3*rate:alnidx+4*rate];
            relwin=intersect(eewin,extwin)-extwin(1)+1;
            datawin=-temp(:,extwin); %flip so original sign
            datawin=deglitchinterpsimple(datawin,3.5,glitchwidth);
            basedata=-temp(:,basewin);   %normalizing to pre-fix and targ window allows normalization trial to trial to look at faster changes within 1 sec 06/04 without phase trend eg. 67/68/69
            basedata=-basedata;
            for ix=1:size(basedata,1)
                [xx, xb]=rmoutliers(basedata(ix,:));
                basedata(ix,xb)=nan;
            end
            baseline=nanstd(basedata,[],2);          %normalized to each trial so independent of phase changes
            dataz=(datawin-nanmean(basedata,2))./baseline;
            datawin=dataz;
            eyewin=dataz(:,relwin);     %user selected windo
            eye=nanmean(eyewin,2);
        end
        %max velocity within 1 s of targ
        %eye distance integral of velocity
        targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
            strcmp({xbinfos.sitelfp},'eyev') &...
            strcmp({xbinfos.sessionid},sessid) &...
            strcmp({xbinfos.event},event)); 
       if ~isempty(targrows)
            etrialsall=xbinfos(targrows(1)).daall.trials;   
            temp=xbinfos(targrows(1)).daall.lfptracesaln;
            alnidx=xbinfos(targrows(1)).daall.mididx;
            eewin=[alnidx+wine(1)*rate:alnidx+wine(2)*rate];
            eyewin=-temp(:,eewin); 
            eyevmax=max(abs(eyewin),[],2);    
            eyevcum=cumtrapz(abs(eyewin),2);
            eyedist=max(eyevcum,[],2);
        end
end
for ida=1:length(etrials)
    tid=find(ismember(maplist(itype).trialnums,etrials(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).eye=eye(ida);
    end
end
for ida=1:length(etrialsall)
    tid=find(ismember(maplist(itype).trialnums,etrialsall(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).eyev=eyevmax(ida);
        trlist(maplist(itype).trlist(tid)).eyedist=eyedist(ida);
    end
end

%get signal values 
    dtm{itype}=[];    
    %datm/betatm values
    %get ch # for corresponding session, based on site label
    ich=0;
    lch=0;
    if ~isempty(siteda)
        %targeted da values from datm
    targsi=find([sites.sessnum]==sessnum & ...
    strcmp({sites.probeid},siteda));
    ich=sites(targsi(1)).ch;        %ch # defined by fscv recording ch
    dtm{itype}=datms{targses}{itype}{ich};           %get all da targeted values for current trial type big/small
    targtrials=dtm{itype}.trialnums;
            if ~isempty(xinfo)
                %compare trial #'s, get trial #'s from xinfo
                 targrows=find((strcmp({xinfo.sitelfp},sitelfp) & ...
                        strcmp({xinfo.siteda},siteda) & ...
                        strcmp({xinfo.event},event)) & ...
                        strcmp({xinfo.sessionid},sessid) &...
                        contains({xinfo.sessiontype},sessiontypes{itype})==1);   
                 xtrials=xinfo(targrows(1)).daall.trials;
                 targtrials=intersect(dtm{itype}.trialnums,xtrials);
            end
   % targtrials=find(ismember(dtm{itype}.trialnums,maplist(itype).trialnums));
    datatemp=getfield(dtm{itype},metric);
    targtrialids=find(ismember(dtm{itype}.trialnums,targtrials));
    dadata=datatemp(targtrialids);  
            for ida=1:length(targtrials)
                tid=find(ismember(maplist(itype).trialnums,targtrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).da=dadata(ida);
                end
            end
    end

    if ~isempty(sitelfp) 
        %targeted beta values from betatm
        if ~relfp
            %get precalc values in lfptms
            lfptm=lfptms{targses}{itype};
            for ii=1:length(lfptm)
                if strcmp(lfptm{ii}.site,sitelfp)
                    lch=ii;
                end
            end
            ltm{itype}=lfptm{lch};           %get all da targeted values for current trial type big/small
            targtrials=ltm{itype}.trialnums;
            if ~isempty(xinfo)
                %compare trial #'s, get trial #'s from xinfo
                 targrows=find((strcmp({xinfo.sitelfp},sitelfp) & ...
                        strcmp({xinfo.siteda},siteda) & ...
                        strcmp({xinfo.event},event)) & ...
                        strcmp({xinfo.sessionid},sessid) &...
                        contains({xinfo.sessiontype},sessiontypes{itype})==1);   
                 xtrials=xinfo(targrows(1)).daall.trials;
                 targtrials=intersect(ltm{itype}.trialnums,xtrials);
            end
          %  targtrials=find(ismember(ltm{itype}.trialnums,maplist(itype).trialnums));
            ltemp=getfield(ltm{itype},metriclfp);
                targtrialids=find(ismember(ltm{itype}.trialnums,targtrials));
            ldata=ltemp(targtrialids);
            for ida=1:length(targtrials)
                tid=find(ismember(maplist(itype).trialnums,targtrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).lfp=ldata(ida);
                end
            end
        elseif relfp && ~isempty(xinfo)
            %else recalculate from xinfos and defined win
            emptyflag=0;
            targrows=find((strcmp({xinfo.sitelfp},sitelfp) & ...
                contains({xinfo.siteda},siteda) & ...
            strcmp({xinfo.event},event)) & ...
            strcmp({xinfo.sessionid},sessid) &...
            contains({xinfo.sessiontype},sessiontypes{itype})==1);   
            if isempty(targrows)
                emptyflag=1;
                disp(['ch ' sitelfp ' not found in xinfo']);
            end
            if length(targrows)>1
                disp(['ch ' sitelfp ' & event ' event ' more than one instance in xinfo']);
                if isempty(siteda)
                    %choose lfp channel with more good trials, don't care
                    %about da pairing, only lfp signal used
                    trialsintargrows={xinfo(targrows).goodtrials};
                    maxid=targrows(1);
                    for itt=2:length(trialsintargrows)
                        if length(trialsintargrows{itt})>length(trialsintargrows{itt-1})
                            maxid=targrows(itt);
                        end
                    end
                    targrows=maxid;
                else                    
                    dbstop in plotcortrials2 at 563 if length(targrows)>1
                end
            end
            xtrials=xinfo(targrows(1)).daall.trials;
            alnidx=xinfo(targrows(1)).daall.mididx;
            lfpvals=[];
         %   targtrials=find(ismember(xtrials,maplist(itype).trialnums));
            if any(win~=0)
                %recalculate average based on user supplied window
            xdata=xinfo(targrows(1)).daall.lfptracesaln;
           % lfpwinids=[alnidx+win(1)*rate+offset(1)*rate:alnidx+win(2)*rate+offset(2)*rate];
           % %windowing used previously < 4/2020 for manual [0 1] wins
                     %   lfpwinids=[alnidx+win(1)*rate+1:alnidx+win(2)*rate+1];      %windowing used in targimwin generated in setTrialAx in plotsession
            lfpwinids=[alnidx+win(1)*rate+offset(1)*rate:alnidx+win(2)*rate+offset(2)*rate];
            lfpvals=nanmean(xdata(:,lfpwinids),2)';
            end
            if xvarflag
                %open xinfos variable as defined by metriclfp
                lfpvals=getfield(xinfo(targrows(1)).daall,metriclfp);
            end
            for ida=1:length(xtrials)
                tid=find(ismember(maplist(itype).trialnums,xtrials(ida)));
                if ~isempty(tid)
                    trlist(maplist(itype).trlist(tid)).lfp=lfpvals(ida);
                end
            end
        end
    end   
    
end
bigtype=double(strcmp({trlist.type},'bigreward'));
bigtype(bigtype==0)=nan;
smalltype=double(strcmp({trlist.type},'smallreward'));
smalltype(smalltype==0)=nan;
targtype=double(strcmp({trlist.type},'targetbreak'));
targtype(targtype==0)=nan;
fixtype=double(strcmp({trlist.type},'fixbreak'));
fixtype(fixtype==0)=nan;
sidetype=[trlist.side];
%find switch trials
smallside=smalltype+sidetype;
bigside=bigtype+sidetype;
smallside2=fillmissing(smallside,'previous');
bigside2=fillmissing(bigside,'previous');
rewsides=smallside2-bigside2;
swtrials=rewsides==0;

iax=1;
nonnan=[];
for ip=1:length(plotvars)
    curx=plotvars{ip};      %for each plotvars assign as x variable and look at all other plotvars as y in relation below
    remplots=numplots-ip;       %subtract remaining plots after each cycle of unique correlations
    for iip=ip+1:length(plotvars)
        if ~noplot
        cla(axa{iax});
        end
        datax=[];
        datay=[];
        cury=plotvars{iip};     %next other variable in plotvars assigned to y
        if strcmp(curx,'trialnum')
            labx='trial #';
            datax=1:length(trlist);
        else
            datax=[trlist.(curx)];
        end
        if strcmp(cury,'trialnum')
            laby='trial #';
            datay=1:length(trlist);
        else
            datay=[trlist.(cury)];
        end
        switch curx
            case 'da'
                labx=['\DeltaDA (nM) ' metric ' ' siteda];
            case 'lfp'
                labx=['\beta-LFP (\muV^2) ' metriclfp ' ' sitelfp];
            case 'hr'
                labx=['hr (bpm)'];
            case 'trt'
                labx=['target reaction time (s)'];
            case 'frt'
                labx=['center reaction time (s)'];
            case 'lick'
                labx='lick';
                datax=(datax-nanmean(datax))./nanstd(datax);
            case 'eye'
                labx='eye';
                normval=abs(nanmedian(datax));
                datax=datax./normval;
           case 'eyedist'
                labx='eyedist';
                normval=abs(nanmedian(datax));
                datax=datax./normval;
           case 'eyev'
                labx='eyev';
                normval=abs(nanmedian(datax));
                datax=datax./normval;
            case 'rrstd'
                labx='RRstd (ms)';
                datax=datax.*1e3;
            otherwise
                labx=curx;
        end
        switch cury
            case 'da'
                laby=['\DeltaDA (nM) ' metric ' ' siteda];
            case 'lfp'
                laby=['\beta-LFP (\muV^2) ' metriclfp ' ' sitelfp];
            case 'hr'
                laby=['hr (bpm)'];
            case 'trt'
                laby=['target reaction time (s)'];
            case 'frt'
                laby=['center reaction time (s)'];
            case 'lick'
                laby='lick';
                datay=(datay-nanmean(datay))./nanstd(datay);
            case 'eye'
                laby='eye';
                normval=abs(nanmedian(datay));
                datay=datay./normval;
           case 'eyedist'
                laby='eyedist';
                normval=abs(nanmedian(datax));
                datax=datax./normval;
           case 'eyev'
                laby='eyev';
                normval=abs(nanmedian(datax));
                datax=datax./normval;
           case 'rrstd'
                laby='RRstd (ms)';
                datay=datay.*1e3;
            otherwise
                laby=cury;
        end
        if plothist(iip)<0
            %future
            labx=[curx ' next trial'];
        elseif plothist(iip)>0
            %past
            labx=[curx ' past trial'];
        end
        if plothist(ip)<0
            %future
            laby=[cury ' next trial'];
        elseif plothist(ip)>0
            %past
            laby=[cury ' past trial'];
        end
       % datax=datax([maplist.trlist]);
        datax=circshift(datax,[0 plothist(iip)]);       %shift -1 for future trial, +1 for past trial
        
       % datay=datay([maplist.trlist]);
        datay=circshift(datay,[0 plothist(ip)]);
         databigy=bigtype;
         databigx=bigtype;
         datasmally=smalltype;
         datasmallx=smalltype;
         if hfsametype
             %constrict to same type for future/past trial
             if plothist(iip)~=0
                %future/past trial type
                labx=[labx ' (same condition)'];
             end
             if plothist(ip)~=0
                %future/past
                laby=[laby ' (same condition)'];
             end
             databigy=circshift(bigtype,[0 plothist(ip)]);
             databigx=circshift(bigtype,[0 plothist(iip)]);
             datasmally=circshift(smalltype,[0 plothist(ip)]);
             datasmallx=circshift(smalltype,[0 plothist(iip)]);
         end   
        xoutthres=nanmean(datax)+3*nanstd(datax);
        youtthres=nanmean(datay)+3*nanstd(datay);
        remx=[];        %remove outliers
        remy=[];
        if ~strcmp(curx,'eye')
        remx=find(abs(datax)>xoutthres | abs(datax)<nanmean(datax)-4*nanstd(datax));
        end
        if ~strcmp(cury,'eye')
        remy=find(abs(datay)>youtthres );
        end
        if strcmp(curx,'frt')
            remx2=find(datax<.05);
            remx=unique([remx remx2]);
        end
        if strcmp(cury,'frt')
            remy2=find(datay<.05);
            remy=unique([remy remy2]);
        end
        if ~isempty(remx)
            datax(remx)=nan;
        end
        if ~isempty(remy)
            datay(remy)=nan;
        end
       % nonnan=find(~isnan(datax) & ~isnan(datay));
        nonnan=find(~isnan(datax) & ~isnan(datay) & ~isoutlier(datax) & ~isoutlier(datay)); %based on median outlier removal rather than mean 04/25/2020
        if hfsametype
            nonnan=find(~isnan(datax) & ~isnan(datay) & (~isnan(datasmallx) | ~isnan(databigx)) & (~isnan(datasmally) | ~isnan(databigy)));
        end
        origdatax=datax;
        origdatay=datay;
        datax=datax(nonnan);
        datay=datay(nonnan);
        databig=bigtype(nonnan);
        datasmall=smalltype(nonnan);
        if strcmp(plottypes{iip},'all')
            if ~noplot
            scatter(axa{iax},databig.*datax,databig.*datay,20,'o','markeredgecolor',[1 0 0],'markerfacecolor',[1 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);
            scatter(axa{iax},datasmall.*datax,datasmall.*datay,20,'o','markeredgecolor',[0 0 0],'markerfacecolor',[0 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);
            end
        elseif strcmp(plottypes{iip},'big')
             datax=databig.*datax;
             datay=databig.*datay;
             xoutthres=nanmean(datax)+3*nanstd(datax);
            youtthres=nanmean(datay)+3*nanstd(datay);
            if ~strcmp(curx,'eye')
                remx=find(abs(datax)>xoutthres);
            end
            if ~strcmp(cury,'eye')
                remy=find(abs(datay)>youtthres);
            end
            if ~isempty(remx)
                datax(remx)=nan;
            end
            if ~isempty(remy)
                datay(remy)=nan;
            end
             %nonnan=find(~isnan(datax) & ~isnan(datay));
             nonnan=find(~isnan(datax) & ~isnan(datay) & ~isoutlier(datax) & ~isoutlier(datay)); %based on median outlier removal rather than mean 04/25/2020
             datax=datax(nonnan);
             datay=datay(nonnan);
             if ~noplot
             scatter(axa{iax},datax,datay,20,'o','markeredgecolor',[1 0 0],'markerfacecolor',[1 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);
             end
        elseif strcmp(plottypes{iip},'small') 
             datax=datasmall.*datax;
             datay=datasmall.*datay;
            xoutthres=nanmean(datax)+3*nanstd(datax);
            youtthres=nanmean(datay)+3*nanstd(datay);
            if ~strcmp(curx,'eye')
            remx=find(abs(datax)>xoutthres);
            end
            if ~strcmp(cury,'eye')
                                        remy=find(abs(datay)>youtthres);
            end
            if ~isempty(remx)
            datax(remx)=nan;
            end
            if ~isempty(remy)
            datay(remy)=nan;
            end
            % nonnan=find(~isnan(datax) & ~isnan(datay));
            nonnan=find(~isnan(datax) & ~isnan(datay) & ~isoutlier(datax) & ~isoutlier(datay)); %based on median outlier removal rather than mean 04/25/2020
             datax=datax(nonnan);
             datay=datay(nonnan);
             if ~noplot
                         scatter(axa{iax},datax,datay,20,'o','markeredgecolor',[0 0 0],'markerfacecolor',[0 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);  
             end
        end
        if ~noplot
        ylabel(axa{iax},laby);
        xlabel(axa{iax},labx);
        end
        [r,psig]=corr(datax',datay');
        [p,s]=polyfit(datax,datay,1);
        slope=p(1);
        intercep=p(2);
        yfit=slope*datax+intercep;
        yresid=datay-yfit;
        ssresid=sum(yresid.^2);
        sstotal=(length(datay)-1)*var(datay);
        rsq=1-ssresid/sstotal;
        if ~noplot
        fitline=plot(axa{iax},datax,yfit,'-','color',[0 0 0]);
        fitline.Color(4)=0.5;
        xsiz=get(axa{iax},'position');
        ylims=get(axa{iax},'ylim');
        text(axa{iax},xsiz(3)-120,xsiz(4)-50,...
            {['r: ' num2str(r)],...
            ['p: ' num2str(psig)]},'color',[0 0 0],'units','pixels','FontSize',fontsize)
                set(findall(axa{iax},'-property','FontSize'),'FontSize',fontsize)
        end
        if strcmp(curx,'da') && strcmp(cury,'lfp') && iax==1
            dataout(1).sess=sessnum;
            dataout(1).x=curx;
            dataout(1).sitechda=siteda;
            dataout(1).siteda=siteda;
            dataout(1).y=cury;
            dataout(1).sitechlfp=sitelfp;
            dataout(1).sitelfp=sitelfp;
            dataout(1).event=event;
            dataout(1).type=plottypes{1};
            dataout(1).metricda=metric;
            dataout(1).metriclfp=metriclfp;
            if ~nodata
                    %dataout(1).datax=datax;
                   % dataout(1).datay=datay;
                    dataout(1).datax=origdatax;     %keep original trial ids 05/20/2020 for subseq sorting of groups, bigtype. etc
                    dataout(1).datay=origdatay;
                                        dataout(1).dataxused=datax;     %keep original trial ids 05/20/2020 for subseq sorting of groups, bigtype. etc
                    dataout(1).datayused=datay;
                    dataout(1).nonnan=nonnan;
                else
                    dataout(1).trials=length(origdatax); %keep original trial ids 05/20/2020 for subseq sorting of groups
                    
            end
            if relfp && ~xvarflag
                dataout(1).metriclfp=win;
            end
            dataout(1).r=r;
            dataout(1).p=psig;
            bigtype(isnan(bigtype))=0;
            bigtype=logical(bigtype);
            dataout(1).big=bigtype;
            smalltype(isnan(smalltype))=0;
            smalltype=logical(smalltype);
            sidetype(isnan(sidetype))=0;
            sidetype=logical(sidetype);
            dataout(1).small=smalltype;
            dataout(1).side=sidetype;
            targtype(isnan(targtype))=0;
            targtype=logical(targtype);
            dataout(1).targbreak=targtype;
            fixtype(isnan(fixtype))=0;
            fixtype=logical(fixtype);
            dataout(1).fixbreak=fixtype;
            dataout(1).switch=swtrials;
        else            
            %other physiologic correlations
            %dataout(1).y2{iax}=cury;
           % dataout(1).datay2{iax}=datay;
           % dataout(1).r2{iax}=r;
           % dataout(1).p2{iax}=psig;
           if iax==1
                dataout(1).sess=sessnum;
                 dataout(1).x=curx;
                dataout(1).sitechda=siteda;
                dataout(1).siteda=siteda;
                dataout(1).y=cury;
                dataout(1).sitechlfp=sitelfp;
                dataout(1).sitelfp=sitelfp;
                dataout(1).event=event;
                dataout(1).type=plottypes{1};
                dataout(1).metricda=metric;
                dataout(1).metriclfp=metriclfp;
                if ~nodata
                    %dataout(1).datax=datax;
                   % dataout(1).datay=datay;
                    dataout(1).datax=origdatax;     %keep original trial ids 05/20/2020 for subseq sorting of groups, bigtype. etc
                    dataout(1).datay=origdatay;
                    dataout(1).dataxused=datax;     %keep original trial ids 05/20/2020 for subseq sorting of groups, bigtype. etc
                    dataout(1).datayused=datay;
                    dataout(1).nonnan=nonnan;
                else
                    dataout(1).trials=length(origdatax); %keep original trial ids 05/20/2020 for subseq sorting of groups
                end
                 dataout(1).r=r;
                dataout(1).p=psig;
           else           
                dataout=setfield(dataout,{1},['x' num2str(iax)],curx);
                dataout=setfield(dataout,{1},['y' num2str(iax)],cury);
                if ~nodata
                    dataout=setfield(dataout,{1},['datay' num2str(iax)],datax);
                    dataout=setfield(dataout,{1},['datay' num2str(iax)],datay);
                else
                  dataout=setfield(dataout,{1},['trials' num2str(iax)],length(datax));
                end
                dataout=setfield(dataout,{1},['r' num2str(iax)],r);
                dataout=setfield(dataout,{1},['p' num2str(iax)],psig);
           end
        end
        
        iax=iax+1;
    end
      
end

if ~noplot
text(axa{1},-50,axpos{1}(4)+50,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');
end
savename=[savepath '_' labels];
outdata=trlist;
%savefig(hf,savename);
if ~noplot
saveas(hf,savename,'jpg')
print(hf,savename,'-painters','-depsc');
end
end


