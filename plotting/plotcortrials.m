function outdata=plotcortrials(sessnum,trlists,datms,lfptms,binfos,plotparam,varargin)
%4/2020, output relevant parameters for tabulating

sessiontypes={'bigreward','smallreward'};
condtype='all';
sessid=num2str(sessnum);
event='targ';
trlistid=find(strcmp({trlists.sessid},sessid));
trlist=trlists(trlistid).list;
numtrials=length(trlist);
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,['trialcorr_' sessid]);
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
     trlist(ix).lick=nan;   
     trlist(ix).rr=nan;   
     trlist(ix).rrstd=nan;   
     trlist(ix).rmssd=nan;   
end

figpos=[50,50,1400,900];
hf=figure('visible','off');     %figure for each channel
if ispc
hf=figure('visible','on');     %figure for each channel
end
set(hf,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',hf);    %set figure handle to current figure
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
xinfos={};
xbinfos={};
hfsametype=0;
metric='targpeak';
metriclfp='targwin';
win = [0 4];
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
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dax'
            %values from xinfo
            argnum=argnum+1;
            siteda=varargin{argnum};
            argnum=argnum+1;
            xinfos=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'datm'
            %values from datm, supplied instead of xinfos
            argnum=argnum+1;
            siteda=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'lfpx'
            argnum=argnum+1;
            sitelfp=varargin{argnum};            
            argnum=argnum+1;
            xinfos=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='lfp';
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
            bflag=1;
            ploteye=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='eye';
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
        case 'lwin'
            argnum=argnum+1;
            winl=varargin{argnum};
        case 'hrwin'
            argnum=argnum+1;
            winhr=varargin{argnum};
        case 'ewin'
            argnum=argnum+1;
            wine=varargin{argnum};
    end
    argnum=argnum+1;
end     

%set up plots
if isempty(plothist)
plothist=zeros(1,length(plotvars));       %1 means 1 trial in future value
end
if isempty(plottypes)
plottypes=repmat({'all'},1,length(plotvars));
end
clf(hf,'reset');
set(hf,'color',[1 1 1]);
numplots=sum([length(plotvars)-1:-1:1]);
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
titletext=['sess ' sessid ' | ' event ' | ' siteda ' ' metric ' | ' sitelfp ' ' metriclfp ];  
labels=[];
for ix=1:length(plotvars)
    if ix<length(plotvars)
    labels=[labels plotvars{ix} '_'];
    else
    labels=[labels plotvars{ix}];
    end        
end
labels=[labels '_' event];
if ~isempty(siteda)
    labels=[labels '_' metric '_' siteda];
end
if ~isempty(sitelfp)
    labels=[labels '_' metriclfp '_' sitelfp];
end
if ~isempty(winhr)
    labels=[labels '_winhr_' num2str(winhr(1)*rate) '_' num2str(winhr(2)*rate)];
end
if ~isempty(wine)
    labels=[labels '_wine_' num2str(wine(1)*rate) '_' num2str(wine(2)*rate)];
end
if ~isempty(winl)
    labels=[labels '_winl_' num2str(winl(1)*rate) '_' num2str(winl(2)*rate)];
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
    incswgrp=find(strcmp(trialinfo(itype).names,'switch')==1);
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
end
for ida=1:length(etrials)
    tid=find(ismember(maplist(itype).trialnums,etrials(ida)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).eye=eye(ida);
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
    if ~isequal(dtm{itype}.trialnums,maplist(itype).trialnums)
        %fewer trials in data variable than listed good trials, reset to match
        %measured good values
        maplist(itype).trlist=find(ismember([trlist.id],dtm{itype}.trialnums+99) & ...
            strcmp({trlist.type},sessiontypes{itype}));
        maplist(itype).trialnums=dtm{itype}.trialnums;
    end
    targtrials=find(ismember(dtm{itype}.trialnums,maplist(itype).trialnums));
    datatemp=getfield(dtm{itype},metric);
    dadata=datatemp(targtrials);
    for ida=1:length(dadata)
        trlist(maplist(itype).trlist(ida)).da=dadata(ida);
    end
    end

    if ~isempty(sitelfp) 
        %targeted beta values from betatm
    lfptm=lfptms{targses}{itype};
    for ii=1:length(lfptm)
        if strcmp(lfptm{ii}.site,sitelfp)
            lch=ii;
        end
    end
    ltm{itype}=lfptm{lch};           %get all da targeted values for current trial type big/small
    if ~isequal(ltm{itype}.trialnums,maplist(itype).trialnums)
        %fewer trials in data variable than listed good trials, reset to match
        %measured good values
        maplist(itype).trlist=find(ismember([trlist.id],ltm{itype}.trialnums+99) & ...
            strcmp({trlist.type},sessiontypes{itype}));
        maplist(itype).trialnums=ltm{itype}.trialnums;
    end
    targtrials=find(ismember(ltm{itype}.trialnums,maplist(itype).trialnums));
    ltemp=getfield(ltm{itype},metriclfp);
    ldata=ltemp(targtrials);
    for ida=1:length(ldata)
        trlist(maplist(itype).trlist(ida)).lfp=ldata(ida);
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
iax=1;
for ip=1:length(plotvars)
    cury=plotvars{ip};    
    remplots=numplots-ip;
    for iip=ip+1:length(plotvars)
        cla(axa{iax});
        datax=[];
        datay=[];
        curx=plotvars{iip};
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
        xoutthres=nanmean(datax)+4*nanstd(datax);
        youtthres=nanmean(datay)+4*nanstd(datay);
        remx=[];
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
        nonnan=find(~isnan(datax) & ~isnan(datay));
        if hfsametype
            nonnan=find(~isnan(datax) & ~isnan(datay) & (~isnan(datasmallx) | ~isnan(databigx)) & (~isnan(datasmally) | ~isnan(databigy)));
        end
        datax=datax(nonnan);
        datay=datay(nonnan);
        databig=bigtype(nonnan);
        datasmall=smalltype(nonnan);
        if strcmp(plottypes{iip},'all')
            scatter(axa{iax},databig.*datax,databig.*datay,20,'o','markeredgecolor',[1 0 0],'markerfacecolor',[1 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);
            scatter(axa{iax},datasmall.*datax,datasmall.*datay,20,'o','markeredgecolor',[0 0 0],'markerfacecolor',[0 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);
        elseif strcmp(plottypes{iip},'big')

             datax=databig.*datax;
             datay=databig.*datay;
             xoutthres=nanmean(datax)+4*nanstd(datax);
            youtthres=nanmean(datay)+4*nanstd(datay);
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
             nonnan=find(~isnan(datax) & ~isnan(datay));
             datax=datax(nonnan);
             datay=datay(nonnan);
                         scatter(axa{iax},datax,datay,20,'o','markeredgecolor',[1 0 0],'markerfacecolor',[1 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);
        elseif strcmp(plottypes{iip},'small')
 
             datax=datasmall.*datax;
             datay=datasmall.*datay;
            xoutthres=nanmean(datax)+4*nanstd(datax);
            youtthres=nanmean(datay)+4*nanstd(datay);
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
             nonnan=find(~isnan(datax) & ~isnan(datay));
             datax=datax(nonnan);
             datay=datay(nonnan);
                         scatter(axa{iax},datax,datay,20,'o','markeredgecolor',[0 0 0],'markerfacecolor',[0 0 0],...
                'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.3,'linewidth',1);  
        end
        ylabel(axa{iax},laby);
        xlabel(axa{iax},labx);
        [r,psig]=corr(datax',datay');
        [p,s]=polyfit(datax,datay,1);
        slope=p(1);
        intercep=p(2);
        yfit=slope*datax+intercep;
        yresid=datay-yfit;
        ssresid=sum(yresid.^2);
        sstotal=(length(datay)-1)*var(datay);
        rsq=1-ssresid/sstotal;
        fitline=plot(axa{iax},datax,yfit,'-','color',[0 0 0]);
        fitline.Color(4)=0.5;
        xsiz=get(axa{iax},'position');
        ylims=get(axa{iax},'ylim');
        text(axa{iax},xsiz(3)-120,xsiz(4)-50,...
            {['r: ' num2str(r)],...
            ['p: ' num2str(psig)]},'color',[0 0 0],'units','pixels','FontSize',fontsize)
                set(findall(axa{iax},'-property','FontSize'),'FontSize',fontsize)
        
        dataout.x{ix}=curx;
        dataout.y{ix}=cury;

        iax=iax+1;
    end
      
end
text(axa{1},-50,axpos{1}(4)+50,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');

savename=[savepath '_' labels];
outdata=trlist;
%savename=[savepath '_' event];
savefig(hf,savename);
saveas(hf,savename,'jpg')
print(hf,savename,'-painters','-depsc');
end


