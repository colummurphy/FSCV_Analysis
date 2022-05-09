function outdata=plotmultxtraces(sessnum,xinfo,plotparam,varargin)
%called from plotmutliple
%trialinfo must be given in same order as sesstypes
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%same as plotx, but plot big/small/targ on same plots
%for plotting multiple sessions
%1/3/2018, udpates as in plotx/plotxallsess
%PLOT TRACES for DA GROUPS
sessiontypeids={'bigreward','smallreward','targetbreak','fixbreak'};
outdata={};
outdata(1).sessnum=sessnum;
fontsize=15;
figpos=[50,50,1000,800];
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-2 2];     %+/-2 s from aln idx
interval=1;       %in seconds
pad=.2;         %.2 sec pad
ext=.8;         
argnum=1;
sortda=0;
tstraces={};
evt={};
lfpsites={};
sorttype='damax';
plotmarks=0;
sortrts=0;
ratelfp=1000;
datype='goodtrials';
sesstypes={};
event='targ';
tstype='goodtrials';
plotz=0;
sessid=num2str(sessnum);
bplot=0;
btarg=[];
binfo={};
plotmax=0;
plotmin=0;
plotmean=0;
plotimmean=0;
plotabs=0;
datms={};
smoothlfp=0;
smoothpoints=5;
winb=[];
plotrrstd=0;
zbehav=0;
numplots=2;
xflag=0;
tcovflag=0;
dasites={};
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'xcov'
            xflag=1;
            figpos=[50,50,1400,800];
            numplots=3;
        case 'tcov'
            %cross covariance of da and beta across all trials for pairs of
            %time points
            tcovflag=1;
            figpos=[50,50,1400,800];
            numplots=3;
        case 'daneg'
            datype='negtrials';
            tstype='negtrials';
        case 'dapos'
            datype='postrials';
            tstype='postrials';
        case 'daall'
            datype='goodtrials';
            tstype='goodtrials';
        case 'event'
            argnum=argnum+1;    
            event=varargin{argnum};         %event period, eg. fix, targ, outcome
        case 'datm'
            argnum=argnum+1;
            datms=varargin{argnum};
        case 'sort'
            %sort da increasing
            sortda=1;
            argnum=argnum+1;
            sorttype=varargin{argnum};
            if contains(sorttype,'rts')
                sortrts=1;
            end
        case 'ext'
            argnum=argnum+1;
            ext=varargin{argnum};
        case 'win'
            %user specified window +/- s
            argnum=argnum+1;
            win=varargin{argnum};
        case 'plotz'
            plotz=1;           
        case 'plotmarks'
            plotmarks=1;
            tstraces{1}.tstrace={'damaxts'};
            tstraces{2}.tstrace={'lfppostmaxts','lfpmints'};
        case 'lfpmints'
            tstraces{2}.tstrace={'','lfpmints'};
        
        case 'lfps'
            argnum=argnum+1;
            lfpsites=varargin{argnum};      %user provides lfp sites
        case 'dasites'
            argnum=argnum+1;
            dasites=varargin{argnum};
        case 'ttypes'
            argnum=argnum+1;
            ttypes=varargin{argnum};        %condition types grouped eg. {{'big','left'},{'small','left','aftersm'}}
        case 'bplot'
            %plot beh var instead of lfp, 2 arguments needed, xbinfo, and
            %sitelfp
            argnum=argnum+1;
            binfo=varargin{argnum};
            bplot=1;
            argnum=argnum+1;
            btarg=varargin{argnum};
        case 'plotmax'
            plotmax=1;
        case 'plotmin'
            plotmin=1;
        case 'plotmean'
            plotmean=1;
        case 'plotimmean'
            plotimmean=1;
        case 'plotabs'
            plotabs=1;
        case 'pad'
            argnum=argnum+1;
            pad=varargin{argnum};
        case 'smoothlfp'
            smoothlfp=1;
            argnum=argnum+1;
            smoothpoints=varargin{argnum};
        case 'winb'
            %win for beh
            argnum=argnum+1;
            winb=varargin{argnum};
        case 'rrstd'
            plotrrstd=1;
        case 'zbeh'
            zbehav=1;
    end
    argnum=argnum+1;
end
savepath=plotparam.savepath;

mark=[0 0 0; 1 0 0; 0 .7 0;.7 0 0; 0.4 0 .7];
trialinfo=plotparam.trialgrps(strcmp({plotparam.trialgrps.sessid},sessid)).trialinfo;
if ~iscell(win)
    wintemp=win;
    win={};
    win{1}=wintemp;
    win{2}=wintemp;
end
tsx=[win{1}(1):1/rate:win{1}(2)];
tsxlfp=[win{1}(1):1/ratelfp:win{1}(2)];

condlabels=[];
condfull=[];
for ic=1:length(ttypes)
    condlabels=[condlabels '_' ttypes{ic}];
    condfull=[condfull ', ' ttypes{ic}];
    if ic==1
        condlabels=ttypes{ic};
        condfull=ttypes{ic};
    end
end
targdasites=plotparam.dasites;

sites=getsites(sessnum,targdasites);
targdasites=unique({sites.probeid});
if ~isempty(dasites)
    targdasites=dasites;
    sites=getsites(sessnum,targdasites);

end
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));
[dapair,lfppair]=getsitepairs(targdasites);
targses=find(strcmp({plotparam.trialgrps.sessid},sessid));

lfpchs=plotparam.lfpchs;
userlfp=0;
if ~isempty(lfpsites)
    lfpchs=lfpsites;
    userlfp=1;
end
targlsites={};
lsites={};
%get corresponding lfp paired site for each da site for xinfo data
targlsites=lfpchs;
lsites=getlfpsites(sessnum,targlsites);
sessiontypes=unique({xinfo(1:end).sessiontype});
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);

set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axpos={};
axsiz=[750,320];
axoff=100;
targrow=[];
lfptarg=[];
lfpsite=[];
xcovdata=[];
xdata=[];
tsx=[];
psigdata=[];
for ida=1:length(targdasites)
    outdata(ida).dasite=targdasites(ida);
    outdata(ida).sessnum=sessnum;
    %each da separate figure  
    lfpsites={};
    daregion=contains(targdasites(ida),'c');     %1=='c'
     %get ch # for current sess and probe id
    targsi=find([sites.sessnum]==sessnum & ...
    strcmp({sites.probeid},targdasites(ida)));
    ich=sites(targsi(1)).ch;        %ch # defined by fscv recording ch       
    %get lfp pair for da site, for xinfo data
    daid=find(strcmp(dapair,targdasites(ida)));
    if ~bplot
        if userlfp
            %user provided
            targl=find([lsites.sessnum]==sessnum & ...
                strcmp({lsites.probeid},lfpchs));
        if ~isempty(targl)
            lfptarg=lfpchs;
        end
        else                
            lfptarg=lfppair{daid(1)};           %JUST FIRST GOOD LFP PAIR
         %get lfp site pair name iiregardles of lfp type value, if
            %exits in xinfo
        end        
            targl=find([lsites.sessnum]==sessnum & ...
                strcmp({lsites.probeid},lfptarg));
            %check if in xinfo
         targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
            contains({xinfo.sitelfp},lfptarg) & ...
            strcmp({xinfo.event},event)) & ...
            contains({xinfo.sessionid},sessid)==1);    
        if isempty(targrow)
            %target paired lfp not in xinfo, find another
            targrows=find((contains({xinfo.siteda},targdasites(ida))  & ...
                strcmp({xinfo.event},event)) & ...
                contains({xinfo.sessionid},sessid)==1); 
            xinfolfps=unique({xinfo(targrows).sitelfp});
            xinfolfps=intersect({lsites.probeid},xinfolfps);
                        lfptarg=xinfolfps{1};
            targl=find([lsites.sessnum]==sessnum & ... 
                strcmp({lsites.probeid},lfptarg));
            lfpsite=lsites(targl).site; 
        end
    
     if isempty(targl)
          % find another lfppair lfptarg
       lfpsitesinsess=unique({xinfo(find(contains({xinfo.sessionid},sessid))).sitelfp});
       lfptarg=intersect({lfppair{daid}},lfpsitesinsess);
       if length(lfptarg)>1
           lfptarg=lfptarg(1);
       end
       disp(['new lfp targ for sess # ' sessid ' : ' lfptarg{:}]); 
       targl=find([lsites.sessnum]==sessnum & ...
        strcmp({lsites.probeid},lfptarg));
     end
    lfpsite=lsites(targl).site; 
    else
        lfptarg=btarg;
        lfpsite=btarg;
    end
   
    outdata(ida).lfpsite=lfpsite;

da=[];
lfp=[];
tsdata={};
sortdata=[];
    tsdata{1}=[];
    tsdata{2}=[];
    curdata={};
    
%get trials corresponding to groups of trial conditions in current group ttypes
%get trial types (big/small/break) for current group
targtrialtypes=[any(strcmp('big',ttypes))  ...
    any(strcmp('small',ttypes)) ...
    any(strcmp('targetbreak',ttypes)) ...
    any(strcmp('fixbreak',ttypes)) ]; %get logical array for if want big/small/targ/fix
targtrialtypes=find(targtrialtypes==1);
nottypes=find(~strcmp(ttypes,'big') & ...
    ~strcmp(ttypes,'small') &...
    ~strcmp(ttypes,'targetbreak') &...
   ~strcmp(ttypes,'fixbreak'));
%get trial #'s for current condition, ie. not trial type identifier
%get da values for current group condition
for itype=targtrialtypes
    itypeid=find(contains(sessiontypeids,sessiontypes(itype)));
    trialnums=[];
    davals=[];
    targrowb=[];
   curbdata={};
    %get xinfo data
    targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfptarg) & ...
        strcmp({xinfo.event},event)) & ...
        contains({xinfo.sessionid},sessid) &...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);    
    if bplot
        targrowb=find((contains({binfo.siteda},targdasites(ida)) & ...
            contains({binfo.sitelfp},btarg) & ...
            strcmp({binfo.event},event)) & ...
            contains({binfo.sessionid},sessid) &...
            contains({binfo.sessiontype},sessiontypes(itype))==1);  
        curbdata=getfield(binfo(targrowb),'daall'); 
        targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
            strcmp({xinfo.event},event)) & ...
            contains({xinfo.sessionid},sessid) &...
            contains({xinfo.sessiontype},sessiontypes(itype))==1);
        targrow=targrow(1);         %probably multiple since lfp not specified
        curdata=getfield(xinfo(targrow),'daall');    %'dapos' or 'daneg' types get below
    else
        if ~isempty(targrow)
        curdata=getfield(xinfo(targrow),'daall');    %'dapos' or 'daneg' types get below
        end
    end
for tt=nottypes
    curcond=ttypes{tt};
    targt=find(contains(trialinfo(itype).names,curcond)==1);
    trialnums=[trialnums trialinfo(itype).nums{targt}];  
end
trialnums=unique(trialnums);        %trial nums for current trial type big/small
if ~isempty(curdata)
    trialids=[];
    trialidsb=[];
    sortdatat=[];
    trialnums=intersect(trialnums,getfield(xinfo(targrow),datype));     %get dapos/daall/neg trials
    if ~isempty(datms)
        datm=datms{targses}{itypeid}{ich};
        trialnums=intersect(trialnums,datm.trialnums);
    end
    if bplot
        trialnums=intersect(trialnums,getfield(binfo(targrowb),datype));
        trialidsb=find(ismember(curbdata.trials,trialnums)==1); %trial ids for targeted metric
    end
    trialnums=intersect(curdata.trials,trialnums);
    trialids=find(ismember(curdata.trials,trialnums)); %trial ids for targeted metric
   % davals=getfield(datarg,targval);    %all  values for good trials for sess/type
    alnidx=curdata.mididx;
    wins=[alnidx+win{1}(1)*rate:alnidx+win{1}(2)*rate];       %time window to plot
    dat=curdata.datracesaln(trialids,wins);                             %signals aligned to aln idx
    winslfp=[alnidx+win{2}(1)*rate:alnidx+win{2}(2)*rate];       %time window to plot
    lfpt=curdata.lfptracesaln(trialids,winslfp); 
    if bplot
        lfpt=curbdata.lfptracesaln(trialidsb,winslfp); 
        if ~isempty(winb)
            %alternative windwo for beh
            winslfp=[alnidx+winb(1)*rate:alnidx+winb(2)*rate];       %time window to plot
            lfpt=curbdata.lfptracesaln(trialidsb,winslfp); 
            if plotrrstd
                hr=curbdata.lfptracesaln(trialidsb,:); 
                hrwin=hr(:,winslfp);
                hrtemp=nanmean(hrwin,2);
                rrtemp=1./hr.*60.*1e3; 
                rrdata=rrtemp(:,winslfp);
                lfpt=rrdata;
            end
        end
    end
    if sortrts
        sortdatat=getfield(xinfo(targrow),'trt');
        if strcmp(sorttype,'fixrts')
        sortdatat=getfield(xinfo(targrow),'frt');
        end
        origseltrialids=intersect(curdata.trials,origtrialids);
        rttrials=intersect(xinfo(targrow).goodtrials,origseltrialids);
        sortdatat=sortdatat(rttrials);
    else
        if isempty(datms)
        sortdatat=getfield(curdata,sorttype);   %sort by metric acquired in xinfos variable
        sortdatat=sortdatat(trialids);
        else
            sortdatat=getfield(datm,sorttype);
            trialids2=find(ismember(datm.trialnums,trialnums));
            sortdatat=sortdatat(trialids2);
        end
    end   
    if ~isempty(tstraces)
        for id=1:2
            %da & lfp
            tstrace=tstraces{id}.tstrace;
            tsdatat={};
            %plot overlying ts parameter
            for it=1:length(tstrace)
                %max first row, min 2nd row
                if ~isempty(tstrace{it})
                dataforts=getfield(xinfo(targrow),'daall');
                tsdatatemp=getfield(dataforts,tstrace{it})-alnidx-win{ida}(1)*rate+1;
                tsdatat{id}(it,:)=tsdatatemp(trialids);
                end
            end
            tsdata{id}=[tsdata{id} tsdatat{id}];
        end
    end
    sortdata=[sortdata sortdatat];
    da=[da; dat];
    lfp=[lfp; lfpt];
   % yvals{pid}{sitenum}{icond}=[yvals{pid}{sitenum}{icond} davals(trialids)];         %da values for specific group of conditions (ie. trial nums)
end
end
if ~isempty(targrow)
if strcmp(lfpsite,'eye')
    %flip sign so values match actual condition of eye
    %normalize to median of all data
        lfp=-lfp;
     normval=abs(nanmedian(lfp));        %median of all values for all plot conditions for current probe
     lfp=lfp./normval;
end
%complete cycle through all conditions fo values da/lfp/sortdata
if sortda
    [maxvals,sortids]=sort(sortdata);
    sortids=sortids(~isnan(maxvals));       %remove nan values (usually placed as max)
    da=da(sortids,:);
    lfp=lfp(sortids,:);
    tsdatatemp=tsdata;
    tsdata={};
    if ~isempty(tstraces)
        for id=1:2
            tstrace=tstraces{id}.tstrace;
            for it=1:length(tstrace)
                %max first row, min 2nd row
                tsdata{id}(it,:)=tsdatatemp{id}(it,sortids);
            end
        end
    end
end
if plotz
    %da=zscore(da,[],2);            %ZSCORE ONLY NON-NAN IDS
    for it=1:size(da,1)
        %one by one to ignore nan's
        nonnan=~isnan(da(it,:));
        daz=zscore(da(it,nonnan));
        da(it,nonnan)=daz;
    end
    lfp=zscore(lfp,[],2);
end
if xflag    
        xdata={lfp da};
        xrates=[10 10];      %lfp rate & da rate, both at 10 in xinfos
        plotparam.alnts=0;
        vartype='none';  
     xcovdata=xvardata(xdata,xrates,plotparam,[],[],...
        'type',vartype,'noplot','xwin',[0 abs(win{1}(1))+abs(win{1}(2))]);    
    xdata=xcovdata.xcovda;
    tsx=xcovdata.tsx;
    win{3}=[-win{1}(end) win{1}(end)];
end
if tcovflag    
            xdata={lfp da};
        xrates=[10 10];      %lfp rate & da rate, both at 10 in xinfos
        plotparam.alnts=0;
        vartype='none';  
     xcovdata=xvardata(xdata,xrates,plotparam,[],[],...
        'type',vartype,'noplot','xwin',[0 abs(win{1}(1))+abs(win{1}(2))],'tcov');    
    xdata=xcovdata.xcovda;
    tsx=xcovdata.tsx;
    psigdata=xcovdata.pdata;
    win{3}=win{1};
end
freqband=getfield(xinfo(targrow),'freq');
clf(figsess,'reset');
set(figsess,'color',[1 1 1]);
axpos={};
pdata=[];
clabel=[];
for ip=1:numplots       
    axa{ip}=subplot(1,numplots,ip);
    hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
    if size(da,1)>1
    if ip==1
        pdata=da;
        outdata(ida).da=pdata;
        title(axa{ip},[datype ' | ' condfull ' | ' event ' | da ' targdasites{ida}]);
        clabel='\DeltaDA (nM)';
        if plotz
            clabel='z-score DA';
        end
        ylabel(axa{ip},'trial #')
        if sortda
            ylabel(axa{ip}, [sorttype ' sorted trials']);
        end
    end
    if ip==2
        pdata=lfp;
        if ~bplot
         title(axa{ip},['\beta-lfp ' lfpsite]);
        else
            title(axa{ip},[lfpsite]);
        end
         clabel=['\beta-LFP ' num2str(freqband(1)) ...
             '-' num2str(freqband(2)) ' hz'];
    end  
    if ip==3
        pdata=xdata;       
            title(axa{ip},'xcov');        
    end  
    if ip==2
        if plotmax
            %plot max in trilas
            pdata=max(pdata,[],2);
        end
        if plotmean
            pdata=nanmean(pdata,2);
        end
        if plotimmean
            pdata=nanmean(pdata(:,pad*rate:ext*rate),2);
        end
        if plotmin
            pdata=min(pdata,[],2);
        end
        if plotabs
            pdata=max(abs(pdata,[],2));
        end
        if zbehav
            %z score behav only eg. lick, across trials instead of time as done
            %above
           pdata=(pdata-nanmean(pdata))./nanstd(pdata);
        end
        if plotrrstd
               % hr=pdata;            
               % hrtemp=nanmean(hr,2);
                %rrtemp=1./hr.*60.*1e3; 
               % rrdata=rrtemp(:,rrwin(1)*fs+alnidx:rrwin(2)*fs+alnidx);
                pdata=nanstd(pdata,[],2);
        end    
        if smoothlfp
            if smoothpoints(1)~=0
            pdata=smoothdata(pdata,2,'gaussian',smoothpoints(1),'omitnan');      %smooth 5 data points along trial axis
            end
            pdata=smoothdata(pdata,1,'gaussian',smoothpoints(2),'omitnan');      %smooth 5 data points along trial axis
        end
    end
    if size(pdata,2)>1
        %more than one time point
        imagetrials=image(axa{ip}, pdata,'cdatamapping','scaled');
        set(axa{ip},'YDir','reverse')        %flip y-axis values so first trial on top
        artTime=isnan(pdata);   %find artifact points (nan periods)
        artTime=abs(artTime-1);         %make alpha data mask by inverting 1/0's
        artTime2=artTime;
        maskGray=artTime2==0;             %find Zero indices representing artifact mask
        maskGray=maskGray*.15;            %make gray rather than white default by making non-zero
        artTime=artTime+maskGray;
        set(imagetrials, 'AlphaData', artTime);    
        if ip==3
            %mask insig data
          
        end
        if ip<3
            xticklabels=min(win{ip}):interval:max(win{ip});
            xticklabels=round(xticklabels.*rate)./rate;
            xticklabels=num2str(xticklabels');
            xticks=1:round(interval*rate):size(pdata,2);
            set(axa{ip},'xtick',xticks,'xticklabel',xticklabels);
            set(axa{ip},'xlim',[0 size(pdata,2)]);
        else
            xticklabels=min(win{ip}):interval:max(win{ip});
            xticklabels=round(xticklabels.*rate)./rate;
            xticklabels=num2str(xticklabels');
            xticks=1:round(interval*rate):size(pdata,2);
            set(axa{ip},'xtick',xticks,'xticklabel',xticklabels);
            set(axa{ip},'xlim',[0 size(pdata,2)]);
        end
        yints=round(size(pdata,1)/10);      %row intervals
        trialsort=sort(1:size(pdata,1));
        yticks=trialsort(1:yints:end);
        set(axa{ip},'tickdir','out','box','off')
        xlabel(axa{ip},'time (s)')
        set(axa{ip},'ylim',[1 size(pdata,1)]);
        origpos=getpixelposition(axa{ip});      %get original position 
        set(axa{ip},'ytick',yticks);
        h1=colorbar(axa{ip},'southoutside');  
        set(h1,'units','pixels');
        set(axa{ip},'Units','Pixels','Position', [axpos{ip}(1) axpos{ip}(2) axpos{ip}(3) axpos{ip}(4)]);
        set(h1,'position',[axpos{ip}(1) 30 100 10]);
        %climsmax=nanmean(pdata(mididx,:))+nanmean(nanstd(pdata,[],2))*3;
        %climsmin=nanmean(pdata(mididx,:))-nanmean(nanstd(pdata,[],2))*1.5;
        scaleslfp=nanmean(nanstd(lfp,[],2));
        scalesda=nanmean(nanstd(da,[],2));
        climsmax=scalesda*4;
        climsmin=-scalesda*1.5;
    if ip==2
       % climsmax=nanmean(pdata(mididx,:))+nanmean(nanstd(pdata,[],2))*4;
        %climsmin=nanmean(pdata(mididx,:))-nanmean(nanstd(pdata,[],2))*1;
        climsmax=scaleslfp*5;
        climsmin=scaleslfp*1;
        set(axa{ip},'ytick',[]);
        if ~bplot
        outdata(ida).beta=pdata;
        end
    end
    if bplot && ip==2
        climsmax=max(nanmean(pdata,1))+median(nanstd(pdata,[],1))*3;
            climsmin=min(nanmean(pdata,1))-median(nanstd(pdata,[],1))*3;
            outdata(ida).beh=pdata;
    end
    %climsmax=max(nanmean(pdata,1))+median(nanstd(pdata,[],1))*3;
    %climsmin=min(nanmean(pdata,1))-median(nanstd(pdata,[],1))*3;
    if ip<3
    if ~plotz 
    if ~isempty(climsmin) && ~isempty(climsmax) && ~isnan(climsmax) && ~isnan(climsmin)
        set(axa{ip},'clim',[climsmin climsmax])
    end
    else
        if ip==2 && ~bplot
       climsmin=-1;
        climsmax=3;        
        set(axa{ip},'clim',[climsmin climsmax])
        end
    end
    end
    if ~isempty(tstraces) && plotmarks
        tstrace=tstraces{ip}.tstrace;
        markc=[.75 .15 .75; 0 0 0];
        markc=[1 .75 1; 0 0 0];
        marka=[1 .75];
        markt={'s','o'};
        marks=[20,15];
        markl=[2 1.5];
        %plot overlying ts parameter
        for it=1:length(tstrace)
            if ~isempty(tstrace{it})
            %dataforts=getfield(xinfo(targrow),tstype);
           % tsdatatemp=getfield(dataforts,tstrace{it})-alnidx-win(1)*rate+1;
            tsdataplot=tsdata{ip}(it,:);
            meants=nanmean(tsdataplot);
            stdts=nanstd(tsdataplot);
            tsci=stdts./sqrt(length(tsdataplot))*1.96; 
            for jj=1:size(da,1)
                %ama=plot(axa{ip},[tsdata(jj)],jj,'Color',markc(it,:),'LineWidth',2,'LineStyle',':','marker','s','markersize',4);           
                scatter(axa{ip},tsdataplot(jj),jj,marks(it),markt{it},'markeredgecolor',markc(it,:),'MarkerEdgeAlpha',marka(it),'linewidth',markl(it));      
            end
            end

        end
    end
    else 
        %single time point averaged/max
       cla(axa{ip})
        set(axa{ip},'YDir','reverse')        %flip y-axis values so first trial on top 
    set(axa{ip},'tickdir','out','box','off')
    xlims=[nanmean(pdata)-nanstd(pdata)*4 nanmean(pdata)+nanstd(pdata)*4];
    plot(axa{ip},[nanmedian(pdata) nanmedian(pdata)], [1 length(pdata)],'-k','linewidth',.75);
  %  if ~strcmp(lfpsite,'lick')
    set(axa{ip},'xlim',xlims);
   % end
    set(axa{ip},'ylim',[1 length(pdata)]);
    origpos=getpixelposition(axa{ip});      %get original position 
    set(axa{ip},'ytick',[]);
    set(axa{ip},'ycolor','none');
    set(axa{ip},'Units','Pixels','Position', [axpos{ip}(1)-75 axpos{ip}(2) 300 axpos{ip}(4)]);
     scatter(axa{ip}, pdata,1:length(pdata),30,'o','markeredgecolor',[0 0 0],'markerfacecolor',[0 0 0],'markeredgealpha',.2,'markerfacealpha',.2,'linewidth',1);
    end
    set(findall(axa{ip},'-property','FontSize'),'FontSize',fontsize)
    set(h1,'fontsize',8);
    end
    if ip==3
        outdata(ida).xvar=pdata;
        outdata(ida).p=psigdata;
    end
end
%savedir=[savepath 'allsess_traces'  filesep sessid filesep];
savename=[savepath  'traces_sess' sessid '_' datype '_' condlabels '_' event ...
    '_' targdasites{ida} 'x' lfpsite];
if plotz
    savename=[savename '_z'];
end
if sortda    
    savename=[savename '_sorted_' sorttype];
end
if plotmarks
    savename=[savename '_marks'];
end
if plotmean || plotmax || plotmin || plotabs
    savename=[savename '_sints'];
end
if plotimmean
    savename=[savename '_immean_' num2str(pad*rate) '_' num2str(ext*rate) 'win'];
end
if smoothlfp
    savename=[savename '_smooth_' num2str(smoothpoints(1)) '_' num2str(smoothpoints(2)) 'pts' ];
end
if ~isempty(winb)
      savename=[savename '_winb_' num2str(winb(1)*rate) '_' num2str(winb(2)*rate) 's' ];
end 
if plotrrstd
     savename=[savename '_rrstd'];
end  
if xflag
    savename=[savename '_xcov'];
end
savefig(figsess,savename);
saveas(figsess,savename,'tif')
print(figsess,savename,'-painters','-depsc');
delete(findall(figsess,'type','text')) 
end
end
end
