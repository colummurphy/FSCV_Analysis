function outdata=plottrialseries(sessnum,trlists,xinfos,xbinfos,binfos,plotparam,varargin)

sessiontypes={'bigreward','smallreward'};

sessid=num2str(sessnum);
event='targ';
trlistid=find(strcmp({trlists.sessid},sessid));
trlist=trlists(trlistid).list;
numtrials=length(trlist);
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,['trialseries_' sessid '_']);
sessnums=plotparam.sessnums;
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites);
uniquesites=unique({sites(1:end).site});
[dapair,lfppair]=getsitepairs(targdasites);
targses=find(strcmp({trialgrps.sessid},sessid));
trialinfo=trialgrps(targses).trialinfo;

%initialize trlist
for ix=1:length(trlist)
    trlist(ix).side=nan;
    trlist(ix).rts=nan;
    trlist(ix).da=nan;
    trlist(ix).lfp=nan;    
    trlist(ix).hr=nan;
end

ax(1).type='type';
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
axsiz=[1000,100];
axsizbinary=[1000,100];       %binary plots
fontsize=14;
mark={'o','o','.','x'};
marksize=[50,10,10,20];
axpos={};
axa={};
triallim=[];
xdata=1:numtrials;
bigtype=[];
smalltype=[];
targtype=[];
fixtype=[];
marksize=[150,100,20,30];
mark={'.','.','o','x'};

argnum=1;
siteda='';
sitelfp='';
bflag=0;
plotrt=0;
plothr=0;
ploteye=0;
plotnum=length(ax);
datms={};
lfptms={};
linfos={};
metric='targpeak';
metriclfp='targwin';
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'da'
            argnum=argnum+1;
            siteda=varargin{argnum};
            plotnum=plotnum+1;
            ax(plotnum).type='da';
        case 'datm'
            %values from datm, supplied instead of xinfos
            argnum=argnum+1;
            siteda=varargin{argnum};
            argnum=argnum+1;
            datms=varargin{argnum};
            plotnum=plotnum+1;
            ax(plotnum).type='da';
        case 'lfp'
            argnum=argnum+1;
            sitelfp=varargin{argnum};
            plotnum=plotnum+1;
            ax(plotnum).type='lfp';
        case 'lfptm'
            %values from betatm
            argnum=argnum+1;
            sitelfp=varargin{argnum};
            argnum=argnum+1;
            lfptms=varargin{argnum};
            plotnum=plotnum+1;
            ax(plotnum).type='lfp';
        case 'metric'
            %for datm values
            argnum=argnum+1;
            metric=varargin{argnum};
        case 'metriclfp'
            %for betatm values
            argnum=argnum+1;
            metriclfp=varargin{argnum};
        case 'rt'
            plotrt=1;
            plotnum=plotnum+1;
            ax(plotnum).type='rt';
        case 'hr'
            plothr=1;
            plotnum=plotnum+1;
            ax(plotnum).type='hr';
        case 'eye'
            bflag=1;
            ploteye=1;
            plotnum=plotnum+1;
            ax(plotnum).type='eye';
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};     %alignment event/window for values
        case 'tlim'
            argnum=argnum+1;
            triallim=varargin{argnum};
    end
    argnum=argnum+1;
end     
%set up plots
numplots=length(ax);
clf(hf,'reset');
set(hf,'color',[1 1 1]);
for ip=1:numplots
    axa{ip}=subplot(numplots,1,ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    if strcmp(ax(ip).type,'reward') || strcmp(ax(ip).type,'side') || strcmp(ax(ip).type,'error')
        %binary type
        axsize=axsizbinary;
    else
        axsize=axsiz;
    end
    set(axa{ip},'position',[xoff,figpos(4)-axsize(2)*ip-yoff-mar*(ip-1),axsize(1),axsize(2)])
    if ip>1
        set(axa{ip},'position',[xoff,axpos{ip-1}(2)-axsize(2)-mar,axsize(1),axsize(2)])
    end
    axpos{ip}=get(axa{ip},'position');    
    %ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
end
titletext=['sess ' sessid ' | ' event ];  
text(axa{1},-50,axpos{1}(4)+75,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');

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
    
%get hr  from binfos
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    hr=binfos(btarg).pulse;
    hrtrials=binfos(btarg).seltrials;
    for ida=1:length(hrtrials)
        tid=find(ismember(maplist(itype).trialnums,hrtrials(ida)));
        if ~isempty(tid)
            trlist(maplist(itype).trlist(tid)).hr=hr(ida);
        end
    end
    
%get hr  from binfos
%{
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    eye=binfos(btarg).eyed;
    etrials=binfos(btarg).seltrials;
    for ida=1:length(etrials)
        tid=find(ismember(maplist(itype).trialnums,etrials(ida)));
        %{
        if ~isequal(etrials(ida),maplist(itype).trialnums)
            %fewer trials in data variable than listed good trials, reset to match
            %measured good values
            maplist(itype).trlist=find(ismember([trlist.id],etrials(ida)+99) & ...
                strcmp({trlist.type},sessiontypes{itype}));
            maplist(itype).trialnums=etrials(ida);
        end
        %}
        if ~isempty(tid)
            trlist(maplist(itype).trlist(tid)).eye=eye(ida);
        end
    end    
    %}
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


for ip=1:length(ax)
    if strcmp(ax(ip).type,'type')
        bigtype=double(strcmp({trlist.type},'bigreward'));
        bigtype(bigtype==0)=nan;
        smalltype=double(strcmp({trlist.type},'smallreward'));
        smalltype(smalltype==0)=nan;
        targtype=double(strcmp({trlist.type},'targetbreak'));
        targtype(targtype==0)=nan;
        fixtype=double(strcmp({trlist.type},'fixbreak'));
        fixtype(fixtype==0)=nan;
        %plot y ax values according to side, & marker type to distinguish
        %type
        sidetype=[trlist.side];
        scatter(axa{ip},xdata,bigtype.*sidetype,marksize(1),mark{1},'markeredgecolor',[1 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,smalltype.*sidetype,marksize(2),mark{2},'markeredgecolor',[0 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,fixtype.*.5,marksize(3),mark{3},'markeredgecolor',[.7 0 1],...
            'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,targtype.*sidetype.*.5+.25,marksize(4),mark{4},'markeredgecolor',[0 .3 1],...
            'MarkerEdgeAlpha',.5,'linewidth',1);
        ylim(axa{ip},[-0.5 1.5]);
        yticks=[0 1]; yticklabel=['L';'R'];
        set(axa{ip},'ytick',yticks);
        set(axa{ip},'yticklabel',yticklabel);
        ylabel(axa{ip},'trial type');
        if ~isempty(triallim)
            set(axa{ip},'xlim',triallim);
        end
    end
    if strcmp(ax(ip).type,'da')
         %scatter(axa{ip},xdata,[trlist.da],marksize(4),mark{1},'markeredgecolor',[0 0 0],...
        %    'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,bigtype.*[trlist.da],marksize(1),mark{1},'markeredgecolor',[1 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,smalltype.*[trlist.da],marksize(2),mark{2},'markeredgecolor',[0 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        ylabel(axa{ip},['da - ' siteda]);
        if ~isempty(triallim)
            set(axa{ip},'xlim',triallim);
        end
    end
    if strcmp(ax(ip).type,'lfp')
       %  scatter(axa{ip},xdata,[trlist.lfp],marksize(4),mark{1},'markeredgecolor',[0 0 0],...
         %   'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.5,'linewidth',1);
         scatter(axa{ip},xdata,bigtype.*[trlist.lfp],marksize(1),mark{1},'markeredgecolor',[1 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,smalltype.*[trlist.lfp],marksize(2),mark{2},'markeredgecolor',[0 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
       %  aa=plot(axa{ip},xdata,[trlist.lfp],'linestyle','-','color',[1 0 0],...
         %   'linewidth',1);
       % aa.Color(4)=.4;
        ylabel(axa{ip},['beta - ' sitelfp]);
        if ~isempty(triallim)
            set(axa{ip},'xlim',triallim);
        end
    end
    if strcmp(ax(ip).type,'rt')
       %  scatter(axa{ip},xdata,[trlist.rts].*1e3,marksize(4),mark{1},'markeredgecolor',[0 0 0],...
        %    'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.5,'linewidth',1);
         scatter(axa{ip},xdata,bigtype.*[trlist.rts],marksize(1),mark{1},'markeredgecolor',[1 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,smalltype.*[trlist.rts],marksize(2),mark{2},'markeredgecolor',[0 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        ylabel(axa{ip},['rts (ms)']);
        if ~isempty(triallim)
            set(axa{ip},'xlim',triallim);
        end
    end    
    if strcmp(ax(ip).type,'hr')
         scatter(axa{ip},xdata,bigtype.*[trlist.hr],marksize(1),mark{1},'markeredgecolor',[1 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        scatter(axa{ip},xdata,smalltype.*[trlist.hr],marksize(2),mark{2},'markeredgecolor',[0 0 0],...
            'MarkerfaceAlpha',.5,'MarkerEdgeAlpha',.5,'linewidth',1);
        ylabel(axa{ip},['hr (bpm)']);
        if ~isempty(triallim)
            set(axa{ip},'xlim',triallim);
        end
    end    
    if strcmp(ax(ip).type,'eye')
         scatter(axa{ip},xdata,[trlist.eye],marksize(4),mark{1},'markeredgecolor',[0 0 0],...
            'MarkerfaceAlpha',.3,'MarkerEdgeAlpha',.5,'linewidth',1);
        ylabel(axa{ip},['pupil diam']);
        if ~isempty(triallim)
            set(axa{ip},'xlim',triallim);
        end
    end        
end


outdata=trlist;
savename=[savepath '_' event];
if ~isempty(triallim)
    savename=[savename '_tlim' num2str(triallim(1)) '_' num2str(triallim(2))];
end

savefig(hf,savename);
saveas(hf,savename,'jpg')
%print(hf,savename,'-painters','-depsc');
end


