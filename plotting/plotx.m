function plotdata=plotx(xinfo,trialinfo,plotparam,varargin)
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%rev 1/3/2019 to include new timing parameters from updated also xclust
argnum=1;
datype='dapos';
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'daneg'
            datype='daneg';
        case 'dareb'
            datype='dareb';
        case 'daall'
            datype='daall';
    end
    argnum=argnum+1;
end
fontsize=9;
savepath=plotparam.savepath;
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
dasites={xinfo(1:end).siteda};
lfprows={xinfo(1:end).sitelfp};
%lfpchs=unique(lfprows);
targdasites=unique(dasites);
figpos=[50,50,1200,800];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);

set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};

eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);

sessionrows={xinfo(1:end).sessiontype};
sessiontypes=unique({xinfo(1:end).sessiontype});

params={'damax','damaxts','damin','damints','dafallts','darisets',...
    'lfpmax','lfpmaxts','lfpmin','lfpmints','lfprisets','lfpfallts',...
    'lfppostmax','lfppostmaxts'};
params_mag={'damax','damin',...
    'lfpmax','lfpmin','lfppostmax','zlagcoef'};
params_xcov={'maxprelagts',...
    'minprelagts','maxpostlagts','minpostlagts','maxprecoef',...
    'minprecoef','maxpostcoef','minpostcoef'};
params_timing_damax={'delt_lfpmax_damax','delt_lfpfall_damax'...
    'delt_lfpmin_damax','delt_lfprise_damax',...
    'delt_lfppostmax_damax'};
params_timing_darise={'delt_lfpmax_darise','delt_lfpfall_darise'...
    'delt_lfpmin_darise','delt_lfprise_darise',...
    'delt_lfppostmax_darise'};   
params_timing_dafall={'delt_lfpmax_dafall','delt_lfpfall_dafall'...
    'delt_lfpmin_dafall','delt_lfprise_dafall',...
    'delt_lfppostmax_dafall'};
    
params=[params_mag params_xcov params_timing_damax...
   params_timing_darise params_timing_dafall];   

ygrp{1}=params_mag;
ygrp{2}=params_xcov;
ygrp{3}=params_timing_damax;
ygrp{4}=params_timing_darise;
ygrp{5}=params_timing_dafall;

ynam{1}='mag';
ynam{2}='xcov';
ynam{3}='timing_damax';
ynam{4}='timing_darise';
ynam{5}='timing_damin';
%plot for all groups, max da, max beta, etc.

%plot for all groups, max da, max beta, etc.
%single figure has single da-lfp pair (~6 pairs per session)
%separate figure for event tpe (2 types)
%single figure for groups within session (3 sessions)
%all(but switch), phase1, phase4, switch, postswitch, left/right,
%success/fail
%separate ax for magnitudes/timing/xcov (~36 figs each w/ 3 ax)

mark=[0 .2 .5; 0 .3 .4; 0 .5 .2];
markr=[.2 .5 .0; .3 .4 0; .5 .2 0];
%markr=[.2 .7 .0; .4 .3 0; .7 .1 0];

for ida=1:length(targdasites)
    %each da separate subplot  
    lfpsites={};
    daregion=contains(targdasites(ida),'c');     %1=='c'
    if daregion==1
        lfpsites=lfpchs(cgroup);
    else
        lfpsites=lfpchs(pgroup);
    end
for ilfp=1:length(lfpsites)    
for ievent=1:length(eventtypes)
for itype=1:length(sessiontypes)
    %each da-lfp pair separate figure
    %& each window event type & session type separate figure   
    targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);    %'dapos' or 'daneg' types
    %trial ids for this data group for all/phase1/phase4..
    trialtypes=trialinfo(itype).trialtypes; %get trial groups for session
    grptids={};     %trial nums selected to index for this group of data
    ranj={};    %ran jitter, depends on trial nums each group
    for itrial=1:length(trialtypes.names)
        grptids(itrial).name=trialtypes.names{itrial};
        grptids(itrial).nums=...
            find(ismember(curdata.trials,trialtypes.nums{itrial})==1);
        %remove switch trials for all trial conditions except switch itself
        %1/3/2018
        if ~strcmp(trialtypes.names{itrial},'switch')
            %take out switch trials
            swtrials=trialtypes.nums{strcmp(trialtypes.names,'switch')};
            allbutsw=curdata.trials(ismember(curdata.trials,trialtypes.nums{itrial})...
                & ~ismember(curdata.trials,swtrials));
            grptids(itrial).nums=...
                find(ismember(curdata.trials,allbutsw)==1);
        end
        ranj{itrial}=rand(1,length(grptids(itrial).nums))*.25;
    end
for yid=1:length(ygrp)
    %timing/magnitudes/xcov parameters to plot in diff figs    
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
    for ip=1:4
        axa{ip}=subplot(2,2,ip);   hold(axa{ip},'on');
        set(axa{ip},'xtick',1:length(grptids),'xticklabel',{grptids(:).name});
        set(axa{ip},'xlim',[0 length(grptids)+1]);
        set(axa{ip},'xTickLabelRotation',90)
    end
    labels={};
    titletext=[sessiontypes{itype} ' | ' eventtypes{ievent} ' | ' targdasites{ida} ' x ' lfpsites{ilfp}];  
    set(axa{1},'units','pixels')
    axpos=get(axa{1},'position');
    text(axa{1},-30,axpos(4)+30,titletext,'units','pixels','fontweight','bold');
    curgrp=ygrp{yid};
    cntp=1;
    count=1;
    plotdata=[];
for did=1:length(curgrp)
    %scroll parameters for group (ie. all timing variables)
    pid=ceil(cntp/2);        %plot axes, pairs
    grpdata=getfield(curdata,curgrp{did});
    xdata=[];
    xdata={grptids(:).name};
    gids={grptids(:).nums};
    ydata=[];
    xoff=0;
    if mod(cntp,2)
        %odd
        yyaxis(axa{pid},'left');
        xoff=-0.25;
        cmark=mark(1,:);
    else
        %even
        yyaxis(axa{pid},'right');
        xoff=0.25;
        cmark=markr(1,:);
    end
    set(axa{pid},'ycolor',cmark);
    for ix=1:length(xdata)
        plotdata(count).x=grptids(ix).name;
        plotdata(count).yvar=curgrp{did};
        plotdata(count).yavg=nan;
        plotdata(count).yci=nan;
        plotdata(count).ydata=nan;
        if ~isempty(grpdata(gids{ix}))
            ydata{ix}=grpdata(gids{ix});
            xplot=repmat(ix,1,length(ydata{ix}))+xoff+ranj{ix};
            scatter(axa{pid},xplot,...
                ydata{ix},15,'o','markeredgecolor',cmark,...
                'MarkerEdgeAlpha',.3,'linewidth',.5);
            ciy=nanstd(ydata{ix})./sqrt(length(ydata{ix}(~isnan(ydata{ix}))))*1.96;
            avgy=nanmean(ydata{ix});
            ciline=[avgy-ciy  avgy+ciy];
            plot(axa{pid},[ix ix]+xoff-.05,ciline,'-','linewidth',1,'color',cmark);
            plotdata(count).yavg=avgy;
            plotdata(count).yci=ciy;
            plotdata(count).ydata=ydata{ix};
        end
        count=count+1;
    end
    hy=ylabel(axa{pid},curgrp{did},'interpreter','none');
    if ~mod(cntp,2)
        %right axis        
        set(hy,'rotation',270,'units','pixels');
        hypos=get(hy,'position');
        set(hy,'position',[hypos(1)+15,hypos(2),hypos(3)]);  
    end
    cntp=cntp+1;
    set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)
end
savedir=[savepath targdasites{ida} 'x' lfpsites{ilfp} filesep];
if ~isdir(savedir)
    mkdir(savedir);
end
save([savedir ynam{yid} '_' datype '_' sessiontypes{itype} '_' eventtypes{ievent} ...
    '_' targdasites{ida} 'x' lfpsites{ilfp}],'plotdata');
savefig(figsess,[savedir ynam{yid} '_' datype '_' sessiontypes{itype} '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}]);
saveas(figsess,[savedir ynam{yid} '_' datype '_' sessiontypes{itype} '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}],'tif')
delete(findall(figsess,'type','text')) 
        
end
    end
end
end
end
end

end