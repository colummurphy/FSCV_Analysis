function plotdata=plotxallsess(xinfo,trialinfo,plotparam,varargin)
%2/17/2019, just do selective metrics
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%same as plotx, but plot big/small/targ on same plots
%1/3/2018 updated as with plotx same date
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
axb={};

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
    'delt_lfppostmax_damax','lfppostmaxts','damaxts','lfpmints'}; 
params_timing_darise={'delt_lfpmax_darise','delt_lfpfall_darise'...
    'delt_lfpmin_darise','delt_lfprise_darise',...
    'delt_lfppostmax_darise','lfppostmaxts','darisets','lfprisets'};   
params_timing_dafall={'delt_lfpmax_dafall','delt_lfpfall_dafall'...
    'delt_lfpmin_dafall','delt_lfprise_dafall',...
    'delt_lfppostmax_dafall','lfppostmaxts','dafallts','lfpfallts'};
 params_timing_damin={'delt_lfpmax_damin','delt_lfpfall_damin'...
    'delt_lfpmin_damin','delt_lfprise_damin',...
    'delt_lfppostmax_damin','lfppostmaxts','damints','lfpmints'};

params=[params_mag params_xcov params_timing_damax...
   params_timing_darise params_timing_dafall params_timing_damin];   

ygrp{1}=params_mag;
ygrp{2}=params_xcov;
ygrp{3}=params_timing_damax;
%ygrp{4}=params_timing_darise;
%ygrp{5}=params_timing_dafall;
%ygrp{6}=params_timing_damin;

ynam{1}='mag';
ynam{2}='xcov';
ynam{3}='timing_damax';
%ynam{4}='timing_darise';
%ynam{5}='timing_dafall';
%ynam{6}='timing_damin';
%plot for all groups, max da, max beta, etc.
%single figure has single da-lfp pair (~6 pairs per session)
%separate figure for event tpe (2 types)
%single figure for groups within session (3 sessions)
%all(but switch), phase1, phase4, switch, postswitch, left/right,
%success/fail
%separate ax for magnitudes/timing/xcov (~36 figs each w/ 3 ax)

mark=[0 .2 .7; 0 .5 .3; 0 .7 .1];
markr=[.2 .7 .0; .4 .3 0; .7 .1 0];
typeoff=[0.0 0.1 0.2];
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
for yid=1:length(ygrp)    
    %timing/magnitudes/xcov parameters to plot in diff figs    
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
    trialtypes=trialinfo(1); %get trial groups for session
    for ip=1:4
        axa{ip}=subplot(2,2,ip);   hold(axa{ip},'on');
        set(axa{ip},'xtick',1:length(trialtypes.names),'xticklabel',trialtypes.names);
        set(axa{ip},'xlim',[0 length(trialtypes.names)+1]);
        set(axa{ip},'xTickLabelRotation',90)
    end
    labels={};
    titletext=[eventtypes{ievent} ' | ' targdasites{ida} ' x ' lfpsites{ilfp}];  
    set(axa{1},'units','pixels')
    axpos=get(axa{1},'position');
    text(axa{1},-30,axpos(4)+30,titletext,'units','pixels','fontweight','bold');
    curgrp=ygrp{yid};
    count=1;
    plotdata=[];    
for itype=1:length(sessiontypes)
    %each da-lfp pair separate figure
    %& each window event type & session type same figure   
    targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);    %'dapos' or 'daneg' types
    %trial ids for this data group for all/phase1/phase4..
    trialtypes=trialinfo(itype); %get trial groups for session
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
        ranj{itrial}=rand(1,length(grptids(itrial).nums))*.15;
    end
    
    cntp=1;
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
        xoff=-0.25+typeoff(itype);
        cmark=mark(itype,:);
    else
        %even
        yyaxis(axa{pid},'right');
        xoff=0.25+typeoff(itype);
        cmark=markr(itype,:);
    end
    set(axa{pid},'ycolor',cmark);
    for ix=1:length(xdata)
        plotdata(count).sesstype=sessiontypes{itype};
        plotdata(count).x=grptids(ix).name;
        plotdata(count).yvar=curgrp{did};
        plotdata(count).yavg=nan;
        plotdata(count).yci=nan;
        plotdata(count).ydata=nan;
        if ~isempty(grpdata(gids{ix}))
            ydata{ix}=grpdata(gids{ix});
            xplot=repmat(ix,1,length(ydata{ix}))+xoff+ranj{ix};
            scatter(axa{pid},xplot,...
                ydata{ix},15,'.','markeredgecolor',cmark,...
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
    if itype==1
        hy=ylabel(axa{pid},curgrp{did},'interpreter','none');
        if ~mod(cntp,2)
            %right axis        
            set(hy,'rotation',270,'units','pixels');
            hypos=get(hy,'position');
            set(hy,'position',[hypos(1)+15,hypos(2),hypos(3)]);  
        end
    end
    if itype==length(sessiontypes)       
        targys=find(contains({plotdata.yvar},curgrp{did})==1 ...
         & contains({plotdata.x},'all')==1);
        getys=[plotdata(targys).ydata];    %all data for curren setup
        ystd=nanstd(getys);
        ymean=nanmean(getys);
        if ~isnan(ystd) && ~isnan(ymean) && ~isempty(ymean) && ~isempty(ystd) 
            if ystd>0
            set(axa{pid},'ylim',[ymean-ystd*2 ymean+ystd*2]);  
            end
        end
    end
    cntp=cntp+1;
    set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)
end
end
end
savedir=[savepath 'allsess_' targdasites{ida} 'x' lfpsites{ilfp} filesep];
if ~isdir(savedir)
    mkdir(savedir);
end
save([savedir ynam{yid} '_' datype '_'  eventtypes{ievent} ...
    '_' targdasites{ida} 'x' lfpsites{ilfp}],'plotdata');
savefig(figsess,[savedir ynam{yid} '_' datype '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}]);
saveas(figsess,[savedir ynam{yid} '_' datype '_'  eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}],'tif')
delete(findall(figsess,'type','text')) 
        
end
end
end
end

end