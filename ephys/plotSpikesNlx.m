function plotSpikesNlx(filename,path,unit,ev,win,binw,cont)
%fromp lx matlab output:
%ch
%unit
    %TS
    %PC1
    %PC2
    %PC3
%plot spikes aligned to events from neuralynx
%data imported from plexon offline sorter as exported in mat format using
%export-per-waveform data with ch, unit, timestamp, pc1, pc2, pc3
%indicate filename (e.g. csc12.mat) and path
%indicate unit in argument (i.e. 1 is usually A,, 0 unsorted, 2 b, etc.)
%event to align to (ev), e.g. 'display_target' for display target
%win [-4 8] around event
%binw size, default 0.01 s
%cont is contingency e.g. BR eventcode 45 within +4.5 s of ev,
%{EVcode,4.5},cont={'big',7}; or big reward within 7 s
plotmean=0; %plot everything =1
markers={'display_fix','start_target','reward_big','reward_small'};
event_codes={
    '4'     'display_fix' ...
    '5'     'start_fix' ...
    '6'     'break_fix' ...
    '10'    'display_target' ...
    '11'    'start_target'  ...
    '12'    'break_target'  ...
    '14'    'error' ...
    '29'    'left_condition'    ...
    '30'    'right_condition'   ...
    '45'    'reward_big'    ...
    '46'    'reward_small'  ...
    };

figpos=[50,50,800,900];
figsess=figure('position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
ax={};
axsiz=[350 300];
ax.raster=subplot(2,1,1,'parent',figsess);
ylabel(ax.raster,'Trial #'); xlabel(ax.raster,'Time (s) Relative to Evt at 0 s');
hold(ax.raster,'on');
ax.hist=subplot(2,1,2,'parent',figsess);
ylabel(ax.hist,'Spike Rate (S/s)'); xlabel(ax.hist,'Time (s) Relative to Evt at 0 s');
hold(ax.hist,'on');
cla(ax.hist);
cla(ax.raster);

pathparts=split(path,filesep);
targfolder=pathparts{find(contains(pathparts,'ephys'))+1};
titlename={['Unit ' num2str(unit) ' | Event = ' ev ' | Bin Width = ' num2str(binw) ' s']; [filename]; targfolder};
title(ax.raster,titlename,'interpreter','none')


%get events from neuralynx events file
direvents=dir([path '*.nev']);
fileevents=[path direvents.name];       %get events file name
%load neuralynx events file (mat converted from dg_nlx2mat.m
%load dg_Nlx2Mat_EventStrings, dg_Nlx2Mat_Timestamps, dg_Nlx2Mat_TTL
dg_Nlx2Mat(fileevents); 
load([[fileevents(1:end-4)] '.mat']);
nlx_ts=dg_Nlx2Mat_Timestamps.*1e-6;      %convert to sec (units for sorted units)

%get alignment event code
evid=find(contains(event_codes,ev))-1;
evid=str2num(event_codes{evid});

%get ts for alignment event code
aln_ids=find(dg_Nlx2Mat_TTL==evid);
aln_ts=nlx_ts(aln_ids);     %row vector of ts for targ ev


%get spike times for indicated unit
spikes=load([path filename]);
aname=fieldnames(spikes);
spk=getfield(spikes,aname{1});
spk_ids=find(spk(:,2)==unit);
spk_ts=spk(spk_ids,3);          %timestamp in column 3 if exported as indicated above in header
spk_wf=spk(spk_ids,7:end);      %wf in column 7+ if exported as indicated above in header
%figure; plot(spk(spk_ids,7:end)')
%figure; plot(nanmean(spk(spk_ids,7:end))')

%check plot of iniital spike times and task events
%figure;
%plot(spk_ts(1:300),repmat(1,1,300),'.');
%line([aln_ts(1:15); aln_ts(1:15)],repmat([0.8 1.2],15,1)','color',[0 0 0]);

%%%%%bin spikes based on time windows aligned to specified task event, aln_ts

%first organize spikes into matrix where each row is a trial, containing
%all the spikes recorded aligned to the targeted alignment event
spktsreps=repmat(spk_ts',length(aln_ts),1);     %repeat matrix for spike ts' for every alignment event
alnspkts=spktsreps-aln_ts';     %spike ts's now relative to aln ts (each row is different aln event, each column is aspike time)
[spksinwinX, spksinwinY,v]=find(alnspkts>=win(1) & alnspkts <= win(2));
%make matrix of spikes within window
mostspikesforawin=mode(spksinwinX);         %find most frequent values in array to determine largest dimension
nummaxtrialspks=length(find(spksinwinX==mostspikesforawin));
spksinwin=nan(length(aln_ts),nummaxtrialspks);    %this will be spike time stamps relative to aln event, i.e. raster
spksinwinabs=nan(length(aln_ts),nummaxtrialspks); %absolute time stamps
%uniquespksX=unique(spksinwinX); %get unique values, i.e. remove repeats for X (that is aln event)
for alnevid=1:length(aln_ts)
    currtrialspks=find(spksinwinX==alnevid);
    if currtrialspks>0
        spksinwin(alnevid,1:length(currtrialspks))=alnspkts(alnevid,spksinwinY(currtrialspks));
        spksinwinabs(alnevid,1:length(currtrialspks))=spktsreps(alnevid,spksinwinY(currtrialspks));
    end
end
%raster plot for first 5 trials
%figure; plot(spksinwin(1:5,:),repmat([1:5]',1,size(spksinwin(1:5,:),2)),'k.')

%%%%%%%histogram and bin data
edges=win(1):binw:win(2);
%spksinwinvect=reshape(spksinwin,1,[]);      %make a single vector since we're just counting now it does not matter what trial the spikes belong to
%spksvect=spksinwinvect(~isnan(spksinwinvect));      %without nan's padding values left over
[spkN,edgesh,bin] = histcounts(spksinwin,edges);
spkrate=spkN./binw./length(aln_ts);
%[spkN2,edgesh,bin] = histcounts(spksvect,edges);       %no need to vectorize to run hist binning




%figure; plot(edges(1:end-1)+binw/2,spkN,'k.');      %plot averaged binned spike data



%%%%condition indicated contingency (i.e. within big or small reward trial)
%get cont event code (if supplied)
cla(ax.hist);
cla(ax.raster);
contid='';
contts=[];
cont_trials={};         %trial ids, ie aln evt trials with the cont evt
cont_relts={};          %timestamps relative to aln event
otrials=[];
ospkrate=[];
if ~isempty(cont)
    numconts=length(cont)/2;
    for ic=1:numconts
        contid=find(contains(event_codes,cont{2*ic-1}))-1;
        contid=str2num(event_codes{contid});
        %get ts for cont event code
        contids=find(dg_Nlx2Mat_TTL==contid);
        contts{ic}=nlx_ts(contids);     %row vector of ts for targ ev
        alntsrepscont=repmat(aln_ts',1,length(contts{ic}));     %repeat matrix for aln ts' for every cont ts
        relts=contts{ic}-alntsrepscont;        %ts relative to alignment event, each row is still aln evt trial
        [continalnX, continalnY,v]=find(relts<=cont{2} & relts >=0);
        cont_trials{ic}=continalnX;      %trials with contingent event to alignment event
        for id=1:length(cont_trials{ic})
            cont_relts{ic}(id)=relts(cont_trials{ic}(id),continalnY(id));        
        end

        %%%%%%%histogram and bin data with contingency
        [cspkN,~,~] = histcounts(spksinwin(cont_trials{ic},:),edges);
        cspkrate=cspkN./binw./length(cont_trials{ic});
            
        plotcolors=[1 0 0; 0 0 1; 0 1 0];
       
        scatter(ax.hist,edges(1:end-1)+binw/2,cspkrate,10,plotcolors(ic,:),'filled','MarkerFaceAlpha',0.5);  
        cspkplot=plot(ax.hist,edges(1:end-1)+binw/2,cspkrate,'-','LineWidth',1,'Color',[plotcolors(ic,:) .2]);
        cscatterdataX=reshape(spksinwin(cont_trials{ic},:),1,[]);
        cscatterdataY=reshape(repmat([cont_trials{ic}]',1,size(spksinwin,2)),1,[]);
        cscatterrast=scatter(ax.raster,cscatterdataX,cscatterdataY,5,plotcolors(ic,:),'filled','MarkerFaceAlpha',0.2);
    end
    %%%%%%%Other remaining spikes histogram and bin data with contingency
    otrials=find(~ismember(1:length(aln_ts),vertcat(cont_trials{:})));
    [ospkN,~,~] = histcounts(spksinwin(otrials,:),edges);
    ospkrate=cspkN./binw./length(otrials);
end



%plot everything
if plotmean || isempty(cont)
scatterspk=scatter(ax.hist,edges(1:end-1)+binw/2,spkrate,10,[ 0 0 0],'filled','MarkerFaceAlpha',0.5); 
spkplot=plot(ax.hist,edges(1:end-1)+binw/2,spkrate,'-','LineWidth',1,'Color',[0 0 0 .2]);
scatterdataX=reshape(spksinwin(1:end,:),1,[]);
scatterdataY=reshape(repmat([1:length(aln_ts)]',1,size(spksinwin,2)),1,[]);
%plot(ax.raster,spksinwin(1:end,:),repmat([1:length(aln_ts)]',1,size(spksinwin,2)),'k.','markersize',0.2) 
scatterrast=scatter(ax.raster,scatterdataX,scatterdataY,5,[0 0 0],'filled','MarkerFaceAlpha',0.2);
end
%plot everything remaining
if ~isempty(ospkrate)
scatterspk=scatter(ax.hist,edges(1:end-1)+binw/2,ospkrate,10,[ .3 .3 .3],'filled','MarkerFaceAlpha',0.5); 
spkplot=plot(ax.hist,edges(1:end-1)+binw/2,ospkrate,'-','LineWidth',1,'Color',[ .3 .3 .3 .2]);
scatterdataX=reshape(spksinwin(otrials,:),1,[]);
scatterdataY=reshape(repmat(otrials,1,size(spksinwin,2)),1,[]);
%plot(ax.raster,spksinwin(1:end,:),repmat([1:length(aln_ts)]',1,size(spksinwin,2)),'k.','markersize',0.2) 
scatterrast=scatter(ax.raster,scatterdataX,scatterdataY,5,[ .3 .3 .3],'filled','MarkerFaceAlpha',0.2);
end

%plot markers
m_relts={};
m_trials={};
if ~isempty(markers)
    marka=markers(~contains(markers,ev));       %all except current aln evt
    for ic=1:length(marka)
        mid=find(contains(event_codes,marka{ic}))-1;    %get event TTL ID
        mid=str2num(event_codes{mid});
        mids=find(dg_Nlx2Mat_TTL==mid);        %get ts for  event code
        mts{ic}=nlx_ts(mids);     %row vector of ts for targ ev
        alntsrepscont=repmat(aln_ts',1,length(mts{ic}));     %repeat matrix for aln ts' for every cont ts
        relts=mts{ic}-alntsrepscont;        %ts relative to alignment event, each row is still aln evt trial
        [minalnX, minalnY,v]=find(relts<=win(2) & relts >=win(1));    %find those markers within window
        m_trials{ic}=minalnX;      %trials with contingent event to alignment event
        m_relts{ic}=[];
        for id=1:length(m_trials{ic})
            m_relts{ic}(id)=relts(m_trials{ic}(id),minalnY(id));        
        end
        
        %plot markers
        if ~isempty(m_relts{ic})
        markerts=nanmean(m_relts{ic});
        ymarker=get(ax.hist,'ylim'); yscale=ymarker(2)-ymarker(1);
        plot(ax.hist,[markerts markerts],ymarker,'k-','Color',[0 0 0 .3]); 
        text(ax.hist,[markerts],max(ymarker)-yscale*.1,marka{ic},'interpreter','none','rotation',90);
        scatter(ax.raster,m_relts{ic},m_trials{ic},25,[ .3 1 .3],'MarkerEdgeAlpha',0.1);
        end
    end
end

