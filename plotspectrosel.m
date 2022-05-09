function maxfft=plotspectrosel(sessnum,site,trials,event,plotparam,varargin)
%plot average spectrogram of supplied trial numbers across different trial
%types for a given session #, function called from plotmultiple
%not merged file, get lfps separately based on trial num 100 etc.
fs=1000;        %sampling freq lfp downconverted
argnum=1;
subject='patra';
maxfft={};
win=[-5 5];
tapers=[3 5];      %nw (time-bandiwth rpdocut, tapers)
%ffttapers=[1.8 1];      %nw (time-bandiwth rpdocut, tapers)
fftwin=[.75 0.15];
freq=[5 55];
label='';
sub=0;
trials2={};
trtypes={};
if isstruct(trials)
trtypes=unique({trials.trialgrps.type});
label=['trls_' trials.cat '_' trials.site '_' [trtypes{:}]];
end
condition='all';
xinfos={};
binfos={};
types={};
fixedcond={};
getmaxfft=0;
plotboth=0;
%{
tx=trialnums(1).cat;
trtypes=unique({trialnums(1).trialgrps.type});
if length(trialnums)>1
    for iix=2:length({trialnums.cat})
        tx=[tx '_vs_' trialnums(iix).cat];
    end
end
label=['trls_' tx '_' trialnums(1).site '_' [trtypes{:}]];
%}

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'cleo'
            subject='cleo';
        case 'win'
            argnum=argnum+1;
            win=varargin{argnum};     
        case 'label'
            argnum=argnum+1;
            label2=varargin{argnum};
            label=[label '_' label2];
        case 'freq'
            argnum=argnum+1;
            freq=varargin{argnum};
        case 'fftwin'
            argnum=argnum+1;
            fftwin=varargin{argnum};
        case 'tapers'
            argnum=argnum+1;
            tapers=varargin{argnum};
        case 'markers'
            argnum=argnum+1;
            markers=varargin{argnum};
        case 'condition'
            argnum=argnum+1;
            condition=varargin{argnum};
        case 'sub'
            %subtract spectra from 2 groups of trials, supply additional
            %group of trials here
            argnum=argnum+1;
            trials2=varargin{argnum};
            sub=1;
            label=[label '_sub'];
        case 'xinfos'
            argnum=argnum+1;
            xinfos=varargin{argnum};
        case 'binfos'
            argnum=argnum+1;
            binfos=varargin{argnum};
        case 'types'
            argnum=argnum+1;
            types=varargin{argnum}; %all conditions to include i.e. {'big','small'}
        case 'fixedcond'
            argnum=argnum+1;
            fixedcond=varargin{argnum}; %i.e. only left or right {'left'}
        case 'nomaxfft'
            getmaxfft=0;
        case 'plotboth'
            getmaxfft=1;
            plotboth=1;
    end
    argnum=argnum+1;
end
homed=fullfile('Z:', 'inj-monkey2', 'analysis', subject, filesep);
%2022 PITT
homed=fullfile('Y:', 'data_MIT', 'analysis', 'patra', filesep);
homed2=fullfile('Y:', 'data_MIT', 'analysis', 'cleo', filesep);

pathname=fullfile(homed, ['chronic' num2str(sessnum)], filesep);
[chs,cscmap,paths]=getconfig(subject,sessnum);
ncsid=cscmap{find(contains(cscmap,site))-1};
%get reconverted path ie XX_pro directory to retrieve downconverted lfps

if isempty(trtypes)
    %get trials here
    target={sessnum,site,types,event};
    if ~isempty(fixedcond)
        target={sessnum,site,types,event,'fixedcond',fixedcond};
    end
    trials=getmulttrials(plotparam,xinfos,binfos,target);
    trtypes=unique({trials.trialgrps.type});
    label=['trls_' trials.cat '_' trials.site '_' [trtypes{:}]];
end


trials=trials.trialgrps;
trialinfo=plotparam.trialgrps(contains({plotparam.trialgrps.sessid},num2str(sessnum))).trialinfo;
infonames={'bigreward','smallreward','targetbreak','fixbreak'};

tslfp=[];
lfp=[];
alnts=[];
count=1;
otherts={};
count2=1;
otherts2={};
alnts2=[];
lfp2=[];
for itype=1:length(trtypes)
    trialtype=trtypes{itype};
      trtype=find(contains(infonames,trialtype));
trialtypes=trialinfo(trtype); %get trial groups for session   
ttype=find(contains(trialtypes.names,condition)==1);
typename=trialtypes.names{ttype};
    pathlfp=fullfile(paths{1},'matlab',[trialtype '_pro'],filesep);
    %load trialbytrial data with ts for alignment event
    load([pathname trialtype filesep 'trials_chronic' num2str(sessnum) '_' trialtype '.mat'],'trialbytrial');
    eventts=getfield(trialbytrial(1),['t' event]);
    otherevts={};
    if ~isempty(markers)
        for ix=1:length(markers)
            if ~strcmp(markers{ix},'outcome')
                otherevts{ix}=getfield(trialbytrial(1),['t' markers{ix}]);
            else
                otherevts{ix}=repmat(30,1,length(trialbytrial(1).alignts));
            end
        end
    end
    targtrials=intersect(trials(itype).trials,trialtypes.nums{ttype});
for it=1:length(targtrials)
    tnum=targtrials(it)+99;
    tempload=load([pathlfp 'csc' ncsid '_' num2str(tnum) '.mat']);      %load csc data for specifieid trial
    lfp(count,:)=tempload.dg_Nlx2Mat_Samples;
    if count==1
        tslfp=(tempload.dg_Nlx2Mat_Timestamps-tempload.dg_Nlx2Mat_Timestamps(1)).*1e-6;
    end
    alnts(count)=eventts(targtrials(it));
    if ~isempty(otherevts)
    for ix=1:length(otherevts)
        otherts{ix}(count)=otherevts{ix}(targtrials(it));
    end
    end
    count=count+1;
end
if sub
    targtrials2=intersect(trials2(itype).trials,trialtypes.nums{ttype});

    for it2=1:length(targtrials2)
        tnum=targtrials2(it2)+99;
        tempload=load([pathlfp 'csc' ncsid '_' num2str(tnum) '.mat']);      %load csc data for specifieid trial
        lfp2(count2,:)=tempload.dg_Nlx2Mat_Samples;
        alnts2(count2)=eventts(targtrials2(it2));
        if ~isempty(otherevts)
        for ix=1:length(otherevts)
            otherts2{ix}(count2)=otherevts{ix}(targtrials2(it2));
        end
        end
        count2=count2+1;
    end
end
end
figffts=figure;
set(0,'CurrentFigure',figffts);    %set figure handle to current figure
set(figffts, 'Color', [1 1 1]);
set(figffts,'Position',[300,150,1200,500]);
hax=axes;
hold(hax,'on'); 
cla(hax);
label2=label;
%getavgspec(lfp,fs,tslfp,alnts,hax,'win',win,'freq',freq,'tapers',tapers,'fftwin',fftwin,'pinkbl');
if ~sub
    if getmaxfft
        [x,maxfft]=getavgspec(lfp,fs,tslfp,alnts,hax,'win',win,'freq',freq,'tapers',tapers,'fftwin',fftwin,'marks',otherts,'getmaxfft');
        label=[label '_sub_maxfft'];
    else
        [x,maxfft]=getavgspec(lfp,fs,tslfp,alnts,hax,'win',win,'freq',freq,'tapers',tapers,'fftwin',fftwin,'marks',otherts);
    end
else
getavgspec(lfp,fs,tslfp,alnts,hax,'win',win,'freq',freq,'tapers',tapers,'fftwin',fftwin,'marks',otherts,'sub',lfp2,alnts2,'marks2',otherts2);
end
title(hax,['sessnum | ' num2str(sessnum) ' | lfp | ' site ' | ' ...
    condition ' | tapers ' num2str(tapers(1)) ' / ' num2str(tapers(2)) ' | win/step (s) ' ...
    num2str(fftwin(1)) ' / ' num2str(fftwin(2)) ' | ' label ],'interpreter','none')

savename=['fft_' num2str(sessnum) '_' site '_' event '_' condition];
if ~isempty(label)
    savename=['fft_' num2str(sessnum) '_' site '_' event '_' condition '_' label];
end
if getmaxfft
        savename=['fft_' num2str(sessnum) '_' site '_' event '_' condition '_' label '_sub_maxfft'];
end
savepath=[homed 'fft' filesep];
if ~isdir(savepath)
    mkdir(savepath);
end
saveas(figffts,[savepath savename '.jpg'],'jpeg')
savefig(figffts,[savepath savename])

if plotboth
    %both sub/maxfft and non sub
    %generate plot for nonsub here, do not save data, just fig
    getavgspec(lfp,fs,tslfp,alnts,hax,'win',win,'freq',freq,'tapers',tapers,'fftwin',fftwin,'marks',otherts);
    title(hax,['sessnum | ' num2str(sessnum) ' | lfp | ' site ' | ' ...
        condition ' | tapers ' num2str(tapers(1)) ' / ' num2str(tapers(2)) ' | win/step (s) ' ...
        num2str(fftwin(1)) ' / ' num2str(fftwin(2)) ' | ' label2 ],'interpreter','none')
    savename=['fft_' num2str(sessnum) '_' site '_' event '_' condition];
    if ~isempty(label2)
        savename=['fft_' num2str(sessnum) '_' site '_' event '_' condition '_' label2];
    end
    savepath=[homed 'fft' filesep];
    if ~isdir(savepath)
        mkdir(savepath);
    end
    saveas(figffts,[savepath savename '.jpg'],'jpeg')
    savefig(figffts,[savepath savename])
end
end