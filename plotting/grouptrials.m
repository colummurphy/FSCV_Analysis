function trialtypes=grouptrials(goodtrials,trialdata,plotparam, varargin)
%add arg for getting previous trial history
%01/28/2019, all is now everything but switch trials, also delete switch
%trials from all others
%trialdata is trialbytrial(1) created by compiletrialsdir
%exports trialtypes.names and trialtypes.nums
prev=0;
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'prev'
            prev=1;
    end
    argnum=argnum+1;
end
%trialdata.rewardside accurate for what side is big reward for
%calculating switch
for ich=1:length(goodtrials)
    if isfield(plotparam,'sites')        
    trialtypes(ich).site=plotparam.sites(ich);
    else
            trialtypes(ich).site='';
    end
alltrials=1:length(trialdata.alignts);
findswitch=diff(trialdata.rewardside);
findswitch=[findswitch(1) findswitch];
switchtrialsall=find(abs(findswitch)>0);
curtrials=goodtrials{ich};
switchtrials=curtrials(find(ismember(curtrials,switchtrialsall)));
postswitch=switchtrials+1;        %trial after switch
postswitchtrials=intersect(goodtrials{ich},postswitch);

%after big/smal rew
prevbig=find(trialdata.bigside & trialdata.prevside & ~abs(findswitch));   %find previous trials that have same side as current big rew side & not after switch
afterbig=intersect(goodtrials{ich},prevbig);
prevsmall=find(trialdata.bigside & ~trialdata.prevside & ~abs(findswitch));     %find prev side different from current big side & not switch
aftersmall=intersect(goodtrials{ich},prevsmall);

%all good trials now everything but switch
goodtrials2=goodtrials{ich}(~ismember(goodtrials{ich},switchtrials));
goodtrials{ich}=goodtrials2;
afterfailtrials=find(trialdata.prevsuccess==0);
afterfail=goodtrials{ich}(find(ismember(goodtrials{ich},afterfailtrials)));
aftersuccess=goodtrials{ich}(find(~ismember(goodtrials{ich},afterfailtrials)));

nextfailtrials=find(trialdata.nextsuccess==0);
beforesuccess=goodtrials{ich}(find(~ismember(goodtrials{ich},nextfailtrials)));
beforefail=goodtrials{ich}(find(ismember(goodtrials{ich},nextfailtrials)));

rightsidetrials=find(trialdata.rewardside==1);
lsidetrials=find(trialdata.rewardside==0);
goodrightsidetrials=goodtrials{ich}(find(ismember(goodtrials{ich},rightsidetrials)==1));
goodleftsidetrials=goodtrials{ich}(find(ismember(goodtrials{ich},lsidetrials)==1));

quart=floor(length(alltrials)/4);
if quart<10
    phase1trials=[];
    phase4trials=[];
else
    %quarttrials=1:quart:length(alltrials);      %1/28/2019 FIX ALL TRIALS AS IDS NOT JUST GOOD
    %phase1trials=goodtrials(1:quarttrials(2)-1);
    phase1trials=intersect(goodtrials{ich},1:quart+1);
    %phase2trials=goodtrials(quarttrials(2):quarttrials(3)-1);
    %phase3trials=goodtrials(quarttrials(3):quarttrials(4)-1);
    %phase4trials=goodtrials(quarttrials(4):end);
    phase4trials=intersect(goodtrials{ich},length(alltrials)-quart-1:length(alltrials));
   % phase1rtrials=intersect(goodrightsidetrials,phase1trials);
   % phase1ltrials=intersect(goodleftsidetrials,phase1trials);
   % phase4rtrials=intersect(goodrightsidetrials,phase4trials);
    %phase4ltrials=intersect(goodleftsidetrials,phase4trials);
end

%{
trialtypes.names={'all','switch','sleft','sright','aftersuccess','afterfail',...
    'beforesuccess','beforefail','phase1','phase2','phase3','phase4',...
    'srightphase1','sleftphase1','srightphase4','sleftphase4'};
trialtypes.nums={goodtrials{ich},switchtrials,goodleftsidetrials,goodrightsidetrials,...
    aftersuccess,afterfail,beforesuccess,beforefail,phase1trials,phase2trials,phase3trials,phase4trials,...
    phase1rtrials,phase1ltrials,...
    phase4rtrials,phase4ltrials};
%}
trialtypes(ich).names={'all','switch','postswitch','sleft','sright','aftersuccess','afterfail',...
    'phase1','phase4'};

trialtypes(ich).nums={goodtrials{ich},switchtrials,postswitchtrials,goodleftsidetrials,goodrightsidetrials,...
    aftersuccess,afterfail,phase1trials,phase4trials};

if prev
    trialtypes(ich).names={'all','switch','postswitch','sleft','sright','aftersuccess','afterfail',...
    'phase1','phase4','aftersm','afterbig'};
    trialtypes(ich).nums={goodtrials{ich},switchtrials,postswitchtrials,goodleftsidetrials,goodrightsidetrials,...
        aftersuccess,afterfail,phase1trials,phase4trials,aftersmall,afterbig};
end
%fix break trials, only certain events useful
%{
if strcmp(plotparam.trialtype,'fixbreak')
    %just one type of trial, all
    trialtypes={};
    trialtypes.names{1}='all';
    trialtypes.nums{1}=goodtrials{ich};
    nanrts=find(isnan(trialdata.fix_rt)==1);
    numnan=length(nanrts);
    if numnan>2
        %use estimate provided by markers instead
        trialdata.fix_rt=(trialdata.samplesfixeye{1}-trialdata.samplesfix{1})./plotparam.samplespersec;
    end
    %make trial types based on rt's, engaged vs not engaged
    fastrttrials=find(trialdata.fix_rt<=0.2);
    slowrttrials=find(trialdata.fix_rt>0.2);
    trialtypes.names={'all','fastrt','slowrt'};
    trialtypes.nums={goodtrials,fastrttrials,slowrttrials};
end
%}
end
end