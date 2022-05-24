function trids=getTrialIDs(trlists,ttypes,signaltype,varargin)
%Get trialids given the trlists structure containing all trial information
%and variables, as well as signaltype = FSCV or Nlx to find the "good
%trials" excel sheet for FSCV channels or agglomeration of good trials from
%FSCV for ephys
%ttypes = Condition types first arg must be main condition (big, small, targetbreak,fixbreak),
% the rest can be target side (left or right) {'big'},{'left'},{ 
%Other categories:
% 'nosw' - remove trials after block switch/reversal
%   {'post','fixbreak'} - trials after fixbreak
%   {'post','fixbreak','fixbreak'}  - after 2 fixbreaks
%   {'post','reward'} - Trials only after any reward


dach=0;
argnum=1;
while length(varargin)>argnum
    switch varargin{argnum}
        case 'dach'
            argnum=argnum+1;
            dach=varargin{argnum};
    end
    argnum=argnum+1;
end

%Get main condition trials (e.g. big/smal or err)
trids=[];
for ii=1:length(ttypes{1})
    %May be big and small or err combined
        trids=[trids find(contains({trlists.trlist.type},ttypes{1}{ii}))];
end
trids=sort(trids);

side=[];    %Target side condition
preTypes={};%Trial types preceding current trial
numtypes=2;
while length(ttypes)>=numtypes
    if ~iscell(ttypes{numtypes})
        %Just string provided, e.g. 'left' or 'right'
        %{
        switch ttypes{numtypes}
            case 'left'
                side=0; %Leftside trials
            case 'right'
                side=1;
        end
        %}
    else
        %It is a cell containing multiple strings, e.g.
        %{'post','reward'}
        if any(contains(ttypes{numtypes},'left'))
            side=0; %Leftside trials
        end
          if any(contains(ttypes{numtypes},'right'))
            side=1; %Leftside trials
        end      
        if any(contains(ttypes{numtypes},'post'))
            %After certain type of trials only
            posttypes=1;
            for pt=2:length(ttypes{numtypes})
                preTypes{pt-1}=ttypes{numtypes}{pt};    %Get pre-trial conditions
            end
        end
    end
    numtypes=numtypes+1;    %Go to next provided trial condition
end


%Retain only good/valid trials (no artifact) depending on signaltype
if strcmp(signaltype,'fscv')
    goodtrs=find([trlists.fscv{dach}.good]);
    trids=intersect(trids,goodtrs);%Good trials only
else
    %Any other data variable take good trials from fscv (or of all chs)
    %Get union of all chs
    goodtrslist=[];
    for ii=1:length(trlists.fscv)
        goodtrslist(ii,:)=[trlists.fscv{ii}.good];
    end
    goodtrs2=sum(goodtrslist);
    goodtrs=find(goodtrs2~=0);
    trids=intersect(trids,goodtrs);
end

%Get target side condition, if provided as one of the types
if ~isempty(side)
    trids=intersect(trids,find([trlists.trlist.side]==side));
end
%Get pre-trial conditions
if ~isempty(preTypes)
    seltrids=[];
    trids2=trids;
    for n=1:length(preTypes)
        pretrids=trids-n;
        pretrids=pretrids(pretrids>0 & pretrids<length(trlists.trlist));
        pretypes={trlists.trlist(pretrids).type};   %Find nth previous trial type
        pretypeids=find(contains(pretypes,preTypes{n}));    %Get pretrids #
        seltrids=pretrids(pretypeids);  %Get pretrids for those containing indicated pre-trial event at nth pretrial 
        seltrids=seltrids+n;    %Get actual trid for current trial post the pre-trial event
        seltrids=intersect(seltrids,trids);%Only get trids that intersect with current good list--THIS SHOULD STAY SAME since preselected from trids list.
        trids=seltrids;
    end
end