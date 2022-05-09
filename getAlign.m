function [fscvSyncTTL, nlxSyncID, alignmentIDs]=getAlign(setup)
%get parameters for synchronizing Nlx and FSCV based on TTL input in FSCV,
%nlx output ID to send TTL to FSCV, and targeted nlx vcortex behavior IDs
%for designing trials

alignmentIDs=[];
fscvSyncTTL=241;          %port 2 in tarheel when bit 0 goes high 240 changes to 241
%set up in vcortex to output 3 TTL to Neuralynx for trial start in 1-DR
%and this is what triggers neuralynx to send TTL to FSCV, 45 in Ken's
%system (Yuri recording)
nlxSyncID=3;            %trial start indicated in Neuralynx to align to FSCV 241 TTL in all files
nextTrialMarksEnd=0;            %want trial Index to mark trial end for each data clip or instead if zero use duration below
alt_smallRewardID=46;
alt_bigRewardID=45;
correctID=13;
smallRewardedIDs=[22 23 alt_smallRewardID];
bigRewardedIDs=[18 19 alt_bigRewardID];
errorIDs=[14 12];       %error, target fixation break
fixBreakIDs=[6];        %fixation break initial fix
targetBreakID=[12];         %target break error
omissionIDs=[15];        %omission

switch setup
    case {'bigReward','bigreward'}
        alignmentIDs=bigRewardedIDs;          %chosen one for alignment
    case {'smallReward','smallreward'}
        alignmentIDs=smallRewardedIDs;          %chosen one for alignment
    case 'error'
        alignmentIDs=errorIDs;
    case {'fixBreak','fixbreak'}
        alignmentIDs=fixBreakIDs;
    case {'targetBreak','targetbreak'}
        alignmentIDs=targetBreakID;
    case 'omission'
        alignmentIDs=omissionIDs;
    otherwise
        warning('not specified');
end

end