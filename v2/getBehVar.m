function behsig=getBehVar(trlists,behtype,trids,varargin)
%Get behavioral signal given trids and type of beh (behtype)
%for time series data, must provide evids which are the samples for each
%trial for beg/end of signal in which to compute variable average or peak
behsig=[];

evids=[];   %Event ids/sample points to calc data, e.g. lick, given trt data
argnum=1;
while length(varargin)>argnum
    switch varargin{argnum}
        case 'evids'
            argnum=argnum+1;
            evids=varargin{argnum};
    end
    argnum=argnum+1;
end

switch behtype
    case 'rt_fix'
        behsig=[trlists.trlist(trids).frt];
    case 'rt_target'
        behsig=[trlists.trlist(trids).trt];
    case 'lick'
        lickch=find(contains(trlists.nlxnames,'lick'));
        lickdata=vertcat(trlists.nlx{lickch}{trids});
        licksig;


end