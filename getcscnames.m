function csclabels = getcscnames(targets,varargin)
%convert map (varargin) with string array of 'csc#' 'name' repeated. 
%for those targeted in another array target channels targets
csc_map={};
targetids={};
csclabels={};
argnum=1;
while argnum <= length(varargin)
    switch argnum
        case 1
            csc_map=varargin{argnum}; %normalize signals (no extra arg needed)
    end
    argnum = argnum + 1;
end
%convert double target array to string array
for ii=1:length(targets)
    targetids{ii}=num2str(targets(ii));
    csclabels{ii}=['csc' targetids{ii}];
    if ~isempty(csc_map)
        match=strcmp(csc_map,targetids{ii});
        id=find(match==1);
        if isempty(id)
            warning(['missing label for csc' targetids{ii}]);
        else
            csclabels{ii}=csc_map{id+1};        %get label, cell after identifier
        end
    end
end
end
