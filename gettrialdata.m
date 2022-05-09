function data=gettrialdata(chids,selpath,varargin)

argnum=1;
selch=chids;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'selch'
            argnum=argnum+1;
            selch=varargin{argnum};
    end
    argnum=argnum+1;
end
 %chids= param.NCSchannels;
 data={};
%get reconverted path ie XX_pro directory
pathid=strfind(selpath,'_pro');
pathlfp=selpath(1:pathid+4);
dirlfp=dir(pathlfp);
filenums=find(contains({dirlfp.name},['csc' num2str(selch(1)) ])==1);
%get csc signals defined in NCSchannels
if isempty(selch)
    error('no ncs channels defined in param');
end
for it=1:length(filenums)
    tlab=it+100-1;
    for idx=1:length(selch)
        %load in order of appearing in NCSchannels list
        matchcsc=regexpi({dirlfp.name},['csc' num2str(selch(idx)) '_']);
        idscsc=~cellfun('isempty',matchcsc);       %make array of zeros/1 for matching
        matchtrial=regexpi({dirlfp.name},['_' num2str(tlab)]);
        idstrial=~cellfun('isempty',matchtrial);  
        fileid=find((and(idscsc,idstrial))==1);     %get file id that matches both target trial/cscnum
        tempload=load([pathlfp,dirlfp(fileid).name]);
        if it>1 && length(tempload.dg_Nlx2Mat_Samples)~=size(data.lfp{idx},2)
            data.lfp{idx}(it,:)=repmat(nan,1,size(data.lfp{idx},2));
        else
        data.lfp{idx}(it,:)=tempload.dg_Nlx2Mat_Samples;  %store samples, 
        end
        %timestamps in tempload are in microseconds and should match tsLFP already loaded 
        %get/store csc name
        cscnameb=strfind(dirlfp(fileid).name,'_');
        cscname=[dirlfp(fileid).name(1:cscnameb-1) '.ncs'];                
        data.cscNames{idx}=cscname;
    end          
end

end