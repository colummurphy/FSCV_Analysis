function result=zscorenan(data,dim)
%if dim ==1 then take mean across columns (ie returns values in single row)
%if dim==2 then across rows 
invrc=1;
if dim==1
    invrc=2;
end
meandata=nanmean(data,dim);
stddata=nanstd(data,0,dim);
repp=repmat(meandata,invrc,size(data,dim));
subdata=data-repp;
for ii=1:size(data,invrc)
    result(ii,:)=subdata(ii,:)./stddata(ii);
end

end