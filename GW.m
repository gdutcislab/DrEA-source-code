function weight = GW(val,num)
     [dim,siz]         = size(val);
     noval             = val./repmat(sqrt(sum(val.^2)),[dim,1]);
     dis               = noval'*noval;
     weight            = zeros(dim,num);
     loc               = floor(rand*siz)+1;     
     pdis              = dis(loc,:);
     [ser,loc]         = min(pdis);
     weight(:,1)       = noval(:,loc);
    for i=2:num
        [ser,loc]      = min(pdis);
        weight(:,i)    = noval(:,loc);
        pdis           = max([pdis;dis(loc,:)]);
    end
end