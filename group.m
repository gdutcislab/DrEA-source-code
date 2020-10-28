%·Group
function [team,inds] = group(params, pop, nornal_obj, inds)
    [Dim,num_p]       = size(nornal_obj);   
    team              = cell(params.num_class,1);      
    center            = [pop.center];
    dis               = center'*nornal_obj;  
    [minval,minindex] = max(dis,[],1); 
    for i = 1:num_p
        team{minindex(i)}  = [team{minindex(i)},i];
        inds(i).rdist      =  minval(i);
%         cdist      =  sqrt(sum(nornal_obj(:,i).^2)-minval(i)^2);  
        inds(i).cdist      = sqrt(sum((nornal_obj(:,i)-center(:,minindex(i)).*inds(i).rdist).^2));
    end   
end
