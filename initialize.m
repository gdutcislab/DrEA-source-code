function [params,mop,pop,state]= initialize(params,mop,pop,state)
    num_nod             = params.popsize;
    inds                = randompoint(mop, num_nod);    
    %% 评价每一个个体
    [inds,state]        = evaluate(inds,mop,state);
    val                 = [inds.objective];
    params.idealpoint   = min([params.idealpoint,val],[],2);  
    params              = updatenadir(val,params);   
    noval               = (val-repmat(params.idealpoint,[1,num_nod]))./(repmat(params.nadirpoint-params.idealpoint,[1,num_nod]));
    %% group   
    [team,inds]         = group(params,pop,noval,inds);
    for i=1:params.num_class
        pop(i).num_ind  = length(team{i});
        pop(i).inter    = inds(team{i});
%         if num_p<=pop(i).num_ind
%             tst                      = floor(num_nod*rand(1,pop(i).num_ind-num_p))+1;
%             selind                   = inds([team{i},tst]);
%             pop(i).inter             = selind;
%         else
%             selind                   = inds(team{i});
%             pop(i).inter             = selection(params,selind,pop(i));
%         end
    end 
end
%% 产生个体
function ind = randompoint(prob, n)
    if (nargin==1)
        n=1;
    end
    randarray = rand(prob.pd, n);
    lowend = prob.domain(:,1);
    span = prob.domain(:,2)-lowend;
    point = randarray.*(span(:,ones(1, n)))+ lowend(:,ones(1,n));
    cellpoints = num2cell(point, 1);
    indiv = get_structure('individual');
    ind = repmat(indiv, [1, n]);
    [ind.parameter] = cellpoints{:};
end
function params=updatenadir(valf,params)
          dim = size(valf,1);                
          W   = diag(ones(1,dim))+10e-5-diag(10e-5*ones(1,dim));
%% FindExtremePoints
          params.nadirpoint   = max(valf,[],2); 
          z            = ones(dim,dim);           
          for j=1:dim             
              z(:,j)=argASF(valf,W(:,j),params.idealpoint,params.nadirpoint);               
          end 
          params.extrempoint=z;         
          z1=z-repmat(params.idealpoint,[1,dim]);          
         if rank(z1)==dim          
             b=ones(dim,1);
             a=z1'\b;
%              a=lsqnonneg(z1',b);
             intercepts=1./a;
             if ~any(intercepts<10e-5)&&~any(intercepts>10e+5)
                 params.nadirpoint=intercepts+params.idealpoint;
             end
         end
 
end
%%
function y=argASF(val_f,weight,idealpoint,nadirpoint)
         num_nod      = size(val_f,2);        
         noval        = abs(val_f-repmat(idealpoint,[1,num_nod]))./repmat(nadirpoint-idealpoint,[1,num_nod]);
         max_evo      = max(noval./(repmat(weight,[1,num_nod])));
         [cfval,cfloc]= min(max_evo);
         y            = val_f(:,cfloc);
end
