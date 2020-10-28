function pop=selection(params,pop,oldpop,state)
     if(nargin < 3)
         num = params.popsize;   
     else
         num = oldpop.num_ind;
     end
     if ( strcmpi(params.selectMethod, 'NS'))
         pop            = ndsort(pop);
         pop            = selectNS(pop,num);     
     elseif ( strcmpi(params.selectMethod, 'WS'))
         % Todo      
     elseif ( strcmpi(params.selectMethod, 'MM'))
         pop = selectMM(pop,params,oldpop);
     elseif ( strcmpi(params.selectMethod, 'KD')) 
         pop            = kdsort(pop,params,state);
         pop            = selectKD(pop,params,num);
     else
         error( 'No defined the select method');
     end
end

function pop = KDistance(pop,params,num)
         obj    = [pop.objective];
         cdist  = [pop.cdist];
         rdist  = [pop.rdist];
         [dim,n]= size(obj);        
         obj    = (obj-repmat(params.idealpoint,[1,n]))./(repmat(params.nadirpoint-params.idealpoint,[1,n]));
         r      = (max(cdist)-min(cdist))/(num);
         for i = 1:n-1
             for j = i+1:n  
                  tt     = sqrt(sum((obj(:,i)-obj(:,j)).^2)-(rdist(i)-rdist(j))^2);
%                   tt     =  abs(cdist(i)-cdist(j));
%                    d = (sum(abs(obj(:,i)-obj(:,j)).^(1/dim)))^dim;
                 if tt<r
                     sh(i,j) = 1.0/(1+tt);
                     sh(j,i) = 1.0/(1+tt); 
%                  elseif tt<=r&&rdist(i)>rdist(j)
%                      sh(i,j) = 1.5/(1+tt);
%                      sh(j,i) = 0.5/(1+tt);
%                  elseif tt<=r&&rdist(i)==rdist(j)
%                      sh(i,j) = 1.0/(1+tt);
%                      sh(j,i) = sh(i,j);
                 else
                     sh(i,j) = 0;
                     sh(j,i) = 0;
                 end
             end
             sh(i,i) = 0;
         end
         sh(n,n) = 0;
         for i=1:n
             pop(i).distance= sum(sh(i,:));
%              pop(i).distance= cdist(i);
         end
         [val,ind] = min(cdist);
         pop(ind).distance = -Inf;
%          [val,ind]=max(cdist);
%          pop(ind).distance=-Inf;
end

function nextpop = selectKD(pop,params,num)
    nextpop = pop(1:num);    %just for initializing
    rankVector = [pop.rank];
    n = 0;                  % individuals number of next population
    rank = 1;               % current rank number
    idx = find(rankVector==rank);
    numInd = length(idx);   % number of individuals in current front
    while(n + numInd <= num)
        nextpop(n+1:n+numInd) = pop(idx);
        n = n + numInd;
        rank = rank + 1;
        idx = find(rankVector==rank);
        numInd = length(idx);
    end
    if( n < num )
        tpop= pop(idx);
        tpop = KDistance(tpop,params,num-n);
%         tpop = CDistance(tpop);
        distance   = [tpop.distance];
        [val,idx]  = sort(distance);
        nextpop(n+1:num) = tpop(idx(1:num-n));
    end
end

function nextpop = selectMM(pop,params,oldpop)
     val                 = [pop.objective];
     n                   = size(val,2);
     nextpop             = pop(1:oldpop.num_ind);
     noval  = (val-repmat(params.idealpoint,[1,n]))./(repmat(params.nadirpoint-params.idealpoint,[1,n]));
     weight              = oldpop.weight;
    for j=1:oldpop.num_ind
        max_evo          = max(repmat(weight(:,j),[1,n]).*noval);
        [cfval,cfloc]    = min(max_evo);
        nextpop(j)       = pop(cfloc);
        pop(cfloc)       = [];
        noval(:,cfloc)   = [];
        n                = n-1;
    end
end
function y=fitness(pop)
      obj    = [pop.objective];
      cdist  = [pop.cdist];       
      [dim,n]= size(obj);
      dom    = false(n,n);
      for i=1:n
          Q(i).S=0;
      end
      for i=1:n
        for j=i+1:n            
            if sum(obj(:,i)<obj(:,j))>sum(obj(:,i)>obj(:,j))/0.7              
%                 all(obj(:,i)<=obj(:,j)) && any(obj(:,i)<obj(:,j))
                Q(i).S=Q(i).S+1;
                dom(i,j)=true;                
            elseif sum(obj(:,j)<obj(:,i))>sum(obj(:,j)>obj(:,i))/0.7               
%                 all(obj(:,j)<=obj(:,i)) && any(obj(:,j)<obj(:,i))
                Q(j).S=Q(j).S+1;
                dom(j,i)=true;                
            end            
        end
      end
        S=[Q.S];
    for i=1:n
        y(i)=sum(S(dom(:,i)));
    end
end
