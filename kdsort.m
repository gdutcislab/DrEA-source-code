function pop = kdsort(pop,params,state)   
    N = length(pop);    %popsize
    for i = 1:N
        pop(i).rank = 0;
        pop(i).distance = 0;
    end    
    
    domMat      = calcDominationMatrix(pop,params,state); % domination matrix for efficiency    
    [pop,front] = calcRank(pop,domMat);                  % 计算个体的秩和每一层的个体指标

    %*************************************************************************
    % 3. Calculate the distance
    %*************************************************************************

%     pop = calcCrowdingDistance( pop, front);
end

function domMat = calcDominationMatrix(pop,params,state)    
    obj     = [pop(:).objective];
    [dim,N] = size(obj); 
    obj    = (obj-repmat(params.idealpoint,[1,N]))./(repmat(params.nadirpoint-params.idealpoint,[1,N]));
    domMat  = zeros(N, N);
    t  = state.currentGen/(params.iteration/params.popsize);
    inter = pi*((1- 1/(1+exp(-10*t))));
%     inter = pi/2*((1- sqrt(t^3-3*t^2+3*t)));
    temp = cot(inter);
    for p = 1:N-1
        for q = p+1:N
            pdomq = false;
            qdomp = false;
            tt     = sqrt(sum((obj(:,p)-obj(:,q)).^2)-(pop(p).rdist-pop(q).rdist)^2);
%              tt     = sqrt(sum((obj(:,p)-obj(:,q)).^2)-(sum(obj(:,p)-obj(:,q)))^2/dim);
%             for i = 1:numObj
%                 if( obj(i,p) < obj(i,q) )         % objective function is minimization!
%                     pdomq = true;
%                 elseif(obj(i,p) > obj(i,q))
%                     qdomp = true;
%                 end
%             end
%             if  sum(obj(:,p))/sqrt(dim)+0.3*tt<sum(obj(:,q))/sqrt(dim)
%                 pdomq = true;
%             elseif sum(obj(:,q))/sqrt(dim)+0.3*tt<sum(obj(:,p))/sqrt(dim)
%                 qdomp = true;
%             end 

            if  pop(p).rdist + temp*tt < pop(q).rdist
                pdomq = true;
            elseif pop(q).rdist+ temp*tt <pop(p).rdist
                qdomp = true;
            end 
            if( pdomq && ~qdomp )
                domMat(p, q) = 1;
            elseif(~pdomq && qdomp )
                domMat(p, q) = -1;
            end
        end
    end
    domMat = domMat - domMat';
end

function [pop,front] = calcRank(pop,domMat)
    N = length(pop);    %popsize
    ind = repmat(struct('np',0, 'sp', []),[1,N]);    
    
    % Compute np and sp of each indivudal
    for p = 1:N-1
        for q = p+1:N
            if(domMat(p, q) == 1)          % p dominate q
                ind(q).np = ind(q).np + 1;
                ind(p).sp = [ind(p).sp , q];
            elseif(domMat(p, q) == -1)     % q dominate p
                ind(p).np = ind(p).np + 1;
                ind(q).sp = [ind(q).sp , p];
            end
        end
    end

    % The first front(rank = 1)
    front(1).f = [];    % There are only one field 'f' in structure 'front'.
                        % This is intentional because the number of individuals
                        % in the front is difference.
    for i = 1:N
        if( ind(i).np == 0 )
            pop(i).rank = 1;
            front(1).f = [front(1).f, i];
        end
    end  
    
    % Calculate pareto rank of each individuals, viz., pop(:).rank 
    fid = 1;        %pareto front ID
    while( ~isempty(front(fid).f) )
        Q = [];
        for p = front(fid).f
            for q = ind(p).sp
                ind(q).np = ind(q).np -1;
                if( ind(q).np == 0 )
                    pop(q).rank = fid+1;
                    Q = [Q, q];
                end
            end
        end
        fid = fid + 1;
        front(fid).f = Q;
    end
    front(fid) = [];    % delete the last empty front set
end


function pop = calcCrowdingDistance(pop, front)
    numObj = length( pop(1).objective );  % number of objectives
    for fid = 1:length(front)
        idx = front(fid).f;
        frontPop = pop(idx);        % frontPop : individuals in front fid

        numInd = length(idx);       % nInd : number of individuals in current front

        obj = [frontPop.objective]';
        obj = [obj, idx'];          % objctive values are sorted with individual ID
        for m = 1:numObj
            obj = sortrows(obj, m);

            colIdx = numObj+1;
            pop( obj(1, colIdx) ).distance = Inf;         % the first one
            pop( obj(numInd, colIdx) ).distance = Inf;    % the last one

            minobj = obj(1, m);         % the maximum of objective m
            maxobj = obj(numInd, m);    % the minimum of objective m

            for i = 2:(numInd-1)
                id = obj(i, colIdx);
                pop(id).distance = pop(id).distance + (obj(i+1, m) - obj(i-1, m)) / (maxobj - minobj);
            end
        end
    end
end
