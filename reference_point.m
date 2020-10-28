function point = reference_point(mop)
objDim = mop.od;
format long             
% p = 3 or 2
% n = objDim+p-1，
% sample_size = C(n,p) 排列组合数

switch objDim
    case 3
        p1 = 6;            
        n1 = objDim+p1-1;       
        sample_size1 = nchoosek(n1,p1);      
        point  = initweight(objDim, sample_size1); 
    case 5
        p1 = 6;            
        n1 = objDim+p1-1;       
        sample_size1 = nchoosek(n1,p1);      
        point  = initweight(objDim, sample_size1); 
%         point = [point,ones(objDim,1)/objDim,(point+1/objDim)./2];
%         p1 = 2;
%         layer1 = 1;
%         n1 = objDim+p1-1;
%         sample_size1 = nchoosek(n1,p1);
%         point = samplingIGD(objDim, sample_size1, layer1);
% %         point =eye(objDim)+10e-5-diag(10e-5*ones(1,objDim));
%         point = [point;ones(1,objDim)/objDim;(point+1/objDim)./2];
    case 8
        p1 = 3;
        p2 = 2;
        n1 = objDim+p1-1;
        n2 = objDim+p2-1;
        sample_size1 = nchoosek(n1,p1);
        sample_size2 = nchoosek(n2,p2);
        W1 = initweight(objDim, sample_size1);
        W2 = initweight(objDim, sample_size2);
        W2= (1 - 0.5) / objDim * ones(objDim, sample_size2) + 0.5 * W2;
        point  = [W1,W2];
    case 10
        p1 = 3;
        p2 = 2;
        n1 = objDim+p1-1;
        n2 = objDim+p2-1;
        sample_size1 = nchoosek(n1,p1);
        sample_size2 = nchoosek(n2,p2);
        W1 = initweight(objDim, sample_size1);
        W2 = initweight(objDim, sample_size2);
        W2= (1 - 0.5) / objDim * ones(objDim, sample_size2) + 0.5 * W2;
        point  = [W1,W2];
    case 15
        p1 = 2;
        p2 = 3;
        n1 = objDim+p1-1;
        n2 = objDim+p2-1;
        sample_size1 = nchoosek(n1,p1);
        sample_size2 = nchoosek(n2,p2);
        W1 = initweight(objDim, sample_size1);
        W2 = initweight(objDim, sample_size2);
        W2= (1 - 0.5) / objDim * ones(objDim, sample_size2) + 0.5 * W2;
        point  = [W1,W2];
    otherwise 
        error('Undefined test problem name');                
end  

end

% This function is written by Dr. Aimin Zhou for generating any number of weight vectors

function W = initweight(objDim, N)

    U = floor(N^(1/(objDim-1)))-2;
    M = 0;
    while M<N
        U = U+1;
        M = noweight(U, 0, objDim); 
    end

    W      = zeros(objDim, M);
    C      = 0;
    V      = zeros(objDim, 1);
    [W, C] = setweight(W, C, V, U, 0, objDim, objDim);
    W      = W / (U + 0.0);

%     pos     = (W < 1.0E-5);
%     W(pos)  = 1.0E-5;

end

%%
function M = noweight(unit, sum, dim)

    M = 0;

    if dim == 1
        M = 1; 
        return;
    end

    for i = 0 : 1 : (unit - sum)
        M = M + noweight(unit, sum + i, dim - 1);
    end

end

%%
function [w, c] = setweight(w, c, v, unit, sum, objdim, dim)

    if dim == objdim
        v = zeros(objdim, 1);
    end

    if dim == 1
        c       = c + 1;
        v(1)    = unit - sum;
        w(:, c)  = v;
        return;
    end

    for i = 0 : 1 : (unit - sum)
        v(dim)  = i;
        [w, c]  = setweight(w, c, v, unit, sum + i, objdim, dim - 1);
    end

end

%% This function is used to sample points on the PF of DTLZ instance

function P = samplingIGD(objDim, sample_size, layer)
% layer is 1 / 2
    % generate the reference vectors
    if objDim < 6
        W = initweight(objDim, sample_size);
        W = W';
    else
%         W = extract_weight(objDim, no_gaps1, no_gaps2);
%         W = extract_weight(objDim);
        W = initweight(objDim, sample_size);
        W = W';
    end  
   
        denominator = sum(W, 2);
        deMatrix = denominator(:, ones(objDim, 1));
        P = W ./ (layer*deMatrix);   
    
end

function x = dominated_relationship(a, b, m)
% 支配关系判断
% 入口参数：a -个体1，b -个体2，m -目标维数，
% 出口参数：x -1,a支配b; 2,b支配a; 3,ab相等; 4,ab互不支配;
    t = 0;
    q = 0;
    p = 0;

    e = 0.00001;
    for i = 1 : m
        if a(i) <= b(i)
            t = t + 1;
        end
        if  a(i) >= b(i)
            q = q + 1;
        end
        if  a(i) == b(i)
            p = p + 1;
        end
    end
    % for i=1:m
    %     if fix(a(1,c+i)/e)<=fix(b(1,c+i)/e)
    %         t=t+1;
    %     end
    %     if  fix(a(1,c+i)/e)>= fix(b(1,c+i)/e)
    %         q=q+1;
    %     end
    %     if  fix(a(1,c+i)/e)== fix(b(1,c+i)/e)
    %         p=p+1;
    %     end
    % end
    if t == m & p ~= m
        x = 1;
    elseif q == m & p ~= m
        x = 2;
    elseif p == m
        x = 3;
    else
        x = 4;
    end
end

function NPOP = find_nondominated(POP, m)
%寻找非支配解
%入口参数：POP-种群,c-变量个数,m-目标个数
%出口参数：NPOP-非支配解集

i = 1;
while i <= size(POP, 1)
    flag = 0;
    j = i + 1;
    while j <= size(POP, 1)
        x = dominated_relationship(POP(i, :), POP(j, :), m);
        if x == 2
            flag = 1;
            break;
        elseif x == 3
            POP(j, :) = [];
        elseif x == 1
            POP(j, :) = [];
        else
            j = j + 1;
        end
    end
    if flag == 1
        POP(i, :) = [];
    else
        i = i + 1;
    end
end
NPOP = POP;
end
function val_w = Weight(popsize,objDim)
    if objDim==2
        start      = 1/(popsize*100000);
        val_w(1,:) = linspace(start,1-start,popsize);
        val_w(2,:) = ones(1,popsize)-val_w(1,:);
    elseif objDim==3
        val_w      = assign3(popsize,1);
    else
        val_w = lhsdesign(popsize, objDim, 'criterion','maximin', 'iterations', 1000)';
    end
end
function val_w = assign3(popsize,sideL) 
     val_w          = zeros(3,popsize);   
     [numline,numpoint]=comline2(popsize);
     start      = sideL/(numline*10);
     val_one    = linspace(start,sideL-start,numline);
     reg        = 1;
     for i=1:numline
         nump   = numpoint(numline-i+1);
         if nump>1
             val    = zeros(3,nump);
             Lline  = sideL-val_one(i);
             lstart = Lline/(10*nump);
             val(1,:)=linspace(lstart,Lline-lstart,nump);
             val(2,:)=Lline-val(1,:);
             val(3,:)=val_one(i);
         else
             val    = zeros(3,nump);
             val(1)=(sideL-val_one(i))/2;
             val(2)=(sideL-val_one(i))/2;
             val(3)=val_one(i);
         end
         val_w(:,reg:reg+nump-1)=val;
         reg=reg+nump;
     end
end

function [numline,numpoint]=comline2(popsize)
    numline    = round((sqrt(popsize*8+1)-1)/2);
    numpoint   = 1:numline;
    vp         = (numline+1)*numline/2;
    disnum     = popsize-vp;
    numpoint(numline-abs(disnum)+1:numline)=numpoint(numline-abs(disnum)+1:numline)+sign(disnum)*ones(1,abs(disnum));  
end