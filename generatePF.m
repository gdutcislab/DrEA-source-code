function generatePF(name,M)
         input = 10000;
         switch name
             case{'WFG1'}
                 h = UniformPoint(input,M);
                 for i = 1 : size(h,1)
                     c = ones(1,M);
                     k = find(h(i,:)~=0,1);
                     for j = k+1 : M
                         temp     = h(i,j)/h(i,k)*prod(1-c(M-j+2:M-k));
                         c(M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
                     end
                     for j = 1 : M
                         h(i,j) = prod(1-c(1:M-j)).*(1-sqrt(1-c(M-j+1)^2));
                     end
                     temp   = acos(c(1))*2/pi;
                     h(i,M) = 1 - temp - cos(10*pi*temp+pi/2)/10/pi;
                 end
                 h = repmat(2:2:2*M,size(h,1),1).*h;
             case{'WFG2'}
                 x = ReplicatePoint(input,M-1);
                 x(:,1) = sqrt(x(:,1));
                 x = [x,zeros(size(x,1),1)];
                 h = convex(x);
                 h(:,M) = disc(x);
                 h = repmat(2:2:2*M,size(h,1),1).*h;
                 h = h(NDSort(h,1)==1,:);
             case{'WFG3'}
                 PopDec    = (0:1/(input-1):1)';
                 PopDec    = [PopDec,zeros(input,M-2)+0.5,zeros(input,1)];
                 h         = linear(PopDec);
                 h         = repmat(2:2:2*M,size(h,1),1).*h;
             case{'WFG4'}
                 h = UniformPoint(input,M);
                 h = h./repmat(sqrt(sum(h.^2,2)),1,M);
                 h = repmat(2:2:2*M,size(h,1),1).*h;
             case{'WFG5','WFG6','WFG7','WFG8'}
                 h = UniformPoint(input,M);
                 h = h./repmat(sqrt(sum(h.^2,2)),1,M);
                 h = repmat(2:2:2*M,size(h,1),1).*h;
             case{'WFG9'}
                 h = UniformPoint(input,M);
                 h = h./repmat(sqrt(sum(h.^2,2)),1,M);
                 h = repmat(2:2:2*M,size(h,1),1).*h;
         end       
        file      = h';
        filename  = strcat('pf_',name,'-',num2str(M)); 
        filetype  = 'dat';
        folder    = 'PFStar';
        mySave(file,filename,filetype,folder)
end

function mySave(file,filename,filetype,folder)
    name      = strcat(filename,'.',filetype);
    if ( strcmpi(filetype, 'dat'))
        fid       = fopen(name,'w');
        [numline,numi]   = size(file);
        for i=1:numline
            for j=1:numi
                fprintf(fid,'%8.6f ',file(i,j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
    elseif  ( strcmpi(filetype, 'fig'))
        saveas(file,name);
    end
    judge     = exist(folder);
    if judge ~= 7
        system(['mkdir ', folder]);
    end
    file_path = strcat(cd,'\',folder);
    movefile(name,file_path); 
end

function [W,N] = UniformPoint(N,M)
%UniformPoint - Generate a set of uniformly distributed points on the unit
%hyperplane
%
%   [W,N] = UniformPoint(N,M) returns approximate N uniformly distributed
%   points with M objectives.
%
%   Example:
%       [W,N] = UniformPoint(275,10)

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group
% You are free to use the PlatEMO for research purposes. All publications
% which use this platform or any code in the platform should acknowledge
% the use of "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and
% Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, IEEE Computational Intelligence Magazine, 2017, in press".
%--------------------------------------------------------------------------

    H1 = 1;
    while nchoosek(H1+M,M-1) <= N
        H1 = H1 + 1;
    end
    W = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
    W = ([W,zeros(size(W,1),1)+H1]-[zeros(size(W,1),1),W])/H1;
    if H1 < M
        H2 = 0;
        while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
            H2 = H2 + 1;
        end
        if H2 > 0
            W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
            W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
            W  = [W;W2/2+1/(2*M)];
        end
    end
    W = max(W,1e-6);
    N = size(W,1);
end
function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end
function Output = disc(x)
    Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end
function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end