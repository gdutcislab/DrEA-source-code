%用于分析
function F = TestHyperVolume5( points, bounds)
    [nrP, dim] = size(points);
    points = repmat(bounds,nrP,1) - points;
    points(points<0)=0;
    %%均衡化
    Aa = prod(max(points));
    points = points ./ repmat(max(points),nrP,1);
    Bb = prod(max(points));
    %%end均衡化
    F = (Recursion(points,1) * Aa) / Bb;
end

function F  = Recursion(points,Iter)
    global tRP5;
    if Iter > numel(tRP5)
        tRP5 = [tRP5 0];
    end
    [nrP, dim] = size(points);
    tRP5(Iter) = tRP5(Iter) + nrP;
    

%     disp(Iter);
    if nrP <= 3
        if nrP == 2
            Min = min(points);
            F = sum(prod(points,2)) - prod(Min);
        elseif nrP == 1
            F = sum(prod(points));
        else
            Min = min(points);
            M1 = min(points([1,2],:));M2 = min(points([1,3],:));M3 = min(points([2,3],:));
            F = sum(prod(points,2)) - prod(M1) - prod(M2) - prod(M3) + prod(Min);
        end
    else
        %方法1，最大最小值，有点接近，但是不是最好
%         [~ , MIndex] = max(min(points,[],2)); 
        

        %原始最小分割点，原始最小分割点,最小方差点
%         Num_DD = zeros(1,nrP);
%         for i  = 1:nrP
%            temp_DD = repmat(points(i,:) ,nrP,1 ) < points;
% %            Num_DD(i) = std(sum(temp_DD));      %最小方差点
%            Num_DD(i) = sum(sum(temp_DD));      %最小分割点
%         end
% 
% 
%         [~,MIndex] = min(Num_DD);

        

        %改进1：
        %考虑两个因素：最小分割点，和最小方差点
        Num_DV = zeros(nrP,1);Num_DD = zeros(nrP,1);
        for i  = 1:nrP
           temp_DD = repmat(points(i,:) ,nrP,1 ) < points;
           aa = sum(temp_DD);
           Num_DV(i) = sum(aa);
%            Num_DD(i) = std(aa);
        end
        [t,MIndex_T] = min(Num_DV);
%         [aa,MIndex_T] = Dominated_Ref([Num_DV,Num_DD],nrP);
%         if numel(MIndex_T)~= 1
%             %对于所有的非支配点，应该选最小方差点？还是说选最小分割点？
%             [~,MIndex] = min(aa(:,2));       %先试着选最小方差点
%             MIndex = MIndex_T(MIndex);
%         else
            MIndex = MIndex_T;
%         end
       if Iter == 1
            %AA = sum(Num_DD);
            BB = Num_DD(MIndex);
%             disp(BB);
        end


%        %改进2：最小分割数最大的作做ref
%        Num_DD = zeros(1,nrP);
%         for i  = 1:nrP
%            temp_DD = repmat(points(i,:) ,nrP,1 ) < points;
%            aa = sum(temp_DD);
%            Num_DD(i) = max(sum(temp_DD));
%         end
%         [~,MIndex] = min(Num_DD);
        
        
        %改进3:
%         temp = 1:nrP;
%         Num_DV = zeros(nrP,1);Num_DD = zeros(nrP,1);
%         for i  = 1:nrP
%            temp_DD = repmat(points(i,:) ,nrP,1 ) < points;
%            aa = sum(temp_DD);
%            Num_DV(i) = sum(aa);
%            Num_DD(i) = sum(aa > 0);
%         end
%         %aa = max(Num_DD);
%         aa = find(Num_DD - max(Num_DD) == 0);
%         Num_DV = Num_DV(aa);Num_DD = Num_DD(aa);
%         temp = temp(aa); 
%         
%         [~,MIndex] = min(Num_DV);
%         MIndex = temp(MIndex);
        
        
        
        
       %clear Num_D Aa temp_D;
       RP = points(MIndex,:);

        SubPoint = [];
        tempPoint = points;
 
        for i = 1:dim
            Aa = find(tempPoint(:,i) - RP(i) > 0);
            if ~isempty(Aa)
                DimSubPoint = tempPoint(Aa,:);
                DimSubPoint(:,i) =  DimSubPoint(:,i) - RP(i);
                SubPoint = [SubPoint cell(1,1)];
                SubPoint{end} = [SubPoint{end}; DimSubPoint];
                tempPoint(Aa,i) = RP(i);
            end
        end
        F = prod(RP);
        Lenght = length(SubPoint);
%         clear tempPoint RP;
        
        for i = 1 :Lenght
          F = F + Recursion(SubPoint{i},Iter+1);
        end
    end
end

%对所有参照点做的非支配排序
function [points,PIndex] = Dominated_Ref(points,nrP)
    %第一维比较最小，第二维比较最大
    dominated = [];PIndex = 1:nrP;
    for i=1:nrP - 1
        for j = i + 1:nrP
            if sum( (points(i,:) - points(j,:) ) >= 0, 2) == 2
                dominated = [dominated i];
            elseif sum( (points(i,:) - points(j,:) ) <= 0, 2) == 2
                dominated = [dominated j];
            end
        end
    end
    points(dominated,:) = [];
    PIndex(dominated) = [];
end
