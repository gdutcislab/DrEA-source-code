function p=DTLZ1(p)
 p.name='DTLZ1';
 p.od = 3;
 p.pd = 14;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %DTLZ2  evaluation function.
    function y = evaluate(x)
%            y = zeros(3,1);
%         g =100*( (length(x)-2) +  sum((x(3:end)-0.5).^2-cos(20*pi*(x(3:end)-0.5))));
%         y(1) = x(1)*x(2)*(1+g);
%         y(2) = x(1)*(1-x(2))*(1+g);
%         y(3) = (1-x(1))*(1+g);
        M=3;
        [row,col] = size(x);
        if row>col
            x = x';
        end
        Xm = x(M:end);
        gXm =100*(length(Xm)+ sum((Xm-0.5).^2-cos(20*pi*(Xm-0.5))));
        cofficentMatrix = tril(ones(M,M),-1);
        X = repmat(x(1:M),M,1);
        X = X.*cofficentMatrix;
        for i = 1:M-1
            X(i,i) = 1-x(i);
        end
        X(M,M) = 1;
        X = X+(tril(ones(M,M),-1))';
        y = cumprod(X,2);
        y = y(:,end)*(1+gXm)/2;
    end
end