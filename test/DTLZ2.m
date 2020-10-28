function p=DTLZ2(p)
 p.name='DTLZ2';
 p.od = 3;
 p.pd = 12;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %DTLZ2  evaluation function.
    function y = evaluate(x)
%          temp         = (2*x(3:p.pd,:)-1).^2;
%          temp1        = sum(temp);
% %          y(1,:)       = (1+temp1).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2).*cos(x(3,:)*pi/2).*cos(x(4,:)*pi/2);
% %          y(2,:)       = (1+temp1).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2).*cos(x(3,:)*pi/2).*sin(x(4,:)*pi/2);
% %          y(3,:)       = (1+temp1).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2).*sin(x(3,:)*pi/2);
% %          y(4,:)       = (1+temp1).*cos(x(1,:)*pi/2).*sin(x(2,:)*pi/2);
% %          y(5,:)       = (1+temp1).*sin(x(1,:)*pi/2);
%          y(1,:)       = (1+temp1).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2);
%          y(2,:)       = (1+temp1).*cos(x(1,:)*pi/2).*sin(x(2,:)*pi/2);
%          y(3,:)       = (1+temp1).*sin(x(1,:)*pi/2);
        M=3;
        [row,col] = size(x);
        if row>col
            x = x';
        end
        Xm = x(M:end);
        gXm =  sum((Xm-0.5).^2);
        cofficentMatrix = tril(ones(M,M),-1);
        X =cos( pi/2*repmat(x(1:M),M,1));
        X = X.*cofficentMatrix;
        for i = 1:M-1
            X(i,i) = sin(pi/2*x(i));
        end
        X(M,M) = 1;
        X = X+(tril(ones(M,M),-1))';
        y = cumprod(X,2);
        y = y(:,end)*(1+gXm);
    end
end