function p=DTLZ3(p)
 p.name='DTLZ3';
 p.od = 10;
 p.pd = 30;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %DTLZ1  evaluation function.
    function y = evaluate(x)
%          temp         = sum((x(3:p.pd,:)-0.5).^2 - cos(1*pi*(x(3:p.pd,:)-0.5)));
%          temp1        = 1*(p.pd-2)+1*temp;
%          y(1,:)       = (1+temp1).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2).*cos(x(3,:)*pi/2).*cos(x(4,:)*pi/2);
%          y(2,:)       = (1+temp1).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2).*cos(x(3,:)*pi/2).*sin(x(4,:)*pi/2);
%          y(3,:)       = (1+temp1).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2).*sin(x(3,:)*pi/2);
%          y(4,:)       = (1+temp1).*cos(x(1,:)*pi/2).*sin(x(2,:)*pi/2);
%          y(5,:)       = (1+temp1).*sin(x(1,:)*pi/2);
        M=10;
        [row,col] = size(x);
        if row>col
            x = x';
        end
        Xm = x(M:end);
%         gXm =  sum((Xm-0.5).^2);
        gXm = 100*(length(Xm)+ sum((Xm-0.5).^2-cos(20*pi*(Xm-0.5))));
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