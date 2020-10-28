function p=DTLZ5(p)
 p.name = 'DTLZ5';
 p.od      = 10;
 p.pd      = 30;
 range     = ones(p.pd,2); 
 range(:,1) =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %DTLZ2  evaluation function.
    function y = evaluate(x)
%          temp         = (2*x(p.od:p.pd,:)-1).^2;
%          temp1        = sum(temp);
%          alpha        = 100;
%          y(1,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*cos(x(8,:).^alpha*pi/2).*cos(x(9,:).^alpha*pi/2)*cos(x(10,:).^alpha*pi/2).*cos(x(11,:).^alpha*pi/2).*cos(x(12,:).^alpha*pi/2).*cos(x(13,:).^alpha*pi/2).*cos(x(14,:).^alpha*pi/2);
%          y(2,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*cos(x(8,:).^alpha*pi/2).*cos(x(9,:).^alpha*pi/2)*cos(x(10,:).^alpha*pi/2).*cos(x(11,:).^alpha*pi/2).*cos(x(12,:).^alpha*pi/2).*cos(x(13,:).^alpha*pi/2).*sin(x(14,:).^alpha*pi/2);
%          y(3,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*cos(x(8,:).^alpha*pi/2).*cos(x(9,:).^alpha*pi/2)*cos(x(10,:).^alpha*pi/2).*cos(x(11,:).^alpha*pi/2).*cos(x(12,:).^alpha*pi/2).*sin(x(13,:).^alpha*pi/2);
%          y(4,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*cos(x(8,:).^alpha*pi/2).*cos(x(9,:).^alpha*pi/2)*cos(x(10,:).^alpha*pi/2).*cos(x(11,:).^alpha*pi/2).*sin(x(12,:).^alpha*pi/2);
%          y(5,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*cos(x(8,:).^alpha*pi/2).*cos(x(9,:).^alpha*pi/2)*cos(x(10,:).^alpha*pi/2).*sin(x(11,:).^alpha*pi/2);
%          y(6,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*cos(x(8,:).^alpha*pi/2).*cos(x(9,:).^alpha*pi/2)*sin(x(10,:).^alpha*pi/2);
%          y(7,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*cos(x(8,:).^alpha*pi/2).*sin(x(9,:).^alpha*pi/2);
%          y(8,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*cos(x(7,:).^alpha*pi/2).*sin(x(8,:).^alpha*pi/2);
%          y(9,:)       = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*cos(x(6,:).^alpha*pi/2).*sin(x(7,:).^alpha*pi/2);
%          y(10,:)      = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*cos(x(5,:).^alpha*pi/2).*sin(x(6,:).^alpha*pi/2);
%          y(11,:)      = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*cos(x(4,:).^alpha*pi/2).*sin(x(5,:).^alpha*pi/2);
%          y(12,:)      = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*cos(x(3,:).^alpha*pi/2).*sin(x(4,:).^alpha*pi/2);
%          y(13,:)      = (1+temp1).*cos(x(1,:).^alpha*pi/2).*cos(x(2,:).^alpha*pi/2).*sin(x(3,:).^alpha*pi/2);
%          y(14,:)      = (1+temp1).*cos(x(1,:).^alpha*pi/2).*sin(x(2,:).^alpha*pi/2);
%          y(15,:)      = (1+temp1).*sin(x(1,:).^alpha*pi/2);
        M=10;
        [row,col] = size(x);
        if row>col
            x = x';
        end
        Xm = x(M:end);
        gXm =  sum((Xm-0.5).^2);
        cofficentMatrix = tril(ones(M,M),-1);
        X =cos( pi/2*repmat((1+2*gXm*x(1:M))./(4*(1+gXm)),M,1).^alphalParameter);
        X = X.*cofficentMatrix;
        for i = 1:M-1
            X(i,i) = sin(pi/2*x(i).^alphalParameter);
        end
        X(M,M) = 1;
        X = X+(tril(ones(M,M),-1))';
        y = cumprod(X,2);
        y = y(:,end)*(1+gXm);
        
    end    
end