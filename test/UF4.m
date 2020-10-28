function p=UF4(p)
 p.name='UF4';
 p.od = 2;
 p.pd = 10;
 range = 2*ones(p.pd,2); 
 range(:,1)     =  -2; 
 range(1,1)     =  0; 
 range(1,2)     =  1; 
 p.domain=  range;
 p.func = @evaluate;
    %UF1  evaluation function.
   function y = evaluate(x)      
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    H           = zeros(dim,num);
    H(2:dim,:)  = abs(Y(2:dim,:))./(1.0+exp(2.0*abs(Y(2:dim,:))));
    tmp1        = sum(H(3:2:dim,:));  % odd index
    tmp2        = sum(H(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)          + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - x(1,:).^2 + 2.0*tmp2/size(2:2:dim,2);
    clear Y H;
   end
end