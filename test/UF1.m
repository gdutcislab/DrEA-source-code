function p=UF1(p)
 p.name='UF1';
 p.od = 2;
 p.pd = 30;
 range = ones(p.pd,2); 
 range(:,1)     =  -1; 
 range(1,1)     =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %UF1  evaluation function.
   function y = evaluate(x)
       [dim, num]  = size(x);
       tmp         = zeros(dim,num);
       tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
       tmp1        = sum(tmp(3:2:dim,:));  % odd index
       tmp2        = sum(tmp(2:2:dim,:));  % even index
       y(1,:)      = x(1,:)             + 2.0*tmp1/size(3:2:dim,2);
       y(2,:)      = 1.0 - sqrt(x(1,:)) + 2.0*tmp2/size(2:2:dim,2);
       clear tmp;
   end
end