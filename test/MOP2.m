function p=MOP2(p)
 p.name='MOP2';
 p.od = 2;
 p.pd = 10;
 range = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %LZDT2  evaluation function.
    function y = evaluate(x)
      [dim, num]   = size(x);
      temp         = x(2:dim,:) - sin(pi*repmat(x(1,:),[dim-1,1])/2);
      Y            = zeros(dim,num);
      Y(2:dim,:)   = abs(temp)./(1+exp(5*abs(temp)));
      g            = sum(Y);
      g            = 1+10*g*abs(sin(x(1,:)*pi));
      y(1,:)       = g.*x(1,:);
      y(2,:)       = g.*(1-x(1,:).^2);
      clear Y;
    end
end