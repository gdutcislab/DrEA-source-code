function p=MOP4(p)
 p.name='MOP4';
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
      g            = 1+10*sin(pi*x(1,:)).*sum(Y);
      y(1,:)       = g.*x(1,:);
      y(2,:)       = g.*(1-(x(1,:).^0.5).*cos(2*pi*x(1,:)).^2);
      clear Y;
    end
end
