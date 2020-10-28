function p=MOP5(p)
 p.name='MOP5';
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
      Y(2:dim,:)   = -0.9*temp.^2+abs(temp).^0.6;
      g            = 1+2*abs(cos(pi*x(1,:))).*sum(Y);
      y(1,:)       = g*x(1,:);
      y(2,:)       = g.*(1-sqrt(x(1,:)));
      clear Y;
    end
end