function p=MOP7(p)
 p.name='MOP7';
 p.od = 3;
 p.pd = 10;
 range = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %LDTLZ2  evaluation function.
    function y = evaluate(x)
      [dim, num]   = size(x);
      temp         = x(3:dim,:) - x(1,:).*x(2,:);
      Y            = zeros(dim,num);
      Y(3:dim,:)   = -0.9*temp.^2+abs(temp).^0.6;
      g            = 2*sum(Y).*sin(pi*x(1,:));
      y(1,:)       = (1+g).*cos(x(1,:)*pi/2).*cos(x(2,:)*pi/2);
      y(2,:)       = (1+g).*cos(x(1,:)*pi/2).*sin(x(2,:)*pi/2);
      y(3,:)       = (1+g).*sin(x(1,:)*pi/2);
      clear Y;
    end
end