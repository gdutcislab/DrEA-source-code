function p=MOP6(p)
 p.name='MOP6';
 p.od = 3;
 p.pd = 10;
 range = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %LDTLZ1  evaluation function.
   
    function y = evaluate(x)
        [dim, num]   = size(x);
        temp         = x(3:dim,:) - x(1,:).*x(2,:);
        Y            = zeros(dim,num); 
        Y(3:dim,:)   = -0.9*temp.^2+abs(temp).^0.6;
        g            = 2*sum(Y).*sin(pi*x(1,:));
        y(1,:)       = (1+g).*x(1,:).*x(2,:);
        y(2,:)       = (1+g).*(x(1,:).*(1-x(2,:)));
        y(3,:)       = (1+g).*(1-x(1,:));
        clear Y;
    end
end