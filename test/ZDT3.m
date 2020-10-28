function p=ZDT3(p)
 p.name='ZDT3';
 p.od = 2;
 p.pd = 30;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %ZDT3  evaluation function.
    function y = evaluate(x)
        temp         = sum(x(2:p.od,:));
        g            = 1+9*temp/(p.od-1);
        y(1,:)       = x(1,:);
        y(2,:)       = g.*(1-sqrt(y(1,:)./g)-y(1,:)./g.*sin(10*pi*x(1,:)));
        clear Y;
    end
end