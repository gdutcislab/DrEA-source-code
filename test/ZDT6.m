function p=ZDT6(p)
 p.name='ZDT6';
 p.od = 2;
 p.pd = 30;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %ZDT4  evaluation function.
    function y = evaluate(x)
           temp         = sum(x(2:p.pd,:));
           g            = 1+9*(temp/(p.pd-1)).^0.25;
           y(1,:)       = 1-exp(-4*x(1,:)).*sin(6*pi*x(1,:)).^6;
           y(2,:)       = g.*(1-(y(1,:)./g).^2);
    end
end