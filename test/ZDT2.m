function p=ZDT2(p)
 p.name='ZDT2';
 p.od = 2;
 p.pd = 30;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %ZDT2  evaluation function.
    function y = evaluate(x)
       temp         = sum(x(2:p.pd,:));
       g            = 1+9*temp/(p.pd-1);
       y(1,:)       = x(1,:);
       y(2,:)       = g.*(1-(y(1,:)./g).^2);
       clear Y;
    end
end