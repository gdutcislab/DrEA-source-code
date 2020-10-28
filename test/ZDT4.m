function p=ZDT4(p)
 p.name='ZDT4';
 p.od = 2;
 p.pd = 30;
 range           = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %ZDT4  evaluation function.
    function y = evaluate(x)
           temp         = 2*x(2:p.pd,:)-1;      
           g            = 1+10*(p.pd-1)+sum(temp.^2-10*cos(pi*temp));
           y(1,:)       = x(1,:);
           y(2,:)       = g.*(1-sqrt(y(1,:)./g));
    end
end