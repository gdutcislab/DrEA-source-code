function p=ZDT1(p)
 p.name='ZDT1';
 p.od = 2;
 p.pd = 30;
 range = ones(p.pd,2); 
 range(:,1)      =  0; 
 p.domain=  range;
 p.func = @evaluate;
    %ZDT1  evaluation function.
    function y = evaluate(x)
      temp         = sum(x(2:p.pd,:));
      g            = 1+9*temp./(p.pd-1);
      y(1,:)       = x(1,:);
      y(2,:)       = g.*(1-sqrt(y(1,:)./g));
      clear Y; 
    end
end