function p=wfg1(p)
p.name='WFG1';
p.od     = 5;
p.k      = (p.od-1);
p.l      = 20;
p.pd     = p.k+p.l;
range           = ones(p.pd,2);
range(:,2)      = 2*(1:p.pd);
range(:,1)      =  0;
p.domain=  range;
p.func=@evaluate;
function y = evaluate(x)
seq=1;
y = wfg_c(x',p.od,p.k,p.l,seq);
% y1 = wfg(x',p.od,p.k,p.l,seq);
y = y';
end
end