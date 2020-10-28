function p=wfg4(p)
p.name='WFG4';
p.od     = 10;
p.k      = (p.od-1);
p.l      = 20;
p.pd     = p.k+p.l;
range           = ones(p.pd,2);
range(:,2)      = 2*(1:p.pd);
range(:,1)      =  0;
p.domain=  range;
p.func=@evaluate;
function y = evaluate(x)
seq=4;
y = wfg(x',p.od,p.k,p.l,seq);
y = y';
end
end