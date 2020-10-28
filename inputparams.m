function [params,mop,pop]=inputparams(name_f,params)   
    
    % 1)初始化种群规模，最大进化代数，子种群的数目
    mop           = testmop(name_f);
    [params,mop]  = inputparameter(params,mop);
    
    % 2）初始化真是的有效界面
%     if( strcmpi(params.isCauIGD, 'yes'))
%         PFStar   = load(strcat('PFStar/pf_',name_f,'_',num2str(mop.od),'.dat'));
%         params.PFStar   =  PFStar';
%     end
    
    % 4)初始化权重
    if( strcmpi(params.useWeight, 'yes'))
        val_w                  = reference_point(mop);
        val_w                  = val_w./repmat(sqrt(sum(val_w.^2)),[mop.od,1]);
        params.weight          = val_w;
        params.popsize         = size(val_w,2);
    end
    
    % 5)初始化中心点和子种群的中心规模
    val_cp                 = reference_point(mop);
%     val_cp                 = ones(mop.od,1);
    center                 = val_cp./repmat(sqrt(sum(val_cp.^2)),[mop.od,1]);
    params.num_class       = size(center,2);
    sub                    = get_structure('subclass');
    pop                    = repmat(sub, [1,params.num_class]);
    for i=1:params.num_class
        pop(i).center         = center(:,i);
    end

    if( strcmpi(params.useWeight, 'yes'))
        team                      = group(params,pop,val_w);
        for i=1:params.num_class
            pop(i).weight         = 1./val_w(:,team{i});
            pop(i).num_ind        = length(team{i});
        end
    else
        temv1     = mod(params.popsize,params.num_class);
        temv2     = floor(params.popsize/params.num_class);
        for i=1:temv1
            pop(i).num_ind=temv2+1;           
        end
        for i=temv1+1:params.num_class
            pop(i).num_ind=temv2;          
        end
    end   
    params.pmuta           = 1/mop.pd;
    params.idealpoint      = Inf*ones(mop.od,1);     
end

function [params,mop]=inputparameter(params,mop)
%     switch upper(mop.name)
   switch upper(mop.od)
        case {3}
            params.iteration   = 400*92;
            params.popsize     = 92;
        case {5}
            params.iteration   = 750*210;
            params.popsize     = 210;
        case {8}
            params.iteration   = 1500*156;
            params.popsize     = 156;
        case {10}
            params.iteration   = 2000*275;
            params.popsize     = 156;
        case {'MOP1','MOP2','MOP3','MOP4','MOP5'}
            params.iteration   = 300000;
            params.popsize     = 100;
        case {'MOP6','MOP7'}
            params.iteration   = 900000;
            params.popsize     = 300;
        case {'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5'}         
            params.iteration   = 300*1000;
            params.popsize     = 300;
        case {'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'}           
            params.iteration   = 750*210;
            params.popsize     = 210;
    end
end

function val_w = Weight(popsize,objDim)
    if objDim==2
        start      = 1/(popsize*100000);
        val_w(1,:) = linspace(start,1-start,popsize);
        val_w(2,:) = ones(1,popsize)-val_w(1,:);
    elseif objDim==3
        val_w      = assign3(popsize,1);
    else        
        val_w = lhsdesign(popsize-objDim, objDim, 'criterion','maximin', 'iterations', 1000)';
        val_w = [val_w,eye(objDim)];
    end
end
function val_w = assign3(popsize,sideL) 
     val_w          = zeros(3,popsize);   
     [numline,numpoint]=comline2(popsize);
     start      = sideL/(numline*10);
     val_one    = linspace(start,sideL-start,numline);
     reg        = 1;
     for i=1:numline
         nump   = numpoint(numline-i+1);
         if nump>1
             val    = zeros(3,nump);
             Lline  = sideL-val_one(i);
             lstart = Lline/(10*nump);
             val(1,:)=linspace(lstart,Lline-lstart,nump);
             val(2,:)=Lline-val(1,:);
             val(3,:)=val_one(i);
         else
             val    = zeros(3,nump);
             val(1)=(sideL-val_one(i))/2;
             val(2)=(sideL-val_one(i))/2;
             val(3)=val_one(i);
         end
         val_w(:,reg:reg+nump-1)=val;
         reg=reg+nump;
     end
end

function [numline,numpoint]=comline2(popsize)
    numline    = round((sqrt(popsize*8+1)-1)/2);
    numpoint   = 1:numline;
    vp         = (numline+1)*numline/2;
    disnum     = popsize-vp;
    numpoint(numline-abs(disnum)+1:numline)=numpoint(numline-abs(disnum)+1:numline)+sign(disnum)*ones(1,abs(disnum));  
end

