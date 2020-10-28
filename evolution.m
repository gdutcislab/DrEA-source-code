function [params,mop,pop,state]=evolution(params,mop,pop,state)

    % 1)crossover&mutation
    params.delta          = (1- state.reg_obj/params.iteration)^0.7;
    [newpop,state]        = CMOp(params,mop,pop,state);    
    % 2)selection
    [pop,params,state]    = extractPop(pop,newpop,params,state);   
end
%% update population
function [pop,params,state] = extractPop(pop,newpop,params,state)  
          allind              = [pop.inter,newpop.inter];
          val                 = [allind.objective];            
          F                   = non_domination(val,size(val,1));
          f=[];
          i=1;
          while size(f,2) < params.popsize
              f = [f,F(i).f];
              i = i+1;
          end
          params.idealpoint   = min([params.idealpoint,val(:,[F(i).f])],[],2);
          params              = updatenadir(val(:,[F(1).f]),params);
          allind              = allind(f);
          val                 = val(:,f);
          [dim,num_nod]       = size(val);          
          noval               = (val-repmat(params.idealpoint,[1,num_nod]))./(repmat(params.nadirpoint-params.idealpoint,[1,num_nod]));
          [team,allind]       = group(params, pop, noval, allind);       
          num_p               = zeros(1,params.num_class);         
          for i=1:params.num_class
                 num_p(i)        = length(team{i}); 
                 pop(i).num_ind  = 0;
%               if num_p(i)<params.popsize/params.num_class
%                  pop(i).num_ind  = num_p(i);
%                  num_p(i)        = 0;
%               else
%                  pop(i).num_ind  = floor(params.popsize/params.num_class);
%                  num_p(i)        = num_p(i)-pop(i).num_ind; 
%               end            
          end
          while sum([pop.num_ind])<params.popsize              
               [res,index] = sort(num_p);                      
                for i=1:params.num_class
                   if num_p(index(i))>0
                       pop(index(i)).num_ind=pop(index(i)).num_ind+1;                      
                       num_p(index(i))=num_p(index(i))-1;
                       if sum([pop.num_ind])==params.popsize 
                           break;
                       end
                   end
                end           
          end    
          for i=1:params.num_class
              num_p            = length(team{i});
              if num_p==pop(i).num_ind   
                  pop(i).inter                = allind([team{i}]);  
              else
                  selind                      = allind(team{i});
                  pop(i).inter                = selection(params,selind,pop(i),state);
              end
          end
end
%%
function params=updatenadir(valf,params)
          dim = size(valf,1);                
          W   = diag(ones(1,dim))+1e-5-diag(10e-5*ones(1,dim));
%% FindExtremePoints         
          extrempoint         = params.extrempoint;
          if any((params.nadirpoint-params.idealpoint)==0)
              disp('error');
          end
          z                   = ones(dim,dim);           
          for j=1:dim             
              z(:,j) = argASF([valf,extrempoint(:,j)],W(:,j),params.idealpoint,params.nadirpoint);               
          end 
          params.extrempoint = z;         
          z1=z-repmat(params.idealpoint,[1,dim]);          
         if rank(z1)==dim           
            b=ones(dim,1);
            a=z1'\b;
            intercepts=1./a;           
             if ~any(intercepts<0)&&~any(~isfinite(intercepts))&&~any(isnan(intercepts))
                 params.nadirpoint = intercepts+params.idealpoint;
             else
                 params.nadirpoint = max(valf,[],2); 
             end
         else
             params.nadirpoint = max(valf,[],2); 
         end
 
end
%%
function y = argASF(val_f,weight,idealpoint,nadirpoint)
         num_nod      = size(val_f,2);        
         noval        = abs(val_f-repmat(idealpoint,[1,num_nod]));
         max_evo      = max(noval./(repmat(weight,[1,num_nod])));
         [cfval,cfloc]= min(max_evo);
         y            = val_f(:,cfloc);
end
%% generate next offsprings
function [newpop,state] = CMOp(params,mop,pop,state)
         PM       = [pop.inter];
         newpop   = pop;
         x_min    = mop.domain(:,1);
         x_max    = mop.domain(:,2);
         for i=1:params.num_class
             newinter  = pop(i).inter;
             for j=1:pop(i).num_ind
                 parent1=pop(i).inter(j);
                 if rand>params.selectPro
                     parent2  = PM(floor(rand*params.popsize)+1);
                 else
                     loc                 = floor(rand*pop(i).num_ind)+1;
                     parent2             = pop(i).inter(loc);                   
                 end
                 while isequal(parent1,parent2)
                     parent2 = PM(floor(rand*params.popsize)+1);
                 end
                 newinter(j).parameter   = SBXCM(parent1.parameter,parent2.parameter,x_min,x_max,1/mop.pd);
             end
             [newinter,state]            = evaluate(newinter,mop,state);
             newpop(i).inter             = newinter;
         end
end
