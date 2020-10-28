 clear; close; clc;                           
name_func  = {'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9',...
               'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7'}; %  1-7 
    digits(20)  
HV = zeros(9,10);
    %% 算法初始化
    for seq =2:9
        for nrun=1:15
            name_f   = char(name_func(seq));
            seed     = 60+nrun;randn('state',seed);rand('state',seed);
            state    = get_structure('state');
            params   = get_structure('parameter');
            [params,mop,pop]=inputparams(name_f,params);
            [params,mop,pop,state]=initialize(params,mop,pop,state);        
            while state.stopCriterion 
                [params,mop,pop,state]=evolution(params,mop,pop,state);
                % 当前状态输出
%                 state=stateOutput(state,params,pop,mop,nrun);
%                 t1=[pop.inter];
%                 t2=[t1.objective];
%                 figure(1); 
%                 drawnow 
%                 plot(t2);               
%                 f =  Hypervolume( (t2'./repmat(2*(1:mop.od),[size(t2,2),1]))', ones(mop.od,1)+0.1)
%                 fprintf('gen =   %d\n',state.currentGen);           
                state.currentGen=state.currentGen+1;
                % 检测是否终止
                state  = checkstop(params,state);
            end
            t1   = [pop.inter];
            valf = [t1.objective];
%             figure(1);
%             plot(t2)
%             drawnow
%              plot(t2(1,:),t2(2,:),'bo');
%             plot3(t2(1,:),t2(2,:),t2(3,:),'bo');
%               HV(seq,nrun) =  Hypervolume_MEX( valf'./repmat(2*(1:mop.od),[size(valf,2),1]), ones(1,mop.od)+0.1);
%               HV(seq,nrun) = TestHyperVolume5(valf', 2*(1:mop.od)+0.1)/prod(2*(1:mop.od)+0.1); 
            state=stateOutput(state,params,pop,mop,nrun);
        end
    end
  