function state=stateOutput(state,params,pop,mop,nrun)
    if params.isDebug==1||~state.stopCriterion
        gen  = state.currentGen;
        if (gen>0&&gen<11)||(mod(gen,10)==0&&gen<100)||(mod(gen,100)==0)||~state.stopCriterion
        % 当前种群的自变量与目标函数值
            if isempty(state.archive)
                individual    = [pop.inter];
                valf          = [individual.objective];
                valx          = [individual.parameter];
            else
                individual    = [state.archive];
                valf          = [individual.objective];
                valx          = [individual.parameter];
            end
            %计算当前种群的IGD值和H值
%             if strcmpi(computer('arch'),'win64')
%                 state.val_H=          Hypervolume( (valf'./repmat(2*(1:mop.od),[size(valf,2),1]))', ones(mop.od,1)+0.1);
%             else
%                 state.val_H=          Hypervolume_MEX( valf'./repmat(2*(1:mop.od),[size(valf,2),1]), ones(1,mop.od)+0.1);
%             end
%               f =  Hypervolume( (valf'./repmat(2*(1:mop.od),[size(valf,2),1]))', ones(mop.od,1)+0.1)
%               f =  Hypervolume_MEX( valf'./repmat(2*(1:mop.od),[size(valf,2),1]), ones(1,mop.od)+0.1)
%               state.val_IGD=Hypervolume_MEX(valf',2*(1:mop.od)+1)./prod(2*(1:mop.od)+1);
%             f =  hypeSampling( valf'./repmat(2*(1:mop.od),[size(valf,2),1]), ones(1,mop.od)+0.1)
%             IGD_test(params.PFStar,valf)
%             state.val_IGD       = [state.val_IGD;state.currentGen, IGD_test(params.PFStar,valf)];
%             name      = mop.name;
%             file      = state.val_IGD';
%             filename  = strcat('IGD_',name,'_R',num2str(nrun)); 
%             filetype  = 'dat';
%             folder    = 'IGD_data';
%             mySave(file,filename,filetype,folder)
%             state.val_H         = [state.val_H;state.currentGen,Hypervolume_MEX(valf',max(params.PFStar,[],2))];
            %输出当前状态
            output(valf,valx,nrun,params,mop,state);
        end 
    end
end




function output(val_f,val_x,nrun,params,mop,state)
    name         = mop.name;

%% 在算法结束时保存IGD
%     if ~state.stopCriterion
%         file      = state.val_IGD';
%         filename  = strcat('IGD_',name,'_R',num2str(nrun)); 
%         filetype  = 'dat';
%         folder    = 'IGD_data';
%         mySave(file,filename,filetype,folder)
%     end
    
%% 在算法结束时保存H值
%     if ~state.stopCriterion
%         file      = state.val_H';
%         filename  = strcat('H_',name,'_R',num2str(nrun)); 
%         filetype  = 'dat';
%         folder    = 'H_data';
%         mySave(file,filename,filetype,folder)
%     end
    
%% 保存当前的和最后的种群自变量的值
    if ~state.stopCriterion ||  ( strcmpi(params.resultOut, 'save'))
        file      = val_x';
        filename  = strcat('PS_',name,'_R',num2str(nrun),'gen',num2str(state.currentGen)); 
        filetype  = 'dat';
        folder    = 'PS_data';
        mySave(file,filename,filetype,folder)
    end
%% 保存当前的和最后的种群自变量的值
    if ~state.stopCriterion ||  ( strcmpi(params.resultOut, 'save'))
        file      = val_f';
        filename  = strcat('PF_',name,'_R',num2str(nrun),'gen',num2str(state.currentGen)); 
        filetype  = 'dat';
        folder    = 'PF_data';
        mySave(file,filename,filetype,folder)
    end

    
%% 最后作图时删除劣解
    if ~state.stopCriterion
        [val_dim,num_p]  = size(val_f);
        rank             = zeros(1,num_p);
        for k=1:num_p
            for j=k+1:num_p
                    if all(val_f(:,j)<=val_f(:,k))
                        rank(k)  = rank(k)+1;
                    elseif all(val_f(:,k)<=val_f(:,j))
                        rank(j)  = rank(j)+1;
                    end
            end
        end
        val_f(:,rank~=0)=[];
    end

%% 做有效界面图像
%     close;
%     H=figure(1);
%     switch name
%         case {'UF4'}
%             pf          = zeros(2,100);
%             t           = linspace(0,1,100);
%             pf(1,:)     = t;
%             pf(2,:)     = 1-t.^2;
%             h1=plot(pf(1,:),pf(2,:),'-','Color',[0.3137 0.3137 0.3137]);
%             hold on; 
%             clear pf;
%        case {'MOP3'}
%             pf          = zeros(2,100);
%             t           = linspace(0,1,100);
%             pf(1,:)     = cos(0.5*pi*t);
%             pf(2,:)     = sin(0.5*pi*t);
%             h1=plot(pf(1,:),pf(2,:),'-','Color',[0.3137 0.3137 0.3137]);
%             hold on; 
%             clear pf;
%         case {'MOP1','MOP5'}
%             pf          = zeros(2,100);
%             t           = linspace(0,1,100);
%             pf(1,:)     = t;
%             pf(2,:)     = 1-sqrt(t);
%             h1=plot(pf(1,:),pf(2,:),'-','Color',[0.3137 0.3137 0.3137]);
%             hold on; 
%             clear pf;
%         case {'MOP7'}
%             [X,Y]=meshgrid(0:0.1:1);
%             A1=cos(X*pi/2).*cos(Y*pi/2);
%             A2=cos(X*pi/2).*sin(Y*pi/2);
%             A3=sin(X*pi/2);
%             h1=mesh(A1,A2,A3);
%             colormap([0.3137 0.3137 0.3137]);
%             hold on
%        case {'MOP4'}
%             h1=plot(params.PFStar(1,1:93),params.PFStar(2,1:93),'-','Color',[0.3137 0.3137 0.3137]);
%             hold on
%             plot(params.PFStar(1,94:344),params.PFStar(2,94:344),'-','Color',[0.3137 0.3137 0.3137]);
%             plot(params.PFStar(1,345:500),params.PFStar(2,345:500),'-','Color',[0.3137 0.3137 0.3137]);
%             hold on
%         case {'MOP6'}
%             h1=ezmesh('x1-x1*x2','x1*x2','1-x1',[0,1],[0,1],10);
%             colormap([0.3137 0.3137 0.3137]);
%             hold on
%         otherwise
%             close;
%     end
%     
% %% 做有效点的图像  
%     
%     dim_f     = size(val_f,1);
%     
%     if dim_f == 2 
%         plot(val_f(1,:),val_f(2,:),'MarkerFaceColor',[0.502 0.502 0.502],...
%                'MarkerEdgeColor',[0.3137 0.3137 0.3137],...
%                'MarkerSize',10,...
%                'Marker','o','LineStyle','none'); 
%         hold off;
%         xlabel('f1','FontWeight','bold','FontSize',20,'FontName','High Tower Text','Rotation',0);
%         ylabel('f2','FontWeight','bold','FontSize',20,'FontName','High Tower Text');
%         title(name,'FontWeight','bold','FontSize',20);
%         axis([0,1.01,0,1.01]);
%         box off;
%     else
%         plot3(val_f(1,:),val_f(2,:),val_f(3,:),'MarkerFaceColor',[0.502 0.502 0.502],...
%                'MarkerEdgeColor',[0.3137 0.3137 0.3137],...
%                'MarkerSize',10,...
%                'Marker','o','LineStyle','none'); 
%         axis([0,1.01,0,1.01,0,1.01]);
%         grid off
%         hold off;
%         xlabel('f1','FontWeight','bold','FontSize',20,'FontName','High Tower Text');
%         ylabel('f2','FontWeight','bold','FontSize',20,'FontName','High Tower Text');
%         zlabel('f3','FontWeight','bold','FontSize',20,'FontName','High Tower Text');
%         title(name,'FontWeight','bold','FontSize',20);
%     end
%     drawnow
%     %% 保存图像
%     if ~state.stopCriterion ||  ( strcmpi(params.resultOut, 'save'))
%         file = H;
%         filename = strcat('FIG_',name,'_R',num2str(nrun),'gen',num2str(state.currentGen));
%         filetype = 'fig';
%         folder   = 'fig_data';
%         mySave(file,filename,filetype,folder)
%     end
end

function mySave(file,filename,filetype,folder)
    name      = strcat(filename,'.',filetype);
    if ( strcmpi(filetype, 'dat'))
        fid       = fopen(name,'w');
        [numline,numi]   = size(file);
        for i=1:numline
            for j=1:numi
                fprintf(fid,'%8.6f ',file(i,j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
    elseif  ( strcmpi(filetype, 'fig'))
        saveas(file,name);
    end
    judge     = exist(folder);
    if judge ~= 7
        system(['mkdir ', folder]);
    end
    file_path = strcat(cd,'\',folder);
    movefile(name,file_path); 
end