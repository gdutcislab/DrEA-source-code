function SP=n_sort(val_f)            
    [val_dim,num_p]  = size(val_f);
    rank             = zeros(1,num_p);
    prank            = zeros(1,num_p);
    ploc             = 1:num_p;
    count            = zeros(1,num_p);
    for i=1:num_p
        SP{i}        = [];
    end
    
%% count the rank
    for k=1:num_p
        for j=k+1:num_p
                if any(0.1<abs(val_f(:,k)-val_f(:,j)))
                    de=-0.0001;
                else
                    de=0;
                end
                if all(de<val_f(:,k)-val_f(:,j))
                    rank(k)  = rank(k)+1;
                    count(j) = count(j)+1;
                    SP{j}    = [SP{j},k];
                elseif all(de<val_f(:,j)-val_f(:,k))
                    rank(j)  = rank(j)+1;
                    count(k) = count(k)+1;
                    SP{k}    = [SP{k},j];
                end
            end
    end
   pind       = ploc(rank==0);  
      SP=val_f(:,pind);
  
end
