function value=IGD_test(val_cp,val_p)
    [dim_f1,num_cp] = size(val_cp);
    [dim_f2,num_p]  = size(val_p);
    value=0;
    if dim_f1 == dim_f2
        for i=1:num_cp
            dis=sqrt(sum((repmat(val_cp(:,i),[1,num_p])-val_p).^2,1));
            val=min(dis);
            value=value+val;
        end
        value=value/num_cp;
    else
        disp('the dimension of two stes is not equal!');
    end
end
