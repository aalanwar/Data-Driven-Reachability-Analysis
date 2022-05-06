function [gamma,L]= compLipConst(fun,U,R0,steps,initpoints,dim_x)
%compute Lipschitz constant

totalsamples = initpoints*steps;
%input random sample points
for i=1:totalsamples
    u(:,i) = randPointExtreme(U);
end

%get state trajectories
index=1;
for j=1:dim_x:initpoints*dim_x
    x_free(j:j+dim_x-1,1) = randPoint(R0);
    for i=1:steps
        x_free(j:j+dim_x-1,i+1) = fun(x_free(j:j+dim_x-1,i),u(:,index));
        index=index+1;
    end
end


%combine trajectories
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1        
        x_free_vec_1(:,index_1) = x_free(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_free_vec_0(:,index_0) = x_free(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end
% x_free_vec_1= normalize(x_free_vec_1);
% x_free_vec_0 = normalize(x_free_vec_0);
% u = normalize(u);
normType=2;
for idx=1:dim_x
    L(idx) = [0];
    gamma(idx) = 0;
    for i=1:totalsamples
        z1= [x_free_vec_0(idx,i);u(:,i)];
        f1= x_free_vec_1(idx,i);
        for j=1:totalsamples
            z2= [x_free_vec_0(idx,j);u(:,j)];
            f2= x_free_vec_1(idx,j);
            if z1~=z2
                newnorm = norm(f1-f2,normType)./norm(z1-z2,normType);
                newgamma = norm(z1-z2,normType);
                if newnorm > L(idx)
                    L(idx) = newnorm;
                end
                if newgamma > gamma(idx)
                    gamma(idx) = newgamma;
                end
            end
        end
    end

end
end

% L = 0;
% gamma = 0;
% for i=1:totalsamples
%     z1= [x_free_vec_0(:,i);u(:,i)];
%     f1= x_free_vec_1(:,i);
%     for j=1:totalsamples
%         z2= [x_free_vec_0(:,j);u(:,j)];    
%         f2= x_free_vec_1(:,j);
%         if z1~=z2
%             newnorm = norm(f1-f2)./norm(z1-z2);
%             newgamma = norm(z1-z2);
%             if newnorm(1) > L(1)
%                 L(1) = newnorm(1);
%             end
% %             if newnorm(2) > L(2)
% %                 L(2) = newnorm(2);
% %             end            
%             if newgamma > gamma
%                 gamma(1) = newgamma(1);
%             end
% %             if newgamma(2) > gamma
% %                 gamma(2) = newgamma(2);
% %            end            
%         end
%     end
% end
% eps(2)= L .* gamma;
% end
