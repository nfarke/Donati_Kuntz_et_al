function [p,par]=Sample(ParSize,gg)

%parameter and variable names
p={'kcat1','kcat2','kcat3','kcat4','Km2','Kmorn','Kmcbp','Km4','Kmu','Kmu2','alpha1','alpha2','mumax','mumax2','orn','cbp','arg','utp','e2'};

%parameters average flux
if gg == 1 %called by create_SS_solutions_par -> for figure4b
    par_lb=[1, 1, 1, 1, 1, 32*3.3, 0.36*3.3, 0.028*3.3, 7E-5, 0.1*3.3, 95.8, 46.7, 0.01, 0.01, 1, 1, 1, 1, 1]';
    par_ub=[1, 1, 1, 1, 4, 32/3.3, 0.36/3.3, 0.028/3.3, 7E-5, 0.1/3.3, 95.8, 46.7, 0.01, 0.01, 1, 1, 1, 1, 1]';
end

par=zeros(length(par_lb),ParSize);
%random sampling
for i=1:ParSize
    par(:,i)= 10.^((log10(par_lb)-log10(par_ub)).*rand(length(par_lb), 1)+log10(par_ub));
end
end