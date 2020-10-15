function create_SS_solutions_par
ParSize = 1000; %adjust number of Parameter sets
tspan = 0:1:300;
par_end = zeros(19,ParSize);
sol_wt = zeros(5,ParSize);
sol_af = zeros(5,ParSize);
opts = odeset('RelTol',1e-07,'AbsTol',1E-7);
load flux_ss
%syms kcat1 kcat2 kcat3 kcat4 Km2 Kmorn Kmcbp Km4 Kmu Kmu2 alpha1 alpha2 mumax mumax2 orn cbp arg utp e2

for i1 = 1:ParSize %parfor is possible
    disp(i1)
    value = 0; %repeat when stability criteria are not met
    while value == 0
        [p,par] = Sample(1,1); %Get random parameter sets
        y0 = par(end-4:end); %initial conditions
        for af = 0:1  %0 = dysregulated, 1 = regulated
            [par]=calc_assign_par(p,par,flux_ss,af);
            [~,y]  =  ode23s(@(t,c) odemodel(t,c,p,par,af),tspan,y0,opts);
            if af == 0
               orn_noreg = y(:,1);
               noreg(:,i1) = orn_noreg;
            else
              orn_reg = y(:,1);
              reg(:,i1) = orn_reg;
            end
            %Parameters;
            %Variables at t(end)
            orn     =       y(end,1);
            cbp     =       y(end,2);
            arg     =       y(end,3);
            utp     =       y(end,4);
            e2      =       y(end,5);

            Km2     =       par(find(strcmp(p,'Km2')),1);
            Kmorn   =       par(find(strcmp(p,'Kmorn')),1);
            Kmcbp   =       par(find(strcmp(p,'Kmcbp')),1);
            Km4     =       par(find(strcmp(p,'Km4')),1);
            Kmu     =       par(find(strcmp(p,'Kmu')),1);
            Kmu2    =       par(find(strcmp(p,'Kmu2')),1);
            alpha1  =       par(find(strcmp(p,'alpha1')),1);
            alpha2  =       par(find(strcmp(p,'alpha2')),1);
            kcat1   =       par(find(strcmp(p,'kcat1')),1);             
            kcat2   =       par(find(strcmp(p,'kcat2')),1);             
            kcat3   =       par(find(strcmp(p,'kcat3')),1);             
            kcat4   =       par(find(strcmp(p,'kcat4')),1);             
            mumax   =       par(find(strcmp(p,'mumax')),1);             
            mumax2  =       par(find(strcmp(p,'mumax2')),1);             

            
            %Definition of the growth rate
            mue1    =       mumax * arg/(arg + Kmu);
            mue2    =       mumax2 * utp/(utp + Kmu2);
            %rates
            r1      =       kcat1;
            if af == 0
                    r2      =       kcat2 * e2;
            elseif af == 1
                    r2      =       kcat2 * e2 * orn^(Km2);
            end
            r3      =       kcat3 * 1/(1 + Kmorn*Kmcbp/(orn*cbp) + Kmorn/orn + Kmcbp/cbp);
            r4      =       kcat4 * cbp/(cbp + Km4);
            r5      =       alpha1 * mue1;
            r6      =       alpha2 * mue2;
            
            % mass balance
            dorn_dt =       r1 - r3;
            dcbp_dt =       r2 - r3 - r4;
            darg_dt =       r3 - r5;
            dutp_dt =       r4 - r6;
            de2_dt  =       0;

            %check if solutions are in steady state
            F       =       [dorn_dt;dcbp_dt; darg_dt;dutp_dt; de2_dt];
            if af == 0
                SS_noreg    =       sum(abs(F));
            elseif af == 1
                SS_reg      =       sum(abs(F));
            end
            
            %calculation of the Jacobian Matrix
            if af == 0
            J = [ -(kcat3*(Kmorn/orn^2 + (Kmcbp*Kmorn)/(cbp*orn^2)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2,                                                -(kcat3*(Kmcbp/cbp^2 + (Kmcbp*Kmorn)/(cbp^2*orn)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2,                                                             0,                                                0,     0
                  -(kcat3*(Kmorn/orn^2 + (Kmcbp*Kmorn)/(cbp*orn^2)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2, (cbp*kcat4)/(Km4 + cbp)^2 - (kcat3*(Kmcbp/cbp^2 + (Kmcbp*Kmorn)/(cbp^2*orn)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2 - kcat4/(Km4 + cbp),                                                             0,                                                0, kcat2
                  (kcat3*(Kmorn/orn^2 + (Kmcbp*Kmorn)/(cbp*orn^2)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2,                                                 (kcat3*(Kmcbp/cbp^2 + (Kmcbp*Kmorn)/(cbp^2*orn)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2, (alpha1*arg*mumax)/(Kmu + arg)^2 - (alpha1*mumax)/(Kmu + arg),                                                 0,     0
                                                                                                          0,                                                                                                             kcat4/(Km4 + cbp) - (cbp*kcat4)/(Km4 + cbp)^2,                                                             0, (alpha2*mumax2*utp)/(Kmu2 + utp)^2 - (alpha2*mumax2)/(Kmu2 + utp),     0
                                                                                                          0,                                                                                                                                                         0,                                                             0,                                                                 0,     0];
            elseif af == 1
            J = [-(kcat3*(Kmorn/orn^2 + (Kmcbp*Kmorn)/(cbp*orn^2)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2,                                                -(kcat3*(Kmcbp/cbp^2 + (Kmcbp*Kmorn)/(cbp^2*orn)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2,                                                             0,                                                                 0,             0
                 Km2*e2*kcat2*orn^(Km2 - 1) - (kcat3*(Kmorn/orn^2 + (Kmcbp*Kmorn)/(cbp*orn^2)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2, (cbp*kcat4)/(Km4 + cbp)^2 - (kcat3*(Kmcbp/cbp^2 + (Kmcbp*Kmorn)/(cbp^2*orn)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2 - kcat4/(Km4 + cbp),                                                             0,                                     0, kcat2*orn^Km2
                 (kcat3*(Kmorn/orn^2 + (Kmcbp*Kmorn)/(cbp*orn^2)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2,                                                 (kcat3*(Kmcbp/cbp^2 + (Kmcbp*Kmorn)/(cbp^2*orn)))/(Kmcbp/cbp + Kmorn/orn + (Kmcbp*Kmorn)/(cbp*orn) + 1)^2, (alpha1*arg*mumax)/(Kmu + arg)^2 - (alpha1*mumax)/(Kmu + arg),                                                                 0,              0
                                                                                                                                      0,                                                                                                             kcat4/(Km4 + cbp) - (cbp*kcat4)/(Km4 + cbp)^2,                                                 0, (alpha2*mumax2*utp)/(Kmu2 + utp)^2 - (alpha2*mumax2)/(Kmu2 + utp),             0
                                                                                                                                      0,                                                                                                                                                         0,                                                 0,                                                                 0,             0];
           end
            
            %get eigenvalues and steady state solutions
            if af == 0
                val         =       real(eig(J));
                ew_noreg    =       max(val(1:4));
                sol_noreg   =       [orn;cbp;arg;utp;e2];
                parx_noreg  =       par;
            elseif af == 1
                val         =       real(eig(J));
                ew_reg      =       max(val(1:4));
                sol_reg     =       [orn;cbp;arg;utp;e2];
                parx_reg    =       par;
            end
            
        end %af
        %check for steady state, eigenvalues (stability)
        if SS_noreg < 1E-08 && SS_reg < 1E-08 && ew_noreg < -1E-06 && ew_reg < -1E-06
            %if true, exit
            value = 1;
            par_end_noreg(:,i1) =   parx_noreg;
            par_end_reg(:,i1)   =   parx_reg;
            sol1(:,i1)          =   sol_noreg;
            sol2(:,i1)          =   sol_reg;
        end
    end
end

%steady state solutions
par_end_reg(end-4:end,:) = sol2;
par_end_noreg(end-4:end,:) = sol1;

save('PAR.mat','par_end_noreg','par_end_reg')
end

function dcdt  =  odemodel(~,c,p,par,af,~)

par(end-4:end,1) = c;

%Parameters;
kcat1       =       par(find(strcmp(p,'kcat1')),1);
kcat2       =       par(find(strcmp(p,'kcat2')),1);
kcat3       =       par(find(strcmp(p,'kcat3')),1);
kcat4       =       par(find(strcmp(p,'kcat4')),1);
Km2         =       par(find(strcmp(p,'Km2')),1);
Kmorn       =       par(find(strcmp(p,'Kmorn')),1);
Kmcbp       =       par(find(strcmp(p,'Kmcbp')),1);
Km4         =       par(find(strcmp(p,'Km4')),1);
Kmu         =       par(find(strcmp(p,'Kmu')),1);
Kmu2        =       par(find(strcmp(p,'Kmu2')),1);
alpha1      =       par(find(strcmp(p,'alpha1')),1);
alpha2      =       par(find(strcmp(p,'alpha2')),1);
mumax       =       par(find(strcmp(p,'mumax')),1);
mumax2      =       par(find(strcmp(p,'mumax2')),1);

%Variables at t(end)
orn         =       par(end-4,1);
cbp         =       par(end-3,1);
arg         =       par(end-2,1);
utp         =       par(end-1,1);
e2          =       par(end,1);

%Definition of the growth rate
mue1        =       mumax * arg/(arg + Kmu);
mue2        =       mumax2 * utp/(utp + Kmu2);

%rates
r1 = kcat1;
if af == 0
    r2      =       kcat2 * e2;
elseif af == 1
    r2      =       kcat2 * e2 * orn^(Km2);
end
r3          =       kcat3 * 1/(1 + Kmorn*Kmcbp/(orn*cbp) + Kmorn/orn + Kmcbp/cbp);
r4          =       kcat4 * cbp/(cbp + Km4);
r5          =       alpha1 * mue1;
r6          =       alpha2 * mue2;

% mass balance
dcdt(1,1)   =       r1 - r3;
dcdt(2,1)   =       r2 - r3 - r4;
dcdt(3,1)   =       r3 - r5;
dcdt(4,1)   =       r4 - r6;
dcdt(5,1)   =       0;

end

