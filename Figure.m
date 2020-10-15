%this file takes the steady state data and simulates the knock-down
%of carAB over a time window of tspan.

function Figure
    tspan   = 0:1:100;
    opts    = odeset('RelTol',1e-07,'AbsTol',1E-7);
    load PAR %load parameter sets and concentrations
    ParSize = size(par_end_reg,2);
    
    parfor i1 = 1:ParSize %parfor is possible
        
        disp(i1) %show iteration
        [p,~] = Sample(1,1);%get parameter names

        for af =0:1  %0 = no reg, 1= reg
           
           if af == 0
           par      =     par_end_noreg(:,i1);  
           y0       =     par(end-4:end);
           else
           par      =     par_end_reg(:,i1);   
           y0       =     par(end-4:end);
           end
           
           %solve differential equation
           [~,y]    =     ode23s(@(t,c) odemodel(t,c,p,par,af),tspan,y0,opts);

           %Assign parameter values and concentrations of variables
           kcat2    =     par(find(strcmp(p,'kcat2')),1);
           Km2      =     par(find(strcmp(p,'Km2')),1);
           orn      =     y(:,1);
           e2       =     y(:,5);

           if af == 0 %Results with no carAB regulation
                R2_noreg(:,i1)      =      kcat2 * e2; %r2 flux
                orn_no_reg(:,i1)    =      y(:,1);
                cbp_no_reg(:,i1)    =      y(:,2);
                arg_no_reg(:,i1)    =      y(:,3);
                utp_no_reg(:,i1)    =      y(:,4);
           elseif af == 1 %Results with carAB Regulation
                R2_reg(:,i1)        =      kcat2.*e2.*power(orn,Km2); %r2 flux
                orn_reg(:,i1)       =      y(:,1);
                cbp_reg(:,i1)       =      y(:,2);
                arg_reg(:,i1)       =      y(:,3);
                utp_reg(:,i1)       =      y(:,4);
           end
        end      
    end
    
    %create figures
    %% Figure 1
    figure(1)
    subplot(1,2,1)
    title('no_reg')
    h1      =   plot(orn_no_reg,'k');
    hold on
    plot(mean(orn_no_reg,2),'c');
    hold on
    h2      =   plot(arg_no_reg,'b');
    hold on
    plot(mean(arg_no_reg,2),'c');
    hold on
    h3      =	plot(cbp_no_reg,'m');
    hold on
    plot(mean(cbp_no_reg,2),'c');
    hold on
    h4      =   plot(utp_no_reg,'r');
    hold on
    plot(mean(utp_no_reg,2),'c');
    xlim([0 100])
    ylim([0 4])
    title('dysregulated model');

    legend([h1(1) h2(1) h3(1) h4(1)],'orn','arg','cbp','utp')
    ylabel('metabolite concentration')
    xlabel('time, min')

    subplot(1,2,2)
    title('reg')
    h1      =   plot(orn_reg,'k');
    hold on
    plot(mean(orn_reg,2),'c');
    hold on
    h2      =   plot(arg_reg,'b');
    hold on
    plot(mean(arg_reg,2),'c');
    hold on
    h3      =   plot(cbp_reg,'m');
    hold on
    plot(mean(cbp_reg,2),'c');
    hold on
    h4      =   plot(utp_reg,'r');
    hold on
    plot(mean(utp_reg,2),'c');

    ylim([0 4])
    xlim([0 100])
    legend([h1(1) h2(1) h3(1) h4(1)],'orn','arg','cbp','utp')
    xlabel('time, min')
    title('regulated model');

    %% Figure 2
    figure(2)
    subplot(1,2,1)
    plot(R2_noreg/1.425);
    xlim([0 100])
    ylim([0.5 1])
    xlabel('time, min')
    ylabel('Relative carAB flux')
    
    subplot(1,2,2)
    plot(R2_reg/1.425);
    xlim([0 100])
    ylim([0.5 1])
    ylabel('Relative carAB flux')
    xlabel('time, min')
end

function dcdt  =  odemodel(~,c,p,par,af,~)

par(end-4:end,1) =  c;

%Parameters;
kcat1     =    par(find(strcmp(p,'kcat1')),1);
kcat2     =    par(find(strcmp(p,'kcat2')),1);
kcat3     =    par(find(strcmp(p,'kcat3')),1);
kcat4     =    par(find(strcmp(p,'kcat4')),1);
Km2       =    par(find(strcmp(p,'Km2')),1);
Kmorn     =    par(find(strcmp(p,'Kmorn')),1);
Kmcbp     =    par(find(strcmp(p,'Kmcbp')),1);
Km4       =    par(find(strcmp(p,'Km4')),1);
Kmu       =    par(find(strcmp(p,'Kmu')),1);
Kmu2      =    par(find(strcmp(p,'Kmu2')),1);
alpha1    =    par(find(strcmp(p,'alpha1')),1);
alpha2    =    par(find(strcmp(p,'alpha2')),1);
mumax     =    par(find(strcmp(p,'mumax')),1);
mumax2    =    par(find(strcmp(p,'mumax2')),1);

%Variables at t(end)
orn       =    par(end-4,1);
cbp       =    par(end-3,1);
arg       =    par(end-2,1);
utp       =    par(end-1,1);
e2        =    par(end,1);

%Definition of the growth rate
mue1      =    mumax * arg/(arg + Kmu);
mue2      =    mumax2 * utp/(utp + Kmu2);
mue       =    mean([mue1 mue2]);

%rates
r1        =    kcat1;
if     af == 0
      r2  =    kcat2 * e2;
elseif af == 1
      r2  =    kcat2 * e2 * orn^Km2;
end
r3        =    kcat3 * 1/(1 + Kmorn*Kmcbp/(orn*cbp) + Kmorn/orn + Kmcbp/cbp);
r4        =    kcat4 * cbp/(cbp + Km4);
r5        =    alpha1 * mue1;
r6        =    alpha2 * mue2;

% mass balance
dcdt(1,1) =    r1 - r3;
dcdt(2,1) =    r2 - r3 - r4;
dcdt(3,1) =    r3 - r5;
dcdt(4,1) =    r4 - r6;
dcdt(5,1) =    -mue*e2; %dilution of e2 by growth simulates CRISPRi
end
