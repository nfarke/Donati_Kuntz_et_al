function [par]=calc_assign_par(p,par,flux_ss,af)
        Km2   = par(find(strcmp(p,'Km2')),1);
        Kmorn = par(find(strcmp(p,'Kmorn')),1);
        Kmcbp = par(find(strcmp(p,'Kmcbp')),1);
        Km4   = par(find(strcmp(p,'Km4')),1);
        Kmu   = par(find(strcmp(p,'Kmu')),1);
        Kmu2  = par(find(strcmp(p,'Kmu2')),1);

        orn   =  par(end-4);
        cbp   =  par(end-3);
        arg   =  par(end-2);
        utp   =  par(end-1);
        e2    =  par(end);
        
        %calculate and assign kcat1
        par(find(strcmp(p,'kcat1')),1) = flux_ss(1);

        %calculate and assign kcat2
        if af == 0
            kcat2 = flux_ss(2)/(e2);
        elseif af == 1
            kcat2 = flux_ss(2)/(e2*orn^(Km2));
        end
        par(find(strcmp(p,'kcat2')),1) = kcat2;

        %calculate and assign kcat3
        kcat3 = flux_ss(3)*(1 + Kmorn*Kmcbp/(orn*cbp) + Kmcbp/cbp + Kmorn/orn);
        par(find(strcmp(p,'kcat3')),1) = kcat3; 

        %calculate and assign kcat4
        kcat4 = flux_ss(4)/(cbp/(cbp + Km4));
        par(find(strcmp(p,'kcat4')),1) = kcat4; 

        %calculate and assign mumax
        mumax = 0.01/(arg/(arg + Kmu));
        par(find(strcmp(p,'mumax')),1) = mumax;            

        %calculate and assign mumax2
        mumax2 = 0.01/(utp/(utp + Kmu2));
        par(find(strcmp(p,'mumax2')),1) = mumax2;
end
