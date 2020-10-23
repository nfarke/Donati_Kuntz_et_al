function ppc_kinetic_curve

Kmal = 0.23;
mal_wt = 1.68;
mal_new = 0.1050;
mal = 0.01:0.01:10;

Kasp = 0.31;
asp_wt = 4.34;
asp_new = 0.3338;
asp = 0.01:0.01:10;


[ASP,MAL] = meshgrid(asp,mal);

Term_asp = Kasp./(ASP + Kasp);
Term_mal = Kmal./(MAL + Kmal);
Term = Term_asp.*Term_mal;
c = ASP.*MAL;

CO(:,:,1) = zeros(1000); % red
CO(:,:,2) = ones(1000).*linspace(0,1,1000); % green
CO(:,:,3) = ones(1000).*linspace(0,1,1000); % blue

%plot
contour(ASP,MAL,Term,500)
set(gca,'xscale','log');
set(gca,'yscale','log');
colorbar()

hold on
scatter(4.34,1.68,'k','filled')
hold on
scatter(0.33,0.11,'b','filled')
hold on
scatter(0.31,0.01,'r','filled') %Ki value
hold on
scatter(0.01,0.23,'r','filled') %Ki value

xlabel('L-aspartate, mM');
ylabel('Malate, mM');
end

