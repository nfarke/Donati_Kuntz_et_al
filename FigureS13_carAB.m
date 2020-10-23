
function carAB_kinetic_curve
Kmorn = 0.37;
orn_wt = 0.01;
orn_new = 5;
orn = 0.0001:0.0001:10;
Term = orn./(orn + Kmorn);
plot(orn,Term,'k')
hold on
Term2 = 0.01/(0.01 +Kmorn);
Term3 = 5/(5+Kmorn);
scatter(0.01,Term2,'b','filled');
hold on
scatter(5,Term3,'r','filled');
hold on

%plot km
yval = [0,0.5];
xval = [Kmorn,Kmorn];
plot(xval,yval,'k--');

yval = [0.5,0.5];
xval = [0,Kmorn];
plot(xval,yval,'k--');

xlabel('orn, mM')
ylabel('% activation')
end