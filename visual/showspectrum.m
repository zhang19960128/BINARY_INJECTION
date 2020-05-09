data=importdata("injection.dat");
xgrid=data(1:end,1);
xxz=data(1:end,2)*10^(6);
yxz=data(1:end,3)*10^(6);
zxz=data(1:end,4)*10^(6);
figure('units','normalized','outerposition',[0 0 1 1])
plot(xgrid,xxz,'-r+',xgrid,yxz,'-kp',xgrid,zxz,'-m+','LineWidth',2)
legend("xxz","yxz","zxz")
xlabel("Photon Energy(ev)");
xlim([0.2,1.2])
ylabel("Susceptibility($\mu$A/$V^{2}$)",'Interpreter','latex')
set(gca,"FontSize",26)
saveas(gcf,'RhSi.eps','epsc')
