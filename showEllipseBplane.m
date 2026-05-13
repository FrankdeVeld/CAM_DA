function [] = showEllipseBplane(PB,lim,traj,trajVal,metric,Lsc)
%relTrajEllipsoids plots the relative trajectories with ellipsoid
%visualization
%
% Author: Zeno Pavanello, 2022
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

figure()
plot(traj(:,1)*Lsc,traj(:,2)*Lsc,'k')
hold on
plot(trajVal(end,1)*Lsc,trajVal(end,2)*Lsc,'ko','HandleVisibility','off')
plot(0,0,'k*','HandleVisibility','off')
plot(trajVal(:,1)*Lsc,trajVal(:,2)*Lsc,'r')
t          = 0:0.001:2*pi;
if metric == 1
    plot(sqrt(lim)*Lsc*sin(t),sqrt(lim)*Lsc*cos(t),'k','HandleVisibility','off')
else
    [semiaxes,cov2b] = defineEllipsoid(PB,lim);
    a          = semiaxes(1)*Lsc;
    b          = semiaxes(2)*Lsc;
    x          = a*cos(t);
    y          = b*sin(t);
    ellCov     = [x; y];
    ellB       = nan(2,length(t));
    for k = 1:length(t)
        ellB(:,k) = cov2b*ellCov(:,k);
    end
    plot(ellB(1,:),ellB(2,:),'k','HandleVisibility','off');
end
grid on
box on
xlabel('\xi [km]'); ylabel('\zeta [km]')
legend('Optimisation','Validation')
axis equal
hold off
end