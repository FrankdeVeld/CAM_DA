function [] = plotBplane(Pp,Ps,xp,xs,smdLim)
%    plotBplane plots the B-plane configuration
% 
% INPUT: 
%        Pp:     (3,3) double, covariance matrix of the primary, ECI
%        Ps:     (3,3) double, covariance matrix of the secondary, ECI
%        xp:     (6,1) double, state of the primary, ECI
%        xs:     (6,1) double, state of the secondary, ECI
%        smdLim: (1,1) SMD (or MD) threshold
%        
% OUTPUT:
%        
% Author: Zeno Pavanello, 2024
% E-mail: zpav176@aucklanduni.ac.nz
%--------------------------------------------------------------------------

r2ep   = rtn2eci(xp(1:3),xp(4:6));
r2es   = rtn2eci(xs(1:3),xs(4:6));
P      = r2ep*Pp*r2ep' + r2es*Ps*r2es';
for i = 1:length(xp)
    [PB,p(:,i),smd(i)] = Bplane(xp(i,:)',xs(i,:)',P);
end
[semiaxes,cov2b] = defineEllipsoid(PB,smdLim);
a          = semiaxes(1);
b          = semiaxes(2);
tt         = 0:0.001:2*pi;
xx         = a*cos(tt);
yy         = b*sin(tt);
ellCov     = [xx; yy];
ellB       = nan(2,length(tt));
for j = 1:length(tt)
    ellB(:,j) = cov2b*ellCov(:,j);
end
figure()
hold on    
plot(ellB(2,:),ellB(1,:),'k');
plot(p(2,:),p(1,:),'LineWidth',2);
grid on 
xlabel('$\zeta$ [km]','Interpreter','latex')
ylabel('$\xi$ [km]','Interpreter','latex')
hold off
axis equal
box on
end

