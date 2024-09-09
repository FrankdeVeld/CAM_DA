FileName = 'N100_Euc';

xp     = load(['.\write_read\xp_',FileName,'.dat']);
xpadj  = load(['.\write_read\xpadj_',FileName,'.dat']);
xs     = load(['.\write_read\xs_',FileName,'.dat']);
DM     = load(['.\write_read\DM_',FileName,'.dat']);

xpn        = load(['.\write_read\xnp1_Val_',FileName,'.dat']);
xfull_Val  = load(['.\write_read\xfull_Val_',FileName,'.dat']);
xfull_ValR  = load(['.\write_read\xfull_ValR_',FileName,'.dat']);
xfull_ValT  = load(['.\write_read\xfull_ValT_',FileName,'.dat']);

%% 
Disxpadjxpn = sqrt( (xpadj(:,1) - xpn(2:end,1)).^2 +  (xpadj(:,2) - xpn(2:end,2)).^2 +  (xpadj(:,3) - xpn(2:end,3)).^2 );
t = [1:1:length(Disxpadjxpn)]';
figure()
grid on
plot(t,Disxpadjxpn)
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Distance between states (km)', 'FontSize', 12, 'FontWeight', 'bold');

title('Distance xnp1 (DA) and validation propagation (N=100)')

%% 
DisFVal = sqrt( (xfull_Val(end,1) - xs(end,1)).^2 +  (xfull_Val(end,2) - xs(end,2)).^2 +  (xfull_Val(end,3) - xs(end,3)).^2 );
DisPrimVal = sqrt( (xfull_Val(end,1) - xp(end,1)).^2 +  (xfull_Val(end,2) - xp(end,2)).^2 +  (xfull_Val(end,3) - xp(end,3)).^2 );
InitDis = sqrt( (xs(end,1) - xp(end,1)).^2 + (xs(end,2) - xp(end,2)).^2 +  (xs(end,3) - xp(end,3)).^2 );

DMFinal = DM(1);

error = abs(DisFVal - DMFinal)

%% 
norm_nom = sqrt(xp(:,1).^2 + xp(:,1).^2 + xp(:,1).^2 ) ;
norm_Val = sqrt(xfull_Val(:,1).^2 + xfull_Val(:,1).^2 + xfull_Val(:,1).^2 ) ;
norm_ValR = sqrt(xfull_ValR(:,1).^2 + xfull_ValR(:,1).^2 + xfull_ValR(:,1).^2 ) ;
norm_ValT = sqrt(xfull_ValT(:,1).^2 + xfull_ValT(:,1).^2 + xfull_ValT(:,1).^2 ) ;

t = [1:1:length(norm_nom)]';

figure()
grid on
hold on
plot(t,norm_Val-norm_nom, 'LineWidth', 2)
plot(t,norm_ValR-norm_nom, 'LineWidth', 2)
plot(t,norm_ValT-norm_nom, 'LineWidth', 2)
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Norm difference with nominal [km]', 'FontSize', 12, 'FontWeight', 'bold');
legend('Control','Pure Radial','Pure Tangential', 'Location', 'northwest')
title('Norm difference with primary for various control strategies')

% Save the figure as a PNG file
saveas(gcf, ['ControlVal', FileName, '.png']);

% Save the figure as a FIG file
saveas(gcf, ['ControlVal', FileName, '.fig']);

%%
dot(xfull_Val(end,1:3)-xs(end,1:3),xfull_Val(end,4:6)-xs(end,4:6))
dot(xp(end,1:3)-xs(end,1:3),xp(end,4:6)-xs(end,4:6))