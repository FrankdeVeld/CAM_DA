FileName = 'N50_Euc';

xp     = load(['.\write_read\xp_',FileName,'.dat']);
xpadj  = load(['.\write_read\xpadj_',FileName,'.dat']);
xs     = load(['.\write_read\xs_',FileName,'.dat']);
xs_tf_new = load(['.\write_read\xsfnew_',FileName,'.dat']); % with adjusted tca
DM     = load(['.\write_read\DM_',FileName,'.dat']);
DeltaRB = load(['.\write_read\DeltaRB_',FileName,'.dat']);

xpn        = load(['.\write_read\xnp1_Val_',FileName,'.dat']);
xfull_Val  = load(['.\write_read\xfull_Val_',FileName,'.dat']);
xfull_ValR  = load(['.\write_read\xfull_ValR_',FileName,'.dat']);
xfull_ValT  = load(['.\write_read\xfull_ValT_',FileName,'.dat']);

%% Validating whether xp_tnp1 is correctly calculated
Disxpadjxpn = sqrt( (xpadj(:,1) - xpn(2:end,1)).^2 +  (xpadj(:,2) - xpn(2:end,2)).^2 +  (xpadj(:,3) - xpn(2:end,3)).^2 );
t = [1:1:length(Disxpadjxpn)]';
figure()
grid on
plot(t,Disxpadjxpn)
xlabel('Node number (-)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Error (km)', 'FontSize', 12, 'FontWeight', 'bold');

title('Error between adjusted states xnp1, DA and validation')

%% Validating whether the distance metric tracks the right thing
DM_Val       = sqrt( (xfull_Val(end,1) - xs_tf_new(1)).^2 +  (xfull_Val(end,2) - xs_tf_new(2)).^2 +  (xfull_Val(end,3) - xs_tf_new(3)).^2 );
DM_DeltaRB   = sqrt( (DeltaRB(1,1)).^2 + (DeltaRB(1,2)).^2 +  (DeltaRB(1,3)).^2);
DM_Final     = DM(1);

fprintf('Difference DeltaRB DA propagation and validation: %.4f\n', abs(DM_Val-DM_DeltaRB));
fprintf('Difference DM DA propagation and validation: %.4f\n', abs(DM_Val-DM_Final));

%% Distance metric evaluation - figure doesn't say much
DM_Nom = sqrt((xp(end,1)-xs(end,1)).^2 + (xp(end,2)-xs(end,2)).^2 + (xp(end,3)-xs(end,3)).^2 ) ;
fprintf('Nominal final Distance Metric: %.4f\n', DM_Nom);
DM_Val = sqrt((xfull_Val(end,1)-xs_tf_new(1)).^2 + (xfull_Val(end,2)-xs_tf_new(2)).^2 + (xfull_Val(end,3)-xs_tf_new(3)).^2 ) ;
fprintf('First-order control final Distance Metric: %.4f\n', DM_Val);
DM_ValR = sqrt((xfull_ValR(end,1)-xs_tf_new(1)).^2 + (xfull_ValR(end,2)-xs_tf_new(2)).^2 + (xfull_ValR(end,3)-xs_tf_new(3)).^2 ) ;
fprintf('Pure R thrust Distance Metric: %.4f\n', DM_ValR);
DM_ValT = sqrt((xfull_ValT(end,1)-xs_tf_new(1)).^2 + (xfull_ValT(end,2)-xs_tf_new(2)).^2 + (xfull_ValT(end,3)-xs_tf_new(3)).^2 ) ;
fprintf('Pure T final Distance Metric: %.4f\n', DM_ValT);

%% Check if the assessment is truly at the closest approach
IP_Nom = dot(xp(end,1:3)+DeltaRB(end,1:3)-xs(end,1:3),xp(end,4:6)-xs(end,4:6));
fprintf('Nominal final inner product: %.4f\n', IP_Nom);
IP_Val = dot(xfull_Val(end,1:3)-xs_tf_new(1:3),xfull_Val(end,4:6)-xs_tf_new(4:6));
fprintf('First-order control final inner product: %.4f\n', IP_Val);
IP_ValR = dot(xfull_ValR(end,1:3)-xs_tf_new(1:3),xfull_ValR(end,4:6)-xs_tf_new(4:6));
fprintf('Pure R thrust control final inner product: %.4f\n', IP_ValR);
IP_ValT = dot(xfull_ValT(end,1:3)-xs_tf_new(1:3),xfull_ValT(end,4:6)-xs_tf_new(4:6));
fprintf('Pure T thrust control final inner product: %.4f\n', IP_ValT);
% If not zero: the comparison is not at the closest approach