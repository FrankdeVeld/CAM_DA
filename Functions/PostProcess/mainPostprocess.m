function [] = mainPostprocess(outSim, validation, scenario, input, params)

t          = outSim.t_nodes*scenario.Tsc;
md_abs_err = sqrt(outSim.m_d)-validation.m_d;
md_err     = 100*abs(1-sqrt(outSim.m_d)./validation.m_d);

e2b = validation.e2b;
Pb = e2b*input.P*e2b';
showEllipseBplane(Pb,input.lim,outSim.rB,validation.rB, ...
                    input.metric_case,input.Lsc);

figure
plot(t,outSim.control)
legend('R','T','N')
grid on
xlabel('Time before TCA [s]')
ylabel('Normalized control')

figure
plot(t,outSim.deltaTca*scenario.Tsc)
grid on
xlabel('Time before TCA [s]')
ylabel('TCA shift [s]')

if params.metric_case == 2
    poc_abs_err = outSim.poc-validation.poc;
    poc_err = 100*abs(1-outSim.poc./validation.poc);
    smd_err = 100*abs(1-outSim.smd./validation.smd);
    % figure
    % plot(t,outSim.smd)
    % hold on
    % plot(t,validation.smd)
    % plot([t(1),t(end)],scenario.smdLim*ones(2,1),'k--')
    % hold off
    % grid on
    % xlabel('Time before TCA [s]')
    % ylabel('SMD')
    % legend('Optimization','Validation')

    figure
    plot(t,smd_err)
    grid on
    xlabel('Time before TCA [s]')
    ylabel('SMD relative error [$\%$]')
    % 
    figure
    semilogy(t,outSim.poc)
    hold on
    semilogy(t,validation.poc)
    ylim([1e-10,1])

    figure
    semilogy(t,poc_err)
    grid on
    xlabel('Time before TCA [s]')
    ylabel('PoC relative error [$\%$]')
else
    % figure
    % plot(t,sqrt(outSim.m_d)*input.Lsc)
    % hold on
    % plot(t,validation.m_d*input.Lsc)
    % plot([t(1),t(end)],scenario.md_lim*scenario.Lsc*ones(2,1),'k--')
    % hold off
    % grid on
    % xlabel('Time before TCA [s]')
    % ylabel('Miss distance [km]')
    % legend('Optimization','Validation')
    % 
    figure
    plot(t,md_err)
    grid on
    xlabel('Time before TCA [s]')
    ylabel('MD relative error [$\%$]')
end
dvTot = normOfVec(outSim.control(1:end-1,:)')*diff(-t)'*scenario.Asc*scenario.ctrlMax*1000; % m/s
disp(['Total Delta v = ', num2str(dvTot), ' m/s'])
disp(['Maneuver starts ' num2str(outSim.t_start*scenario.Tsc), ' s before TCA'])
if params.metric_case == 2
    disp(['Maximum PoC absolute error ' num2str(max(abs(poc_abs_err)))])
    disp(['Maximum SMD relative error = ' num2str(max(smd_err)), ' %'])
else
    disp(['Maximum miss distance absolute error ' num2str(max(abs(md_abs_err))*scenario.Lsc*1e3), ' m'])
    disp(['Maximum miss distance relative error ' num2str(max(md_err)), ' %'])
end
end

