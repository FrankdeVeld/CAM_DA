function [] = mainPostprocess(simOutput, validation, scenario, input, params)

t = scenario.t;

figure
% plot(tca_shift*scenario.Tsc)
% hold on
% plot(dtca*Tsc)
% hold off
showEllipseBplane(scenario.Pb,input.lim,simOutput.rB,validation.rB, ...
                    input.metric_case,input.Lsc);

figure
plot(t,simOutput.control)
legend('R','T','N')
grid on
xlabel('Orbits before TCA')
ylabel('Normalized control')

figure
plot(t,sqrt(simOutput.m_d)*input.Lsc)
hold on
plot(t,validation.m_d*input.Lsc)
plot([t(1),t(end)],scenario.md_lim*scenario.Lsc*ones(2,1),'k--')
hold off
grid on
xlabel('Orbits before TCA')
ylabel('Miss distance [km]')
legend('Optimization','Validation')

figure
plot(t,1-simOutput.m_d/validation.m_d)
grid on
xlabel('Orbits before TCA')
ylabel('MD relative error [%]')


if params.metric_case == 2
    figure
    plot(t,simOutput.smd)
    hold on
    plot(t,validation.smd)
    plot([t(1),t(end)],scenario.smdLim*ones(2,1),'k--')
    hold off
    grid on
    xlabel('Orbits before TCA')
    ylabel('SMD')
    legend('Optimization','Validation')

    figure
    plot(t,1-simOutput.smd./validation.smd)
    grid on
    xlabel('Orbits before TCA')
    ylabel('SMD relative error [%]')

    figure
    semilogy(t,validation.poc)
    ylim([1e-10,1e-2])
end
end