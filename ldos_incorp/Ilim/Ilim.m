clear all

cd '/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/ldos_incorp/Ilim'

iv_file = 'iv_data.txt';
T = readtable(iv_file);
a_s = T{:,1};
i_lim = T{:,3};
[a_s, ids] = sort(a_s);
i_lim = i_lim(ids);


load('Ilim_exp.mat')
Iexp = Ilim_exp(:,2);
as_exp = Ilim_exp(:,1);

figure()
plot(a_s, i_lim*1e12 ,'-','LineWidth',2, 'MarkerSize', 9);
hold on
plot(as_exp, Iexp ,'o','LineWidth',2, 'MarkerSize', 9);
set(gcf,'Position',[358,301,570,477])
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1,'GridLineStyle','--');
%fig.MarkerSize = 9;
xlabel('Radius (nm)','FontSize',24,'interpreter', 'latex') 
ylabel('Current ($\textnormal{p\AA}$)','FontSize',24,'interpreter', 'latex')
legend('Sim','Exp','Location','best','Interpreter','latex')
legend('boxoff') 
saveas(gcf,'fig1b.png')
hold off

%% Fit

ft = fittype( 'a*x + b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.0 0.0];

[fit_sim, gof] = fit( a_s, i_lim*1e12, ft, opts )

[fit_exp, gof] = fit( as_exp, Iexp, ft, opts )


fit_sim(5)
fit_exp(5)