%clear all;

cd '/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/ldos_incorp/'


%% Plot k rates

load('2.0/k_data_Vapp_0.0.mat');

theta_list = theta;
knum = 24;
theta = string(theta_list);
nsamps = length(rscx);
%Vappl=0;
directory = './';
%load(append(directory,'ldos-th_',theta,'_nsamps_',string(nsamps),'_knum_',string(knum),'.mat')) %Original ldos data (for rscx and rscy)

alpha = 2.47;
sc_alpha = alpha/(2*sind(theta_list(1)/2));



figure(1)
surf(sc_alpha*rscx, sc_alpha*rscy, kox_data)
shading flat
view(2)
axis equal
colorbar
title(['$k_{ox}$ for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
xlabel('$x$ (SC units)','FontSize',20,'interpreter', 'latex')
ylabel('$y$ (SC units)','FontSize',20,'interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',20);

figure(2)
surf(sc_alpha*rscx, sc_alpha*rscy, kred_data)
shading flat
view(2)
axis equal
colorbar
title(['$k_{red}$ for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
xlabel('$x$ (SC units)','FontSize',20,'interpreter', 'latex')
ylabel('$y$ (SC units)','FontSize',20,'interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',20);

figure(3)
surf(sc_alpha*rscx, sc_alpha*rscy, Vdl_data)
shading flat
view(2)
axis equal
colorbar
title(['$V_{dl}$ for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
xlabel('$x$ (SC units)','FontSize',20,'interpreter', 'latex')
ylabel('$y$ (SC units)','FontSize',20,'interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',20);

%% Plot ldos

tar_E = 5500;

figure(4)
surf(sc_alpha*rscx, sc_alpha*rscy, squeeze(data(tar_E,:,:)))
shading flat
view(2)
axis equal
colorbar
%caxis([2.4e-5 15.7e-5])
title(['Local DOS for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
xlabel('$x$ (SC units)','FontSize',20,'interpreter', 'latex')
ylabel('$y$ (SC units)','FontSize',20,'interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize',20);

%% Current map scg

clear all

% Parameters to tweak before running : idx_abn (l:103) and caxis(l:191)


% Load domain areas
load('/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/prefactor/area_domains.mat')
% Load I-V data
iv_file = append('/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/ldos_incorp/Current_map/', ...
'5nm/eta_0.0/1.6/iv_data.txt') ; %6mM_Ru3+/
spl = split(iv_file, '/');

T = readtable(iv_file);
I_arr = T{:,2};
c_x = T{:,4}; % in Ang
c_y = T{:,5}; % in Ang

% Interpolate zero values where solution did not converge
non_z_el = find(I_arr);
z_el = find(~I_arr);
interpolant = scatteredInterpolant(c_x(non_z_el),c_y(non_z_el),I_arr(non_z_el));
for i=1:length(z_el)
   I_arr(z_el(i)) = interpolant(c_x(i), c_y(i)); 
end

idx_ab = 785;
idx_aa = 1;

%Reduce overshooting of currents in map
idx_abn = idx_aa;
z_shoot = find(abs(I_arr) > 1.05*abs(I_arr(idx_abn)));
z_normal = find(abs(I_arr) < 1.05*abs(I_arr(idx_abn)));
interpolant = scatteredInterpolant(c_x(z_normal),c_y(z_normal),I_arr(z_normal));
for i=1:length(z_shoot)
   I_arr(z_shoot(i)) = interpolant(c_x(i), c_y(i)); 
end

% Anomalous points cleaned up in next section

z_ana = find(abs(I_arr) > 1.05*abs(I_arr(idx_abn)));
z_normal = find(abs(I_arr) < 1.05*abs(I_arr(idx_abn)));
interpolant = scatteredInterpolant(c_x(z_normal),c_y(z_normal),I_arr(z_normal));
for i=1:length(z_ana)
   I_arr(z_ana(i)) = interpolant(c_x(i), c_y(i)); 
end

nr = sqrt(length(c_x));
I = reshape(I_arr, [nr, nr]);
X = reshape(c_x, [nr, nr]);
Y = reshape(c_y, [nr, nr]);

alpha = 2.47;
theta = str2double(string(spl(end-1)));
sc_alpha = alpha/(2*sind(theta/2)); % in Ang
sc_area = (sc_alpha^2)*sind(60)*1e-2; % in nm^2
moir_a = 60; %deg
bvec1 = [sc_alpha*cosd(moir_a/2), sc_alpha*sind(moir_a/2)].*0.1; % basis vector 1 in nm                                                   
bvec2 = [sc_alpha*cosd(moir_a/2), -sc_alpha*sind(moir_a/2)].*0.1; % basis vector 2 in nm
center_aa = [0,0];
center_ab = [2*bvec1(1)/3, 0];
center_sp = [bvec1(1), 0];
[~, idx_aa] = min(abs(vecnorm([c_x c_y]*0.1 - center_aa, 2, 2)));
[~, idx_ab] = min(abs(vecnorm([c_x c_y]*0.1 - center_ab, 2, 2)));

circle_orig = center_aa; %max(c_x)/2
circle_rad = str2double(string(regexp(string(spl(12)),'\d*','Match'))); % in nm
th = 0:pi/50:2*pi;
xunit = circle_rad * cos(th) + circle_orig(1);
yunit = circle_rad * sin(th) + circle_orig(2);

dx = linspace(0,1,nr+1);
dx = dx(1:end-1);
[x,y] = meshgrid(dx,dx);

A = [sqrt(3)/2 sqrt(3)/2; -0.5 0.5];
r_sc = A(:,1).*x(:)' + A(:,2).*y(:)';

rx = reshape(r_sc(1,:),nr,nr);
ry = reshape(r_sc(2,:),nr,nr);

scg = 2; % supercell grid scaling
nrsc = (2*scg+1)*nr;
rscx = zeros(nrsc,nrsc);
rscy = rscx;
I_sc = rscx;

for scx = -scg:scg
    for scy = -scg:scg
        tar_idx_x = [1:nr]+nr*(scx+scg);
        tar_idx_y = [1:nr]+nr*(scy+scg);
        sc_h = A(:,1)*scy + A(:,2)*scx;
        
        rscx(tar_idx_x,tar_idx_y) = rx + sc_h(1);
        rscy(tar_idx_x,tar_idx_y) = ry + sc_h(2);
        I_sc(tar_idx_x,tar_idx_y) = I;
        

    end
end

%figure('Renderer','Painters');
figure()
surf(sc_alpha*rscx.*0.1, sc_alpha*rscy.*0.1, I_sc*1e15) % in nm
shading flat
hold on
plot(xunit, yunit, 'r-', 'LineWidth', 2);
view(2)
axis equal
colorbar
%title(['$k_{red}$ for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
xlabel('$x$ (nm)','FontSize',20,'interpreter', 'latex')
ylabel('$y$ (nm)','FontSize',20,'interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize', 20,'LineWidth', 1.4);
%shading(gca,'interp')
c = colorbar;
%w = c.LineWidth;
c.LineWidth = 1.5;
caxis([-70, -50]) %caxis([-75, -40]) %For 5nm, 1.1 deg -> caxis([-75, -40])% For 5nm, 2 deg -> caxis([-75 -60])
box on
grid off
%saveas(gcf,'/Users/mbabar/Desktop/PNP_plots/fig5a.png')
%saveas(gcf,'/Users/mbabar/Desktop/PNP_plots/fig5a.svg');
hold off

%-100*(max(I_arr) - min(I_arr))./min(I_arr)

AA_area = sc_area * 0.01*interp1(angle1, AA, theta);

disp(append('Theta : ', string(theta)), ' deg')
disp(append('Moire sc area: ', string(sc_area), ' nm^2, AA area: ', string(AA_area), ' nm^2'))
disp(append('Basis vector length: ', string(norm(bvec1)), ' nm'))
disp(append('Current AA/AB diff: ', string((I_arr(idx_aa) - I_arr(idx_ab))*1e15), ' fA'))
disp(append('AA current : ', string(I_arr(idx_aa)*1e15), ' fA', ', AB current', string(I_arr(idx_ab)*1e15), ' fA'))

%% Clean current map (contd.)

% Need variables from last section's run
pos_mat = [c_x, c_y]*0.1 ;
I_temp = I_arr;

c1 = [0,0];
c2 = bvec1;
c3 = bvec2;
c4 = bvec1+bvec2;
idcs = ones(length(I_arr), 1);

for i=1:length(pos_mat)
    if norm(pos_mat(i,:) - c1) < 5 % 5 chosen same a_s = 5nm, tweak for different a_s
        %I_temp(i) = 0;
    elseif norm(pos_mat(i,:) - c2) < 5
        %I_temp(i) = 0;
    elseif norm(pos_mat(i,:) - c3) < 5
        %I_temp(i) = 0;
    elseif norm(pos_mat(i,:) - c4) < 5
        %I_temp(i) = 0;
    else
        idcs(i) = 0;
    end
end
idcs = logical(idcs);
I_ab = I_arr(~idcs);
I_ab(abs(I_ab*1e15) > 60) = 0.7*I_ab(abs(I_ab*1e15) > 60);
I_temp(~idcs) = I_ab;

I_aa = I_arr(idcs);
I_aa(abs(I_aa*1e15) > 70) = 0.8*I_aa(abs(I_aa*1e15) > 70);
I_temp(idcs) = I_aa;

I = reshape(I_temp, [nr, nr]);


A = [sqrt(3)/2 sqrt(3)/2; -0.5 0.5];
r_sc = A(:,1).*x(:)' + A(:,2).*y(:)';

rx = reshape(r_sc(1,:),nr,nr);
ry = reshape(r_sc(2,:),nr,nr);
scg = 2; % supercell grid scaling
nrsc = (2*scg+1)*nr;
rscx = zeros(nrsc,nrsc);
rscy = rscx;
I_sc = rscx;

for scx = -scg:scg
    for scy = -scg:scg
        tar_idx_x = [1:nr]+nr*(scx+scg);
        tar_idx_y = [1:nr]+nr*(scy+scg);
        sc_h = A(:,1)*scy + A(:,2)*scx;
        
        rscx(tar_idx_x,tar_idx_y) = rx + sc_h(1);
        rscy(tar_idx_x,tar_idx_y) = ry + sc_h(2);
        I_sc(tar_idx_x,tar_idx_y) = I;
        

    end
end

%figure('Renderer','Painters');
figure()
surf(sc_alpha*rscx.*0.1, sc_alpha*rscy.*0.1, I_sc*1e15) % in nm
shading flat
hold on
plot(xunit, yunit, 'k-', 'LineWidth', 2);
view(2)
axis equal
colorbar
%title(['$k_{red}$ for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
xlabel('$x$ (nm)','FontSize',20,'interpreter', 'latex')
ylabel('$y$ (nm)','FontSize',20,'interpreter', 'latex')
set(gca,'FontName','Times New Roman','FontSize', 20,'LineWidth', 1);
%shading(gca,'interp')
c = colorbar;
%w = c.LineWidth;
c.LineWidth = 1.5;
%caxis([-75, -40]) %For 5nm, 1.1 deg -> caxis([-75, -40])% For 5nm, 2 deg -> caxis([-75 -60])
box on
grid off
%saveas(gcf,'/Users/mbabar/Desktop/PNP_plots/fig5b.png')
%saveas(gcf,'/Users/mbabar/Desktop/PNP_plots/fig5b.svg');
hold off

%% AA/AB CVs

%clear all
iv_file = '/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/ldos_incorp/CV/5nm/2mM/1.6/iv_data.txt' ; 

T = readtable(iv_file);
domain_arr = T{:,1};
I_arr = T{:,2};
Vapp_arr = T{:,3};

ind_aa = find(string(domain_arr) == "AA");
ind_ab = find(string(domain_arr) == "AB");
I_aa = I_arr(ind_aa);
I_ab = I_arr(ind_ab);
V_aa = Vapp_arr(ind_aa);
V_ab = Vapp_arr(ind_ab);
[V_aa, sort_aa] = sort(V_aa);
[V_ab, sort_ab] = sort(V_ab);
I_aa = I_aa(sort_aa);
I_ab = I_ab(sort_ab);

figure(1)
plot(V_aa, I_aa*1e15, 'o-', 'LineWidth', 2)
hold on
plot(V_ab, I_ab*1e15, 'o-', 'LineWidth', 2)
%plot(V_aa0, I_aa0*1e15, 'o-', 'LineWidth', 2)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.4,'GridLineStyle','--'); % ,'yscale','log'
xlabel('$V_{app}$ (V)','interpreter','latex')
ylabel('$I$ (fA)','interpreter','latex')
ylim([min(I_arr*1e15), 0])
%title(['$k_{red}$ for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
lh = legend('AA','AB','Location','best','FontSize',20,'interpreter','latex');
legend boxoff
hold off


%% Fit CVs with Sigmoid (contd.)

ft = fittype( 'c - (a/(1+exp(b*(x-d))))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [170 0.485375648722841 0.0 0.0];

[fit_aa, gof] = fit( V_aa, I_aa*1e15, ft, opts)

[fit_ab, gof] = fit( V_ab, I_ab*1e15, ft, opts)

figure()
plot(V_aa, I_aa*1e15, 'ro', 'LineWidth', 2)
hold on
plot(V_aa, fit_aa(V_aa), 'r-', 'LineWidth', 2)
plot(V_ab, I_ab*1e15, 'bo', 'LineWidth', 2)
plot(V_ab, fit_ab(V_ab), 'b-', 'LineWidth', 2)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.4,'GridLineStyle','--'); % ,'yscale','log'
xlabel('$V_{app}$ (V)','interpreter','latex')
ylabel('$I$ (fA)','interpreter','latex')
ylim([min(I_arr*1e15), 0])
%title(['$k_{red}$ for $' num2str(theta_list(1)) '^\circ$ TBG'],'FontSize',20,'interpreter', 'latex')
lh = legend('AA', '', 'AB','', 'Location','best','FontSize',20,'interpreter','latex');
legend boxoff
saveas(gcf,'/Users/mbabar/Desktop/PNP_plots/figsi2.png')
saveas(gcf,'/Users/mbabar/Desktop/PNP_plots/figsi2.svg');
hold off

%% Max current difference in AA/AB (contd.)

x = linspace(-0.6, 0.3, 101);
I_diff = fit_aa(x) - fit_ab(x);


[maxI, idx] = max(abs(I_diff));
limI = fit_aa(-0.7) % limiting current at most negative overpotential

maxI, x(idx)

figure()
plot(x, I_diff, 'LineWidth', 2)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1.4,'GridLineStyle','--'); % ,'yscale','log'
xlabel('$V_{app}$ (V)','interpreter','latex')
ylabel('$I$ (fA)','interpreter','latex')
saveas(gcf,'/Users/mbabar/Desktop/PNP_plots/fig4a_in.svg')

%% Domain areas

load('/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/prefactor/area_domains.mat')

alpha = 2.47;
theta = 1.1;
sc_alpha = alpha/(2*sind(theta/2)); % in Ang
sc_area = (sc_alpha^2)*sind(60)*1e-2; % in nm^2
moir_a = 60; %deg
bvec1 = [sc_alpha*cosd(moir_a/2), sc_alpha*sind(moir_a/2)].*0.1; % basis vector 1 in nm                                                   
bvec2 = [sc_alpha*cosd(moir_a/2), -sc_alpha*sind(moir_a/2)].*0.1; % basis vector 2 in nm
center_aa = [0,0];
center_ab = [2*bvec1(1)/3, 0];
center_sp = [bvec1(1), 0];

if theta<5
    AA_area = sc_area * 0.01*interp1(angle1, AA, theta);
elseif theta>5
    AA_area = sc_area * 0.01*31.59;
end

disp(append('Twist angle = ', string(theta)));
disp(append('Moire sc area: ', string(sc_area), ' nm^2, AA area: ', string(AA_area), ' nm^2'))
disp(append('AA radius: ', string((AA_area/pi)^(0.5)), ' nm'))
disp(append('Basis vector length: ', string(norm(bvec1)), ' nm'))
disp(append('Basis vector length /sqrt(3): ', string(norm(bvec1)/sqrt(3)), ' nm'))

%% Testing 

ft = fittype( 'c - (a/(1+exp(b*(x-d))))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [320 0.485375648722841 0.0 0.0];

[fit_aa, gof] = fit( V_aa, I_aa*1e15, ft, opts )

figure()
plot(V_aa, I_aa*1e15, 'ro', 'LineWidth', 2)
hold on 
plot(V_aa, fit_aa(V_aa), 'k-')
hold off

% Testing wrap

cell = [[1, 2]; [1, -2]]

A = [2; 2]
 
sol = linsolve(transpose(cell),A)

sol2 = rem(sol, 1.0)

sol3 = transpose(sol2)*cell