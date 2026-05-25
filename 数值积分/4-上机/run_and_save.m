%% run_and_save.m  --  运行四个题目并保存图片（内存安全版）
% 对 problem1 的 h 范围加限制避免内存溢出（n 最大 1e7）

clear; close all; clc;
% 非交互模式：batch 下必须关闭可见性才能正常保存图片
set(groot, 'DefaultFigureVisible', 'off');

outdir = fileparts(mfilename('fullpath'));

%% ============================================================
%  题目1：计算圆周率 π（修改版：h 最小到 1e-7 避免内存溢出）
%% ============================================================

f = @(x) 4 ./ (1 + x.^2);
a = 0; b = 1;
pi_exact = pi;

fprintf('=== 题目1：数值积分计算圆周率 π ===\n\n');

% 第一部分：误差随步长变化（限制 n <= 5e6 以节省内存）
h_vals = logspace(0, -7, 80);
n_vals = ceil((b - a) ./ h_vals);
n_vals = min(n_vals, 5e6);  % 限制最大节点数
err_trap = zeros(1, length(h_vals));
err_simp = zeros(1, length(h_vals));
actual_h = zeros(1, length(h_vals));

for k = 1:length(h_vals)
    n = n_vals(k);
    m = max(1, ceil(n/2));
    h_trap = (b - a) / n;
    h_simp = (b - a) / (2*m);
    actual_h(k) = h_trap;
    x_trap = linspace(a, b, n + 1);
    fx_trap = f(x_trap);
    T_n = h_trap / 2 * (fx_trap(1) + fx_trap(end) + 2 * sum(fx_trap(2:end-1)));
    err_trap(k) = abs(T_n - pi_exact);
    x_simp = linspace(a, b, 2*m + 1);
    fx_simp = f(x_simp);
    S_m = h_simp / 3 * (fx_simp(1) + fx_simp(end) + ...
          4 * sum(fx_simp(2:2:end-1)) + 2 * sum(fx_simp(3:2:end-2)));
    err_simp(k) = abs(S_m - pi_exact);
    clear x_trap fx_trap x_simp fx_simp;
end

[min_err_trap, idx_trap] = min(err_trap);
[min_err_simp, idx_simp] = min(err_simp);
fprintf('梯形最优 h=%.2e，误差=%.4e\n', actual_h(idx_trap), min_err_trap);
fprintf('Simpson最优 h=%.2e，误差=%.4e\n', actual_h(idx_simp), min_err_simp);

% Romberg 积分
max_level = 10; tol_romberg = 1e-12;
T_rom = zeros(max_level, max_level);
T_rom(1,1) = (b-a)/2*(f(a)+f(b));
romberg_result = T_rom(1,1); converged_level = 1;
fprintf('\nRomberg 表：\n');
fprintf('m=1: T=%.15f  误差=%.4e\n', T_rom(1,1), abs(T_rom(1,1)-pi_exact));
for m = 2:max_level
    h_prev = (b-a)/2^(m-2);
    h_cur  = (b-a)/2^(m-1);
    k_v = 1:2^(m-2);
    x_new = a + (k_v-0.5)*h_prev;
    T_rom(m,1) = 0.5*T_rom(m-1,1) + h_cur*sum(f(x_new));
    for j = 2:m
        c = 4^(j-1);
        T_rom(m,j) = (c*T_rom(m,j-1)-T_rom(m-1,j-1))/(c-1);
    end
    fprintf('m=%d: T=%.15f  误差=%.4e\n', m, T_rom(m,m), abs(T_rom(m,m)-pi_exact));
    if m>=2 && abs(T_rom(m,m)-T_rom(m-1,m-1))<tol_romberg
        romberg_result=T_rom(m,m); converged_level=m;
        fprintf('  收敛于第%d层\n',m); break;
    end
    romberg_result=T_rom(m,m); converged_level=m;
end
fprintf('\nRomberg结果: %.15f，误差: %.4e\n', romberg_result, abs(romberg_result-pi_exact));

% 自适应Simpson
tol_adaptive = 1e-10;
S_init = simp1_local(f,a,b);
[pi_adaptive, func_count] = adapt_simp_local(f,a,b,tol_adaptive,S_init);
fprintf('自适应Simpson: %.15f，误差: %.4e，求值次数: %d\n', ...
    pi_adaptive, abs(pi_adaptive-pi_exact), func_count);

%% 绘图1：误差随h变化
fig1 = figure('Position',[50 50 900 500]);
loglog(actual_h, err_trap, 'b-', 'LineWidth',1.5, 'DisplayName','复合梯形误差');
hold on;
loglog(actual_h, err_simp, 'r-', 'LineWidth',1.5, 'DisplayName','复合Simpson误差');
h_ref = logspace(0,-5,50);
loglog(h_ref, 0.3*h_ref.^2, 'b--', 'LineWidth',1.0, 'DisplayName','O(h²)');
loglog(h_ref, 0.05*h_ref.^4, 'r--', 'LineWidth',1.0, 'DisplayName','O(h⁴)');
xline(actual_h(idx_trap),'b:','梯形最优h','LineWidth',1.5,'FontSize',9);
xline(actual_h(idx_simp),'r:','Simpson最优h','LineWidth',1.5,'FontSize',9);
xlabel('步长 h（对数刻度）');
ylabel('|近似值 - π|（对数刻度）');
title('题目1：复合梯形 vs 复合Simpson：误差随步长变化');
legend('Location','NorthEast','FontSize',9);
grid on;
print(fig1, fullfile(outdir,'fig1_error_vs_h'), '-dpng', '-r150');
fprintf('图1已保存\n');

%% 绘图2：Romberg收敛 + 方法精度对比
fig2 = figure('Position',[100 100 900 400]);
subplot(1,2,1);
romberg_diag = zeros(1,converged_level);
for m=1:converged_level; romberg_diag(m)=abs(T_rom(m,m)-pi_exact); end
semilogy(1:converged_level, romberg_diag, 'bo-', 'LineWidth',1.5, 'MarkerSize',7);
xlabel('对分次数 m'); ylabel('|T_{m,m} - π|');
title('Romberg 对角线收敛过程'); grid on;

subplot(1,2,2);
h0=1e-4; n0=round((b-a)/h0); m0=ceil(n0/2);
x0t=linspace(a,b,n0+1); fx0t=f(x0t);
T0=h0/2*(fx0t(1)+fx0t(end)+2*sum(fx0t(2:end-1)));
x0s=linspace(a,b,2*m0+1); fx0s=f(x0s);
S0=(h0/2)/3*(fx0s(1)+fx0s(end)+4*sum(fx0s(2:2:end-1))+2*sum(fx0s(3:2:end-2)));
errs=[abs(T0-pi_exact),abs(S0-pi_exact),abs(romberg_result-pi_exact),abs(pi_adaptive-pi_exact)];
bar(1:4,errs);
set(gca,'YScale','log');
set(gca,'XTickLabel',{'梯形(h=1e-4)','Simpson(h=1e-4)','Romberg','自适应Simpson'});
xtickangle(20); ylabel('误差（对数刻度）');
title('各方法精度比较（题目1）'); grid on;
sgtitle('题目1：数值积分计算圆周率 π');
print(fig2, fullfile(outdir,'fig1_romberg_adaptive'), '-dpng', '-r150');
fprintf('图2已保存\n\n');
close all;

%% ============================================================
%  题目2：Planck 黑体辐射积分
%% ============================================================
fprintf('=== 题目2：Planck 黑体辐射积分 ===\n\n');

f_p = @(x) planck_local(x);
I_exact2 = 7*pi^4/120;
x_min = 1e-10; x_max = 50;
fprintf('精确值 7π⁴/120 = %.10f\n\n', I_exact2);

% 复合梯形
n_trap2 = 5000;
h_trap2 = (x_max-x_min)/n_trap2;
x_trap2 = linspace(x_min,x_max,n_trap2+1);
fx_trap2 = arrayfun(f_p,x_trap2);
I_trap2 = h_trap2/2*(fx_trap2(1)+fx_trap2(end)+2*sum(fx_trap2(2:end-1)));
fprintf('复合梯形(n=%d): I=%.10f，误差=%.4e\n',n_trap2,I_trap2,abs(I_trap2-I_exact2));

% 复合Simpson
m_simp2 = 2500;
h_simp2 = (x_max-x_min)/(2*m_simp2);
x_simp2 = linspace(x_min,x_max,2*m_simp2+1);
fx_simp2 = arrayfun(f_p,x_simp2);
I_simp2 = h_simp2/3*(fx_simp2(1)+fx_simp2(end)+...
          4*sum(fx_simp2(2:2:end-1))+2*sum(fx_simp2(3:2:end-2)));
fprintf('复合Simpson(m=%d): I=%.10f，误差=%.4e\n',m_simp2,I_simp2,abs(I_simp2-I_exact2));

% Romberg
T_r2 = zeros(12,12); max_lv2=12; tol_rom2=1e-8;
T_r2(1,1)=(x_max-x_min)/2*(f_p(x_min)+f_p(x_max));
I_romberg2=T_r2(1,1); rom_level2=1;
fprintf('\nRomberg对角线:\n');
fprintf('m=1: %.12f  误差=%.4e\n',T_r2(1,1),abs(T_r2(1,1)-I_exact2));
for m=2:max_lv2
    h_rp=(x_max-x_min)/2^(m-2); h_rc=(x_max-x_min)/2^(m-1);
    k_r=1:2^(m-2); xnr=x_min+(k_r-0.5)*h_rp;
    T_r2(m,1)=0.5*T_r2(m-1,1)+h_rc*sum(arrayfun(f_p,xnr));
    for j=2:m; c=4^(j-1); T_r2(m,j)=(c*T_r2(m,j-1)-T_r2(m-1,j-1))/(c-1); end
    fprintf('m=%d: %.12f  误差=%.4e\n',m,T_r2(m,m),abs(T_r2(m,m)-I_exact2));
    if abs(T_r2(m,m)-T_r2(m-1,m-1))<tol_rom2
        I_romberg2=T_r2(m,m); rom_level2=m;
        fprintf('  收敛于第%d层\n',m); break;
    end
    I_romberg2=T_r2(m,m); rom_level2=m;
end

% 自适应Simpson
S_init2=simp_one_local(f_p,x_min,x_max);
[I_adapt2,fcnt2]=adaptive_simp_local(f_p,x_min,x_max,1e-8,S_init2);
fprintf('\nRomberg: %.10f，误差=%.4e\n',I_romberg2,abs(I_romberg2-I_exact2));
fprintf('自适应Simpson: %.10f，误差=%.4e\n',I_adapt2,abs(I_adapt2-I_exact2));

%% 绘图3：Planck 积分
fig3 = figure('Position',[50 50 1000 600]);
subplot(2,2,1);
x_plot2=linspace(0.01,15,500); y_plot2=arrayfun(f_p,x_plot2);
fill([x_plot2,fliplr(x_plot2)],[y_plot2,zeros(1,500)],[0.7 0.85 1],'EdgeColor','b','LineWidth',1.2);
xlabel('x'); ylabel('f(x)=x³/(eˣ+1)');
title('Planck 被积函数（积分区域阴影）'); grid on;

subplot(2,2,2);
errs2=abs([I_trap2,I_simp2,I_romberg2,I_adapt2]-I_exact2);
bar(errs2); set(gca,'YScale','log');
set(gca,'XTickLabel',{'梯形','Simpson','Romberg','自适应'}); xtickangle(30);
ylabel('误差（对数刻度）'); title('四种方法精度比较'); grid on;

subplot(2,2,3);
rd2=zeros(1,rom_level2);
for m=1:rom_level2; rd2(m)=abs(T_r2(m,m)-I_exact2); end
semilogy(1:rom_level2,rd2,'ro-','LineWidth',1.5,'MarkerSize',7);
xlabel('对分层数 m'); ylabel('|T_{m,m}-I_{exact}|');
title('Romberg 积分收敛过程'); grid on;

subplot(2,2,4);
n_rng=round(logspace(1,4,25));
err_rng=zeros(1,length(n_rng));
for k=1:length(n_rng)
    nk=n_rng(k); hk=(x_max-x_min)/nk;
    xk=linspace(x_min,x_max,nk+1);
    fxk=arrayfun(f_p,xk);
    Tk=hk/2*(fxk(1)+fxk(end)+2*sum(fxk(2:end-1)));
    err_rng(k)=abs(Tk-I_exact2);
end
loglog(n_rng,err_rng,'b-o','MarkerSize',4,'LineWidth',1.2);
hold on;
loglog(n_rng,10./n_rng.^2,'r--','LineWidth',1.2,'DisplayName','O(n^{-2})');
xlabel('子区间数 n'); ylabel('误差'); title('梯形公式收敛阶验证');
legend({'实际误差','O(n^{-2})'},'Location','NorthEast'); grid on;
sgtitle('题目2：Planck 黑体辐射积分 ∫₀^∞ x³/(eˣ+1)dx = 7π⁴/120');
print(fig3, fullfile(outdir,'fig2_planck'), '-dpng', '-r150');
fprintf('\n图3已保存\n\n');
close all;

%% ============================================================
%  题目3：蒲丰投针
%% ============================================================
fprintf('=== 题目3：蒲丰投针估计 π ===\n\n');
rng(42);
l=1.0; d=2.0;
N_total3=1e6;
y_c=d*rand(1,N_total3); theta_c=pi*rand(1,N_total3);
hp=0.5*l*sin(theta_c);
cross_flag=(y_c<hp)|(y_c>d-hp);
N_cross=sum(cross_flag);
pi_est3=2*l*N_total3/(d*N_cross);
fprintf('N=%.0e，相交%d次，π估计=%.10f，误差=%.4e\n',N_total3,N_cross,pi_est3,abs(pi_est3-pi));

N_seq3=round(logspace(2,6,50));
pi_seq3=zeros(1,length(N_seq3)); err_seq3=zeros(1,length(N_seq3));
cumcross3=cumsum(cross_flag);
for k=1:length(N_seq3)
    Nk=N_seq3(k);
    ck=cumcross3(Nk);
    if ck>0
        pi_seq3(k)=2*l*Nk/(d*ck);
    else
        pi_seq3(k)=Inf;
    end
    err_seq3(k)=abs(pi_seq3(k)-pi);
end
p0=2*l/(pi*d);
theory_std3=pi*sqrt((1-p0)./(N_seq3*p0));

% 重复实验
M_rep3=300; N_each3=3000;
pi_rep3=zeros(1,M_rep3);
for m=1:M_rep3
    ym=d*rand(1,N_each3); tm=pi*rand(1,N_each3);
    hpm=(l/2)*sin(tm); cm=sum((ym<hpm)|(ym>d-hpm));
    if cm>0
        pi_rep3(m)=2*l*N_each3/(d*cm);
    else
        pi_rep3(m)=NaN;
    end
end
pi_rep3=pi_rep3(~isnan(pi_rep3));
fprintf('重复实验均值=%.8f，std=%.6f，理论std=%.6f\n',mean(pi_rep3),std(pi_rep3),pi*sqrt((1-p0)/(N_each3*p0)));

%% 绘图4：蒲丰投针
fig4=figure('Position',[50 50 1100 700]);
subplot(2,3,1);
n_show=200; y_sh=d*rand(1,n_show); th_sh=pi*rand(1,n_show);
hp_sh=(l/2)*sin(th_sh); cr_sh=(y_sh<hp_sh)|(y_sh>d-hp_sh);
hold on;
for lk=0:3; plot([0,l*3],[lk*d,lk*d],'k-','LineWidth',1.5); end
for k=1:n_show
    yy=mod(y_sh(k)+floor(rand()*4)*d,4*d);
    xc=l*rand(); dx=(l/2)*cos(th_sh(k)); dy=(l/2)*sin(th_sh(k));
    if cr_sh(k); plot([xc-dx,xc+dx],[yy-dy,yy+dy],'r-','LineWidth',1.2);
    else;        plot([xc-dx,xc+dx],[yy-dy,yy+dy],'b-','LineWidth',0.8); end
end
axis equal; xlim([0,l*3]); ylim([0,4*d]);
title(sprintf('投针示意（红=相交，蓝=不相交）\nn=%d',n_show)); xlabel('x'); ylabel('y');

subplot(2,3,2);
semilogx(N_seq3,pi_seq3,'b-','LineWidth',1.2); hold on;
yline(pi,'r--','LineWidth',1.5);
fill([N_seq3,fliplr(N_seq3)],[pi+2*theory_std3,fliplr(pi-2*theory_std3)],...
     [1 0.7 0.7],'FaceAlpha',0.3,'EdgeColor','none');
xlabel('投针次数 N'); ylabel('π 估计值');
title('π 估计值随实验次数的收敛'); grid on;

subplot(2,3,3);
loglog(N_seq3,err_seq3,'b-','LineWidth',1.2,'DisplayName','实际误差');
hold on;
loglog(N_seq3,theory_std3,'r--','LineWidth',1.5,'DisplayName','理论σ=C/√N');
xlabel('N'); ylabel('|π̂-π|'); title('收敛速度验证：误差∝1/√N');
legend('Location','NorthEast'); grid on;

subplot(2,3,4);
histogram(pi_rep3,25,'Normalization','pdf','FaceColor',[0.4 0.6 0.8],'EdgeColor','white');
hold on;
mu3=pi; sig3=pi*sqrt((1-p0)/(N_each3*p0));
xn3=linspace(pi-4*sig3,pi+4*sig3,200);
yn3=(1/(sig3*sqrt(2*pi)))*exp(-0.5*((xn3-mu3)/sig3).^2);
plot(xn3,yn3,'r-','LineWidth',2); xline(pi,'k--','LineWidth',1.5);
xlabel('π 估计值'); ylabel('概率密度');
title(sprintf('π 估计值分布\n均值=%.4f，std=%.4f',mean(pi_rep3),std(pi_rep3)));

subplot(2,3,5);
th_ax=linspace(0,pi,200); P_th=(l/d)*sin(th_ax);
fill([th_ax,fliplr(th_ax)],[P_th,zeros(1,200)],[0.8 0.9 0.7],'EdgeColor','g','LineWidth',1.2);
xlabel('θ（弧度）'); ylabel('P(相交|θ)');
title(sprintf('条件相交概率\n面积=2l/(πd)=%.4f',2*l/(pi*d)));
xticks([0 pi/4 pi/2 3*pi/4 pi]); xticklabels({'0','π/4','π/2','3π/4','π'}); grid on;

subplot(2,3,6);
ns6=1:min(10000,N_total3);
cum_pi6=2*l*ns6./(d*max(1,cumcross3(ns6)));
plot(ns6,cum_pi6,'b-','LineWidth',0.8); hold on;
yline(pi,'r--','LineWidth',1.5);
xlabel('投针次数 N'); ylabel('π 累积估计值');
title('π 估计值实时收敛（前10000次）'); grid on;
sgtitle('题目3：蒲丰投针 Monte Carlo 模拟估计 π');
print(fig4, fullfile(outdir,'fig3_buffon'), '-dpng', '-r150');
fprintf('图4已保存\n\n');
close all;

%% ============================================================
%  题目4：Monte Carlo 计算 ∫₀¹ e^{-x} dx
%% ============================================================
fprintf('=== 题目4：Monte Carlo 计算 I=∫₀¹e^{-x}dx ===\n\n');
rng(2024);
f4=@(x)exp(-x);
I_exact4=1-exp(-1);
N4=5e5;
fprintf('精确值: %.15f\n\n',I_exact4);

x_hit=rand(1,N4); y_hit=rand(1,N4);
hits=sum(y_hit<=f4(x_hit));
I_hm=hits/N4;
fprintf('随机投点法: I=%.10f，误差=%.4e\n',I_hm,abs(I_hm-I_exact4));

x_smp=rand(1,N4);
I_sm=mean(f4(x_smp));
fprintf('样本均值法: I=%.10f，误差=%.4e\n\n',I_sm,abs(I_sm-I_exact4));

% 收敛分析
N_list4=round(logspace(2,5,40));
err_hit4=zeros(1,length(N_list4)); err_smp4=zeros(1,length(N_list4));
hit_all4=(y_hit<=f4(x_hit)); fval_all4=f4(x_smp);
cumhit4=cumsum(hit_all4); cumfv4=cumsum(fval_all4);
for k=1:length(N_list4)
    Nk=N_list4(k);
    err_hit4(k)=abs(cumhit4(Nk)/Nk-I_exact4);
    err_smp4(k)=abs(cumfv4(Nk)/Nk-I_exact4);
end
std_hit4=sqrt(I_exact4*(1-I_exact4)./N_list4);
var_f4=(1-exp(-2))/2-I_exact4^2;
std_smp4=sqrt(var_f4./N_list4);

%% 绘图5：Monte Carlo 积分
fig5=figure('Position',[50 50 1100 800]);
subplot(2,3,1);
xc4=linspace(0,1,300); yc4=f4(xc4);
fill([xc4,fliplr(xc4)],[yc4,zeros(1,300)],[0.65 0.85 0.95],'EdgeColor','b','LineWidth',2);
hold on; plot(xc4,yc4,'b-','LineWidth',2);
text(0.35,0.3,sprintf('面积=%.6f',I_exact4),'FontSize',12,'Color','b','FontWeight','bold');
xlabel('x'); ylabel('y=e^{-x}'); title('积分区域（曲线下阴影）');
xlim([0,1]); ylim([0,1.1]); grid on;

subplot(2,3,2);
Ns2=2000;
fill([xc4,fliplr(xc4)],[yc4,zeros(1,300)],[0.9 0.95 0.98],'EdgeColor','none');
hold on;
xsh=x_hit(1:Ns2); ysh=y_hit(1:Ns2); hsh=hit_all4(1:Ns2);
plot(xsh(~hsh),ysh(~hsh),'r.','MarkerSize',2,'DisplayName','未命中');
plot(xsh(hsh), ysh(hsh), 'b.','MarkerSize',2,'DisplayName','命中');
plot(xc4,yc4,'k-','LineWidth',1.5);
xlabel('x'); ylabel('y'); title(sprintf('随机投点法示意(N=%d)',Ns2));
xlim([0,1]); ylim([0,1]);

subplot(2,3,3);
histogram(f4(x_smp(1:5000)),25,'Normalization','pdf','FaceColor',[0.4 0.7 0.4],'EdgeColor','white');
hold on;
tr4=linspace(exp(-1),1,200); pdf4=1./tr4;
plot(tr4,pdf4,'r-','LineWidth',2);
xline(I_exact4,'b--','LineWidth',1.5);
xlabel('f(x)=e^{-x}的取值'); ylabel('概率密度');
title('样本均值法：f(X)的分布');

subplot(2,3,4);
loglog(N_list4,err_hit4,'r-','LineWidth',1.5,'DisplayName','投点法误差');
hold on;
loglog(N_list4,err_smp4,'b-','LineWidth',1.5,'DisplayName','样本均值法误差');
loglog(N_list4,std_hit4,'r--','LineWidth',1.0,'DisplayName','投点法理论σ');
loglog(N_list4,std_smp4,'b--','LineWidth',1.0,'DisplayName','样本均值法理论σ');
Nrf=[1e2,1e5]; loglog(Nrf,0.5./sqrt(Nrf),'k:','LineWidth',1.2,'DisplayName','O(N^{-1/2})');
xlabel('样本数 N'); ylabel('误差'); title('两种方法收敛速度比较');
legend('Location','SouthWest','FontSize',7); grid on;

subplot(2,3,5);
Nlv=min(N4,5e4); na5=1:Nlv;
lh5=cumhit4(1:Nlv)./na5; ls5=cumfv4(1:Nlv)./na5;
plot(na5,lh5,'r-','LineWidth',0.8,'DisplayName','投点法');
hold on; plot(na5,ls5,'b-','LineWidth',0.8,'DisplayName','样本均值法');
yline(I_exact4,'k--','LineWidth',1.5,DisplayName=sprintf('精确值=%.6f',I_exact4));
xlabel('样本数 N'); ylabel('I 累积估计值');
title('实时收敛（前5万次）');
legend('Location','NorthEast','FontSize',8); grid on;
ylim([I_exact4-0.05, I_exact4+0.05]);

subplot(2,3,6);
Mr6=500; Ne6=500;
ph6=zeros(1,Mr6); ps6=zeros(1,Mr6);
for m=1:Mr6
    xh6=rand(1,Ne6); yh6=rand(1,Ne6);
    ph6(m)=sum(yh6<=f4(xh6))/Ne6;
    xs6=rand(1,Ne6); ps6(m)=mean(f4(xs6));
end
histogram(ph6,25,'Normalization','pdf','FaceColor',[1 0.6 0.6],'EdgeColor','white','FaceAlpha',0.6,...
    'DisplayName',sprintf('投点法(σ=%.4f)',std(ph6)));
hold on;
histogram(ps6,25,'Normalization','pdf','FaceColor',[0.6 0.7 1],'EdgeColor','white','FaceAlpha',0.6,...
    'DisplayName',sprintf('样本均值法(σ=%.4f)',std(ps6)));
xline(I_exact4,'k--','LineWidth',2);
xlabel('I 估计值'); ylabel('概率密度');
title(sprintf('估计值分布对比(M=%d,N=%d)',Mr6,Ne6));
legend('Location','NorthWest','FontSize',7);
sgtitle('题目4：Monte Carlo 计算 I=∫₀¹e^{-x}dx');
print(fig5, fullfile(outdir,'fig4_monte_carlo'), '-dpng', '-r150');
fprintf('图5已保存\n\n');
close all;

fprintf('=== 所有图片已保存至 %s ===\n', outdir);

%% ---- 局部辅助函数 ----
function S=simp1_local(f,a,b); m=(a+b)/2; S=(b-a)/6*(f(a)+4*f(m)+f(b)); end
function [I,cnt]=adapt_simp_local(f,a,b,tol,S_ab)
    m=(a+b)/2; Sam=simp1_local(f,a,m); Smb=simp1_local(f,m,b);
    S2=Sam+Smb; cnt=3;
    if abs(S2-S_ab)/15<=tol; I=S2+(S2-S_ab)/15;
    else
        [Il,cl]=adapt_simp_local(f,a,m,tol/2,Sam);
        [Ir,cr]=adapt_simp_local(f,m,b,tol/2,Smb);
        I=Il+Ir; cnt=cnt+cl+cr;
    end
end
function y=planck_local(x)
    if x<=0; y=0; elseif x>500; y=0; elseif x>37; y=x^3*exp(-x);
    else; y=x^3/(exp(x)+1); end
end
function S=simp_one_local(f,a,b); m=(a+b)/2; S=(b-a)/6*(f(a)+4*f(m)+f(b)); end
function [I,cnt]=adaptive_simp_local(f,a,b,tol,S_ab)
    m=(a+b)/2; Sam=simp_one_local(f,a,m); Smb=simp_one_local(f,m,b);
    S2=Sam+Smb; cnt=3;
    if abs(S2-S_ab)/15<=tol; I=S2+(S2-S_ab)/15;
    else
        [Il,cl]=adaptive_simp_local(f,a,m,tol/2,Sam);
        [Ir,cr]=adaptive_simp_local(f,m,b,tol/2,Smb);
        I=Il+Ir; cnt=cnt+cl+cr;
    end
end
