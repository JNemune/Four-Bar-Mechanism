% init
clear;
clc;

% data
m1 = 1.2;
m2 = 7;
m3 = 1.8;
g = 9.80665;
l1 = 80e-3;
l2 = 240e-3;
l3 = 200e-3;
cx = -190e-3;
cy = -70e-3;
theta_ = 2 * pi;

% I
i1o = 1/3 * m1 * l1 ^ 2;
i2g = 1/12 * m2 * l2 ^ 2;
i3c = 1/3 * m3 * l3 ^ 2;

% alpha beta and thete derivations
syms theta alpha beta alpha_ beta_ alpha__ beta__;
theta_vec = [0 0 theta_];
alpha_vec = [0 0 alpha_];
beta_vec = [0 0 beta_];
alpha__vec = [0 0 alpha__];
beta__vec = [0 0 beta__];

% solve for alpha and beta
roa = l1 .* [cos(theta) sin(theta) 0];
rab = l2 .* [cos(alpha) sin(alpha) 0];
rbc = l3 .* [cos(beta) sin(beta) 0];
roc = [cx cy 0];

slv1 = solve(roa + rab + rbc == roc, [alpha beta]);
alpha = simplify(slv1.alpha(2));
beta = simplify(slv1.beta(2));
rab = l2 .* [cos(alpha) sin(alpha) 0];
rbc = l3 .* [cos(beta) sin(beta) 0];

% solve for omega
va = cross(theta_vec, roa);
vb = va + cross(alpha_vec, rab);

slv2 = solve([0 0 0] == vb + cross(beta_vec, rbc), [alpha_, beta_]);
alpha_ = simplify(slv2.alpha_);
beta_ = simplify(slv2.beta_);
alpha_vec = [0 0 alpha_];
beta_vec = [0 0 beta_];
vb = va + cross(alpha_vec, rab);

% solve for alpha
aa = cross(theta_vec, cross(theta_vec, roa));
ab = aa + cross(alpha_vec, cross(alpha_vec, rab)) + cross(alpha__vec, rab);

slv3 = solve([0 0 0] == ab + cross(beta_vec, cross(beta_vec, rbc)) + cross(beta__vec, rbc), [alpha__, beta__]);
alpha__ = simplify(slv3.alpha__);
beta__ = simplify(slv3.beta__);
alpha__vec = [0 0 alpha__];
beta__vec = [0 0 beta__];
ab = aa + cross(alpha_vec, cross(alpha_vec, rab)) + cross(alpha__vec, rab);

% acceleration of G1, G2, G3
ag1 = cross(theta_vec, cross(theta_vec, roa ./ 2));
ag2 = aa + cross(alpha_vec, cross(alpha_vec, rab ./ 2)) + cross(alpha__vec, rab ./ 2);
ag3 = ab + cross(beta_vec, cross(beta_vec, rbc ./2)) + cross(beta__vec, rbc ./ 2);

syms  fox foy f21x f21y f32x f32y fcx fcy M;

fo = [fox foy 0];
f21 = [f21x f21y 0];
f12 = -f21;
f32 = [f32x f32y 0];
f23 = -f32;
fc = [fcx fcy 0];

w1 = [0 -m1 * g 0];
w2 = [0 -m2 * g 0];
w3 = [0 -m3 * g 0];
mo = [0 0 M];

eqns = [
    fo + w1 + f21 == m1 .* ag1;
    mo + cross(roa ./ 2, w1) + cross(roa, f21) == [0 0 0];
    
    f12 + w2 + f32 == m2 .* ag2;
    cross(-rab ./ 2, f12) + cross(rab ./ 2, f32) == i2g .* alpha__vec;
    
    f23 + w3 + fc == m3 .* ag3;
    cross(rbc ./ 2, w3) + cross(rbc, fc) == i3c .* beta__vec;
];

slv4 = solve(eqns, [fox foy f21x f21y f32x f32y fcx fcy M]);

M_sim = simplify(slv4.M);
fo_sim = sqrt(simplify(slv4.fox)^2 + simplify(slv4.foy)^2);
fa_sim = sqrt(simplify(slv4.f21x)^2 + simplify(slv4.f21y)^2);
fb_sim = sqrt(simplify(slv4.f32x)^2 + simplify(slv4.f32y)^2);
fc_sim = sqrt(simplify(slv4.fcx)^2 + simplify(slv4.fcy)^2);

% output and plot
theta_values = linspace(0, 2 * pi, 100);
M_values = double(subs(M_sim, theta, theta_values));

figure;
plot(rad2deg(theta_values), M_values);
xlabel('Theta');
ylabel('M');
title('M vs. Theta');
grid on;

fo = double(subs(fo_sim, theta, theta_values));
fa = double(subs(fa_sim, theta, theta_values));
fb = double(subs(fb_sim, theta, theta_values));
fc = double(subs(fc_sim, theta, theta_values));

figure;
plot(rad2deg(theta_values), fo);
hold on;
plot(rad2deg(theta_values), fa);
plot(rad2deg(theta_values), fb);
plot(rad2deg(theta_values), fc);
hold off;
legend('O', 'A', 'B', 'C');
xlabel('Theta');
ylabel('Forces Magnitude');
title('Forces Magnitude vs. Theta');
grid on;

theta_m_max = vpasolve(diff(M_sim) == 0, theta, 2.5);
theta_o_max = vpasolve(diff(fo_sim) == 0, theta, 3);
theta_a_max = vpasolve(diff(fa_sim) == 0, theta, 3);
theta_b_max = vpasolve(diff(fb_sim) == 0, theta, 3);
theta_c_max = vpasolve(diff(fc_sim) == 0, theta, 3);

subs_theta_m_max = subs(M_sim, theta, theta_m_max);
subs_theta_o_max = subs(fo_sim, theta, theta_o_max);
subs_theta_a_max = subs(fa_sim, theta, theta_a_max);
subs_theta_b_max = subs(fb_sim, theta, theta_b_max);
subs_theta_c_max = subs(fc_sim, theta, theta_c_max);

fprintf('Theta_m_max: %.2f\n', rad2deg(theta_m_max));
fprintf('M(theta_m_max): %.2f\n', double(subs_theta_m_max));

fprintf('Theta_o_max: %.2f\n', rad2deg(theta_o_max));
fprintf('Fo(theta_o_max): %.2f\n', double(subs_theta_o_max));

fprintf('Theta_a_max: %.2f\n', rad2deg(theta_a_max));
fprintf('Fa(theta_a_max): %.2f\n', double(subs_theta_a_max));

fprintf('Theta_b_max: %.2f\n', rad2deg(theta_b_max));
fprintf('Fb(theta_b_max): %.2f\n', double(subs_theta_b_max));

fprintf('Theta_c_max: %.2f\n', rad2deg(theta_c_max));
fprintf('Fc(theta_c_max): %.2f\n', double(subs_theta_c_max));
