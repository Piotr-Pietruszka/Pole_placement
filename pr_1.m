
% Control system model
Ts = 2;
Gp_c = tf(1, [15, 1], 'InputDelay', 2)
Gp_d = c2d(Gp_c, Ts, 'zoh')


% Reference model
ks=12; Tn=12;
Gref = tf(ks, [Tn/2.2, 1])
Gref_d = c2d(Gref, Ts, 'zoh')

% plot zeros and poles
% plot__pz(Gp_c, Gp_d, Gref, Gref_d)


% Polynomial's coefficients
k = 2;
[B, A] = tfdata(Gp_d, 'v');
[Bm, Am] = tfdata(Gref_d, 'v');

% special case
a1 = A(2);
am1 = Am(2);
b0 = B(2);
bm0 = Bm(2);
f1 = am1 - a1;
g0 = -a1*f1;


S = g0; % 0.1248
T = bm0; % 3.684
r1 = b0*f1; % 0.0227
r0 = b0; % 0.1248


% sim_d(S, T, r1, r0, Ts, k, a1, b0, Tn);
% step(Gref_d)

% figure(1)
% hold on
sim_c(S, T, r1, r0, Ts, k, Gp_c, Tn)


function [] = sim_d(S, T, r1, r0, Ts, k, a1, b0, Tn)
    % Discrete-time simulation
    N = 500;
    t_sim = 1:Ts:Ts*N;    
    u = zeros(1, N+1);
    y = zeros(1, N+1);
    yr = gen_square_wave(20*Tn/Ts, N+1);
    for t = 2:1:(N-k)
        u(t) = (-r1*u(t-1) + T*yr(t) - S*y(t)) / r0;
        y(t+2) = -a1*y(t+1) + b0*u(t);
    end
    
    % remove element 0
    u = u(2:end);
    y = y(2:end);
    yr = yr(2:end);
    
    % plot
    figure(1)
    stairs(t_sim,u)
    title("u(t)")
    xlabel("t [s]")
    
    figure(2)
    hold on
    stairs(t_sim, y)
    title("y(t)")
    stairs(t_sim, yr)
    legend("y(t)", "yr(t)")
    xlabel("t [s]")

end


function sim_c(S, T, r1, r0, Ts, k, Gp_c, Tn)
    
    N = 500;
    total_time = N*Ts;
    step = 0.01; % simulation step
    M = ceil(total_time/step);
    Ts_steps = round(Ts/step) % no time-steps in one sample period
    
    % Discrete state-space model of Gp
    gd = c2d(Gp_c, step, 'zoh')
    [a, b, c, d] = ssdata(gd);
    delay = gd.InputDelay;
    % Signals
    u = zeros(1, M);
    y = zeros(1, M);
    x = zeros(1, M);
    t_sim = (1:1:M-1)*step;
    yr = gen_square_wave(round(20*Tn/step), M+1);
    % temp? - discrete simulation ---------------
    ud = zeros(1, M);
    yd = zeros(1, M);
    a1 = -0.8752;
    b0 = 0.1248;
    counter = 0;
    %     ------------------------
    % Time loop
    for t = delay+1:1:M-1
        % change control signal every Ts
        if mod(t, Ts_steps)== 0
           u(t:t+Ts_steps) = (-r1*u(t-Ts_steps) + T*yr(t) - S*y(t)) / r0;
           %     ------------------------
           ud(t:t+Ts_steps) = (-r1*ud(t-Ts_steps) + T*yr(t) - S*yd(t)) / r0;
           yd(t+Ts_steps:t+2*Ts_steps) = -a1*yd(t+Ts_steps) + b0*ud(t);
           if counter ==1 % temp2
               % print something
           end
           counter = counter + 1;
           %     ------------------------
        end
        % calculate new system response very step
        x(t+1) = a*x(t) + b*u(t-delay);
        y(t) = c*x(t);
    end
    %     ------------------------
    ud = ud(1:end-1);
    yd = yd(1:end-Ts_steps-1);
%     size(yd)
%     size(t_sim)
%     plot(t_sim, ud)
    plot(t_sim, yd)
    hold on
    %     ------------------------
    
    y = y(1:end-1);
    u = u(1:end-1);
%     plot(t_sim, u)
    plot(t_sim, y)
%     plot(yr)
    
end


% genereate square wave
% period - period in samples
% M - length
function [sq] = gen_square_wave(period, M)
    sq = zeros(1, M);
    for i = 1:2*period:(M-period)
        sq(i:i+period) = 1;
    end
end


% Plot and print zeros and poles of given systems
function [] = plot__pz(Gp_c, Gp_d, Gref_c, Gref_d)
    figure(1)
    
    % continous Gp
    zeros_c = zero(Gp_c)
    poles_c = pole(Gp_c)
    subplot(2, 2, 1)
    pzmap(Gp_c)
    title("Poles and zeros of cont Gp") 
    
    % discrete Gd
    zeros_d = zero(Gp_d)
    poles_d = pole(Gp_d)
    subplot(2, 2, 2)
    pzmap(Gp_d)
    hold on
    plot(0, 0, 'o') % z^-1 delay isn't normally  plotted 
    title("Poles and zeros of discrete Gp") 
    
    % continous Gref
    zeros_ref_c = zero(Gref_c)
    poles_ref_c = pole(Gref_c)
    subplot(2, 2, 3)
    pzmap(Gref_c)
    title("Poles and zeros of cont Gref") 
    
    % discrete Gref
    zeros_ref_d = zero(Gref_d)
    poles_ref_d = pole(Gref_d)
    subplot(2, 2, 4)
    pzmap(Gref_d)
    title("Poles and zeros of discrete Gref") 

end
