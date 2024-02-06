close all; clear; clc;

%% Parametri
k = 5e4;
beta = 0.1;
m = 5000;
b = 250;
n = 3;
x_e = [0.5;0];


%solo per visualizzione, pulsazione minima e massima
omega_plot_min = 1e-2;
omega_plot_max = 1e5;

%% Funzione di Trasferimento
A = [0,1;-k/m,-b/m];
B = [0;1/m];
C = [1,0];
D = 0;
I = eye(2);

s = tf('s');
GG=C*inv(s* I -A)*B+D;




%% Specifiche

% ampiezze gradini
WW = 2;
DD = 2;

% errore a regime
e_star = 0.01;

% Margine di fase
Mf_esp = 45;

% attenuazione disturbo sull'uscita
A_d = 60;
% non si puo` impostare omega_d_min esattamette a 0 perche` non 
% puo` essere rappresentato in scala logarimtica
omega_d_min = 0.00001;
omega_d_MAX = 0.5;


% attenuazione disturbo di misura
A_n = 50;
omega_n_min = 5*1e4;
omega_n_MAX = 5*1e6;


% Sovraelongazione massima e tempo d'assestamento all'1%
S_star = 15; %(S%<=15%)
T_star = 4*1e-3;


%% Diagramma di Bode

figure(1);
bode(GG,{omega_plot_min,omega_plot_max});
title("G(s)")
grid on, zoom on;


%% Diagrammi di Bode di G con specifiche
figure(2);
hold on;

% Calcolo specifiche S% => Margine di fase
xi_star = abs(log(S_star/100))/sqrt(pi^2 + log(S_star/100)^2);
Mf      = max(xi_star*100,Mf_esp);

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione critica)
omega_Ta_min = 1e-4; % lower bound per il plot
omega_Ta_max = 460/(Mf*T_star); % omega_c >= 460/(Mf*T^*) ~ 4.6
Bnd_Ta_x = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG,{omega_plot_min,omega_plot_max});
title("G(s) con specifiche")
grid on; zoom on;



% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_max;
omega_c_max = omega_n_min;

phi_up = Mf - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y = [phi_up; phi_up; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);



%% Regolatore statico

% valore minimo prescritto per L(0)
mu_s_error = (DD+WW)/e_star;
mu_s_dist  = 10^(A_d/20);

% guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(GG,0));
G_omega_d_MAX = abs(evalfr(GG,1j*omega_d_MAX));


RR_s = max(mu_s_error/G_0,mu_s_dist/G_omega_d_MAX);

% Sistema esteso
GG_e = RR_s*GG;


%% Diagrammi di Bode di Ge con specifiche
figure(3);
hold on;

% Calcolo specifiche S% => Margine di fase
xi_star = abs(log(S_star/100))/sqrt(pi^2 + log(S_star/100)^2);
Mf      = max(xi_star*100,Mf_esp);

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione critica)
omega_Ta_min = 1e-4; % lower bound per il plot
omega_Ta_max = 460/(Mf*T_star); % omega_c >= 460/(Mf*T^*) ~ 4.6
Bnd_Ta_x = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
title("G(s) con regolatore statico")
grid on; zoom on;



% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_max;
omega_c_max = omega_n_min;

phi_up = Mf - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y = [phi_up; phi_up; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);


%% Design del regolatore dinamico


Mf_star = Mf+5;
omega_c_star = 2400;

mag_omega_c_star_dB = abs(evalfr(GG_e,j*omega_c_star));
arg_omega_c_star    = rad2deg(angle(evalfr(GG_e,j*omega_c_star)));

M_star = 1/mag_omega_c_star_dB;
phi_star = Mf_star - 180 - arg_omega_c_star;

tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;

if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi');
    return;
end


RR_d = (1 + tau*s)/(1 + alpha*tau*s);
RR = RR_s*RR_d;

LL = RR*GG; % funzione di anello

figure(4);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
title("G(s) con rete anticipatrice")
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);


%% Aggiunta di poli

T_p = 1/4.81e4;
polo_reale = 1/(1+s * T_p);

RR_d = RR_d * polo_reale;

RR = RR_s*RR_d;

LL = RR*GG; % funzione di anello

figure(5)
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
title("G(s) con rete anticipatrice e polo in alta frequenza")
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);

% STOP qui per sistema con controllore dinamico + specifiche
if 0
    return;
end


%% Check prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);

% Risposta al gradino
figure(6);

T_simulation = 2*T_star;
[y_step,t_step] = step(WW*FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

LV = evalfr(WW*FF,0);

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure(7);

% Simulazione disturbo in uscita a pulsazione 0.125
omega_d = 0.125;
tt = 0:1e-2:2e2;

dd = 0;
for k = 1:1:4
    dd = dd + DD*sin(omega_d*k*tt) ;
end

y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('d(t)','y_d(t)')

%% Check disturbo di misura

% Funzione di sensitività complementare
FF = LL/(1+LL);
figure(8);

% Simulazione disturbo di misura a pulsazione 50000
omega_n = 5e4;

NN      = 0.2;

tt = 0:1e-5:2*1e-3;

nn = 0;

for k = 1:1:4
    nn = nn + NN*sin(omega_n*k*tt) ;
end


y_n = lsim(FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('n(t)','y_n(t)')


%% Animazione
%% parametri del disegno

interv = 0:0.1:30; % da 0 a 30 secondi con passo 0.1


% condizione iniziale della sospensione
pos_init = x_e(1); % [m]
vel_init = 0; % [m/s]
x0 = [pos_init; vel_init];

uu = zeros(length(interv),1);

% state-space model
modello = ss(A, B, C, D);

% evoluzione sistema
[YY, TT, XX] = lsim(modello, uu, interv, x0);

yy=YY;%variabile altezza
tt=TT;%periodo di simulazione



width_sos = 1;
hight_sos = 1;
% inizializza figura
figure(9);
clf;
min_T=0;
max_T=max(tt);
min_Y=min(yy);
max_Y=max(yy);
subplot(2,1,1)
axis([min_T max_T (min_Y-0.2) (max_Y+0.2)])
subplot(2,1,2)
axis([-10 10  0 5])

for i=1:5:length(tt)
    figure(9);
    clf;
    
    % aggiorna figura superiore
    subplot(2,1,1)
    plot(tt(1:i), yy(1:i),'b');
    axis([min_T max_T (min_Y-0.2) (max_Y+0.2)])
    title(['Evoluzione del sistema al tempo: ',num2str(tt(i)),' [s]']);
    legend('Uscita y(t)')
    grid on; box on;
    
    % aggiorna figura inferiore
    subplot(2,1,2)
    axis([-1.5 1.5  0 3])
    
    patch([hight_sos hight_sos -hight_sos -hight_sos],yy(i)+1.6+[-width_sos/2 width_sos/2 width_sos/2 -width_sos/2] ,'y') % cart
    hold on;
 
    plot([-10 -10 20],[0 0 0],'k','LineWidth',4) % ground.
    plot(-0.5+[0 0 0 .05 -.05 .05 -.05 .05 -.05 .05 -.05 0 0 0], ...
        10+[-10, -9.9, -9.9:((+9.9 +yy(i)-9)/9.9):yy(i)-9, yy(i)-9, yy(i)-8.9],'r','LineWidth',1) % spring
    title('Risposta del sistema')
    
    % pausa prima del prossimo frame
    pause(0.0001);
end