

%clc;
%clear all;

delete(gcp('nocreate'));
%poolobj = parpool ;
%parpool local 2
%parpool(12);
%spmd

%Pauli_matrix
i = sqrt(-1);   Sigma_x = [0  1; 1  0]/2;   Sigma_y = [0 -i; i  0]/2;   Sigma_z = [1  0; 0 -1]/2;   
Identity = eye(2);

%Control_Parameters
B_0 = 1; % 1T
gamma_e =28000;% 28000;%28024.95164;%B*10^6(Hz/T);
W_0 = B_0*gamma_e ; 

A_1 = 80;%1440;(MHZ)
A_2 = 120;%2420;
%J = interaction(t,t_f) ; %10;%(80)

omega_x1_arr = [-0.000620588238213570,0.000443967356049757,0.000884029395718997,0.000222778367359921,-0.000812400612925980,-0.000313002137147762,-0.000747436871576532,0.000902960093355495];
omega_y1_arr = [-0.000998515957636616,-0.000236882078744332,0.000107814066130186,0.000533518618328375,-2.54633067364100e-06,0.000887237930345645,0.000547025270128294,-7.71545698550659e-05];

Parameters = [omega_x1_arr  omega_y1_arr ];% omega_x2_arr  omega_y2_arr ];

t_i = 0;
t_f = 0.05 ;%0.6;
U_i = eye(4);   U_i_r = reshape(U_i , [] , 1);
U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);    
if isempty(odeget(opts,'Stats'))
    odeset(opts,'Stats','on');
end
    
[Tsol,Usol] = ode45(@(t,U) Schrodinger_H_p_rf(t,U,Parameters), [t_i  t_f] , U_i_r , opts);
    
Uf = Usol(end,:);
Uf = reshape(Uf,[],4)
    
    %Uf = reshape(Uf,[],round(sqrt(length(Uf))));
    %    finalState = U_0'*Uf;
    %Uf = Uf'*Uf;
    %U_t = U_t'*U_t;
    %F_ele = trace(sqrtm(sqrtm(U_t)*Uf*sqrtm(U_t)));
    %infidelity = 1-sqrt(conj(F_ele) *F_ele)/16;
    
U_t'
trace(U_t'*Uf)
infidelity = 1-(abs(trace((U_t'*Uf)))^2/16)
 
 


function dUdt = Schrodinger_H_p_rf(t,U,Parameters)
    i = sqrt(-1);
    B_0 = 1;    gamma_e = 28000;
    t_f = 0.05 ; 
    
    A_1 = 80;    A_2 = 120;
    delta_A = (A_2 - A_1)/2;    ave_A = (A_2 + A_1)/2;
    J = interaction(t,t_f) ;%10;
    
    omega_x1_arr = Parameters(1:length(Parameters)/2);
    omega_y1_arr = Parameters(length(Parameters)/2+1:length(Parameters));
    %omega_x2_arr = Parameters(13:18);
    %omega_y2_arr = Parameters(19:24);
    B_x_s_1 = Control_MF_x_1(omega_x1_arr, t, t_f,B_0, gamma_e, delta_A);
    B_y_s_1 = Control_MF_y_1(omega_y1_arr, t, t_f,B_0, gamma_e, delta_A);
    B_x_s_2 = Control_MF_x_2(omega_x1_arr, t, t_f,B_0, gamma_e, delta_A);
    B_y_s_2 = Control_MF_y_2(omega_y1_arr, t, t_f,B_0, gamma_e, delta_A);
    
    H = zeros(4);
    
    H(1,1) = gamma_e*B_0 + delta_A/2 + J/4;     H(4,4) = -gamma_e*B_0 - delta_A/2 + J/4;    
    H(2,2) = ave_A/2 - J/4;     H(3,3) = -ave_A/2 - J/4;    
    H(2,3) =  J/2;    H(3,2) =  J/2;
    
    H(1,2) = 1/2*gamma_e*(B_x_s_2 - (B_y_s_2*i));%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x_s_1 - (B_y_s_1*i));%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x_s_2 + (B_y_s_2*i));%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x_s_1 - (B_y_s_1*i));%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x_s_1 + (B_y_s_1*i));%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x_s_2 - (B_y_s_2*i));%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x_s_1 + (B_y_s_1*i));%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x_s_2 + (B_y_s_2*i));%*RF_com_n;
    
    H_p_rf = kron(eye(4),H);
    
    %U_m = reshape(U,[],4);
    dUdt = -i * 2*pi*  H_p_rf * U ;
    %dUdt_m = -i * 2*pi*  H_p_rf * U_m ;
    
    %dUdt   = reshape(dUdt_m,[],1);
 
end

function J = interaction(t,t_f) 
    J = 60*sin(pi*t/t_f) ;    
end

function B_x_s_1 = Control_MF_x_1 (omega_x1_arr, t, t_f, B_0, gamma_e, delta_A)
    i = sqrt(-1);  
    B_x_s_dr_1 = 0 ;
        for n = 1:length(omega_x1_arr)
            B_x_dr_1 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f) ;
            B_x_s_dr_1 = B_x_s_dr_1 + B_x_dr_1;
        end
    B_x_s_1 = B_x_s_dr_1* (exp((i*2*pi*(B_0*gamma_e + delta_A/2))*t) + exp((-i*2*pi*(B_0*gamma_e + delta_A/2))*t))/2; 
end


function B_y_s_1 = Control_MF_y_1 (omega_y1_arr, t, t_f, B_0, gamma_e, delta_A)
    i = sqrt(-1);  
    B_y_s_dr_1 = 0 ;
        for n = 1:length(omega_y1_arr)
            B_y_dr_1 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_1 = B_y_s_dr_1 + B_y_dr_1;
        end
    B_y_s_1 = B_y_s_dr_1*(exp((i*2*pi*(B_0*gamma_e + delta_A/2))*t) + exp((-i*2*pi*(B_0*gamma_e + delta_A/2))*t))/2 ; 
end

function B_x_s_2 = Control_MF_x_2 (omega_x1_arr, t, t_f, B_0, gamma_e, delta_A)   
    i = sqrt(-1);
    B_x_s_dr_2 = 0 ;
        for n = 1:length(omega_x1_arr)
            B_x_dr_2 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f);
            B_x_s_dr_2 = B_x_s_dr_2 + B_x_dr_2;
        end
    B_x_s_2 = B_x_s_dr_2*(exp((i*2*pi*(B_0*gamma_e + delta_A/2))*t) + exp((-i*2*pi*(B_0*gamma_e + delta_A/2))*t))/2; 
end

function B_y_s_2 = Control_MF_y_2 (omega_y1_arr, t, t_f, B_0, gamma_e, delta_A)
    i = sqrt(-1);
    B_y_s_dr_2 = 0 ;
        for n = 1:length(omega_y1_arr)
            B_y_dr_2 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_2 = B_y_s_dr_2 + B_y_dr_2;
        end
    B_y_s_2 = B_y_s_dr_2*(exp((i*2*pi*(B_0*gamma_e + delta_A/2))*t) + exp((-i*2*pi*(B_0*gamma_e + delta_A/2))*t))/2;
    %cos((i*2*pi*(B_0*gamma_e + delta_A/2))*t)   ; 
end


%{
clc;
clear all;
i = sqrt(-1);
Sigma_x = [0  1; 1  0]/2;   Sigma_y = [0 -i; i  0]/2;   Sigma_z = [1  0; 0 -1]/2;   Identity = eye(2);

omega_x1_arr = [-0.000620588238213570,0.000443967356049757,0.000884029395718997,0.000222778367359921,-0.000812400612925980,-0.000313002137147762,-0.000747436871576532,0.000902960093355495];
omega_y1_arr = [-0.000998515957636616,-0.000236882078744332,0.000107814066130186,0.000533518618328375,-2.54633067364100e-06,0.000887237930345645,0.000547025270128294,-7.71545698550659e-05];
Parameters = [omega_x1_arr omega_y1_arr ];%omega_x2_arr omega_y2_arr];
%Control_Parameters
B_0 = 1; % 1T
gamma_e =28000;% 28000;%28024.95164;%B*10^6(Hz/T);
W_0 = B_0*gamma_e ; 

A_1 = 80;%1440;(MHZ)
A_2 = 120;%2420;
%Parameter_for_schrodinger
t_i = 0;    t_f = 0.05 ;
U_i = eye(4);   U_i_r = reshape(U_i , [] , 1);
U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
   
    opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end
    
    [Tsol,Usol] = ode45(@(t,U) Schrodinger_H_p_rf(t,U,Parameters), [t_i  t_f] , U_i_r , opts);
    
    Uf = Usol(end,:);
    Uf = reshape(Uf,[],4);
    %[V,D,W] = eig(Uf)
    
    %Uf = reshape(Uf,[],round(sqrt(length(Uf))));
    %    finalState = U_0'*Uf;
    
    %F_ele = trace(sqrtm(sqrtm(U_t)*Uf*sqrtm(U_t)));
    %infidelity = 1-(F_ele *F_ele')/16
    
    infidelity_J1 = 1-(abs(trace(Uf*U_t')))^2/16
    %infidelity = 1-(abs(trace(Uf*obj.U_t')))^2/16;
    
function dUdt = Schrodinger_H_p_rf(t,U,Parameters)
    i = sqrt(-1);
    B_0 = 1;    W_0 = 28000;    gamma_e = 28000;
    t_f = 0.05;
    A_1 = 80;%1440;
    A_2 = 120;%2420;
    J = interaction(t,t_f);
    
    %RF_com = exp(i*2*pi*(W_0+(A_1-A_2)/4)*t);
    omega_x1_arr = Parameters(1:length(Parameters)/2);
    omega_y1_arr = Parameters(length(Parameters)/2+1:length(Parameters));
    %omega_x2_arr = Parameters(11:15);
    %omega_y2_arr = Parameters(16:20);
    B_x_s_1 = Control_MF_x_1(B_0,gamma_e,A_1, A_2, omega_x1_arr, t, t_f);
    B_y_s_1 = Control_MF_y_1(B_0,gamma_e,A_1, A_2, omega_y1_arr, t, t_f);
    B_x_s_2 = Control_MF_x_2(B_0,gamma_e,A_1, A_2, omega_x1_arr, t, t_f);
    B_y_s_2 = Control_MF_y_2(B_0,gamma_e,A_1, A_2, omega_y1_arr, t, t_f);
    
    H = zeros(4);
    
    H(1,1) = gamma_e*B_0 + (A_2 - A_1)/4 + J/4;
    H(4,4) = -gamma_e*B_0 - (A_2 - A_1)/4 + J/4;
    
    H(2,2) = (A_2 + A_1)/4 - J/4;
    H(3,3) = -(A_2 + A_1)/4 - J/4;
    
    H(2,3) =  J/2;
    H(3,2) =  J/2;
    
    H(1,2) = 1/2*gamma_e*(B_x_s_2 - (B_y_s_2*i));%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x_s_1 - (B_y_s_1*i));%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x_s_2 + (B_y_s_2*i));%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x_s_1 - (B_y_s_1*i));%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x_s_1 + (B_y_s_1*i));%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x_s_2 - (B_y_s_2*i));%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x_s_1 + (B_y_s_1*i));%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x_s_2 + (B_y_s_2*i));%*RF_com_n;
    
    %{
    H(1,2) = 1/2*RF_com*gamma_e*(B_x_s_2/2 - (B_y_s_2*i)/2);
    H(1,3) = 1/2*RF_com*gamma_e*(B_x_s_1/2 - (B_y_s_1*i)/2);
    H(2,1) = 1/2*RF_com'*gamma_e*(B_x_s_2/2 + (B_y_s_2*i)/2);
    H(2,4) = 1/2*RF_com*gamma_e*(B_x_s_1/2 - (B_y_s_1*i)/2);
    H(3,1) = 1/2*RF_com'*gamma_e*(B_x_s_1/2 + (B_y_s_1*i)/2);
    H(3,4) = 1/2*RF_com*gamma_e*(B_x_s_2/2 - (B_y_s_2*i)/2);
    H(4,2) = 1/2*RF_com'*gamma_e*(B_x_s_1/2 + (B_y_s_1*i)/2);
    H(4,3) = 1/2*RF_com'*gamma_e*(B_x_s_2/2 + (B_y_s_2*i)/2);
    %}    
    H_p_rf = kron(eye(4),H);
    dUdt   = -i * 2*pi*  H_p_rf * U ; 
    %U_m = reshape(U,[],4);
    %dUdt_m = -i * H_p_rf * U_m ;
    %dUdt   = reshape(dUdt_m,[],1);
 
end


function J = interaction(t,t_f) 
    J = 60*sin(pi*t/t_f) ;    
end

function B_x_s_1 = Control_MF_x_1 (B_0,gamma_e,A_1, A_2, omega_x1_arr, t, t_f)
     
    w_dr = B_0*gamma_e + (A_2 - A_1)/4;           
    B_x_s_dr_1 = 0 ;
        for n = 1:length(omega_x1_arr)
            B_x_dr_1 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f) ;
            B_x_s_dr_1 = B_x_s_dr_1 + B_x_dr_1;
        end
    B_x_s_1 = B_x_s_dr_1* cos(w_dr*2*pi*t)  ; 
end


function B_y_s_1 = Control_MF_y_1 (B_0,gamma_e,A_1, A_2, omega_y1_arr, t, t_f)
    
    w_dr = B_0*gamma_e + (A_2 - A_1)/4;          
    B_y_s_dr_1 = 0 ;
        for n = 1:length(omega_y1_arr)
            B_y_dr_1 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_1 = B_y_s_dr_1 + B_y_dr_1;
        end
    B_y_s_1 = B_y_s_dr_1*cos(w_dr*2*pi*t)  ; 
end

function B_x_s_2 = Control_MF_x_2 (B_0,gamma_e,A_1, A_2, omega_x1_arr, t, t_f)   
 
    w_ur = B_0*gamma_e + (A_2 - A_1)/4;
            
    B_x_s_dr_2 = 0 ;
        for n = 1:length(omega_x1_arr)
            B_x_dr_2 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f);
            B_x_s_dr_2 = B_x_s_dr_2 + B_x_dr_2;
        end
    B_x_s_2 = B_x_s_dr_2*cos(w_ur*2*pi*t)   ; 
end

function B_y_s_2 = Control_MF_y_2 (B_0,gamma_e,A_1, A_2, omega_y1_arr, t, t_f)
    
    w_ur = B_0*gamma_e + (A_2 - A_1)/4;
            
    B_y_s_dr_2 = 0 ;
        for n = 1:length(omega_y1_arr)
            B_y_dr_2 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_2 = B_y_s_dr_2 + B_y_dr_2;
        end
    B_y_s_2 = B_y_s_dr_2*cos(w_ur*2*pi*t)   ; 
    
end

%}
