

clc;
clear all;

delete(gcp('nocreate'));
%poolobj = parpool ;
%parpool local 2
%parpool(12);
%spmd

%Pauli_matrix
i = sqrt(-1);
Sigma_x = [0  1; 1  0]/2;
Sigma_y = [0 -i; i  0]/2;
Sigma_z = [1  0; 0 -1]/2;
Identity = eye(2);


%Control_Parameters
B_0 = 1; % 1T
gamma_e =28000;% 28000;%28024.95164;%B*10^6(Hz/T);
W_0 = B_0*gamma_e ; 

A_1 = 80;%1440;(MHZ)
A_2 = 120;%2420;
J  = 10%10;%(80)

%omega_x1_arr =1.0e-03 *[-0.226788535320793  -0.568599398997114  0.345220753516151   0.346271411320385 0.401016362667913  -0.254732852361272 0.822719501693017  -0.326517997259468 -0.557467980777712  -0.101131625389798 ];
%omega_y1_arr =1.0e-03 *[-0.212557410469814   0.193607681766298 0.345135442154544  -0.062986785135234  0.433495473477444  -0.270821713888270 0.982505583726377   0.290043569306544 -0.318220437082291  -0.360434923205752];
omega_x1_arr =1.0e-03 *[0.547563286161170   0.003779638026392  -0.021972895457217  -0.546738856577191  -0.612487262727591   0.134752210069026 -0.138123865955881  -0.736599330170021 ];
omega_y1_arr =1.0e-03 *[ 0.168138875336672   0.919665294156572  0.111571589945191  -0.572553897844082 0.791702880444480   0.829473591689472 -0.990928127352689  -0.190937505825745];
Parameters = [omega_x1_arr  omega_y1_arr ];% omega_x2_arr  omega_y2_arr ];

t_i = 0;
t_f = 0.1%0.6;
U_i = eye(4);
U_t = [ 1 0 0 0 ;0 1 0 0 ;0 0 0 1 ; 0 0 1 0 ];
   
    
U_i_r = reshape(U_i , [] , 1);
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
infidelity = 1-(abs(trace(Uf*U_t')))^2/16
 
 


function dUdt = Schrodinger_H_p_rf(t,U,Parameters)

    B_0 = 1;
    i = sqrt(-1);
    gamma_e = 28000;
    A_1 = 80;
    A_2 = 120;
    delta_A = (A_2 - A_1)/2;
    ave_A = (A_2 + A_1)/2;
    J  = 10;%10;
    t_f = 0.1;%0.6;
    
  
    
    omega_x1_arr = Parameters(1:8);
    omega_y1_arr = Parameters(9:16);
    %omega_x2_arr = Parameters(13:18);
    %omega_y2_arr = Parameters(19:24);
    B_x_s_1 = Control_MF_x_1(omega_x1_arr, t, t_f,B_0, gamma_e, delta_A);
    B_y_s_1 = Control_MF_y_1(omega_y1_arr, t, t_f,B_0, gamma_e, delta_A);
    B_x_s_2 = Control_MF_x_2(omega_x1_arr, t, t_f,B_0, gamma_e, delta_A);
    B_y_s_2 = Control_MF_y_2(omega_y1_arr, t, t_f,B_0, gamma_e, delta_A);
    
    H = zeros(4);
    
    H(1,1) = gamma_e*B_0 + delta_A/2 + J/4;
    H(4,4) = -gamma_e*B_0 - delta_A/2 + J/4;
    
    H(2,2) = ave_A/2 - J/4;
    H(3,3) = -ave_A/2 - J/4;
    
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
    
    H_p_rf = kron(eye(4),H);
    
    %U_m = reshape(U,[],4);
    dUdt = -i * 2*pi*  H_p_rf * U ;
    %dUdt_m = -i * 2*pi*  H_p_rf * U_m ;
    
    %dUdt   = reshape(dUdt_m,[],1);
 
end


function B_x_s_1 = Control_MF_x_1 (omega_x1_arr, t, t_f, B_0, gamma_e, delta_A)
       
    B_x_s_dr_1 = 0 ;
        for n = 1:8
            B_x_dr_1 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f) ;
            B_x_s_dr_1 = B_x_s_dr_1 + B_x_dr_1;
        end
    B_x_s_1 = B_x_s_dr_1* cos((2*pi*(B_0*gamma_e + delta_A/2))*t)  ; 
end


function B_y_s_1 = Control_MF_y_1 (omega_y1_arr, t, t_f, B_0, gamma_e, delta_A)
      
    B_y_s_dr_1 = 0 ;
        for n = 1:8
            B_y_dr_1 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_1 = B_y_s_dr_1 + B_y_dr_1;
        end
    B_y_s_1 = B_y_s_dr_1*cos((2*pi*(B_0*gamma_e + delta_A/2))*t)  ; 
end

function B_x_s_2 = Control_MF_x_2 (omega_x1_arr, t, t_f, B_0, gamma_e, delta_A)   
    
    B_x_s_dr_2 = 0 ;
        for n = 1:8
            B_x_dr_2 = omega_x1_arr(n) * sin((2*n-1)*pi*t/t_f);
            B_x_s_dr_2 = B_x_s_dr_2 + B_x_dr_2;
        end
    B_x_s_2 = B_x_s_dr_2*cos((2*pi*(B_0*gamma_e + delta_A/2))*t)   ; 
end

function B_y_s_2 = Control_MF_y_2 (omega_y1_arr, t, t_f, B_0, gamma_e, delta_A)
    
    B_y_s_dr_2 = 0 ;
        for n = 1:8
            B_y_dr_2 = omega_y1_arr(n) * sin((2*n)*pi*t/t_f) ;
            B_y_s_dr_2 = B_y_s_dr_2 + B_y_dr_2;
        end
    B_y_s_2 = B_y_s_dr_2*cos((2*pi*(B_0*gamma_e + delta_A/2))*t)   ; 
end


