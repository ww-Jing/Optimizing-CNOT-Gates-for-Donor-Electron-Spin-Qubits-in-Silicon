clear vars;
clear all; 
clc;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
    if isempty(odeget(opts,'Stats'))
        odeset(opts,'Stats','on');
    end
omega = 184;

delta_E = 200;
Parameters =[0.03, sqrt((2*pi*omega)^2 - delta_E^2)]; %565.7 . 1142.63
t_i = 0;    t_f =Parameters(1);
%beta_t = linspace(t_i,t_f,100);
%beta_s = beta_arr(beta_t);

t = linspace(t_i,t_f,3000);     
Usol = zeros(16,length(t)); 
h=(t_f - t_i)/3000;  
Uud_ud_arr_num = zeros(1,length(t));

k = 1000;
for tol = 1: k
U_i = eye(4);   
U_i_r = reshape(U_i , [] , 1);    Usol(:,1) = U_i_r; 

gamma_ou_arr = linspace(10^-8,1,25);
for ga = 1: length(gamma_ou_arr)
    
    gamma_ou_t = gamma_ou_arr(ga);
    beta_s = beta_arr(t, gamma_ou_t);
    %U_n+1 = f ( t_n, U_n )
    dU = @(t,U, beta_s) schrodinger(t,U, beta_s, Parameters);
for m=1:(length(t)-1)     
   
    k_1 = dU(t(m),Usol(:,m),beta_s(m));
    k_2 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_1),beta_s(m));
    k_3 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_2),beta_s(m));
    k_4 = dU((t(m)+h),(Usol(:,m)+h*k_3 ),beta_s(m));
    
    Usol(:,m+1) = Usol(:,m) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*h;  % main equation
end

%{
for m=1:(length(t)-1)     
    beta_s = beta_s(m);
    dU = @(t,U) schrodinger(t,U, beta_s, Parameters);

    k_1 = dU(t(m),Usol(:,m));
    k_2 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_1 ));
    k_3 = dU((t(m)+0.5*h),(Usol(:,m)+0.5*h*k_2 ));
    k_4 = dU((t(m)+h),(Usol(:,m)+h*k_3 ));
    
    Usol(:,m+1) = Usol(:,m) + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*h;  % main equation
end
%}
Uud_ud_arr = [];
for r = 1: length(t)
    Uf_i = Usol(:,r);   
    Uf = reshape(Uf_i,[],4);
    Uud_ud = real(Uf(2,3)*conj(Uf(2,3)));
    Uud_ud_arr = [Uud_ud_arr Uud_ud] ;   
end
Uud_ud_arr_num = Uud_ud_arr_num + Uud_ud_arr ;
end
end
plot(t,Uud_ud_arr_num/k/length(gamma_ou_arr));

function dU = schrodinger(t,U, beta_s, Parameters )
t_f = Parameters(1);%360/28000/3;%0.0628;%pi/50;
%B_x = 0; B_y = 0;
J_arr = Parameters(2:length(Parameters));
J_0 = interaction(J_arr ,t, t_f);
J = J_0  *( 1 + beta_s ) ;
B_x = 0;         B_y = 0;
gamma_e = 28000; B_0 = 2; E_0 = gamma_e*B_0;
Delta_E = 200;%MHz

H = zeros(4);
i = sqrt(-1);
    H(1,1) = Delta_E/2; %Delta_E/2;%E_0 + (Delta_E/2) ;    
    H(4,4) = -Delta_E/2; %-Delta_E/2;%-E_0 - (Delta_E/2) ;    
    H(2,2) =  Delta_E/2 - J/2 ; %+ Delta_E/2 - J/2 ; %-J/2;    
    H(3,3) =  - Delta_E/2 - J/2; %- Delta_E/2  - J/2 ; %- J/2;
    H(2,3) =  J/2;    
    H(3,2) =  J/2;
 %{
    H(1,2) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(1,3) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(2,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(2,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(3,1) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(3,4) = 1/2*gamma_e*(B_x - (B_y)*i)/2;%*RF_com_p;
    H(4,2) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n;
    H(4,3) = 1/2*gamma_e*(B_x + (B_y)*i)/2;%*RF_com_n; 
%}          
H_un = kron(eye(4),H);    
dU = -i*H_un*U;
end

function beta_s = beta_arr(t, gamma_ou_t )
%beta_s = 1*ones(1,length(t));; %zeros(1,length(t));%
%beta_s = randn()*0.1*ones(1,length(t));

sigma_ou_t = 0.145;         
beta_ss_t = [randn()*sigma_ou_t];
for m = 1 : length(t)-1
    dW = randn()*sqrt(t(m+1)-t(m));
    beta_ss_t(m+1) = (1-gamma_ou_t*(t(m+1)-t(m))) * beta_ss_t(m) ...
        + sigma_ou_t*sqrt(2*gamma_ou_t)*dW;
end


beta_s  =  beta_ss_t ;

end

function J = interaction(J_arr , t, t_f)
%tunnel = 4300;
ep = 0 ; 
for n = 1:length(J_arr)
ep_1 = J_arr(n);%*(sin((2*n-1)*pi*t/t_f)).^3 ; 
ep = ep +ep_1;
end
J = ep; %ep/2 + sqrt((ep/2)^2 + tunnel^2);
end

%/opt/matlab2019a/bin/matlab
