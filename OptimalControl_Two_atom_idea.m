
%{
function [xmin,fmin,ct]=OptimalControl_Two_atom_idea(myfunction,tol,max_feval,ain,bin,cin)
global xinit;
global B0 gamma_e gamma_n A10 A20 t0 t X theda mu dalta_A ave_A E_ud_w E_du_w mu_S mu_T E_uu E_ud E_du E_dd B_1 ;
B0 = 1;
gamma_e = 1.760859644*10^5;%1.760859644*10^11
gamma_n = 10.8291;%108.291*10^5
i = sqrt(-1);
%%coefficient
A10 = 1.8199*10^2 ; %1.8199*10^8
%117*10^5
A20 = 1.8199*10^2 ; 
t0 = 10^(-3); 
t = 1; 


% input argument testing
myfunction = @solve_ode2;


a11in = 25*10^3; %25*10^5;
a12in = 0.32;
%a12in = 28; %28*10^3;
%a13in = 0.32; %320;
%a21in = -20*10^3; %-20*10^5;
%a22in = -25; %-25*10^3;
%a23in = -0.3; %-300;

b11in = 25*10^3; %25*10^5;
b12in = 0.32;
%b12in = 28; %28*10^3;
%b13in = 0.32; %320;
%b21in = -20*10^2; %-20*10^5;
%b22in = -25; %-25*10^3;
%b23in = -0.3; %-300;

c11in = 25*10^3; %25*10^6;
c12in = 0.32;
%c12in = 28; %28*10^3;
%c13in = 0.032; %32;
%c21in = -20*10^3; %-20*10^6;
%c22in = -25; %-25*13^3;
%c23in = -0.03; %-30;

tol = 10^-6;
max_feval = 10000;

% Order according to the values at the vertices
% Check convergency
%{
xinit = [a11in,a12in,a13in,a21in,a22in,a23in,...
    b11in,b12in,b13in,b21in,b22in,b23in,...
    c11in,c12in,c13in,c21in,c22in,c23in];
%}
xinit = [a11in,a12in,b11in,b12in,c11in,c12in];
x0 = xinit(:)';  % x0 is a row vector.

myfunction = fcnchk(myfunction);

dim = size(x0); % dimension of the problem
dim = dim(2); %取size of dim's second term.
X_arr = zeros(dim+1,dim);
F_arr = zeros(dim+1);
% set up adaptive parameters
alpha=1; beta=1+2/dim; gamma=0.75-0.5/dim; delta=1-1/dim;

% Construct the initial simplex: Large initial simplex is used.
scalefactor = min(max(max(abs(x0)),1),10);

D0 = eye(dim);
D0(dim+1,:)=(1-sqrt(dim+1))/dim*ones(1,dim);
for i=1:dim+1
   X = x0 + scalefactor * D0(i,:);
   X_arr(i,:) = X;
   F_arr(i) = feval(myfunction);
end

ct = dim+1;

% Create start point
[F_arr,Index]=sort(F_arr);
X_arr = X_arr(Index,:);

% Main iteration

Infidelity_Container = [];

while max(max(abs(X_arr(2:dim+1,:)-X_arr(1:dim,:)))) >= scalefactor*tol 
    if ct>max_feval
        break;
    end   
    M=mean(X_arr(1:dim,:));  % Centroid of the dim best vertices
    % FM=mean(FX(1:dim)); 
    X = (1+alpha)*M - alpha * X_arr(dim+1,:);
    xref = X;
    Fref = feval(myfunction);
    ct=ct+1;
    if Fref<F_arr(1)
        % expansion
        X = (1+alpha*beta)* - alpha*beta*X_arr(dim+1,:);
        Fexp = feval(myfunction);
        ct=ct+1;
        if Fexp < Fref
          X_arr(dim+1,:) = X;
          F_arr(dim+1)=Fexp;
        else
          X_arr(dim+1,:) = xref;
          F_arr(dim+1) = Fref;
        end
    else
        if Fref<F_arr(dim)
            % accept reflection point
            X_arr(dim+1,:) = xref;
            F_arr(dim+1)=Fref;
        else 
            if Fref<F_arr(dim+1)
             % Outside contraction
                xoc=(1+alpha*gamma)*M-alpha*gamma*X_arr(dim+1,:);
                X_arr = xoc;
                Foc=feval(myfunction);
                ct=ct+1;           
              if Foc<=Fref
                X_arr(dim+1,:)=xoc;
                F_arr(dim+1)=Foc;
              else
                % shrink
                for i=2:dim+1
                    X = X_arr(1,:)+ delta*(X_arr(i,:)-X_arr(1,:));
                    X_arr(i,:) = X;
                    F_arr(i)=feval(myfunction);
                end
                ct=ct+dim;
              end
            else
              %inside contraction
              xic=(1-gamma)*M+gamma*X_arr(dim+1,:);
              X = xic;
              Fic=feval(myfunction);
              ct=ct+1;
              if Fic<F_arr(dim+1)
                X_arr(dim+1,:)=xic;
                F_arr(dim+1)=Fic;
              else
                % shrink
                for i=2:dim+1
                    X = X_arr(1,:)+ delta*(X_arr(i,:)-X_arr(1,:));
                    X_arr(i,:) = X;
                    F_arr(i)=feval(myfunction);
                    
                end
                ct=ct+dim;
              end
            end
        end
    end
    [F_arr,Index] = sort(F_arr);
    X_arr = X_arr(Index,:);
    
    Infidelity_Container = [Infidelity_Container F_arr(1)];
    plot(Infidelity_Container,'o-')
end
xmin=X_arr(1,:)
fmin=F_arr(1)

plot(Fidelity_Container)

Fidelity = [];
Fidelity = [Fidelity new_value];
plot(Fidelity)

hold on
disp(myfunction)
plot(x,y,'o')  
end
%}
function Infidelity=OptimalControl_Two_atom_idea()
%function Infidelity = solve_ode2()
global B0 gamma_e gamma_n A10 A20 t t0;
B0 = 10^2;%10^5
gamma_e =  1.760859644*10^5;%1.760859644*10^11
gamma_n = 10.8291;%108.291*10^5
i = sqrt(-1);
%%coefficient
%A10 = 1.8199*10^2 ; %1.8199*10^8
%117*10^5
%A20 = 1.8199*10^2 ; 
%U_epsilon = 10.591*10^6;%10.591*10^12 Hz %/h
t0 = 10^(-3); 
t = 1; 

U0_e1 = [1;1];
U0_e2 = [0;1];
PSI0 = kron(U0_e1,U0_e2);    

[Tsol,PSIsol] = ode45(@Schrodinger,[0,10*t0],PSI0);

Infidelity = 1-(abs( PSIsol(end,:) * conj(PSI0) ))^2/16;
fprintf('Infidelity =');
disp(Infidelity);


end

%%

function H_2 = Hamitonian()
global B0 gamma_e A10 A20 t t0 theda mu dalta_A ave_A E_ud_w E_du_w mu_S mu_T E_uu E_ud E_du E_dd B_1 ;
i = sqrt(-1);
%A1 = ElectricFieldA1(t);
%A2 = ElectricFieldA2(t);
%J = ElectricFieldtunnel(t);
B_1 = 10^2;
A1 =117 ;%117*10^3
A2 =117 ;
J =10^3 ; %10^6

theda = 1/2*atan(J./ave_A);

dalta_A = (A1 - A2)/2;
ave_A = (A1 + A2)/2;
E_uu = gamma_e*B0 + J/4 +dalta_A/2;
E_ud = -J/4 - ave_A/2;
E_du = -J/4 + ave_A/2;
E_dd = -gamma_e*B0 + J/4 -dalta_A/2;


E_ud_w = -J/4 - (ave_A^2 + J^2)^(1/2)/2;
E_du_w = -J/4 + (ave_A^2 + J^2)^(1/2)/2;
mu_S = (cos(theda)-sin(theda))/2;
mu_T = (cos(theda)+sin(theda))/2;

mu = E_ud_w - E_dd;

H_2 = [(E_uu - mu) gamma_e*B_1*mu_S gamma_e*B_1*mu_T 0 ;...
    
    gamma_e*B_1*mu_S E_ud_w 0 gamma_e*B_1*mu_S;...
    gamma_e*B_1*mu_T 0 E_du_w gamma_e*B_1*mu_T;...
    0 gamma_e*B_1*mu_S gamma_e*B_1*mu_T (E_dd + mu)] ;
end

function dPSIdt = Schrodinger(t,PSI,H_sw,Un0,Un0sol)
global B0 gamma_e A10 A20 t0 theda mu dalta_A ave_A E_ud_w E_du_w mu_S mu_T E_uu E_ud E_du E_dd B_1 ;

H_2 = Hamitonian();
I = [1 0 ; 0 1];
U_rwa_l = [exp(-i*mu*t) 0;0 exp(i*mu*t)];
U_rwa_l4  = kron(U_rwa_l,I);
U_rwa_r = [exp(i*mu*t) 0; 0 exp(-i*mu*t)];
U_rwa_r4  = kron(U_rwa_r,I);

H_rwa = U_rwa_l4 * H_2 * U_rwa_r4;
dPSIdt = i * H_rwa * PSI;
end



%%
%{
function A1 = ElectricFieldA1(t)
   global A10 t0 X;
   a11 = X(1);
   a12 = X(2);
   %a13 = X(3);
   %a21 = X(4);
   %a22 = X(5);
   %a23 = X(6);
   A1 = A10 +(a11*cos(pi()/100/t0 * t)+a12*cos(2*pi()/100/t0 * t));
end

function A2 = ElectricFieldA2(t)
   global A20 t0 X;
   b11 = X(3);
   b12 = X(4);
   %b13 = X(9);
   %b21 = X(10);
   %b22 = X(11);
   %b23 = X(12);
   A2 = A20 + (b11*cos(pi/100/t0 * t)+b12*cos(2*pi/100/t0 * t));
end

function J = ElectricFieldtunnel(t)
   global X t0;
   c11 = X(5);
   c12 = X(6);
   %c13 = X(15);
   %c21 = X(16);
   %c22 = X(17);
   %c23 = X(18);
   J = (c11*cos(pi/100/t0 * t)+c12*cos(2*pi/100/t0 * t));
end
%}