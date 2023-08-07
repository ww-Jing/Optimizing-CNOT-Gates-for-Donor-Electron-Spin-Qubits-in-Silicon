

function [xmin,fmin,ct]=OptimalControl_leackage_Two_atom(myfunction,tol,max_feval,ain,bin,cin)
global xinit;
global B0 gamma_e gamma_n A10 A20 tunnel10 tunnel20 U_epsilon t0 t X ;
B0 = 1;
gamma_e = 1.760859644*10^5;%1.760859644*10^11
gamma_n = 10.8291;%108.291*10^5
i = sqrt(-1);
%%coefficient
A10 = 1.8199*10^2 ; %1.8199*10^8
%117*10^5
A20 = 1.8199*10^2 ; 
U_epsilon = 10.591*10^6;%10.591*10^12 Hz %/h
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

tol = 10^-9;
max_feval = 1000000;

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

Fidelity_Container = [];

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
    
    Fidelity_Container = [Fidelity_Container F_arr(1)];
    plot(Fidelity_Container)
end
xmin=X_arr(1,:)
fmin=F_arr(1)
%plot(Fidelity_Container)

%Fidelity = [];
%Fidelity = [Fidelity new_value];
%plot(Fidelity)

%hold on
%disp(myfunction)
%plot(x,y,'o')  
end

%{
function fidelity(Tsol,PSIsol)
Tsol = transpose(Tsol);
PSIsol = transpose(PSIsol);
norm_fact = sqrt(sum(PSIsol.*conj(PSIsol),1));
PSIsol = PSIsol ./ norm_fact;
PSIf = sqrt((PSI0'*PSIsol).*conj(PSI0'*PSIsol));
PSIf
figure(1);
plot(Tsol,PSIf,'b');

end
%}

function Fidelity = solve_ode2()
global B0 gamma_e gamma_n A10 A20 tunnel10 tunnel20 U_epsilon dEz  t0 aEz t;
B0 = 1;
gamma_e =  1.760859644*10^5;%1.760859644*10^11
gamma_n = 10.8291;%108.291*10^5
i = sqrt(-1);
%%coefficient
A10 = 1.8199*10^2 ; %1.8199*10^8
%117*10^5
A20 = 1.8199*10^2 ; 
U_epsilon = 10.591*10^6;%10.591*10^12 Hz %/h
t0 = 10^(-3); 
t = 1; 

U0_e1 = [1;0];
U0_e2 = [0;1];
U0_n1 = [1;0];
U0_n2 = [0;1];
PSI0 = Cross_4_Basis(U0_e1,U0_e2,U0_n1,U0_n2);    
[Tsol,PSIsol] = ode45(@Schrodinger,[0,10*t0],PSI0);

Fidelity = 1-(abs( PSIsol(end,:) * conj(PSI0) ))^2/16;

disp(Fidelity)

%disp(Fidelity);

end

function dPSIdt = Schrodinger(t,PSI)
global B0 gamma_e gamma_n A10 A20 tunnel10 tunnel20 U_epsilon t0 ;

i = sqrt(-1);
A1 = ElectricFieldA1(t);
A2 = ElectricFieldA2(t);
tunnel = ElectricFieldtunnel(t);

Sigma_x = [0 1;  1  0] / 2;
Sigma_y = [0 i; -i  0] / 2;
Sigma_z = [1 0;  0 -1] / 2;
I2 = eye(2);

S1_x = Cross_4_Basis(Sigma_x ,I2 ,I2 ,I2);
S1_y = Cross_4_Basis(Sigma_y ,I2 ,I2 ,I2);
S1_z = Cross_4_Basis(Sigma_z ,I2 ,I2 ,I2);

S2_x = Cross_4_Basis(I2 ,Sigma_x ,I2 ,I2);
S2_y = Cross_4_Basis(I2 ,Sigma_y ,I2 ,I2);
S2_z = Cross_4_Basis(I2 ,Sigma_z ,I2 ,I2);

I1_x = Cross_4_Basis(I2 ,I2 ,Sigma_x ,I2);
I1_y = Cross_4_Basis(I2 ,I2 ,Sigma_y ,I2);
I1_z = Cross_4_Basis(I2 ,I2 ,Sigma_z ,I2);

I2_x = Cross_4_Basis(I2 ,I2 ,I2 ,Sigma_x);
I2_y = Cross_4_Basis(I2 ,I2 ,I2 ,Sigma_y);
I2_z = Cross_4_Basis(I2 ,I2 ,I2 ,Sigma_z);

S1I1 = S1_x*I1_x + S1_y*I1_y + S1_z*I1_z;
S2I2 = S2_x*I2_x + S2_y*I2_y + S2_z*I2_z;
%S1S2 = S1_x*S2_x + S1_y*S2_y + S1_z*S2_z;

H = gamma_e*B0*(S1_z+S2_z) + gamma_n*B0*(I1_z+I2_z) + A1*S1I1 + A2*S2I2 ;


%extend matrix
H_leackage = zeros(20,20);
%%column
Sc = zeros(16,4);
Sc_ele = eye(4);
Sc_ele_pt = tunnel*Sc_ele;
Sc_ele_nt = (-1)*tunnel*Sc_ele;

Sc(5:8,1:4) = Sc_ele_pt;
Sc(9:12,1:4) = Sc_ele_nt;
%%row
Sr = zeros([4,20]);
Sr_ele1 = eye(4);
Sr_ele1_pt = tunnel*Sr_ele1;
Sr_ele1_nt = (-tunnel)*Sr_ele1;
Sr_ele2 = (U_epsilon)*(eye(4));

Sr(1:4,5:8) = Sr_ele1_pt;
Sr(1:4,9:12) = Sr_ele1_nt;
Sr(1:4,17:20) = Sr_ele2;

H_leackage(1:16,17:20) = Sc ;
H_leackage(17:20,1:20) = Sr ;
H_leackage(1:16,1:16) = H ;

%matrix S
pgamma = (tunnel/(U_epsilon-(A1-A2)/2)) ;
ngamma = (tunnel/(U_epsilon+(A1-A2)/2)) ;
psigma = (tunnel/(U_epsilon-(A1-A2)/2)) ;
nsigma = (tunnel/(U_epsilon+(A1-A2)/2)) ;

S = zeros([20,20]);
S(4,5:8) = -ngamma ;
S(5,5:8) = -nsigma ;
S(4,9:12) =  ngamma ;
S(5,9:12) =  nsigma ;

S(5:8,4) =  ngamma ;
S(5:8,5) =  nsigma ;
S(9:12,4) = -ngamma ;
S(9:12,5) = -nsigma ;

H_e = expm(S) * H_leackage * expm(-S);

SW_C = H_e(:,17:20);
SW_L = H_e(17:20,:);

H_e(:,17:20) = [];
H_e(17:20,:) = [];

%disp(H_e);
%disp(SW_C);
%disp(SW_L);

%RWA H = U+ * H * U - i * U+ * dU
Un0 = zeros(16,16);
Un0_1 = exp(-i*(A1+A2)/2*2*pi*t) * eye(4);
Un0_23 = eye(8);
Un0_4 = exp(i*(A1+A2)/2*2*pi*t) * eye(4);
Un0(1:4,1:4) = Un0_1;
Un0(5:12,5:12) = Un0_23;
Un0(13:16,13:16) = Un0_4;

Un0_0 = eye(4);
[tsol,Un0sol] = ode45(@diffU0,[0,10*t0],Un0_0);
Un0sol = Un0sol(26:41,1:16);

H_sw_rf = ctranspose(Un0) * H_e * Un0 - i  * ctranspose(Un0) * Un0sol;

%disp(H_sw_rf);

dPSIdt = i * H_sw_rf * PSI;
%disp(H_sw_rf);

end


function dUn0dt = diffU0(t,Un0)
t=1;
i = sqrt(-1);
A1 = ElectricFieldA1(t);
A2 = ElectricFieldA2(t);

Un0 = zeros(16,16);
Un0_1 = exp(-i*(A1+A2)/2*2*pi*t) * eye(4);
Un0_23 = eye(8);
Un0_4 = exp(i*(A1+A2)/2*2*pi*t) * eye(4);                 
Un0(1:4,1:4) = Un0_1;
Un0(5:12,5:12) = Un0_23;
Un0(13:16,13:16) = Un0_4;

dUn0dt = zeros(16,1);

end

function A1 = ElectricFieldA1(t)
   global A10 t0 X;
   a11 = X(1);
   a12 = X(2);
   %a13 = X(3);
   %a21 = X(4);
   %a22 = X(5);
   %a23 = X(6);
   A1 = A10 +(a11*cos(pi()/100/t0 * t)+a12*cos(2*pi()/100/t0 * t));
   %+a13*cos(3*pi()/100/t0 * t) )-...
   %(a21*cos(4*pi()/100/t0 * t)+a22*cos(5*pi()/100/t0 * t)+a23*cos(6*pi()/100/t0 * t) );
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
   %+b13*cos(3*pi/100/t0 * t))-...
   %(b21*cos(4*pi/100/t0 * t)+b22*cos(5*pi/100/t0 * t)+b23*cos(6*pi/100/t0 * t));
end

function tunnel = ElectricFieldtunnel(t)
   global t0 X;
   c11 = X(5);
   c12 = X(6);
   %c13 = X(15);
   %c21 = X(16);
   %c22 = X(17);
   %c23 = X(18);
   tunnel = (c11*cos(pi/100/t0 * t)+c12*cos(2*pi/100/t0 * t));
   %+c13*cos(3*pi/100/t0 * t))-...
   %(c21*cos(4*pi/100/t0 * t)+c22*cos(5*2*pi/100/t0 * t)+c23*cos(6*pi/100/t0 * t));
end

function  mat_full_basis = Cross_4_Basis(mat_e1 ,mat_e2 ,mat_n1 ,mat_n2)
    mat_full_basis = kron(kron(kron(mat_e1 ,mat_e2 ),mat_n1 ),mat_n2);
end




