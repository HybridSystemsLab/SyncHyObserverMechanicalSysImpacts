%% Parameters

lambda = 0.8;
rho = 0.1;

Ac = [0,1;0,-rho];
Ad = [-1,0;0,-lambda];
Hd = [1,0];

%% Polytopic embedding

taumin = 0;
taumax = 5;

Lambda = eig(Ac);
res = compute_residues(Ac,eye(2),eye(2));
F1 = res(1:2,:);
F2 = res(3:4,:);
betabornes = zeros(2,4);
for k=1:2
    l = Lambda(k);
    if isreal(l)
        betabornes(k,:) = [exp(l*taumin),exp(l*taumax),0,0]; % if multiplicity=1
    else
        betabornes(k,:) = [2*exp(real(l)*taumin)*cos(imag(l)*taumin),2*exp(real(l)*taumax)*cos(imag(l)*taumax),-2*exp(real(l)*taumin)*sin(imag(l)*taumin),-2*exp(real(l)*taumax)*sin(imag(l)*taumax)];
    end
end
% if 2 real eigenvalues
M = zeros(2,2,4);
ind = 1;
for k = 1:2
    for j=1:2
        M(:,:,ind) = F1*betabornes(1,k)+F2*betabornes(2,j);
        ind = ind+1;
    end
end
sizeM = size(M);
P = sdpvar(2,2);
Ltilde = sdpvar(2,1);
constraints = [P>=0];
for k=1:sizeM(3)
    cons = [P,M(:,:,k)'*(P*Ad-Ltilde*Hd)';(P*Ad-Ltilde*Hd)*M(:,:,k),P];
    constraints = [constraints,cons<=0];
end
sol = optimize(constraints);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 Ld = inv(value(P))*value(Ltilde)
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
