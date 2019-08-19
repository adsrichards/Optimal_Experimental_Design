Nt = 1000;  % Number of time points
Ns = 4;     % Number of state variables
No = Ns-1;  % Number of observable state variables
Np = 6;     % Number of parameters
Nm = 6;     % Number of measurements
Nk = 40;    % Number of measurement intervals

%------ TIME GRID ------%
tmin = 0;
tmax = 360;
t = linspace(tmin, tmax, Nt);
h = (tmax - tmin)/Nt;

%------ SENSITIVITY ANALYSIS ------%
S = load('~/Documents/Project/OED Code/Conjugation/exp2_opt/sensitivities2.txt');

% average sensitivities on Nk intervals
Sk = zeros(No*Nk, Np);

kn = 0;
for k=1:Nt
    for s=1:No
        Sk(s + No*kn,:) = Sk(s + No*kn,:) + S(s + No*(k-1),:);
    end
    if mod(k,ceil(Nt/Nk)) == 0
        kn = kn + 1;  
    end
end
Sk = Sk/ceil(Nt/Nk);

function J = Objective(m)
    %------ CALCULATING J ------%
    J = zeros(1,Nc);
    for k=1:Nk
        if ismember(k, C(c,:))
            m(k) = 1;
        end
    end
    mS = zeros(No*Nm, Np);
    mk = 1;
    for k=1:Nk
        if m(k) == 1
            for s=1:No
                for p=1:Np
                    mS(s + No*(mk-1),p) = Sk(s + No*(k-1), p);            
                end
            end
        end
        mk = mk + 1;
    end

    J(c) = log(det(mS'*mS));
end
