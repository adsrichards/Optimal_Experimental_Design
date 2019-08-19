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
S = load('~/Documents/Project/OED Code/Conjugation/OptimalExperiments/exp6_opt/sensitivities6.txt');

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

%------ CALCULATING J ------%
C = combnk(1:Nk,Np);
Nc = size(C,1);
J = zeros(1,Nc);
for c=1:Nc
    m = zeros(1,Nk);
    for k=1:Nk
        if ismember(k, C(c,:))
            m(k) = 1;
        end
    end
    mS = zeros((Ns-1)*Nm, Np);
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

%------ SORTING ------%
[Js,I] = sort(J, 'descend');

disp("optimal schedule is");
disp(C(I(1),:));
%------ PLOTTING ------%
plotting(Js,I,Nc);

function plotting(Js,I,Nc)

    figure(1)

    clf;
    hold on
    cgrid = linspace(1,Nc,Nc);
    %cgrid2 = 10.^(linspace(0,7,50));
    %A = 445;
    %b = 0.35;
    %c = 60000;
    %strexp_fit = A*exp(-(log(cgrid)/c));
    %strexp_fit = A*exp(-((cgrid2-1)/c).^b);
    plot(cgrid, exp(Js), '-', 'Color', [36, 7, 133]/256,'LineWidth',2);
    %plot(cgrid2, strexp_fit, 'o', 'Color', 'k','LineWidth',1.5);
    ax = gca;
    ax.FontSize = 32; 
    xlabel('sampling schedule configuration','FontSize',32)
    ylabel('J(c)','FontSize',32)
    set(gca, 'XScale', 'log');
    %axis([1 Nc 0 600]);
    %xticks([1e0 1e2 1e4 1e6 1e8 1e10])
    %xticklabels({'10^0' '10^2' '10^4' '10^6' '10^8' '10^10'})
    %yticks([0 100 200 300])
    %yticklabels({'0' '100' '200' '300'})
    %ylims = ylim;
    text(2,60, strcat('c*=', num2str(I(1))), 'FontSize',32, 'fontweight','bold');
    text(2,30, strcat('J(c*)=', num2str(exp(Js(1)))), 'FontSize',32, 'fontweight','bold');
    legend('exact','Location','northeast')
    box on
end
