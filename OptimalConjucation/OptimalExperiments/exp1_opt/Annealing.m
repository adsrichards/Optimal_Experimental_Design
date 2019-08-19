tic

NT = 64;    % Number of temperature points
Ns = 4;     % Number of state variables
No = Ns-1;  % Number of observable state variables
Np = 6;     % Number of parameters
Nm = 6;     % Number of measurements
Nk = 40;    % Number of measurement intervals

%------ MONTE CARLO SWEEPS (1 sweep = Nm steps) ------%
eqmcs = 1e4;  % equilibration Monte Carlo sweeps
msmcs = 1e5; % measurement Monte Carlo sweeps
grmcs = 1e4;  % greedy Monte Carlo sweeps

%------ TEMPERATURE GRID ------%
Tmin = 0.1;
Tmax = 120.0;
T = linspace(Tmax, Tmin, NT); % annealing grid (high temp -> low temp)
hT = (Tmax - Tmin)/NT;

%------ LOAD SENSITIVITIES ------%
disp("Loading sensitivities file (./Sensitivities/sensitivities1.txt)...");
S = load('Sensitivities/sensitivities1.txt');
Nt = size(S,2); % Number of time points
S = reshape(S,Ns,Np,Nt); % unflatten S
% computing average sensitivities on each of the Nk intervals
Sk = zeros(No,Np,Nk);
for k=1:Nk
    iL = ceil(Nt/Nk); % interval length
    Sk(:,:,k) = sum( S(1:No, :, iL*(k-1)+1:iL*k), 3)/iL;
end

%------ INITIALIZE STATE ------%
disp("Initializing state...");
[m,J] = Initialize(Sk,Nk,Nm,Np);

%------ ANNEALING ------%
disp(strcat("Starting annealing. Initial temperature = ", num2str(T(1)), "..."));
% measured quantities
J_ave = zeros(1,NT);
J2_ave = zeros(1,NT);
for Ti=1:NT
    tic
    %------ EQUILIBRATION LOOP ------%
    for i=1:eqmcs
        for mi=1:Nm
            [m,J] = MCStep(m,J,Sk,T(Ti),Nk,Nm,Np);
        end
    end

    %------ MEASUREMENT LOOP ------%
    disp(strcat("Taking measurements for temperature point ", num2str(Ti), "..."));
    for i=1:msmcs
        for mi=1:Nm
            [m,J] = MCStep(m,J,Sk,T(Ti),Nk,Nm,Np);
        end 
        J_ave(Ti) = J_ave(Ti) + J;
        J2_ave(Ti) = J2_ave(Ti) + J*J;
    end
    % normalizing
    J_ave(Ti) = -J_ave(Ti)/msmcs;
    J2_ave(Ti) = J2_ave(Ti)/msmcs;
    disp(strcat("Measurements finished for temperature point ", num2str(Ti)));
    toc
end
disp("Running greedy loop (temperature = 0)...")
%------ GREEDY LOOP ------%
for i=1:grmcs
    for mi=1:Nm
        T0=1e-8; [m,J] = MCStep(m,J,Sk,T0,Nk,Nm,Np);
    end
end

%------ DISPLAY OPTIMAL SCHEDULE ------%
disp("Annealing finished");
disp("The optimal samping schedule is")
disp(find(m==1));

%------ OUTPUT THERMODYNAMICS ------%
cV = (J2_ave - J_ave.*J_ave)./(T.*T);
dlmwrite(strcat("Thermodynamics/J_Nk=",num2str(Nk),"_msmcs=",num2str(msmcs),".txt"), J_ave, 'delimiter', ' ', 'precision', 6);
dlmwrite(strcat("Thermodynamics/cV_Nk=",num2str(Nk),"_msmcs=",num2str(msmcs),".txt"), cV, 'delimiter', ' ', 'precision', 6);
disp("Plotting thermodynamics - check thermodynamics plots for bad statistics.");

%------ PLOTTING THERMODYNAMICS ------%
plotting(J_ave,J2_ave,T);
disp("done");
toc

%====== FUNCTIONS ======%

%------ INITIALIZE STATE ------%
function [m,J] = Initialize(Sk,Nk,Nm,Np)
    m = zeros(1,Nk);
    m(1:Nm) = 1;
    mS = Sk(:,:,(find(m==1)));
    FIM = zeros(Np,Np);
    for mi=1:Nm
       FIM = FIM + mS(:,:,mi)'*mS(:,:,mi);
    end
    J = -det(FIM);
end

%------ MC STEP ------% 
%This is where all of the computational cost comes from
function [m,J] = MCStep(mo,Jo,Sk,T,Nk,Nm,Np)
    m1 = find(mo==1);  r1 = ceil(Nm*rand);       k1 = m1(r1); % random one index
    m0 = find(mo==0);  r0 = ceil((Nk-Nm)*rand);  k0 = m0(r0); % random zero index
    
    n = mo;  n(k0) = 1; n(k1) = 0; % update proposal
   
    nS = Sk(:,:,(find(n==1)));
    FIM = zeros(Np,Np);
    for i=1:Nm
       FIM = FIM + nS(:,:,i)'*nS(:,:,i);
    end
    
    delJ = -det(FIM) - Jo;
    
    if delJ < 0
        J = Jo + delJ;
        m = n;
        
    elseif rand < exp(-delJ/T)
        J = Jo + delJ;
        m = n;
    else
        J = Jo;
        m = mo;
    end  
end

%------ PLOTTING ------%
function plotting(J_ave, cV, T)
    figure(1)
    clf;
    hold on

    plot(-T, J_ave, '-o', 'Color', [36, 7, 133]/256,'LineWidth',1.5);
    %plot(, , 'o', 'Color', 'k','LineWidth',1.5);
    ax = gca;
    ax.FontSize = 48; 
    xlabel('T','FontSize',48)
    ylabel('\langle{J}\rangle','FontSize',48)
    %set(gca, 'XScale', 'log');
    axis([-120 0  0 600]);
    xticks([-120 -80 -40 0])
    xticklabels({'120' '80' '40' '0'})
    yticks([0 200 400 600])
    yticklabels({'0' '200' '400' '600'})
    %ylims = ylim;
    %text(2,60, strcat('c*=', num2str(I(1))), 'FontSize',32, 'fontweight','bold');
    %text(2,30, strcat('J(c*)=', num2str(253.4)), 'FontSize',32, 'fontweight','bold');
    %legend('exact','fit','Location','northeast')
    box on
    
    figure(2)
    clf;
    hold on
    
    plot(-T, cV, '-o', 'Color', [36, 7, 133]/256,'LineWidth',1.5);
    %plot(, , 'o', 'Color', 'k','LineWidth',1.5);
    ax = gca;
    ax.FontSize = 48; 
    ylabel('d\langle{J}\rangle/dT','FontSize',48)
    xlabel('T','FontSize',48)
    %set(gca, 'XScale', 'log');
    %axis([1 Nc 0 600]);
    xticks([-120 -80 -40 0])
    xticklabels({'120' '80' '40' '0'})
    %yticks([0 100 200 300])
    %yticklabels({'0' '100' '200' '300'})
    %ylims = ylim;
    %text(2,60, strcat('c*=', num2str(I(1))), 'FontSize',32, 'fontweight','bold');
    %text(2,30, strcat('J(c*)=', num2str(253.4)), 'FontSize',32, 'fontweight','bold');
    %legend('exact','fit','Location','northeast')
    box on
    
    figure(3)
    clf;
    hold on
    
    T1 = fliplr(T);
    cV1 = fliplr(cV);
    
    NT = size(T,2);
    hT = T(1)-T(2);
    ST = zeros(1,NT);
    for i=2:size(T,2)
       ST(i) = ST(i-1) + hT*cV1(i-1)/T1(i-1);
    end
    ST = fliplr(ST);
    
    plot(-T, ST, '-o', 'Color', [36, 7, 133]/256,'LineWidth',1.5);
    %plot(, , 'o', 'Color', 'k','LineWidth',1.5);
    ax = gca;
    ax.FontSize = 48; 
    ylabel('S(T)','FontSize',48)
    xlabel('T','FontSize',48)
    %set(gca, 'XScale', 'log');
    axis([-120 0 0 16]);
    xticks([-120 -80 -40 0])
    xticklabels({'120' '80' '40' '0'})
    yticks([0 5 10 15])
    yticklabels({'0' '5' '10' '15'})
    %ylims = ylim;
    %text(2,60, strcat('c*=', num2str(I(1))), 'FontSize',32, 'fontweight','bold');
    %text(2,30, strcat('J(c*)=', num2str(253.4)), 'FontSize',32, 'fontweight','bold');
    %legend('exact','fit','Location','northeast')
    box on
end
