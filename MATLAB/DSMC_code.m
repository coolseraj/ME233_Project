%% DSMC CODE
clear all
clc
%%
%read data
und(1) = 100;
und(2) = 150;
und(3) = 250;
und(4) = 94;
MNM = sum(und);


%inputs for each location
mu = 0;
gym = 0.2;
library = 0.5; 
dining = 1;
sigma.gym = gym;
sigma.dining = dining;
sigma.library = library; 
num_p = round(num_location(sigma,mu, MNM));
groups = each_year(und, num_p);

I_matrix = [0.5 0.78 0.15 0.36; 0.78 0.3 0.98 0.56; 0.15 0.98 0.5 0.6; 0.36 0.56 0.6 0.1];
%%

%{
MNM: max # of molecules
MNC: max # of cells
MNSC: max # of sub-cells
%}
%MNM = 500;
MNC = 1;
MNSC = 3;

%{
NM: number of molecules
PP(MNM): x coordinate of molecule M
PV(1:3,MNM): u,v,w velocity comp of molecule M
IP(MNM): sub-cell number of molecule M
         Input: molecule #
         Output: sub-cell it is located at
IR(MNM): cross-reference array, molecule number in order of sub-cells
         ex/ the 7th molecule in the order of sub-cells is molecule #98
%}
NM = 0;
PP = zeros(1,MNM);
PV = ones(1,MNM)*1.7;
IP = zeros(1,MNM);
IR = zeros(1,MNM);
infected  = zeros(1,MNM);
susceptible = ones(1,MNM);
for i = 1:MNM
    % assign initial infected people distribution 
    if rand > 0.5
        infected(i) = 1;
    end
end

%{
IC(2,MNC)
    1: start address-1 of molecule # in array IR
    2: # of molecules in the cell
ISC(MNSC): the cell in which the sub-cell lies
ISCG(2,MNSC): index info on subcell
    1: start address-1 of molecule
    2: # of molecules in the sub-cell
%}
IC = zeros(2,MNC);
ISC = zeros(1,MNSC);
ISCG = zeros(2,MNSC);

%{
CS(5,MNC)
    1: # in the sample
    2,3,4: sum of u,v,w
    5: sum of u^2+v^2+w^2
FND: # density
FNUM: weight factor
DTM: time step
%}
DTM = 1e-1;

%{
CW: cell width
NSC: # of sub-cells per cell
XF,XR: min and max x-coordinate
%}

CW = 1; %CW = (XR - XF)/MNC;
XF = 0;
XR = CW*MNC + XF;
NSC = MNSC/MNC;


CG = zeros(2,MNC);

% INITIALIZATION
for M = 1:MNC
    if M > 1
        CG(1,M) = CG(2,M-1);
    end
    CG(2,M) = CG(1,M) + CW;
%     CCG(2,M) = 0;
%     CCG(1,M) = 1.25e-16; %initialize it to a reasonable low value
end

% % SET SUBCELLS 
for N = 1:MNC
    for M = 1:NSC
        L = (N-1)*NSC + M;
        ISC(L) = N;
    end
end

% GENERATE INITIAL GAS IN EQUILIBRIUM AT REF_TEMP
A = MNM/MNC; %FND*CW/FNUM; %# of simulated molecules in each cell
%NM = 0; it is already initialized to zero
for N = 1:MNC
    for M = 1:A
        R = rand; %random # between 0,1
        NM = NM + 1;
        PP(NM) = CG(1,N) + R*(CG(2,N) - CG(1,N));
        IP(NM) = NSC*(N-1) + ceil(R*NSC);
    end
end

c = 1;

while c<500
% MOVE MOLECULES - no infections happen here
for N = 1:NM
     %at each time step, only people in one SC can infect each other and
     %no infection happens when they move from one SC to another
     %next time step(next day, we initialize the people either randomely or from a dist
     %and assign new SC to them)
    MSC = IP(N);
    MC = ISC(MSC); %initial cell number
    XI = PP(N);
    DX = PV(1,N)*DTM;
    X = XI + DX;
    if X<XF || X>XR
    [X, PV_new] = BC(XF,XR,X,PV(1,N));
    PV(1,N) = PV_new;
    end
    if X<CG(1,MC) || X>CG(2,MC) %if molecule has moved from its initial cell
        MC = ceil((X - XF)/CW); %new cell #
    end
    IP(N) = NSC*(MC-1) + ceil((X - CG(1,MC))/CW*NSC); %new sub-cell #
    PP(N) = X;
end
% INDEXING 
% molecule # are arranged in the order of cells and within the cells, in order of sub-cells
for N = 1:NM
    MSC = IP(N);
    ISCG(2,MSC) = ISCG(2,MSC) + 1;
    MC = ISC(MSC);
    IC(2,MC) = IC(2,MC) + 1;
end
M = 0;
for N = 1:MNC %start address - 1 has been set for the cells
    IC(1,N) = M;
    M = M + IC(2,N);
end
M = 0;
for N = 1:MNSC %start address - 1 has been set for sub-cells
    ISCG(1,N) = M;
    M = M + ISCG(2,N);
    ISCG(2,N) = 0;
end
for N = 1:NM
    MSC = IP(N);
    ISCG(2,MSC) = ISCG(2,MSC) + 1;
    k = ISCG(1,MSC) + ISCG(2,MSC);
    IR(k) = N;
end

% COLLISION
cell = 1;
for N = 1:MNSC %goes through each collision cell
    AVN = ISCG(2,N); %# of molecules in a sub-cell
    if AVN>1
    ASEL = AVN*(AVN - 1)/2; %# of pairs to be selected in a sub-cell eq.11.3
    for ISEL = 1:ASEL
        kk = ceil(rand*ISCG(2,N)) + ISCG(1,N);
        L = IR(kk); %molecule L from collision cell N has been randomely chosen
        choose = 0;
        while choose == 0
        % choose the second molecule at ronadom from molecules that are in sub-cell MSC
        kk = ceil(rand*ISCG(2,N)) + ISCG(1,N);
        M = IR(kk);
        if M == L
            choose = 0;
        else
            choose = 1;
        end
        end
        if abs(PP(L) - PP(M)) <  0.001 %collision is accepted
            [year_L, year_M] = check_year(L, M, und);
            if infected(L) == 1 && infected(M) == 0
                infected(M) = round(I_matrix(year_L,year_M)); %it's either zero or one based on the probability
            end
            if infected(L) == 0 && infected(M) == 1
                infected(L) = round(I_matrix(year_L,year_M)); %it's either zero or one based on the probability
            end
        end
    end
    end
%      if rem(c,2)==0 %sample for the cell
%          for i = 1:MNM
%              if infected(i) == 1 && rand < accuracy 
%                  susceptible(i) = 0;
%              end
%          end
%      end
%      if rem(c,14)==0 %sample for the cell
%          for i = 1:MNM
%              if susceptible(i) == 0
%                  susceptible(i) = 1;
%              end
%          end
%      end


end
infected_tot(c) = sum(infected);
c = c + 1
IC = zeros(2,MNC);
ISCG = zeros(2,MNSC);
end
plot(infected_tot)
