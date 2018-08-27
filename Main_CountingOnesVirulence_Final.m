%Main Code for Final Project
%Counting Ones Domain
%Simple Coevolutionary Algorithm
%Kristina Harper


clear;
NoRuns = 10;

%Initialize Variables
NoHost = 25; 
NoPar = 25;
NoGenes = 100;
NoGenerations = 600;    	
mutationRate= 0.03; 
TournSize = 5;
BPar = 0.9;%0.5;% 0.9; %Introduce parasite mutation bias. Varies (0-1). Cartlidge tested 0.9


StoreTotalOnesHost = zeros(NoRuns, NoGenerations); %To store for average across runs
StoreTotalOnesPar = zeros(NoRuns, NoGenerations);
StoreRelativeFitH = zeros(NoRuns,NoGenerations);
StoreRelativeFitP = zeros(NoRuns,NoGenerations);
MeanOnesH = zeros(1,NoGenerations);
MeanOnesP = zeros(1,NoGenerations);
MeanRelHostFit = zeros(1,NoGenerations);
MeanRelParFit = zeros(1,NoGenerations);
BestRelativeFitHost = zeros(1,NoGenerations);
BestRelativeFitPar = zeros(1,NoGenerations);

%Virulence 

%lambda = 1.0;    %Maximum Virulence 
%lambda = 0.75;   %Moderate Virulence
lambda = 0.5;    %Null Virulence
%lambda = 0.3;


%Initialize target string - genotype of ones
TargetString = ones(1,NoGenes);

%Start loop for multiple runs
for runs = 1:NoRuns
    TotalOnesHost = zeros(1,NoGenerations); %Moved from initialization to loop for 10 runs. 
    TotalOnesPar = zeros(1,NoGenerations);

%Need 2 populations of genotype zeros to start
    HostPop = zeros(NoHost,NoGenes);
    ParPop = zeros(NoPar,NoGenes);
    HostFitness = zeros(1,NoHost);  %Initialize array to hold ranked & normalized host fitness
    HostFit = zeros(1,NoHost);      %Initialize array to hold simple host fitness
    [HostFit,HostFitness] = SimpleHostFitness(HostPop,NoHost,HostFitness); %Call on host fitness function

    XScore = zeros(1,NoPar);        %Initialize xscore array each run
    ParasiteFit = zeros(1,NoPar);
    [ParFit] = SimpleParFitness(ParPop,NoPar); %Find simple parasite fitness by counting ones
    [ParasiteFit] = ParasiteFitness(NoPar,ParFit,lambda); %Set initial parasite 

    TournComp = zeros(1,TournSize);     %Tournament size 5


    %Main Tournament Loop
    for c = 1:NoGenerations
        RelHostFitness = zeros(1,NoHost);   %Initialize relative host fitness 
        RelParFitness = zeros(1,NoPar);     %Initialize relative parasite fitness
        
        for i = 1:NoHost
            for d = 1:TournSize
                TournComp(d) = floor(rand* NoPar + 1); %Select opponent indeces from parasite population
            end

            [maxFitCalc,maxIndex] = max(ParasiteFit(TournComp)); %Find fittest opponent
            if maxFitCalc > HostFit(i)      % Select parasite as winner
                Win = TournComp(maxIndex);  % if Winner is Parasite, win index is from tournament competitors
                HostWin = 0;                % Boolean value for application of mutation bias 
            elseif maxFitCalc == HostFit(i) % If equal fitnesses, randomly select Par or Host as winner
                if rand > 0.5
                    Win = TournComp(maxIndex); 
                    HostWin = 0;
                else
                    Win = i;                % If winner is host, win index from host population
                    HostWin = 1;
                end
            else
                Win = i; 
                HostWin = 1;
            end 

            for e = 1:TournSize             % Assign relative fitness to parasite and host individuals 
                if HostFit(i) > ParasiteFit(TournComp(e)) 
                    RelHostFitness(i) = RelHostFitness(i)+1; % one point if wins tournament
                elseif HostFit(i) == ParasiteFit(TournComp(e))
                    RelHostFitness(i) = RelHostFitness(i)+0.5; % 0.5 point to host and par if draw
                    RelParFitness(TournComp(e)) = RelParFitness(TournComp(e))+0.5;
                else
                    RelParFitness(TournComp(e)) = RelParFitness(TournComp(e))+1;
                end
            end

        %Mutation of winner given Mutation Bias
        %"Given mutation at particular parasite locus, substitution of 1 occurs with
        %BPar, substitution of 0 occurs with (1-BPar). Coevolving Host substitutes
        %with equal likelihood whenever mutation occurs." 
            for g = 1:NoGenes       % Determine if mutation occurs in each gene of the winning individual
                if HostWin == 1 && (rand < mutationRate) 
                    HostPop(Win,g) = not(HostPop(Win,g)); 
                elseif HostWin == 0 && (rand < mutationRate) 
                    if ParPop(Win,g) == 0 && (rand < BPar)
                        ParPop(Win,g) = not(ParPop(Win,g)); 
                    elseif ParPop(Win,g) == 1 && (rand < (1-BPar))
                        ParPop(Win,g) = not(ParPop(Win,g));
                    end
                end
            end

            %Recalculate fitness with Host Fitness functions 
            [HostFit,HostFitness] = SimpleHostFitness(HostPop,NoHost,HostFitness);

            %Recalculate fitness with Parasite Fitness functions
            [ParFit] = SimpleParFitness(ParPop,NoPar); %Find simple parasite fitness by counting ones
            [ParasiteFit] = ParasiteFitness(NoPar,ParFit,lambda); 
        end

        for i = 1:NoPar         % Same process as above with host individuals, but with parasite population. 
            for d = 1:TournSize
                TournComp(d) = floor(rand* NoHost + 1); %Select opponent indeces from host population
            end

            for e = 1:TournSize
                if ParasiteFit(i) > HostFit(TournComp(e)) 
                    RelParFitness(i) = RelParFitness(i)+1;
                elseif ParasiteFit(i) == HostFit(TournComp(e))
                    RelParFitness(i) = RelParFitness(i)+ 0.5;
                    RelHostFitness(TournComp(e)) = RelHostFitness(TournComp(e))+0.5;
                else
                    RelHostFitness(TournComp(e)) = RelHostFitness(TournComp(e))+1;
                end
            end

            [maxFitCalc,maxIndex] = max(HostFit(TournComp));
            if maxFitCalc > ParasiteFit(i) %select host as winner
                Win = maxIndex; 
                HostWin = 1;
            elseif maxFitCalc == ParasiteFit(i)
                if rand > 0.5
                    Win = i; %If equal fitnesses, randomly select par or host as winner
                    HostWin = 0;
                else
                    Win = maxIndex; 
                    HostWin = 1;
                end
            else
                Win = i; 
                HostWin = 0;
            end  

            for g = 1:NoGenes
                if HostWin == 1 && (rand < mutationRate)
                   HostPop(TournComp(Win),g) = not(HostPop(TournComp(Win),g)); 
                end 

                if HostWin == 0 && (rand < mutationRate) %if  Winner's in Parasite
                    if ParPop(Win,g) == 0 && (rand < BPar)
                        ParPop(Win,g) = not(ParPop(Win,g)); 
                    elseif ParPop(Win,g) == 1 && (rand < (1-BPar))
                        ParPop(Win,g) = not(ParPop(Win,g));
                    end
                end
            end
            %Find new fitness Host Pop
            [HostFit,HostFitness] = SimpleHostFitness(HostPop,NoHost,HostFitness);

            %Recalculate par's fitness %Call on ParasiteFitness function
            [ParFit] = SimpleParFitness(ParPop,NoPar); %Find simple parasite fitness by counting ones
            [ParasiteFit] = ParasiteFitness(NoPar,ParFit,lambda); 

        end 
        
        % Set up statistical tests 
        relHostFit(c) = (sum(RelHostFitness))/NoHost;   % Track relative host fitness across generations
        relParFit(c) = (sum(RelParFitness))/NoPar;      % Track relative parasite fitness across generations
        StoreRelativeFitH(runs,c) = relHostFit(c);      % Store relHostFit for averaging later
        StoreRelativeFitP(runs,c) = relParFit(c);       % Store relParFit for averaging later
%         relHostFit(c) = sum(RelHostFitness);
%         relParFit(c) = sum(RelParFitness);
%         meanHostFit(c) = mean(HostFit);
%         meanParFit(c) = mean(ParasiteFit); 
        MeanRelHostFit(c) = mean(RelHostFitness); %Same for host and par. 
        MeanRelParFit(c) = mean(RelHostFitness);
        BestRelativeFitHost(c) = max(RelHostFitness);
        BestRelativeFitPar(c) = max(RelParFitness);
        
        bestHostFitnessSim(c) = max(HostFitness);
        bestParFitnessSim(c) = max(ParFit);
            
        PartSum = sum(HostPop,2);        % Determine sum of all ones for each indiv in HostPop    
        %PartSumH(:,c) = sum(HostPop,2);
        TotalOnesHost(c) = sum(PartSum); % Total ones for entire population of Host
        StoreTotalOnesHost(runs,c) = TotalOnesHost(c);
       % TotalOnesHostIndivs(runs,c) = PartSumH(:,c);
        MeanOnesH(c) = mean(PartSum);   % Calculate mean number of ones in Host Population
        
        PartSum = sum(ParPop,2);    
        %PartSumP(:,c) = sum(ParPop,2);
        TotalOnesPar(c) = sum(PartSum); 
        StoreTotalOnesPar(runs,c) = TotalOnesPar(c);
        %TotalOnesParIndivs(runs,c) = PartSumP;
        MeanOnesP(c) = mean(PartSum);
    end      

% Set up figures

    figure(2)
    plot(bestHostFitnessSim,'b')
    hold on;
    plot(bestParFitnessSim,'r')
    hold on;
    title('Best Individual Fitness Per Generation')
    xlabel('Generations');ylabel('Fitness');
    legend('Host','Parasite')
    
%     figure(2)
%     plot(relHostFit,'b.')
%     hold on;
%     plot(relParFit,'r.')
%     hold on;
%     title('Relative Fitness')
%     xlabel('Generations');ylabel('Relative Fitness');
%     legend('Host', 'Parasite')
%     
%     figure(3)
%     plot(MeanRelHostFit, 'b')
%     hold on;
%     plot(MeanRelParFit,'r')
%     title('Mean Relative Fitness')
%     ylabel('relative fitness');xlabel('Relative Fitness')
%     legend('Host', 'Parasite')

%     figure(3)
%     plot(meanHostFit, 'b')
%     hold on;
%     plot(meanParFit,'r')
%     figure(3)
%     title('Average Fitness')
%     xlabel('Generations');ylabel('Average Fitness');
%     %ylim([0 1])
%     legend('Host', 'Parasite')%, 'Average Host', 'Average Parasite')

%     figure(4)
%     plot(TotalOnesHost,'b')
%     hold on;
%     plot(TotalOnesPar, 'r') 
%     hold on;
%     title('Total Ones Per Generation')
%     xlabel('Generations'); ylabel('Number of Ones');
%     xlim([0 NoGenerations])
%     legend('Host', 'Parasite')
%     
%     figure(5)
%     plot(MeanOnesH, 'b')
%     hold on;
%     plot(MeanOnesP,'r')
%     hold on;
%     title('Mean Ones Per Individual Per Generation')
%     legend('Host', 'Parasite')
%     xlabel('Generations');ylabel('Mean Ones');
%     ylim([0 NoGenes]);
end

% [AvgHostOnes] = AverageHostOnes(StoreTotalOnesHost,NoRuns);
% figure(4)
% LengthGen = 1:1:NoGenerations;
% plot(LengthGen,AvgHostOnes, 'k*', 'LineWidth', 2, 'MarkerSize', 2);
% hold on;
% 
% [AvgParOnes] = AverageParasiteOnes(StoreTotalOnesPar,NoRuns);
% figure(4)
% plot(LengthGen,AvgParOnes, 'y*', 'LineWidth', 2, 'MarkerSize', 2);
% hold on;
% 
% [AvgRelFitH,AvgRelFitP] = AverageRelativeFitHost(StoreRelativeFitP,StoreRelativeFitH,NoRuns);
% figure(2)
% plot(LengthGen, AvgRelFitH, 'k-*', 'MarkerSize', 2, 'LineWidth', 2)
% hold on;
% plot(LengthGen,AvgRelFitP,'y-*','Markersize',2, 'LineWidth', 2)
% 
% [MeanRelFitHost, MeanRelFitPar] = AverageRelativeFit(BestRelativeFitHost,BestRelativeFitPar,...
%     NoGenerations);
% figure(1)
% plot(LengthGen,MeanRelFitHost, 'k*', 'MarkerSize', 2)
% hold on;
% plot(LengthGen,MeanRelFitPar, 'y*', 'MarkerSize', 2)
% hold on;



%StandardDev(StoreTotalOnesHost,StoreTotalOnesPar,NoRuns, LengthGen, AvgHostOnes,AvgParOnes)

function [MeanRelFitHost, MeanRelFitPar] = AverageRelativeFit(BestRelativeFitHost,BestRelativeFitPar,...
    NoGenerations)
    SumFit = sum(BestRelativeFitHost);
    MeanRelFitHost = SumFit/NoGenerations;
    SumFit = sum(BestRelativeFitPar);
    MeanRelFitPar = SumFit/NoGenerations;
end

function [AvgHostOnes] = AverageHostOnes(StoreTotalOnesHost,NoRuns) %NoHost,NoGenerations)
    SumofOnes = sum(StoreTotalOnesHost);
    AvgHostOnes = SumofOnes/NoRuns;
end

function [AvgParOnes] = AverageParasiteOnes(StoreTotalOnesPar,NoRuns)
    SumofOnes = sum(StoreTotalOnesPar);
    AvgParOnes = SumofOnes/NoRuns;
end

function [AvgRelFitH,AvgRelFitP] = AverageRelativeFitHost(StoreRelativeFitP,StoreRelativeFitH,NoRuns)
    SumRelFit = sum(StoreRelativeFitH);
    AvgRelFitH = SumRelFit/NoRuns;
    SumRelFit = sum(StoreRelativeFitP);
    AvgRelFitP = SumRelFit/NoRuns;
end



function StandardDev(StoreTotalOnesHost,StoreTotalOnesPar,NoRuns, LengthGen, AvgHostOnes,AvgParOnes)
    sHost = std(StoreTotalOnesHost);
    sPar = std(StoreTotalOnesPar);
    figure(4)
    errorbar(LengthGen,AvgHostOnes,sHost,'m') 
    errorbar(LengthGen,AvgParOnes,sPar,'m')
end

function [HostFit,HostFitness] = SimpleHostFitness(HostPop,NoHost,HostFitness)
    for i = 1:NoHost
       HostFitness(i) = sum(HostPop(i,:));
    end
    
    [RankId,RankFit]= sort(HostFitness);
    
    %Setting Host Fitness based on number of ones, scaling from 0-1 
    
    for v = 1:NoHost
       HostFitA(RankFit(v)) = v/NoHost;
    end
    %Virulence function for host but lambda always 1
    for v = 1:NoHost
       HostFit(v) = ((2*HostFitA(v))/1) - ((HostFitA(v)^2)/(1^2));
    end
end

function [ParFit] = SimpleParFitness(ParPop,NoPar)
    for i = 1:NoPar
        ParFit(i) = sum(ParPop(i,:));
    end 
end

function [ParasiteFit] = ParasiteFitness(NoPar,ParFit,lambda)
    
    [rankId,rankFit] = sort(ParFit); 

    for p = 1:NoPar
        XScore(rankFit(p)) = p/NoPar;
    end
    %Virulence comes into play here. Parasite fitness with modified virulence
    for b = 1:NoPar
        ParasiteFit(b) = ((2*XScore(b))/lambda) - ((XScore(b)^2)/(lambda^2)); 
    end
end
