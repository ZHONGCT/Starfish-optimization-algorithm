function [xposbest,fvalbest,Curve] = SFOA(Npop,Max_it,lb,ub,nD,fobj)
% Starfish optimization algorithm (SFOA)
% Created by Dr. Changting Zhong (Email: zhongct@hainanu.edu.cn)
% Paper: Changting Zhong, Gang Li, Zeng Meng, Haijiang Li, Ali Riza Yildiz, Seyedali Mirjalili. 
% Starfish Optimization Algorithm (SFOA): A bio-inspired metaheuristic algorithm for global optimization compared with 100 optimizers
% Neural Computing and Applications, 2025, 37: 3641-3683.

%% Initialization
GP=0.5;     % parameter 

if size(ub,2) == 1
    lb = lb*ones(1,nD); ub = ub*ones(1,nD);
end

fvalbest = inf;
Curve = zeros(1,Max_it);
Xpos = rand(Npop,nD).*(ub-lb)+lb;
for i = 1:Npop
    Fitness(i) = feval(fobj,Xpos(i,:));
end
[fvalbest,order] = min(Fitness);    % global best fitness
xposbest = Xpos(order,:);        % global best position

newX = zeros(Npop,nD);
%% Evolution
T = 1;
while T <= Max_it
    theta = pi/2*T./Max_it;
    tEO = (Max_it-T)/Max_it*cos(theta);
    if rand < GP %  exploration of starfish
        for i = 1:Npop
            if nD > 5
                % for nD is larger than 5
                jp1 = randperm(nD,5);
                for j = 1:5
                    pm = (2*rand-1)*pi;
                    if rand < GP
                        newX(i,jp1(j)) = Xpos(i,jp1(j)) + pm*(xposbest(jp1(j))-Xpos(i,jp1(j)))*cos(theta);
                    else
                        newX(i,jp1(j)) = Xpos(i,jp1(j)) - pm*(xposbest(jp1(j))-Xpos(i,jp1(j)))*sin(theta);
                    end
                    if newX(i,jp1(j))>ub(jp1(j)) || newX(i,jp1(j))<lb(jp1(j))
                        newX(i,jp1(j)) = Xpos(i,jp1(j));
                    end
                end
            else
                % for nD is not larger than 5
                jp2 = ceil(nD*rand);
                im = randperm(Npop);
                rand1 = 2*rand-1;
                rand2 = 2*rand-1;
                newX(i,jp2) = tEO*Xpos(i,jp2) + rand1*(Xpos(im(1),jp2)-Xpos(i,jp2))+rand2*(Xpos(im(2),jp2)-Xpos(i,jp2));
                if newX(i,jp2)>ub(jp2) || newX(i,jp2)<lb(jp2)
                    newX(i,jp2) = Xpos(i,jp2);
                end  
            end
            newX(i,:) = max(min(newX(i,:),ub),lb);  % boundary check
        end
    else % exploitation of starfish
        df = randperm(Npop,5);
        dm(1,:) = xposbest - Xpos(df(1),:);
        dm(2,:) = xposbest - Xpos(df(2),:);
        dm(3,:) = xposbest - Xpos(df(3),:);
        dm(4,:) = xposbest - Xpos(df(4),:);
        dm(5,:) = xposbest - Xpos(df(5),:);  % five arms of starfish
        for i = 1:Npop
            r1 = rand; r2 = rand;
            kp = randperm(length(df),2);
            newX(i,:) = Xpos(i,:) + r1*dm(kp(1),:) + r2*dm(kp(2),:);   % exploitation
            if i == Npop
                newX(i,:) = exp(-T*Npop/Max_it).*Xpos(i,:);  % regeneration of starfish
            end
            newX(i,:) = max(min(newX(i,:),ub),lb);  % boundary check
        end
    end
    
    % Fitness evaluation
    for i = 1:Npop
        newFit = feval(fobj,newX(i,:));
        if newFit < Fitness(i)
            Fitness(i) = newFit;
            Xpos(i,:) = newX(i,:);
            if newFit < fvalbest
                fvalbest = Fitness(i);
                xposbest = Xpos(i,:);
            end
        end
    end
    
    Curve(T) = fvalbest;
    T = T+1;
end