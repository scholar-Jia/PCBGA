%% ACO
function [luxian,T] = ACO(mission,dij,record,Tmax)
%% Get the coordinates of the city to be processed and make the initial map
a = record(1:end-1);
citys = mission(a,1:2);
%% Calculate the distance between each city
n = size(citys,1);
D = dij(a,a);
D(logical(eye(size(D)))) = eps;

%% Population initialization, setting of parameters
m = 75;                             
alpha = 1;                           
beta = 5;                            
vol = 0.2;                           
Q = 10;                              
Heu_F = 1./D;                        
Tau = ones(n,n);                     
Table = zeros(m,n);                  
Table(:,1) = 1;                         
iter_max = 100;                     
omega = 5;    
sudu = 0.02;   

Route_best = zeros(iter_max,n);      
Length_best = zeros(iter_max,1);     
Length_ave = zeros(iter_max,1);      

%% Iterate to find the best path
for iter = 1:iter_max
    for i=1:m
        Table(i,2) = randi([2,n],1,1);  
    end
   
    citys_index = 1:n;
    for i = 1:m
        for j = 3:n
            has_visited = Table(i,1:(j - 1)); 
            allow_index = ~ismember(citys_index,has_visited);
            allow = citys_index(allow_index);  
            P = allow;
            
            for k = 1:length(allow)
                P(k) = Tau(has_visited(end),allow(k))^alpha * (Heu_F(has_visited(end),allow(k)))^beta;
            end
            P = P/sum(P);
            Pc = cumsum(P); 
            target_index = find(Pc >= rand);
            target = allow(target_index(1));
            Table(i,j) = target;
        end
    end
    
    %%  First record the best route after the last iteration in the first position
    if iter>=2
        Table(1,:) = Route_best(iter-1,:);
    end
    
    %% Calculate the distance of each ant's path
    Length = zeros(m,1);
    for i = 1:m
        Route = Table(i,:); 
        for j = 1:(n - 1)
            Length(i) = Length(i) + D(Route(j),Route(j + 1)) ;
        end
        Length(i) = Length(i) + D(Route(n),Route(1)) ;
    end
    
    %% Calculate shortest path distance and average distance
    [min_Length,min_index] = min(Length);   
    Length_best(iter,:)= min_Length;           
    Route_best(iter,:) = Table(min_index,:);  
    Length_ave(iter) = mean(Length);       
    
    %% Update pheromone
    Delta_Tau = zeros(n,n);

    for i = 1:m
        for j = 1:(n - 1)
            Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
        end    
        Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/Length(i);
    end
    
    Tau = (1-vol)* Tau + Delta_Tau;
    Table = zeros(m,n);  
    Table(:,1) = 1;
end

%% The result of the command window shows
[~,index] = min(Length_best); 
Shortest_Route = Route_best(index,:);
Shortest_Route(end+1) = 1;

TACO=0;
Route = Shortest_Route; 
for j = 1:n
    if  j >= 2
        TACO= TACO + D(Route(j),Route(j + 1))/sudu +...
            turn_angle(citys(Route(j-1),:),citys(Route(j),:),citys(Route(j+1),:))/omega;
    else
        TACO = TACO+ D(Route(j),Route(j + 1))/sudu ;
    end
end

if TACO<Tmax
    luxian = record(Shortest_Route);
    T = TACO;
else
    luxian = record;
    T = Tmax;
end
end