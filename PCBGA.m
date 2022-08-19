clc;clear;close all;

%% Define
num_missionpoints=100;
num_uav = 5;
v = 0.02;
omega = 5;
Tmax = 5400;
alpha = 0.5; %alpha = (0, 1)

[real_missionpoints,anco] = Hotpoints(num_missionpoints);
missionpoints = real_missionpoints(:,1:2);
missionpoints(1,1:4) = [25,25,0,0];
weight_missionpoint = real_missionpoints(:,3);
weight_missionpoint(1,1) = 0;

theta = 0;
theta_base= 0;
camera_distance = 1;
T_Residual = 35;

P_missionpoints = missionpoints;
num_missionpoints = length(missionpoints);
sum_of_weights = sum(weight_missionpoint);
base = missionpoints(1,:);
%% Calculate the distance between mission points
distance = zeros(size(missionpoints,1));
for i = 1 : num_missionpoints
    for j = 1 : num_missionpoints
        distance(i,j) = sqrt( ( P_missionpoints(i,1) - P_missionpoints(j,1) ) ^ 2 + ( P_missionpoints(i,2) - P_missionpoints(j,2) ) ^ 2);
    end
end

%% Algorithm start
Sorting = 1:1:length(missionpoints);
remaining_mission_points = 1:1:length(missionpoints);
number = 0;
cnm = [];
xline = [];
Meancover2 = [];
k = 0;
everyluxian=cell(1,num_uav);
flight_path = 0;
cover_weight = zeros(1,num_uav);
Cover_temp = [];
secondary_advantage_Record = [];
T = zeros(num_uav,1);
Cover1_temp = [];
Cover2_temp = [];

while num_missionpoints >= 1 && k < num_uav
    k = k+1;
    flight_path = 0;
    flight_path(1) = 1;
    Tk_all = 0;
    T(k,1)= 0;
    T_judge =0;
    num = 1;
    theta_base = 0;
    Cover = [];Cover1 = [];Cover2 = [];
    Tcheck = zeros(26,1);
    
    for a = 1:(num_missionpoints-1-length(Cover))
        
        turning_Point_father(a,:) = missionpoints(Sorting(:,1),:);
        Sorting(:,1)=[];
        dis = zeros(1,(num_missionpoints-a-length(Cover)));
        cost = zeros(1,num_missionpoints-a-length(Cover));
        Timef = 0;
        theta = 0;
        
        if num_missionpoints - a -length(Cover) <= 0
            break
        end
        
        for b = 1:(num_missionpoints-a-length(Cover))
            dis(b) = distance(num, Sorting(b));
            if a >= 2
                theta(b) = acosd(dot([turning_Point_father(a,1)-turning_Point_father((a-1),1),turning_Point_father(a,2)-turning_Point_father((a-1),2)],[missionpoints(Sorting(b),1)-turning_Point_father(a,1),missionpoints(Sorting(b),2)-turning_Point_father(a,2)])...
                    /(norm([turning_Point_father(a,1)-turning_Point_father((a-1),1),turning_Point_father(a,2)-turning_Point_father((a-1),2)])*norm([missionpoints(Sorting(b),1)-turning_Point_father(a,1),missionpoints(Sorting(b),2)-turning_Point_father(a,2)])));
            else
                theta(b) = 0;
            end
            Timef(b) = dis(b)/v + theta(b)/omega;
            cost(b) = alpha * Timef(b) + (1-alpha)/weight_missionpoint(Sorting(b));
        end
        
        best_point = find( cost == min(cost));
        if length(best_point) >=2
            best_point = best_point(1);
        end
        
        theta_base = acosd(dot([missionpoints(Sorting(best_point),1)-turning_Point_father(a,1),missionpoints(Sorting(best_point),2)-turning_Point_father(a,2)],[base(1)-missionpoints(Sorting(best_point),1),base(2)-missionpoints(Sorting(best_point),2)])...
            /(norm([missionpoints(Sorting(best_point),1)-turning_Point_father(a,1),missionpoints(Sorting(best_point),2)-turning_Point_father(a,2)])*norm([base(1)-missionpoints(Sorting(best_point),1),base(2)-missionpoints(Sorting(best_point),2)])));
        T_judge= Tk_all + Timef(best_point)+distance(Sorting(best_point),1)/v + theta_base/omega;
        %% Seek secondary advantage point
        Waitjudge=length(cost)-1;
        T_judge_2 = T_judge;
        
        while T_judge_2>Tmax  && Waitjudge ~= 0
            Waitjudge = Waitjudge-1;
            delete = find(cost == min(cost));
            cost(delete) = inf;
            secondary_advantage = find(cost == min(cost));
            theta_base = turn_angle(turning_Point_father(a,:),missionpoints(Sorting(secondary_advantage),:),base);
            T_judge_2 = Tk_all + Timef(secondary_advantage)+distance(Sorting(secondary_advantage),1)/v + theta_base/omega;
        end
        
        if T_judge > Tmax && Waitjudge == 0
            break;
        elseif T_judge > Tmax &&  Waitjudge~=0
            secondary_advantage_Record = [secondary_advantage_Record,Sorting(secondary_advantage)];
            best_point = secondary_advantage;
        end
        
        temp = Sorting(best_point);
        father = P_missionpoints(num,:);
        son = P_missionpoints(Sorting(best_point),:);
        
        %% Points that can be covered by the camera at the same time
        best_point_temp = Sorting(best_point);
        Meancover1 = [];
        for i = 1:length(Sorting)
            if distance(best_point_temp,Sorting(i))<camera_distance && distance(best_point_temp,Sorting(i))>0
                Meancover1 = [Meancover1,Sorting(i)];
            end
        end
        Meancover1 = setdiff(Meancover1,Cover);
        if isempty(Meancover1) == 0
            [~,pos] = ismember(Meancover1,Sorting);
            Cover1 = [Cover1,Meancover1];
            Sorting(pos) = [];
            Timef(pos) = [];
            theta(pos) = [];
        end
        
        Tk_all = Tk_all + Timef(Sorting == best_point_temp);
        Tcheck(a,1)=Tk_all;
        Cover2 = [];
        Cover = [Cover1,Cover2];
        t = Sorting(1);
        Sorting(Sorting == best_point_temp) = t;
        Sorting(1) = best_point_temp;
        num = Sorting(1);
        flight_path(a+1) = Sorting(1);
    end
    
    %% End the path of UAV k and complete the return to base point
    number = length(flight_path);
    if  number >= 2
        done = [flight_path(2:number),Cover];
        remaining_mission_points = setdiff(remaining_mission_points,done);
        Sorting = remaining_mission_points;
    else
        break
    end
    
    cover_weight(k) = sum(weight_missionpoint(flight_path(2:number)));
    flight_path(number+1) = 1;
    father_of_father = P_missionpoints(flight_path(number-1),:);
    father = P_missionpoints(flight_path(number),:);
    son = base(1,:);
    theta_back = turn_angle(father_of_father,father,son);
    T(k)=Tk_all +  distance(flight_path(number),1)/v + theta_back/omega;
    [flight_path(1,:),T(k)] = ACO(missionpoints,distance,flight_path,T(k));
    
    everyluxian(1,k) = {flight_path};
    Cover_temp = [Cover_temp,Cover];
    Cover1_temp = [Cover1_temp,Cover1];
    Cover2_temp = [Cover2_temp,Cover2];
    num_missionpoints = length(remaining_mission_points);
end


%% Resulting parameters
T = T/60;
Tasktime = max(T);
Quan_T = max(T);
Quan_quanzhong = sum(cover_weight)+sum(weight_missionpoint(Cover_temp));
Zyita = Quan_quanzhong/sum_of_weights;

%% Draw Picture
hold on ;
suiji = 'wwwww';
xiantiao={'-','-','-','-','-'};
for i = 1:num_missionpoints
    scatter(P_missionpoints(i,1),P_missionpoints(i,2),20,[1,1,1],'filled');
end
ll=0;
for i = 1:num_uav
    huatu = everyluxian{1,i};
    for j = 1:length(huatu)-1
        line([P_missionpoints(huatu(j),1),P_missionpoints(huatu(j+1),1)],[P_missionpoints(huatu(j),2),P_missionpoints(huatu(j+1),2)],'Color',suiji(1,i),'linestyle',xiantiao{i},'linewidth',1.5);
    end
end
scatter(P_missionpoints(1,1),P_missionpoints(1,2),'sk','filled');
title('PCBGA','FontName','Times New Roman','fontsize',14);
xlabel('X / km','FontName','Times New Roman','fontsize',13);
ylabel('Y / km','FontName','Times New Roman','fontsize',13);
axis square;
grid on