function [real_taskpoints,anco] = Hotpoints(inputArg1)

num_taskpoints = inputArg1;
num_point = 100;
taskpoint = 50*rand(num_point,4);
taskpoint(1,:) = [25,25,0,0];
taskpoint(:,3:4) = 0;
taskpoint = round(taskpoint*1e1)/1e1;
weight = 10*rand(num_point,1)';
weight(1,1) = 0;


anco = taskpoint;
anco(:,3)=weight';
anco = sortrows(anco,1);
anco_x=anco(:,1);
anco_y=anco(:,2);
anco_z=anco(:,3);

X=anco_x;
Y=anco_y';
Z=anco_z;
[x,y]=meshgrid(0:1:50,0:1:50);
z=griddata(X,Y,Z,x,y,'v4');
[C,h] = contourf(x,y,z,6);
% close

nn = zeros;
i = 1;
ii = 0;
while i < size(C,2)
    ii = ii+1;
    nn(ii,1) = C(1,i);
    nn(ii,2) = C(2,i);
    nn(ii,3) = i;
    i = i + C(2,i)+1;
end

step_1 = nn;
step_1((step_1(:,1)<0),:) = [];
step_2 = unique(step_1(:,1));
point_and_weight = [];
segment_num = zeros(size(step_2,1)+1,1);
for i = 1: size(step_2,1)+1
    if i == 1
        [row,col] = find(0<z & z<step_2(i,1));
        v_num = z((0<z & z<step_2(i,1)));
    elseif i == size(step_2,1)+1
        [row,col] = find(z>step_2(i-1,1));
        v_num = z(z>step_2(i-1,1));
    else
        [row,col] = find(step_2(i-1,1)<z & z<step_2(i,1));
        v_num = z(step_2(i-1,1)<z & z<step_2(i,1));
    end
    M=[y(col,1),x(1,row)'];
    segment_num(i,1) = size(v_num,1);
    point_and_weight = [point_and_weight;M,v_num];
end

segment_num(:,2) = round(segment_num(:,1) * num_taskpoints / sum(segment_num(:,1)));
real_taskpoints = [];
ii = 1;
for i = 1 :size(segment_num,1)
    num = segment_num(i,1);
    temp = point_and_weight(ii:ii+num-1,:);
    step_3 = temp(randperm(num,segment_num(i,2)),:);
    
    real_taskpoints = [real_taskpoints;step_3];
    temp = zeros;
    ii = ii + num;
end
end

