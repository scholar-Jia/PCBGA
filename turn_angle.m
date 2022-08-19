%½Ç¶È
function s = turn_angle(fu,now,future)
x = fu;
y = now;
z = future;
s = acosd(dot([y(1,1)-x(1,1),y(1,2)-x(1,2)],[z(1,1)-y(1,1),z(1,2)-y(1,2)])/(norm([y(1,1)-x(1,1),y(1,2)-x(1,2)])*norm([z(1,1)-y(1,1),z(1,2)-y(1,2)])));
end