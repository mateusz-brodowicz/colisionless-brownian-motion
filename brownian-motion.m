close all
clear all
clc

%plotting the box
grid on
hold on
x_lim = [-500, -500, 500, 500, -500];
y_lim = [500, -500, -500, 500, 500];
xlim([-500,500]);
ylim([-500,500]);
plot(x_lim,y_lim);

%creating N particles in random places within R range from the centre
%defining N and the radius of the boundary
N = 20;
R_lim = 250;

%creating points in polar coordinate system
polar_ang = rand(N,1)*360;
polar_r = rand(N,1)*R_lim;

%converting the polar coordiante system to carthesian coordinates and
%plotting
[x, y] = pol2cart(polar_ang,polar_r);
plot_pts = scatter(x,y,'filled','blue');

%create trace arrays
trace_len = 40;
trace_x = x;
trace_y = y;
for i=1:trace_len-1
    trace_x = [trace_x,x];
    trace_y = [trace_y,y];
end

%create initial trace plot
for j=1:N
    trace_plot(j) = plot(trace_x(j,:),trace_y(j,:),"--k");
end

%defining the maximum and minimum possible velocity
v_max = 80;
v_min = 30;

%multiply random samples with difference of vmax and vmin then add/subtract
%v_min(depending on the sign of a value in a vector) to assure each value is in range+/-(vmin,vmax)
v_x = -(v_max-v_min) + 2*(v_max-v_min)*rand(N,1);
v_x = v_x + sign(v_x)*v_min;

v_y = -(v_max-v_min) + 2*(v_max-v_min)*rand(N,1);
v_y = v_y + sign(v_y)*v_min;

%defining starting time t and delta t dt
t = 0;
dt = 0.1;

%defining gravity
g = -9.8;

%defining hit velocity loss (%)
loss = 0.1;

%define indexing variable i which will help with adding values to vector
%used to plot energy graphs
i = 1;

%run the function
while isgraphics(plot_pts)
    
    %move the particles
    v_y = v_y + g*dt;
    x = x+v_x*dt;
    y = y+v_y*dt;
    
    %update first column in trace matrix and shift columns to the right
    trace_x(:,1)=x;
    trace_y(:,1)=y;
    trace_x(:,2:end)=trace_x(:,1:end-1);
    trace_y(:,2:end)=trace_y(:,1:end-1);
    
    %check if particle crossed boundaries on x or y
    crossed_y_lim=find(y<=-500 | y>=500);
    crossed_x_lim=find(x<=-500 | x>=500);
      
    %if particle bounced off the ceiling/floor, find nearest boundary(-500 or 500) and
    %set y value to it
    y(crossed_y_lim)=500*sign(y(find(abs(y)>500)));
    
    %apply energy loss if bounced off the ceiling/floor
    v_y(crossed_y_lim)=v_y(crossed_y_lim)*(1-loss);
    v_x(crossed_y_lim)=v_x(crossed_y_lim)*(1-loss);
    
    % flip vertical direction of velocity
    v_y(crossed_y_lim)=-v_y(crossed_y_lim);
     
    %if particle bounced off a wall, find nearest boundary(-500 or 500) and
    %set x value to it
    x(crossed_x_lim)=500*sign(x(find(abs(x)>500)));
    
    %apply energy loss if bounced off the wall
    v_y(crossed_x_lim)=v_y(crossed_x_lim)*(1-loss);
    v_x(crossed_x_lim)=v_x(crossed_x_lim)*(1-loss);
    
    % flip horizontal direction of velocity
    v_x(crossed_x_lim)=-v_x(crossed_x_lim);
    
    %plot updated locations of particles
    plot_pts.XData=x;
    plot_pts.YData=y;
    
    for j=1:N
        trace_plot(j).XData=trace_x(j,:);
        trace_plot(j).YData=trace_y(j,:);
    end
    
    %add time stamp and energy stamp to the vectors, mass set as 1 for
    %simplicty
    t=t+dt;
    e_k_stamps(i) = sum([1/2*(v_x.^2+v_y.^2).^0.5]); 
    t_stamps(i)=t;
    
    %saving each 20th-snapshot (iteration) within first 500 iterations
    %if (mod(i,10)==0) && (i<=600) % if iter is multiple of 10 and is not greater than 500
        %filename=['gravityloss' num2str(i) '.png'];
        %saveas(gcf,filename);
    %end
    i=i+1;
    drawnow;
end

plot(t_stamps,e_k_stamps)
title('Kinetic energy in time')
xlabel('time(s)')
ylabel('energy(j)')