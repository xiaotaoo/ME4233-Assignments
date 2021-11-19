%ME4233 Part 2 Assignment 2
%Ho Jun Yang  A0242277A
%Ong Jun Xian A0197294M
%Yeo Huai Zhe A0206224U
%Solving a lid-driven cavity flow problem with a rectangular domain

clear;
close All;

%Reynolds number is calculated by taking the average of the first three
%digits of our matriculation number
Re=(1/3)*(24+19+20);

%Setting up the grid points and spaces
%Number of grid points is 51 means number of grid intervals is 50
Nx=50;Ny=20;
Lx=3;Ly=1;
%Create the grids and dx dy.
x=linspace(0,Lx,Nx+1);dx=x(2)-x(1);
y=linspace(0,Ly,Ny+1);dy=y(2)-y(1);

[x2,y2]=meshgrid(x,y); %Obtain the mesh to plot the velocity fields

% specify the value of dt
%dt = 0.001; %Time step of 0.001s until the time reaches 2s              
dt = 0.02; %Max explicit dt at 0.0154918
%dt = 0.02; %dt at value > Max explicit dt.
t_end = 2;
NMax=t_end/dt;   % maximum no. of time steps

% assemble A matrix to solve the Poisson equation Au=b
[A]=assembleA(Nx,Ny,dx,dy);
% initial condition for vorticity field
vort = zeros(Nx-1,Ny-1);
% initial time
t=0;
% specify method
method = 2; %1 -> Explicit Euler, 2 -> Implicit Euler
vel_arr=zeros(NMax,1);
% time evoluion 
for i=1:NMax
    streamfunc = Poisson_Solver(vort,A,Nx,Ny);  % solve Poisson equation for streamfunction
    vort = adv_vort(method,streamfunc,vort,Nx,Ny,dx,dy,dt,Re,t);  %advance the vorticity equation using explict/implict discretization
    disp(['Finish timestep' num2str(i)])
    t=t+dt;

    % post-processing   
    % plot the streamfunction with boundary conditions
    streamfunc_plot=[0             zeros(1,Ny-1)   0
                     zeros(Nx-1,1) streamfunc         zeros(Nx-1,1)
                     0             zeros(1,Ny-1)   0 ];   
    % for plotting (vector plot)
    fig1 = figure(1);
    movegui(fig1, [350,1080/2 + 10]);
    [u,v]=construct(streamfunc_plot,Nx,Ny,dx,dy,t);
    quiver(x2,y2,u',v')
    xlim([0 3]);ylim([0 1.01]);
    title(['At time step = ' num2str(i)])
    xlabel('x');ylabel('y')
    pause(0.2)
    
    % Streamfunction Plot
    fig2 = figure(2);
    movegui(fig2, [1920/2 + 40,1080/2 + 10]);
    contourf(x,y,streamfunc_plot')
    title('Streamfunction')
    xlabel('x');ylabel('y')
    
    % Part (c)
    vel_arr(i) = u(26,11); %Getting value of u at x=1.5 and y=0.5
end    

% spectral analysis
vel_time=linspace(0,dt*NMax,NMax);
fig3 = figure(3);
movegui(fig3, [350,10]);
plot(vel_time,vel_arr)
xlim([0 2]);ylim([-0.06 0.06]);
title('U at C(1.5,0.5) over time')

vel_Fourier = fft(vel_arr);
vel_Magnitude = abs(vel_Fourier);
vel_FreqN = length(vel_Fourier);
vel_Fbins=((0: 1/vel_FreqN: 1-1/vel_FreqN)*(1/dt)).';

fig4 = figure(4);
movegui(fig4, [1920/2 + 40,10]);
plot(vel_Fbins, vel_Magnitude)
xlim([0 10]);
title('Fourier Analysis of U at C')
xlabel('Frequency(Hz)')