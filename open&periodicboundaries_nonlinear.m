clear all; close all;

%% initialise variables - NON-LINEAR SYSTEM, OPEN/Periodic BOUNDARIES
U = 5.0;                    % background flow (m/s)
H = 10.0;                   % average height of fluid layer (m)
g = 9.81;                   % gravitational acceleration (m/s^2)
L = 100.0;                  % length of channel (m)
h0 = 0.5;                   % amplitude of initial distortion (m)
CFL = 0.4;                  % CFL criterion
dx = 0.25;                  % spatial distance (m)
dt = 0.001;                 % time step
ttot = 15.0;                % total time (s)
steps = round(ttot./dt);    % number of time steps

%% define initial distortion and speed at t = 0 s
for m = 1:4*L+3
    En=0;
    x(m) = (m-2) * dx;
    if abs(x(m)./L-0.5) <= 1.0./9.0
        hprime(m,1) = h0.*9.^4.*(((x(m)./L-0.5)).^2-(1.0./9.0).^2).^2;
    else
        hprime(m,1) = 0;
    end
    uprime(m,1) = 0;
    Enonlin=((H+hprime(m,1))*(((uprime(m,1))^2)/2))+(g*((hprime(m,1))^2)/2); %Energy (nonlinear)
    En=En+Enonlin;
end
E(1)=En;
c=sqrt(g*H);

%% first integration step using Euler's method
for n=2:2
    En=0;
    for m = 3:4*L+1
        uprime(m,n)=uprime(m,n-1)-(dt./dx)*((U+(uprime(m,n-1)))*0.5.*(uprime(m+1,n-1)-uprime(m-1,n-1))+g.*(hprime(m,n-1)-hprime(m-1,n-1))); 
        hprime(m,n)=hprime(m,n-1)-(dt./dx)*((U+(uprime(m,n-1)))*0.5.*(hprime(m+1,n-1)-hprime(m-1,n-1))+(H+(hprime(m,n-1)))*(uprime(m+1,n-1)-uprime(m,n-1)));
        Enonlin=((H+hprime(m,n))*(((uprime(m,n))^2)/2))+(g*((hprime(m,n))^2)/2); %Energy (nonlinear)
        En=En+Enonlin;
    end
    Enonlin=((H+hprime(2,n))*(((uprime(2,n))^2)/2))+(g*((hprime(2,n))^2)/2); %Energy (nonlinear)
    En=En+Enonlin;
    Enonlin=((H+hprime(402,n))*(((uprime(402,n))^2)/2))+(g*((hprime(402,n))^2)/2); %Energy (nonlinear)
    En=En+Enonlin;
    %%open boundaries
    %uprime(2,n)=uprime(2,n-1)+(c-U)*dt*((uprime(3,n-1)-uprime(2,n-1))/dx);
    %hprime(2,n)=hprime(2,n-1)+(c-U)*dt*((hprime(3,n-1)-hprime(2,n-1))/dx);
    %uprime(402,n)=uprime(402,n-1)-(U+c)*dt*((uprime(402,n-1)-uprime(401,n-1))/dx);
    %hprime(402,n)=hprime(402,n-1)-(U+c)*dt*((hprime(402,n-1)-hprime(401,n-1))/dx);
    E(n)=En;
end
%periodic boundaries
uprime(1,n) = uprime(402,n);
uprime(403,n) = uprime(2,n);
hprime(1,n) = hprime(402,n);
hprime(403,n) = hprime(2,n);

%% further integration steps with leapfrog method
for n=3:steps
    En=0;
    for m = 3:4*L+1
        uprime(m,n)=uprime(m,n-2)-(2.0.*dt./dx)*((U+(uprime(m,n-1))).*0.5.*(uprime(m+1,n-1)-uprime(m-1,n-1))+g.*(hprime(m,n-1)-hprime(m-1,n-1))); 
        hprime(m,n)=hprime(m,n-2)-(2.0.*dt./dx)*((U+(uprime(m,n-1))).*0.5.*(hprime(m+1,n-1)-hprime(m-1,n-1))+(H+(hprime(m,n-1))).*(uprime(m+1,n-1)-uprime(m,n-1)));
        Enonlin=((H+hprime(m,n))*(((uprime(m,n))^2)/2))+(g*((hprime(m,n))^2)/2); %Energy (nonlinear)
        En=En+Enonlin;
    end
  
    %%open boundaries
    %uprime(2,n)=uprime(2,n-1)+(c-U)*dt*((uprime(3,n-1)-uprime(2,n-1))/dx);
    %hprime(2,n)=hprime(2,n-1)+(c-U)*dt*((hprime(3,n-1)-hprime(2,n-1))/dx);
    %uprime(402,n)=uprime(402,n-1)-(U+c)*dt*((uprime(402,n-1)-uprime(401,n-1))/dx);
    %hprime(402,n)=hprime(402,n-1)-(U+c)*dt*((hprime(402,n-1)-hprime(401,n-1))/dx);  
    Enonlin=((H+hprime(2,n))*(((uprime(2,n))^2)/2))+(g*((hprime(2,n))^2)/2); %Energy (nonlinear)
    En=En+Enonlin;
    Enonlin=((H+hprime(402,n))*(((uprime(402,n))^2)/2))+(g*((hprime(402,n))^2)/2); %Energy (nonlinear)
    En=En+Enonlin;
    %periodic boundaries
    uprime(1,n) = uprime(402,n);
    uprime(403,n) = uprime(2,n);
    hprime(1,n) = hprime(402,n);
    hprime(403,n) = hprime(2,n);
    E(n)=En;
end

%% plot the results depending on time step
g = figure; 
for n = 1:50:steps
    time = n*dt;
subplot(3,1,1)       
    plot(x(:),uprime(:,n))
    axis([0 100 -0.5 0.5])
    xlabel('Position in region (m)')
    ylabel('u^\prime (m/s)')
    title(['Velocity fluctuation, t = ' num2str(time) ' s'])

subplot(3,1,2)       
    plot(x(:),hprime(:,n))
    axis([0 100 -0.5 0.5])
    xlabel('Position in region (m)')
    ylabel('h^\prime (m)')
    title(['Height fluctuation, t = ' num2str(time) ' s'])

subplot(3,1,3)
    plot(E(1:n))
    axis([0 15000 0 70])
    xlabel('Time (s)')
    ylabel('Energy (J)')
    title(['Energy of the non-linear system, t = ' num2str(time) ' s'])
pause(0.0001)   
end
