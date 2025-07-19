clear, close all
format longE
addpath 'tools/'

% ------------------------------------------------------------ %

% declaration of global variables:
global dy nu Dx Dy Kx Ky K2 Lx Ly Nx Ny Nxy X Y

% ------------------------------------------------------------ %

%%% Physical parameters:
nu = 1E-4;	% 1/Re or viscosity

% ------------------------------------------------------------ %

%%% Domain definition:
Lx = pi;    % Domain half-length in x-direction
Ly = pi;    % Domain half-length in y-direction

%%% Numerical parameters:
Nx = 1024;  % number of Fourier modes in discrete solution x-dir
Ny = 1024;	% number of Fourier modes in discrete solution y-dir
Nxy = Nx*Ny;

dx = 2*Lx/Nx;   		% distance between two physical points
x = (1-Nx/2:Nx/2)'*dx;% physical space discretization
%x=(0:Nx-1)'*(2*Lx/Nx);  		
%y=(0:Ny-1)'*(2*Ly/Ny); 

dy = 2*Ly/Ny;   		% distance between two physical points
y = (1-Ny/2:Ny/2)'*dy;  % physical space discretization

[X,Y] = meshgrid(x,y);	% 2D composed grid

% ------------------------------------------------------------ %

% vectors of wavenumbers in the transformed space:
 kx = [0:Nx/2 1-Nx/2:-1]'*pi/Lx;
 ky = [0:Ny/2 1-Ny/2:-1]'*pi/Ly;

%kx = [Nx/2:Nx 1-Nx:(-Nx/2-1)]'*pi/Lx;
%ky = [Ny/2:Ny 1-Ny:(-Ny/2-1)]'*pi/Ly;

% antialising treatment
jx = (Nx/4+2:Nx/4*3);  % the frequencies we sacrify
kx(jx) = 0;

jy = (Ny/4+2:Ny/4*3);  % the frequencies we sacrify
ky(jy) = 0;

% ------------------------------------------------------------ %

%%% Some operators arising in NS equations:
[Kx, Ky] = meshgrid(kx,ky);
K2 = Kx.^2 + Ky.^2;     % to compute the Laplace operator

K2inv = zeros(size(K2));
K2inv(K2 ~= 0) = 1./K2(K2 ~= 0);

Dx = 1i*Ky.*K2inv;		% u velocity component reconstruction from the vorticity
Dy = 1i*Kx.*K2inv;		% v velocity component reconstruction from the vorticity

% ------------------------------------------------------------ %

fftw('planner', 'hybrid');

%%% Set random number generator (for the initial condition)
s = RandStream('dsfmt19937','Seed',2);
RandStream.setGlobalStream(s);

% ------------------------------------------------------------ %

%%% Time-stepping parameters:
t = 0.0;           	% the discrete time variable
Tf = 50;          	% final simulation time
ds = 0.1;			% write time of the results

ops = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'OutputFcn', @odetpbar);

% ------------------------------------------------------------ %

Om_hat = zeros(Nx, Ny);

	%We take a random distribution of the vorticity:
	Om_hat(1,5) = randn(1) + 1i*randn(1);
	Om_hat(2,2) =randn(1) + 1i*randn(1);
	Om_hat(4,1) = randn(1) + 1i*randn(1);

	
Om = real(ifft2(Om_hat));
Om = Om/max(max(Om));	% normalize to O(1)
Om = InitCond();
% Om = zeros(Nx, Ny);
 FigHandle = figure(1);
 set(FigHandle, 'Position', [100, 100, 600, 550]);
% 
 frame = 0;
 Plot(Om, t, frame);	% we plot the initial condition
 saveas(figure(1),'first','jpg')
% 
% % ------------------------------------------------------------ %
% 


Om_vec = reshape(Om, Nxy, 1);
%vid=VideoWriter('video','MPEG-4');
%open(vid);
while (t < Tf) % main loop in time
	[~, v] = ode113(@RHS, [0:ds:ds], Om_vec, ops);

	Om_vec = v(end,:);
	Om = reshape(Om_vec, Nx, Ny);

	t = t + ds; frame = frame + 1;
    FigHandle = figure(1);
    set(FigHandle, 'Position', [100, 100, 600, 550]);
	Plot(Om, t, frame); % and we plot the solution
    %F=getframe(figure(2));
    %writeVideo(vid, F);
end % while (t)
%close(vid);


% calculating streamfunction
% solving poisson equation

% in x-direction
 for p=1:Nx
     Q_x=fft(-Om);
 end
 Q_x=Q_x/Nx;

 for i =1:Nx
    Q_x_shifted(1:Nx,i)=fftshift(Q_x(1:Nx,i));
 end
%in y-direction
 for m=1:Ny
         Q_y(m,:)=fft(Q_x_shifted(m,:));
 end
 Q_xy=Q_y/Ny;

 for i =1:Ny
   Q_xy_shifted(i,1:Ny)=fftshift(Q_xy(i,1:Ny));
 end

% phi_cap(m,n)

 for  m=1:Nx
     for n=1:Ny
         k(m)=(2.*pi.*(m-(Nx/2+1)))./(2*Lx);
         l(n)=(2.*pi.*(n-(Ny/2+1)))./(2*Ly);
         phi_cap(m,n)=-Q_xy_shifted(m,n)./((k(m).^2)+(l(n).^2));
     end
 end


 phi_cap(Nx/2+1,Ny/2+1)=0;

 
 for i =1:Ny
     phi_cap(i,1:Ny)=ifftshift(phi_cap(i,1:Ny));
 end


 %inverse 2D FFT for phi(j,p)

% 2D inverse FFT

 %inverse in y dir:
for m=1:Ny
          phi_y(m,1:Ny)=(ifft(phi_cap(m,1:Ny).*Ny));
end


for i= 1:Ny
      phi_y(1:Ny,i)=fftshift(phi_y(1:Ny,i));
  end

% in x-direction
 for p=1:Ny
     phi_x=real((ifft(phi_y.*Ny)));
 end

 % plotting the streamfunction
figure(2)
[X,Y] = meshgrid(x,y);   % 2D mesh grid
surf(X, Y, phi_x);
view([0 90]);
shading flat;
colormap("jet")
cc = colorbar;
title ('Stream Function Distribution at t=50');
xlabel('x');
xlim([-Lx Lx]); 
ylabel('y');
ylim([-Ly+dy Ly]);
set(gca,'XTickLabel',[0 1 2 3 4 5 6])
set(gca,'YTickLabel',[0 1 2 3 4 5 6])
xlabel(cc, '\psi(x,y,t)');

[u,v] = gradient(phi_x) ;
figure(3)
surf(X, Y, u);
view([0 90]);
shading flat;
colormap("jet")
cc = colorbar;
title ('u velocity Distribution at t=50');
xlabel('x');
xlim([-Lx Lx]); 
ylabel('y');
ylim([-Ly+dy Ly]);
set(gca,'XTickLabel',[0 1 2 3 4 5 6])
set(gca,'YTickLabel',[0 1 2 3 4 5 6])
xlabel(cc, 'u(x,y,t)');

figure(4)
surf(X, Y, v);
view([0 90]);
shading flat;
colormap("jet")
cc = colorbar;
title ('v velocity Distribution at t=50');
xlabel('x');
xlim([-Lx Lx]); 
ylabel('y');
ylim([-Ly+dy Ly]);
set(gca,'XTickLabel',[0 1 2 3 4 5 6])
set(gca,'YTickLabel',[0 1 2 3 4 5 6])
xlabel(cc, 'v(x,y,t)');