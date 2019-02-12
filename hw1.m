clear all; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks); %frequency
%matrix to store average frequency signal, which will contain the
%central signal since the noise get canceled out
Uave = zeros(n,n,n);

for j=1:20
    %Average out for each time step
    Un(:,:,:)=reshape(Undata(j,:),n,n,n); %this is the data we received each time
    Ut(:,:,:) = fftn(Un); %convert to frequency space
    Uave = Uave + Ut;
end
Uave = fftshift(abs(Uave)) / 20; %take average after each time step
Umax = max(max(max(abs(Uave)))); %spectral signature
close all, isosurface(Kx,Ky,Kz, Uave / Umax,0.4);

axis([-20 20 -20 20 -20 20]), grid on, drawnow
title('Fig 1: Average spectral profiles for the 20 realizations')
xlabel('Frequency x direction')
ylabel('Frequency y direction')
zlabel('Frequency z direction')
[freq_max_x, freq_max_y, freq_max_z] = ind2sub([n,n,n], find(Uave == Umax));


width = 0.2;
filter = exp(-width * (Kx - Kx(freq_max_x,freq_max_y,freq_max_z)).^2 - width * (Ky -Ky(freq_max_x,freq_max_y,freq_max_z)).^2 - width * (Kz -Kz(freq_max_x,freq_max_y,freq_max_z)).^2);
filter = fftshift(filter);
position_over_time = zeros(20,3);
%%
for k = 1:20
    Un(:,:,:)=reshape(Undata(k,:),n,n,n); %this is the data we received each time
    Ut(:,:,:) = fftn(Un);
    Utf(:,:,:) = Ut(:,:,:) .* filter;
    Utif(:,:,:) = ifftn(Utf);
    isosurface(X,Y,Z,(abs(Utif)),0.4)
    [index_x, index_y, index_z] = ind2sub([n,n,n], find(abs(Utif) == max(max(max(abs(Utif))))));
    position_x = x(index_x);
    position_y = y(index_y);
    position_z = z(index_z);
    position_over_time(k, 1) = position_x;
    position_over_time(k, 2) = position_y;
    position_over_time(k, 3) = position_z;
      
    title('Fig 2: Spatial position of object in animal body for the 20 realizations')
    xlabel('x ')
    ylabel('y ')
    zlabel('z')
     axis([-20 20 -20 20 -20 20]), grid on, 
     drawnow
     pause(.1)
end

%%
close all, isosurface(X,Y,Z,abs(Un),0.4)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
title('Figure 1: Noisy ultrasound data of small area in intestines where the marble is expected to be at 20th realization')
xlabel('x ')
ylabel('y ')
zlabel('z')
%%
figure
plot3(position_over_time(:,1), position_over_time(:, 2), position_over_time(:,3));
xlabel('x ')
ylabel('y ')
zlabel('z ')


