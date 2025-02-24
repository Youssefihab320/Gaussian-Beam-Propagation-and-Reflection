% Parameters
lambda = 3.2e-6;              % Wavelength in meters
w0 = 20e-6;                   % Beam waist in meters
k = 2 * pi / lambda;          % Wave number
z0 = pi * w0^2 / lambda;      % Rayleigh range

% Spatial domain
x = linspace(-5*w0, 5*w0, 1024);

y = x;
[X, Y] = meshgrid(x, y);
rho = sqrt(X.^2 + Y.^2);

% Initial Gaussian beam at z=0
U0 = exp(-rho.^2 / w0^2);

% Fourier Transform propagation
dk = 2 * pi / (x(end) - x(1));
kx = linspace(-dk*length(x)/2, dk*length(x)/2 - dk, length(x));
ky = kx;
[KX, KY] = meshgrid(kx, ky);
kz_squared = k^2 - KX.^2 - KY.^2;
kz = sqrt(max(0, kz_squared));         % Avoid imaginary kz

% Propagate beam to 0.5z0 and z0
Uz_0_5 = ifft2(ifftshift(fftshift(fft2(U0)) .* exp(-1j * kz * (0.5 * z0))));
Uz_1 = ifft2(ifftshift(fftshift(fft2(U0)) .* exp(-1j * kz * z0)));

% Reflect from mirror
mirror = exp(-1j * k * (X.^2 + Y.^2) / (-2 * z0));
U_mirror = Uz_1 .* mirror;

% Propagate after reflection
Uz_ref_0_5 = ifft2(ifftshift(fftshift(fft2(U_mirror)) .* exp(-1j * kz * (0.5 * z0))));
Uz_ref_1 = ifft2(ifftshift(fftshift(fft2(U_mirror)) .* exp(-1j * kz * z0)));
Uz_ref_2 = ifft2(ifftshift(fftshift(fft2(U_mirror)) .* exp(-1j * kz * (2 * z0))));

% Plot intensity at z = 0
figure;
imagesc(x * 1e6, y * 1e6, abs(U0).^2);
title('Intensity at z = 0');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% Plot intensity at z = 0.5z0
figure;
imagesc(x * 1e6, y * 1e6, abs(Uz_0_5).^2);
title('Intensity at z = 0.5z0');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% Plot intensity at z = z0
figure;
imagesc(x * 1e6, y * 1e6, abs(Uz_1).^2);
title('Intensity at z = z0');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% Plot intensity after reflection at z = 0.5z0
figure;
imagesc(x * 1e6, y * 1e6, abs(Uz_ref_0_5).^2);
title('Intensity after reflection at z = 0.5z0');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% Plot intensity after reflection at z = z0
figure;
imagesc(x * 1e6, y * 1e6, abs(Uz_ref_1).^2);
title('Intensity after reflection at z = z0');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% Plot intensity after reflection at z = 2z0
figure;
imagesc(x * 1e6, y * 1e6, abs(Uz_ref_2).^2);
title('Intensity after reflection at z = 2z0');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

figure;

% Plot intensity at z = 0
subplot(3, 2, 1);
plot(x * 1e6, abs(U0(round(end/2), :)).^2, 'b');
title('z = 0');
xlabel('x (\mum)');
ylabel('Intensity');
grid on;

% Plot intensity at z = 0.5z0
subplot(3, 2, 3);
plot(x * 1e6, abs(Uz_0_5(round(end/2), :)).^2, 'r');
title('z = 0.5z0');
xlabel('x (\mum)');
ylabel('Intensity');
grid on;

% Plot intensity at z = z0
subplot(3, 2, 5);
plot(x * 1e6, abs(Uz_1(round(end/2), :)).^2, 'g');
title('z = z0');
xlabel('x (\mum)');
ylabel('Intensity');
grid on;

% Plot intensity after reflection at z = 0.5z0
subplot(3, 2, 2);
plot(x * 1e6, abs(Uz_ref_0_5(round(end/2), :)).^2, 'm');
title('After reflection at z = 0.5z0');
xlabel('x (\mum)');
ylabel('Intensity');
grid on;

% Plot intensity after reflection at z = z0
subplot(3, 2, 4);
plot(x * 1e6, abs(Uz_ref_1(round(end/2), :)).^2, 'c');
title('After reflection at z = z0');
xlabel('x (\mum)');
ylabel('Intensity');
grid on;

% Plot intensity after reflection at z = 2z0
subplot(3, 2, 6);
plot(x * 1e6, abs(Uz_ref_2(round(end/2), :)).^2, 'k');
title('After reflection at z = 2z0');
xlabel('x (\mum)');
ylabel('Intensity');
grid on;

