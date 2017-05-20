clear all;
close all;

Ch = load('CHANNEL.mat');

% mx
figure;
subplot(2,1,1);
plot(Ch.f/1000,1e6*abs(Ch.Hxmx),'k');
hold on;
grid on;
plot(Ch.f/1000,1e6*abs(Ch.Hymx),'r');
plot(Ch.f/1000,1e6*abs(Ch.Hzmx),'b');
title('Amplitude of underground RX field from x directed magnetic dipole above ground');
xlabel('f[kHz]');
ylabel('|H|[uA/m]');
legend('Hx','Hy','Hz');

subplot(2,1,2);
plot(Ch.f/1000,phase(Ch.Hxmx)*180/pi,'k');
grid on;
title('Phase of underground RX field from x directed magnetic dipole above ground');
xlabel('f[kHz]');
ylabel('<H[deg]');
legend('Hx');

% my
figure;
subplot(2,1,1);
plot(Ch.f/1000,1e6*abs(Ch.Hxmy),'k');
hold on;
grid on;
plot(Ch.f/1000,1e6*abs(Ch.Hymy),'r');
plot(Ch.f/1000,1e6*abs(Ch.Hzmy),'b');
title('Amplitude of underground RX field from y directed magnetic dipole above ground');
xlabel('f[kHz]');
ylabel('|H|[uA/m]');
legend('Hx','Hy','Hz');

subplot(2,1,2);
plot(Ch.f/1000,phase(Ch.Hymy)*180/pi,'r');
grid on;
title('Phase of underground RX field from y directed magnetic dipole above ground');
xlabel('f[kHz]');
ylabel('<H[deg]');
legend('Hy');

% mz
figure;
subplot(2,1,1);
plot(Ch.f/1000,1e6*abs(Ch.Hxmz),'k');
hold on;
grid on;
plot(Ch.f/1000,1e6*abs(Ch.Hymz),'r');
plot(Ch.f/1000,1e6*abs(Ch.Hzmz),'b');
title('Amplitude of underground RX field from z directed magnetic dipole above ground');
xlabel('f[kHz]');
ylabel('|H|[uA/m]');
legend('Hx','Hy','Hz');

subplot(2,1,2);
plot(Ch.f/1000,phase(Ch.Hzmz)*180/pi,'b');
grid on;
title('Phase of underground RX field from z directed magnetic dipole above ground');
xlabel('f[kHz]');
ylabel('<H[deg]');
legend('Hz');