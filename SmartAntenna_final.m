clear all
clc
tic
%% plot array pattern
f=1*10^9;
c=physconst('LightSpeed');
lamda=c/f;
d=lamda/2;
L=17;
r=[-4*d -3*d -2*d -d d 2*d 3*d 4*d];
R=[r zeros(1,8) 0;zeros(1,8) r 0];
theta=linspace(pi,-pi,1000);
v=[cos(theta);sin(theta)];
tao=[];
s=[];
for l=1:L
    tao(l,:)=R(:,l)'*v/c;
    s(l,:)=exp(1j*2*pi*f*tao(l,:));
end
polarplot(theta,abs(s'*ones(L,1)),'m')
figure(1)
title('Array Pattern of Cross Array')

%% setting noises with SNR
SNR=10;%dB
SNR_SONI1 = 30*rand;
SNR_SONI2 = 30*rand;
SNR_SONI3 = 30*rand;
N0=1;
Eb=N0*10^(SNR/10);
Eb_SONI1=N0*10^(SNR_SONI1/10);
Eb_SONI2=N0*10^(SNR_SONI2/10);
Eb_SONI3=N0*10^(SNR_SONI3/10);

%% input an image for TX data
imageData = imread('DSC_8468.tif');
size=[1000 1500]/2;
imageData = imresize(imageData,size); 
% transform image into binary as a vector
imageData_bit = de2bi(imageData);
imageData_bit_col=double(  reshape(imageData_bit,[size(1)*size(2)*8 1])  );


%% BPSK mod.& AWGN

bpsk_signal= sqrt(Eb)*(2*imageData_bit_col-1);

binary_stream = [zeros( length(bpsk_signal)/2 ,1)  ; ones(length(bpsk_signal)/2,1)];
SONI1 = binary_stream(randperm(length(bpsk_signal)));
SONI2 = binary_stream(randperm(length(bpsk_signal)));
SONI3 = binary_stream(randperm(length(bpsk_signal)));

bpsk_SONI1 = sqrt(Eb_SONI1)*(2*SONI1-1);
bpsk_SONI2 = sqrt(Eb_SONI2)*(2*SONI2-1);
bpsk_SONI3 = sqrt(Eb_SONI3)*(2*SONI3-1);

%% steering vector of SOI & SONI
theta_SOI = rand*2*pi;
theta_SONI1 = rand*2*pi;
theta_SONI2 =  rand*2*pi;
theta_SONI3 =  rand*2*pi;
s_SOI=[];
s1=[];
s2=[];
s3=[];

v_SOI=[cos(theta_SOI);sin(theta_SOI)];
v1=[cos(theta_SONI1);sin(theta_SONI1)];
v2=[cos(theta_SONI2);sin(theta_SONI2)];
v3=[cos(theta_SONI3);sin(theta_SONI3)];
tao=[];
for l=1:L
    tao(l,:)=R(:,l)'*v1/c;
    s1(l,:)=exp(1j*2*pi*f*tao(l,:));
end
for l=1:L
    tao(l,:)=R(:,l)'*v2/c;
    s2(l,:)=exp(1j*2*pi*f*tao(l,:));
end
for l=1:L
    tao(l,:)=R(:,l)'*v3/c;
    s3(l,:)=exp(1j*2*pi*f*tao(l,:));
end
for l=1:L
    tao(l,:)=R(:,l)'*v_SOI/c;
    s_SOI(l,:)=exp(1j*2*pi*f*tao(l,:));
end

%% LMS Algorithm
weight=zeros(L,1);
set_0f_SNR=[Eb Eb_SONI1 Eb_SONI2 Eb_SONI3];
mu=0.005%2/L/max(set_0f_SNR);
for bitnumber=1:length(bpsk_signal)
    noise = N0*randn(length(bpsk_signal),1);
    desire_signal = bpsk_signal(bitnumber);
    u = desire_signal*s_SOI + bpsk_SONI1(bitnumber)*s1 + bpsk_SONI2(bitnumber)*s2 ...
        + bpsk_SONI3(bitnumber)*s3 + noise(bitnumber);
        
    y = weight'*u;%desire_signal*(weight'*s_SOI) + weight'*ones(L,1)*noise (bitnumber);
    e = desire_signal-y;
    weight = weight + mu/(u'*u+10e-8)*u*e';

SINR = Eb*abs(weight'*s_SOI)/(Eb_SONI1*abs(weight'*s1)+Eb_SONI2*abs(weight'*s2)...
       +Eb_SONI3*abs(weight'*s3)+ N0*sqrt(weight'*weight) ) 
if SINR>0.7
    break;
end
%% Recieving 

bpsk_signal_AWGN = bpsk_signal*(weight'*s_SOI) +bpsk_SONI1*(weight'*s1) ...
                   + bpsk_SONI2*(weight'*s2) + bpsk_SONI3*(weight'*s3)...
                   + noise*(weight'*ones(L,1));

%% Decision & Demod.


Demod_signal=zeros(length(bpsk_signal_AWGN),1);
for index=1:length(bpsk_signal_AWGN)
    if bpsk_signal_AWGN(index)> (0)
        Demod_signal(index)=1;
    else
        Demod_signal(index)=0;
    end
        
end


polarplot(theta,abs(s'*weight),'m')
figure(1)
title('Array Pattern of Cross Array')
%% show image
img=Demod_signal;
img=uint8( reshape(img,[size(1)*size(2) 8]) );
img=reshape(bi2de(img),size);
imshow(img)
figure(2)
end

toc
