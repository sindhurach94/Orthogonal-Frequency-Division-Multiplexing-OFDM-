% 	Code to generate the plot of Symbol Error Rate (SER) Vs Signal to Noise Ratio (SNR) f or LS and MMSE methods:
%plot of snr vs ser
clc;
close all;
clear all;
%Generation of a naive training sequence.
%Assuming BPSK modulation.
X=zeros(64,64);
d=rand(64,1);
      for i=1:64
       if(d(i)>=0.5)
           d(i)=+1;
       else
           d(i)=-1;
       end
    end
 for i=1:64
     X(i,i)=d(i);
 end
%Calculation of G[The channel Matrix]
 %The channnel is:
  tau=[0.5 3.5];%The fractionally spaced taps..
%Generation of the G matrix...
for k=1:64
      s=0;
      for m=1:2
         s=s+(exp(-j*pi*(1/64)*(k+63*tau(m))) * (( sin(pi*tau(m)) / sin(pi*(1/64)*(tau(m)-k)))));
         
      end
g(k)=s/sqrt(64);
end
G=g';%The channel vector is evaluated
H=fft(G);% In the freq domain..
XFG=X*H;
n1=ones(64,1);
n1=n1*0.000000000000000001i;   % to ensure that the function AWGN adds 'complex Gaussian noise'.
noise=awgn(n1,8);   %Assuming the 'channel learning' is happening at 8db.
variance=var(noise);
N=fft(noise);
Y=XFG+N;
 
% Evaluation of the autocovariance matrix of G-Rgg
 
gg=zeros(64,64);
for i=1:64
    gg(i,i)=G(i);
end
gg_myu = sum(gg, 1)/64;                    
gg_mid = gg - gg_myu(ones(64,1),:);        
sum_gg_mid= sum(gg_mid, 1);
Rgg = (gg_mid' * gg_mid- (sum_gg_mid'  * sum_gg_mid) / 64) / (64 - 1);
 
 
%use of the LS and the MMSE algorithms.
 
%EVALUATION OF Hls
%ls=inv(X)*Y;
 
H_ls=(inv(X)) * Y;
Hls=zeros(64,64);
for i=1:64
    Hls(i,i)=H_ls(i);
end
 
 
%EVALUATION OF H_MMSE
%Hmmse=F*Rgg*inv(Rgy)*Y;
 
u=rand(64,64);
F=fft(u)*inv(u);   %The 64 X 64 twiddle factor matrix.
I=eye(64,64);
Rgy=Rgg * F'* X';
Ryy=X * F * Rgg * F' *X' + variance * I;
for i=1:64
    yy(i,i)=Y(i);
end
Gmmse=Rgy * inv(Ryy)* Y;
H_mmse=fft(Gmmse);
for i=1:64
  Hmmse(i,i)=H_mmse(i); 
end
 
% real time simulations.
for n=1:6
 
SNR_send=5*n;
error_count_ls=0;%Clear the error_count..
error_count_mmse=0;%Clear the error_count.
 
%Sending around 1000 data vectors through the channel
%Roughly like 1000 simulations per SNR reading.
for c=1:1000
%Generate Random Data[i/p matrix]
X=zeros(64,64);
d=rand(64,1);
      for i=1:64
       if(d(i)>=0.5)
           d(i)=+1;
       else
           d(i)=-1;
       end
    end
 for i=1:64
     X(i,i)=d(i);
 end
XFG=X*H;%Let it go through the actual channel.
n1=ones(64,1);
n1=n1*0.000000000000000001i;%Just to ensure that the function awgn adds 'complex gaussian noise'.
noise=awgn(n1,SNR_send);
variance=var(noise);
N=fft(noise);
Y=XFG+N;%o/p got by the receiver.
%The receiver begins.
 
% I:LS ESTIMATOR BASED RECEIVER:
 
    %I(k) represents the decision matrix.
    I=inv(Hls)* Y;
     for k=1:64
       
        if(real(I(k))>0)%Putting it through a slicer
            I(k)=1;
         else
            I(k)=-1;
         end
     end 
   for k=1:64
        if(I(k)~=d(k))
            error_count_ls=error_count_ls+1;
        end
    end
 
% I:MMSE ESTIMATOR BASED RECEIVER:
 %I(k) represents the decision matrix.
    I=inv(Hmmse)* Y;
     for k=1:64
       
        if(real(I(k))>0)%Putting it through a slicer
            I(k)=1;
         else
            I(k)=-1;
         end
     end 
   for k=1:64
        if(I(k)~=d(k))
            error_count_mmse=error_count_mmse+1;
        end
    end
 
end%End of the 1000 run simulation.
 
ser_ls(n)=error_count_ls/64000;
ser_mmse(n)=error_count_mmse/64000;
ser_ls
ser_mmse
SNR(n)=SNR_send;
 
end;
 
%Now just the display part.
semilogy(SNR,ser_mmse,'k-');
grid on;
xlabel('SNR in DB');
ylabel('Symbol Error Rate');
title('PLOT OF SNR V/S SER FOR AN OFDM SYSTEM WITH MMSE/LS ESTIMATOR BASED RECEIVERS');
 
hold on;
semilogy(SNR,ser_ls,'b*');
semilogy(SNR,ser_ls,'b-');
semilogy(SNR,ser_mmse,'kv');
grid on;
xlabel('SNR in DB');
ylabel('Symbol Error Rate');
title('PLOT OF SNR V/S SER FOR AN OFDM SYSTEM WITH MMSE/LS ESTIMATOR BASED RECEIVERS');
legend ('MMSE channel estimation','LS channel estimation');

