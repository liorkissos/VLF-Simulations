
clc
clear
close all

%Fs_Audio=88.2e3;
Fs_Audio=44.1e3;

N_samples_frame=2^10;
q_bits=16;

Trecord=2.5;

%%
Mic=dsp.AudioRecorder('SampleRate',Fs_Audio,'NumChannels',1,'SamplesPerFrame',N_samples_frame);
Speaker=dsp.AudioPlayer('SampleRate',Fs_Audio); % must have the same Fs as Microphone so that the voice does not distort


%% Option 2

WAV_write = dsp.AudioFileWriter('Filename','Lior.wav','FileFormat','WAV'...
    ,'SampleRate',Fs_Audio,'DataType','int16');
% WAV_write = dsp.AudioFileWriter('Filename','Lior.wav','FileFormat','WAV'...
%     ,'SampleRate',Fs_Audio,'DataType','inherit');

WAV_read=dsp.AudioFileReader('Filename','Lior_recovered.wav'...
    ,'OutputDataType','int16');

%% 1)Record from Mic
disp(['Sing into microphone for ',num2str(Trecord),' seconds!'])

pause(1)
A=0;
tic
while toc<Trecord
    a=step(Mic);
    A=[A;a];
end
toc


    %A=step(Mic);


disp('Stop singing!')
release(Mic)


%% 2) Write to wav file
%audiowrite('Lior.wav',A,Fs_Audio,'BitsPerSample',q_bits)
step(WAV_write,A)

%% 3) Encode: wav->bin
[a,b]=system('cp_amrwb_encoder 8 Lior.wav Lior.bin');

C=fread(fopen('Lior.bin'),'uint16');

%% 4) decode: bin->wav
fwrite(fopen('Lior_r_w.bin','w'),C,'uint16')

[a,b]=system('cp_amrwb_decoder Lior_r_w.bin Lior_recovered.wav');

%% 5) Read from wav file
%B=audioread('Lior_recovered.wav');
B=0;
while ~isDone(WAV_read)
    b=step(WAV_read);
    B=[B;b];
end

%% 6) Send to Speaker
step(Speaker,B)
pause(3)
release(Speaker)

%% compression calvulation

Compression_factor=length(A)/length(C)
N_symbols_64QAM=length(C)*16/log2(64)
