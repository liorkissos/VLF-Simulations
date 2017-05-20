function [ Symbol_stream_Tx,data_PreCoded_Tx,data_Tx,OFDM_config,N_bits,t1_tx ] = Modulator(bit_stream,Coding_config,OFDM_config,Simulation_config,Audio_config)
%GENERATOR Summary of this function goes here
%   Detailed explanation goes here

Debug_flag=0;

%% Inputs
Coding_flag=Coding_config.Coding_flag;
K_coding=Coding_config.K_coding;
N_coding=Coding_config.N_coding;
Interleave_flag=Coding_config.Interleave_flag;

M=OFDM_config.M;
N_data=OFDM_config.N_data;
P_data=OFDM_config.P_data;
N_preamble_CE=OFDM_config.N_preamble_CE;
N_preamble_synch=OFDM_config.N_preamble_synch;
% IF_chain_config
% Simulation_config


Voice_flag=Audio_config.Voice_flag;


N_symbols=Simulation_config.N_symbols;

Identical_symbols=Simulation_config.Identical_symbols;


%% Generation & QAM Modulation
hMod_data = comm.RectangularQAMModulator('ModulationOrder',M,'NormalizationMethod','Average Power','BitInput',1);


if Voice_flag
    q_bits=Audio_config.q_bits;
    Fs_Audio=Audio_config.Fs_Audio;
    
    % q=quantizer([q_bits q_bits-1]);
    q=quantizer([q_bits 0],'mode','ufixed');
    Mic= dsp.AudioRecorder('SampleRate',Fs_Audio,'NumChannels',1);
    % WAV_write = dsp.AudioFileWriter('Filename','Lior.wav','FileFormat','WAV','SampleRate',Fs_Audio,'DataType','inherit','Compressor','none');
    WAV_write = dsp.AudioFileWriter('Filename','Lior.wav','FileFormat','WAV'...
        ,'SampleRate',Fs_Audio,'DataType','int16');
    
    
    % else
    %     N_samples_frame=[];
end

N_data_bits_sSymbol=N_data*log2(M); % data bit per super symbol (=OFDM symbol)

if Coding_flag
    %%% data bits generation
    if Voice_flag
        
        Rcompression=35.39; % the ratio of compressor 0 is ~35.39 and we wish to create less than required
        
        
        N_data_bits_sSymbol=N_data*log2(M); % data bit per super symbol (=OFDM symbol)
        % Nbits_post_compression=3*double(lcm(sym([N_data_bits_sSymbol 8 K_coding]))); % the number of bits to be transmitted; needs to divide by N bits per super symbol(for an integer number of super symbols), by the number of bits in a codeword (for an integer number of codewords),
        % Nbits_post_compression=300312; % in case of reference configuration, in the absence of sybolic toolbox
        % Nbits_post_compression=100104; % in case of reference configuration, in the absence of sybolic toolbox
        Nbits_post_compression=lcm(N_data_bits_sSymbol,lcm(8,K_coding));
        
        
        Nbits_pre_compression=floor(Rcompression*Nbits_post_compression);
        
        if double(Nbits_post_compression/log2(M))>150e3;
            error('too many symbols per packet')
        end
        
        N_samples_frame=8;
        N_frames=floor(Nbits_pre_compression/(N_samples_frame*q_bits));
        
        Mic.SamplesPerFrame=N_samples_frame;
        
        Signal_voice=zeros(N_frames*N_samples_frame,1);
        
        disp(['Sing for ',num2str(N_frames*N_samples_frame*(1/Fs_Audio)),' seconds!']);
        
        tic
        for kk=1:N_frames % +1 in order to hav always more than needed
            Signal_voice_frame=step(Mic);
            step(WAV_write,Signal_voice_frame);
        end
        
        toc
        disp('stop singing now!')
        release(Mic)
        release(WAV_write)
        
        t1_tx=clock; % comm stsrts here and ends after full demodulation of bits and handing them to speaker
        
        [a,b]=system('cp_amrwb_encoder 0 Lior.wav Lior.bin'); % vocoder encoding
        Signal_voice_compressed=fread(fopen('Lior.bin'),'uint8');% we parse with 8 bits word so to enable lower LCM (low common multiplier) than 16 bits words for example, with no damage to quality of signal
        
        if Debug_flag
            fwrite(fopen('Lior.bin','w'),Signal_voice_compressed,'uint8')
            [a,b]=system('cp_amrwb_decoder Lior.bin Lior_recovered.wav');
            
            B=audioread('Lior_recovered.wav');
            Speaker=dsp.AudioPlayer('SampleRate',Fs_Audio);
            step(Speaker,B)
        end
        
        
        data_PreCoded_Tx1=de2bi(Signal_voice_compressed,8); % 8 because we parse the bin file (fread function) as containing uint8 words
        data_PreCoded_Tx=reshape(data_PreCoded_Tx1',[],1);
        N_bits=length(data_PreCoded_Tx);
        
        tail=N_data_bits_sSymbol-mod(length(data_PreCoded_Tx),N_data_bits_sSymbol); % to accomplish to an integer number of super symbols
        data_PreCoded_Tx=[data_PreCoded_Tx;zeros(tail,1)]; % TEMP: chop to a length fitting an OFDM frame
        N_QAM_symbols=N_bits/log2(M)
        
    else % Not voice signal
        if isempty(bit_stream)
            N_bits=ceil(N_symbols/(K_coding*N_data))*(N_data*K_coding)*log2(M);
            data_PreCoded_Tx = randi([0 1],N_bits,1); %
        else
            LCM=lcm(N_data_bits_sSymbol,lcm(8,K_coding)); % the number of bits to be transmitted; needs to divide by N bits per super symbol(for an integer number of super symbols), by the number of bits in a codeword (for an integer number of codewords),
            
            if LCM>=length(bit_stream)
                tail=LCM-length(bit_stream);
            else
                tail=LCM-mod(length(bit_stream),LCM);
            end
            
            data_PreCoded_Tx=[bit_stream;zeros(tail,1)];
            N_bits=length(data_PreCoded_Tx);
            
        end
        t1_tx=clock; % comm starts here and ends after full demodulation of bits and handing them to speaker
        
    end % Voice
    
    
    %%% Encoding
    enc = comm.RSEncoder('BitInput',true,'CodewordLength',N_coding,'MessageLength',K_coding);
    data_Tx = step(enc,data_PreCoded_Tx);
    %%% Interleaving
    
    if Interleave_flag
        data_Tx=matintrlv(data_Tx,length(data_Tx)/N_coding,N_coding);
    end
    %%% QAM modulating
    Symbol_stream_Tx = step(hMod_data, data_Tx); % QAM modulation
    
    %%%%% power ratio
    Symbol_stream_PreCoded_Tx=step(hMod_data,data_PreCoded_Tx); % generated only for the sake of correct addition of noise
    
    E00=sum(sum((abs(Symbol_stream_PreCoded_Tx)).^2)); % Only payload subcarriers Energy, before concatenatin of DC zeros
    E01=sum(sum((abs(Symbol_stream_Tx)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
    R0=E00/E01; % should result in 1, aince guard bands do not add any power
    
else % no coding
    
    if Voice_flag
        
        Rcompression=35.39; % the ratio of compressor 0 is ~35.39 and we wish to create less than required
        N_samples_frame=8;
        
        N_data_bits_sSymbol=N_data*log2(M); % data bit per super symbol (=OFDM symbol)
        Nbits_post_compression=70*double(lcm(sym([N_data_bits_sSymbol 8])));
        Nbits_pre_compression=floor(Rcompression*Nbits_post_compression);
        
        if double(Nbits_post_compression/log2(M))>150e3;
            error('too many symbols per packet')
        end
        
        
        N_frames=floor(Nbits_pre_compression/(N_samples_frame*q_bits));
        
        Mic.SamplesPerFrame=N_samples_frame;
        
        Signal_voice=zeros(N_frames*N_samples_frame,1);
        
        disp(['Sing for ',num2str(N_frames*N_samples_frame*(1/Fs_Audio)),' seconds!']);
        
        tic
        for kk=1:N_frames % +1 in order to hav always more than needed
            Signal_voice_frame=step(Mic);
            step(WAV_write,Signal_voice_frame);
        end
        
        toc
        disp('stop singing now!')
        release(Mic)
        release(WAV_write)
        
        t1_tx=clock; % comm stsrts here and ends after full demodulation of bits and handing them to speaker
        
        [a,b]=system('cp_amrwb_encoder 0 Lior.wav Lior.bin'); % vocoder encoding
        Signal_voice_compressed=fread(fopen('Lior.bin'),'uint8');
        
        if Debug_flag
            fwrite(fopen('Lior.bin','w'),Signal_voice_compressed,'uint8')
            [a,b]=system('cp_amrwb_decoder Lior.bin Lior_recovered.wav');
            
            B=audioread('Lior_recovered.wav');
            Speaker=dsp.AudioPlayer('SampleRate',Fs_Audio);
            step(Speaker,B)
        end
        
        
        data_Tx1=de2bi(Signal_voice_compressed,8);
        data_Tx=reshape(data_Tx1',[],1);
        N_bits=length(data_Tx);
        
        tail=N_data_bits_sSymbol-mod(length(data_Tx),N_data_bits_sSymbol);
        data_Tx=[data_Tx;zeros(tail,1)]; % TEMP: chop to a length fitting an OFDM frame
        N_QAM_symbols=N_bits/log2(M)
        
    else
        if isempty(bit_stream)
            N_bits=round(N_symbols/N_data)*N_data*log2(M);
            data_Tx = randi([0 1],N_bits,1); % N_data is the # of sc's per OFDM symbol. Thus,round(N_symbols/N_data)*N_data is a number of QAM symbols that can be transmitted within an integer number of OFDM (!) symbols
        else
            LCM=lcm(N_data_bits_sSymbol,8); % the number of bits to be transmitted; needs to divide by N bits per super symbol(for an integer number of super symbols), by the number of bits in a codeword (for an integer number of codewords),
            
            if LCM>length(bit_stream)
                tail=LCM-length(bit_stream);
            else
                tail=LCM-mod(length(bit_stream),LCM);
            end
            
            data_Tx=[bit_stream;zeros(tail,1)];
            N_bits=length(data_Tx);
            
        end
        
        
    end
    
    t1_tx=clock; % comm stsrts here and ends after full demodulation of bits and handing them to speaker
    
    Symbol_stream_Tx = step(hMod_data, data_Tx); % QAM modulation
    
    data_PreCoded_Tx=[];
    R0=1;
    
end


OFDM_config.R0=R0; % passed to Channel_SW function ofr a correct addition of noise

if Identical_symbols
    data_Tx=repmat(data_Tx(1:N_data*log2(M)),length(data_Tx)/(N_data*log2(M)),1);
end

Symbol_stream_Tx=Symbol_stream_Tx*(sqrt(P_data)/sqrt(var(Symbol_stream_Tx))); % enforcement of the required average power. derived from the P_total parameter

%N_frames_Tx=length(Symbol_stream_Tx)/N_data;
%N_frames_Rx=N_frames_Tx+N_preamble_CE+N_preamble_synch;
end

