How to use speech coder:

Open the "command" prompt and cd to the directory where you have extracted the files.

1. The program "cp_amrwb_encoder.exe" encodes a .wav file

	For example: 
			 8 the_voice_of_peace.wav the_voice_of_peace.bin 
		encodes the wav file to a binary file at rate of 23.85Kbit/Sec. 
		Encoding rate can be as low as 6.6 Kbit/Sec with reasonable quality.
	
	To get all parameters simply type: cp_amrwb_encoder

2. The program "cp_amrwb_decoder.exe" decodes a binary file.

	For example:
	cp_amrwb_decoder  the_voice_of_peace.bin the_voice_of_peace_recovered.bin