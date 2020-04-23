all:
	LD_LIBRARY_PATH="." gcc -o sandbox sandbox.c cyperus_filter_varslope_lowpass.c dsp_math_utils.c rtqueue.c -lm -ljack -lpthread -lm 
