all:
	LD_LIBRARY_PATH="." gcc -o sandbox sandbox.c modules/dsp/filter/varslope_lowpass/modules_dsp_filter_varslope_lowpass.c modules/dsp/math_utils.c rtqueue.c -lm -ljack -lpthread -lm 
