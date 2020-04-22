all:
	LD_LIBRARY_PATH="." gcc -o cyperus_filter_varslope_lowpass cyperus_filter_varslope_lowpass.c rt_nonfinite.c rtGetInf.c rtGetNaN.c rtqueue.c -lm -ljack -lpthread -lm 
