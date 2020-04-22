all:
	LD_LIBRARY_PATH="." gcc -o cyperus_lowpass_module cyperus_lowpass_module.c rt_nonfinite.c rtGetInf.c rtGetNaN.c -lm
