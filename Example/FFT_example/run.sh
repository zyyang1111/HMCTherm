#!/bin/sh
$(your multi2sim directory)/bin/m2s --x86-sim detailed --x86-report ./pipeline.out --mem-report ./mem.out --x86-config ../sim_info/cpuconfig --net-config ../sim_info/netconfig --mem-config ../sim_info/memconfig --net-report net.out --x86-max-inst 1000000000 ../splash2/codes/kernels/fft/FFT -m18 -p16 -n65536 -l4 > mem_trace.tr
