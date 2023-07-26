#!/usr/bin/gnuplot
set terminal postscript color eps
set output "jacobi_kernel_titane.eps"
set logscale x
set grid
set xtics nomirror
set ytics nomirror
set xrange[1000:20000*36]
set xlabel "N cells"
set ylabel "GFlops/core"
set title "Best SSE Jacobi kernel (10flops/cell) on X5570 (Titane, gcc 4.5.0)"
set size 0.8
plot "jacobi_kernel_nehalem.txt" using ($1*36):2 title "8 threads, double precision" w lp lw 2,\
     "jacobi_kernel_nehalem.txt" using ($1*36):3 title "8 threads, single precision" w lp lw 2,\
"jacobi_kernel_nehalem_1thread.txt" using ($1*36):2 title "1 thread, double precision" w lp lw 2,\
"jacobi_kernel_nehalem_1thread.txt" using ($1*36):3 title "1 thread, single precision" w lp lw 2

set output "jacobi_kernel_titane_bw.eps"
set title "SSE Jacobi kernel (10flops/cell) on X5570: extra L1 bandwidth"
set ylabel "GB/s"
plot "jacobi_kernel_nehalem.txt" using ($1*36):($2/10*6*8*8) title "8 threads, double precision" w lp lw 2,\
     "jacobi_kernel_nehalem.txt" using ($1*36):($3/10*6*4*8) title "8 threads, single precision" w lp lw 2

!epstopdf jacobi_kernel_titane.eps
!epstopdf jacobi_kernel_titane_bw.eps