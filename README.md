# fftiny
A zig implementation of 'A Depth-First Iterative Algorithm for the Conjugate Pair Fast Fourier Transform' by Becoulet et al.
Here, I am not aiming to be the fastest FFT algorithm, but just to provide a portable FFT I can dogfeed in my own projects.
To this end I aim to provide only one dimensional `cfft` and `rfft` interfaces (thus keeping it 'tiny', so to say).

NOTE: HEAVY WIP DONT USE

## Benchmarking
I've created a benchmark vs FFTW, a well known FFT algorithm. Just run 
```zig
zig build bench -Doptimize=ReleaseFast # optimization recommended...!
```

### Latest
For my own comparision purposes as I slowly improve things (with luck), my latest benchmark yielded me:
```
Starting size 8...
->       fwFFTIterative:  +53509188 ns, mean:       5.35 ns
->       fwFFTRecursive:  +53165311 ns, mean:       5.32 ns
->               fwFFTW: +102397641 ns, mean:      10.24 ns

Starting size 16...
->       fwFFTIterative: +196957426 ns, mean:      19.70 ns
->       fwFFTRecursive: +197684084 ns, mean:      19.77 ns
->               fwFFTW: +339632590 ns, mean:      33.96 ns

Starting size 32...
->       fwFFTIterative: +490274206 ns, mean:      49.03 ns
->       fwFFTRecursive: +487848452 ns, mean:      48.78 ns
->               fwFFTW: +365264862 ns, mean:      36.53 ns

Starting size 64...
->       fwFFTIterative: +1051167220 ns, mean:     105.12 ns
->       fwFFTRecursive: +1064958684 ns, mean:     106.50 ns
->               fwFFTW: +340410747 ns, mean:      34.04 ns
```
My system:
```
niels@pc:~/fftiny.zig$ lscpu
Architecture:             x86_64
  CPU op-mode(s):         32-bit, 64-bit
  Address sizes:          48 bits physical, 48 bits virtual
  Byte Order:             Little Endian
CPU(s):                   16
  On-line CPU(s) list:    0-15
Vendor ID:                AuthenticAMD
  Model name:             AMD Ryzen 7 7840HS w/ Radeon 780M Graphics
    CPU family:           25
    Model:                116
    Thread(s) per core:   2
    Core(s) per socket:   8
    Socket(s):            1
    Stepping:             1
    BogoMIPS:             7585.29
```

