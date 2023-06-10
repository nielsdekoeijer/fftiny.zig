# fftiny.zig

A tiny functional fft library written in pure zig. Took some inspiration from RustFFT, PFFFT and [this](https://gist.github.com/rygorous/500e48a94c64c4d83c7d) implementation of split radix I found. 
By no means optimal: mainly an intellectual exercise for me, but by all means use it for whatever you'd like. 
It run's about 2 times as slow as RustFFT for me on x86 and 1.5 to 2 times slower than RustFFT on a raspberry pi for size 256 FFTs.
Heavy disclaimer: these are micro-benchmarks so take this with a grain of salt.
