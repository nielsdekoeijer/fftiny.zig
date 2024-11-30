# fftiny
A zig implementation of 'A Depth-First Iterative Algorithm for the Conjugate Pair Fast Fourier Transform' by Becoulet et al.
Here, I am not aiming to be the fastest FFT algorithm, but just to provide a portable FFT I can dogfeed in my own projects.
To this end I aim to provide only one dimensional `cfft` and `rfft` interfaces (thus keeping it 'tiny', so to say).
