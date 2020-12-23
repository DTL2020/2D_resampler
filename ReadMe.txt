Sinc and Jinc (EWA Lanczos and other) kernels processed by 2D convolution resampler. It is upsample only for now. It accepts 2 command line params:

Have optional performance test part. All options are switchable by comment in and out parts of code before compile time only.

It loads in.tif as input and outputs out.tif. 8 bit uncompessed. Processes only 1 (red) channel now.

It requires VS2019 for build or any other C compiler with C17.0 math library with bessel_j function.

