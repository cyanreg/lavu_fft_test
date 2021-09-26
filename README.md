Development, testing and benchmarking code for FFmpeg's libavutil/tx
--------------------------------------------------------------------

Tests and benchmarks the new FFT in FFmpeg.
To configure FFmpeg:
``` bash
./configure --disable-doc --disable-shared --enable-static --disable-ffplay --disable-ffmpeg --disable-ffprobe --disable-everything --disable-swresample --disable-avcodec --disable-avfilter --disable-avdevice --disable-swresample --enable-avutil --disable-swscale --disable-autodetect --disable-iconv --disable-stripping --extra-cflags='-mtune=native -march=native' --optflags=-Og --cc='clang-13 -fuse-ld=lld' --nm=llvm-nm-13 --ar=llvm-ar-13 --ranlib=llvm-ranlib-13
```

To compile FFmpeg, the program, and run it:
``` bash
FFMPEG_PATH="path_to_your_configured_ffmpeg_folder make"
```

You want the RMSE value to remain **unchanged** from the C version
