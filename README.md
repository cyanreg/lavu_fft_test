Development, testing and benchmarking code for FFmpeg's libavutil/tx
--------------------------------------------------------------------

Directory structure:
```
./some_directory/
    ./ffmpeg/
    ./lavu_fft_test/
```

Apply the patch `0001-lavu-tx-add-placeholder-aarch64-double-precision-opt.patch` to ffmpeg.

Run ./configure in ./ffmpeg/ with the following arguments:
```
./configure --disable-doc --disable-shared --enable-static --disable-ffplay --disable-ffmpeg --disable-ffprobe --disable-everything --disable-swresample --disable-avcodec --disable-avfilter --disable-avdevice --disable-swresample --enable-avutil --disable-swscale --disable-autodetect --disable-iconv --disable-stripping --extra-cflags='-mtune=native -march=native' --optflags=-Og --disable-pthreads --enable-linux-perf
```

Run make in `./ffmpeg/`

To run the test program:
```
cd ../ffmpeg ; make ; cd ../lavu_fft_test ; cc fft_test.c ../ffmpeg/libavutil/libavutil.a -lm -I../ffmpeg/ && ./a.out
```

You want the RMSE value to remain **unchanged** from the C version
