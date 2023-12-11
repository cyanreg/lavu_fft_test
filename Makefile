#FFMPEG_PATH = "use your own path here or as an argument to make"
INCLUDES = -I "$(FFMPEG_PATH)"
LAVU_LIB = "$(FFMPEG_PATH)/libavutil/libavutil.a"
#LAVC_LIB = "$(FFMPEG_PATH)/libavcodec/libavcodec.a"
CFLAGS = -O0 -g -Wall -lpthread
LIBS = $(LAVU_LIB) -lm -lfftw3 -lfftw3f -lpthread
OBJ = fft_test.c

lavu:
	$(MAKE) -C $(FFMPEG_PATH)

.DEFAULT_GOAL :=
fft_test: lavu $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJ) -o $@ $(LIBS)

clean:
	-rm -f fft_test
