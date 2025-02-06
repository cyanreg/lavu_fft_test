#FFMPEG_PATH = "use your own path here or as an argument to make"
INCLUDES = -I "$(FFMPEG_PATH)"
LAVU_LIB = "$(FFMPEG_PATH)/libavutil/libavutil.a"
#LAVC_LIB = "$(FFMPEG_PATH)/libavcodec/libavcodec.a"
CFLAGS = -O2 -g -Wall -lpthread
DEFINES = -DUSE_ALLOCA -DCUSTOM_MODES
LIBS = $(LAVU_LIB) -lm -lfftw3 -lfftw3f -lpthread
OBJ = fft_test.c kiss_fft/kiss_fft.c kiss_fft/mdct.c

lavu:
	$(MAKE) -C $(FFMPEG_PATH)

ne10:
	$(CC) $(DEFINES) -DHAVE_ARM_NE10 $(CFLAGS) $(INCLUDES) $(OBJ) kiss_fft/arm/celt_fft_ne10.c kiss_fft/arm/celt_mdct_ne10.c -o $@ $(LIBS) -lNE10

.DEFAULT_GOAL :=
fft_test: lavu $(OBJ)
	$(CC) $(DEFINES) $(CFLAGS) $(INCLUDES) $(OBJ) -o $@ $(LIBS)

clean:
	-rm -f fft_test
