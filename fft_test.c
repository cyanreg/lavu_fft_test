#include <libavutil/cpu.h>
#include <libavutil/mem.h>
#include <libavutil/lfg.h>
#include <libavutil/tx.h>
#include <libavutil/random_seed.h>

#if 1
#include <libavutil/timer.h>
#else
#define START_TIMER
#define STOP_TIMER(a)
#endif

#define AVFFT     0
#define FFTW      0

#define FFT_LEN   8
#define DOUBLE    1
#define REPS     (1 << 0)
#define IN_PLACE  0
#define NO_SIMD   0
#define INVERSE   0

#if DOUBLE == 1
#undef AVFFT
#define AVFFT 0
#endif

#if AVFFT
#include <libavcodec/avfft.h>
#endif

#if FFTW
#include <fftw3.h>
#endif

#if DOUBLE
#undef AVFFT
#define AVFFT 0
typedef AVComplexDouble TXComplex;
#else
typedef AVComplexFloat TXComplex;
#endif

TXComplex cmult(TXComplex z1, TXComplex z2)
    { TXComplex res; res.re = z1.re*z2.re - z1.im*z2.im; res.im = z2.re*z1.im + z1.re*z2.im; return res; }
TXComplex cadd(TXComplex z1, TXComplex z2)
    { TXComplex res; res.re = z1.re + z2.re; res.im = z1.im + z2.im; return res; }

void do_naive_tx(int inv, TXComplex *output, TXComplex *input, int N)
{
    double phase = inv ? 2*M_PI : -2*M_PI;

    for(int i = 0; i < N; i++) {
        TXComplex tmp = { 0 };
        for(int j = 0; j < N; j++) {
            TXComplex twiddle = { cos(phase*j*i/N), sin(phase*j*i/N), };
        	tmp = cadd(tmp, cmult(input[j], twiddle));
        }
        output[i] = tmp;
    }
}

av_tx_fn tx;
void do_avfft_tx(AVTXContext *s, TXComplex *output, TXComplex *input, int len)
{
    for (int i = 0; i < REPS; i++) {
        if (IN_PLACE)
            memcpy(output, input, len*sizeof(TXComplex));

        START_TIMER
        tx(s, output, IN_PLACE ? output : input, sizeof(TXComplex));
        STOP_TIMER("             av_tx_fn");
    }
}

#if AVFFT
void do_lavc_tx(FFTContext *avfft, TXComplex *output, TXComplex *input, int len)
{
    for (int i = 0; i < REPS; i++) {
        START_TIMER
        memcpy(output, input, len*sizeof(TXComplex));
        av_fft_permute(avfft, (FFTComplex *)output);
        av_fft_calc(avfft, (FFTComplex *)output);
        STOP_TIMER("  av_fft_permute+calc");
    }
}
#endif

#if FFTW
#if DOUBLE
void do_fftw_tx(fftw_plan fftw_plan, TXComplex *output, TXComplex *input, int len)
#else
void do_fftw_tx(fftwf_plan fftw_plan, TXComplex *output, TXComplex *input, int len)
#endif
{
    for (int i = 0; i < REPS; i++) {
        if (IN_PLACE)
            memcpy(output, input, len*sizeof(TXComplex));
        START_TIMER
#if DOUBLE
        fftw_execute(fftw_plan);
#else
        fftwf_execute(fftw_plan);
#endif
        STOP_TIMER("        fftwf_execute");
    }
}
#endif

void compare_results(TXComplex *ref, TXComplex *dut, int len, const char *str)
{
    double err = 0.0;
    for (int i = 0; i < len; i++) {
        err += (ref[i].re - dut[i].re)*(ref[i].re - dut[i].re);
        err += (ref[i].im - dut[i].im)*(ref[i].im - dut[i].im);
    }

    printf("RMSE %s = %f\n", str, sqrtf(err / (len*2)));
}

int main(void)
{
    int inv = INVERSE;
    int tx_len = FFT_LEN;

    if (NO_SIMD)
        av_force_cpu_flags(0);

    AVTXContext *avfftctx;
    TXComplex scale = { 1.0f };
    int ret = av_tx_init(&avfftctx, &tx, DOUBLE ? AV_TX_DOUBLE_FFT : AV_TX_FLOAT_FFT,
                         inv, tx_len, &scale.re,
                         (IN_PLACE ? AV_TX_INPLACE : 0x0) |
                         (NO_SIMD ? AV_TX_UNALIGNED : 0x0));
    if (ret) {
        av_log(NULL, AV_LOG_WARNING, "ERRROR = %i\n", ret);
        return 1;
    }

    TXComplex *input = av_mallocz(sizeof(TXComplex)*tx_len);
#if AVFFT
    TXComplex *output_lavc = av_mallocz(sizeof(TXComplex)*tx_len);
#endif
#if FFTW
    TXComplex *output_fftw = av_mallocz(sizeof(TXComplex)*tx_len);
#endif

    TXComplex *output_new = av_mallocz(sizeof(TXComplex)*tx_len);
    TXComplex *output_naive = av_mallocz(sizeof(TXComplex)*tx_len);

#if AVFFT
    FFTContext *avfft = av_fft_init(av_log2(tx_len), inv);
#endif
#if FFTW
#if DOUBLE
    fftw_plan fftw_plan = fftw_plan_dft_1d  (tx_len, IN_PLACE ?
                                                     (fftw_complex *)output_fftw :
                                                     (fftw_complex *)input,
                                             (fftw_complex *)output_fftw,
                                             inv ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_MEASURE |
                                             NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftwf_plan fftw_plan = fftwf_plan_dft_1d(tx_len, IN_PLACE ?
                                                     (fftwf_complex *)output_fftw :
                                                     (fftwf_complex *)input,
                                             (fftwf_complex *)output_fftw,
                                             inv ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_MEASURE |
                                             NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#endif

//    srand(av_get_random_seed());

    for (int i = 0; i < tx_len; i++) {
        TXComplex gen = {
            (rand()/(double)RAND_MAX),
            (rand()/(double)RAND_MAX),
        };
        input[i] = gen;
    }

    do_avfft_tx(avfftctx, output_new, input, tx_len);
#if AVFFT
    if (!(tx_len & (tx_len - 1)))
        do_lavc_tx(avfft, output_lavc, input, tx_len);
#endif
#if FFTW
    do_fftw_tx(fftw_plan, output_fftw, input, tx_len);
#endif
    do_naive_tx(inv, output_naive, input, tx_len);

    compare_results(output_naive, output_new, tx_len, " avtx");
#if FFTW
    compare_results(output_naive, output_fftw, tx_len, " fftw");
#endif
#if AVFFT
    if (!(tx_len & (tx_len - 1)))
        compare_results(output_naive, output_lavc, tx_len, "avfft");
#endif

#if FFTW
#if DOUBLE
    fftw_destroy_plan(fftw_plan);
#else
    fftwf_destroy_plan(fftw_plan);
#endif
#endif

#if AVFFT
    av_fft_end(avfft);
#endif
    av_tx_uninit(&avfftctx);

    av_free(input);
#if AVFFT
    av_free(output_lavc);
#endif
    av_free(output_naive);
#if FFTW
    av_free(output_fftw);
#endif
    av_free(output_new);

    return 0;
}
