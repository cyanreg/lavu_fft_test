#include <libavutil/cpu.h>
#include <libavutil/mem.h>
#include <libavutil/lfg.h>
#include <libavutil/tx.h>
#include <libavutil/time.h>
#include <libavutil/random_seed.h>

#if 1
#include <libavutil/timer.h>
#else
#define START_TIMER
#define STOP_TIMER(x)
#endif

#define FFT       0 // Standard complex-to-complex FFT
#define MDCT      1 // Standard MDCT
#define R2C       2
#define C2R       3

#define MODE C2R
#define INVERSE   1 // Transform direction (only for MDCT and FFT)

/* Double-precision instead of float */
#define DOUBLE    0

#define FFTW      1
#define AVFFT     1

#define FFT_LEN   32
#define REPS     (1 << 0)
#define IN_PLACE  0
#define NO_SIMD   1
#define IMDCT_F   0
#define SEED      0x8063
#define PRINTOUT  1

#if MODE == MDCT
#undef FFTW
#define FFTW 0
#endif

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
typedef double TXSample;
typedef AVComplexDouble TXComplex;
#else
typedef float TXSample;
typedef AVComplexFloat TXComplex;
#endif

static TXComplex cmult(TXComplex z1, TXComplex z2)
{
    TXComplex res;
    res.re = z1.re*z2.re - z1.im*z2.im;
    res.im = z2.re*z1.im + z1.re*z2.im;
    return res;
}

static TXComplex cadd(TXComplex z1, TXComplex z2)
{
    TXComplex res;
    res.re = z1.re + z2.re;
    res.im = z1.im + z2.im;
    return res;
}

static void naive_fft(TXComplex *output, TXComplex *input, int len)
{
    double phase = INVERSE ? 2*M_PI : -2*M_PI;

    for (int i = 0; i < len; i++) {
        TXComplex tmp = { 0 };
        for(int j = 0; j < len; j++) {
            TXComplex twiddle = { cos(phase*j*i/len), sin(phase*j*i/len) };
            tmp = cadd(tmp, cmult(input[j], twiddle));
        }
        output[i] = tmp;
    }
}

static void naive_rdft_r2c(TXComplex *output, TXComplex *input, int len)
{
    TXSample *in = (TXSample *)input;
    double phase = -2*M_PI;

    for (int i = 0; i < (len + 1); i++) {
        TXComplex tmp = { 0 };
        for(int j = 0; j < len*2; j++) {
            TXSample n = in[j];
            TXComplex twiddle = { cos(phase*j*i/(len*2)), sin(phase*j*i/(len*2)) };
            TXComplex res = { n * twiddle.re, n * twiddle.im };
        	tmp = cadd(tmp, res);
        }
        output[i] = tmp;
    }
}

static void naive_rdft_c2r(TXComplex *output, TXComplex *input, int len)
{
    TXSample *out = (TXSample *)output;
    double phase = 2*M_PI;
    int len_complex = len + 1;

    len *= 2;

    for (int i = 0; i < len; i++) {
        /* DC value */
        out[i] = input[0].re;

        int end = len_complex;

        /* Even lengths are special cased */
        if (!(len & 1)) {
            double tw = cos((phase*i*(len_complex - 1))/(len));
            out[i] += input[len_complex - 1].re * tw;
            end = len_complex - 1;
        }

        for (int j = 1; j < end; j++) {
            double phi = phase*j*i/(len);
            TXComplex tw = { cos(phi), sin(phi) };

            /* No need for a full complex multiply. */
            double val = input[j].re*tw.re - input[j].im*tw.im;
            out[i] += 2.0 * val;
        }
    }
}

static void naive_imdct(TXSample *dst, TXSample *src, int len)
{
    len >>= 1;

    int len2 = len*2;
    const double phase = M_PI/(4.0*len2);

    for (int i = 0; i < len; i++) {
        double sum_d = 0.0;
        double sum_u = 0.0;
        double i_d = phase * (4*len  - 2*i - 1);
        double i_u = phase * (3*len2 + 2*i + 1);
        for (int j = 0; j < len2; j++) {
            double a = (2 * j + 1);
            double a_d = cos(a * i_d);
            double a_u = cos(a * i_u);
            double val = src[j];
            sum_d += a_d * val;
            sum_u += a_u * val;
        }
        dst[i +   0] =  sum_d;
        dst[i + len] = -sum_u;
    }
}

static void naive_imdct_full(TXSample *dst, TXSample *src, int len)
{
    const double phase = M_PI/(4.0*len);

    for (int i = 0; i < len*2; i++) {
        double sum = 0.0;
        for (int j = 0; j < len; j++) {
            int a = (2 * i + 1 + len) * (2 * j + 1);
            sum += src[j] * cos(M_PI * a / (double) (4 * len));
        }
        dst[i] = -sum;
    }
}

static void naive_mdct(TXSample *dst, TXSample *src, int len)
{
    const double phase = M_PI/(4.0*len);

    for (int i = 0; i < len; i++) {
        double sum = 0.0;
        for (int j = 0; j < len*2; j++) {
            int a = (2*j + 1 + len) * (2*i + 1);
            sum += src[j] * cos(a * phase);
        }
        dst[i] = sum;
    }
}

#if 0
static void expand_output(TXComplex *buffer, int len)
{
    TXSample *input = (TXSample *)(buffer + len/2);
    input--;

    for (int i = (len - 1); i >= 0; i--) {
        buffer[i] = (TXComplex){*input, 0};
        input--;
    }
}
#endif

av_tx_fn tx;
void do_avfft_tx(AVTXContext *s, TXComplex *output, TXComplex *input, int len)
{
#if IN_PLACE
    memcpy(output, input, len*sizeof(TXComplex));
#endif

    int64_t start = av_gettime_relative();

    for (int i = 0; i < REPS; i++) {
        START_TIMER
        tx(s, output, IN_PLACE ? output : input, sizeof(TXComplex) >> MDCT);
#if MODE == MDCT
#if INVERSE
        STOP_TIMER("        av_tx (imdct)");
#else
        STOP_TIMER("         av_tx (mdct)");
#endif
#elif MODE == C2R
        STOP_TIMER("          av_tx (c2r)");
#elif MODE == R2C
        STOP_TIMER("          av_tx (r2c)");
#else
#if INVERSE
        STOP_TIMER("         av_tx (ifft)");
#else
        STOP_TIMER("          av_tx (fft)");
#endif
#endif
    }

    if (REPS > 1)
        printf("Total for len %i reps %i = %f s\n", len, REPS,
               (av_gettime_relative() - (double)start)/1000000.0);
}

#if AVFFT
#if MODE == C2R || MODE == R2C
void do_lavc_tx(RDFTContext *avfft, TXComplex *output, TXComplex *input, int len)
#else
void do_lavc_tx(FFTContext *avfft, TXComplex *output, TXComplex *input, int len)
#endif
{
#if MODE == C2R
    TXSample tmp = input[len].re;
    input[len].re = input[0].im;
    input[0].im = tmp;
#endif

#if (MODE != MDCT) && IN_PLACE
    memcpy(output, input, len*sizeof(TXComplex));
#endif

    int64_t start = av_gettime_relative();

    for (int i = 0; i < REPS; i++) {
        START_TIMER
#if (MODE != MDCT) && !IN_PLACE
        memcpy(output, input, len*sizeof(TXComplex));
#endif

#if MODE == MDCT
#if INVERSE
#if IMDCT_F
        av_imdct_calc(avfft, (FFTSample *)output, (const FFTSample *)input);
        STOP_TIMER("        av_imdct_calc");
#else
        av_imdct_half(avfft, (FFTSample *)output, (const FFTSample *)input);
        STOP_TIMER("        av_imdct_half");
#endif
#else
        av_mdct_calc(avfft, (FFTSample *)output, (const FFTSample *)input);
        STOP_TIMER("         av_mdct_calc");
#endif
#elif MODE == C2R || MODE == R2C
        av_rdft_calc(avfft, (FFTSample *)output);
        STOP_TIMER("         av_rdft_calc");
#else
        av_fft_permute(avfft, (FFTComplex *)output);
        av_fft_calc(avfft, (FFTComplex *)output);
        STOP_TIMER("  av_fft_permute+calc");
#endif
    }

    if (REPS > 1)
        printf("Total for len %i reps %i = %f s\n", len, REPS,
               (av_gettime_relative() - (double)start)/1000000.0);

#if MODE == C2R
    for (int i = 0; i < len; i++) {
        output[i].re *= 2;
        output[i].im *= 2;
    }
#elif MODE == R2C
    TXSample tmp = output[len].re;
    output[len].re = output[0].im;
    output[0].im = tmp;
#endif
}
#endif

#if FFTW
#if DOUBLE
void do_fftw_tx(fftw_plan fftw_plan, TXComplex *output, TXComplex *input, int len)
#else
void do_fftw_tx(fftwf_plan fftw_plan, TXComplex *output, TXComplex *input, int len)
#endif
{
#if IN_PLACE
    memcpy(output, input, len*sizeof(TXComplex));
#endif

    int64_t start = av_gettime_relative();

    for (int i = 0; i < REPS; i++) {
        START_TIMER
#if DOUBLE
        fftw_execute(fftw_plan);
        STOP_TIMER("         fftw_execute");
#else
        fftwf_execute(fftw_plan);
        STOP_TIMER("        fftwf_execute");
#endif
    }

    if (REPS > 1)
        printf("Total for len %i reps %i = %f s\n", len, REPS,
               (av_gettime_relative() - (double)start)/1000000.0);
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

void compare_all(TXComplex *all[], int len, const char *names[], int num_arr)
{
    printf("Dump:\n\t\t\t");
    for (int j = 0; j < num_arr; j++) {
        printf("%s%s", names[j], j == (num_arr - 1) ? "\n" : "\t\t\t");
    }

    for (int i = 0; i < len; i++) {
        printf("\t%02i. ", i + 1);
        for (int j = 0; j < num_arr; j++) {
            printf("{%f, %f}%s", all[j][i].re, all[j][i].im,
                   j == (num_arr - 1) ? "\n" : "\t");
        }
    }
}

int main(void)
{
    int inv = INVERSE;
    int tx_len = FFT_LEN;

    if (NO_SIMD)
        av_force_cpu_flags(0);

#if MODE == C2R || MODE == R2C
    printf("Mode: %s %s DFT %s\n",
           DOUBLE ? "double" : "float",
           MODE == C2R ? "complex-to-real" : "real-to-complex",
           NO_SIMD ? "(no SIMD)" : "(with SIMD)");
#else
    printf("Mode: %s %s %s %s\n",
           DOUBLE ? "double" : "float",
           INVERSE ? "inverse" : "forward",
           MDCT ? IMDCT_F ? "mdct (half)" : "mdct" : "fft",
           NO_SIMD ? "(no SIMD)" : "(with SIMD)");
#endif

    uint32_t seed = !SEED ? av_get_random_seed() : SEED;
    printf("Seed: 0x%x\n", seed);

    AVTXContext *avfftctx;
    TXSample scale = 1.0;
    int ret = av_tx_init(&avfftctx, &tx,
#if MODE == MDCT
                         DOUBLE ? AV_TX_DOUBLE_MDCT : AV_TX_FLOAT_MDCT,
#elif MODE == R2C
                         DOUBLE ? AV_TX_FLOAT_R2C_DFT : AV_TX_DOUBLE_R2C_DFT,
#elif MODE == C2R
                         DOUBLE ? AV_TX_FLOAT_C2R_DFT : AV_TX_DOUBLE_R2C_DFT,
#else
                         DOUBLE ? AV_TX_DOUBLE_FFT : AV_TX_FLOAT_FFT,
#endif
                         inv, tx_len, &scale,
                         (IN_PLACE ? AV_TX_INPLACE : 0x0) |
                         (IMDCT_F ? AV_TX_FULL_IMDCT : 0x0) |
                         (NO_SIMD ? AV_TX_UNALIGNED : 0x0));
    if (ret) {
        av_log(NULL, AV_LOG_WARNING, "ERROR = %s\n", av_err2str(ret));
        return 1;
    }

    printf("Length: %i\n", tx_len);
    printf("Function: %p\n", tx);

    int num_in = tx_len + (MODE == C2R);
    int num_out = tx_len/(1 + ((MODE == MDCT) && INVERSE && !IMDCT_F)) + (MODE == R2C);

    size_t alloc_in = sizeof(TXComplex)*num_in;
    size_t alloc_out = sizeof(TXComplex)*num_out;

    TXComplex *input = av_mallocz(alloc_in);

#if AVFFT
    TXComplex *input_lavc = av_mallocz(alloc_in);
    TXComplex *output_lavc = av_mallocz(alloc_out);
#endif
#if FFTW
    TXComplex *input_fftw = av_mallocz(alloc_in);
    TXComplex *output_fftw = av_mallocz(alloc_out);
#endif

    TXComplex *output_new = av_mallocz(alloc_out);
    TXComplex *output_naive = av_mallocz(alloc_out);

#if AVFFT
#if MODE == MDCT
    FFTContext *avfft = av_mdct_init(av_log2(tx_len) + 1, inv, scale);
#elif MODE == R2C
    RDFTContext *avfft = av_rdft_init(av_log2(tx_len) + 1, DFT_R2C);
#elif MODE == C2R
    RDFTContext *avfft = av_rdft_init(av_log2(tx_len) + 1, IDFT_C2R);
#else
    FFTContext *avfft = av_fft_init(av_log2(tx_len), inv);
#endif
#endif
#if FFTW
#if MODE == FFT
#if DOUBLE
    fftw_plan fftw_plan = fftw_plan_dft_1d  (tx_len, IN_PLACE ?
                                                     (fftw_complex *)output_fftw :
                                                     (fftw_complex *)input_fftw,
                                             (fftw_complex *)output_fftw,
                                             inv ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_MEASURE |
                                             NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftwf_plan fftw_plan = fftwf_plan_dft_1d(tx_len, IN_PLACE ?
                                                     (fftwf_complex *)output_fftw :
                                                     (fftwf_complex *)input_fftw,
                                             (fftwf_complex *)output_fftw,
                                             inv ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_MEASURE |
                                             NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#elif MODE == R2C
#if DOUBLE
    fftw_plan fftw_plan = fftw_plan_dft_r2c_1d  (tx_len*2, IN_PLACE ?
                                                         (double *)output_fftw :
                                                         (double *)input_fftw,
                                                 (fftw_complex *)output_fftw,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftwf_plan fftw_plan = fftwf_plan_dft_r2c_1d(tx_len*2, IN_PLACE ?
                                                         (float *)output_fftw :
                                                         (float *)input_fftw,
                                                 (fftwf_complex *)output_fftw,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#elif MODE == C2R
#if DOUBLE
    fftw_plan fftw_plan = fftw_plan_dft_c2r_1d  (tx_len*2, IN_PLACE ?
                                                         (fftw_complex *)output_fftw :
                                                         (fftw_complex *)input_fftw,
                                                 (double *)output_fftw,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftwf_plan fftw_plan = fftwf_plan_dft_c2r_1d(tx_len*2, IN_PLACE ?
                                                         (fftwf_complex *)output_fftw :
                                                         (fftwf_complex *)input_fftw,
                                                 (float *)output_fftw,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#endif
#endif

    srand(seed);

    for (int i = 0; i < num_in; i++) {
        TXComplex gen = {
            (rand()/(double)RAND_MAX)/tx_len,
            (rand()/(double)RAND_MAX)/tx_len,
        };

        input[i] = gen;
    }

#if FFTW
    memcpy(input_fftw, input, alloc_in);
#endif
#if AVFFT
    memcpy(input_lavc, input, alloc_in);
#endif

    do_avfft_tx(avfftctx, output_new, input, tx_len);
#if AVFFT
    if (!(tx_len & (tx_len - 1)))
        do_lavc_tx(avfft, output_lavc, input_lavc, tx_len);
#endif
#if FFTW
    do_fftw_tx(fftw_plan, output_fftw, input, tx_len);
#endif
#if MODE == MDCT
#if INVERSE
#if IMDCT_F
    naive_imdct_full((TXSample *)output_naive, (TXSample *)input, tx_len);
#else
    naive_imdct((TXSample *)output_naive, (TXSample *)input, tx_len);
#endif
#else
    naive_mdct((TXSample *)output_naive, (TXSample *)input, tx_len);
#endif
#elif MODE == R2C
    naive_rdft_r2c(output_naive, input, tx_len);
#elif MODE == C2R
    naive_rdft_c2r(output_naive, input, tx_len);
#else
    naive_fft(output_naive, input, tx_len);
#endif

#if PRINTOUT
     compare_all(
        (TXComplex *[]){
            output_naive,
            output_new,
#if FFTW
            output_fftw,
#endif
#if AVFFT
            output_lavc,
#endif
        },
        num_out,
        (const char *[]){
            "naive",
            "lavu",
#if FFTW
            "fftw",
#endif
#if AVFFT
            "lavc",
#endif
        },
        2 + !!FFTW + !!AVFFT);
#endif

    compare_results(output_naive, output_new, num_out,      "  av_tx");

#if FFTW
#if DOUBLE
    compare_results(output_naive, output_fftw, num_out,     "  fftw3");
#else
    compare_results(output_naive, output_fftw, num_out,     " fftw3f");
#endif
#endif

#if AVFFT
    if (!(tx_len & (tx_len - 1)))
        compare_results(output_naive, output_lavc, num_out, "  avfft");
#endif

#if FFTW
#if DOUBLE
    fftw_destroy_plan(fftw_plan);
#else
    fftwf_destroy_plan(fftw_plan);
#endif
#endif

#if AVFFT
#if MODE == MDCT
    av_mdct_end(avfft);
#elif MODE == R2C || MODE == C2R
    av_rdft_end(avfft);
#else
    av_fft_end(avfft);
#endif
#endif
    av_tx_uninit(&avfftctx);

    av_free(input);
#if AVFFT
    av_free(input_lavc);
    av_free(output_lavc);
#endif
    av_free(output_naive);
#if FFTW
    av_free(input_fftw);
    av_free(output_fftw);
#endif
    av_free(output_new);

    return 0;
}
