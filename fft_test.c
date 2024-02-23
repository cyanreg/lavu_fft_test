#include <float.h>

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
#define DCT23     4
#define DCT1      5
#define DST1      6

#define MODE FFT
#define INVERSE   0 // Transform direction (only for MDCT and FFT)

/* Double-precision instead of float */
#define DOUBLE    0

#define FFTW      1

#define FFT_LEN   2097152
#define REPS     (1 << 0)
#define IN_PLACE  0
#define NOSHUF    0
#define IMDCT_F   0
#define REAL_TO_REAL 0
#define REAL_TO_IMAGINARY 0
#define SEED      0x8083
#define PRINTOUT  0
#define PRINT_IN  0

#ifndef NO_SIMD
#define NO_SIMD   0
#endif

#if MODE == MDCT
#undef FFTW
#define FFTW 0
#endif

#if FFTW
#include <fftw3.h>
#endif

#if DOUBLE
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

static av_unused void naive_fft(TXComplex *output, TXComplex *input, int len)
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

static av_unused void naive_rdft_r2c(TXComplex *output, TXComplex *input, int len)
{
    TXSample *in = (TXSample *)input;
    double phase = -2*M_PI;

#if REAL_TO_REAL || REAL_TO_IMAGINARY
    TXSample *out = (TXSample *)output;
    int cnt = 0;
#endif

    for (int i = 0; i < (len + 1); i++) {
        TXComplex tmp = { 0 };
        for(int j = 0; j < len*2; j++) {
            TXSample n = in[j];
            TXComplex twiddle = { cos(phase*j*i/(len*2)), sin(phase*j*i/(len*2)) };
            TXComplex res = { n * twiddle.re, n * twiddle.im };
        	tmp = cadd(tmp, res);
        }
#if REAL_TO_REAL
        out[cnt++] = tmp.re;
#elif REAL_TO_IMAGINARY
        out[cnt++ - 1] = tmp.im;
#else
        output[i] = tmp;
#endif
    }
}

static av_unused void naive_rdft_c2r(TXComplex *output, TXComplex *input, int len)
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

static av_unused void naive_idct(TXSample *output, TXSample *input, int n)
{
    for (int i = 0; i < n; i++) {
        double s = 0.5 * input[0];
        for (int k = 1; k < n; k++) {
            double a = M_PI * k * (i + 0.5) / n;
            s += input[k] * cos(a);
        }
        output[i] = 2*s;
    }
}

static av_unused void naive_dct(TXSample *output, TXSample *input, int n)
{
    for (int k = 0; k < n; k++) {
        double s = 0;
        for (int i = 0; i < n; i++) {
            double a = M_PI * k * (i + 0.5) / n;
            s += input[i] * cos(a);
        }
        output[k] = 2*s;
    }
}

static av_unused void naive_dct1(TXSample *output, TXSample *input, int n)
{
    for (int k = 0; k < n; k++) {
        double s = 0;
        for (int i = 1; i < (n - 1); i++) {
            double a = M_PI * k * i / (n - 1);
            s += input[i] * cos(a);
        }
        output[k] = input[0] + ((k % 2) ? -1 : 1)*input[n - 1] + 2*s;
    }
}

static av_unused void naive_dst1(TXSample *output, TXSample *input, int n)
{
    for (int k = 0; k < n; k++) {
        double s = 0;
        for (int i = 0; i < n; i++) {
            double a = M_PI * (k + 1)*(i + 1) / (n + 1);
            s += input[i] * sin(a);
        }
        output[k] = 2*s;
    }
}

static av_unused void naive_imdct(TXSample *dst, TXSample *src, int len)
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

static av_unused void naive_imdct_full(TXSample *dst, TXSample *src, int len)
{
    for (int i = 0; i < len*2; i++) {
        double sum = 0.0;
        for (int j = 0; j < len; j++) {
            int a = (2 * i + 1 + len) * (2 * j + 1);
            sum += src[j] * cos(M_PI * a / (double) (4 * len));
        }
        dst[i] = -sum;
    }
}

static av_unused void naive_mdct(TXSample *dst, TXSample *src, int len)
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

av_tx_fn tx;
void do_avfft_tx(AVTXContext *s, TXComplex *output, TXComplex *_input, int len)
{
    TXComplex *input = _input;
#if IN_PLACE || MODE == C2R || MODE == R2C || MODE == DCT23
    input = av_mallocz(200*(len+(MODE == C2R))*sizeof(*input));
    memcpy(input, _input, (len+(MODE == C2R))*sizeof(*input));
#endif

    int64_t start = av_gettime_relative();

    for (int i = 0; i < REPS; i++) {
        START_TIMER
        tx(s, output, input,
#if MODE == MDCT || MODE == DCT1 || MODE == DST1 || MODE == R2C
           sizeof(TXSample));
#else
           sizeof(TXComplex));
#endif
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
#elif MODE == DCT23
        STOP_TIMER("       av_tx (dct2/3)");
#else
#if INVERSE
        STOP_TIMER("         av_tx (ifft)");
#else
        STOP_TIMER("          av_tx (fft)");
#endif
#endif

#if IN_PLACE || MODE == C2R || MODE == R2C || MODE == DCT23
        memcpy(input, _input, (len+(MODE == C2R))*sizeof(*input));
#endif
    }

#if IN_PLACE || MODE == C2R || MODE == R2C || MODE == DCT23
    av_freep(&input);
#endif

    if (REPS > 1)
        printf("Total for len %i reps %i = %f s\n", len, REPS,
               (av_gettime_relative() - (double)start)/1000000.0);
}

#if FFTW
#if DOUBLE
void do_fftw_tx(fftw_plan fftw_plan, TXComplex *output, TXComplex *_input, int len)
#else
void do_fftw_tx(fftwf_plan fftw_plan, TXComplex *output, TXComplex *_input, int len)
#endif
{
    TXComplex *input = _input;
#if MODE == C2R || MODE == R2C
    input = av_mallocz((len+(MODE == C2R))*sizeof(*input));
    memcpy(input, _input, (len+(MODE == C2R))*sizeof(*input));
#endif

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

#if MODE == C2R || MODE == R2C
    av_freep(&input);
#endif

    if (REPS > 1)
        printf("Total for len %i reps %i = %f s\n", len, REPS,
               (av_gettime_relative() - (double)start)/1000000.0);
}
#endif

#include <libavutil/intreadwrite.h>

void compare_results(TXComplex *ref, TXComplex *dut, int len, const char *str)
{
    int matches = 0;
    double err = 0.0;
    int first_mismatch_idx = -1;
    for (int i = 0; i < len; i++) {
        int nb;
        int m_before = matches;
        if (1) {
            if (((MODE == R2C) && (i == (len - 1)))) {
                matches += (fabs(ref[i].re - dut[i].re) < FLT_EPSILON);
                nb = 1;
            } else {
                matches += (fabs(ref[i].re - dut[i].re) < FLT_EPSILON) +
                           (fabs(ref[i].im - dut[i].im) < FLT_EPSILON);
                nb = 2;
            }
        } else {
            if (((MODE == R2C) && (i == (len - 1)))) {
                for (int j = 0; j < len; j++) {
                    matches += (fabs(ref[j].re - dut[i].re) < FLT_EPSILON);
                    nb = 1;
                }
            } else {
                for (int j = 0; j < len; j++) {
                    matches += (fabs(ref[j].re - dut[i].re) < FLT_EPSILON) +
                               (fabs(ref[j].im - dut[i].im) < FLT_EPSILON);
                    nb = 2;
                }
            }
        }

        if (((m_before == matches) || (m_before + nb != matches)) && (first_mismatch_idx < 0))
            first_mismatch_idx = i;

        err += (ref[i].re - dut[i].re)*(ref[i].re - dut[i].re);
        if (!((MODE == R2C) && (i == (len - 1))))
            err += (ref[i].im - dut[i].im)*(ref[i].im - dut[i].im);
    }

    printf("RMSE %s = %f (%i matches, first mismatch at %i)\n",
           str, sqrtf(err / (len*2)), matches, first_mismatch_idx);
}

void compare_all(TXComplex *all[], int len, const char *names[], int num_arr)
{
    printf("Dump:\n\t\t\t");
    for (int j = 0; j < num_arr; j++) {
        printf("%s%s", names[j], j == (num_arr - 1) ? "\n" : "\t\t\t");
    }

    for (int i = 0; i < len; i++) {
        printf("\t%02i. ", i);
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

    av_log_set_level(AV_LOG_DEBUG);

    if (NO_SIMD)
        av_force_cpu_flags(0);

#if MODE == DCT23 || MODE == DCT1
    printf("Mode: %s %s DCT %s\n",
           DOUBLE ? "double" : "float",
           (MODE == DCT1) ? "DCT-I" : (!INVERSE ? "DCT-II" : "DCT-III"),
           NO_SIMD ? "(no SIMD)" : "(with SIMD)");
#elif MODE == DST1
    printf("Mode: %s %s DST %s\n",
           DOUBLE ? "double" : "float",
           "DST-I",
           NO_SIMD ? "(no SIMD)" : "(with SIMD)");
#elif MODE == C2R || MODE == R2C
    printf("Mode: %s %s DFT %s\n",
           DOUBLE ? "double" : "float",
           MODE == C2R ? "complex-to-real" : "real-to-complex",
           NO_SIMD ? "(no SIMD)" : "(with SIMD)");
#else
    printf("Mode: %s %s %s %s\n",
           DOUBLE ? "double" : "float",
           INVERSE ? "inverse" : "forward",
           (MODE == MDCT ? (IMDCT_F ? "mdct (half)" : "mdct") : "fft"),
           NO_SIMD ? "(no SIMD)" : "(with SIMD)");
#endif

    uint32_t seed = !SEED ? av_get_random_seed() : SEED;
    printf("Seed: 0x%x\n", seed);

    AVTXContext *avfftctx;
    TXSample scale = 1.0;
    int ret = av_tx_init(&avfftctx, &tx,
#if MODE == MDCT
                         DOUBLE ? AV_TX_DOUBLE_MDCT : AV_TX_FLOAT_MDCT,
#elif MODE == R2C || MODE == C2R
                         DOUBLE ? AV_TX_DOUBLE_RDFT : AV_TX_FLOAT_RDFT,
#elif MODE == DCT23
                         DOUBLE ? AV_TX_DOUBLE_DCT : AV_TX_FLOAT_DCT,
#elif MODE == DCT1
                         DOUBLE ? AV_TX_DOUBLE_DCT_I : AV_TX_FLOAT_DCT_I,
#elif MODE == DST1
                         DOUBLE ? AV_TX_DOUBLE_DST_I : AV_TX_FLOAT_DST_I,
#else
                         DOUBLE ? AV_TX_DOUBLE_FFT : AV_TX_FLOAT_FFT,
#endif
#if MODE == DCT23
                         INVERSE,
#else
                         MODE == C2R ? 1 : inv,
#endif
#if MODE == R2C || MODE == C2R || MODE == DCT23 || (MODE == DCT1) || (MODE == DST1)
                         2*tx_len, &scale,
#else
                         tx_len, &scale,
#endif
                         (IN_PLACE ? AV_TX_INPLACE : 0x0) |
                         (IMDCT_F ? AV_TX_FULL_IMDCT : 0x0) |
                         (NO_SIMD ? AV_TX_UNALIGNED : 0x0) |
                         (REAL_TO_REAL ? AV_TX_REAL_TO_REAL : 0x0) |
                         (REAL_TO_IMAGINARY ? AV_TX_REAL_TO_IMAGINARY : 0x0) |
                         (NOSHUF ? (1ULL << 61) : 0x0));
    if (ret) {
        av_log(NULL, AV_LOG_WARNING, "ERROR = %s\n", av_err2str(ret));
        return 1;
    }

    printf("Length: %i\n", tx_len);
    printf("Function: %p\n", tx);

    int num_in = tx_len + (MODE == C2R);
    int num_out = tx_len/(1 + ((MODE == MDCT) && INVERSE && !IMDCT_F)) + (MODE == R2C);

    size_t alloc_in = sizeof(TXComplex)*num_in*20;
    size_t alloc_out = sizeof(TXComplex)*num_out*20;

#if PRINTOUT && PRINT_IN
    alloc_in = FFMAX(alloc_in, alloc_out);
    alloc_out = FFMAX(alloc_in, alloc_out);
#endif

    TXComplex *input = av_mallocz(alloc_in);
    TXComplex *output_new = av_mallocz(alloc_out);
    TXComplex *output_naive = av_mallocz(alloc_out);

#if FFTW
    TXComplex *input_fftw = av_mallocz(alloc_in);
    TXComplex *output_fftw = av_mallocz(alloc_out);

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
#elif MODE == DCT23
#if DOUBLE
    fftw_plan fftw_plan = fftw_plan_r2r_1d  (tx_len*2, IN_PLACE ?
                                                         (double *)output_fftw :
                                                         (double *)input_fftw,
                                                         (double *)output_fftw,
                                                 !INVERSE ? FFTW_REDFT10 : FFTW_REDFT01,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftwf_plan fftw_plan = fftwf_plan_r2r_1d(tx_len*2, IN_PLACE ?
                                                         (float *)output_fftw :
                                                         (float *)input_fftw,
                                                         (float *)output_fftw,
                                                 !INVERSE ? FFTW_REDFT10 : FFTW_REDFT01,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#elif MODE == DCT1
#if DOUBLE
    fftw_plan fftw_plan = fftw_plan_r2r_1d  (tx_len*2, IN_PLACE ?
                                                         (double *)output_fftw :
                                                         (double *)input_fftw,
                                                         (double *)output_fftw,
                                                 FFTW_REDFT00,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftwf_plan fftw_plan = fftwf_plan_r2r_1d(tx_len*2, IN_PLACE ?
                                                         (float *)output_fftw :
                                                         (float *)input_fftw,
                                                         (float *)output_fftw,
                                                 FFTW_REDFT00,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#elif MODE == DST1
#if DOUBLE
    fftw_plan fftw_plan = fftw_plan_r2r_1d  (tx_len*2, IN_PLACE ?
                                                         (double *)output_fftw :
                                                         (double *)input_fftw,
                                                         (double *)output_fftw,
                                                 FFTW_RODFT00,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftwf_plan fftw_plan = fftwf_plan_r2r_1d(tx_len*2, IN_PLACE ?
                                                         (float *)output_fftw :
                                                         (float *)input_fftw,
                                                         (float *)output_fftw,
                                                 FFTW_RODFT00,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#elif MODE == R2C
#if DOUBLE
#if REAL_TO_REAL || REAL_TO_IMAGINARY
    fftw_plan fftw_plan = fftw_plan_r2r_1d(tx_len*2, IN_PLACE ?
                                                         (double *)output_fftw :
                                                         (double *)input_fftw,
                                                         (double *)output_fftw,
                                                 FFTW_R2HC,
                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#else
    fftw_plan fftw_plan = fftw_plan_dft_r2c_1d  (tx_len*2, IN_PLACE ?
                                                         (double *)output_fftw :
                                                         (double *)input_fftw,
                                                 (fftw_complex *)output_fftw,

                                                 FFTW_MEASURE |
                                                 NO_SIMD ? FFTW_UNALIGNED : 0);
#endif
#else
#if REAL_TO_REAL || REAL_TO_IMAGINARY
    fftwf_plan fftw_plan = fftwf_plan_r2r_1d(tx_len*2, IN_PLACE ?
                                                         (float *)output_fftw :
                                                         (float *)input_fftw,
                                                         (float *)output_fftw,
                                                 FFTW_R2HC,
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
    av_log_set_level(AV_LOG_INFO);

    do_avfft_tx(avfftctx, output_new, input, tx_len);
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
#elif (MODE == DCT23 && INVERSE)
    naive_idct((TXSample *)output_naive, (TXSample *)input, 2*tx_len);
#elif (MODE == DCT23 && !INVERSE)
    naive_dct(output_naive, input, 2*tx_len);
#elif (MODE == DCT1)
    naive_dct1(output_naive, input, 2*tx_len);
#elif (MODE == DST1)
    naive_dst1(output_naive, input, 2*tx_len);
#else
    naive_fft(output_naive, input, tx_len);
#endif

#if REAL_TO_REAL || REAL_TO_IMAGINARY
    num_out = 1 + num_out/2;
#endif

#if PRINTOUT
     compare_all(
        (TXComplex *[]){
#if PRINT_IN
            input,
#endif
            output_naive,
            output_new,
#if FFTW
            output_fftw,
#endif
        },
        num_out,
        (const char *[]){
#if PRINT_IN
            "input",
#endif
            "naive",
            "lavu",
#if FFTW
            "fftw",
#endif
        },
        2 + !!FFTW + !!PRINT_IN);
#endif

    compare_results(output_naive, output_new, num_out,      "  av_tx");

#if FFTW
#if DOUBLE
    compare_results(output_naive, output_fftw, num_out,     "  fftw3");
#else
    compare_results(output_naive, output_fftw, num_out,     " fftw3f");
#endif
#endif

#if FFTW
#if DOUBLE
    fftw_destroy_plan(fftw_plan);
#else
    fftwf_destroy_plan(fftw_plan);
#endif
#endif

#if FFTW
    av_free(input_fftw);
    av_free(output_fftw);
#endif

    av_free(input);
    av_free(output_naive);
    av_free(output_new);
    av_tx_uninit(&avfftctx);

    return 0;
}
