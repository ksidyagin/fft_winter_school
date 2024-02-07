#include "dft.h"

extern void dft2Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    // dst:  0.re + 1.re | 0.im + 1.im || 0.re - 1.re | 0.im - 1.im
    // =                                  =
    // src : 0.re        | 0.im        || 0.re        | 0.im
    // +                                  -
    // src : 1.re        | 1.im        || 1.re        | 1.im
    vfloat32m1_t src0 = vle32_v_f32m1((float*)(pSrc + 0), 2);
    vfloat32m1_t src1 = vle32_v_f32m1((float*)(pSrc + 1), 2);
    //   dst1    =   src0     -    src1
    //pDst[1].re = pSrc[0].re - pSrc[1].re;
    //pDst[1].im = pSrc[0].im - pSrc[1].im;
    vfloat32m1_t dst1 = vfsub_vv_f32m1(src0, src1, 2);
    //   dst1    =   src0     +    src1
    //pDst[0].re = pSrc[0].re + pSrc[1].re;
    //pDst[0].im = pSrc[0].im + pSrc[1].im;
    vfloat32m1_t dst0 = vfadd_vv_f32m1(src0, src1, 2);
    vse32_v_f32m1((float*)(pDst + 0), dst0,2);
    vse32_v_f32m1((float*)(pDst + 1), dst1,2);
}
