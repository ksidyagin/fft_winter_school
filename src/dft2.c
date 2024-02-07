#include "dft.h"

extern void dft2Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    pDst[1].re = pSrc[0].re - pSrc[1].re;
    pDst[1].im = pSrc[0].im - pSrc[1].im;
    pDst[0].re = pSrc[0].re + pSrc[1].re;
    pDst[0].im = pSrc[0].im + pSrc[1].im;
}
