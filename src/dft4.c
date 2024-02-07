#include "dft.h"

extern void dft4Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    cfloat32_t tmpDst[4];
    tmpDst[0].re = pSrc[0].re   + pSrc[2].re;
    tmpDst[0].im = pSrc[0].im   + pSrc[2].im;
    tmpDst[1].re = pSrc[0].re   - pSrc[2].re;
    tmpDst[1].im = pSrc[0].im   - pSrc[2].im;
    tmpDst[2].re = pSrc[1].re   + pSrc[3].re;
    tmpDst[2].im = pSrc[1].im   + pSrc[3].im;
    tmpDst[3].re = pSrc[1].re   - pSrc[3].re;
    tmpDst[3].im = pSrc[1].im   - pSrc[3].im;
    pDst[0].re   = tmpDst[0].re + tmpDst[2].re;
    pDst[0].im   = tmpDst[0].im + tmpDst[2].im;
    pDst[1].re   = tmpDst[1].re + tmpDst[3].im;
    pDst[1].im   = tmpDst[1].im - tmpDst[3].re;
    pDst[2].re   = tmpDst[0].re - tmpDst[2].re;
    pDst[2].im   = tmpDst[0].im - tmpDst[2].im;
    pDst[3].re   = tmpDst[1].re - tmpDst[3].im;
    pDst[3].im   = tmpDst[1].im + tmpDst[3].re;
}
