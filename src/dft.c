#include "dft.h"

extern void refDftFwd(const cfloat32_t* pSrc, cfloat32_t* pDst, const int n)
{
    for (int k = 0; k < n; k++)
    {
        float yRe = 0.f;
        float yIm = 0.f;
        for (int j = 0; j < n; j++)
        {
            yRe += pSrc[j].re * cos((OWN_2PI * j * k) / n) + pSrc[j].im * sin((OWN_2PI * j * k) / n);
            yIm +=-pSrc[j].re * sin((OWN_2PI * j * k) / n) + pSrc[j].im * cos((OWN_2PI * j * k) / n);
        }
        pDst[k].re = yRe;
        pDst[k].im = yIm;
    }
}
