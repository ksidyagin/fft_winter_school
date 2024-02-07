#include "dft.h"

void dft32Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    refDftFwd(pSrc, pDst, 32);
}
