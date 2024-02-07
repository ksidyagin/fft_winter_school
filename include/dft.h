#ifndef DFT_H
#define DFT_H

#include "common_defs.h"

extern void refDftFwd(const cfloat32_t* pSrc, cfloat32_t* pDst, const int len);

extern void dft2Fwd (const cfloat32_t *pSrc, cfloat32_t *pDst);
extern void dft4Fwd (const cfloat32_t *pSrc, cfloat32_t *pDst);
extern void dft8Fwd (const cfloat32_t *pSrc, cfloat32_t *pDst);
extern void dft16Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst);
extern void dft32Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst);

#endif // DFT_H
