#include "dft.h"


extern void dft8Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{

    cfloat32_t tmpDst8[8];


// **Шаг 2** - вычисляем $n_1$ $n_2$-точечных DFT. _Напоминаю, что мы их уже написали_

        tmpDst8[1].re = pSrc[0].re - pSrc[4].re;
        tmpDst8[1].im = pSrc[0].im - pSrc[4].im;
        tmpDst8[0].re = pSrc[0].re + pSrc[4].re;
        tmpDst8[0].im = pSrc[0].im + pSrc[4].im;

        tmpDst8[3].re = pSrc[1].re - pSrc[5].re;
        tmpDst8[3].im = pSrc[1].im - pSrc[5].im;
        
        tmpDst8[2].re = pSrc[1].re + pSrc[5].re;
        tmpDst8[2].im = pSrc[1].im + pSrc[5].im;

        tmpDst8[5].re = pSrc[2].im - pSrc[6].im;
        tmpDst8[5].im = -(pSrc[2].re - pSrc[6].re);

        tmpDst8[4].re = pSrc[2].re + pSrc[6].re;
        tmpDst8[4].im = pSrc[2].im + pSrc[6].im;

        tmpDst8[7].re = pSrc[3].re - pSrc[7].re;
        tmpDst8[7].im = pSrc[3].im - pSrc[7].im;
        tmpDst8[6].re = pSrc[3].re + pSrc[7].re;
        tmpDst8[6].im = pSrc[3].im + pSrc[7].im;
// **Шаг 3** - производим поэлементное комплексное умножение на тригонометрические константы. _Упражнение: догадайтесь, почему `i1` начинается с единицы, а не с нуля_

    // for (uint32_t i1 = 1; i1 < n1; i1++)
    // {
    //     for (uint32_t i2 = 0; i2 < n2; i2++)
    //     {
    //         float twdRe = cos(OWN_2PI * i2 * i1 / length);
    //         float twdIm =-sin(OWN_2PI * i2 * i1 / length);
    //         float re = tmpDst8[i1 * n2 + i2].re * twdRe - tmpDst8[i1 * n2 + i2].im * twdIm;
    //         float im = tmpDst8[i1 * n2 + i2].im * twdRe + tmpDst8[i1 * n2 + i2].re * twdIm;
    //         tmpDst8[i1 * n2 + i2].re = re;
    //         tmpDst8[i1 * n2 + i2].im = im;
    //     }
    // }
    float re;
    float im;

            //tmpDst8[2].re = tmpDst8[2].re * 1 - tmpDst8[2].im * 0;
            //tmpDst8[2].im = tmpDst8[2].im * 1 + tmpDst8[2].re * 0;

            re = tmpDst8[3].re * 0.707107 - tmpDst8[3].im * (-0.707107);
            im = tmpDst8[3].im * 0.707107 + tmpDst8[3].re * (-0.707107);
            tmpDst8[3].re = re;
            tmpDst8[3].im = im;


            //tmpDst8[4].re = tmpDst8[4].re * 1 - tmpDst8[4].im * 0;
            //tmpDst8[4].im = tmpDst8[4].im * 1 + tmpDst8[4].re * 0;

            //re = tmpDst8[5].re * 0 - tmpDst8[5].im * (-1);
            //im = tmpDst8[5].im * 0 + tmpDst8[5].re * (-1);
            //tmpDst8[5].re = re;
            //tmpDst8[5].im = im;

            //re = tmpDst8[6].re * 1 - tmpDst8[6].im * 0;
            //im = tmpDst8[6].im * 1 + tmpDst8[6].re * 0;
            //tmpDst8[6].re = re;
            //tmpDst8[6].im = im;

            re = tmpDst8[7].re * (-0.707107) - tmpDst8[7].im * (-0.707107);
            im = tmpDst8[7].im * (-0.707107) + tmpDst8[7].re * (-0.707107);
            tmpDst8[7].re = re;
            tmpDst8[7].im = im;

// **Шаг 5** - вычисляем $n_2$ $n_1$-точечных DFT

    // for (uint32_t i2 = 0; i2 < n2; i2++)
    // {
    //     dft4Fwd(&tmpSrc8[i2 * n1], &tmpDst8[i2 * n1]);

    // }

        pDst[0].re  = tmpDst8[0].re   + tmpDst8[4].re + tmpDst8[2].re   + tmpDst8[6].re;
        pDst[0].im  = tmpDst8[0].im   + tmpDst8[4].im + tmpDst8[2].im   + tmpDst8[6].im;
        pDst[2].re  = (tmpDst8[0].re   - tmpDst8[4].re) + (tmpDst8[2].im   - tmpDst8[6].im);
        pDst[2].im  = (tmpDst8[0].im   - tmpDst8[4].im) - (tmpDst8[2].re   - tmpDst8[6].re);
        pDst[4].re  = (tmpDst8[0].re   + tmpDst8[4].re) - (tmpDst8[2].re   + tmpDst8[6].re);
        pDst[4].im  = (tmpDst8[0].im   + tmpDst8[4].im) - (tmpDst8[2].im   + tmpDst8[6].im);
        pDst[6].re  = (tmpDst8[0].re   - tmpDst8[4].re) - (tmpDst8[2].im   - tmpDst8[6].im);
        pDst[6].im  = (tmpDst8[0].im   - tmpDst8[4].im) + (tmpDst8[2].re   - tmpDst8[6].re);


        pDst[1].re = tmpDst8[1].re   + tmpDst8[5].re + tmpDst8[3].re   + tmpDst8[7].re;
        pDst[1].im = tmpDst8[1].im   + tmpDst8[5].im + tmpDst8[3].im   + tmpDst8[7].im;
        pDst[3].re = (tmpDst8[1].re   - tmpDst8[5].re) + (tmpDst8[3].im   - tmpDst8[7].im);
        pDst[3].im = (tmpDst8[1].im   - tmpDst8[5].im) - (tmpDst8[3].re   - tmpDst8[7].re);
        pDst[5].re = (tmpDst8[1].re   + tmpDst8[5].re) - (tmpDst8[3].re   + tmpDst8[7].re);
        pDst[5].im = (tmpDst8[1].im   + tmpDst8[5].im) - (tmpDst8[3].im   + tmpDst8[7].im);
        pDst[7].re = (tmpDst8[1].re   - tmpDst8[5].re) - (tmpDst8[3].im   - tmpDst8[7].im);
        pDst[7].im = (tmpDst8[1].im   - tmpDst8[5].im) + (tmpDst8[3].re   - tmpDst8[7].re);


    return;
}
