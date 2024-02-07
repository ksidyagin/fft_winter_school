#include "dft.h"

extern void dft16Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    cfloat32_t tmpSrc[16];
    cfloat32_t tmpDst[16];
    //n1=8
    //n2=2

// **Шаг 1** - транспонируем входные данные по столбцам, по которым будем вычислять $n_2$-точечные DFT

// **Шаг 2** - вычисляем $n_1$ $n_2$-точечных DFT. _Напоминаю, что мы их уже написали_

    // for (uint32_t i1 = 0; i1 < 8; i1++)
    // {
    //     dft2Fwd(&tmpSrc[i1 * 2], &tmpDst[i1 * 2]);

    // }
        tmpDst[1].re = pSrc[0].re - pSrc[8].re;
        tmpDst[1].im = pSrc[0].im - pSrc[8].im;
        tmpDst[0].re = pSrc[0].re + pSrc[8].re;
        tmpDst[0].im = pSrc[0].im + pSrc[8].im;

        tmpDst[3].re = pSrc[1].re - pSrc[9].re;
        tmpDst[3].im = pSrc[1].im - pSrc[9].im;
        tmpDst[2].re = pSrc[1].re + pSrc[9].re;
        tmpDst[2].im = pSrc[1].im + pSrc[9].im;

        tmpDst[5].re = pSrc[2].re - pSrc[10].re;
        tmpDst[5].im = pSrc[2].im - pSrc[10].im;
        tmpDst[4].re = pSrc[2].re + pSrc[10].re;
        tmpDst[4].im = pSrc[2].im + pSrc[10].im;

        tmpDst[7].re = pSrc[3].re - pSrc[11].re;
        tmpDst[7].im = pSrc[3].im - pSrc[11].im;
        tmpDst[6].re = pSrc[3].re + pSrc[11].re;
        tmpDst[6].im = pSrc[3].im + pSrc[11].im;

        tmpDst[9].re = pSrc[4].re - pSrc[12].re;
        tmpDst[9].im = pSrc[4].im - pSrc[12].im;
        tmpDst[8].re = pSrc[4].re + pSrc[12].re;
        tmpDst[8].im = pSrc[4].im + pSrc[12].im;

        tmpDst[11].re = pSrc[5].re - pSrc[13].re;
        tmpDst[11].im = pSrc[5].im - pSrc[13].im;
        tmpDst[10].re = pSrc[5].re + pSrc[13].re;
        tmpDst[10].im = pSrc[5].im + pSrc[13].im;

        tmpDst[13].re = pSrc[6].re - pSrc[14].re;
        tmpDst[13].im = pSrc[6].im - pSrc[14].im;
        tmpDst[12].re = pSrc[6].re + pSrc[14].re;
        tmpDst[12].im = pSrc[6].im + pSrc[14].im;

        tmpDst[15].re = pSrc[7].re - pSrc[15].re;
        tmpDst[15].im = pSrc[7].im - pSrc[15].im;
        tmpDst[14].re = pSrc[7].re + pSrc[15].re;
        tmpDst[14].im = pSrc[7].im + pSrc[15].im;


// **Шаг 3** - производим поэлементное комплексное умножение на тригонометрические константы. _Упражнение: догадайтесь, почему `i1` начинается с единицы, а не с нуля_

    float twdRe;
    float twdIm;
    float re;
    float im;

    //i1=1        
            // re = tmpDst[1 * 2 + 0].re * 1 - tmpDst[1 * 2 + 0].im * 0;
            // im = tmpDst[1 * 2 + 0].im * 1 + tmpDst[1 * 2 + 0].re * 0;
            // tmpDst[1 * 2 + 0].re = re;
            // tmpDst[1 * 2 + 0].im = im;

            re = tmpDst[1 * 2 + 1].re * 0.92388 - tmpDst[1 * 2 + 1].im * (-0.382683);
            im = tmpDst[1 * 2 + 1].im * 0.92388 + tmpDst[1 * 2 + 1].re * (-0.382683);
            tmpDst[1 * 2 + 1].re = re;
            tmpDst[1 * 2 + 1].im = im;

    //i1=2


            re = tmpDst[2 * 2 + 1].re * 0.707107 - tmpDst[2 * 2 + 1].im * (-0.707107);
            im = tmpDst[2 * 2 + 1].im * 0.707107 + tmpDst[2 * 2 + 1].re * (-0.707107);
            tmpDst[2 * 2 + 1].re = re;
            tmpDst[2 * 2 + 1].im = im;


    //i1=3


            re = tmpDst[3 * 2 + 1].re * 0.382683 - tmpDst[3 * 2 + 1].im * (-0.92388);
            im = tmpDst[3 * 2 + 1].im * 0.382683 + tmpDst[3 * 2 + 1].re * (-0.92388);
            tmpDst[3 * 2 + 1].re = re;
            tmpDst[3 * 2 + 1].im = im;   

    //i1=4

            re = tmpDst[4 * 2 + 1].re * 0 - tmpDst[4 * 2 + 1].im * (-1);
            im = tmpDst[4 * 2 + 1].im * 0 + tmpDst[4 * 2 + 1].re * (-1);
            tmpDst[4 * 2 + 1].re = re;
            tmpDst[4 * 2 + 1].im = im;     

    //i1=5

            re = tmpDst[5 * 2 + 1].re * (-0.382683) - tmpDst[5 * 2 + 1].im * (-0.92388);
            im = tmpDst[5 * 2 + 1].im * (-0.382683) + tmpDst[5 * 2 + 1].re * (-0.92388);
            tmpDst[5 * 2 + 1].re = re;
            tmpDst[5 * 2 + 1].im = im;  
    
    //i1=6

            re = tmpDst[6 * 2 + 1].re * (-0.707107) - tmpDst[6 * 2 + 1].im * (-0.707107);
            im = tmpDst[6 * 2 + 1].im * (-0.707107) + tmpDst[6 * 2 + 1].re * (-0.707107);
            tmpDst[6 * 2 + 1].re = re;
            tmpDst[6 * 2 + 1].im = im;  

            
    //i1=7
;
            re = tmpDst[7 * 2 + 1].re * (-0.92388) - tmpDst[7 * 2 + 1].im * (-0.382683);
            im = tmpDst[7 * 2 + 1].im * (-0.92388) + tmpDst[7 * 2 + 1].re * (-0.382683);
            tmpDst[7 * 2 + 1].re = re;
            tmpDst[7 * 2 + 1].im = im;  


    //i1=8   

            re = tmpDst[8 * 2 + 1].re * (-1) - tmpDst[8 * 2 + 1].im * 0;
            im = tmpDst[8 * 2 + 1].im * (-1) + tmpDst[8 * 2 + 1].re * 0;
            tmpDst[8 * 2 + 1].re = re;
            tmpDst[8 * 2 + 1].im = im;   


    //i1=9

            twdRe = cos(OWN_2PI * 1 * 9 / 16);
            twdIm =-sin(OWN_2PI * 1 * 9 / 16);
            re = tmpDst[9 * 2 + 1].re * twdRe - tmpDst[9 * 2 + 1].im * twdIm;
            im = tmpDst[9 * 2 + 1].im * twdRe + tmpDst[9 * 2 + 1].re * twdIm;
            tmpDst[9 * 2 + 1].re = re;
            tmpDst[9 * 2 + 1].im = im;  

    //i1=10

            twdRe = cos(OWN_2PI * 1 * 9 / 16);
            twdIm =-sin(OWN_2PI * 1 * 9 / 16);
            re = tmpDst[10 * 2 + 1].re * twdRe - tmpDst[10 * 2 + 1].im * twdIm;
            im = tmpDst[10 * 2 + 1].im * twdRe + tmpDst[10 * 2 + 1].re * twdIm;
            tmpDst[10 * 2 + 1].re = re;
            tmpDst[10 * 2 + 1].im = im;  

    //i1=11

            twdRe = cos(OWN_2PI * 1 * 9 / 16);
            twdIm =-sin(OWN_2PI * 1 * 9 / 16);
            re = tmpDst[11 * 2 + 1].re * twdRe - tmpDst[11 * 2 + 1].im * twdIm;
            im = tmpDst[11 * 2 + 1].im * twdRe + tmpDst[11 * 2 + 1].re * twdIm;
            tmpDst[11 * 2 + 1].re = re;
            tmpDst[11 * 2 + 1].im = im; 

    //i1=12

            twdRe = cos(OWN_2PI * 1 * 9 / 16);
            twdIm =-sin(OWN_2PI * 1 * 9 / 16);
            re = tmpDst[12 * 2 + 1].re * twdRe - tmpDst[12 * 2 + 1].im * twdIm;
            im = tmpDst[12 * 2 + 1].im * twdRe + tmpDst[12 * 2 + 1].re * twdIm;
            tmpDst[12 * 2 + 1].re = re;
            tmpDst[12 * 2 + 1].im = im; 

    //i1=13

            twdRe = cos(OWN_2PI * 1 * 9 / 16);
            twdIm =-sin(OWN_2PI * 1 * 9 / 16);
            re = tmpDst[13 * 2 + 1].re * twdRe - tmpDst[13 * 2 + 1].im * twdIm;
            im = tmpDst[13 * 2 + 1].im * twdRe + tmpDst[13 * 2 + 1].re * twdIm;
            tmpDst[13 * 2 + 1].re = re;
            tmpDst[13 * 2 + 1].im = im; 

    //i1=14

            twdRe = cos(OWN_2PI * 1 * 9 / 16);
            twdIm =-sin(OWN_2PI * 1 * 9 / 16);
            re = tmpDst[14 * 2 + 1].re * twdRe - tmpDst[14 * 2 + 1].im * twdIm;
            im = tmpDst[14 * 2 + 1].im * twdRe + tmpDst[14 * 2 + 1].re * twdIm;
            tmpDst[14 * 2 + 1].re = re;
            tmpDst[14 * 2 + 1].im = im; 

    //i1=15

            twdRe = cos(OWN_2PI * 1 * 9 / 16);
            twdIm =-sin(OWN_2PI * 1 * 9 / 16);
            re = tmpDst[15 * 2 + 1].re * twdRe - tmpDst[15 * 2 + 1].im * twdIm;
            im = tmpDst[15 * 2 + 1].im * twdRe + tmpDst[15 * 2 + 1].re * twdIm;
            tmpDst[15 * 2 + 1].re = re;
            tmpDst[15 * 2 + 1].im = im; 
// **Шаг 4** - транспонируем промежуточные результаты по строкам для $n_1$-точечных DFT

    // for (uint32_t i1 = 0; i1 < 8; i1++)
    // {
    //     for (uint32_t i2 = 0; i2 < 2; i2++)
    //     {
    //         tmpSrc[8 * i2 + i1].re = tmpDst[2 * i1 + i2].re;
    //         tmpSrc[8 * i2 + i1].im = tmpDst[2 * i1 + i2].im;
    //     }
    // }

        tmpSrc[0].re = tmpDst[0].re;
        tmpSrc[0].im = tmpDst[0].im;
        tmpSrc[8].re = tmpDst[1].re;
        tmpSrc[8].im = tmpDst[1].im;

        tmpSrc[1].re = tmpDst[2].re;
        tmpSrc[1].im = tmpDst[2].im;
        tmpSrc[9].re = tmpDst[3].re;
        tmpSrc[9].im = tmpDst[3].im;

        tmpSrc[2].re = tmpDst[4].re;
        tmpSrc[2].im = tmpDst[4].im;
        tmpSrc[10].re = tmpDst[5].re;
        tmpSrc[10].im = tmpDst[5].im;

        tmpSrc[3].re = tmpDst[6].re;
        tmpSrc[3].im = tmpDst[6].im;
        tmpSrc[11].re = tmpDst[7].re;
        tmpSrc[11].im = tmpDst[7].im;

        tmpSrc[4].re = tmpDst[8].re;
        tmpSrc[4].im = tmpDst[8].im;
        tmpSrc[12].re = tmpDst[9].re;
        tmpSrc[12].im = tmpDst[9].im;

        tmpSrc[5].re = tmpDst[10].re;
        tmpSrc[5].im = tmpDst[10].im;
        tmpSrc[13].re = tmpDst[11].re;
        tmpSrc[13].im = tmpDst[11].im;

        tmpSrc[6].re = tmpDst[12].re;
        tmpSrc[6].im = tmpDst[12].im;
        tmpSrc[14].re = tmpDst[13].re;
        tmpSrc[14].im = tmpDst[13].im;

        tmpSrc[7].re = tmpDst[14].re;
        tmpSrc[7].im = tmpDst[14].im;
        tmpSrc[15].re = tmpDst[15].re;
        tmpSrc[15].im = tmpDst[15].im;

// **Шаг 5** - вычисляем $n_2$ $n_1$-точечных DFT

    // for (uint32_t i2 = 0; i2 < 2; i2++)
    // {
    //     dft8Fwd(&tmpSrc[i2 * 8], &tmpDst[i2 * 8]);
    // }

    //i2=0
            
    cfloat32_t tmpDst8[16];

        tmpDst8[1].re = tmpSrc[0].re - tmpSrc[4].re;
        tmpDst8[1].im = tmpSrc[0].im - tmpSrc[4].im;
        tmpDst8[0].re = tmpSrc[0].re + tmpSrc[4].re;
        tmpDst8[0].im = tmpSrc[0].im + tmpSrc[4].im;

        tmpDst8[3].re = tmpSrc[1].re - tmpSrc[5].re;
        tmpDst8[3].im = tmpSrc[1].im - tmpSrc[5].im;
        
        tmpDst8[2].re = tmpSrc[1].re + tmpSrc[5].re;
        tmpDst8[2].im = tmpSrc[1].im + tmpSrc[5].im;

        tmpDst8[5].re = tmpSrc[2].im - tmpSrc[6].im;
        tmpDst8[5].im = -(tmpSrc[2].re - tmpSrc[6].re);

        tmpDst8[4].re = tmpSrc[2].re + tmpSrc[6].re;
        tmpDst8[4].im = tmpSrc[2].im + tmpSrc[6].im;

        tmpDst8[7].re = tmpSrc[3].re - tmpSrc[7].re;
        tmpDst8[7].im = tmpSrc[3].im - tmpSrc[7].im;
        tmpDst8[6].re = tmpSrc[3].re + tmpSrc[7].re;
        tmpDst8[6].im = tmpSrc[3].im + tmpSrc[7].im;

            re = tmpDst8[3].re * 0.707107 - tmpDst8[3].im * (-0.707107);
            im = tmpDst8[3].im * 0.707107 + tmpDst8[3].re * (-0.707107);
            tmpDst8[3].re = re;
            tmpDst8[3].im = im;

            re = tmpDst8[7].re * (-0.707107) - tmpDst8[7].im * (-0.707107);
            im = tmpDst8[7].im * (-0.707107) + tmpDst8[7].re * (-0.707107);
            tmpDst8[7].re = re;
            tmpDst8[7].im = im;


        pDst[0].re  = tmpDst8[0].re   + tmpDst8[4].re + tmpDst8[2].re   + tmpDst8[6].re;
        pDst[0].im  = tmpDst8[0].im   + tmpDst8[4].im + tmpDst8[2].im   + tmpDst8[6].im;
        pDst[4].re  = (tmpDst8[0].re   - tmpDst8[4].re) + (tmpDst8[2].im   - tmpDst8[6].im);
        pDst[4].im  = (tmpDst8[0].im   - tmpDst8[4].im) - (tmpDst8[2].re   - tmpDst8[6].re);
        pDst[8].re  = (tmpDst8[0].re   + tmpDst8[4].re) - (tmpDst8[2].re   + tmpDst8[6].re);
        pDst[8].im  = (tmpDst8[0].im   + tmpDst8[4].im) - (tmpDst8[2].im   + tmpDst8[6].im);
        pDst[12].re  = (tmpDst8[0].re   - tmpDst8[4].re) - (tmpDst8[2].im   - tmpDst8[6].im);
        pDst[12].im  = (tmpDst8[0].im   - tmpDst8[4].im) + (tmpDst8[2].re   - tmpDst8[6].re);


        pDst[2].re = tmpDst8[1].re   + tmpDst8[5].re + tmpDst8[3].re   + tmpDst8[7].re;
        pDst[2].im = tmpDst8[1].im   + tmpDst8[5].im + tmpDst8[3].im   + tmpDst8[7].im;
        pDst[6].re = (tmpDst8[1].re   - tmpDst8[5].re) + (tmpDst8[3].im   - tmpDst8[7].im);
        pDst[6].im = (tmpDst8[1].im   - tmpDst8[5].im) - (tmpDst8[3].re   - tmpDst8[7].re);
        pDst[10].re = (tmpDst8[1].re   + tmpDst8[5].re) - (tmpDst8[3].re   + tmpDst8[7].re);
        pDst[10].im = (tmpDst8[1].im   + tmpDst8[5].im) - (tmpDst8[3].im   + tmpDst8[7].im);
        pDst[14].re = (tmpDst8[1].re   - tmpDst8[5].re) - (tmpDst8[3].im   - tmpDst8[7].im);
        pDst[14].im = (tmpDst8[1].im   - tmpDst8[5].im) + (tmpDst8[3].re   - tmpDst8[7].re);



    //i2=1


        tmpDst8[9].re = tmpSrc[8].re - tmpSrc[12].re;
        tmpDst8[9].im = tmpSrc[8].im - tmpSrc[12].im;
        tmpDst8[8].re = tmpSrc[8].re + tmpSrc[12].re;
        tmpDst8[8].im = tmpSrc[8].im + tmpSrc[12].im;

        tmpDst8[11].re = tmpSrc[9].re - tmpSrc[13].re;
        tmpDst8[11].im = tmpSrc[9].im - tmpSrc[13].im;
        
        tmpDst8[10].re = tmpSrc[9].re + tmpSrc[13].re;
        tmpDst8[10].im = tmpSrc[9].im + tmpSrc[13].im;

        tmpDst8[13].re = tmpSrc[10].im - tmpSrc[14].im;
        tmpDst8[13].im = -(tmpSrc[10].re - tmpSrc[14].re);

        tmpDst8[12].re = tmpSrc[10].re + tmpSrc[14].re;
        tmpDst8[12].im = tmpSrc[10].im + tmpSrc[14].im;

        tmpDst8[15].re = tmpSrc[11].re - tmpSrc[15].re;
        tmpDst8[15].im = tmpSrc[11].im - tmpSrc[15].im;
        tmpDst8[14].re = tmpSrc[11].re + tmpSrc[15].re;
        tmpDst8[14].im = tmpSrc[11].im + tmpSrc[15].im;


            re = tmpDst8[11].re * 0.707107 - tmpDst8[11].im * (-0.707107);
            im = tmpDst8[11].im * 0.707107 + tmpDst8[11].re * (-0.707107);
            tmpDst8[11].re = re;
            tmpDst8[11].im = im;

            re = tmpDst8[15].re * (-0.707107) - tmpDst8[15].im * (-0.707107);
            im = tmpDst8[15].im * (-0.707107) + tmpDst8[15].re * (-0.707107);
            tmpDst8[15].re = re;
            tmpDst8[15].im = im;


        pDst[1].re  = tmpDst8[8].re   + tmpDst8[12].re + tmpDst8[10].re   + tmpDst8[14].re;
        pDst[1].im  = tmpDst8[8].im   + tmpDst8[12].im + tmpDst8[10].im   + tmpDst8[14].im;
        pDst[5].re  = (tmpDst8[8].re   - tmpDst8[12].re) + (tmpDst8[10].im   - tmpDst8[14].im);
        pDst[5].im  = (tmpDst8[8].im   - tmpDst8[12].im) - (tmpDst8[10].re   - tmpDst8[14].re);
        pDst[9].re  = (tmpDst8[8].re   + tmpDst8[12].re) - (tmpDst8[10].re   + tmpDst8[14].re);
        pDst[9].im  = (tmpDst8[8].im   + tmpDst8[12].im) - (tmpDst8[10].im   + tmpDst8[14].im);
        pDst[13].re  = (tmpDst8[8].re   - tmpDst8[12].re) - (tmpDst8[10].im   - tmpDst8[14].im);
        pDst[13].im  = (tmpDst8[8].im   - tmpDst8[12].im) + (tmpDst8[10].re   - tmpDst8[14].re);


        pDst[3].re = tmpDst8[9].re   + tmpDst8[13].re + tmpDst8[11].re   + tmpDst8[15].re;
        pDst[3].im = tmpDst8[9].im   + tmpDst8[13].im + tmpDst8[11].im   + tmpDst8[15].im;
        pDst[7].re = (tmpDst8[9].re   - tmpDst8[13].re) + (tmpDst8[11].im   - tmpDst8[15].im);
        pDst[7].im = (tmpDst8[9].im   - tmpDst8[13].im) - (tmpDst8[11].re   - tmpDst8[15].re);
        pDst[11].re = (tmpDst8[9].re   + tmpDst8[13].re) - (tmpDst8[11].re   + tmpDst8[15].re);
        pDst[11].im = (tmpDst8[9].im   + tmpDst8[13].im) - (tmpDst8[11].im   + tmpDst8[15].im);
        pDst[15].re = (tmpDst8[9].re   - tmpDst8[13].re) - (tmpDst8[11].im   - tmpDst8[15].im);
        pDst[15].im = (tmpDst8[9].im   - tmpDst8[13].im) + (tmpDst8[11].re   - tmpDst8[15].re);




// **Шаг 6** - правильно расставляем результаты вычислений.

    // for (uint32_t i2 = 0; i2 < 2; i2++)
    // {
    //     for (uint32_t i1 = 0; i1 < 8; i1++)
    //     {
    //         pDst[2 * i1 + i2].re = tmpDst[i1 + i2 * 8].re;
    //         pDst[2 * i1 + i2].im = tmpDst[i1 + i2 * 8].im;
    //     }
    // }
    // //i2=0
    //         pDst[0].re = tmpDst[0].re;
    //         pDst[0].im = tmpDst[0].im;

    //         pDst[2].re = tmpDst[1].re;
    //         pDst[2].im = tmpDst[1].im;

    //         pDst[4].re = tmpDst[2].re;
    //         pDst[4].im = tmpDst[2].im;

    //         pDst[6].re = tmpDst[3].re;
    //         pDst[6].im = tmpDst[3].im;

    //         pDst[8].re = tmpDst[4].re;
    //         pDst[8].im = tmpDst[4].im;

    //         pDst[10].re = tmpDst[5].re;
    //         pDst[10].im = tmpDst[5].im;

    //         pDst[12].re = tmpDst[6].re;
    //         pDst[12].im = tmpDst[6].im;

    //         pDst[14].re = tmpDst[7].re;
    //         pDst[14].im = tmpDst[7].im;

    // //i2=1
    //         pDst[1].re = tmpDst[8].re;
    //         pDst[1].im = tmpDst[8].im;

    //         pDst[3].re = tmpDst[9].re;
    //         pDst[3].im = tmpDst[9].im;

    //         pDst[5].re = tmpDst[10].re;
    //         pDst[5].im = tmpDst[10].im;

    //         pDst[7].re = tmpDst[11].re;
    //         pDst[7].im = tmpDst[11].im;

    //         pDst[9].re = tmpDst[12].re;
    //         pDst[9].im = tmpDst[12].im;

    //         pDst[11].re = tmpDst[13].re;
    //         pDst[11].im = tmpDst[13].im;

    //         pDst[13].re = tmpDst[14].re;
    //         pDst[13].im = tmpDst[14].im;

    //         pDst[15].re = tmpDst[15].re;
    //         pDst[15].im = tmpDst[15].im;
    return;
}
