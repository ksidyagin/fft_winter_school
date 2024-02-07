#include "dft.h"

extern void dft16Fwd(const cfloat32_t *pSrc, cfloat32_t *pDst)
{
    //cfloat32_t tmpSrc[16];
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

    float re;
    float im;

    //i1=1        
            // re = tmpDst[1 * 2 + 0].re * 1 - tmpDst[1 * 2 + 0].im * 0;
            // im = tmpDst[1 * 2 + 0].im * 1 + tmpDst[1 * 2 + 0].re * 0;
            // tmpDst[1 * 2 + 0].re = re;
            // tmpDst[1 * 2 + 0].im = im;

            re = tmpDst[3].re * 0.92388 - tmpDst[3].im * (-0.382683);
            im = tmpDst[3].im * 0.92388 + tmpDst[3].re * (-0.382683);
            tmpDst[3].re = re;
            tmpDst[3].im = im;

    //i1=2


            re = tmpDst[5].re * 0.707107 - tmpDst[5].im * (-0.707107);
            im = tmpDst[5].im * 0.707107 + tmpDst[5].re * (-0.707107);
            tmpDst[5].re = re;
            tmpDst[5].im = im;


    //i1=3


            re = tmpDst[7].re * 0.382683 - tmpDst[7].im * (-0.92388);
            im = tmpDst[7].im * 0.382683 + tmpDst[7].re * (-0.92388);
            tmpDst[7].re = re;
            tmpDst[7].im = im;   

    //i1=4

            re = tmpDst[9].re * 0 - tmpDst[9].im * (-1);
            im = tmpDst[9].im * 0 + tmpDst[9].re * (-1);
            tmpDst[9].re = re;
            tmpDst[9].im = im;     

    //i1=5

            re = tmpDst[11].re * (-0.382683) - tmpDst[11].im * (-0.92388);
            im = tmpDst[11].im * (-0.382683) + tmpDst[11].re * (-0.92388);
            tmpDst[11].re = re;
            tmpDst[11].im = im;  
    
    //i1=6

            re = tmpDst[13].re * (-0.707107) - tmpDst[13].im * (-0.707107);
            im = tmpDst[13].im * (-0.707107) + tmpDst[13].re * (-0.707107);
            tmpDst[13].re = re;
            tmpDst[13].im = im;  

            
    //i1=7

            re = tmpDst[15].re * (-0.92388) - tmpDst[15].im * (-0.382683);
            im = tmpDst[15].im * (-0.92388) + tmpDst[15].re * (-0.382683);
            tmpDst[15].re = re;
            tmpDst[15].im = im;  



// **Шаг 4** - транспонируем промежуточные результаты по строкам для $n_1$-точечных DFT

    // for (uint32_t i1 = 0; i1 < 8; i1++)
    // {
    //     for (uint32_t i2 = 0; i2 < 2; i2++)
    //     {
    //         tmpSrc[8 * i2 + i1].re = tmpDst[2 * i1 + i2].re;
    //         tmpSrc[8 * i2 + i1].im = tmpDst[2 * i1 + i2].im;
    //     }
    // }



// **Шаг 5** - вычисляем $n_2$ $n_1$-точечных DFT

    // for (uint32_t i2 = 0; i2 < 2; i2++)
    // {
    //     dft8Fwd(&tmpSrc[i2 * 8], &tmpDst[i2 * 8]);
    // }

    //i2=0
            
    cfloat32_t tmpDst8[16];

        tmpDst8[1].re = tmpDst[0].re - tmpDst[8].re;
        tmpDst8[1].im = tmpDst[0].im - tmpDst[8].im;
        tmpDst8[0].re = tmpDst[0].re + tmpDst[8].re;
        tmpDst8[0].im = tmpDst[0].im + tmpDst[8].im;

        tmpDst8[3].re = tmpDst[2].re - tmpDst[10].re;
        tmpDst8[3].im = tmpDst[2].im - tmpDst[10].im;
        
        tmpDst8[2].re = tmpDst[2].re + tmpDst[10].re;
        tmpDst8[2].im = tmpDst[2].im + tmpDst[10].im;

        tmpDst8[5].re = tmpDst[4].im - tmpDst[12].im;
        tmpDst8[5].im = -(tmpDst[4].re - tmpDst[12].re);

        tmpDst8[4].re = tmpDst[4].re + tmpDst[12].re;
        tmpDst8[4].im = tmpDst[4].im + tmpDst[12].im;

        tmpDst8[7].re = tmpDst[6].re - tmpDst[14].re;
        tmpDst8[7].im = tmpDst[6].im - tmpDst[14].im;
        tmpDst8[6].re = tmpDst[6].re + tmpDst[14].re;
        tmpDst8[6].im = tmpDst[6].im + tmpDst[14].im;

            re = tmpDst8[3].re * 0.707107;
            im= tmpDst8[3].im * 0.707107;

            tmpDst8[3].re = re + im;
            tmpDst8[3].im = im - re;


            re = - tmpDst8[7].re * 0.707107;
            im =   tmpDst8[7].im * 0.707107;
            tmpDst8[7].re = re+im;
            tmpDst8[7].im = re-im;
        
        cfloat32_t tmpRes8[16];
        tmpRes8[0].re = tmpDst8[0].re   + tmpDst8[4].re;
        tmpRes8[0].im = tmpDst8[0].im   + tmpDst8[4].im;
        tmpRes8[1].re = tmpDst8[2].re   + tmpDst8[6].re;
        tmpRes8[1].im = tmpDst8[2].im   + tmpDst8[6].im;
        tmpRes8[2].re = tmpDst8[1].re   + tmpDst8[5].re;
        tmpRes8[2].im = tmpDst8[1].im   + tmpDst8[5].im;
        tmpRes8[3].re = tmpDst8[3].re   + tmpDst8[7].re;
        tmpRes8[3].im = tmpDst8[3].im   + tmpDst8[7].im;
        tmpRes8[4].re = tmpDst8[0].re   - tmpDst8[4].re;
        tmpRes8[4].im = tmpDst8[0].im   - tmpDst8[4].im;
        tmpRes8[5].re = tmpDst8[2].re   - tmpDst8[6].re;
        tmpRes8[5].im = tmpDst8[2].im   - tmpDst8[6].im;
        tmpRes8[6].re = tmpDst8[1].re   - tmpDst8[5].re;
        tmpRes8[6].im = tmpDst8[1].im   - tmpDst8[5].im;
        tmpRes8[7].re = tmpDst8[3].re   - tmpDst8[7].re;
        tmpRes8[7].im = tmpDst8[3].im   - tmpDst8[7].im;
        
        pDst[0].re  = tmpRes8[0].re + tmpRes8[1].re;
        pDst[0].im  = tmpRes8[0].im + tmpRes8[1].im;
        pDst[4].re  = tmpRes8[4].re + tmpRes8[5].im;
        pDst[4].im  = tmpRes8[4].im - tmpRes8[5].re;
        pDst[8].re  = tmpRes8[0].re - tmpRes8[1].re;
        pDst[8].im  = tmpRes8[0].im - tmpRes8[1].im;
        pDst[12].re  = tmpRes8[4].re - tmpRes8[5].im;
        pDst[12].im  = tmpRes8[4].im + tmpRes8[5].re;


        pDst[2].re = tmpRes8[2].re + tmpRes8[3].re;
        pDst[2].im = tmpRes8[2].im + tmpRes8[3].im;
        pDst[6].re = tmpRes8[6].re + tmpRes8[7].im;
        pDst[6].im = tmpRes8[6].im - tmpRes8[7].re;
        pDst[10].re = tmpRes8[2].re - tmpRes8[3].re;
        pDst[10].im = tmpRes8[2].im - tmpRes8[3].im;
        pDst[14].re = tmpRes8[6].re - tmpRes8[7].im;
        pDst[14].im = tmpRes8[6].im + tmpRes8[7].re;


    //i2=1


        tmpDst8[9].re = tmpDst[1].re - tmpDst[9].re;
        tmpDst8[9].im = tmpDst[1].im - tmpDst[9].im;
        tmpDst8[8].re = tmpDst[1].re + tmpDst[9].re;
        tmpDst8[8].im = tmpDst[1].im + tmpDst[9].im;

        tmpDst8[11].re = tmpDst[3].re - tmpDst[11].re;
        tmpDst8[11].im = tmpDst[3].im - tmpDst[11].im;
        
        tmpDst8[10].re = tmpDst[3].re + tmpDst[11].re;
        tmpDst8[10].im = tmpDst[3].im + tmpDst[11].im;

        tmpDst8[13].re = tmpDst[5].im - tmpDst[13].im;
        tmpDst8[13].im = -(tmpDst[5].re - tmpDst[13].re);

        tmpDst8[12].re = tmpDst[5].re + tmpDst[13].re;
        tmpDst8[12].im = tmpDst[5].im + tmpDst[13].im;

        tmpDst8[15].re = tmpDst[7].re - tmpDst[15].re;
        tmpDst8[15].im = tmpDst[7].im - tmpDst[15].im;
        tmpDst8[14].re = tmpDst[7].re + tmpDst[15].re;
        tmpDst8[14].im = tmpDst[7].im + tmpDst[15].im;


            re = tmpDst8[11].re * 0.707107;
            im = tmpDst8[11].im * 0.707107;
            tmpDst8[11].re = re+im;
            tmpDst8[11].im = im-re;

            re = -tmpDst8[15].re * 0.707107;
            im = tmpDst8[15].im * 0.707107;
            tmpDst8[15].re = re+im;
            tmpDst8[15].im = re-im;


        tmpRes8[8].re = tmpDst8[8].re   + tmpDst8[12].re;
        tmpRes8[8].im = tmpDst8[8].im   + tmpDst8[12].im;
        tmpRes8[9].re = tmpDst8[10].re   + tmpDst8[14].re;
        tmpRes8[9].im = tmpDst8[10].im   + tmpDst8[14].im;
        tmpRes8[10].re = tmpDst8[9].re   + tmpDst8[13].re;
        tmpRes8[10].im = tmpDst8[9].im   + tmpDst8[13].im;
        tmpRes8[11].re = tmpDst8[11].re   + tmpDst8[15].re;
        tmpRes8[11].im = tmpDst8[11].im   + tmpDst8[15].im;
        tmpRes8[12].re = tmpDst8[8].re   - tmpDst8[12].re;
        tmpRes8[12].im = tmpDst8[8].im   - tmpDst8[12].im;
        tmpRes8[13].re = tmpDst8[10].re   - tmpDst8[14].re;
        tmpRes8[13].im = tmpDst8[10].im   - tmpDst8[14].im;
        tmpRes8[14].re = tmpDst8[9].re   - tmpDst8[13].re;
        tmpRes8[14].im = tmpDst8[9].im   - tmpDst8[13].im;
        tmpRes8[15].re = tmpDst8[11].re   - tmpDst8[15].re;
        tmpRes8[15].im = tmpDst8[11].im   - tmpDst8[15].im;

        pDst[1].re  = tmpRes8[8].re + tmpRes8[9].re;
        pDst[1].im  = tmpRes8[8].im + tmpRes8[9].im;
        pDst[5].re  = tmpRes8[12].re + tmpRes8[13].im;
        pDst[5].im  = tmpRes8[12].im - tmpRes8[13].re;
        pDst[9].re  = tmpRes8[8].re - tmpRes8[9].re;
        pDst[9].im  = tmpRes8[8].im - tmpRes8[9].im;
        pDst[13].re  = tmpRes8[12].re - tmpRes8[13].im;
        pDst[13].im  = tmpRes8[12].im + tmpRes8[13].re;


        pDst[3].re = tmpRes8[10].re + tmpRes8[11].re;
        pDst[3].im = tmpRes8[10].im + tmpRes8[11].im;
        pDst[7].re = tmpRes8[14].re + tmpRes8[15].im;
        pDst[7].im = tmpRes8[14].im - tmpRes8[15].re;
        pDst[11].re = tmpRes8[10].re - tmpRes8[11].re;
        pDst[11].im = tmpRes8[10].im - tmpRes8[11].im;
        pDst[15].re = tmpRes8[14].re - tmpRes8[15].im;
        pDst[15].im = tmpRes8[14].im + tmpRes8[15].re;




    return;
}
