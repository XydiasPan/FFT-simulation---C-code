/*
Last update : Thursday 18 August 2022

Here is a demonstration of FFT algorithm applied on parallel real data stored 
in a single buffer. This file contains the initialization function, complex type 
of elements and the main FFT algorithm. Also, a time variable is provided to verify
time passed when FFT is applied. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define SAMPLES 2048  //define number of samples that equals to FFT's size
#define PI 3.1415926 //define number of PI 

/*
   cplx -> A struct that holds two float values, real corresponds to complex numbers' real part, 
   while img corresponds to their imaginary part with respect
*/
typedef struct complex{
    float real;
    float img;
}cplx;

/*
   init_fft -> This method initializes array W of type cplx, which holds cosine/sine values on each index.
   Each element is computed to the index returned by the product k*n. In case that product is bigger than N, the cosine/sine value is 
   actually a sum of L*2*PI (where L is integer) + phase.The "phase" part should only be evaluated for that index, but it will not 
   result a new number in comparisson with the current ones stored in the array. 
   That being said, indexes 1*6 and 2*3 should return the same results for the cosines/sines evaluated. 
   By using the methodology described above, we manage to reduce our total RAM space needed to store all the elements in array W.
   Recall that 1 float = 4 Bytes, so 1 cplx item = 8 Bytes.

   fft -> This method applies fft algorithm in a specified buffer and returns a new buffer of type cplx. 
   Each element of the returning buffer will now have a real and img value.

*/

cplx *fft(float inputBuf[SAMPLES]);
void init_fft(int N);
cplx W[SAMPLES];


int main(void){

    /*
       scale_factor -> A typical value used to custom-normalise data visualized in printing commands below

       t -> A clock_t element used to identify the time consumed by our PC to apply fft and return the results. Notice that PCs usually run at 3+ GHz clock, 
       so our mcu board will probably scale the total time up by a factor equal to 3+GHz/fclock MHz.

       inputBuf -> A typical buffer of SAMPLES elements which stores the input data. 
       This buffer was pre-constructed to check if fft works well (e.g, there is a sum of three different sines). In normal cases, this buffer will be updaited on each 
       while loop of the mcu board, and once there are SAMPLE number of elements, fft will be applied.

    */
    int scale_factor = 30;
    clock_t t;
    float phase = 0.0f;
    float inputBuf[SAMPLES];
    for(int k=0;k<SAMPLES;k++){
        phase += 2*PI*250/44100;
        if(phase>= 2*PI) phase -= 2*PI;
        inputBuf[k] = (4*sin(phase) + 2*sin(2*phase) + 1*sin(3*phase));
    }
    
    init_fft(SAMPLES);
    t = clock();
    cplx *outputBuf = fft(inputBuf);
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("time elapsed = %.20f\n",time_taken);
    for(int k=0;k<SAMPLES/8;k++){
        printf("bin %d f = %f",k,(float)k*44100/SAMPLES);
        int mag = ((int)sqrt(outputBuf[k].real*outputBuf[k].real + outputBuf[k].img*outputBuf[k].img))/scale_factor;
        for(int r=0;r<mag;r++) printf("|");
        printf("\n");
    }
    return 0;
}


void init_fft(int N){
    for(int k=0;k<N;k++){
        for(int m=0;m<N/2;m++){
            int idx;
            if(m*k >= N) idx = m*k - N*((int)m*k/N);
            else idx = m*k;
            W[idx].real = cos(2*PI*idx/N);
            W[idx].img = sin(2*PI*idx/N);
        }
    }
}


cplx *fft(float inputBuf[SAMPLES]){
    static cplx out[SAMPLES];
    for(int k=0;k<SAMPLES/2-1;k++){
        cplx xeven;
        cplx xodd;
        xeven.real = 0.0f;
        xeven.img = 0.0f;
        xodd.real = 0.0f;
        xodd.img = 0.0f;
        for(int n=0;n<SAMPLES/2;n++){
            int idx;
            if(k*2*n >= SAMPLES) idx = 2*n*k - SAMPLES*((int)2*n*k/SAMPLES);
            else idx = 2*n*k;
            xeven.real += inputBuf[2*n]*W[idx].real;
            xeven.img += inputBuf[2*n]*W[idx].img;
            if(k*(2*n+1) >= SAMPLES) idx = (2*n+1)*k - SAMPLES*((int)(2*n+1)*k/SAMPLES);
            else idx = (2*n+1)*k;
            xodd.real += inputBuf[2*n+1]*W[idx].real;
            xodd.img += inputBuf[2*n+1]*W[idx].img;
        }
        out[k].real = xeven.real + W[k].real*xodd.real - W[k].img*xodd.img;
        out[k].img = xeven.img + W[k].real*xodd.img + W[k].img*xodd.real;
    }
    return out;
}
