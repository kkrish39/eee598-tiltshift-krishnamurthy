#include <iostream>
#include <pthread.h>
#include <cmath>


/*Function to construct the gaussian kernel with the given sigma value*/
double* constructGaussianKernel(double sigma){
    if(sigma < 0.6)
        return NULL;

    double radius = ceil(2*sigma);

    int kernelSize = (int)(ceil(radius)*2) + 1;
    jdouble * kernelVector = new double[kernelSize] ;

    int wholeRadius = (int)ceil(radius);
    double sigmaSquare = sigma * sigma;
    double twoPiSigmaSquare = 2 * M_PI * sigmaSquare;
    double sqrtTwoPiSigmaSquare = sqrt(twoPiSigmaSquare);

    double firstTerm = 1/(sqrtTwoPiSigmaSquare);

    for(int k=-wholeRadius ;k<=wholeRadius;k++){
        double secondTerm = -1* k*k/(2*sigmaSquare);
        double weight = exp(secondTerm)*firstTerm;

        kernelVector[k+wholeRadius] = weight;
    }

    return kernelVector;
}

/*Sigma Calculation function
 * Based of whether the sigma region is far or near, the sigma values is calculated accordingly
 * */
double calculateSigma(int low, int high, int y, float sigma, bool isFarSigma){
    return  (double) sigma * (isFarSigma ? (high - y) : (y - low))/(high - low);
}

/*Function to perform vertical convolution with constant sigma for every row*/
static void performVerticalConvolutionWithGivenSigma(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, double gaussianKernelVector[]){
    for (int i=top; i<top+width; i++){
        for(int j = i; j < bottom; j=j+width){
            float bluePixel = 0;
            float greenPixel = 0;
            float redPixel = 0;

            int count = -1;
            int pixelVal;

            int rangeToBeConvoluted = radius * width;
            for(int k = j - rangeToBeConvoluted; k <= j + rangeToBeConvoluted; k = k + width) {
                count++;
                if(k < totalPixels) {
                    pixelVal = pixelsIn[k];
                }else{
                    pixelVal = 0;
                }

                int blue = pixelVal & 0xff;
                int green = (pixelVal >> 8) & 0xff;
                int red = (pixelVal >> 16) & 0xff;

                redPixel += (gaussianKernelVector[count] * red );
                greenPixel += (gaussianKernelVector[count] * green );
                bluePixel += (gaussianKernelVector[count] * blue );
            }

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;

            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
        }
    }
}

/*Function to perform horizontal convolution with constant sigma for every row*/
static void performHorizontalConvolutionWithGivenSigma(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, double gaussianKernelVector[]){
    for(int i=top;i<bottom;i=i+width){
        int pixelRight = i + width;
        for(int j = i; j < pixelRight; j++) {
            float bluePixel = 0;
            float greenPixel = 0;
            float redPixel = 0;

            int pixelVal;
            for(int k = -radius; k <= radius; k++) {
                int pixelIndex = j + k;
                int vectorIndex = k + radius;
                if(pixelIndex >= 0 && pixelIndex < pixelRight && pixelIndex < totalPixels){
                    pixelVal = pixelsIn[pixelIndex];
                }else{
                    pixelVal = 0;
                }

                int blue = pixelVal & 0xff;
                int green = (pixelVal >> 8) & 0xff;
                int red = (pixelVal >> 16) & 0xff;

                redPixel += (gaussianKernelVector[vectorIndex] * red);
                greenPixel += (gaussianKernelVector[vectorIndex] * green);
                bluePixel += (gaussianKernelVector[vectorIndex] * blue);
            }

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;


            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
        }
    }
}
/*Function to perform vertical convolution with varying sigma for every row*/
static void performVerticalConvolution(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
    int lowIndex = top/width;
    int highIndex = bottom/width;
    for (int i=top; i<=top+width; i++){
        for(int j = i; j < bottom; j=j+width){
            float bluePixel = 0;
            float greenPixel = 0;
            float redPixel = 0;

            int count = -1;
            int pixelVal;

            double sigmaVal = calculateSigma(lowIndex,highIndex,j/width,sigma,isSigmaFar);
            if(sigmaVal < 0.6){
                pixelsOut[j] = pixelsIn[j];
                continue;
            }
            double radiusVal = ceil(2*sigmaVal);
            int kernelSize = (int)(ceil(radiusVal)*2) + 1;
            double *gaussianKernelVector;
            gaussianKernelVector = constructGaussianKernel(sigmaVal);


            if(gaussianKernelVector == NULL){
                pixelsOut[j] = pixelsIn[j];
                continue;
            }
            int radius = kernelSize/2;
            int rangeToBeConvoluted = radius * width;

            for(int k = j - rangeToBeConvoluted; k <= j + rangeToBeConvoluted; k = k + width) {
                count++;
                if(k < totalPixels){
                    pixelVal = pixelsIn[k];
                }else{
                    pixelVal = 0;
                }

                int blue = pixelVal & 0xff;
                int green = (pixelVal >> 8) & 0xff;
                int red = (pixelVal >> 16) & 0xff;


                redPixel += (gaussianKernelVector[count] * red );
                greenPixel += (gaussianKernelVector[count] * green );
                bluePixel += (gaussianKernelVector[count] * blue );

            }
            /*de allocating the Gaussian kernel vector*/
            delete [] gaussianKernelVector;

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;

            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
        }
    }
}
/*Function to perform horizontal convolution with varying sigma for every row*/
static void performHorizontalConvolution(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
    int lowIndex = top/width;
    int highIndex = bottom/width;
    for(int i=top;i<bottom;i=i+width){
        int pixelRight = i + width - 1;
        double sigmaVal = calculateSigma(lowIndex,highIndex,i/width,sigma,isSigmaFar);

        if(sigmaVal < 0.6){
            for(int x=i;x<=pixelRight-i+1;x++){
                pixelsOut[x] = pixelsIn[x];
            }
            continue;
        }

        double radiusVal = ceil(2*sigmaVal);
        int kernelSize = (int)(ceil(radiusVal)*2) + 1;
        double *gaussianKernelVector;
        gaussianKernelVector = constructGaussianKernel(sigmaVal);

        int radius = kernelSize/2;

        for(int j = i; j <= pixelRight; j++) {
            float bluePixel = 0;
            float greenPixel = 0;
            float redPixel = 0;

            int pixelVal;

            for(int k = -radius; k <= radius; k++) {
                int pixelIndex = j + k;
                int vectorIndex = k + radius;

                if(pixelIndex >= 0 && pixelIndex < pixelRight && pixelIndex < totalPixels){
                    pixelVal = pixelsIn[pixelIndex];
                }else{
                    pixelVal = 0;
                }

                int blue = pixelVal & 0xff;
                int green = (pixelVal >> 8) & 0xff;
                int red = (pixelVal >> 16) & 0xff;

                redPixel += (gaussianKernelVector[vectorIndex] * red);
                greenPixel += (gaussianKernelVector[vectorIndex] * green);
                bluePixel += (gaussianKernelVector[vectorIndex] * blue);
            }
            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;
            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);

        }
        /*de allocating the Gaussian kernel vector*/
        delete [] gaussianKernelVector;
    }
}

/*Function to be called by the thread with the @struct threadArgs based on region of interest*/
void *performConvolution(void *threadarg){
    struct threadArgs *args;
    args = (struct threadArgs *) threadarg;

    int top = args->top;
    int bottom = args->bottom;
    int *pixelsIn = args->pixelsIn;
    int *pixelsIntermediate = args ->pixelsIntermediate;
    int *pixelsOut  = args->pixelsOut;
    int width = args ->width;
    int totalPixels = args -> totalPixels;
    float sigma = args ->sigma;
    bool isSigmaFar = args->isSigmaFar;
    bool singleSigma = args->singleSigma;

    if(sigma <= 0){
        for(int i=top;i<bottom;i++){
            pixelsOut[i] = pixelsIn[i];
        }
        pthread_exit(NULL);
    }

    if(singleSigma){
        double radius = ceil(2*sigma);
        int kernelSize = (int)(ceil(radius)*2) + 1;
        double *gaussianKernelVector = constructGaussianKernel(sigma);

        performVerticalConvolutionWithGivenSigma(top, bottom, pixelsIn, pixelsIntermediate, width, totalPixels, kernelSize/2, gaussianKernelVector);
        performHorizontalConvolutionWithGivenSigma(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, kernelSize/2,gaussianKernelVector);
        delete [] gaussianKernelVector;
    }else{
        performVerticalConvolution(top, bottom, pixelsIn, pixelsIntermediate,width, sigma, isSigmaFar, totalPixels);
        performHorizontalConvolution(top, bottom, pixelsIntermediate, pixelsOut, width, sigma, isSigmaFar, totalPixels);
    }

    pthread_exit(NULL);
}