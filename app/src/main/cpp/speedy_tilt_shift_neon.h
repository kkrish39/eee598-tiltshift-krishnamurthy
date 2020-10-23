#include <iostream>
#include <pthread.h>
#include <cmath>
#include <arm_neon.h>

double* constructGaussianKernelNeon(double sigma){
    if(sigma < 0.6)
        return NULL;

    double radius = ceil(2*sigma);

    int kernelSize = (int)(ceil(radius)*2) + 1;
    double * kernelVector = new double[kernelSize] ;

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


double calculateSigmaNeon(int low, int high, int y, float sigma, bool isFarSigma){
    return  (double) sigma * (isFarSigma ? (high - y) : (y - low))/(high - low);
}

static void performVerticalConvolutionWithGivenSigmaNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, double gaussianKernelVector[]){
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

static void performHorizontalConvolutionWithGivenSigmaNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, double gaussianKernelVector[]){
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

static void performVerticalConvolutionNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
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

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;

            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
        }
    }
}

static void performHorizontalConvolutionNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
    int lowIndex = top/width;
    int highIndex = bottom/width;
    for(int i=top;i<bottom;i=i+width){
        int pixelRight = i + width - 1;
        double sigmaVal = calculateSigmaNeon(lowIndex,highIndex,i/width,sigma,isSigmaFar);

        if(sigmaVal < 0.6){
            for(int x=i;x<=pixelRight-i+1;x++){
                pixelsOut[x] = pixelsIn[x];
            }
            continue;
        }

        double radiusVal = ceil(2*sigmaVal);
        int kernelSize = (int)(ceil(radiusVal)*2) + 1;
        double *gaussianKernelVector;
        gaussianKernelVector = constructGaussianKernelNeon(sigmaVal);

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
    }
}

void performConvolutionNeon(int top, int bottom, int *pixelsIn, int *pixelsIntermediate, int *pixelsOut, int height, int width, int totalPixels, float sigma, bool isSigmaFar, bool singleSigma){

    if(sigma <= 0){
        for(int i=top;i<bottom;i++){
            pixelsOut[i] = pixelsIn[i];
        }
        return;
    }

    if(singleSigma){
        double radius = ceil(2*sigma);
        int kernelSize = (int)(ceil(radius)*2) + 1;
        double *gaussianKernelVector = constructGaussianKernelNeon(sigma);

        performVerticalConvolutionWithGivenSigmaNeon(top, bottom, pixelsIn, pixelsIntermediate, width, totalPixels, kernelSize/2, gaussianKernelVector);
        performHorizontalConvolutionWithGivenSigmaNeon(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, kernelSize/2,gaussianKernelVector);
    }else{
        performVerticalConvolutionNeon(top, bottom, pixelsIn, pixelsIntermediate,width, sigma, isSigmaFar, totalPixels);
        performHorizontalConvolutionNeon(top, bottom, pixelsIntermediate, pixelsOut, width, sigma, isSigmaFar, totalPixels);
    }

}