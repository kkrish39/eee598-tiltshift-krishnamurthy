#include <iostream>
#include <pthread.h>
#include <cmath>
#include <arm_neon.h>

/*Construct gaussian kernel with the given Sigma Value*/
float* constructGaussianKernelNeon(double sigma){
    if(sigma < 0.6)
        return NULL;
    double radius = ceil(2*sigma);

    int kernelSize = (int)(ceil(radius)*2) + 1;
    float * kernelVector = new float[kernelSize] ;

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

/*Calculate the sigma value based on the current row and the indexes*/
float calculateSigmaNeon(int low, int high, int y, float sigma, bool isFarSigma){
    if(isFarSigma)
        return (sigma * (high-y) + 1.0f*(y-low)) / (high-low);

    return (1.0f * (high-y) + sigma * (y-low)) / (high-low);
}


/*
 * px - pixels to be processed
 * blur_val - the gaussian blur value that needs to be applied
 * */
uint32x4_t processPixel(int *px, float blur_val){
    uint8_t *eight_pxl_arr = (uint8_t *)px;

    /*Reading the first 2 pixels without interleaving*/
    uint8x8_t eight_pxl_vector = vld1_u8(eight_pxl_arr);

    /*Scaling up the value from 8bit to 16bit vectors*/
    uint16x8_t scaled_up_vector = vmovl_u8(eight_pxl_vector);

    /*Reading only the first pixel value*/
    uint16x4_t low_scaled_up_vector = vget_low_u16(scaled_up_vector);

    /*Applying gaussian_blur with the given sigma value after scaling it up to a factor of 64*/
    low_scaled_up_vector = vmul_n_u16(low_scaled_up_vector, (uint16_t)(blur_val*64));

    /*Removing the scale value of 64 shifting right to 6 digits*/
    low_scaled_up_vector = vshr_n_s16(low_scaled_up_vector,6);

    /*scaling the 16bit value to 32 bit and return*/
    return vmovl_u32(low_scaled_up_vector);
}

/*Function to perform vertical convolution with constant sigma for every row*/
static void performVerticalConvolutionWithGivenSigmaNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, float gaussianKernelVector[]){
    for (int i=top; i<top+width; i++){
        for(int j = i; j < bottom; j=j+width){
            float bluePixel = 0;
            float greenPixel = 0;
            float redPixel = 0;

            int count = -1;

            int rangeToBeConvoluted = radius * width;
            for(int k = j - rangeToBeConvoluted; k <= j + rangeToBeConvoluted; k = k + width) {
                count++;
                /*Store the processed Pixel*/
                uint32x4_t processedPixel;
                if(k < totalPixels) {
                    processedPixel = processPixel(&pixelsIn[k],gaussianKernelVector[count]);
                }else{
                    processedPixel = vdupq_n_u32(0);
                }

                /*Getting each channel value of 32 bits from the specific lane*/
                int redLane = vgetq_lane_u32(processedPixel, 2);
                int greenLane = vgetq_lane_u32(processedPixel, 1);
                int blueLane = vgetq_lane_u32(processedPixel, 0);

                /*Aggregating the pixel values over the given gaussian kernel size*/
                redPixel += (redLane & 0xff) << 16;
                greenPixel += (greenLane & 0xff) << 8;
                bluePixel += (blueLane & 0xff);
            }

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;

            /*Re constructing the pixel with the new channel values*/
            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
        }
    }
}

/*Function to perform horizontal convolution with constant sigma for every row*/
static void performHorizontalConvolutionWithGivenSigmaNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, float gaussianKernelVector[]){
    for(int i=top;i<bottom;i=i+width){
        int pixelRight = i + width;
        for(int j = i; j < pixelRight; j++) {
            float bluePixel = 0;
            float greenPixel = 0;
            float redPixel = 0;

            for(int k = -radius; k <= radius; k++) {
                int pixelIndex = j + k;
                int vectorIndex = k + radius;
                uint32x4_t processedPixel;
                if(pixelIndex >= 0 && pixelIndex < pixelRight && pixelIndex < totalPixels){
                    processedPixel = processPixel(&pixelsIn[pixelIndex],gaussianKernelVector[vectorIndex]);
                }else{
                    processedPixel = vdupq_n_u32(0);
                }

                int redLane = vgetq_lane_u32(processedPixel, 2);
                int greenLane = vgetq_lane_u32(processedPixel, 1);
                int blueLane = vgetq_lane_u32(processedPixel, 0);

                redPixel += (redLane & 0xff) << 16;
                greenPixel += (greenLane & 0xff) << 8;
                bluePixel += (blueLane & 0xff);
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
static void performVerticalConvolutionNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
    int lowIndex = top/width;
    int highIndex = bottom/width;
    for (int i=top; i<=top+width; i++){
        for(int j = i; j < bottom; j=j+width){
            float bluePixel = 0;
            float greenPixel = 0;
            float redPixel = 0;

            int count = -1;

            double sigmaVal = calculateSigma(lowIndex,highIndex,j/width,sigma,isSigmaFar);
            if(sigmaVal < 0.6){
                pixelsOut[j] = pixelsIn[j];
                continue;
            }
            double radiusVal = ceil(2*sigmaVal);
            int kernelSize = (int)(ceil(radiusVal)*2) + 1;
            float *gaussianKernelVector;
            gaussianKernelVector = constructGaussianKernelNeon(sigmaVal);


            if(gaussianKernelVector == NULL){
                pixelsOut[j] = pixelsIn[j];
                continue;
            }
            int radius = kernelSize/2;
            int rangeToBeConvoluted = radius * width;

            for(int k = j - rangeToBeConvoluted; k <= j + rangeToBeConvoluted; k = k + width) {
                count++;
                uint32x4_t processedPixel;
                if(k < totalPixels) {
                    processedPixel = processPixel(&pixelsIn[k],gaussianKernelVector[count]);
                }else{
                    processedPixel = vdupq_n_u32(0);
                }

                int redLane = vgetq_lane_u32(processedPixel, 2);
                int greenLane = vgetq_lane_u32(processedPixel, 1);
                int blueLane = vgetq_lane_u32(processedPixel, 0);

                redPixel += (redLane & 0xff);
                greenPixel += (greenLane & 0xff) << 8;
                bluePixel += (blueLane & 0xff)<< 16;
            }
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
static void performHorizontalConvolutionNeon(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
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
        float *gaussianKernelVector;
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
                uint32x4_t processedPixel;
                if(pixelIndex >= 0 && pixelIndex < pixelRight && pixelIndex < totalPixels){
                    processedPixel = processPixel(&pixelsIn[pixelIndex],gaussianKernelVector[vectorIndex]);
                }else{
                    processedPixel = vdupq_n_u32(0);
                }

                int redLane = vgetq_lane_u32(processedPixel, 2);
                int greenLane = vgetq_lane_u32(processedPixel, 1);
                int blueLane = vgetq_lane_u32(processedPixel, 0);

                redPixel += (redLane & 0xff) ;
                greenPixel += (greenLane & 0xff) << 8;
                bluePixel += (blueLane & 0xff)<< 16;
            }

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;
            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
        }
        delete [] gaussianKernelVector;
    }
}

/*Function to be called by the thread with the @struct threadArgs based on region of interest*/
void *performConvolutionNeon(void *threadarg){
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
        float *gaussianKernelVector = constructGaussianKernelNeon(sigma);

        performVerticalConvolutionWithGivenSigmaNeon(top, bottom, pixelsIn, pixelsIntermediate, width, totalPixels, kernelSize/2, gaussianKernelVector);
        performHorizontalConvolutionWithGivenSigmaNeon(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, kernelSize/2,gaussianKernelVector);
    }else{
        performVerticalConvolutionNeon(top, bottom, pixelsIn, pixelsIntermediate,width, sigma, isSigmaFar, totalPixels);
        performHorizontalConvolutionNeon(top, bottom, pixelsIntermediate, pixelsOut, width, sigma, isSigmaFar, totalPixels);
    }

    pthread_exit(NULL);
}