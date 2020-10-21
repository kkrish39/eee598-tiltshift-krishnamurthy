#include <jni.h>
#include <string>
#include <math.h>
#include <cpu-features.h>
#include <android/log.h>


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


double calculateSigma(int low, int high, int y, float sigma, bool isFarSigma){
    return  (double) sigma * (isFarSigma ? (high - y) : (y - low))/(high - low);
}

private static void performVerticalConvolutionWithGivenSigma(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, double gaussianKernelVector[]){
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

private static void performHorizontalConvolutionWithGivenSigma(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, int totalPixels, int radius, double gaussianKernelVector[]){
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

private static void performVerticalConvolution(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
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

private static void performHorizontalConvolution(int top, int bottom, int *pixelsIn, int *pixelsOut, int width, double sigma, bool isSigmaFar, int totalPixels){
    int lowIndex = top/width;
    int highIndex = bottom/width;
    for(int i=top;i<bottom;i=i+width){
        int pixelRight = i + width - 1;
        double sigmaVal = calculateSigma(lowIndex,highIndex,i/width,sigma,isSigmaFar);

        if(sigmaVal < 0.6){
            for(int x=i;x<pixelRight-i+1;x++){
                pixelsOut[x] = pixelsIn[x];
            }
            continue;
        }

        double radiusVal = ceil(2*sigmaVal);
        int kernelSize = (int)(ceil(radiusVal)*2) + 1;
        double *gaussianKernelVector;
        gaussianKernelVector = constructGaussianKernel(sigmaVal);

        if(gaussianKernelVector == NULL){
            for(int x=i;x<pixelRight-i+1;x++){
                pixelsOut[x] = pixelsIn[x];
            }
            continue;
        }

        int radius = kernelSize/2;
        __android_log_print(ANDROID_LOG_ERROR, "RADIUS ", "%d", radius);
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


//                __android_log_print(ANDROID_LOG_ERROR, "PixelIndex", "%d", pixelIndex);
//                __android_log_print(ANDROID_LOG_ERROR, "pixelRight", "%d", pixelRight);
//                __android_log_print(ANDROID_LOG_ERROR, "Total", "%d", totalPixels);
//                __android_log_print(ANDROID_LOG_ERROR, "P in", "%d", pixelIndex > 0 ? pixelsIn[pixelIndex] : 0);

            }

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;
            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
//            __android_log_print(ANDROID_LOG_ERROR, "OUT", "%d", pixelsOut[j]);
        }
    }
}

void performConvolution(int top, int bottom, int *pixelsIn, int *pixelsIntermediate, int *pixelsOut, int height, int width, int totalPixels, float sigma, bool isSigmaFar, bool singleSigma){
    if(singleSigma){
        double radius = ceil(2*sigma);
        int kernelSize = (int)(ceil(radius)*2) + 1;
        double *gaussianKernelVector = constructGaussianKernel(sigma);

        if(kernelSize == 0){return;}

        performVerticalConvolutionWithGivenSigma(top, bottom, pixelsIn, pixelsIntermediate, width, totalPixels, kernelSize/2, gaussianKernelVector);
        performHorizontalConvolutionWithGivenSigma(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, kernelSize/2,gaussianKernelVector);
    }else {

        performVerticalConvolution(top, bottom, pixelsIn, pixelsIntermediate,width, sigma, isSigmaFar, totalPixels);
        performHorizontalConvolution(top, bottom, pixelsIntermediate, pixelsOut, width, sigma, isSigmaFar, totalPixels);
    }
}


extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftcppnative(JNIEnv *env,
                                                                               jobject instance,
                                                                               jintArray inputPixels_,
                                                                               jintArray outputPixels_,
                                                                               jint width,
                                                                               jint height,
                                                                               jfloat sigma_far,
                                                                               jfloat sigma_near,
                                                                               jint a0, jint a1,
                                                                               jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);
    jint *pixelsIntermediate = env->GetIntArrayElements(outputPixels_, NULL);
    int totalPixels = width*height;
    __android_log_print(ANDROID_LOG_ERROR, "Total Pixels ", "%d", width);
    __android_log_print(ANDROID_LOG_ERROR, "Total Pixels ", "%d", height);
    __android_log_print(ANDROID_LOG_ERROR, "Total Pixels ", "%d", totalPixels);
    sigma_far = 5.0;
    a0 = 0;
    a1 = height;

    performConvolution(a0*width,a1*width,pixels,pixelsIntermediate,outputPixels,height,width,totalPixels,sigma_far,true, false);

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
     return 0;
}

extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftneonnative(JNIEnv *env,
                                                                               jclass instance,
                                                                               jintArray inputPixels_,
                                                                               jintArray outputPixels_,
                                                                               jint width,
                                                                               jint height,
                                                                               jfloat sigma_far,
                                                                               jfloat sigma_near,
                                                                               jint a0, jint a1,
                                                                               jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);

    for (int j=0;j<height;j++){
        for (int i=0;i<width;i++) {
            int B = pixels[j*width+i]%0x100;
            int G = (pixels[j*width+i]>>8)%0x100;
            int R = (pixels[j*width+i]>>16)%0x100;
            int A = 0xff;
            R=0;
            int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 8 | (B & 0xff);

            outputPixels[j*width+i]=color;
        }
    }

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}