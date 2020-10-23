#include <jni.h>
#include <string>
#include <math.h>
#include <cpu-features.h>
#include <android/log.h>
#include <arm_neon.h>
#include <iostream>
#include <pthread.h>

#define NUM_THREADS 4

struct threadArgs{
    int top;
    int bottom;
    int *pixelsIn;
    int *pixelsIntermediate;
    int *pixelsOut;
    int height;
    int width;
    int totalPixels;
    float sigma;
    bool isSigmaFar;
    bool singleSigma;
};

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

            int combinedAlpha = 0xff;
            int combinedRed = (int) redPixel;
            int combinedGreen = (int) greenPixel;
            int combinedBlue = (int) bluePixel;

            pixelsOut[j] = (combinedAlpha & 0xff) << 24 | (combinedRed & 0xff) << 16 | (combinedGreen & 0xff) << 8 | (combinedBlue & 0xff);
        }
    }
}

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
//        __android_log_print(ANDROID_LOG_ERROR, "RADIUS ", "%d", radius);
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

        }
    }
}

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

        if(kernelSize == 0){ pthread_exit(NULL);}

        performVerticalConvolutionWithGivenSigma(top, bottom, pixelsIn, pixelsIntermediate, width, totalPixels, kernelSize/2, gaussianKernelVector);
        performHorizontalConvolutionWithGivenSigma(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, kernelSize/2,gaussianKernelVector);
    }else{
        performVerticalConvolution(top, bottom, pixelsIn, pixelsIntermediate,width, sigma, isSigmaFar, totalPixels);
        performHorizontalConvolution(top, bottom, pixelsIntermediate, pixelsOut, width, sigma, isSigmaFar, totalPixels);
    }

    pthread_exit(NULL);
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

    __android_log_print(ANDROID_LOG_ERROR, "a0", "%d", a0);
    __android_log_print(ANDROID_LOG_ERROR, "a1", "%d", a1);
    __android_log_print(ANDROID_LOG_ERROR, "a2", "%d", a2);
    __android_log_print(ANDROID_LOG_ERROR, "a3", "%d", a3);
    __android_log_print(ANDROID_LOG_ERROR, "sigma_far", "%f", sigma_far);
    __android_log_print(ANDROID_LOG_ERROR, "sigma_near", "%f", sigma_near);


    pthread_t threads[NUM_THREADS];
    struct threadArgs th[NUM_THREADS];

    th[0].top = 0;
    th[0].bottom = a0*width;
    th[0].pixelsIn = pixels;
    th[0].pixelsIntermediate = pixelsIntermediate;
    th[0].pixelsOut = outputPixels;
    th[0].height = height;
    th[0].width = width;
    th[0].totalPixels = totalPixels;
    th[0].sigma = sigma_far;
    th[0].isSigmaFar = true;
    th[0].singleSigma = true;

    th[1].top = a0*width;
    th[1].bottom = a1*width;
    th[1].pixelsIn = pixels;
    th[1].pixelsIntermediate = pixelsIntermediate;
    th[1].pixelsOut = outputPixels;
    th[1].height = height;
    th[1].width = width;
    th[1].totalPixels = totalPixels;
    th[1].sigma = sigma_far;
    th[1].isSigmaFar = true;
    th[1].singleSigma = false;

    th[2].top = a2*width;
    th[2].bottom = a3*width;
    th[2].pixelsIn = pixels;
    th[2].pixelsIntermediate = pixelsIntermediate;
    th[2].pixelsOut = outputPixels;
    th[2].height = height;
    th[2].width = width;
    th[2].totalPixels = totalPixels;
    th[2].sigma = sigma_near;
    th[2].isSigmaFar = false;
    th[2].singleSigma = false;

    th[3].top = a3*width;
    th[3].bottom = height*width;
    th[3].pixelsIn = pixels;
    th[3].pixelsIntermediate = pixelsIntermediate;
    th[3].pixelsOut = outputPixels;
    th[3].height = height;
    th[3].width = width;
    th[3].totalPixels = totalPixels;
    th[3].sigma = sigma_near;
    th[3].isSigmaFar = false;
    th[3].singleSigma = true;

    for(int i=0;i<NUM_THREADS;i++){
        pthread_create(&threads[i], NULL, performConvolution, (void *)&th[i]);
    }

    for(int i=(a1*width);i<=(a2*width);i++){
        outputPixels[i] = pixels[i];
    }

    for(int i=0;i<NUM_THREADS;i++){
        pthread_join(threads[i], NULL);
    }
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


    __android_log_print(ANDROID_LOG_ERROR, "a0", "%d", a0);
    __android_log_print(ANDROID_LOG_ERROR, "a1", "%d", a1);
    __android_log_print(ANDROID_LOG_ERROR, "a2", "%d", a2);
    __android_log_print(ANDROID_LOG_ERROR, "a3", "%d", a3);
    __android_log_print(ANDROID_LOG_ERROR, "sigma_far", "%f", sigma_far);
    __android_log_print(ANDROID_LOG_ERROR, "sigma_near", "%f", sigma_near);

    uint8_t * arrayInPtr = (uint8_t *)pixels;
    uint8_t * arrayOutPtr = (uint8_t *)outputPixels;

    int totalPixels  = width*height;
    int numIterations = totalPixels / 16;
    uint8x16x4_t pixelChannels;


    sigma_far = 3.0f;
    double radius = ceil(2*sigma_far);
    int kernelSize = (int)(ceil(radius)*2) + 1;
    double *gaussianKernelVector = constructGaussianKernel(sigma_far);

    if(kernelSize > 0) {
        for(int i=0;i<kernelSize;i++){
            __android_log_print(ANDROID_LOG_ERROR, "GK", "%f", gaussianKernelVector[i]);
        }

        for (int i = 0; i < a0 * width; i++) {

        }
    }


    for(int i=0;i<numIterations;i++){
        pixelChannels = vld4q_u8(arrayInPtr);
        uint8x16_t b_vector = pixelChannels.val[0];
        uint8x16_t g_vector = pixelChannels.val[1];
        uint8x16_t r_vector = pixelChannels.val[2];
        uint8x16_t a_vector = pixelChannels.val[3];

        pixelChannels.val[0] = b_vector;
        pixelChannels.val[1] = g_vector;
        pixelChannels.val[2] = r_vector;
        pixelChannels.val[3] = a_vector;


        vst4q_u8(arrayOutPtr,pixelChannels);
        arrayOutPtr = arrayOutPtr+64;
        arrayInPtr = arrayInPtr+64;
    }

//    outputPixels = (jint *)(arrayOutPtr);

//
//    for(int i=0;i<num8x16; i=i++){
//        pixelChannels = vld4q_u8(arrayInPtr + 4*16*i);
//        uint8x16_t b_vector = pixelChannels.val[0];
//        uint8x16_t g_vector = pixelChannels.val[1];
//        uint8x16_t r_vector = pixelChannels.val[2];
//        uint8x16_t a_vector = pixelChannels.val[3];
//
//        uint32x4_t combinedPixel;
//
//        uint32_t combined_b_vector = vgetq_lane_u8(b_vector, 0);
//        uint32_t combined_g_vector = vgetq_lane_u8(g_vector, 0);
//        uint32_t combined_r_vector = vgetq_lane_u8(r_vector, 0);
//        uint32_t combined_a_vector = vgetq_lane_u8(a_vector, 0);
//
//        vsetq_lane_f32(combined_b_vector,combinedPixel,3);
//        vsetq_lane_f32(combined_g_vector,combinedPixel,2);
//        vsetq_lane_f32(combined_r_vector,combinedPixel,1);
//        vsetq_lane_f32(combined_a_vector,combinedPixel,0);
//
//        float sigma = 0.9f;
//
//        combinedPixel = vmulq_n_u32(combinedPixel,(uint32_t) sigma*64);
//        combinedPixel = vshrq_n_s32(combinedPixel,6);
//
//
//        __android_log_print(ANDROID_LOG_ERROR, "Marker", "%s", "----------------");
//        __android_log_print(ANDROID_LOG_ERROR, "Marker", "%s", "----------------");
//
//
//
//        vst4q_u32(arrayOutPtr,combinedPixel)
//        arrayOutPtr = arrayOutPtr + 4*16;
//    }

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}