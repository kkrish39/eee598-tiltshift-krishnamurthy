#include <jni.h>
#include <string>
#include <cpu-features.h>
#include <android/log.h>
#include <arm_neon.h>
#include "speedy_tilt_shift_cpp.h"
#include "speedy_tilt_shift_neon.h"
#define NUM_THREADS 4

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

    __android_log_print(ANDROID_LOG_DEBUG, "a0", "%d", a0);
    __android_log_print(ANDROID_LOG_DEBUG, "a1", "%d", a1);
    __android_log_print(ANDROID_LOG_DEBUG, "a2", "%d", a2);
    __android_log_print(ANDROID_LOG_DEBUG, "a3", "%d", a3);
    __android_log_print(ANDROID_LOG_DEBUG, "sigma_far", "%f", sigma_far);
    __android_log_print(ANDROID_LOG_DEBUG, "sigma_near", "%f", sigma_near);


    pthread_t threads[NUM_THREADS];
    struct threadArgs th[NUM_THREADS];

    th[0].top = 0;
    th[0].bottom = a0*width;
    th[0].pixelsIn = pixels;
    th[0].pixelsIntermediate = pixelsIntermediate;
    th[0].pixelsOut = outputPixels;
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