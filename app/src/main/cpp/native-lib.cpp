#include <jni.h>
#include <string>
#include <cpu-features.h>
#include <android/log.h>
#include <arm_neon.h>
#include <chrono>
#include "neon_cpp_commom.h"
#include "speedy_tilt_shift_cpp.h"
#include "speedy_tilt_shift_neon.h"
#define NUM_THREADS 4

using namespace std::chrono;

extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftcppnative(JNIEnv *env, jobject instance, jintArray inputPixels_, jintArray outputPixels_, jint width,
        jint height, jfloat sigma_far, jfloat sigma_near, jint a0, jint a1, jint a2, jint a3) {

    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);
    jint *pixelsIntermediate = env->GetIntArrayElements(outputPixels_, NULL);

    int totalPixels = width*height;
    __android_log_print(ANDROID_LOG_DEBUG, "a0", "a0:%d  a1:%d  a2:%d   a3:%d   sigma_far:%f  sigma_near:%f", a0, a1, a2, a3, sigma_far, sigma_near);

    /*Initializing four threads to operate in four region. No blur region is handled in the main thread itself*/
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

    auto start = high_resolution_clock::now();
    for(int i=0;i<NUM_THREADS;i++){
        pthread_create(&threads[i], NULL, performConvolution, (void *)&th[i]);
    }

    for(int i=(a1*width);i<=(a2*width);i++){
        outputPixels[i] = pixels[i];
    }

    for(int i=0;i<NUM_THREADS;i++){
        pthread_join(threads[i], NULL);
    }
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    __android_log_print(ANDROID_LOG_INFO, "Execution time - CPP (microseconds)", "%lld", duration);


    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}

extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftneonnative(JNIEnv *env, jclass instance, jintArray inputPixels_, jintArray outputPixels_, jint width,
        jint height, jfloat sigma_far, jfloat sigma_near, jint a0, jint a1,jint a2, jint a3) {

    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);
    jint *pixelsIntermediate = env->GetIntArrayElements(outputPixels_, NULL);

    __android_log_print(ANDROID_LOG_DEBUG, "a0", "a0:%d  a1:%d  a2:%d   a3:%d   sigma_far:%f  sigma_near:%f", a0, a1, a2, a3, sigma_far, sigma_near);

    uint8_t * arrayInPtr = (uint8_t *)pixels;
    uint8_t * arrayOutPtr = (uint8_t *)outputPixels;

    int totalPixels  = width*height;
    uint8x16x4_t pixelChannels;

/*Initializing four threads to operate in four region. No blur region is handled in the main thread itself*/
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

    auto start = high_resolution_clock::now();
    for(int i=0;i<NUM_THREADS;i++){
        pthread_create(&threads[i], NULL, performConvolutionNeon, (void *)&th[i]);
    }

    for(int i=0;i<NUM_THREADS;i++){
        pthread_join(threads[i], NULL);
    }

    int rangeVal = ((a2-a1)*width)/16;
    arrayInPtr = arrayInPtr + (a1*width*4);
    arrayOutPtr = arrayOutPtr + (a1*width*4);

    /*Straight forward copy from given image to the output using NEON*/
    for(int i=0;i<rangeVal;i++){
        pixelChannels = vld4q_u8(arrayInPtr);
        vst4q_u8(arrayOutPtr,pixelChannels);
        arrayOutPtr = arrayOutPtr+64;
        arrayInPtr = arrayInPtr+64;
    }

    auto stop = high_resolution_clock::now();

    /*Clock duration in micro seconds*/
    auto duration = duration_cast<microseconds>(stop - start);
    __android_log_print(ANDROID_LOG_INFO, "Execution time - NEON (microseconds)", "%lld", duration);

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}