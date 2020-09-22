package edu.asu.ame.meteor.speedytiltshift2018;

import android.graphics.Bitmap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SpeedyTiltShift {
    static SpeedyTiltShift Singleton = new SpeedyTiltShift();

    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }


    private static double[] constructGaussianKernel(double sigma){
        double radius = Math.ceil(2*sigma);

        int kernelSize = (int)(Math.ceil(radius)*2) + 1;
//        System.out.println("Kernel Size ---> "+"  "+radius+"  "+ sigma +"  "+kernelSize);
        double[] kernelVector = new double[kernelSize];

        int wholeRadius = (int)Math.ceil(radius);
        double sigmaSquare = sigma * sigma;
        double twoPiSigmaSquare = 2 * Math.PI * sigmaSquare;
        double sqrtTwoPiSigmaSquare = Math.sqrt(twoPiSigmaSquare);

        double firstTerm = 1/(sqrtTwoPiSigmaSquare);


        for(int k=-wholeRadius ;k<=wholeRadius;k++){
            double secondTerm = -1* k*k/(2*sigmaSquare);
            double weight = Math.exp(secondTerm)*firstTerm;

            kernelVector[k+wholeRadius] = weight;
        }

        return kernelVector;
    }

    private static double calculateSigma(int low, int high, int y, float sigma, boolean isFarSigma){
        double calculatedSigma = sigma * (isFarSigma ? (high - y) : (y - low))/(high - low);

        return calculatedSigma < 0.6 ? sigma : calculatedSigma;
    }

    private static void performVerticalConvolution(int top, int bottom, int[] pixelsIn, int[] pixelsOut, int width, double[][] gaussianKernelVectors){
        for (int i=top; i<top+width; i++){
            for(int j = i; j < bottom; j=j+width) {
                float bluePixel = 0;
                float greenPixel = 0;
                float redPixel = 0;

                int count = -1;
                int pixelVal;
                double[] gaussianKernelVector = gaussianKernelVectors[j/width];
                int radius = gaussianKernelVector.length/2;
                for(int k = j - (radius * width); k <= j + (radius * width); k = k + width) {
                    count++;
                    // TODO: Must check the time taken to process
                    try{
                        pixelVal = pixelsIn[k];
                    }catch (ArrayIndexOutOfBoundsException a){
                        pixelVal = 1;
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
    private static void performHorizontalConvolution(int top, int bottom, int[] pixelsIn, int[] pixelsOut, int width, int totalPixels, double[][] gaussianKernelVectors){
        for(int i=top;i<bottom+width;i=i+width){
            int pixelRight = i + width - 1;

            double[] gaussianKernelVector = gaussianKernelVectors[i/width];
            int radius = gaussianKernelVector.length/2;

            for(int j = i; j< pixelRight; j++) {
                float bluePixel = 0;
                float greenPixel = 0;
                float redPixel = 0;

                int pixelVal;

                for(int k = -radius; k <= radius; k++) {
                    int pixelIndex = j + k;
                    int vectorIndex = k + radius;

                    // TODO: Must check the time taken to process
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
    private static void performVerticalConvolutionWithGivenSigma(int top, int[] pixelsIn, int[] pixelsOut, int width, int totalPixels, int radius, double[] gaussianKernelVector){
        for (int i=top; i<top+width; i++){
            for(int j = i; j < totalPixels; j=j+width) {
                float bluePixel = 0;
                float greenPixel = 0;
                float redPixel = 0;

                int count = -1;
                int pixelVal;

                int rangeToBeConvoluted = radius * width;
                for(int k = j - rangeToBeConvoluted; k <= j + rangeToBeConvoluted; k = k + width) {
                    count++;
                    // TODO: Must check the time taken to process
                    try{
                        pixelVal = pixelsIn[k];
                    }catch (ArrayIndexOutOfBoundsException a){
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
    private static void performHorizontalConvolutionWithGivenSigma(int top, int[] pixelsIn, int[] pixelsOut, int width, int totalPixels, int radius, double[] gaussianKernelVector){
        for(int i=top;i<totalPixels;i=i+width){
            int pixelRight = i + width;
            for(int j = i; j < pixelRight; j++) {
                float bluePixel = 0;
                float greenPixel = 0;
                float redPixel = 0;

                int pixelVal;

                for(int k = -radius; k <= radius; k++) {
                    int pixelIndex = j + k;
                    int vectorIndex = k + radius;

                    // TODO: Must check the time taken to process
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


    private static void performConvolution(int top, int bottom, int[] pixelsIn, int[] pixelsIntermediate, int[] pixelsOut, int height, int width, int totalPixels, float sigma, boolean isSigmaFar, boolean singleSigma){
        if(singleSigma){
            double[] gaussianKernelVector = constructGaussianKernel(sigma);
            performVerticalConvolutionWithGivenSigma(top, pixelsIn, pixelsIntermediate, width, totalPixels, gaussianKernelVector.length/2,gaussianKernelVector);
            performHorizontalConvolutionWithGivenSigma(top, pixelsIntermediate, pixelsOut, width, totalPixels, gaussianKernelVector.length/2,gaussianKernelVector);
        }else {
            double[][] gaussianKernelVector = new double[height][];

            int lowIndex = top/width;
            int highIndex = bottom/width;
            for(int i=lowIndex;i<highIndex;i++){
                gaussianKernelVector[i] = constructGaussianKernel(calculateSigma(lowIndex,highIndex,i,sigma,isSigmaFar));
            }

            performVerticalConvolution(top, bottom, pixelsIn, pixelsIntermediate,width,gaussianKernelVector);
            performHorizontalConvolution(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, gaussianKernelVector);
        }
    }
    public static Bitmap tiltshift_java(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);

        /* Image dimensions */
        int imageWidth = input.getWidth();
        int imageHeight = input.getHeight();
        int totalPixels = input.getWidth()*input.getHeight();

        int[] pixelsIn = new int[totalPixels];
        int[] pixelsIntermediate  = new int[totalPixels];
        int[] pixelsOut = new int[totalPixels];

        System.out.println(a0+"  "+a1+"  "+a2+"   "+a3+" "+sigma_far+"  "+sigma_near);

        input.getPixels(pixelsIn,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        performConvolution(0,(a0-1)*imageWidth,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_far,true, true);
        System.out.println("**** Finished Processing First Part ****");

        performConvolution(a0*imageWidth,a1*imageWidth,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_far,true, false);
        System.out.println("**** Finished Processing Second Part ****");

        for(int i=(a1+1)*imageWidth; i<(a2-1)*imageWidth; i++){
            pixelsOut[i] = pixelsIn[i];
        }
        System.out.println("**** Finished Processing Third Part ****");

        performConvolution(a2*imageWidth,a3*imageWidth,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_near,false, false);
        System.out.println("**** Finished PRocessing Fourth Part ****");

        performConvolution((a3+1)*imageWidth,imageHeight-2,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_near,false, true);
        System.out.println("--------- Processing Done ------------");

        outBmp.setPixels(pixelsOut,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        return outBmp;
    }
    public static Bitmap tiltshift_cpp(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        tiltshiftcppnative(pixels,pixelsOut,input.getWidth(),input.getHeight(),sigma_far,sigma_near,a0,a1,a2,a3);

        outBmp.setPixels(pixelsOut,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());
        return outBmp;
    }
    public static Bitmap tiltshift_neon(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        tiltshiftneonnative(pixels,pixelsOut,input.getWidth(),input.getHeight(),sigma_far,sigma_near,a0,a1,a2,a3);

        outBmp.setPixels(pixelsOut,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());
        return outBmp;
    }


    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
    public static native int tiltshiftcppnative(int[] inputPixels, int[] outputPixels, int width, int height, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3);
    public static native int tiltshiftneonnative(int[] inputPixels, int[] outputPixels, int width, int height, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3);

}
