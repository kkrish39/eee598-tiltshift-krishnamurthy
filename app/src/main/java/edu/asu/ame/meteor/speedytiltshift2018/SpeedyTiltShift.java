package edu.asu.ame.meteor.speedytiltshift2018;

import android.graphics.Bitmap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;

public class SpeedyTiltShift {
    static SpeedyTiltShift Singleton = new SpeedyTiltShift();

    private static final Logger logger = Logger.getLogger(String.valueOf(SpeedyTiltShift.class));

    static {
        System.loadLibrary("native-lib");
    }


    private static double[] constructGaussianKernel(double sigma){
        if(sigma < 0.6)
            return new double[0];
        double radius = Math.ceil(2*sigma);

        int kernelSize = (int)(Math.ceil(radius)*2) + 1;
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
        return  sigma * (isFarSigma ? (high - y) : (y - low))/(high - low);
    }

    private static void performVerticalConvolution(int top, int bottom, int[] pixelsIn, int[] pixelsOut, int width, double[][] gaussianKernelVectors){
        for (int i=top; i<=top+width; i++){
            for(int j = i; j < bottom; j=j+width) {
                float bluePixel = 0;
                float greenPixel = 0;
                float redPixel = 0;

                int count = -1;
                int pixelVal;

                double[] gaussianKernelVector = gaussianKernelVectors[j/width];

                if(gaussianKernelVector.length == 0){
                    pixelsOut[j] = pixelsIn[j];
                    continue;
                }
                int radius = gaussianKernelVector.length/2;
                int rangeToBeConvoluted = radius * width;

                for(int k = j - rangeToBeConvoluted; k <= j + rangeToBeConvoluted; k = k + width) {
                    count++;
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
    private static void performHorizontalConvolution(int top, int bottom, int[] pixelsIn, int[] pixelsOut, int width, int totalPixels, double[][] gaussianKernelVectors){
        for(int i=top;i<bottom;i=i+width){
            int pixelRight = i + width - 1;
            double[] gaussianKernelVector = gaussianKernelVectors[i/width];

            if(gaussianKernelVector.length == 0){
                if (pixelRight - i >= 0)
                    System.arraycopy(pixelsIn, i, pixelsOut, i, pixelRight - i+1);
                continue;
            }

            int radius = gaussianKernelVector.length/2;
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
    private static void performVerticalConvolutionWithGivenSigma(int top, int bottom, int[] pixelsIn, int[] pixelsOut, int width, int totalPixels, int radius, double[] gaussianKernelVector){
        for (int i=top; i<top+width; i++){
            for(int j = i; j < bottom; j=j+width) {
                float bluePixel = 0;
                float greenPixel = 0;
                float redPixel = 0;

                int count = -1;
                int pixelVal;

                int rangeToBeConvoluted = radius * width;
                for(int k = j - rangeToBeConvoluted; k <= j + rangeToBeConvoluted; k = k + width) {
                    count++;
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
    private static void performHorizontalConvolutionWithGivenSigma(int top, int bottom, int[] pixelsIn, int[] pixelsOut, int width, int totalPixels, int radius, double[] gaussianKernelVector){
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


    private static void performConvolution(int top, int bottom, int[] pixelsIn, int[] pixelsIntermediate, int[] pixelsOut, int height, int width, int totalPixels, float sigma, boolean isSigmaFar, boolean singleSigma){
        if(singleSigma){
            double[] gaussianKernelVector = constructGaussianKernel(sigma);
            if(gaussianKernelVector.length == 0){
                logger.log(Level.WARNING, "Gaussian Kernel Vector cannot be empty");
                return;
            }
            performVerticalConvolutionWithGivenSigma(top, bottom, pixelsIn, pixelsIntermediate, width, totalPixels, gaussianKernelVector.length/2,gaussianKernelVector);
            performHorizontalConvolutionWithGivenSigma(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, gaussianKernelVector.length/2,gaussianKernelVector);
        }else {
            double[][] gaussianKernelVector = new double[height+1][];

            int lowIndex = top/width;
            int highIndex = bottom/width;

            for(int i=lowIndex;i<=highIndex;i++){
                gaussianKernelVector[i] = constructGaussianKernel(calculateSigma(lowIndex,highIndex,i,sigma,isSigmaFar));
            }

            performVerticalConvolution(top, bottom, pixelsIn, pixelsIntermediate,width,gaussianKernelVector);
            performHorizontalConvolution(top, bottom, pixelsIntermediate, pixelsOut, width, totalPixels, gaussianKernelVector);
        }
    }
    public static Bitmap tiltshift_java(Bitmap input, final float sigma_far, final float sigma_near, final int a0, final int a1, final int a2, final int a3) throws InterruptedException {
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);

        /* Image dimensions */
        final int imageWidth = input.getWidth();
        final int imageHeight = input.getHeight();
        final int totalPixels = input.getWidth()*input.getHeight();

        /* Arrays to keep track of image pixels */
        final int[] pixelsIn = new int[totalPixels];
        final int[] pixelsIntermediate  = new int[totalPixels];
        final int[] pixelsOut = new int[totalPixels];

        input.getPixels(pixelsIn,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        Callable<Void> farSigmaThread = new Callable<Void>()
        {
            @Override
            public Void call() throws Exception
            {
                performConvolution(0,(a0)*imageWidth,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_far,true, true);
                logger.log(Level.INFO, "farSigmaThread Completed Execution...");
                return null;
            }
        };

        Callable<Void> a0a1Thread = new Callable<Void>()
        {
            @Override
            public Void call() throws Exception
            {
                performConvolution(a0*imageWidth,a1*imageWidth,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_far,true, false);
                logger.log(Level.INFO, "a0a1Thread Completed Execution...");
                return null;
            }
        };

        Callable<Void> noBlurThread = new Callable<Void>()
        {
            @Override
            public Void call() throws Exception
            {

                if ((a2) * imageWidth - (a1) * imageWidth >= 0)
                    System.arraycopy(pixelsIn, (a1) * imageWidth, pixelsOut, (a1) * imageWidth, (a2) * imageWidth - (a1) * imageWidth);

                logger.log(Level.INFO, "noBlurThread Completed Execution...");
                return null;
            }
        };

        Callable<Void> a2a3Thread = new Callable<Void>()
        {
            @Override
            public Void call() throws Exception
            {
                performConvolution(a2*imageWidth,a3*imageWidth,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_near,false, false);
                logger.log(Level.INFO, "a2a3Thread Completed Execution...");
                return null;
            }
        };

        Callable<Void> nearSigmaThread = new Callable<Void>()
        {
            @Override
            public Void call() throws Exception
            {
                performConvolution((a3)*imageWidth,imageHeight*imageWidth,pixelsIn,pixelsIntermediate,pixelsOut,imageHeight,imageWidth,totalPixels,sigma_near,false, true);
                logger.log(Level.INFO, "nearSigmaThread Completed Execution...");
                return null;
            }
        };


        List<Callable<Void>> operationList = new ArrayList<>(Arrays.asList(farSigmaThread,a0a1Thread,noBlurThread, a2a3Thread,nearSigmaThread));

        ExecutorService ex = Executors.newFixedThreadPool(5);

        try{
            ex.invokeAll(operationList);
        }catch (InterruptedException e){
            logger.log(Level.SEVERE, "Thread preempted. Failure!!!");
            Thread.currentThread().interrupt();
        }

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
