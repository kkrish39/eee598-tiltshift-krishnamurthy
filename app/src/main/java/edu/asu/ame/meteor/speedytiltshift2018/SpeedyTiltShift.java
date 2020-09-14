package edu.asu.ame.meteor.speedytiltshift2018;

import android.graphics.Bitmap;

public class SpeedyTiltShift {
    static SpeedyTiltShift Singleton = new SpeedyTiltShift();

    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }


    public static double[] constructGaussianKernel(float radius, float sigma){
        int kernelSize = (int)(Math.ceil(radius)*2) + 1;

        double[] kernelVector = new double[kernelSize];

        int wholeRadius = (int)Math.ceil(radius);
        double sigmaSquare = sigma * sigma;
        double twoPiSigmaSquare = 2 * Math.PI * sigmaSquare;
        double sqrtTwoPiSigmaSquare = Math.sqrt(twoPiSigmaSquare);

        double firstTerm = 1/(sqrtTwoPiSigmaSquare);


        for(int i=-wholeRadius ;i<=wholeRadius;i++){
            double secondTerm = -1* i*i/(2*sigmaSquare);
            double weight = Math.exp(secondTerm)*firstTerm;

            kernelVector[i+wholeRadius] = weight;
        }

        return kernelVector;
    }
    public static Bitmap tiltshift_java(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        //cannot write to input Bitmap, since it may be immutable
        //if you try, you may get a java.lang.IllegalStateException

        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsIntermediate  = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        double[] gaussianKernelVector = constructGaussianKernel(2,0.8F);

        int imageWidth = input.getWidth();
        int imageHeight = input.getHeight();
        int totalPixels = input.getWidth()*input.getHeight();
        int gaussianVectorLength = gaussianKernelVector.length;

        System.out.println("Image Width and Height "+ imageWidth +"   "+imageHeight);
        for (int i=0; i<imageHeight; i++){
            int pixelLeft = i*imageWidth;
            int pixelRight = pixelLeft + imageWidth-1;


            for(int j=pixelLeft; j<pixelRight && j < totalPixels; j++) {
                int color = 0;
                float bluePixel = 0, greenPixel = 0, redPixel = 0;
//                System.out.println("Pixel to be mainpulated ---->" + pixels[j]);
                for(int k=j-2;k<=j+2;k++) {
                    if(k < 0 || k > pixelRight-1 || k > totalPixels) continue;
                    int B = pixels[k] % 0x100;
                    int G = (pixels[k] >> 8) % 0x100;
                    int R = (pixels[k] >> 16) % 0x100;
                    int A = 0xff;

//                    A += (int)(gaussianKernelVector[k % gaussianVectorLength] * A);
                    redPixel += gaussianKernelVector[k % gaussianVectorLength] * ((R >> 16)& 0xff);
                    greenPixel += gaussianKernelVector[k % gaussianVectorLength] * ((G >> 8)& 0xff);
                    bluePixel += (gaussianKernelVector[k % gaussianVectorLength] * (B & 0xff));
                }

                int ap = 0xff;
                int rp = (int) redPixel;
                int gp = (int) greenPixel;
                int bp = (int) bluePixel;
//                System.out.println(rp +"   "+gp+"    "+bp);
                color = (ap  << 24) | (rp << 16) | (gp << 8) | bp;
                pixelsIntermediate[j] = (int)color;
//                System.out.println(pixelsIntermediate[j]);
            }


        }
        outBmp.setPixels(pixelsIntermediate,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

//        for (int i=0; i<imageWidth; i++){
//            int pixelTop = i*imageHeight;
//
//
//            for(int j=pixelTop; j < totalPixels; j=j+imageWidth) {
//                int color = 0;
//                float bluePixel = 0, greenPixel = 0, redPixel = 0;
////                System.out.println("Pixel to be mainpulated ---->" + pixels[j]);
//                for(int k=j-2;k<=j+2;k++) {
//                    if(k < pixelTop || k > pixelBottom-1) continue;
//                    int B = pixelsIntermediate[k] % 0x100;
//                    int G = (pixelsIntermediate[k] >> 8) % 0x100;
//                    int R = (pixelsIntermediate[k] >> 16) % 0x100;
//                    int A = 0xff;
//
////                    A += (int)(gaussianKernelVector[k % gaussianVectorLength] * A);
//                    redPixel += gaussianKernelVector[k % gaussianVectorLength] * ((R >> 16)& 0xff);
//                    greenPixel += gaussianKernelVector[k % gaussianVectorLength] * ((G >> 8)& 0xff);
//                    bluePixel += (gaussianKernelVector[k % gaussianVectorLength] * (B & 0xff));
//                }
//
//                int ap = 0xff;
//                int rp = (int) redPixel;
//                int gp = (int) greenPixel;
//                int bp = (int) bluePixel;
////                System.out.println(rp +"   "+gp+"    "+bp);
//                color = (ap  << 24) | (rp << 16) | (gp << 8) | bp;
//                pixelsOut[j] = (int)color;
////                System.out.println(pixelsIntermediate[j]);
//            }
//
//
//        }
//        outBmp.setPixels(pixelsOut,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());
//
//        System.out.println("Someting");
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
