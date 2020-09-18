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


        for(int k=-wholeRadius ;k<=wholeRadius;k++){
            double secondTerm = -1* k*k/(2*sigmaSquare);
            double weight = Math.exp(secondTerm)*firstTerm;

            kernelVector[k+wholeRadius] = weight;
        }

        return kernelVector;
    }
    public static Bitmap tiltshift_java(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        //cannot write to input Bitmap, since it may be immutable
        //if you try, you may get a java.lang.IllegalStateException
        System.out.println(sigma_far+ " "+sigma_near+" "+a0 + " "+a1+ " "+a2+ " "+a3);

        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsIntermediate  = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());


        float sigma = 3F;
        double radius = Math.ceil(2*sigma);
        int r = (int)radius;
        System.out.println("Radius  " + r);
        double[] gaussianKernelVector = constructGaussianKernel(r,sigma);

        int imageWidth = input.getWidth();
        int imageHeight = input.getHeight();
        int totalPixels = input.getWidth()*input.getHeight();
        int gaussianVectorLength = gaussianKernelVector.length;

        System.out.println("Image Width and Height "+ imageWidth +"   "+imageHeight +" "+gaussianKernelVector.length);

        for(double x: gaussianKernelVector){
            System.out.print(x+" ");
        }

        for (int i=0; i<imageWidth; i++){
            for(int j = i; j < totalPixels; j=j+imageWidth) {
                int color = 0;
                float bluePixel = 0, greenPixel = 0, redPixel = 0, alphaPixel = 0;
                int count = -1;
                int px;
                for(int k = j - (r * imageWidth); k <= j + (r * imageWidth); k = k + imageWidth) {
                    count++;
                    if(k < i || k > totalPixels -1) {
                        px = 1;
                    }else{
                        px = pixels[k];
                    }
                    int B = px % 0xff;
                    int G = (px >> 8) % 0xff;
                    int R = (px >> 16) % 0xff;
                    int A = (px >> 24) % 0xff;

                    redPixel += (gaussianKernelVector[count] * R );
                    greenPixel += (gaussianKernelVector[count] * G );
                    bluePixel += (gaussianKernelVector[count] * B );
                    alphaPixel += (gaussianKernelVector[count] * A);
                }

                int ap = (int) alphaPixel;
                int rp = (int) redPixel;
                int gp = (int) greenPixel;
                int bp = (int) bluePixel;

                color = (ap  << 24) | (rp << 16) | (gp << 8) | bp;
                pixelsIntermediate[j] = (int) color;
            }
        }
//
        for (int i=0; i<imageHeight; i++){
            int pixelLeft = i*imageWidth;
            int pixelRight = pixelLeft + imageWidth-1;

            for(int j=pixelLeft; j<pixelRight && j < totalPixels; j++) {
                int color = 0;
                float bluePixel = 0, greenPixel = 0, redPixel = 0, alphaPixel = 0;
                int px;
                for(int k=-r;k<=r;k++) {
                    if(k < 0 || k > pixelRight-1 || k > totalPixels){
                        px = 1;
                    }else{
                        px = pixelsIntermediate[k];
                    }
                    int B = px % 0xff;
                    int G = (px >> 8) % 0xff;
                    int R = (px >> 16) % 0xff;
                    int A = (px >> 24) % 0xff;

                    redPixel += (gaussianKernelVector[k + r] * R);
                    greenPixel += (gaussianKernelVector[k + r] * G);
                    bluePixel += (gaussianKernelVector[k +r] * B);
                    alphaPixel += (gaussianKernelVector[k + r] * A);
                }

                int ap = (int) alphaPixel;
                int rp = (int) redPixel;
                int gp = (int) greenPixel;
                int bp = (int) bluePixel;

                color = (ap  << 24) | (rp << 16) | (gp << 8) | bp;
                pixelsOut[j] = (int)color;

            }
        }

        outBmp.setPixels(pixelsIntermediate,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        System.out.println("Something");
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
