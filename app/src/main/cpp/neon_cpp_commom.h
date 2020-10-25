

/*structure to send the arguments to the function from the thread*/
struct threadArgs{
    int top;
    int bottom;
    int *pixelsIn;
    int *pixelsIntermediate;
    int *pixelsOut;
    int width;
    int totalPixels;
    float sigma;
    bool isSigmaFar;
    bool singleSigma;
};