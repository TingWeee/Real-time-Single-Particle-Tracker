/*
 * codes containing 2nd algorightm (stand alone)
 */
package RoughLocalisation;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import javax.swing.SwingWorker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import twodgaussfunction.TwoDGaussFunction;
import twodgaussianplugin.TwoDGaussianPlugin;
import twodgaussianplugin.TwoDGaussianPlugin.doTwoFit;
import static twodgaussianplugin.TwoDGaussianPlugin.printlog;


/**
 *
 * @author ting
 */
public class SpotDectector {
    //find local maxima in an image to pixel level accuracy
    //then do 2d guassian fitting and start tracking
    // free parameters: particle size, offset, size of fitting box 
    
    //image to process,image should be bright spots on dark background
    static ImagePlus imp = new ImagePlus(); 
    //if i put static ImagePlus the whole damn thing doesnt work
    
    //the minimum brightness of a pixel that might be local maxima
    static int threshold; 
   
    //if data is noisy (e.g single particle has multiple maxima), then set this optional
    static int size=8; 
    
    int framNum;
   
    //
    static ArrayList<int[]> brightspots = new ArrayList<int[]>();
    static ArrayList<int[]> brightspotsFinal = new ArrayList<int[]>();
    
    //constructor
    public SpotDectector(ImagePlus imp) {
        this.imp= imp;
    }
    
    //finding threshold as the average count value 
    //putting this as a separate function so that users can choose threshold easily 
    public static void findThreshold(ImagePlus imp,int framNum){
        double sum =0;
        int w = imp.getWidth();
        int h = imp.getHeight();
        double count= 0;
        for (int x=0; x<w;x++){
            for (int y=0;y<h;y++){
                count = imp.getStack().getProcessor(framNum).get(x,y);
                sum+=count;
            }
        }
        threshold = (int)(sum/(w*h)*3) ;
    }
        
    //find all pixels above threshold
    //check each pixel to see it its brighter than its neighbours
    //if its not brighter, remove it 
    public static void findPixels(ImagePlus imp, int framNum){
        int w = imp.getWidth();
        int h = imp.getHeight();
        //double[] countArray = new double[w*h];
        
        
        int count;
        int count1;
        int count2;
        int count3;
        int count4;
        //Toremove
        int totalsuccess = 0;
        for (int x=1; x<w-1;x++){
            for (int y=1;y<h-1;y++){
          
                count = imp.getStack().getProcessor(framNum).get(x,y);
                count1 = imp.getStack().getProcessor(framNum).get(x+1,y);
                count2 = imp.getStack().getProcessor(framNum).get(x,y+1);
                count3 = imp.getStack().getProcessor(framNum).get(x-1,y);
                count4 = imp.getStack().getProcessor(framNum).get(x,y-1);

                if (count>threshold && count>count1 && count>count2 && count>count3 && count>count4){
                   
                    totalsuccess++;
                    int[] temp = new int[2];
                    temp[0]=x;
                    temp[1]=y;
                    
                    brightspots.add(temp);
                }
                
            }
        }
        brightspotsFinal.addAll(brightspots);
        removeRedundant();
    }
    //remove doublecounted values
    private static void removeRedundant(){
        int initialNum = brightspots.size();//this give me twice the size, major oof
        initialNum=initialNum/2;
        for (int i =0; i<initialNum-1; i++){
            boolean doubleCounted = false;
            for (int w = i+1; w<initialNum; w++){
                if (Math.sqrt((Math.pow(brightspots.get(w)[0]-brightspots.get(i)[0], 2))+(Math.pow(brightspots.get(w)[1]-brightspots.get(i)[1], 2)))<size)
                    //brightspotsFinal.remove(brightspots.get(w));
                    doubleCounted = true;
                
            }
            if (doubleCounted){
                brightspotsFinal.remove(brightspots.get(i));
            }
        }
    }
    
    static int numParticles;
    static int numFrames;//can change this to the amount of frames that the live prog feeds in
    //global 3D array to store result
    static double[][][] finalResult; 
    //multithreading
    
    public static void loopParticles (){
        
        IJ.log("Main thread: " + Thread.currentThread().getName());
        findThreshold(imp, 1);
        findPixels(imp, 1);
        numFrames = imp.getStackSize();
        numParticles = brightspotsFinal.size();
        
        brightspotsFinal.forEach(i -> System.out.println(Arrays.toString(i)));
        
        /* to print out X and Y coordinates of rough localisation 
        int [] xarr = new int[numParticles];
        int [] yarr = new int[numParticles];
        for (int i=0; i<numParticles; i++){
            xarr[i]=brightspotsFinal.get(i)[0];
            yarr[i]=brightspotsFinal.get(i)[1];
        }
        System.out.println(Arrays.toString(xarr));
        System.out.println(Arrays.toString(yarr));
        */
  
        finalResult = new double[numParticles][numFrames][2];
        
        for (int i = 0; i<numParticles;i++){
            
            startTrack2(brightspotsFinal.get(i),i);
        }
        IJ.log("numP"+numParticles);
        IJ.log("finish starttrack");
        //
        
        
    }
    
    static int ROI_width_2 = 20;
    static int ROI_height_2 = 20;
    
    public static void startTrack2(int[] coor, int index){ //coordinates of(x,y) of object taken in 
        SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>(){
            @Override
            protected Void doInBackground() throws Exception{
       
                IJ.log("StartTrack2 thread ID: " + Thread.currentThread().getName());
                int label = index; //label the particle, labelling started from 0 for convenience
                
                //initialise the roi for first frame
                
                Roi roi1 =new Roi(Math.max(coor[0]-ROI_width_2/2,0), Math.max(coor[1]-ROI_height_2/2,0), ROI_width_2, ROI_height_2);
                imp.setRoi(roi1); //check if this works in multithread..
                
                //other parameters
                int shiftR_2 = 3; //200/110 theen i just round up will include some formulas next time 
                int frameno_2 = 1;
                boolean iterate_2 = false;
                int incrementX_2 = 0;//if increment not zero, go to while loop and set newStart_2[1] =0
                int incrementY_2 = 0;//if increment not zero, go to while loop and set newStart_2[2] =0

                double[] data_2= new double[ROI_width_2*ROI_height_2]; //flattened 2d array 
                double[] newStart_2 = { //set initial parameters
                    100,//amplitude arbitrary, roughly //
                    1,
                    1,
                    1,//depends on camera pixel size and magnification
                    1,
                    400 // find average of all cnts in one frame
                };
                double[] optimalValues_2 = new double[6]; //FOR??? why is it its own class
                int[] optim_param_2;//for fiitting, dont need know
                int data_width_2;//assuming its a sqaure region, width of square in pixel
          
                newStart_2[1] = ROI_width_2/ 2;
                newStart_2[2] = ROI_height_2/ 2;
                
                fitMeDialogue(frameno_2, roi1, data_2, newStart_2, label, optimalValues_2);
                

               
                for (int framNum_ =frameno_2+1; framNum_ <=numFrames; framNum_++){
                    //insert comparator 
                    //if (overlap(framNum_-1, label)){
                    //    System.out.println("stopped tracking particle: "+ label);
                    //    break;
                    //}
                    Integer[] newres;
                    for (int i=0; i<=4; i++){

                        newStart_2[1]= ROI_width_2/2;
                        newStart_2[2]=ROI_height_2/2;
                        if(i==0){
                            //uses meanX and meanY of previous frame to reposition the ROI
                            newres = getNewRoiPosition(optimalValues_2[1], optimalValues_2[2], 0, 0);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi1.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi1.getBounds().getX() + ", roi y: " + roi1.getBounds().getY() + ", roi w: " + roi1.getBounds().getWidth() + ", roi h: " + roi1.getBounds().getHeight());

                            if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                break;
                            else if (iterate_2){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration();
                                for (int[] x:combination){
                                    newStart_2[1]=x[0];
                                    newStart_2[2]=x[1];
                                    if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2));
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }
                        }else if(i==1){
                            newres = getNewRoiPosition(optimalValues_2[1], optimalValues_2[2], shiftR_2, 0);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi1.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi1.getBounds().getX() + ", roi y: " + roi1.getBounds().getY() + ", roi w: " + roi1.getBounds().getWidth() + ", roi h: " + roi1.getBounds().getHeight());

                            if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                break;
                            else if (iterate_2){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration();
                                for (int[] x:combination){
                                    newStart_2[1]=x[0];newStart_2[2]=x[1];
                                    if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2));
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }

                        }else if(i==2){
                            newres = getNewRoiPosition(optimalValues_2[1], optimalValues_2[2], -1*shiftR_2, 0);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi1.setLocation(newres[0], newres[1]);//WHAT DOES THIS DO?
                            printlog("roi x: " + roi1.getBounds().getX() + ", roi y: " + roi1.getBounds().getY() + ", roi w: " + roi1.getBounds().getWidth() + ", roi h: " + roi1.getBounds().getHeight());
                            if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                break;
                            else if (iterate_2){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration();
                                for (int[] x:combination){
                                    newStart_2[1]=x[0];newStart_2[2]=x[1];
                                    if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }
                        }else if(i==3){
                            newres = getNewRoiPosition(optimalValues_2[1], optimalValues_2[2], 0, shiftR_2);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi1.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi1.getBounds().getX() + ", roi y: " + roi1.getBounds().getY() + ", roi w: " + roi1.getBounds().getWidth() + ", roi h: " + roi1.getBounds().getHeight());
                            if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                break;
                            else if (iterate_2){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration();
                                for (int[] x:combination){
                                    newStart_2[1]=x[0];newStart_2[2]=x[1];
                                    if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }
                        }else{
                            newres = getNewRoiPosition(optimalValues_2[1], optimalValues_2[2], 0, -1*shiftR_2);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi1.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi1.getBounds().getX() + ", roi y: " + roi1.getBounds().getY() + ", roi w: " + roi1.getBounds().getWidth() + ", roi h: " + roi1.getBounds().getHeight());
                            if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                break;
                            else if (iterate_2){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration();
                                for (int[] x:combination){
                                    newStart_2[1]=x[0];newStart_2[2]=x[1];
                                    if(fitMeDialogue(framNum_, roi1, data_2, newStart_2, label, optimalValues_2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }

                        }
                    }
                }
            
            return null;  
            }
        };
        
        worker.execute();
    }
    
    private static boolean overlap(int prevFrame, int label){//prevFrame-2
        for (int x=0; x<numParticles;x++){
            if (x !=label){
                if (Math.sqrt((Math.pow((finalResult[x][prevFrame-1][0])-(finalResult[label][prevFrame-1][0]), 2))+(Math.pow((finalResult[x][prevFrame-1][1])-(finalResult[label][prevFrame-1][1]), 2)))<size)
                    return true;
            }
            
        }
        return false;
    }
    
    private static int[][] WayOfIteration() {
        int[][] combination = new int[4][2];
        combination[0][0] = ROI_width_2 / 4;
        combination[0][1] = ROI_height_2 / 4;
        combination[1][0] = -1 * ROI_width_2 / 4;
        combination[1][1] = -1 * ROI_height_2 / 4;
        combination[2][0] = ROI_width_2 / 4;
        combination[2][1] = -1 * ROI_height_2 / 4;
        combination[3][0] = -1 * ROI_width_2 / 4;
        combination[3][1] = ROI_height_2 / 4;
        return combination;
    }

    private static double[] get1Darray(ImagePlus imp, int frameindex, Roi roi) { // index start at 1
        printlog("inside get1Darray framindex: " + frameindex + ", roix: " + roi.getBounds().getX() + ", roiy: " + roi.getBounds().getY() + ", roiwidth: " + roi.getBounds().getWidth() + ", roiheight: " + roi.getBounds().getHeight());

        Rectangle rect = roi.getBounds();
        int l_disp = (int) rect.getX(); //i changed these into global variables
        int t_disp = (int) rect.getY();
        int w_disp = (int) rect.getWidth();
        int h_disp = (int) rect.getHeight();
        int f = frameindex;

        double[] res = new double[w_disp * h_disp];

        for (int y = 0; y < h_disp; y++) {
            for (int x = 0; x < w_disp; x++) {
                res[x + (y * w_disp)] = imp.getStack().getProcessor(f).get(x + l_disp, y + t_disp);
                //GET THE COUNTS
                //have to remap since coordinates are in terms of 
            }
        }
        printlog("2inside get1Darray framindex: " + frameindex + ", roix: " + roi.getBounds().getX() + ", roiy: " + roi.getBounds().getY() + ", roiwidth: " + roi.getBounds().getWidth() + ", roiheight: " + roi.getBounds().getHeight());

        return res;//1D ARRAY
    }

    private static Integer[] getNewRoiPosition(double fittedX, double fittedY, int shiftx, int shifty) {
        // index start 0 fittedX,fittedY
        //index start 1 initW, initH
        Integer[] result = new Integer[2];
        int initH = ROI_height_2;
        int initW = ROI_width_2;

        int new_left = (int) (Math.round(fittedX - (initW / 2) + shiftx));
        int new_top = (int) (Math.round(fittedY - (initH / 2) + shifty));

        if (new_left < 0) {
            new_left = 0;
        }

        if (new_top < 0) {
            new_top = 0;
        }
        result[0] = new_left;
        result[1] = new_top;

        return result;
    }

    private static boolean fitMeDialogue(int framNum, Roi roi1, double[] data_2, double[] newStart_2, int label, double[] optimalValues_2){
        data_2 = get1Darray(imp, framNum, roi1);
        //data_width_2= ROI_width_2;
        doTwoFit_2 x = new doTwoFit_2(data_2, newStart_2, ROI_width_2, new int[]{1000, 100});
        try {//do LevenbergMarquardt optimization and get optimized parameters
            LeastSquaresOptimizer.Optimum opt = x.fit2dGauss(); //result store here

            final double[] optimalValues = opt.getPoint().toArray();
            printlog("optimalalues[1]: " + optimalValues[1]);
            printlog("optimalalues[2]: " + optimalValues[2]);
            printlog("roi x: " + roi1.getBounds().getX() + ", roi y: " + roi1.getBounds().getY() + ", roi w: " + roi1.getBounds().getWidth() + ", roi h: " + roi1.getBounds().getHeight());

            optimalValues[1] = optimalValues[1] + roi1.getBounds().getX(); 
            optimalValues[2] = optimalValues[2] + roi1.getBounds().getY();
            
            optimalValues_2[1]=optimalValues[1];
            optimalValues_2[2]=optimalValues[2];
            

            printlog("amplitude" + optimalValues[0] + "meanx: " + optimalValues[1] + " meany: " + optimalValues[2] + " sigmax: " + optimalValues[3] + " sigma y: " + optimalValues[4] + " offset: " + optimalValues[5]);

            IJ.log("Thread " + Thread.currentThread().getName() + " framNum: " + framNum);
            finalResult[label][framNum-1][0] = optimalValues[1];//changesmade
    //        tfMeanY.setText(Double.toString(optimalValues[2]));
            finalResult[label][framNum-1][1] = optimalValues[2];//changesmade
    //        tfSigmaX.setText(Double.toString(optimalValues[3]));
    //        tfSigmaY.setText(Double.toString(optimalValues[4]));
    //        tfOffset.setText(Double.toString(optimalValues[5]));
            optimalValues_2 = optimalValues;
            printlog("optimalalues_1[1]: " + optimalValues_2[1]);
            printlog("optimalalues_1[2]: " + optimalValues_2[2]);

            //output data
//            System.out.println("v0: " + optimalValues[0]);
//            System.out.println("v1: " + optimalValues[1]);
//            System.out.println("v2: " + optimalValues[2]);
//            System.out.println("v3: " + optimalValues[3]);
//            System.out.println("v4: " + optimalValues[4]);
//            System.out.println("v5: " + optimalValues[5]);
//            System.out.println("Iteration number: " + opt.getIterations());
//            System.out.println("Evaluation number: " + opt.getEvaluations());
            return true;
        } catch (Exception e) {
    //        tfAmplitude.setText(Double.toString(Double.NaN));
    //        tfMeanX.setText(Double.toString(Double.NaN));
    //        tfMeanY.setText(Double.toString(Double.NaN));
    //        tfSigmaX.setText(Double.toString(Double.NaN));
    //        tfSigmaY.setText(Double.toString(Double.NaN));
    //        tfOffset.setText(Double.toString(Double.NaN));
            System.out.println(e.toString());
            printlog("-----First frame fit NOK");
            return false;

        }
    }


    public static class doTwoFit_2 {
        double [] data;
        double[] newStart;
        int data_width;
        int[] optim_param;
        
        private static LeastSquaresProblem lsp;//not sure if putting here is appropriate
        
        public doTwoFit_2(double[] data, double[] newStart, int data_width, int[] optim_param) { //takes 4 arg and stores in global variables
            this.data = data;//this.data = data;
            this.newStart= newStart;
            this.data_width = data_width;
            this.optim_param = optim_param;
            buildlsb();//the one that does the fitting
        }
        /**
         * build LeastSquareProblem by using constructor data
         */
        private void buildlsb() {
            //construct two-dimensional Gaussian function
            TwoDGaussFunction tdgf = new TwoDGaussFunction(data_width, data.length);

            //prepare construction of LeastSquresProblem by builder
            LeastSquaresBuilder lsb = new LeastSquaresBuilder();

            //set model function and its jacobian
            lsb.model(tdgf.retMVF(), tdgf.retMMF());
            //set target data
            lsb.target(data);
            //set initial parameters
            lsb.start(newStart);
            //set upper limit of evaluation time
            lsb.maxEvaluations(optim_param[0]);
            //set upper limit of iteration time
            lsb.maxIterations(optim_param[1]);

            lsp = lsb.build(); // TODOmultithreading: pay attention here, make sure each individual thread has its own lsp
        }

        /**
         * Do two dimensional Gaussian fit
         *
         * @return return the fitted data as Optimum
         */
        public LeastSquaresOptimizer.Optimum fit2dGauss() {
            LevenbergMarquardtOptimizer lmo = new LevenbergMarquardtOptimizer();
            LeastSquaresOptimizer.Optimum lsoo = lmo.optimize(lsp);

            return lsoo;
        }
    }
    
    public static void printResult() throws IOException{
        int i =0;
        
        for(double[][] particle: finalResult){
            i++;
            int f =0;
            for (double[] frame: particle){
                f++;
                System.out.println("particle: "+i+" frame: "+f+" x: "+frame[0]+" y: "+frame[1]);
                writeFile();
            }
        }
        
    }
    
    //to save file in csv format
    public static void writeFile() throws IOException{
        StringBuilder builder = new StringBuilder();
        for(int i =0;i<finalResult.length;i++){
            for (int j =0;j<finalResult[i].length;j++){
                builder.append(Arrays.toString(finalResult[i][j])+"");
                if (j<finalResult.length-1)
                    builder.append(" , ");
            }
            builder.append("\n");
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter("sampledata.txt"));
        writer.write(builder.toString());
        writer.close();
       
    }
    
    
   
   
}

        /*

        findThreshold(imp, 1);
        findPixels(imp, 1, threshold);

        */





