package twodgaussianplugin;

//import com.sun.org.apache.xml.internal.security.utils.XalanXPathAPI;
import twodgaussfunction.TwoDGaussFunction;
import RoughLocalisation.SpotDectector;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Roi;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.naming.spi.DirStateFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;

/**
 *
 * @author TingWei
 */
public class TwoDGaussianPlugin {
   
    static ImagePlus imp_ = new ImagePlus();
    static ImageCanvas imp_can;
    static int totalFrame_;
    private static SpotDectector fml;

    private static LeastSquaresProblem lsp;

    //replace IJ.log--> printlog
    public static void printlog(String msg) {
        if (false) { IJ.log(msg);}}
    
    //method to do rough localisation on frame-> gives an array of location of particles
    //referenced from Binary Spot detector used in TrackMate
    class Particle{
        public double[] positionArr;//[frameNum, positionX, positionY]
        public int index;//each particle object is labelled 1,2,3,4 etc
        public Particle(int idx){
            index = idx;
        }
    }
   
    public static class SPT {

        public SPT() {}
        private static doTwoFit x;
 
        //raw data

        private static Roi roi2;
        private static int ROI_width = 8; //adjust ROI_width and heigh to fit size of particle
        private static int ROI_height = 8;
        private static int shiftR = 3; 
        private static int frameno_ = 1;
        private static double SPTresultarray[][] = new double[totalFrame_][2];

        private static boolean iterate = false;
        int incrementX = 0;//if increment not zero, go to while loop and set newStart_2[1] =0
        int incrementY = 0;//if increment not zero, go to while loop and set newStart_2[2] =0

        private static double[] data_; //flattened 2d array
        //set initial parameters
        private static double[] newStart_ = {
            14000,//amplitude arbitrary//
            1,
            1,
            1,//depends on camera pixel size and magnification
            1,
            2000 // find average of all cnts in one frame
        };

        private static double[] optimalValues_1 = new double[6];
        private static int[] optim_param_;//for fitting
        private static int data_width_;//assuming its a sqaure region, width of square in pixel

        //blackbox that performs the fitting
        //private static double SPTresultarray[][] = new double[totalFrame_][2];
        //Main Method
        public static void StartTrack() throws IOException {
           newStart_[1] = ROI_width / 2;
           newStart_[2] = ROI_width / 2;
            //for display
            GenericDialog gd = new GenericDialog("Fitting 1 frame");
            gd.addMessage("Set Initial fit parameter");
            gd.addNumericField("Init Amplitude: ", newStart_[0], 0);
            gd.addNumericField("Init MeanX: ", newStart_[1], 0);
            gd.addNumericField("Init MeanY: ", newStart_[2], 0);
            gd.addNumericField("Init SigmaX: ", newStart_[3], 0);
            gd.addNumericField("Init SigmaY: ", newStart_[4], 0);
            gd.addNumericField("Init Offset: ", newStart_[5], 0);
            gd.addMessage("Select frame number (index start at 1)");
            gd.addNumericField("Frame no to be fitted (" + "1 to " + totalFrame_ + ")", frameno_, 0);

            gd.showDialog();

            if (gd.wasOKed()) {
                // Set initial fit parameters
                newStart_[0] = (double) gd.getNextNumber();
                newStart_[1] = (double) gd.getNextNumber();
                newStart_[2] = (double) gd.getNextNumber();
                newStart_[3] = (double) gd.getNextNumber();
                newStart_[4] = (double) gd.getNextNumber();
                newStart_[5] = (double) gd.getNextNumber();

                frameno_ = (int) gd.getNextNumber();
                printlog("frame no: " + frameno_);
                IJ.log("totalFrame: " + totalFrame_);
                IJ.log("--------START");
                boolean ffok = FitMeDialogue(frameno_, roi2);
                if (!ffok && iterate) {
                    boolean foundMatch = false;
                    int[][] combination = WayOfIteration(ROI_width, ROI_height);
                    for (int[] x : combination) {
                        newStart_[1] = x[0];
                        newStart_[2] = x[1];
                        if (FitMeDialogue(frameno_, roi2)) {
                            foundMatch = true;
                            break;
                        }
                    }
                }

                IJ.log("--------END");
               
               for (int framNum =frameno_+1; framNum <=totalFrame_; framNum++){
                    Integer[] newres;
                    for (int i=0; i<=4; i++){
                        newStart_[1]= ROI_width/2;
                        newStart_[2]=ROI_height/2;
                        if(i==0){
                            //uses meanX and meanY of previous frame to reposition the ROI
                            newres = getNewRoiPosition(optimalValues_1[1], optimalValues_1[2], 0, 0,ROI_width, ROI_height);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi2.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi2.getBounds().getX() + ", roi y: " + roi2.getBounds().getY() + ", roi w: " + roi2.getBounds().getWidth() + ", roi h: " + roi2.getBounds().getHeight());

                            if(FitMeDialogue(framNum,roi2))
                                break;
                            else if (iterate){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration(ROI_width, ROI_height);
                                for (int[] x:combination){
                                    newStart_[1]=x[0];
                                    newStart_[2]=x[1];
                                    if(FitMeDialogue(framNum,roi2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }
                        }else if(i==1){
                            newres = getNewRoiPosition(optimalValues_1[1], optimalValues_1[2], shiftR, 0, ROI_width, ROI_height);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi2.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi2.getBounds().getX() + ", roi y: " + roi2.getBounds().getY() + ", roi w: " + roi2.getBounds().getWidth() + ", roi h: " + roi2.getBounds().getHeight());

                            if(FitMeDialogue(framNum,roi2))
                                break;
                            else if (iterate){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration(ROI_width, ROI_height);
                                for (int[] x:combination){
                                   newStart_[1]=x[0];
                                   newStart_[2]=x[1];
                                   if(FitMeDialogue(framNum,roi2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }

                        }else if(i==2){
                            newres = getNewRoiPosition(optimalValues_1[1], optimalValues_1[2], -1*shiftR, 0,  ROI_width, ROI_height);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi2.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi2.getBounds().getX() + ", roi y: " + roi2.getBounds().getY() + ", roi w: " + roi2.getBounds().getWidth() + ", roi h: " + roi2.getBounds().getHeight());
                            if(FitMeDialogue(framNum,roi2))
                                break;
                            else if (iterate){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration( ROI_width, ROI_height);
                                for (int[] x:combination){
                                    newStart_[1]=x[0];newStart_[2]=x[1];
                                    if(FitMeDialogue(framNum,roi2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }

                        }else if(i==3){
                            newres = getNewRoiPosition(optimalValues_1[1], optimalValues_1[2], 0, shiftR,  ROI_width, ROI_height);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi2.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi2.getBounds().getX() + ", roi y: " + roi2.getBounds().getY() + ", roi w: " + roi2.getBounds().getWidth() + ", roi h: " + roi2.getBounds().getHeight());
                            if(FitMeDialogue(framNum,roi2))
                                break;
                            else if (iterate){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration(ROI_width, ROI_height);
                                for (int[] x:combination){
                                    newStart_[1]=x[0];newStart_[2]=x[1];
                                    if(FitMeDialogue(framNum,roi2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }
                        }else{
                            newres = getNewRoiPosition(optimalValues_1[1], optimalValues_1[2], 0, -1*shiftR, ROI_width, ROI_height);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi2.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi2.getBounds().getX() + ", roi y: " + roi2.getBounds().getY() + ", roi w: " + roi2.getBounds().getWidth() + ", roi h: " + roi2.getBounds().getHeight());
                            if(FitMeDialogue(framNum,roi2))
                                break;
                            else if (iterate){
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration( ROI_width, ROI_height);
                                for (int[] x:combination){
                                    newStart_[1]=x[0];newStart_[2]=x[1];
                                    if(FitMeDialogue(framNum,roi2))
                                        foundMatch=true;
                                        break;
                                }
                                if (foundMatch)
                                    break;
                            }
                        } else {
                            IJ.log("inside frameNum: " + framNum + ", i: " + i);
                            newres = getNewRoiPosition(optimalValues_1[1], optimalValues_1[2], 0, -1 * shiftR, ROI_width, ROI_height);
                            printlog("newres[0]: " + newres[0] + " newres[1]: " + newres[1]);
                            roi2.setLocation(newres[0], newres[1]);
                            printlog("roi x: " + roi2.getBounds().getX() + ", roi y: " + roi2.getBounds().getY() + ", roi w: " + roi2.getBounds().getWidth() + ", roi h: " + roi2.getBounds().getHeight());
                            if (FitMeDialogue(framNum, roi2)) {
                                break;
                            } else if (iterate) {
                                boolean foundMatch = false;
                                int[][] combination = WayOfIteration(ROI_width, ROI_height);
                                for (int[] x : combination) {
                                    newStart_[1] = x[0];
                                    newStart_[2] = x[1];
                                    IJ.log("inside iterate frameNUM: " + framNum + ", newStart[1]: " + newStart_[1] + ", newStart[2]: " + newStart_[2]);
                                    if (FitMeDialogue(framNum, roi2)) {
                                        foundMatch = true;
                                        break;
                                    }
                                }
                                if (foundMatch) {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            int t = 0;
            for (double[] q : SPTresultarray) {
                t++;
                System.out.println("frame:" + t + " " + Arrays.toString(q));
            }
            writeFile();
        }

      
        //Methods
        private static int[][] WayOfIteration(int boxWidth, int boxHeight) {
            int[][] combination = new int[4][2];
            combination[0][0] = (int) Math.floor(boxWidth * 0.25);
            combination[0][1] = (int) Math.floor(boxHeight * 0.25);
            combination[1][0] = (int) Math.floor(boxWidth * 0.75);
            combination[1][1] = (int) Math.floor(boxHeight * 0.25);
            combination[2][0] = (int) Math.floor(boxWidth * 0.75);
            combination[2][1] = (int) Math.floor(boxHeight * 0.75);
            combination[3][0] = (int) Math.floor(boxWidth * 0.25);
            combination[3][1] = (int) Math.floor(boxHeight * 0.75);
            return combination;
        }

        private static boolean FitMeDialogue(int framNum, Roi roi2) {
            data_ = get1Darray(imp_, framNum, roi2);
            data_width_ = ROI_width;
            printlog("Init Amplitude: " + newStart_[0] + ", Init MeanX: " + newStart_[1] + ", Init MeanY: " + newStart_[2] + ", Init SigmaX: " + newStart_[3] + ", Init SigmaY: " + newStart_[4] + ", Init Offset: " + newStart_[5]);
            printlog("fitting frame no: " + framNum);
            x = new doTwoFit(data_, newStart_, data_width_, new int[]{1000, 100});
           
            try {//do LevenbergMarquardt optimization and get optimized parameters
                LeastSquaresOptimizer.Optimum opt = x.fit2dGauss();

                final double[] optimalValues = opt.getPoint().toArray();//opt after fitting is done 
                printlog("optimalalues[1]: " + optimalValues[1]);
                printlog("optimalalues[2]: " + optimalValues[2]);
                printlog("roi x: " + roi2.getBounds().getX() + ", roi y: " + roi2.getBounds().getY() + ", roi w: " + roi2.getBounds().getWidth() + ", roi h: " + roi2.getBounds().getHeight());

                optimalValues[1] = optimalValues[1] + roi2.getBounds().getX(); //changesmade
                optimalValues[2] = optimalValues[2] + roi2.getBounds().getY();

                printlog("amplitude" + optimalValues[0] + "meanx: " + optimalValues[1] + " meany: " + optimalValues[2] + " sigmax: " + optimalValues[3] + " sigma y: " + optimalValues[4] + " offset: " + optimalValues[5]);

                //fill textfield UI
                tfAmplitude.setText(Double.toString(optimalValues[0]));
                tfMeanX.setText(Double.toString(optimalValues[1]));
                tfMeanY.setText(Double.toString(optimalValues[2]));
                tfSigmaX.setText(Double.toString(optimalValues[3]));
                tfSigmaY.setText(Double.toString(optimalValues[4]));
                tfOffset.setText(Double.toString(optimalValues[5]));
                optimalValues_1 = optimalValues;
                SPT.SPTresultarray[framNum - 1][0] = optimalValues[1];//mean X
                SPT.SPTresultarray[framNum - 1][1] = optimalValues[2];//mean Y
                SPT.SPTresultarray[framNum - 1][2] = optimalValues[3];//sd X
                SPT.SPTresultarray[framNum - 1][3] = optimalValues[4];//sd Y
                SPT.SPTresultarray[framNum - 1][4] = optimalValues[0];//Aeff
                SPT.SPTresultarray[framNum - 1][5] = optimalValues[5];//Offset
                printlog("optimalalues_1[1]: " + optimalValues_1[1]);
                printlog("optimalalues_1[2]: " + optimalValues_1[2]);

                //output data
                System.out.println("v0: " + optimalValues[0]);
                System.out.println("v1: " + optimalValues[1]);
                System.out.println("v2: " + optimalValues[2]);
                System.out.println("v3: " + optimalValues[3]);
                System.out.println("v4: " + optimalValues[4]);
                System.out.println("v5: " + optimalValues[5]);
                System.out.println("Iteration number: " + opt.getIterations());
                System.out.println("Evaluation number: " + opt.getEvaluations());
            } catch (Exception e) {
                tfAmplitude.setText(Double.toString(Double.NaN));
                tfMeanX.setText(Double.toString(Double.NaN));
                tfMeanY.setText(Double.toString(Double.NaN));
                tfSigmaX.setText(Double.toString(Double.NaN));
                tfSigmaY.setText(Double.toString(Double.NaN));
                tfOffset.setText(Double.toString(Double.NaN));
                System.out.println(e.toString());
                printlog("-----First frame fit NOK");
                return false;
            }
            return true;
        }

        private static double[] get1Darray(ImagePlus imp, int frameindex, Roi roi) { // index start at 1
            printlog("inside get1Darray framindex: " + frameindex + ", roix: " + roi.getBounds().getX() + ", roiy: " + roi.getBounds().getY() + ", roiwidth: " + roi.getBounds().getWidth() + ", roiheight: " + roi.getBounds().getHeight());

            Rectangle rect = roi.getBounds();
            int l_disp = (int) rect.getX(); 
            int t_disp = (int) rect.getY();
            int w_disp = (int) rect.getWidth();
            int h_disp = (int) rect.getHeight();
            int f = frameindex;

            double[] res = new double[w_disp * h_disp];

            for (int y = 0; y < h_disp; y++) {
                for (int x = 0; x < w_disp; x++) {
                    res[x + (y * w_disp)] = imp.getStack().getProcessor(f).get(x + l_disp, y + t_disp);
                }
            }
            return res;//1D ARRAY
        }

         public static Integer[] getNewRoiPosition(double fittedX, double fittedY, int shiftx, int shifty, int boxWidth, int boxHeight) {
            // index start 0 fittedX,fittedY
            //index start 1 initW, initH
            IJ.log("inside getNewRoiPosition fittexX: " + fittedX + ", fittedY: " + fittedY + ", shiftx: " + shiftx + ", shifty: " + shifty);
            Integer[] result = new Integer[2];
            int new_left = (int) (Math.round(fittedX - (boxWidth / 2) + shiftx));
            int new_top = (int) (Math.round(fittedY - (boxHeight / 2) + shifty));

            if (new_left < 0) {
                new_left = 0;
            }

            if (new_top < 0) {
                new_top = 0;
            }

            if ((new_left + boxWidth  - 1) > (imp_.getWidth() - 1)) {
                new_left = imp_.getWidth() - boxWidth ;
            }

            if ((new_top + boxHeight - 1) > (imp_.getHeight() - 1)) {
                new_top = imp_.getHeight() - boxHeight;
            }
            result[0] = new_left;
            result[1] = new_top;

            IJ.log("getNewRoiPosition l: " + result[0] + ", t: " + result[1]);
            return result;
        }
        public static void initializeroi(int initposx, int initposy) {
            IJ.log("initposx: " + initposx + ", initposy: " + initposy);
            int tempX = 0;
            int tempY = 0;
            if ((initposx - ROI_width / 2) >= 0) {
                tempX = initposx - ROI_width / 2;
            }
            if ((initposy - ROI_height / 2) >= 0) {
                tempY = initposy - ROI_height / 2;
            }
            roi2 = new Roi(tempX, tempY, ROI_width, ROI_height);
            imp_.setRoi(roi2);
        }

    }
    //End of SPT class

    public static class doTwoFit {
        public doTwoFit(double[] data, double[] newStart, int data_width, int[] optim_param) { //takes 4 arg and stores in global variables
            SPT.data_ = data;
            SPT.newStart_ = newStart;
            SPT.data_width_ = data_width;
            SPT.optim_param_ = optim_param;
            buildlsb();//Performs fitting
        }
        
        /**
         * build LeastSquareProblem by using constructor data
         */
        private void buildlsb() {
            //construct two-dimensional Gaussian function
            TwoDGaussFunction tdgf = new TwoDGaussFunction(SPT.data_width_, SPT.data_.length);

            //prepare construction of LeastSquresProblem by builder
            LeastSquaresBuilder lsb = new LeastSquaresBuilder();

            //set model function and its jacobian
            lsb.model(tdgf.retMVF(), tdgf.retMMF());
            //set target data
            lsb.target(SPT.data_);
            //set initial parameters
            lsb.start(SPT.newStart_);
            //set upper limit of evaluation time
            lsb.maxEvaluations(SPT.optim_param_[0]);
            //set upper limit of iteration time
            lsb.maxIterations(SPT.optim_param_[1]);

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

    public static JFrame TwoDGaussFrame = new JFrame("TRY 2D Gauss SPT");
    private JButton btnLoad;
    private JButton btnFitme1;//first algorithm
    private JButton btnFitme2;//second algorithm
    private JButton btnTest;
    private static JTextField tfAmplitude;
    private static JTextField tfMeanX;
    private static JTextField tfMeanY;
    private static JTextField tfSigmaX;
    private static JTextField tfSigmaY;
    private static JTextField tfOffset;

    public void createPanel() {
        TwoDGaussFrame.setFocusable(true);
        TwoDGaussFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        TwoDGaussFrame.setLayout(new GridLayout(10, 2)); //this is y,x?
        TwoDGaussFrame.setLocation(new Point(390, 125));
        TwoDGaussFrame.setSize(new Dimension(350, 230));
        TwoDGaussFrame.setResizable(false);

        //Init
        btnLoad = new JButton("Load tiff");
        btnFitme1 = new JButton("FitMe1");
        btnFitme2 = new JButton("FitMe2");
        btnTest = new JButton("Test");
        tfMeanX = new JTextField(Double.toString(0), 4);
        tfMeanX.setEditable(false);
        tfMeanY = new JTextField(Double.toString(0), 4);
        tfMeanY.setEditable(false);
        tfAmplitude = new JTextField(Double.toString(0), 4);
        tfAmplitude.setEditable(false);
        tfSigmaX = new JTextField(Double.toString(0), 4);
        tfSigmaX.setEditable(false);
        tfSigmaY = new JTextField(Double.toString(0), 4);
        tfSigmaY.setEditable(false);
        tfOffset = new JTextField(Double.toString(0), 4);
        tfOffset.setEditable(false);

        //row1
        TwoDGaussFrame.add(new JLabel("Step 1 load imagestack"));
        TwoDGaussFrame.add(btnLoad);
        TwoDGaussFrame.add(new JLabel("..."));
        TwoDGaussFrame.add(btnFitme1);
        TwoDGaussFrame.add(new JLabel("..."));
        TwoDGaussFrame.add(btnFitme2);
        TwoDGaussFrame.add(new JLabel("Amplitude"));
        TwoDGaussFrame.add(tfAmplitude);
        TwoDGaussFrame.add(new JLabel("Mean x"));
        TwoDGaussFrame.add(tfMeanX);
        TwoDGaussFrame.add(new JLabel("Mean y"));
        TwoDGaussFrame.add(tfMeanY);
        TwoDGaussFrame.add(new JLabel("Sigma x"));
        TwoDGaussFrame.add(tfSigmaX);
        TwoDGaussFrame.add(new JLabel("Sigma y"));
        TwoDGaussFrame.add(tfSigmaY);
        TwoDGaussFrame.add(new JLabel("Offset"));
        TwoDGaussFrame.add(tfOffset);
        TwoDGaussFrame.add(btnTest);
        TwoDGaussFrame.add(new JLabel("..."));
        

        // add listener
        btnLoad.addActionListener(btnLoadPressed);
        btnFitme1.addActionListener(btnFitmePressed1);
        btnFitme2.addActionListener(btnFitmePressed2);
        btnTest.addActionListener(btnTestPressed);

        TwoDGaussFrame.setVisible(true);

    }

    ActionListener btnLoadPressed = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent arg0) {
            imp_ = IJ.openImage();

            totalFrame_= imp_.getStackSize();
            if (imp_ != null) {
                imp_.show();
                imp_can = imp_.getCanvas();
                imp_can.setFocusable(true);
                imp_can.addMouseListener(impcanMouseListener);
            }
        }
    };

    ActionListener btnFitmePressed1 = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent arg0) {
            try {
                SPT.StartTrack();
            } catch (IOException ex) {
                Logger.getLogger(TwoDGaussianPlugin.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    };
    ActionListener btnFitmePressed2 = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent arg0) {
            fml = new SpotDectector(imp_);
            fml.loopParticles();   
        }
    };
    
    ActionListener btnTestPressed = new ActionListener(){
        @Override
        public void actionPerformed(ActionEvent e) {
            try {
                SpotDectector.printResult();
            } catch (IOException ex) {
                Logger.getLogger(TwoDGaussianPlugin.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    };

    MouseListener impcanMouseListener = new MouseListener() { //declare imp canvas 
        @Override
        public void mouseClicked(MouseEvent e) {
            int px = e.getX();
            int py = e.getY();

            //recalculate this into coordinates in the image
            int initposx = (int) Math.floor(imp_can.offScreenX(px));
            int initposy = (int) Math.floor(imp_can.offScreenY(py));

            SPT.initializeroi(initposx, initposy);
        }

        @Override
        public void mousePressed(MouseEvent e) {
        }

        @Override
        public void mouseReleased(MouseEvent e) {
        }

        @Override
        public void mouseEntered(MouseEvent e) {
            //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public void mouseExited(MouseEvent e) {
            //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
    };
    
    public static void writeFile() throws IOException{
        StringBuilder builder = new StringBuilder();
        for(int i =0;i<SPT.SPTresultarray.length;i++){
            builder.append(Arrays.toString(SPT.SPTresultarray[i])+"");
            builder.append("\n");
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter("sampledata.txt"));
        writer.write(builder.toString());
        writer.close();
       
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        TwoDGaussianPlugin twodobj = new TwoDGaussianPlugin();
        twodobj.createPanel();
        //System.out.println(SPT.SPTresultarray);

    }
}
