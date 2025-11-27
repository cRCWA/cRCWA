package visualization;

import org.jtransforms.fft.DoubleFFT_2D;

import javax.swing.*;
import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.text.*;
import java.util.*;

import storage.*;
import timer.*;

/**
    Panel employed to draw extremely fast data files coming from different
    sources.

    Copyright 2012-2017 by Davide Bucci
    @author Davide Bucci
*/
public class GraphingPanel extends JPanel
{
    // Store the data to be drawn
    MatrixStore sd;
    boolean hasChanged;
    Color[] palette;
    private final static int MAXPALETTECOL = 256;
    // Contains the rendered image
    private BufferedImage bufferedImage;
    private BufferedImage colorScale;
    private String fileFormat;

    private boolean isScatter;
    private Color foreColor;

    private boolean automaticPointSize;
    private int scatterColumnX;
    private int scatterColumnY;
    private int scatterColumnZ;
    private int scatterColumnR;
    private int scatterColumnI;
    private boolean measureTime;

    private int wantedImageSizeX;
    private int wantedImageSizeY;

    private double min;
    private double max;
    private double xmin;
    private double ymin;
    private double xmax;
    private double ymax;
    private double zmin;
    private double zmax;
    private double zquote;
    private double ztol;
    private int xSize;
    private int ySize;

    private boolean calculateMin; // Evaluate the min value in raster drawing
    private boolean calculateMax; // Evaluate the max value in raster drawing

    // Those values are used for the raster drawing if the calculateMin and
    // calculateMax flags are false. If they are true, they are updated to
    // the value to be used, evaluated each time a new raster image is
    // drawn.

    private double crangeMin;     // Min value in raster drawing
    private double crangeMax;     // Max value in raster drawing

    // Visualisation settings
    public final static int VIS_REAL = 0;
    public final static int VIS_IMAG = 1;
    public final static int VIS_MAGNITUDE = 2;
    public final static int VIS_LOGMAG = 4;
    public final static int VIS_PHASE = 3;

    // Preprocessing settings
    public final static int PRE_NO = 0;
    public final static int PRE_FFT = 1;
    public final static int PRE_IFFT = 2;

    private int visType;    // Type of visualisation
    private int preType;    // Type of preprocessing

    /** Constructor for the panel.
    */
    public GraphingPanel()
    {
        foreColor = Color.green.darker().darker();
        sd = null;
        automaticPointSize=true;
        xSize=5;
        ySize=5;

        calculateMin=true;
        calculateMax=true;
        wantedImageSizeX = 500;
        wantedImageSizeY = 400;


        visType = VIS_REAL;
        fileFormat="";
        isScatter=false;
        palette=createGrayscalePalette();

        zquote=0;
        ztol=1e-4;
        hasChanged = true;
        bufferedImage =null;
        scatterColumnX = 0;
        scatterColumnY = 0;
        scatterColumnZ = 0;
        scatterColumnR = 0;
        scatterColumnI = 0;
    }

    /** Set the current palette.
        @param p the palette to be used.
    */
    public void setPalette(Color[] p)
    {
        palette=p;
    }

    /** Set a grey scale palette.
        @return a grayscale palette.
    */
    public Color[] createGrayscalePalette()
    {
        Color[] p = new Color[MAXPALETTECOL];
        // Creat a greyscale palette
        double c = 256.0 / MAXPALETTECOL;
        for(int i=0; i<MAXPALETTECOL; ++i){
            p[i] = new Color((int)(i*c), (int)(i*c), (int)(i*c));
        }
        return p;
    }

     /**
        Decode an HTML color string like '#F567BA;' into a {@link Color}
        @param colorString The string to decode
        @return The decoded color
        @throws IllegalArgumentException if the color sequence is not valid
        from http://www.java2s.com/Tutorial/Java/0120__Development/\
        DecodeanHTMLcolorstringlikeF567BAintoaColor.htm
     */
    public static Color decodeHtmlColorString(String colorString)
    {
        Color color;

        if (colorString.startsWith("#")) {
            colorString = colorString.substring(1);
        }
        if (colorString.endsWith(";")) {
            colorString = colorString.substring(0, colorString.length() - 1);
        }

        int red, green, blue;
        switch (colorString.length()) {
            case 6:
                red = Integer.parseInt(colorString.substring(0, 2), 16);
                green = Integer.parseInt(colorString.substring(2, 4), 16);
                blue = Integer.parseInt(colorString.substring(4, 6), 16);
                color = new Color(red, green, blue);
                break;
            case 3:
                red = Integer.parseInt(colorString.substring(0, 1), 16);
                green = Integer.parseInt(colorString.substring(1, 2), 16);
                blue = Integer.parseInt(colorString.substring(2, 3), 16);
                color = new Color(red, green, blue);
                break;
            case 1:
                red = green = blue = Integer.parseInt(colorString.substring(0,
                    1), 16);
                color = new Color(red, green, blue);
                break;
            default:
                throw new IllegalArgumentException("Invalid color: " +
                    colorString);
        }
        return color;
    }

    /** Creates the Matlab "Jet" palette.
        @return the "Jet" palette.
    */
    public Color[] createJetPalette()
    {
        int n_max=9;
        Color[] pal= new Color[n_max];
        float[] idx= new float[n_max];
        pal[0]= decodeHtmlColorString("#000090");
        pal[1]= decodeHtmlColorString("#000fff");
        pal[2]= decodeHtmlColorString("#0090ff");
        pal[3]= decodeHtmlColorString("#0fffee");
        pal[4]= decodeHtmlColorString("#90ff70");
        pal[5]= decodeHtmlColorString("#ffee00");
        pal[6]= decodeHtmlColorString("#ff7000");
        pal[7]= decodeHtmlColorString("#ee0000");
        pal[8]= decodeHtmlColorString("#7f0000");
        idx[0]=0.0f * MAXPALETTECOL;
        idx[1]=1.0f * MAXPALETTECOL;
        idx[2]=2.0f * MAXPALETTECOL;
        idx[3]=3.0f * MAXPALETTECOL;
        idx[4]=4.0f * MAXPALETTECOL;
        idx[5]=5.0f * MAXPALETTECOL;
        idx[6]=6.0f * MAXPALETTECOL;
        idx[7]=7.0f * MAXPALETTECOL;
        idx[8]=8.0f * MAXPALETTECOL;
        return createGenericPalette(pal, idx, n_max, false);
    }

    /** "Temperature" palette.
        @return the "temp." palette.
    */
    public Color[] createTempPalette()
    {
        int n_max=5;
        Color[] pal= new Color[n_max];
        float[] idx= new float[n_max];
        int i=0;
        pal[i++]= decodeHtmlColorString("#0909FF");
        pal[i++]= decodeHtmlColorString("#1050FF");
        pal[i++]= decodeHtmlColorString("#FF0000");
        pal[i++]= decodeHtmlColorString("#FFFF00");
        pal[i++]= decodeHtmlColorString("#FFFFFF");

        i=0;
        idx[i++]=0.0f * MAXPALETTECOL;
        idx[i++]=1.2f * MAXPALETTECOL;
        idx[i++]=2.6f * MAXPALETTECOL;
        idx[i++]=3.9f * MAXPALETTECOL;
        idx[i++]=5.0f * MAXPALETTECOL;

        return createGenericPalette(pal, idx, n_max, true);
    }

    /** Creates the Blue/White/Red palette. Useful to separate positive
        and negative values, centered in zero.
        @return the palette.
    */
    public Color[] createBlueRedPalette()
    {
        // Creat a Blue/White/Red palette
        double c = 256.0 / MAXPALETTECOL;

        int n_max=3;
        Color[] pal= new Color[n_max];
        float[] idx= new float[n_max];
        pal[0]= decodeHtmlColorString("#ff0000");
        pal[1]= decodeHtmlColorString("#ffffff");
        pal[2]= decodeHtmlColorString("#0000ff");
        idx[0]=0.0f * MAXPALETTECOL;
        idx[1]=1.5f * MAXPALETTECOL;
        idx[2]=3.0f * MAXPALETTECOL;
        return createGenericPalette(pal, idx, n_max, false);
    }

    /** Create a generic palette.
    */
    private Color[] createGenericPalette(Color [] pal, float[] idx, int n_max,
        boolean useGamma)
    {
        float r=0.5f;
        Color color1=Color.GRAY;
        Color color2=Color.GRAY;
        Color [] p = new Color[MAXPALETTECOL];

        for(int i=0; i<MAXPALETTECOL; ++i){
            for (int j=1; j<n_max;++j) {
                if(i*n_max<idx[j]) {
                    color1=pal[j-1];
                    color2=pal[j];
                    r=(idx[j]-i*n_max)/(idx[j]-idx[j-1]);
                    break;
                }
            }
            if(useGamma) {
                p[i] = blendColorsGamma(color1, color2, r);
            } else {
                p[i] = blendColors(color1, color2, r);
            }
        }
        return p;
    }

    /**
        Blend two colors. From
     http://www.java2s.com/Code/Java/2D-Graphics-GUI/Commoncolorutilities.htm

        @param color1  First color to blend.
        @param color2  Second color to blend.
        @param r   Blend ratio. 0.5 will give even blend, 1.0 will return
                   color1, 0.0 will return color2 and so on.
        @return        Blended color.
    */
    public static Color blendColors(Color color1, Color color2, float r)
    {
        float ir = (float) 1.0 - r;

        float rgb1[] = new float[3];
        float rgb2[] = new float[3];

        color1.getColorComponents (rgb1);
        color2.getColorComponents (rgb2);

        Color color = new Color (rgb1[0] * r + rgb2[0] * ir,
                             rgb1[1] * r + rgb2[1] * ir,
                             rgb1[2] * r + rgb2[2] * ir);

        return color;
    }

    /** Takes into account a gamma of 0.45/2.2 when blending colors. This
        allows to obtain much more natural color rendering in some cases.
    */
    public static Color blendColorsGamma
        (Color color1, Color color2, float r)
    {
        float ir = (float) 1.0 - r;

        float rgb1[] = new float[3];
        float rgb2[] = new float[3];

        color1.getColorComponents (rgb1);
        color2.getColorComponents (rgb2);

        float dir=0.45f;
        float inv=1/dir;

        double rr=Math.pow(rgb1[0],inv) * r + Math.pow(rgb2[0],inv) * ir;
        double gg=Math.pow(rgb1[1],inv) * r + Math.pow(rgb2[1],inv) * ir;
        double bb=Math.pow(rgb1[2],inv) * r + Math.pow(rgb2[2],inv) * ir;

        rr=Math.pow(rr,dir);
        gg=Math.pow(gg,dir);
        bb=Math.pow(bb,dir);

        if(rr>1.0) rr=1.0;
        if(gg>1.0) gg=1.0;
        if(bb>1.0) bb=1.0;

        Color color = new Color ((float)rr, (float)gg, (float)bb);
        return color;
    }

    public void setAutomaticPointSize(boolean b)
    {
        automaticPointSize=b;
        forceRedraw();
    }

    public void setSizePoint(int x, int y)
    {
        xSize=x;
        ySize=y;
        forceRedraw();
    }

    public void chooseVisualization(int v)
    {
        visType=v;
        hasChanged=true;
        repaint();
    }

    /** Calculate the visualisation value from the real and imaginary part of
        the datum to be represented.
        @param r the real part
        @param i the imaginary part
        @return the result.
    */
    private double calculateVisualization(double r, double i)
    {
        double v;

        switch(visType){
            case VIS_REAL:      // Real part
                v=r;
                break;
            case VIS_IMAG:      // Imaginary part
                v=i;
                break;
            case VIS_MAGNITUDE: // Magnitude
                v=Math.sqrt(r*r+i*i);
                break;
            case VIS_PHASE:     // Phase
                v=Math.atan(i/r);
                break;
            case VIS_LOGMAG:    // Log(magnitude)
                 v=Math.log(Math.sqrt(r*r+i*i));
                break;
            default:
                v=0;
        }
        return v;
    }

    public void measureTimings(boolean b)
    {
        measureTime =b;
    }

    public void setScatterType(int sx, int sy, int sz, int sr, int si)
    {
        isScatter = true;

        scatterColumnX = sx;
        scatterColumnY = sy;
        scatterColumnZ = sz;
        scatterColumnR = sr;
        scatterColumnI = si;
        if(sd==null) return;

        hasChanged=true;
    }
    public void setRasterType()
    {
        isScatter = false;
    }

    public void setFileFormat(String f)
    {
        fileFormat =f;
    }

    /** Change the data to be drawn
    */
    public void setData(MatrixStore m)
    {
        sd = m;
        hasChanged = true;
    }

    /** Change the z quota to be drawn
    */
    public void setZQuote(double z)
    {
        zquote=z;
    }

    /** Change the z tolerance to be used when drawing.
    */
    public void setZToll(double t)
    {
        ztol=t;
    }

    /** Draw a scattering graph.
    */
    private void drawScatter(Graphics2D gg,
            MatrixStore d, int xsize, int ysize)
    {
        int x;
        int y;
        int xinc, yinc;
        int npoints = d.getNrow();
        if(scatterColumnX<0 || scatterColumnY<0)
            return;
        if(d.getNcol()<=scatterColumnY||d.getNrow()<1)
            return;


        xmin=d.getElement(0,scatterColumnX,0);
        xmax=xmin;
        ymin=d.getElement(0,scatterColumnY,0);
        ymax=ymin;
        double vx;
        double vy;
        double v;

        int index;

        if(scatterColumnR>=0) {
            if(scatterColumnR>=d.getNcol())
                return;
            min=d.getElement(0,scatterColumnR,0);
            max=min;
        }
        // At first, we determine the min and max values
        // of the x and y coordinates
        for(int i=0; i<npoints;++i) {
            vx = d.getElement(i,scatterColumnX,0);
            vy = d.getElement(i,scatterColumnY,0);
            if(scatterColumnR>=0){
                v = d.getElement(i,scatterColumnR,0);
                if(scatterColumnI>0) {
                    v=calculateVisualization(v,
                        d.getElement(i,scatterColumnI,0));
                }
                if(v<min)
                    min=v;
                if(v>max)
                    max=v;
            }
            if(vx<xmin)
                xmin=vx;
            if(vx>xmax)
                xmax=vx;
            if(vy<ymin)
                ymin=vy;
            if(vy>ymax)
                ymax=vy;
        }

        // If the flags calculateXXX are true, update the limits to the
        // min/max values calculated above. If they are false, overwrite
        // the limits employed in the drawing.
        if (calculateMin)
            crangeMin=min;
        else
            min=crangeMin;

        if(calculateMax)
            crangeMax=max;
        else
            max=crangeMax;

        // Then, we calculate the visualization parameters
        int xs = d.getScanline(scatterColumnY, 0);

        int zs=0;
        if(scatterColumnZ>=0)
            zs= d.getScanline(scatterColumnZ,0);

        int ppy;
        if(zs>0)
            ppy=zs;
        else
            ppy=npoints;
        int ys;
        boolean lineconnect=false;
        if(xs>1) {
            if(automaticPointSize) {
                ys=ppy/xs;
                xinc=xsize/xs+1;
                if(ys>0)
                    yinc=ysize/ys+1;
                else
                    yinc=1;
            } else {
                xinc=xSize;
                yinc=ySize;
            }
            gg.setColor(foreColor);
        } else {
            xinc = 1;
            yinc = 1;
            lineconnect=true;
            gg.setColor(Color.white);
        }

        // We draw a background, to visually differentiate where the data is
        // obtained from the data file, and where it is not.
        gg.fillRect(0,0, xsize, ysize);

        double vz;
        // Then, we plot the scatter diagram
        int ox=0, oy=0;
        int writtenPoints=0;
        for(int i=0; i<npoints;++i) {
            vx = d.getElement(i,scatterColumnX,0);
            vy = d.getElement(i,scatterColumnY,0);
            if(scatterColumnZ>=0) {
                vz = d.getElement(i,scatterColumnZ,0);
                if (Math.abs(vz-zquote)>ztol*Math.abs(vz))
                    continue;
            }
            ++writtenPoints;
            x = (int)((vx-xmin)/(xmax-xmin)*(double)xsize);
            y = (int)((vy-ymin)/(ymax-ymin)*(double)ysize);
            if(scatterColumnR>0){
                v = d.getElement(i, scatterColumnR,0);
                if(scatterColumnI>0) {
                    v=calculateVisualization(v,
                        d.getElement(i,scatterColumnI,0));
                }

                index = (int)((v - min)*((double)MAXPALETTECOL-1)/(max-min));
            } else
                index = 0;

            // Check if the index is in the limits (this should happen when)
            // the min and max have been automatically calculated, but
            // it might not be the case when they have been manually entered.

            if(index<0) {
                index=0;
            } else if(index>MAXPALETTECOL-1) {
                index=MAXPALETTECOL-1;
            }

            gg.setColor(palette[index]);
            if(lineconnect) {
                if(i>0) {
                    gg.drawLine(ox, ysize-oy, x,ysize-y);

                }
                ox=x; oy=y;
            } else {
                gg.fillRect(x,ysize-y,xinc,yinc);
            }
        }
    }

    /** Set the z-slice at a given percent of the whole interval of data
        present in memory.
    */
    public void setZpercent(int p)
    {
        if(sd==null||scatterColumnZ<0)
            return;
        double v;
        // Measure max and min of data towards z.
        zmax=zmin=sd.getElement(0,scatterColumnZ,0);
        for(int i=0; i<sd.getNrow();++i) {
            v = sd.getElement(i,scatterColumnZ,0);
            if (v<zmin)
                zmin=v;
            if(v>zmax)
                zmax=v;
        }

        zquote = (double)p/100.0*(zmax-zmin)+zmin;

        double d=1e308;
        double rd=zquote;
        for(int i=0; i<sd.getNrow();++i) {
            v = sd.getElement(i,scatterColumnZ,0);
            if(Math.abs(v-zquote)<d) {
                d=Math.abs(v-zquote);
                rd=v;
            }
        }
        zquote = rd;
        forceRedraw();
    }

    /** Force an in-depth redraw
    */
    public void forceRedraw()
    {

        hasChanged=true;
        repaint();
    }

    /** Return true if the color scale of the picture should be adapted to
        the minimum value contained in the data.
    */
    public boolean getCalculateMin()
    {
        return calculateMin;
    }

    /** Determine if the color scale of the picture should be adapted to
        the minimum value contained in the data.
    */
    public void setCalculateMin(boolean c)
    {
        calculateMin=c;
    }

    /** Return the minimum value of the color.
    */
    public double getCrangeMin()
    {
        return crangeMin;
    }

    /** Sets the maximum value of the color.
    */
    public void setCrangeMin(double c)
    {
        crangeMin=c;
    }

    /** Return true if the color scale of the picture should be adapted to
        the maximum value contained in the data.
    */
    public boolean getCalculateMax()
    {
        return calculateMax;
    }

    /** Determine if the color scale of the picture should be adapted to
        the maximum value contained in the data.
    */
    public void setCalculateMax(boolean c)
    {
        calculateMax=c;
    }

    /** Return the maximum value of the color.
    */
    public double getCrangeMax()
    {
        return crangeMax;
    }

    /** Sets the maximum value of the color.
    */
    public void setCrangeMax(double c)
    {
        crangeMax=c;
    }

    /** Update the min and max values for the visualisation of data contained
        in matrix d. The way the value is calculated depends on the current
        visualisation option. For example, if a "module" visualisation is
        currently active, the module of the complex data contained in the
        matrix will be used for the calculations. If instead the "real" is
        employed, only the real part will be considered and so on.

        @param d the data to be taken for calculating the min/max value.
    */
    private void calculateMinMax(MatrixStore d)
    {
        double v=calculateVisualization(
                d.getElement(0,0,0),
                d.getElement(0,0,1));
        min=max=v;
        // Determine the minimum and maximum values to be used for the
        // color scale.

        if(calculateMin || calculateMax) {
            for(int i=0; i<d.getNrow(); ++i) {
                for(int j=0; j<d.getNcol();++j){
                    v=calculateVisualization(
                        d.getElement(i,j,0),
                        d.getElement(i,j,1));
                    if(v<min) {
                        min=v;
                    }
                    if(v>max) {
                        max=v;
                    }
                }
            }
        }

        // If the flags calculateXXX are true, update the limits to the
        // min/max values calculated above. If they are false, overwrite
        // the limits employed in the drawing.
        if (calculateMin)
            crangeMin=min;
        else
            min=crangeMin;

        if(calculateMax)
            crangeMax=max;
        else
            max=crangeMax;
    }

    /** Draw a raster picture.
    */
    private void drawRaster(Graphics2D gg, MatrixStore d,
            int xsize, int ysize)
    {
        if(d==null)
            return;

        int x;
        int y;
        int xinc=xsize/d.getNcol()+1;
        int yinc=ysize/d.getNrow()+1;
        double v;
        int index;

        if(d.getNPlanes()<2)
            return;

        for(int i=0; i<d.getNrow(); ++i) {
            for(int j=0; j<d.getNcol();++j){
                x = j*xsize/d.getNcol();
                y = i*ysize/d.getNrow();

                v=calculateVisualization(
                    d.getElement(i,j,0),
                    d.getElement(i,j,1));

                // Calculate the index in the palette table.
                index = (int)((v - min)*((double)MAXPALETTECOL-1)/(max-min));

                // Handle a possible saturation of the values.
                if(index<0)
                    index=0;
                if(index>MAXPALETTECOL-1)
                    index=MAXPALETTECOL-1;

                // Draw the point!
                gg.setColor(palette[index]);
                gg.fillRect(x,ysize-y-1,xinc,yinc);
            }
        }
    }

    /** Create an image which can be used for the redraw operation without
        recalculating everything about the visualisation.
        @param d the data to be represented.
        @param xsize the horizontal size in pixels.
        @param ysize the vertical size in pixels.
    */
    private BufferedImage createImage(MatrixStore d, int xsize, int ysize)
    {
        MyTimer mt = new MyTimer();

        BufferedImage bImage = new BufferedImage(xsize, ysize,
                            BufferedImage.TYPE_INT_BGR);

        Graphics2D gg = bImage.createGraphics();

        // Activate antialiasing for any drawing operation.
        gg.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                         RenderingHints.VALUE_ANTIALIAS_ON);

        if(d==null) {
            gg.setColor(Color.white);
            gg.fillRect(0,0, xsize, ysize);
            gg.setColor(foreColor);
            gg.drawLine(0,0,xsize+1,ysize+1);
            gg.drawLine(0,ysize,xsize,0);
            gg.drawRect(0,0, xsize-1, ysize-1);
            return bImage;
        }

        gg.setColor(Color.white);
        gg.fillRect(0,ysize, xsize, 50);
        gg.setColor(foreColor);
        if(isScatter) {
            drawScatter(gg,d, xsize, ysize);
        } else {
            if(preType!=PRE_NO) {
                MatrixStore m;
                if(preType==PRE_FFT) {
                    m=preprocessing(d, true);
                } else {
                    m=preprocessing(d, false);
                }
                calculateMinMax(m);
                drawRaster(gg,m, xsize, ysize);
            } else {
                calculateMinMax(d);
                drawRaster(gg,d, xsize, ysize);
            }
        }

        hasChanged=false;
        if(measureTime)
            System.out.println("Image completed in: "+mt.getElapsed()+" ms.");
        return bImage;
    }

    /** Draw the y axis with the corresponding scale.
        @param g the graphic context on which to draw.
        @param xpos the x coordinate of the leftmost position in pixels.
        @param ypos the y coordinate of the topmost position in pixels.
        @param height the height of the y axis scale in pixels.
    */
    private void createYaxis(Graphics g, int xpos, int ypos, int height)
    {
        int nd = 5;
        int tw ;
        int th = 12;
        int pl, p;
        double vv;
        FontMetrics fm = g.getFontMetrics();
        th=fm.getAscent();
        g.setColor(foreColor);
        for(int i=0; i<=nd; ++i) {
            pl=i*height/nd;
            g.drawLine(xpos, ypos+pl, xpos+3, ypos+pl);
            vv=ymin+(ymax-ymin)*(nd-i)/nd;
            String s=String.format("%3.3g",vv);
            tw=fm.stringWidth(s);

            if(i==0)
                p=pl+th;
            else if(i==nd)
                p=pl;
            else
                p=pl+th/2;
            g.drawString(s,xpos+5,ypos+p);
        }
    }

    /** Draw the y axis with the corresponding scale.
        @param g the graphic context on which to draw.
        @param xpos the x coordinate of the leftmost position in pixels.
        @param ypos the y coordinate of the topmost position in pixels.
        @param height the width of the x axis scale in pixels.
    */
    private void createXaxis(Graphics g, int xpos, int ypos, int width)
    {
        int nd = 4;
        int tw ;
        int th = 12;
        int pl, p;
        double vv;
        FontMetrics fm = g.getFontMetrics();
        th=fm.getAscent();
        g.setColor(foreColor);
        for(int i=0; i<=nd; ++i) {
            pl=i*width/nd;
            g.drawLine(xpos+pl, ypos, xpos+pl, ypos+3);
            vv=xmin+(xmax-xmin)*i/nd;
            String s=String.format("%3.3g",vv);
            tw=fm.stringWidth(s);

            if(i==0)
                p=pl;
            else if(i==nd)
                p=pl-tw;
            else
                p=pl-tw/2;
            g.drawString(s,xpos+p,ypos+th+5);
        }
    }

    /** Create an image containing the colour scale associated to the given
        visualization settings.
        @param colorScaleWidth the with of the colour legend in pixels.
        @param imageHeight the height of the colour legend in pixels.
    */
    private BufferedImage createColorScale(int colorScaleWidth, int imageHeight)
    {
        MatrixStore m=new MatrixStore(MAXPALETTECOL,1,2);

        int i=0;
        for (double dd=min; i<MAXPALETTECOL;dd+=(max-min)/(MAXPALETTECOL+1)) {
            m.setElement(i,0,0,dd);
            m.setElement(i,0,1,dd);
            ++i;
        }
        colorScale  = new BufferedImage(colorScaleWidth+100, imageHeight+20,
                      BufferedImage.TYPE_INT_BGR);

        Graphics2D g = colorScale.createGraphics();
        g.setColor(Color.white);
        g.fillRect(0,0,colorScaleWidth+100, imageHeight+20);
        g.setColor(foreColor);


        FontMetrics fm = g.getFontMetrics();
        int th=fm.getAscent();
        int nd = 5;
        int pl, p;
        double vv;
        for(i=0; i<=nd; ++i) {
            pl=i*imageHeight/nd;
            g.drawLine(colorScaleWidth, pl, colorScaleWidth+3, pl);
            vv=min+(max-min)*(nd-i)/nd;
            if(i==0)
                p=pl+th;
            else if(i==nd)
                p=pl;
            else
                p=pl+th/2;

            String s=String.format("%3.3g",vv);
            g.drawString(s,colorScaleWidth+7,p);
        }

        int old=visType;
        visType=VIS_REAL;
        drawRaster(g, m, colorScaleWidth,imageHeight);
        visType=old;
        g.setColor(foreColor);
        g.drawRect(0,0,colorScaleWidth, imageHeight);

        return colorScale;
    }

    public void setXYRanges(double xm, double xmm, double ym, double ymm)
    {
        xmin=xm;
        xmax=xmm;
        ymin=ym;
        ymax=ymm;
    }

    public void setImageSize(int wx)
    {
        wantedImageSizeX=Math.max(wx, 300);
        wantedImageSizeY=(int)Math.round(wantedImageSizeX*0.8);
    }

    /** The paint callback function for the Swing component.
        @param gg the graphic context.
    */
    protected void paintComponent(Graphics gg)
    {
        int imageWidth = wantedImageSizeX;
        int imageHeight = wantedImageSizeY;
        int colorScaleWidth = 50;

        Graphics2D g2D = (Graphics2D)gg;

        // Activate antialiasing for any drawing operation.
        g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                         RenderingHints.VALUE_ANTIALIAS_ON);

        // Calculate the most adapted separation between the
        // graph and the vertical scale. This correspond roughly to
        // the Y axis plus a small spacing.
        int separatorx = gg.getFontMetrics().stringWidth(" -1.00e-05  ");
        int separatory = gg.getFontMetrics().getAscent()+
            gg.getFontMetrics().getDescent();

        gg.setColor(Color.white);
        gg.fillRect(0,0,2*imageWidth, 2*imageHeight);
        if(hasChanged){
            bufferedImage=createImage(sd, imageWidth, imageHeight);
            if(sd!=null)
                colorScale=createColorScale(colorScaleWidth, imageHeight);
        }
        gg.drawImage(bufferedImage,0,0,null);
        gg.drawImage(colorScale, imageWidth+separatorx, 0, null);
        createXaxis(gg, 0,imageHeight, imageWidth);
        createYaxis(gg, imageWidth,0, imageHeight);

        if(scatterColumnZ>0) {
            String s="Z quote: "+zquote;
            gg.drawString(s, 20, imageHeight+separatory+30);
        }

        String m;
        if(sd!=null) {
            if(isScatter)
                m="Scatter plot: "+sd.getNrow()+" points";
            else
                m="Raster image: "+sd.getNrow()+"x"+sd.getNcol()+" pixels";

            gg.drawString(m, 20, imageHeight+2*separatory+30);
        }
    }

    /** Update the current preprocessing settings.
        @param p the new preprocessing settings.
    */
    public void choosePreprocessing(int p)
    {
        preType=p;
        hasChanged=true;
        repaint();
    }

    /** Apply preprocessing operations to the data to be represented.
        This may include a Fast Fourier Transform operation, direct or inverse.
        @param d the data to be transformed.
        @param direct true if a direct FFT has to be performed. False if an
            inverse FFT has to be done.
    */
    private MatrixStore preprocessing(MatrixStore d, boolean direct)
    {
        double[][] fft2D = new double[d.getNrow()][2*d.getNcol()];

        double r=1.0;
        for(int i=0; i<d.getNrow(); ++i) {
            for(int j=0; j<d.getNcol();++j){
                r=Math.pow(-1.0, i*j);
                fft2D[i][2*j]=d.getElement(i,j,0)*r;
                fft2D[i][2*j+1]=d.getElement(i,j,1)*r;
            }
        }

        DoubleFFT_2D fftDo = new DoubleFFT_2D(d.getNrow(), d.getNcol());

        if(direct) {
            fftDo.complexForward(fft2D);
        } else {
            fftDo.complexInverse(fft2D, false);
        }

        MatrixStore res = new MatrixStore(d.getNrow(),d.getNcol(),2);

        for(int i=0; i<d.getNrow(); ++i) {
            for(int j=0; j<d.getNcol();++j){
                res.setElement(i,j,0,fft2D[i][2*j]);
                res.setElement(i,j,1,fft2D[i][2*j+1]);
            }
        }
        return res;
    }
}
