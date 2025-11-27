package read;

import java.util.*;
import java.io.*;
import java.util.regex.*;
import java.awt.*;

import storage.*;

public class ReadVariousFormats{

    public static String read10lines(String filename)
        throws FileNotFoundException
    {
        String s="";
        try {
            Scanner scanner = new Scanner(new File(filename));
            for(int i=0; i<9;++i)
                s+=scanner.nextLine()+"\n";
            s+=scanner.nextLine();
        } catch (InputMismatchException E) {
            return "";
        } catch (NoSuchElementException F) {
            return "";
        }

        return s;
    }
    public static int getGnuplotColumn(String filename)
        throws FileNotFoundException
    {
        int ncol=-1;
        try {
            Scanner scanner = new Scanner(new File(filename));
            // In any case, an Optiwave file will use the dot as
            // the decimal separator
            scanner.useLocale(Locale.US);

            // Skip two lines
            scanner.nextLine();
            scanner.nextLine();

            // Read a line
            String s=scanner.nextLine();
            // Now the scan is done on the line just read.
            scanner = new Scanner(s);

            scanner.useLocale(Locale.US);

            // Here we know the size and the structure of the file
            // we can thus create the matrix to store the data and read
            // the file.
            scanner.useDelimiter(Pattern.compile(" |\t|\n|\r"));

            double x;
            int i=-1;
            //System.out.println(s);
            // We try to read a line, to see if it is ok.
            try {
                for(i=0; i<10; ++i) {
                    x=scanner.nextDouble(); // x
                }
            } catch (InputMismatchException E) {

            } catch (NoSuchElementException F) {

            }
            ncol = i;
        } catch (InputMismatchException E) {
            return -1;
        } catch (NoSuchElementException F) {
            return -1;
        }
            // All tests have been passed!

        return ncol;

    }

    public static boolean isGnuplot(String filename)
        throws FileNotFoundException
    {
        try {
            BufferedReader r=new BufferedReader(new FileReader(filename));
            CustomTokenizer scanner = new CustomTokenizer(r);
            // Skip two lines
            scanner.nextLine();
            scanner.nextLine();

            // We try to read a line, to see if it has at least six doubles.
            scanner.nextDouble();
            scanner.nextDouble();
            scanner.nextDouble();
            scanner.nextDouble();
            scanner.nextDouble();
            scanner.nextDouble();
        } catch (InputMismatchException E) {
            return false;
        } catch (NoSuchElementException F) {
            return false;
        } catch (IOException G) {
            return false;
        }
            // All tests have been passed!
        return true;

    }

    public static boolean isOptiWave(String fileName)
        throws FileNotFoundException
    {
        try {
            File f=new File(fileName);
            Scanner scanner = new Scanner(f);

            // In any case, an Optiwave file will use the dot as
            // the decimal separator
            scanner.useLocale(Locale.US);

            scanner.nextLine();
            int nx = scanner.nextInt();
            int ny = scanner.nextInt();
            double xmin = scanner.nextDouble();
            double xmax = scanner.nextDouble();
            double ymin = scanner.nextDouble();
            double ymax = scanner.nextDouble();
            // If we got to this point, the header of the file does appear
            // to be not too distant from a standard Optiwave one. We thus
            // need to check if the data has some sense.

            if (nx<=0 || ny<=0)
                return false;
            if (xmin>xmax)
                return false;
            if (ymin>ymax)
                return false;

        } catch (InputMismatchException E) {
            return false;
        } catch (NoSuchElementException F) {
            return false;
        }
            // All tests have been passed!
        return true;
    }

    public static FileInfo readOptiWaveFast(String fileName)
        throws FileNotFoundException,IOException
    {
        Reader r=null;

        MatrixStore dd=null;
        FileInfo f=new FileInfo();
        try {
            r= new BufferedReader(new FileReader(fileName),1<<16);
            CustomTokenizer scanner = new CustomTokenizer(r);
            scanner.nextLine();
            int nx = scanner.nextInt();
            int ny = scanner.nextInt();
            double xmin = scanner.nextDouble();
            double xmax = scanner.nextDouble();
            double ymin = scanner.nextDouble();
            double ymax = scanner.nextDouble();

            dd = new MatrixStore(ny, nx,2);
            double v;

            for(int j=0;j<ny;++j) {
                for(int i=0;i<nx;++i) {
                    v=scanner.nextDouble();
                    dd.setElement(j,i,0,v);
                    v=scanner.nextDouble();
                    dd.setElement(j,i,1,v);
                }
            }
            f.xmin=xmin;
            f.xmax=xmax;
            f.ymin=ymin;
            f.ymax=ymax;
            f.mm=dd;
        } catch (IOException E){

        } finally {
            if (r != null)
                r.close();
        }
        return f;
    }


    private boolean stopImmediately;
    private boolean reading;

    public MatrixStore readGnuplot(String fileName,  int ncol)
        throws FileNotFoundException, IOException
    {
        Reader r=null;
        int i=0;
        if(ncol<0) return null;

        MatrixStore dd=new MatrixStore(100,ncol,1);
        try {
            r = new BufferedReader(new FileReader(fileName),1<<18);
            CustomTokenizer scanner = new CustomTokenizer(r);

            double z;

            // Since we do not know 'a priori' the file size,
            // we gradually increase the size of the matrix used
            // to store the data.

            synchronized (dd){
                stopImmediately = false;
                reading=true;
                while(!stopImmediately) {
                    for(int k=0; k<ncol;++k) {
                        z = scanner.nextDouble();
                        dd.increaseSize(i+1, ncol, 1);
                        dd.setElement(i,k,0, z);
                    }
                    ++i;
                }
                reading=false;
            }
        } catch (IOException E){

        } finally {
            reading=false;
            if (r != null)
                r.close();
        }
        reading=false;
        if(stopImmediately)
            dd=null;
        stopImmediately = false;
        return dd;

    }
    public boolean isReading()
    {
        System.out.println("reading="+reading);
        return reading;
    }
    public void stopReading()
    {
        System.out.println("STOP READING!");
        stopImmediately=true;
    }
}


