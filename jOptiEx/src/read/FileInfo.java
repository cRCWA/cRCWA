package read;
import storage.*;

public class FileInfo {
    public double xmin;
    public double xmax;
    public double ymin;
    public double ymax;
    public MatrixStore mm;

    public FileInfo()
    {
        mm=null;
        xmin=xmax=ymin=ymax=0;
    }
}

