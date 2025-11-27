package storage;

/** The MatrixStore class contains the contents of a matrix, but also the
    data about the x and y range

*/
public class MatrixStore {
    // Store the matrix size
    private int ncol;
    private int nrow;

    private double[][][] mat;
    private int allocatedRows;
    private int allocatedColumns;
    private int allocatedPlanes;

    // Store the x and y range
    private int xmin;
    private int ymin;
    private int xmax;
    private int ymax;

    private int nPlanes;
    private int actualPlane;

    public int getScanline(int r, int k)
    {
        double y=0, oy=0;
        int i;
        for(i=0; i<getNrow();++i) {
            y=getElement(i,r,k);
            if(y!=oy&&i>0) {
                break;
            }
            oy=y;
        }
        return i;
    }


    public double getMax()
    {
        double max = mat[actualPlane][0][0];

        for(int i=0; i<nrow; ++i) {
            for(int j=0; j<ncol; ++j) {
                if(mat[actualPlane][i][j]>max)
                    max=mat[actualPlane][i][j];
            }
        }
        return max;
    }

    public void increaseSize(int N, int M, int K)
    {
        // We check if we need a bigger matrix
        if(N>allocatedRows || M>allocatedColumns || K>allocatedPlanes) {

            int P  = N*160/100;     // We take some margin.

            // We allocate the new space
            double[][][] n_mat = new double[K][P][M];
            allocatedRows = P;
            allocatedColumns = M;
            allocatedPlanes = K;

            // We copy the old values
            for(int k=0; k<nPlanes; ++k) {
                for(int i=0; i<nrow; ++i) {
                    for(int j=0; j<ncol; ++j) {
                        n_mat[k][i][j]=mat[k][i][j];
                    }
                }
            }
            // We discard the old values.
            mat = n_mat;
        }
        nrow=N;
        ncol=M;
        nPlanes=K;
    }

    public double getMin()
    {
        double min = mat[actualPlane][0][0];

        for(int i=0; i<nrow; ++i) {
            for(int j=0; j<ncol; ++j) {
                if(mat[actualPlane][i][j]<min)
                    min=mat[actualPlane][i][j];
            }
        }

        return min;
    }

    /** Standard constructor: create a matrix of the given size
    */
    public MatrixStore(int row, int col, int np)
    {
        ncol = col;
        nrow = row;
        nPlanes = np;
        allocatedRows=nrow;
        allocatedPlanes=nPlanes;
        allocatedColumns=ncol;
        mat = new double[nPlanes][nrow][ncol];
    }

    /** Set the current plane
    */
    public void setPlane(int p)
    {
        if(p<0 || p>=nPlanes)
            return;

        actualPlane=p;
    }

    /** Get the number of available planes.
    */
    public int getNPlanes()
    {
        return nPlanes;
    }

    /** Standard constructor: create a matrix of the given size
    */
    public MatrixStore(int row, int col)
    {
        ncol = col;
        nrow = row;
        nPlanes = 1;
        mat = new double[nPlanes][nrow][ncol];
    }

    /** Get the number of columns stored
    */
    public int getNcol()
    {
        return ncol;
    }
    /** Get the number of rows stored
    */
    public int getNrow()
    {
        return nrow;
    }

    public double getElement(int i, int j)
    {
        return getElement(i,j,actualPlane);
    }
    public double getElement(int i, int j, int np)
    {
        if (i>=nrow || i<0 || j>=ncol || j<0||np<0||np>=nPlanes) {
            throw new RuntimeException("Incorrect matrix indexes "+i+" "+j);
        }

        return mat[np][i][j];
    }

    public void setElement(int i, int j, double v)
    {
        setElement(i,j,actualPlane, v);
    }
    public void setElement(int i, int j, int np, double v)
    {
        /*if (i>=nrow || i<0 || j>=ncol || j<0||np<0||np>=nPlanes) {
            throw new RuntimeException("Incorrect matrix indexes "+i+" "+j);
        }*/

        mat[np][i][j] = v;

    }

    public double getXmin()
    {
        return xmin;
    }

    public double getXmax()
    {
        return xmax;
    }

    public double getYmin()
    {
        return ymin;
    }

    public double getYmax()
    {
        return ymax;
    }
}
