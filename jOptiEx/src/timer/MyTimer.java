package timer;

// import java.awt.*;
// import java.awt.image.*;
// import javax.swing.*;



public class MyTimer {
    private final long start;

    /** Standard constructor. Time measurement begins from here.

    */
    public MyTimer() {
        start = System.currentTimeMillis();
    }

    /** Get the elapsed time from class construction.
        @return the elapsed time in milliseconds. Time resolution will
        depend on your operating system.
    */
    public long getElapsed() {
        return System.currentTimeMillis() - start;
    }
}
