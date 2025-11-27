package dialogs;

import java.awt.event.*;
import java.awt.*;
import javax.swing.*;

/*  Some routines useful with dialog windows.
    Copyright 2007-2016 by Davide Bucci
*/

public class DialogUtil {

    /** Center the frame on the screen and set its height and width to half
        the screen size */
    public static void center(Window frame)
    {
        GraphicsEnvironment ge =
            GraphicsEnvironment.getLocalGraphicsEnvironment();
        Point center = ge.getCenterPoint();
        //Rectangle bounds = ge.getMaximumWindowBounds();

        int w = frame.getWidth();
        int h = frame.getHeight();

        int x = center.x - w/2, y = center.y - h/2;

        frame.setBounds(x, y, w, h);
  /*      if (w == bounds.width && h == bounds.height)
            frame.setExtendedState(Frame.MAXIMIZED_BOTH);*/
        frame.validate();
    }

    /** Center the frame on the screen and set its height and width to the
        given proportion of the screen size */
    public static void center(Window frame, double propX, double propY)
    {
        GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
        Point center = ge.getCenterPoint();
        Rectangle bounds = ge.getMaximumWindowBounds();
        int w = Math.max((int)(bounds.width*propX), Math.min(frame.getWidth(),
                        bounds.width));
        int h = Math.max((int)(bounds.height*propY), Math.min(frame.getHeight(),
                        bounds.height));
        int x = center.x - w/2, y = center.y - h/2;
        frame.setBounds(x, y, w, h);
    /*    if (w == bounds.width && h == bounds.height)
            frame.setExtendedState(Frame.MAXIMIZED_BOTH);*/
        frame.validate();
    }


    /** Center the frame on the screen and set its height and width to the
        given proportion of the screen size. Specify also the minimum
        size in pixels
    */
    public static void center(Window frame, double propX, double propY,
        int minx, int miny)
    {
        GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
        Point center = ge.getCenterPoint();
        Rectangle bounds = ge.getMaximumWindowBounds();
        int w = Math.max((int)(bounds.width*propX), Math.min(frame.getWidth(),
                        bounds.width));

        if(w<minx)
            w=minx;

        int h = Math.max((int)(bounds.height*propY), Math.min(frame.getHeight(),
                        bounds.height));

        if(h<miny)
            h=miny;

        int x = center.x - w/2, y = center.y - h/2;
        frame.setBounds(x, y, w, h);
    /*    if (w == bounds.width && h == bounds.height)
            frame.setExtendedState(Frame.MAXIMIZED_BOTH);*/
        frame.validate();
    }

    /** Maps the escape key with the action given (ideally, it should probably
        close the dialog window). Example follows:

        <pre>
        // Here is an action in which the dialog is closed

        AbstractAction cancelAction = new AbstractAction ()
        {
            public void actionPerformed (ActionEvent e)
            {
                setVisible(false);
            }
        };
        dialogUtil.addCancelEscape (this, cancelAction);
        </pre>
        @param f the dialog window on which the action should be applied
        @param cancelAction the action to be performed in response to the Esc
            key

    */
    public static void addCancelEscape (JDialog f,
        AbstractAction cancelAction)
    {
        // Map the Esc key to the cancel action description in the dialog box's
        // input map.

        String CANCEL_ACTION_KEY = "CANCEL_ACTION_KEY";

        int noModifiers = 0;

        KeyStroke escapeKey =
            KeyStroke.getKeyStroke (KeyEvent.VK_ESCAPE, noModifiers, false);

        InputMap inputMap = f.getRootPane ().getInputMap
            (JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);

        inputMap.put (escapeKey, CANCEL_ACTION_KEY);

        f.getRootPane ().getActionMap ().put (CANCEL_ACTION_KEY, cancelAction);
    }

    /** Set up the constraints for the GridLayout manager
        @param gridx the x position in the grid
        @param gridy the y position in the grid
        @param width the width of the control
        @param height the heigth of the control
        @param weigthx
        @param weigthy
        @param anch
        @param fill
        @param insets

    */
    public static GridBagConstraints createConst(int gridx, int gridy,
        int width, int height, int weightx, int weighty,
        int anch, int fill, Insets insets)
    {
        GridBagConstraints constraints=new GridBagConstraints();
        constraints.gridx=gridx;
        constraints.gridy=gridy;
        constraints.gridwidth=width;
        constraints.gridheight=height;
        constraints.weightx=weightx;
        constraints.weighty=weighty;
        constraints.fill = fill;
        constraints.anchor = anch;
        constraints.insets=insets;

        return constraints;
    }
}
