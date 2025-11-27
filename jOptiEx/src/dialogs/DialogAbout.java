package dialogs;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import java.io.*;
import javax.imageio.*;
import java.net.*;


/**
    Shows a rather standard "About" dialog. Nothing more exotic than showing
    the nice icon of the program, its name as well as three lines of
    description.

    Copyright 2007-2016 by Davide Bucci
    @author Davide Bucci
*/
public class DialogAbout extends JFrame implements ComponentListener
{
    // The minimu size in pixels.
    private static final int MIN_WIDTH=300;
    private static final int MIN_HEIGHT=250;


    /** Required for the implementation of the ComponentListener interface.
        In this case, prevents from resizing the dialog in a size which is
        too small.
    */
    public void componentResized(ComponentEvent e)
    {
        int width = getWidth();
        int height = getHeight();

        boolean resize = false;
        if (width < MIN_WIDTH) {
            resize = true;
            width = MIN_WIDTH;
         }
         if (height < MIN_HEIGHT) {
            resize = true;
            height = MIN_HEIGHT;
         }
         if (resize) {
            setSize(width, height);
         }
    }
    public void componentMoved(ComponentEvent e)
    {
    }
    public void componentShown(ComponentEvent e)
    {
    }
    public void componentHidden(ComponentEvent e)
    {
    }

    /** Standard constructor: it needs the parent frame.

        @param parent the dialog's parent
    */
    public DialogAbout (JFrame parent)
    {
        // super(parent,"", true);
        super("");
        DialogUtil.center(this, .40,.45,450,350);
        setResizable(false);
        addComponentListener(this);

        // Shows the icon of the program and then three lines read from the
        // resources which describe the software and give the credits.

        GridBagLayout bgl=new GridBagLayout();
        GridBagConstraints constraints=new GridBagConstraints();
        Container contentPane=getContentPane();
        contentPane.setLayout(bgl);

        URL url=DialogAbout.class.getResource(
            "program_icons/imep_lahc.png");
        JLabel icon=new JLabel("");
        constraints.weightx=100;
        constraints.weighty=100;
        constraints.gridx=0;
        constraints.gridy=0;
        constraints.gridwidth=1;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.CENTER;
        constraints.insets=new Insets(10,20,0,20);

        if (url != null) icon.setIcon(new ImageIcon(url));
        contentPane.add(icon, constraints);


        JLabel programName=new JLabel("jOptiEx");

        Font f=new Font("Lucida Grande",Font.BOLD,18);

        programName.setFont(f);
        constraints.gridx=0;
        constraints.gridy=1;
        constraints.gridwidth=1;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.CENTER;
        constraints.insets=new Insets(0,20,0,20);
        contentPane.add(programName, constraints);


        JLabel programVersion=new JLabel("Version 1.2.1");

        constraints.gridx=0;
        constraints.gridy=2;
        constraints.gridwidth=1;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.CENTER;
        contentPane.add(programVersion, constraints);

        JLabel programDescription1=new JLabel("A data simulation explorer");
        constraints.gridx=0;
        constraints.gridy=3;
        constraints.gridwidth=1;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.CENTER;
        contentPane.add(programDescription1, constraints);

        JLabel programDescription2=new JLabel(
                "Davide Bucci, June 30, 2016");
        constraints.gridx=0;
        constraints.gridy=4;
        constraints.gridwidth=3;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.CENTER;
        constraints.insets=new Insets(0,20,20,20);
        contentPane.setBackground(Color.white);
        contentPane.add(programDescription2, constraints);

        JLabel programDescription3=new JLabel("bucci@minatec.grenoble-inp.fr");
        constraints.gridx=0;
        constraints.gridy=5;
        constraints.gridwidth=1;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.CENTER;
        constraints.insets=new Insets(0,20,20,20);

        contentPane.add(programDescription3, constraints);
    }
}
