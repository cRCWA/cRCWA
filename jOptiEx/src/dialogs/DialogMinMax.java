package dialogs;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/** Choose min and max ranges for a drawing.
    Davide Bucci 2015-2016
*/
public class DialogMinMax extends JDialog implements ComponentListener
{
    private static final int MIN_WIDTH=400;
    private static final int MIN_HEIGHT=250;

    private final JCheckBox calcMinCB;
    private final JTextField calcMinTF;
    private final JCheckBox calcMaxCB;
    private final JTextField calcMaxTF;

    private boolean isOk=false;

    /** Standard constructor: it needs the parent frame.
        @param parent the dialog's parent
    */
    public DialogMinMax (JFrame parent)
    {
        super(parent,"Use min/max in raster drawings", true);
        addComponentListener(this);

        // Ensure that under MacOSX >= 10.5 Leopard, this dialog will appear
        // as a document modal sheet

        getRootPane().putClientProperty("apple.awt.documentModalSheet",
                Boolean.TRUE);

        GridBagLayout bgl=new GridBagLayout();
        GridBagConstraints constraints=new GridBagConstraints();
        Container contentPane=getContentPane();
        contentPane.setLayout(bgl);

        constraints.insets.right=30;

        JLabel empty=new JLabel("  ");
        constraints.weightx=100;
        constraints.weighty=100;
        constraints.gridx=0;
        constraints.gridy=0;
        constraints.gridwidth=1;
        constraints.gridheight=1;
        constraints.fill=GridBagConstraints.BOTH;

        contentPane.add(empty, constraints);            // Add "   " label

        JLabel empty1=new JLabel("  ");
        constraints.weightx=100;
        constraints.weighty=100;
        constraints.gridx=3;
        constraints.gridy=0;
        constraints.gridwidth=1;
        constraints.gridheight=1;
        constraints.fill=GridBagConstraints.BOTH;
        contentPane.add(empty1, constraints);           // Add "   " label

        calcMinCB=new JCheckBox("Calculate automatically minimum value");
        constraints.gridx=1;
        constraints.gridy=0;
        constraints.gridwidth=2;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.WEST;
        constraints.fill=GridBagConstraints.BOTH;
        contentPane.add(calcMinCB, constraints);
        calcMinCB.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e)
            {
                if (!calcMinCB.isSelected())
                    calcMinTF.setEnabled(true);
                else
                    calcMinTF.setEnabled(false);
            }
        });


        calcMinTF =new JTextField("");
        constraints.gridx=1;
        constraints.gridy=1;
        constraints.gridwidth=2;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.WEST;
        constraints.fill=GridBagConstraints.HORIZONTAL;
        contentPane.add(calcMinTF, constraints);

        calcMaxCB=new JCheckBox("Calculate automatically maximum value");
        constraints.gridx=1;
        constraints.gridy=2;
        constraints.gridwidth=2;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.WEST;
        constraints.fill=GridBagConstraints.BOTH;
        contentPane.add(calcMaxCB, constraints);

        calcMaxCB.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e)
            {
                if (!calcMaxCB.isSelected())
                    calcMaxTF.setEnabled(true);
                else
                    calcMaxTF.setEnabled(false);
            }
        });

        calcMaxTF =new JTextField("");
        constraints.gridx=1;
        constraints.gridy=3;
        constraints.gridwidth=2;
        constraints.gridheight=1;
        constraints.fill=GridBagConstraints.HORIZONTAL;
        constraints.anchor=GridBagConstraints.WEST;
        contentPane.add(calcMaxTF, constraints);

        // Put the OK and Cancel buttons and make them active.
        JButton ok=new JButton("Ok");
        JButton cancel=new JButton("Cancel");

        constraints.gridx=0;
        constraints.gridy=4;
        constraints.gridwidth=4;
        constraints.gridheight=1;
        constraints.anchor=GridBagConstraints.EAST;

        Box b=Box.createHorizontalBox();
        b.add(Box.createHorizontalGlue());
        ok.setPreferredSize(cancel.getPreferredSize());

        b.add(cancel);
        b.add(Box.createHorizontalStrut(12));
        b.add(ok);

        ok.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent evt)
            {
                try{
                    Double.parseDouble(calcMinTF.getText());
                    Double.parseDouble(calcMaxTF.getText());
                    isOk=true;
                    setVisible(false);
                } catch (java.lang.NumberFormatException J) {
                    JOptionPane.showMessageDialog(null,
                        "I could not understand the numbers you wrote.",
                        "Woha!", JOptionPane.ERROR_MESSAGE);
                }

            }
        });
        cancel.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent evt)
            {
                isOk=false;
                setVisible(false);
            }
        });

        // Here is an action in which the dialog is closed
        AbstractAction cancelAction = new AbstractAction ()
        {
            public void actionPerformed (ActionEvent e)
            {
                setVisible(false);
            }
        };
        contentPane.add(b, constraints);        // Add OK/cancel dialog

        DialogUtil.addCancelEscape(this, cancelAction);
        pack();
        DialogUtil.center(this);
        getRootPane().setDefaultButton(ok);
    }

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

    public void setMinMaxValues(double min, double max)
    {
        calcMinTF.setText(""+min);
        calcMaxTF.setText(""+max);
    }

    public boolean getMinCB()
    {
        return calcMinCB.isSelected();
    }

    public boolean getMaxCB()
    {
        return calcMaxCB.isSelected();
    }

    public double getMinValue()
    {
        return Double.parseDouble(calcMinTF.getText());
    }

    public double getMaxValue()
    {
        return Double.parseDouble(calcMaxTF.getText());
    }

    public void setMinMaxCB(boolean d, boolean f)
    {
        calcMinCB.setSelected(d);
        if (!calcMaxCB.isSelected())
            calcMaxTF.setEnabled(true);
        else
            calcMaxTF.setEnabled(false);
        calcMaxCB.setSelected(f);
        if (!calcMaxCB.isSelected())
            calcMaxTF.setEnabled(true);
        else
            calcMaxTF.setEnabled(false);
    }


    public void componentMoved(ComponentEvent e)
    {
        // Nothing to do
    }
    public void componentShown(ComponentEvent e)
    {
        // Nothing to do
    }
    public void componentHidden(ComponentEvent e)
    {
        // Nothing to do
    }

    public boolean getIsOk()
    {
        return isOk;
    }
}