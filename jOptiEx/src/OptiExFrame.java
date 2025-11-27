import javax.swing.*;
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.image.*;
import javax.swing.event.*;
import javax.swing.text.*;
import java.lang.*;
import java.awt.event.*;
import java.lang.reflect.*;
import java.nio.file.*;


import storage.*;
import visualization.*;
import read.*;
import timer.*;
import dialogs.*;

/**
The main frame of jOptiEx

    Copyright 2012-2017 by Davide Bucci

    @author Davide Bucci

*/

class OptiExFrame extends JFrame implements ItemListener,
      ActionListener
{
    boolean isDebug;
    public GraphingPanel gp;
    JList<String> fileJList;
    File currentPath;
    private String filename;
    long date;

    final static String DIR_TAG  = "D  ";
    final static String FILE_TAG = "   ";
    final static int N_MAX_COL = 10;
    final static int MAX_FILE_SIZE = 220000000;

    ListSelectionListener fileListSelectionListener;
    String fileFormat;
    boolean showHiddenFiles;

    JPanel postProcessing;
    JList<String> postColumnX;
    JList<String> postColumnY;
    JList<String> postColumnZ;
    JList<String> postColumnR;
    JList<String> postColumnI;

    JSplitPane splitPaneTexts;

    JScrollPane ts;
    JSlider zquotaSlider;
    JSlider pointSizeSlider;
    JCheckBox automaticPointSize;
    int columnX;
    int columnY;
    int columnZ;
    int columnR;
    int columnI;

    private JTextArea textFile;

    private JRadioButton checkReal;
    private JRadioButton checkImag;
    private JRadioButton checkMagnitude;
    private JRadioButton checkPhase;
    private JRadioButton checkLogMag;

    private JRadioButton checkNoPreproc;
    private JRadioButton checkFFT;
    private JRadioButton checkIFFT;

    private ReadVariousFormats reader;

    private JPanel visualization;
    private JPanel preproc;

    private FileInfo f;

    // Watcher service (in package java.nio): notification of file/dir changes
    private WatchService watcher;
    private Path watchDir;
    private Path watchFile;

    JTextField fileNameField;
    static private int openWindowsNumber;

    static private boolean weAreOnAMac;
    // shortcut key to be used:
    static private int shortcutKey;
    // META (Command) for Macintoshes
    // CTRL elsewhere
    static private boolean useNativeFileDialogs;

    /** Standard constructor
    */
    OptiExFrame()
    {
        super("IMEP-LAHC optical simulation browser");
        DialogUtil.center(this, .75,.75,800,800);
        ++openWindowsNumber;
        setDefaultCloseOperation(
                javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);

        useNativeFileDialogs=false;
        if (System.getProperty("os.name").startsWith("Mac")) {
            weAreOnAMac=true;
            useNativeFileDialogs=true;
            // From what I know, only Mac users expect to use the Command (meta)
            // key for shortcuts, while others will use Control.
            shortcutKey=InputEvent.META_MASK;

            System.setProperty("apple.laf.useScreenMenuBar","true");
            // Here we use the reflection provided by Java to understand
            // if the AppleSpecific class is available on the system.
            // This class should be compiled separately from the main
            // program since the compilation can be successful only on
            // a MacOSX system.
            boolean notCompiledForApple=false;
            try {
                Class a = Class.forName("AppleSpecific");
                Object b = a.newInstance();
                Method m = a.getMethod("answerFinder");
                m.invoke(b);
            } catch (InstantiationException exc) {
                notCompiledForApple=true;
            } catch (IllegalAccessException exc1) {
                notCompiledForApple=true;
            } catch (ClassNotFoundException exc2) {
                notCompiledForApple=true;
            } catch (NoSuchMethodException exc3) {
                notCompiledForApple=true;
            } catch (InvocationTargetException exc4) {
                notCompiledForApple=true;
            }

            if(notCompiledForApple) {
               weAreOnAMac = false;
               System.out.println("It seems that this software has been "+
                       "compiled on a system different from MacOSX. Some nice "+
                       "integrations with MacOSX will therefore be absent. If "+
                       "you have compiled on MacOSX, make sure you used the "+
                       "'compile' or 'rebuild' script along with the 'mac' "+
                       "option.");
            }
        } else {
            weAreOnAMac=false;
        }
        isDebug = true;
        gp = new GraphingPanel();
        fileFormat ="";
        gp.measureTimings(isDebug);
        showHiddenFiles = false;
        filename="";
        watcher=null;
        try{ // TODO: implement a file watcher
            watcher = FileSystems.getDefault().newWatchService();
        } catch (Exception E)
        {}
        //keys = new HashMap<WatchKey,Path>();
    }

    /**
     * Register the given directory with the WatchService
     */
    private void register(Path dir) throws IOException {
        WatchKey key = dir.register(watcher,
            StandardWatchEventKinds.ENTRY_CREATE,
            StandardWatchEventKinds.ENTRY_DELETE,
            StandardWatchEventKinds.ENTRY_MODIFY);
    }

    public void init()
    {
        addWindowListeners();

        fileJList = new JList<String>();

        fileJList.setPrototypeCellValue(
            "   boule_bucci_r6um_resTE_Ey_m_14.mode   ");
        fileJList.setCellRenderer(new ListFilesRenderer());

        fileJList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        fileListSelectionListener =
            new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent listSelectionEvent)
                {
                    String elem=(String)fileJList.getSelectedValue();

                    if(listSelectionEvent.getValueIsAdjusting())
                        return;
                    // Now, we must identify wether the user has selected a file
                    // or a directory.
                    if (elem==null)
                        return;
                    if(elem.startsWith(FILE_TAG)) {
                        // It is a file: we try to read it.
                        OpenFile openf=new OpenFile(OptiExFrame.this,
                                elem.substring(3));
                        Thread thread = new Thread(openf);
                        //thread.setDaemon(true);
                        // Start the thread
                        thread.start();
                        //setFileToBeRead(elem.substring(3));
                    } else if (elem.startsWith(DIR_TAG)) {

                        String s=elem.substring(3);
                        if ("..".equals(s)) {
                            s=currentPath.getAbsolutePath();
                            int p=s.lastIndexOf(File.separatorChar);

                            if(p>=0)
                                currentPath=new File(s.substring(0,p));
                        } else {
                            currentPath=new File(currentPath, s);
                        }
                        updateFileList();
                        fileNameField.setText(currentPath.getAbsolutePath());
                    }
                    gp.repaint();
                }
            };

        currentPath=new File(".");
        String s=currentPath.getAbsolutePath();
        int p=s.lastIndexOf(File.separatorChar);
        if(p>=0)
            currentPath=new File(s.substring(0,p));

        updateFileList();
        fileJList.addListSelectionListener(fileListSelectionListener);
        Container cp = getContentPane();

        fileNameField = new JTextField(80);
        fileNameField.setText(currentPath.getPath());
        //fileNameField.setEditable(false);

        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);

        Dimension windowSize = getSize();
        gp.setPreferredSize(new Dimension(windowSize.width*85/100,
            windowSize.width*60/100));

        JScrollPane jsg = new JScrollPane(createVisualizationPanel());
        JScrollPane jsl = new JScrollPane(fileJList);
        splitPane.setTopComponent(jsg);
        splitPane.setBottomComponent(jsl);
        splitPane.setResizeWeight(.8);

        postProcessing = createPostProcessingPanel();
        splitPaneTexts = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        ts = new JScrollPane(textFile);
        splitPaneTexts.setTopComponent(ts);
        splitPaneTexts.setBottomComponent(new JScrollPane(postProcessing));

        splitPaneTexts.setResizeWeight(.8);

        JSplitPane splitPane3 = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        splitPane3.setResizeWeight(.8);

        splitPane3.setTopComponent(splitPane);

        JPanel bottom=new JPanel(new GridBagLayout());
        GridBagConstraints d = new GridBagConstraints();
        int position=0;
        d.gridx=0;
        d.gridheight=1;
        d.gridwidth=1;
        d.gridy=position++;
        d.weightx=1.0;
        d.weighty=0;
        d.fill=GridBagConstraints.BOTH;

        bottom.add(createModPanel(),d);

        d.gridy=position++;
        d.weightx=1.0;
        d.weighty=1.0;
        bottom.add(splitPaneTexts,d);

        JScrollPane bb=new JScrollPane(bottom);
        bb.setHorizontalScrollBarPolicy(
            ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        splitPane3.setBottomComponent(bottom);

        cp.add(splitPane3);

        setJMenuBar(createMenu());

        updateScatter(0,1,2,3,4);
    }

    public void addWindowListeners()
    {
        /*  Add a window listener to close the application when the frame is
            closed. This behavior is platform dependent, for example a
            Macintosh application can be made run without a visible frame.
            There would anyway the need to customize the menu bar, in order
            to allow the user to open a new FidoFrame, when it has been
            closed once. The easiest solution to implement is therefore to
            make the application close when the user closes the last frame.
         */
        addWindowListener(new WindowAdapter(){
                public void windowClosing(WindowEvent e)
                {
                    setVisible(false);
                    dispose();

                    --openWindowsNumber;

                    if (openWindowsNumber<1)
                        System.exit(0);
                }
        });

        /* Add a focus listener. Mainly to refresh the file list when the
           software is activated.
        */
        addWindowFocusListener(new WindowAdapter() {
            public void windowGainedFocus(WindowEvent e)
            {
                String sname=(String)fileJList.getSelectedValue();

                //updateFileList();
                if(!filename.equals("")) {
                    File file = new File(filename);
                    if(file.exists()) {
                        // Check if the file has changed
                        Long lastModified = file.lastModified();
                        if(lastModified!=date) {
                            gp.forceRedraw();
                        }
                        if(e!=null && sname!=null) {
                            int pos=fileJList.getNextMatch(sname,0,
                                Position.Bias.Forward);
                            if(isDebug && sname!=null)
                                System.out.println(sname+" Select: "+pos);

                            // fileJList.setSelectedIndex(pos);
                        }

                    } else {
                        setFileToBeRead("");
                        if(isDebug)
                            System.out.println("File deleted?");
                    }
                }
            }
        });
    }

    /** Create the menu bar of the program.
    */
    public JMenuBar createMenu()
    {
    // Menu creation
        JMenuBar menuBar=new JMenuBar();

        JMenu fileMenu=new JMenu("File");
        JMenuItem fileNew = new JMenuItem("New");
        JMenuItem fileOpen = new JMenuItem("Open");
        JMenuItem fileClose = new JMenuItem("Close");
        JMenuItem fileQuit = new JMenuItem("Quit");
        JMenuItem fileAbout = new JMenuItem("About");

        fileMenu.add(fileNew);
        fileMenu.add(fileOpen);
        fileMenu.add(fileClose);

        fileNew.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_N,
                    shortcutKey));
        fileOpen.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O,
                    shortcutKey));
        fileClose.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_W,
                    shortcutKey));


        if(!weAreOnAMac) {
            fileMenu.addSeparator();
            fileMenu.add(fileAbout);
            fileAbout.addActionListener((ActionListener)this);
            fileMenu.addSeparator();
            fileMenu.add(fileQuit);
        }
        fileNew.addActionListener((ActionListener)this);
        fileOpen.addActionListener((ActionListener)this);
        fileQuit.addActionListener((ActionListener)this);
        fileClose.addActionListener((ActionListener)this);

        menuBar.add(fileMenu);

        JMenu viewMenu=new JMenu("View");
        JMenuItem grayscalePalette = new JMenuItem("Grayscale palette");
        JMenuItem jetPalette = new JMenuItem("Jet palette");
        JMenuItem blueRedPalette = new JMenuItem("Blue/White/Red palette");
        JMenuItem tempPalette = new JMenuItem("Temperature palette");
        JMenuItem setMinMaxRange = new JMenuItem("Set min/max color range");
        JMenuItem readAgain = new JMenuItem("Refresh");

        viewMenu.add(readAgain);
        viewMenu.addSeparator();
        viewMenu.add(grayscalePalette);
        viewMenu.add(jetPalette);
        viewMenu.add(blueRedPalette);
        viewMenu.add(tempPalette);
        viewMenu.addSeparator();
        viewMenu.add(setMinMaxRange);

        readAgain.addActionListener((ActionListener)this);
        grayscalePalette.addActionListener((ActionListener)this);
        jetPalette.addActionListener((ActionListener)this);
        blueRedPalette.addActionListener((ActionListener)this);
        tempPalette.addActionListener((ActionListener)this);

        setMinMaxRange.addActionListener((ActionListener)this);


        menuBar.add(viewMenu);

        return menuBar;
    }

    /** Create a post-processing panel containing all the choices on how to
        assign the columns.
    */
    public JPanel createPostProcessingPanel()
    {
        JPanel postProcessing = new JPanel();
        textFile = new JTextArea();
        textFile.setEditable(false);

        columnX=0;
        columnY=1;
        columnZ=-1;
        columnR=2;
        columnI=3;
        Vector<String> v=new Vector<String>();

        for(int i=0; i<10;++i)
            v.add("dummy");

        postColumnX=new JList<String>(v);
        postColumnY=new JList<String>();
        postColumnZ=new JList<String>();
        postColumnR=new JList<String>();
        postColumnI=new JList<String>();

        postColumnX.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        postColumnY.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        postColumnZ.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        postColumnR.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        postColumnI.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);


        postColumnX.addListSelectionListener(
            new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent
                    listSelectionEvent)
                {
                    columnX=postColumnX.getSelectedIndex()-1;

                    gp.setScatterType(columnX,columnY,columnZ,
                        columnR,columnI);
                    gp.repaint();
                }
            }
        );

        postColumnY.addListSelectionListener(
            new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent
                    listSelectionEvent)
                {
                    columnY=postColumnY.getSelectedIndex()-1;
                    gp.setScatterType(columnX,columnY,columnZ,
                        columnR,columnI);
                    gp.repaint();
                }
            }
        );
        postColumnZ.addListSelectionListener(
            new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent
                    listSelectionEvent)
                {
                    columnZ=postColumnZ.getSelectedIndex()-1;
                    if(columnZ>0)
                        zquotaSlider.setVisible(true);
                    else
                        zquotaSlider.setVisible(false);
                    gp.setScatterType(columnX,columnY,columnZ,
                        columnR,columnI);
                    gp.repaint();
                }
            }
        );

        postColumnR.addListSelectionListener(
            new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent
                    listSelectionEvent)
                {
                    columnR=postColumnR.getSelectedIndex()-1;
                    gp.setScatterType(columnX,columnY,columnZ,
                        columnR,columnI);
                    gp.repaint();
                }
            }
        );

        postColumnI.addListSelectionListener(
            new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent
                    listSelectionEvent)
                {
                    columnI=postColumnI.getSelectedIndex()-1;
                    if(columnI>0)
                        visualization.setVisible(true);
                    else
                        visualization.setVisible(false);

                    gp.setScatterType(columnX,columnY,columnZ,
                        columnR,columnI);
                    gp.repaint();
                }
            }
        );

        postColumnX.setPrototypeCellValue("1234567890");
        postColumnY.setPrototypeCellValue("1234567890");
        postColumnZ.setPrototypeCellValue("1234567890");
        postColumnR.setPrototypeCellValue("1234567890");
        postColumnI.setPrototypeCellValue("1234567890");

        GridBagConstraints c = new GridBagConstraints();
        c.gridy=0;
        c.anchor=GridBagConstraints.WEST;
        c.insets=new Insets(0,5,5,0);
        postProcessing.setLayout(new GridBagLayout());
        postProcessing.add(new JLabel("X data"),c);
        postProcessing.add(new JLabel("Y data"),c);
        postProcessing.add(new JLabel("Z data"),c);
        postProcessing.add(new JLabel("Real data"),c);
        postProcessing.add(new JLabel("Imag data"),c);
        c.gridy=1;
        postProcessing.add(postColumnX,c);
        postProcessing.add(postColumnY,c);
        postProcessing.add(postColumnZ,c);
        postProcessing.add(postColumnR,c);
        postProcessing.add(postColumnI,c);

        postProcessing.setVisible(false);
        return postProcessing;
    }

    /** Create a visualization panel populated by all those interface elements
        needed for its use (select z slice, point size, real, imaginary...).
    */
    public JPanel createVisualizationPanel()
    {

        return gp;
    }

    public JPanel createModPanel()
    {
        GridBagConstraints d = new GridBagConstraints();
        JPanel pp = new JPanel(new GridBagLayout());
        int position=0;
        d.gridx=0;
        d.gridheight=1;
        d.gridwidth=2;
        d.gridy=position++;
        d.weightx=1;
        d.fill=GridBagConstraints.BOTH;
        pp.add(fileNameField,d);
        zquotaSlider = new JSlider(JSlider.HORIZONTAL,0,100,0);
        d.gridy=position++;
        pp.add(zquotaSlider,d);
        zquotaSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                JSlider source = (JSlider)e.getSource();
                if (!source.getValueIsAdjusting()) {
                    gp.setZpercent((int)source.getValue());
                }
            }
        });
        pointSizeSlider =  new JSlider(JSlider.HORIZONTAL,1,25,5);
        automaticPointSize = new JCheckBox("Automatic point size");
        d.gridy=position++;
        d.gridheight=1;
        d.gridwidth=1;
        d.gridx=0;
        pp.add(automaticPointSize,d);
        automaticPointSize.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                gp.setAutomaticPointSize(automaticPointSize.isSelected());
                pointSizeSlider.setEnabled(!automaticPointSize.isSelected());
            }
        });
        d.gridx=1;
        pp.add(pointSizeSlider,d);
        pointSizeSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                JSlider source = (JSlider)e.getSource();
                if (!source.getValueIsAdjusting()) {
                    int p=(int)source.getValue();
                    gp.setSizePoint(p,p);
                }
            }
        });
        automaticPointSize.setSelected(true);
        pointSizeSlider.setEnabled(false);


        checkReal = new JRadioButton("Real");
        checkImag = new JRadioButton("Imag");
        checkMagnitude = new JRadioButton("Magnitude");
        checkLogMag = new JRadioButton("ln(Magnitude)");
        checkPhase = new JRadioButton("Phase");


        checkReal.setSelected(true);
        ButtonGroup gr = new ButtonGroup();
        visualization = new JPanel();
        visualization.add(checkReal);
        visualization.add(checkImag);
        visualization.add(checkMagnitude);
        visualization.add(checkLogMag);
        visualization.add(checkPhase);

        gr.add(checkReal);
        gr.add(checkImag);
        gr.add(checkMagnitude);
        gr.add(checkLogMag);
        gr.add(checkPhase);

        ButtonGroup gr1 = new ButtonGroup();
        preproc = new JPanel();
        checkNoPreproc=new JRadioButton("RAW");
        checkFFT = new JRadioButton("FFT");
        checkIFFT = new JRadioButton("IFFT");

        checkNoPreproc.setSelected(true);
        preproc.add(checkNoPreproc);
        preproc.add(checkFFT);
        preproc.add(checkIFFT);

        gr1.add(checkNoPreproc);
        gr1.add(checkFFT);
        gr1.add(checkIFFT);

        d.gridy=position++;
        pp.add(visualization,d);
        d.gridy=position++;
        pp.add(preproc,d);

        checkReal.addItemListener(this);
        checkImag.addItemListener(this);
        checkMagnitude.addItemListener(this);
        checkLogMag.addItemListener(this);
        checkPhase.addItemListener(this);

        checkFFT.addItemListener(this);
        checkIFFT.addItemListener(this);
        checkNoPreproc.addItemListener(this);

        return pp;
    }

    /** Event handler for a change of visualisation or preprocessing state.
        @param e the event to be treated.
    */
    public void itemStateChanged(ItemEvent e)
    {
        Object source = e.getItemSelectable();
        if(source.equals(checkReal)){
            gp.chooseVisualization(GraphingPanel.VIS_REAL);
        } else if(source.equals(checkImag)){
            gp.chooseVisualization(GraphingPanel.VIS_IMAG);
        } else if(source.equals(checkMagnitude)){
            gp.chooseVisualization(GraphingPanel.VIS_MAGNITUDE);
        }  else if(source.equals(checkLogMag)) {
            gp.chooseVisualization(GraphingPanel.VIS_LOGMAG);
        } else if(source.equals(checkPhase)){
            gp.chooseVisualization(GraphingPanel.VIS_PHASE);
        } else if(source.equals(checkFFT)) {
            gp.choosePreprocessing(GraphingPanel.PRE_FFT);
        } else if(source.equals(checkIFFT)) {
            gp.choosePreprocessing(GraphingPanel.PRE_IFFT);
        } else if(source.equals(checkNoPreproc)) {
            gp.choosePreprocessing(GraphingPanel.PRE_NO);
        }
    }

    /** Create a column list for the given number of columns.
        @param ncol the number of columns to be used.
    */
    public void activateColumnList(int ncol)
    {
        Vector<String> v=new Vector<String>();
        v.add("None");
        for(int i=0; i<ncol;++i) {
            v.add("Column "+(i+1));
        }
        postColumnX.setListData(v);
        postColumnY.setListData(v);
        postColumnZ.setListData(v);
        postColumnR.setListData(v);
        postColumnI.setListData(v);
    }

    /** Handle menu actions.
        @param e the event to be treated.
    */
    public void actionPerformed(ActionEvent e)
    {
        String arg=e.getActionCommand();
        if("New".equals(arg)){
            createNewInstance();
        } else if ("Open".equals(arg)){
            String fname=null;

            if(useNativeFileDialogs) {
                // File chooser provided by the host system.
                // Vastly better on MacOSX
                FileDialog fd = new FileDialog(this, "Open");
                fd.setDirectory(currentPath.getAbsolutePath());

                fd.setVisible(true);
                currentPath=new File(fd.getDirectory());
                fname =fd.getFile();
            } else {
                // File chooser provided by Swing.
                // Better on Linux

                JFileChooser fc = new JFileChooser();
                fc.setCurrentDirectory(currentPath);
                fc.setDialogTitle("Open");
                if(fc.showOpenDialog(this)!=JFileChooser.APPROVE_OPTION)
                    return;

                currentPath=new File(
                       fc.getSelectedFile().getParentFile().getPath());

                fname =fc.getSelectedFile().getName();
            }
            updateFileList();
            setFileToBeRead(fname);
            repaint();
        } else if ("About".equals(arg)){
            DialogAbout d=new DialogAbout(this);
            d.setVisible(true);
        } else if ("Quit".equals(arg)){
            System.exit(0);
        } else if ("Close".equals(arg)){
            setVisible(false);
            dispose();
            --openWindowsNumber;
            if (openWindowsNumber<1)
                System.exit(0);
        }  else if ("Grayscale palette".equals(arg)){
            gp.setPalette(gp.createGrayscalePalette());
            gp.forceRedraw();
        } else if ("Jet palette".equals(arg)){
            gp.setPalette(gp.createJetPalette());
            gp.forceRedraw();
        }  else if ("Blue/White/Red palette".equals(arg)){
            gp.setPalette(gp.createBlueRedPalette());
            gp.forceRedraw();
        } else if ("Temperature palette".equals(arg)){
            gp.setPalette(gp.createTempPalette());
            gp.forceRedraw();
        } else if ("Set min/max color range".equals(arg)){
            DialogMinMax m=new DialogMinMax(this);
            m.setMinMaxValues(gp.getCrangeMin(), gp.getCrangeMax());
            m.setMinMaxCB(gp.getCalculateMin(), gp.getCalculateMax());

            m.setVisible(true);
            if(m.getIsOk()) {
                try{
                    gp.setCrangeMin(m.getMinValue());
                    gp.setCrangeMax(m.getMaxValue());
                    gp.setCalculateMax(m.getMaxCB());
                    gp.setCalculateMin(m.getMinCB());
                    gp.forceRedraw();
                } catch (java.lang.NumberFormatException J) {
                    JOptionPane.showMessageDialog(this,
                        "I could not understand the numbers you wrote.",
                        "Woha!", JOptionPane.ERROR_MESSAGE);
                }
            } 
        } else if("Refresh".equals(arg)){
            int tx=columnX;
            int ty=columnY;
            int tz=columnZ;
            int tr=columnR;
            int ti=columnI;
            setFileToBeRead(null);
            // That does not work. Synchronization issue?
            gp.setScatterType(tx,ty,tz,tr,ti);
            repaint();
        }
    }

    /** Create a new instance of the window.
      @return the created instance
     */
    public OptiExFrame createNewInstance()
    {
        OptiExFrame popFrame=new OptiExFrame();
        popFrame.init();

        popFrame.setBounds(getX()+30, getY()+30, popFrame.getWidth(),
                popFrame.getHeight());

        popFrame.setVisible(true);

        return popFrame;
    }

    /** Update the list containing the files and the directories
        in the current path. Note that hidden files might be shown
        or not, depending on the value of the flag showHiddenFiles.

    */
    public void updateFileList()
    {
        File folder = currentPath;

        File[] listOfFiles = folder.listFiles();
        Vector<String>  fileList = new Vector<String>();

        if (listOfFiles ==null){
            listOfFiles = File.listRoots();
            currentPath = listOfFiles[0];
        } else
            fileList.add(DIR_TAG+"..");

        // Some operating systems (MacOSX) already give a sorted
        // list, others do not (Linux). I suppose it depends on
        // the file system used.

        Arrays.sort(listOfFiles);
        for (int i = 0; i < listOfFiles.length; ++i) {
            String name = listOfFiles[i].getName();
            if (listOfFiles[i].isFile()) {
                if(showHiddenFiles || !listOfFiles[i].isHidden())
                    fileList.add(FILE_TAG+name);
            } else if (listOfFiles[i].isDirectory()) {
                if(showHiddenFiles || !name.startsWith("."))
                    fileList.add(DIR_TAG+name);
            }
        }

        fileJList.setListData(fileList);
        fileJList.setVisibleRowCount(-1);
    }

    /** Check the file size and ask to the user if a very big file has to
        be opened.
        @param fn the filename.
    */
    public boolean checkFile(String fn)
    {
        filename=fn;
        File file=new File(filename);
        long length=file.length();

        if(length>MAX_FILE_SIZE) {
            String f =String.format("%3.3g",(double)length/(1024*1024));

            Object[] options = {"Yes",
                "No"};
            int answer=JOptionPane.showOptionDialog(this,
                    "This file seems to be BIG: "+f+" MiB\n"
                    +"Are you sure you want to try to open it?",
                    "Woha!",
                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.QUESTION_MESSAGE,
                    null, options, options[1]);
            if(answer==JOptionPane.NO_OPTION)
                return false;
        }
        return true;
    }

    class SetData implements Runnable
    {
        MatrixStore r;
        int ncol;
        String filename;
        final static int GNUPLOT=1;
        final static int OPTIWAVE=2;
        final static int NONE=0;
        int type;

        /** constructor */
        public SetData()
        {
            r=null;
            type=NONE;
        }
        public void run()
        {
            if(type==GNUPLOT) {
                setDataGnuplot();
            } else if (type==OPTIWAVE){
                setDataOptiwave();
            } else {
                gp.setData(null);
                gp.setFileFormat("none");
            }
            gp.repaint();
        }

        public void selectGnuplot(MatrixStore m, int n, String f)
        {
            type=GNUPLOT;
            r=m;
            ncol=n;
            filename=f;
        }

        public void selectOptiwave(FileInfo fi, String fn)
        {
            f=fi;
            type=OPTIWAVE;
            filename=fn;
        }

        /** Set the program for reading Optiwave .f3d or .rid files.
        */
        public void setDataOptiwave()
        {
            fileNameField.setText(filename);
            gp.setRasterType();
            r=f.mm;
            gp.setXYRanges(f.xmin, f.xmax, f.ymin, f.ymax);
            postProcessing.setVisible(false);
            visualization.setVisible(true);
            preproc.setVisible(true);
            splitPaneTexts.setVisible(false);
            textFile.setText("");
            gp.setData(r);
        }

        /** Set the program for reading Gnuplot (scatter) type files.
        */
        public void setDataGnuplot()
        {
            activateColumnList(ncol);
            fileNameField.setText(filename);
            if(ncol==3)
                updateScatter(0,1,-1,2,-1);
            else if(ncol==4)
                updateScatter(0,1,-1,2,3);
            else if(ncol==2)
                updateScatter(0,-1,-1,1,0);
            else
                updateScatter(0,1,2,3,4);

            try{
                textFile.setText(
                    ReadVariousFormats.read10lines(filename));
            } catch(FileNotFoundException E) {
                System.out.printf("uhm... this is weird...");
            }

            ts.getHorizontalScrollBar().setValue(0);
            ts.getVerticalScrollBar().setValue(0);
            postProcessing.setVisible(true);
            preproc.setVisible(false);
            splitPaneTexts.setVisible(true);
            gp.choosePreprocessing(GraphingPanel.PRE_NO);
            checkNoPreproc.setSelected(true);
            if(columnZ>0) {
                zquotaSlider.setValue(0);
            }
            gp.setData(r);
            gp.setFileFormat(fileFormat);
        }
    }

    /** Specify the file to be read.
        @param fn the file to read. If null, read again the same file
            (refresh).
    */
    public ReadVariousFormats setFileToBeRead(String fn)
    {
        if(fn!=null) 
            filename = currentPath.getAbsolutePath()+File.separator+fn;
        if(!checkFile(filename)) {
            synchronized(gp) {gp.setData(null);}
            return reader;
        }

        File file = new File(filename);
        date = file.lastModified();

        // If necessary, performs some thread synchronization.
        if(reader!=null) {
            // If a reading is already in process, we stop it
            // and we wait until everything has finished.

            while(reader.isReading()){
                reader.stopReading();
                try {
                    Thread.sleep(10);
                } catch (InterruptedException E)
                {}
            }
        } else {
            reader = new ReadVariousFormats();
        }
        if(isDebug)
            System.out.println("--> "+filename);

        SetData clearVisualization = new SetData();
        try {
            if(SwingUtilities.isEventDispatchThread()) {
                SwingUtilities.invokeLater(clearVisualization);
            } else {
                SwingUtilities.invokeAndWait(clearVisualization);
            }
        } catch (InterruptedException E) {
        } catch (InvocationTargetException F) {
        }
        MyTimer mt=new MyTimer();
        int ncol=0;
        try{
            MatrixStore r=null;

            if(ReadVariousFormats.isOptiWave(filename)) {
                FileInfo f= ReadVariousFormats.readOptiWaveFast(filename);
                fileFormat = "Optiwave.";
                SetData doWorkRunnable = new SetData();
                doWorkRunnable.selectOptiwave(f,filename);
                SwingUtilities.invokeLater(doWorkRunnable);
            } else if(ReadVariousFormats.isGnuplot(filename)) {
                ncol = ReadVariousFormats.getGnuplotColumn(filename);
                r= reader.readGnuplot(filename,ncol);
                SetData doWorkRunnable = new SetData();
                doWorkRunnable.selectGnuplot(r,ncol,filename);
                SwingUtilities.invokeLater(doWorkRunnable);
            }
        } catch (java.io.FileNotFoundException E) {
            System.out.println("File not found error.");
        } catch (java.io.IOException F) {
            System.out.println("Something  went wrong.");
        } catch (java.lang.OutOfMemoryError G) {
            System.out.println("Out of mermory error!");
            JOptionPane.showMessageDialog(this,
                    "Not enough memory! I give up.");
            try {
                 SwingUtilities.invokeAndWait(clearVisualization);
             } catch (InterruptedException E) {
             } catch (InvocationTargetException F) {
             }
        }
        if(isDebug)
            System.out.println("Read: elapsed time "+mt.getElapsed()+" ms.");

        return reader;
    }

    public void updateScatter(int X, int Y, int Z, int R, int I)
    {
        // setScatterType raises an IOException if the
        gp.setScatterType(X,Y,Z,R,I);

        columnX=X;
        postColumnX.setSelectedIndex(X+1);
        columnY=Y;
        postColumnY.setSelectedIndex(Y+1);
        columnZ=Z;
        postColumnZ.setSelectedIndex(Z+1);
        columnR=R;
        postColumnR.setSelectedIndex(R+1);
        columnI=I;
        postColumnI.setSelectedIndex(I+1);
    }

    public static void main(String[] args)
    {
        OptiExFrame o = new OptiExFrame();

        o.init();
        o.show();
    }

    private class ListFilesRenderer extends JLabel
        implements ListCellRenderer<Object>
    {
        private Icon emptyIcon;

        public ListFilesRenderer()
        {
            setOpaque(true);
        }

        public Component getListCellRendererComponent(JList<?> list,
                                                   Object value,
                                                   int index,
                                                   boolean isSelected,
                                                   boolean cellHasFocus)
        {
            String ss=value.toString();
            if(ss.startsWith(DIR_TAG)) {
                setIcon(UIManager.getIcon("FileView.directoryIcon"));
                ss=ss.substring(DIR_TAG.length());
            } else {
                ss=ss.substring(FILE_TAG.length());
                boolean ic=false;
                try{
                    String fp=ss;
                    if(currentPath!=null)
                        fp=currentPath.getAbsolutePath()+
                            File.separator+ss;
                    ic=ReadVariousFormats.isGnuplot(fp) ||
                        ReadVariousFormats.isOptiWave(fp);
                } catch(FileNotFoundException E)
                {
                    ic=false;
                }
                if(ic) {
                    setIcon(UIManager.getIcon("FileView.fileIcon"));
                } else {
                    if(emptyIcon==null) {
                        createEmptyIcon();
                    }
                    setIcon(emptyIcon);
                }
            }

            if(isSelected) {
                setBackground(Color.BLUE);
                setForeground(UIManager.getColor("text"));
            } else {
                setBackground(UIManager.getColor("control"));
                setForeground(UIManager.getColor("controlText"));
            }

            setText(ss);
            return this;
        }
        
        private void createEmptyIcon()
        {
            int xsize=UIManager.getIcon("FileView.fileIcon").getIconWidth();
            int ysize=UIManager.getIcon("FileView.fileIcon").getIconHeight();
            BufferedImage bufferedImage = new BufferedImage(xsize,ysize,
                BufferedImage.TYPE_INT_ARGB);
            emptyIcon=new ImageIcon(bufferedImage);
        }
    }
}
