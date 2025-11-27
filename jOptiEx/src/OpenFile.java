class OpenFile implements Runnable {
    private final OptiExFrame parent;
    private final String file;


    public OpenFile(OptiExFrame o, String f)
    {
        file=f;
        parent=o;
    }

    public void run()
    {
        parent.setFileToBeRead(file);
        parent.gp.repaint();
    }
}
