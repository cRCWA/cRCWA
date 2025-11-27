import javax.swing.*;
import java.util.*;

import dialogs.*;

import com.apple.eawt.*;


/** The class AppleSpecific implements a few mechanism for interacting with
    the MacOSX operating system. This class will only be used if the program
    detects it is being run on a MacOSX operating system.
    This can be a problem when the program is not compiled on a MacOSX
    operating system, since the com.apple.eawt package is made available
    only under this platform. You should thus remove each reference to the
    AppleSpecific class in the code when compiling under an alternative
    system.

    Copyright 2009-2014 by Davide Bucci
*/
class AppleSpecific implements ApplicationListener
{

    /** Create an application listener able to respond to a few Finder events
    */
    public void answerFinder() {
        Application app = new Application();
        app.setEnabledPreferencesMenu(true);
        app.getApplication().addApplicationListener(this);
    }

    /** Respond to an user clicking on an About menu.
    */
    public void handleAbout(ApplicationEvent evt)
    {
        DialogAbout d=new DialogAbout(null);
        d.setVisible(true);
        evt.setHandled(true);

    }
    /** Respond to an user opening the application.

    */
    public void handleOpenApplication(ApplicationEvent evt)
    {
    }
    /** Respond to an user double clicking on a FCD file

    */
    public void handleOpenFile(ApplicationEvent evt)
    {

    }

    /** Respond to an user clicking on the Preferences menu.

    */
    public void handlePreferences(ApplicationEvent evt)
    {
    }

    /** Respond to an user wanting to print a particular file.
    */
    public void handlePrintFile(ApplicationEvent evt)
    {
        // nothing to do
    }

    /** Ask for confirmation when quitting.

    */
    public void handleQuit(ApplicationEvent evt)
    {
        evt.setHandled(true);
    }

    public void handleReOpenApplication(ApplicationEvent evt)
    {
        // Nothing to do
    }
}
