/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2017, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.demos.autoanalysis;


import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import javax.swing.ToolTipManager;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.MediaTracker;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.lang.reflect.Constructor;
import java.net.URL;

/**
 * This class provides a single GUI to launch the AutoAnalyser GUIs from.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public class AutoAnalyserLauncher extends JFrame
	implements ActionListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	protected JButton[] launcherButtons;
	
	protected String[] buttonLabels = {
			"Entropy",
			"Mutual Info",
			"Conditional Mutual Info",
			"Active Info Storage",
			"Transfer Entropy",
			"Conditional Transfer Entropy"
	};

	@SuppressWarnings("rawtypes")
	protected Class[] launcherClasses = {
			AutoAnalyserEntropy.class,
			AutoAnalyserMI.class,
			AutoAnalyserCMI.class,
			AutoAnalyserAIS.class,
			AutoAnalyserTE.class,
			AutoAnalyserCTE.class
	};
	
	protected String appletTitle = "JIDT AutoAnalyser Launcher";
	
	protected String offsetFromJidtToAutoAnalyserFolder = "demos/AutoAnalyser/";
	
	protected String autoAnalyserFolder = null;
	
	/**
	 * Constructor to generate the application windows
	 */
	public AutoAnalyserLauncher(boolean useCurrentLocation) {
		
		String jidtFolder = null;
		if (!useCurrentLocation) {
			// Figure out the location of the resource this class was run from.
			// Presumably it is the infodynamics.jar file:
			URL resource = getClass().getProtectionDomain().getCodeSource().getLocation();
			if (resource.getPath().endsWith("jar")) {
				// This resource came from a jar
				jidtFolder = resource.getPath().replaceFirst("infodynamics.jar", "");
				autoAnalyserFolder = jidtFolder + offsetFromJidtToAutoAnalyserFolder;
			}
			// System.out.println(resource.getPath());
		}
		if (autoAnalyserFolder == null) {
			// Either useCurrentLocation is true, or the current resource did not
			//  come from the infodynamics.jar file. Either way, we're simply
			//  going to assume the current location is the auto analyser folder
			autoAnalyserFolder = System.getProperty("user.dir") + "/";
			jidtFolder = autoAnalyserFolder + "/../../";
		}

		// Build the swing applet
		
		ImageIcon icon = new ImageIcon(jidtFolder + "JIDT-logo.png"); // Location for distributions
		if (icon.getImageLoadStatus() != MediaTracker.COMPLETE) {
			// Try the alternative image location for git checkouts
			icon = new ImageIcon(jidtFolder + "web/JIDT-logo.png");
		}
		setIconImage(icon.getImage());
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(250,300);
		setTitle(appletTitle);
		// Centre in the middle of the screen
		setLocationRelativeTo(null);

		// Add buttons to launch each AutoAnalyser:
		launcherButtons = new JButton[buttonLabels.length];
		for (int i = 0; i < launcherButtons.length; i++) {
			launcherButtons[i] = new JButton(buttonLabels[i]);
			launcherButtons[i].addActionListener(this);
		}
				
		// Add all the components in:
		/*
		add(calcTypePanel, BorderLayout.NORTH);
		add(dataFileChooserPanel, BorderLayout.EAST);
		add(dataFileDescriptorPanel, BorderLayout.WEST);
		add(computeButton, BorderLayout.SOUTH);
		*/
		// gridbag.setConstraints(calcTypePanel, c);
		// add(calcTypePanel);
		JPanel calcButtonsPanel = new JPanel();
		calcButtonsPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
		GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
        calcButtonsPanel.setLayout(gridbag);
        c.anchor = GridBagConstraints.CENTER; // Not sure what I put EAST for?
        
        // Add the buttons for each AutoAnalyser
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.fill = GridBagConstraints.BOTH;
        c.weightx = 1.0;
		JLabel textLabel = new JLabel("Select AutoAnalyser to launch:");
		textLabel.setSize(10, 10);
		textLabel.setHorizontalAlignment(SwingConstants.CENTER);
        calcButtonsPanel.add(textLabel, c);
		JLabel dummyLabel = new JLabel(" ");
		dummyLabel.setSize(10, 10);
		dummyLabel.setHorizontalAlignment(SwingConstants.CENTER);
        calcButtonsPanel.add(dummyLabel, c);
        for (int i = 0; i < launcherButtons.length; i++) {
        	launcherButtons[i].setHorizontalAlignment(SwingConstants.CENTER);
            calcButtonsPanel.add(launcherButtons[i], c);
    		JLabel dummyLabel1 = new JLabel(" ");
    		dummyLabel1.setSize(10, 10);
            calcButtonsPanel.add(dummyLabel1, c);
        }
        
		// Add all panels into the frame with Border layout
		add(calcButtonsPanel, BorderLayout.WEST);
		
		setVisible(true);
		
		// The default tool tip delay before dismissing was too short to read these, so 
		// I'm setting it to 30 sec.
		ToolTipManager.sharedInstance().setDismissDelay(30000);
		
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		// Check which (if any) button the action came from
		for (int i = 0; i < launcherButtons.length; i++) {
			if (e.getSource() == launcherButtons[i]) {
				try {
					// Call the constructor for the relevant AutoAnalyser
					//  passing in the offset to the AutoAnalyser directory.
					@SuppressWarnings({ "rawtypes", "unchecked" })
					Constructor cons = launcherClasses[i].getConstructor(String.class);
					cons.newInstance(autoAnalyserFolder);
				} catch (Exception ex) {
					JOptionPane.showMessageDialog(null, ex.getMessage(),
							"Error launching " + buttonLabels[i] + " AutoAnalyser",
							JOptionPane.ERROR_MESSAGE);
					ex.printStackTrace();
					System.out.println();
					return;
				}
				// The AutoAnalyser started correctly, so we can close this launcher
				dispose();
			}
		}
		// Else nothing extra to do
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Assume that any command argument means to take the current
		//  directory as the AutoAnalyser folder.
		//  Double clicking jar won't pass any arguments in.
		new AutoAnalyserLauncher(args.length > 0);
	}
}
