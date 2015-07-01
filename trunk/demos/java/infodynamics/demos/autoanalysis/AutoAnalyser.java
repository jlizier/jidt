/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2015, Joseph T. Lizier
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

import infodynamics.measures.continuous.ConditionalMutualInfoMultiVariateCommon;
import infodynamics.measures.continuous.TransferEntropyCalculator;
import infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian;
import infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel;
import infodynamics.measures.continuous.kraskov.ConditionalMutualInfoCalculatorMultiVariateKraskov;
import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.MatrixUtils;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JComboBox;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JFileChooser;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ToolTipManager;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.MediaTracker;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Vector;

public class AutoAnalyser extends JFrame
	implements ActionListener, DocumentListener, MouseListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	// Set options for the TE Calculator type:
	public final static String CALC_TYPE_DISCRETE = "Discrete";
	public final static String CALC_TYPE_GAUSSIAN = "Gaussian";
	public final static String CALC_TYPE_KRASKOV  = "Kraskov (KSG)";
	public final static String CALC_TYPE_KERNEL   = "Kernel";
	protected String[] calcTypes = {
			CALC_TYPE_DISCRETE, CALC_TYPE_GAUSSIAN,
			CALC_TYPE_KRASKOV, CALC_TYPE_KERNEL};
	protected String[] unitsForEachCalc = {"bits", "nats", "nats", "bits"};
	
	// Store the properties for each calculator
	protected String[] commonContPropertyNames = {
			TransferEntropyCalculator.K_PROP_NAME
	};
	protected String[] commonContPropertiesFieldNames = {
			"K_PROP_NAME"
	};
	protected String[] commonContPropertyDescriptions = {
			"Destination history embedding length (k_HISTORY)"
	};
	protected String[] gaussianProperties = {
			TransferEntropyCalculator.K_TAU_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.L_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.L_TAU_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.DELAY_PROP_NAME, // Not common to Kernel
	};
	protected String[] gaussianPropertiesFieldNames = {
			"TransferEntropyCalculator.K_TAU_PROP_NAME", // Not common to Kernel
			"TransferEntropyCalculator.L_PROP_NAME", // Not common to Kernel
			"TransferEntropyCalculator.L_TAU_PROP_NAME", // Not common to Kernel
			"TransferEntropyCalculator.DELAY_PROP_NAME", // Not common to Kernel
	};
	protected String[] gaussianPropertyDescriptions = {
			"Destination history embedding delay (k_TAU)",
			"Source history embedding length (l_HISTORY)",
			"Source history embeding delay (l_TAU)",
			"Delay from source to destination (in time steps)"
	};
	protected String[] kernelProperties = {
			TransferEntropyCalculatorKernel.KERNEL_WIDTH_PROP_NAME,
			TransferEntropyCalculatorKernel.DYN_CORR_EXCL_TIME_NAME,
			TransferEntropyCalculatorKernel.NORMALISE_PROP_NAME,			
	};
	protected String[] kernelPropertiesFieldNames = {
			"KERNEL_WIDTH_PROP_NAME",
			"DYN_CORR_EXCL_TIME_NAME",
			"NORMALISE_PROP_NAME"			
	};
	protected String[] kernelPropertyDescriptions = {
			"Kernel width to be used in the calculation. <br/>If the property " +
					TransferEntropyCalculatorKernel.NORMALISE_PROP_NAME +
					" is set, then this is a number of standard deviations; " +
					"otherwise it is an absolute value.",
			"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
					"0 (default) means no dynamic exclusion window",
			"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not  (recommended)",
	};
	protected String[] kraskovProperties = {
			TransferEntropyCalculator.K_TAU_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.L_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.L_TAU_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.DELAY_PROP_NAME, // Not common to Kernel
			ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE,
			ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K,
			ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE,
			ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME,
			ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE,
			ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS,
			TransferEntropyCalculatorKraskov.PROP_KRASKOV_ALG_NUM,
			TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD,
			TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX,
			TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX,
			TransferEntropyCalculatorKraskov.PROP_RAGWITZ_NUM_NNS
	};
	protected String[] kraskovPropertiesFieldNames = {
			"TransferEntropyCalculator.K_TAU_PROP_NAME", // Not common to Kernel
			"TransferEntropyCalculator.L_PROP_NAME", // Not common to Kernel
			"TransferEntropyCalculator.L_TAU_PROP_NAME", // Not common to Kernel
			"TransferEntropyCalculator.DELAY_PROP_NAME", // Not common to Kernel
			"ConditionalMutualInfoMultiVariateCommon.PROP_NORMALISE",
			"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_K",
			"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_ADD_NOISE",
			"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_DYN_CORR_EXCL_TIME",
			"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NORM_TYPE",
			"ConditionalMutualInfoCalculatorMultiVariateKraskov.PROP_NUM_THREADS",
			"PROP_KRASKOV_ALG_NUM",
			"PROP_AUTO_EMBED_METHOD",
			"PROP_K_SEARCH_MAX",
			"PROP_TAU_SEARCH_MAX",
			"PROP_RAGWITZ_NUM_NNS"
	};
	protected String[] kraskovPropertyDescriptions = {
			"Destination history embedding delay (k_TAU)",
			"Source history embedding length (l)",
			"Source history embeding delay (l_TAU)",
			"Delay from source to destination (in time steps)",
			"(boolean) whether to normalise <br/>each incoming time-series to mean 0, standard deviation 1, or not (recommended)",
			"Number of k nearest neighbours to use <br/>in the full joint kernel space in the KSG algorithm",
			"Standard deviation for an amount <br/>of random Gaussian noise to add to each variable, " +
					"to avoid having neighbourhoods with artificially large counts. <br/>" +
					"(\"false\" may be used to indicate \"0\".). The amount is added in after any normalisation.",
			"Dynamic correlation exclusion time or <br/>Theiler window (see Kantz and Schreiber); " +
					"0 (default) means no dynamic exclusion window",
			"<br/>Norm type to use in KSG algorithm between the points in each marginal space. <br/>Options are: " +
					"\"MAX_NORM\" (default), otherwise \"EUCLIDEAN\" or \"EUCLIDEAN_SQUARED\" (both equivalent here)",
			"Number of parallel threads to use <br/>in computation: an integer > 0 or \"USE_ALL\" " +
					"(default, to indicate to use all available processors)",
			"Which KSG algorithm to use (1 or 2)",
			"Method to automatically determine embedding lengths (k_HISTORY,l_HISTORY)<br/> and delays (k_TAU, l_TAU) for " +
					"destination and potentially source time-series. Default is \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_NONE +
					"\" meaning values are set manually; other values include: <br/>  -- \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ +
					" for use of the Ragwitz criteria for both source and destination (searching up to \"" + TransferEntropyCalculatorKraskov.PROP_K_SEARCH_MAX +
					"\" and \"" + TransferEntropyCalculatorKraskov.PROP_TAU_SEARCH_MAX + "\"); <br/>  -- \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_RAGWITZ_DEST_ONLY +
					"\" for use of the Ragwitz criteria for the destination only. <br/>Use of values other than \"" + TransferEntropyCalculatorKraskov.AUTO_EMBED_METHOD_NONE +
					"\" leads to any previous settings for embedding lengths and delays for the destination and perhaps source to be overwritten after observations are supplied",
			"Max. embedding length to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD + ")",
			"Max. embedding delay to search to <br/>if auto embedding (as determined by " + TransferEntropyCalculatorKraskov.PROP_AUTO_EMBED_METHOD + ")",
			"Number of k nearest neighbours for <br/>Ragwitz auto embedding (if used; defaults to match property \"k\")"
	};
	
	protected JButton computeButton;
	protected JButton openDataButton;
	// Stores the current data file
	protected File dataFile = null;
	// Combo box for selecting the calculator type
	protected JComboBox<String> calcTypeComboBox;
	// Displays the current data file
	protected JTextField dataFileTextField;
	// Descriptor of the data:
	protected JLabel dataFileDescriptorLabel;
	// The data we're working with
	protected double[][] data = null;
	// Number of rows and columns of the data:
	protected int dataRows = 0;
	protected int dataColumns = 0;
	// Source column field
	protected JTextField sourceColTextField;
	// Destination column field
	protected JTextField destColTextField;
	// Table for the properties
	protected JTable propertiesTable;
	// Table model (local class) for the table
	protected PropertiesTableModel propertiesTableModel;
	// Names of the properties
	protected Vector<String> propertyNames;
	// Names of the fields for the properties
	protected Vector<String> propertyFieldNames;
	// Descriptions of the fields for the properties
	protected Vector<String> propertyDescriptions;
	// Values of the properties
	protected HashMap<String,String> propertyValues;
	// Results text
	protected JLabel resultsLabel;
	// Text area for the generated Java code
	protected JTextArea javaCodeTextArea;
	// Text area for the generated Python code
	protected JTextArea pythonCodeTextArea;
	// Text area for the generated Matlab code
	protected JTextArea matlabCodeTextArea;
	
	protected String codeDefaultText = "    ... Awaiting new parameter selection (press compute) ...";
	
	public class TextAreaWithImage extends JTextArea {

	    /**
		 * Default serialVersionUID
		 */
		private static final long serialVersionUID = 1L;
		/**
		 * Image for the background of the text area
		 */
		protected Image image;
		protected boolean rescaled = false;
		protected int xOffset = 0, yOffset = 0;

	    public TextAreaWithImage(String defaultText, Image image) {
	        super(defaultText);
	        setOpaque(false); // I think this is set by default
	        this.image = image;
	    }

	    @Override
	    protected void paintComponent(Graphics g) {
	    	if (!rescaled) {
	    		int width = getWidth();
	    		int height = getHeight();
	    		if (width < height) {
	    			image = image.getScaledInstance(width, -1, 0);
	    			// TODO The following doesn't work because the image height is
	    			//  not returned.
	    			// yOffset = height - image.getHeight(null)/2;
	    		} else {
	    			image = image.getScaledInstance(-1, height, 0);
	    			// The following doesn't work because the image height is
	    			//  not returned.
	    			// xOffset = width - image.getWidth(this)/2;
	    		}
	    		rescaled = true;
	    	}
	    	// Hacking an offset in because I haven't worked out how to
	    	//  access the current width/height. This link might have a solution:
	    	// http://stackoverflow.com/questions/26386422/how-to-set-background-image-to-a-jtextarea-in-a-jpanel
	        // g.drawImage(image,xOffset,yOffset,null);
	        g.drawImage(image,50,0,null);
	        super.paintComponent(g);
	    }
	}
	
	/**
	 * Constructor to generate the application windows
	 */
	public AutoAnalyser() {
		
		// Build the swing applet
		
		ImageIcon icon = new ImageIcon("../../JIDT-logo.png"); // Location for distributions
		if (icon.getImageLoadStatus() != MediaTracker.COMPLETE) {
			// Try the alternative image location for SVN checkouts
			icon = new ImageIcon("../../web/JIDT-logo.png");
		}
		setIconImage(icon.getImage());
		
		Image watermarkImage = (new ImageIcon("JIDT-logo-watermark.png")).getImage();
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(1100,530);
		setTitle("JIDT Transfer Entropy Auto-Analyser");
		// Centre in the middle of the screen
		setLocationRelativeTo(null);

		// Create the fields for the calc type:
		JLabel calcTypeLabel = new JLabel("TE Calculator Type:");
		calcTypeLabel.setToolTipText("Select estimator type. \"Discrete\" is for discrete or pre-binned data; all others for continuous data.");
		calcTypeLabel.setBorder(BorderFactory.createEmptyBorder(0,0,10,0));
		calcTypeComboBox = (JComboBox<String>) new JComboBox<String>(calcTypes);
		calcTypeComboBox.setSelectedIndex(2); // Select Kraskov by default
		calcTypeComboBox.addActionListener(this);
		// Don't set an empty border as it becomes clickable as well,
		//  we'll use an empty label instead
		//calcTypeComboBox.setBorder(BorderFactory.createEmptyBorder(0,0,10,0));

		// Create the fields for the data file
		JLabel fileLabel = new JLabel("Data file:");
		fileLabel.setToolTipText("Must be a text file of time-series values with time increasing in rows, and variables along columns");
		dataFileTextField = new JTextField(30);
		dataFileTextField.setEnabled(false);
		// TODO set action listener for on click of this field
		dataFileTextField.addMouseListener(this);
		// Don't set border around this text field as it doesn't look right
		// Button to open data file
		openDataButton = new JButton("Select");
		openDataButton.addActionListener(this);
		// Description about the data
		dataFileDescriptorLabel = new JLabel("No data file selected yet ...");
		dataFileDescriptorLabel.setHorizontalAlignment(JLabel.RIGHT);
		dataFileDescriptorLabel.setBorder(BorderFactory.createEmptyBorder(0,0,10,0));
		
		// From column:
		JLabel sourceLabel = new JLabel("Source column:");
		sourceLabel.setToolTipText("(first is 0)");
		sourceColTextField = new JTextField(5);
		sourceColTextField.setEnabled(true);
		sourceColTextField.setText("0");
		// Must add document listener, not add action listener
		sourceColTextField.getDocument().addDocumentListener(this);
		// To column:
		JLabel destLabel = new JLabel("Destination column:");
		destLabel.setToolTipText("(first is 0)");
		destColTextField = new JTextField(5);
		destColTextField.setEnabled(true);
		destColTextField.setText("1");
		destColTextField.getDocument().addDocumentListener(this);

		JLabel dummyLabel1 = new JLabel(" ");
		dummyLabel1.setSize(10, 10);
		JLabel dummyLabel2 = new JLabel(" ");
		dummyLabel2.setSize(10, 10);
		JLabel dummyLabel3 = new JLabel(" ");
		dummyLabel3.setSize(10, 10);
		
		putCalcPropertiesInTable();
		propertiesTableModel = new PropertiesTableModel();
		propertiesTable = new TableWithToolTip(propertiesTableModel);
		// Make sure any properties are saved when the compute button is clicked
		propertiesTable.putClientProperty("terminateEditOnFocusLost", Boolean.TRUE);
		Font headerFont = propertiesTable.getTableHeader().getFont();
		propertiesTable.getTableHeader().setFont(headerFont.deriveFont(Font.BOLD));
		TableColumn valueColumn = propertiesTable.getColumn("Property value");
		valueColumn.setMinWidth(130);
		valueColumn.setMaxWidth(130);
		JScrollPane propsTableScrollPane = new JScrollPane(propertiesTable);
		// Set up for ~18 rows maximum (the +6 is exact to fit all props
		//  for KRaskov in without scrollbar)
		Dimension d = propertiesTable.getPreferredSize();
		propsTableScrollPane.setPreferredSize(
		    new Dimension(d.width,propertiesTable.getRowHeight()*17+6));
		
		
		// Button to compute
		computeButton = new JButton("Compute");
		computeButton.addActionListener(this);
		
		// Results label
		resultsLabel = new JLabel(" ");

		// Code panel
		JPanel codePanel = new JPanel();
		codePanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Generated code"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
		JTabbedPane codeTabbedPane = new JTabbedPane();
		// Generate Java code text area
		javaCodeTextArea = new TextAreaWithImage(codeDefaultText, watermarkImage);
		javaCodeTextArea.setOpaque(false);
		javaCodeTextArea.setEditable(false);
		javaCodeTextArea.setBorder(BorderFactory.createCompoundBorder(
				javaCodeTextArea.getBorder(), 
	            BorderFactory.createEmptyBorder(5,5,5,5)));
		//codeTextArea.setLineWrap(true); // This makes it all unreadable
		//codeTextArea.setWrapStyleWord(true);
		JScrollPane javaAreaScrollPane = new JScrollPane(javaCodeTextArea);
        javaAreaScrollPane.setVerticalScrollBarPolicy(
        		JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        javaAreaScrollPane.setHorizontalScrollBarPolicy(
        		JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        int codeTextAreaWidth = 560;
        int codeTextAreaHeight = 460;
        Dimension codeTextAreaDimension = 
        		new Dimension(codeTextAreaWidth, codeTextAreaHeight);
        javaAreaScrollPane.setPreferredSize(codeTextAreaDimension);
        javaAreaScrollPane.setMinimumSize(codeTextAreaDimension);
        javaAreaScrollPane.setMaximumSize(codeTextAreaDimension);
        codeTabbedPane.addTab("Java", javaAreaScrollPane);
		// Generate Python code text area
		pythonCodeTextArea = new TextAreaWithImage(codeDefaultText, watermarkImage);
		pythonCodeTextArea.setOpaque(false);
		pythonCodeTextArea.setEditable(false);
		pythonCodeTextArea.setBorder(BorderFactory.createCompoundBorder(
				pythonCodeTextArea.getBorder(), 
	            BorderFactory.createEmptyBorder(5,5,5,5)));
		JScrollPane pythonAreaScrollPane = new JScrollPane(pythonCodeTextArea);
		pythonAreaScrollPane.setVerticalScrollBarPolicy(
        		JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		pythonAreaScrollPane.setHorizontalScrollBarPolicy(
        		JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		pythonAreaScrollPane.setPreferredSize(codeTextAreaDimension);
		pythonAreaScrollPane.setMinimumSize(codeTextAreaDimension);
		pythonAreaScrollPane.setMaximumSize(codeTextAreaDimension);
        codeTabbedPane.addTab("Python", pythonAreaScrollPane);
		// Generate Matlab code text area
		matlabCodeTextArea = new TextAreaWithImage(codeDefaultText, watermarkImage);
		matlabCodeTextArea.setOpaque(false);
		matlabCodeTextArea.setEditable(false);
		matlabCodeTextArea.setBorder(BorderFactory.createCompoundBorder(
				matlabCodeTextArea.getBorder(), 
	            BorderFactory.createEmptyBorder(5,5,5,5)));
		JScrollPane matlabAreaScrollPane = new JScrollPane(matlabCodeTextArea);
		matlabAreaScrollPane.setVerticalScrollBarPolicy(
        		JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		matlabAreaScrollPane.setHorizontalScrollBarPolicy(
        		JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		matlabAreaScrollPane.setPreferredSize(codeTextAreaDimension);
		matlabAreaScrollPane.setMinimumSize(codeTextAreaDimension);
		matlabAreaScrollPane.setMaximumSize(codeTextAreaDimension);
        codeTabbedPane.addTab("Matlab", matlabAreaScrollPane);
        // Now add the tabbed pane to the panel
		codePanel.add(codeTabbedPane);
		codePanel.setSize(codeTextAreaWidth+10, codeTextAreaHeight+10);
		
		// Add all the components in:
		/*
		add(calcTypePanel, BorderLayout.NORTH);
		add(dataFileChooserPanel, BorderLayout.EAST);
		add(dataFileDescriptorPanel, BorderLayout.WEST);
		add(computeButton, BorderLayout.SOUTH);
		*/
		// gridbag.setConstraints(calcTypePanel, c);
		// add(calcTypePanel);
		JPanel paramsPanel = new JPanel();
		paramsPanel.setBorder(BorderFactory.createCompoundBorder(
                                BorderFactory.createTitledBorder("Calculation parameters"),
                                BorderFactory.createEmptyBorder(5,5,5,5)));
		GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
        paramsPanel.setLayout(gridbag);
        c.anchor = GridBagConstraints.EAST; // Not sure what I put EAST for?
        
        // Add the CalcType label and combobox
        c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
        c.fill = GridBagConstraints.NONE;      //reset to default
        c.weightx = 0.0;                       //reset to default
        paramsPanel.add(calcTypeLabel, c);
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        paramsPanel.add(calcTypeComboBox, c);
        // Add dummy label for spacing
        paramsPanel.add(dummyLabel3, c);

        
        // Add the data file chooser fields
        c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
        c.fill = GridBagConstraints.NONE;      //reset to default
        c.weightx = 0.0;                       //reset to default
        paramsPanel.add(fileLabel, c);
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        paramsPanel.add(dataFileTextField, c);
        c.gridx = 1;
        paramsPanel.add(openDataButton, c);
        c.gridx = -1; // Reset to no indication
        
        paramsPanel.add(dataFileDescriptorLabel, c);

        // Add the source selector
        c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
        c.fill = GridBagConstraints.NONE;      //reset to default
        c.weightx = 0.0;                       //reset to default
        paramsPanel.add(sourceLabel, c);
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        paramsPanel.add(sourceColTextField, c);
        // Add the destination selector
        c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
        c.fill = GridBagConstraints.NONE;      //reset to default
        c.weightx = 0.0;                       //reset to default
        paramsPanel.add(destLabel, c);
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        paramsPanel.add(destColTextField, c);

        // Add dummy label for spacing
        paramsPanel.add(dummyLabel1, c);
        // Add the properties table
        paramsPanel.add(propsTableScrollPane, c);
        // Add dummy label for spacing
        paramsPanel.add(dummyLabel2, c);
        // Add the compute button
        paramsPanel.add(computeButton, c);
        // Add the results text
        paramsPanel.add(resultsLabel, c);
        	
		// Add both panels into the frame with Border layout
		add(paramsPanel, BorderLayout.WEST);
		add(codePanel);
		
		setVisible(true);
		
		// The default tool tip delay before dismissing was too short to read these, so 
		// I'm setting it to 30 sec.
		ToolTipManager.sharedInstance().setDismissDelay(30000);
		
		// Try to force the watermark images to come up
		javaCodeTextArea.repaint();
		pythonCodeTextArea.repaint();
		pythonCodeTextArea.repaint();
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		// Clear text fields for now if anything changes
		resultsLabel.setText(" ");
		javaCodeTextArea.setText(codeDefaultText);
		pythonCodeTextArea.setText(codeDefaultText);
		matlabCodeTextArea.setText(codeDefaultText);
		
		if (e.getSource() == computeButton) {
			computeTE();
		} else if (e.getSource() == openDataButton) {
			selectFileAction();
		} else if (e.getSource() == calcTypeComboBox) {
			putCalcPropertiesInTable();
			propertiesTableModel.fireTableDataChanged(); // Alerts to refresh the table contents
			System.out.println("Added properties for new calculator");
		}
		// Else nothing extra to do
	}

	public void repaint() {
		System.out.println("repainting ...");
		super.repaint();
	}
	
	protected void selectFileAction() {
		System.out.println("Open data button pressed ...");
		// Give user a choice of file, starting from the currently selected
		//  one
		JFileChooser dataFileChooser;
		if (dataFile == null) {
			dataFileChooser = new JFileChooser(System.getProperty("user.dir") + "/../data/");
		} else {
			dataFileChooser = new JFileChooser(dataFile);
		}
		int returnVal = dataFileChooser.showOpenDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			// File selection was approved:
			dataFile = dataFileChooser.getSelectedFile();
			dataFileTextField.setText(dataFile.getAbsolutePath());
			System.out.println("Data file selected: " + dataFile.getAbsolutePath());
			// Now load the file in:
			ArrayFileReader afr = new ArrayFileReader(dataFile);
			try {
				data = afr.getDouble2DMatrix();
				dataRows = data.length;
				if (dataRows > 0) {
					dataColumns = data[0].length;
				} else {
					dataColumns = 0;
				}
				dataFileDescriptorLabel.setText(
						String.format("Valid data file with %d rows and %d columns",
						dataRows, dataColumns));
				System.out.printf("Read in data with %d rows and %d columns\n",
						dataRows, dataColumns);
			} catch (Exception ex) {
				ex.printStackTrace(System.err);
				JOptionPane.showMessageDialog(this,
						ex.getMessage());
				dataFileDescriptorLabel.setText("Invalid data file, please load another");
				data = null;
			}
		}
		// Else do nothing
	}
	
	protected void computeTE() {
		
		System.out.println("Compute button pressed ...");
		resultsLabel.setText("Computing ...");
		if (data == null) {
			JOptionPane.showMessageDialog(this, "No valid data source selected");
			return;
		}
		int sourceColumn = Integer.parseInt(sourceColTextField.getText());
		int destColumn = Integer.parseInt(destColTextField.getText());
		if ((sourceColumn < 0) || (sourceColumn >= dataColumns)) {
			JOptionPane.showMessageDialog(this,
					String.format("Source column must be between 0 and %d for this data set",
							dataColumns-1));
			return;
		}
		if ((destColumn < 0) || (destColumn >= dataColumns)) {
			JOptionPane.showMessageDialog(this,
					String.format("Destination column must be between 0 and %d for this data set",
							dataColumns-1));
			return;
		}

		String selectedCalcType = (String)
				calcTypeComboBox.getSelectedItem();
		String units = unitsForEachCalc[calcTypeComboBox.getSelectedIndex()];
		
		// Generate headers:
		// 1. Java
		StringBuffer javaCode = new StringBuffer();
		javaCode.append("package infodynamics.demos.autoanalysis;\n\n");
		javaCode.append("import infodynamics.utils.ArrayFileReader;\n");
		javaCode.append("import infodynamics.utils.MatrixUtils;\n\n");
		// Cover TransferEntropyCalculator and any common conditional MI classes
		//  used for property names
		javaCode.append("import infodynamics.measures.continuous.*;\n");
		// 2. Python
		StringBuffer pythonCode = new StringBuffer();
		pythonCode.append("from jpype import *\n");
		pythonCode.append("import numpy\n");
		pythonCode.append("# I think this is a bit of a hack, python users will do better on this:\n");
		pythonCode.append("sys.path.append(\"../python\")\n");
		pythonCode.append("import readFloatsFile\n\n");
		pythonCode.append("# Add JIDT jar library to the path\n\n");
		pythonCode.append("jarLocation = \"../../infodynamics.jar\"\n");
		pythonCode.append("# Start the JVM (add the \"-Xmx\" option with say 1024M if you get crashes due to not enough memory space)\n");
		pythonCode.append("startJVM(getDefaultJVMPath(), \"-ea\", \"-Djava.class.path=\" + jarLocation)\n\n");
		// 3. Matlab:
		StringBuffer matlabCode = new StringBuffer();
		matlabCode.append("% Add JIDT jar library to the path\n");
		matlabCode.append("javaaddpath('../../infodynamics.jar');\n");
		matlabCode.append("% Add utilities to the path\n");
		matlabCode.append("addpath('../octave');\n\n");
		
		try{
			// Create a Kraskov TE calculator:
			TransferEntropyCalculator teCalc;
			
			// Construct an instance of the selected calculator:
			String javaConstructorLine = null;
			String pythonPreConstructorLine = null;
			String matlabConstructorLine = null;
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				throw new Exception("Not implemented yet");
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
				teCalc = new TransferEntropyCalculatorGaussian();
				// Cover the TE calculator and any references to conditional MI calculator properties
				javaCode.append("import infodynamics.measures.continuous.gaussian.*;\n");
				javaConstructorLine = "    teCalc = new TransferEntropyCalculatorGaussian();\n";
				pythonPreConstructorLine = "teCalcClass = JPackage(\"infodynamics.measures.continuous.gaussian\").TransferEntropyCalculatorGaussian\n";
				matlabConstructorLine = "teCalc = javaObject('infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian');\n";
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
				teCalc = new TransferEntropyCalculatorKraskov();
				// Cover the TE calculator and any references to conditional MI calculator properties
				javaCode.append("import infodynamics.measures.continuous.kraskov.*;\n");
				javaConstructorLine = "    teCalc = new TransferEntropyCalculatorKraskov();\n";
				pythonPreConstructorLine = "teCalcClass = JPackage(\"infodynamics.measures.continuous.kraskov\").TransferEntropyCalculatorKraskov\n";
				matlabConstructorLine = "teCalc = javaObject('infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov');\n";
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
				teCalc = new TransferEntropyCalculatorKernel();
				// Cover the TE calculator and any references to conditional MI calculator properties
				javaCode.append("import infodynamics.measures.continuous.kernel.*;\n");
				javaConstructorLine = "    teCalc = new TransferEntropyCalculatorKernel();\n";
				pythonPreConstructorLine = "teCalcClass = JPackage(\"infodynamics.measures.continuous.kernel\").TransferEntropyCalculatorKernel\n";
				matlabConstructorLine = "teCalc = javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');\n";
			} else {
				throw new Exception("No recognised calculator selected: " +
						selectedCalcType);
			}
			
			javaCode.append("\npublic class GeneratedTECalculator {\n\n");
			javaCode.append("  public static void main(String[] args) throws Exception {\n\n");
			
			// Code to read in data:
			String loadDataComment = "0. Load/prepare the data:\n";
			// 1. Java
			javaCode.append("    // " + loadDataComment);
			javaCode.append("    String dataFile = \"" + dataFile.getAbsolutePath() + "\";\n");
			javaCode.append("    ArrayFileReader afr = new ArrayFileReader(dataFile);\n");
			javaCode.append("    double[][] data = afr.getDouble2DMatrix();\n");
			javaCode.append("    double[] source = MatrixUtils.selectColumn(data, " + sourceColumn + ");\n");
			javaCode.append("    double[] dest = MatrixUtils.selectColumn(data, " + destColumn + ");\n\n");
			// 2. Python
			pythonCode.append("# " + loadDataComment);
			pythonCode.append("dataRaw = readFloatsFile.readFloatsFile(\"" + dataFile.getAbsolutePath() + "\")\n");
			pythonCode.append("# As numpy array:\n");
			pythonCode.append("data = numpy.array(dataRaw)\n");
			pythonCode.append("source = data[:," + sourceColumn + "]\n");
			pythonCode.append("dest = data[:," + destColumn + "]\n\n");
			// 2. Matlab
			matlabCode.append("% " + loadDataComment);
			matlabCode.append("data = load('" + dataFile.getAbsolutePath() + "');\n");
			matlabCode.append("% Column indices start from 1 in Matlab:\n");
			matlabCode.append("source = octaveToJavaDoubleArray(data(:," + (sourceColumn+1) + "));\n");
			matlabCode.append("dest = octaveToJavaDoubleArray(data(:," + (destColumn+1) + "));\n\n");

			// Construct the calculator:
			String constructComment = "1. Construct the calculator:\n";
			// 1. Java
			javaCode.append("    // " + constructComment);
			javaCode.append("    TransferEntropyCalculator teCalc;\n");
			javaCode.append(javaConstructorLine);
			// 2. Python
			pythonCode.append("# " + constructComment);
			pythonCode.append(pythonPreConstructorLine);
			pythonCode.append("teCalc = teCalcClass()\n");
			// 3. Matlab
			matlabCode.append("% " + constructComment);
			matlabCode.append(matlabConstructorLine);
			
			// Set properties 
			String setPropertiesComment = "2. Set any properties to non-default values:\n";
			javaCode.append("    // " + setPropertiesComment);
			pythonCode.append("# " + setPropertiesComment);
			matlabCode.append("% " + setPropertiesComment);
			int i = 0;
			for (String propName : propertyNames) {
				String propValue = null;
				String propFieldName = propertyFieldNames.get(i++);
				try {
					propValue = propertyValues.get(propName);
				} catch (Exception ex) {
					ex.printStackTrace(System.err);
					JOptionPane.showMessageDialog(this,
							ex.getMessage());
					resultsLabel.setText("Cannot find a value for property " + propName);
				}
				// Check whether this property value is different to the default for
				//  this calculator. This is more for generating the minimal code.
				if (!propValue.equalsIgnoreCase(teCalc.getProperty(propName))) {
					// We need to set this property:
					teCalc.setProperty(propName, propValue);
					// 1. Java Code -- use full field name here
					javaCode.append("    teCalc.setProperty(" + propFieldName +
								",\n        \"" +
								propValue + "\");\n");
					// 2. Python code
					pythonCode.append("teCalc.setProperty(\"" + propName + "\", \"" +
							propValue + "\")\n");
					// 3. Matlab code
					matlabCode.append("teCalc.setProperty('" + propName + "', '" +
							propValue + "');\n");
				}
			}
			
			// Initialise
			teCalc.initialise();
			String initialiseComment = "3. Initialise the calculator for (re-)use:\n";
			javaCode.append("    // " + initialiseComment);
			javaCode.append("    teCalc.initialise();\n");
			pythonCode.append("# " + initialiseComment);
			pythonCode.append("teCalc.initialise()\n");
			matlabCode.append("% " + initialiseComment);
			matlabCode.append("teCalc.initialise();\n");
			
			// Set observations
			teCalc.setObservations(
					MatrixUtils.selectColumn(data, sourceColumn),
					MatrixUtils.selectColumn(data, destColumn));
			String supplyDataComment = "4. Supply the sample data:\n";
			// 1. Java
			javaCode.append("    // " + supplyDataComment);
			javaCode.append("    teCalc.setObservations(source, dest);\n");
			// 2. Python
			pythonCode.append("# " + supplyDataComment);
			pythonCode.append("teCalc.setObservations(source, dest)\n");
			// 3. Matlab
			matlabCode.append("% " + supplyDataComment);
			matlabCode.append("teCalc.setObservations(source, dest);\n");
			
			// Compute TE 
			double teResult = teCalc.computeAverageLocalOfObservations();
			String computeComment = "5. Compute the estimate:\n";
			javaCode.append("    // " + computeComment);
			javaCode.append("    double teResult = teCalc.computeAverageLocalOfObservations();\n");
			pythonCode.append("# " + computeComment);
			pythonCode.append("teResult = teCalc.computeAverageLocalOfObservations()\n");
			matlabCode.append("% " + computeComment);
			matlabCode.append("teResult = teCalc.computeAverageLocalOfObservations();\n");
			String resultsPrefixString = String.format("TE_%s(col_%d -> col_%d) = ",
					selectedCalcType, sourceColumn, destColumn);
			resultsLabel.setText(String.format(resultsPrefixString + "%.4f %s", teResult, units));
			// And generate code to write the results and finalise:
			// 1. Java
			javaCode.append("    System.out.printf(\"" + resultsPrefixString + "%.4f " + units + "\\n\", teResult);\n");
			javaCode.append("  }\n");
			javaCode.append("}\n\n");
			// 2. Python
			pythonCode.append("print(\"" + resultsPrefixString + "%.4f " + units + "\\n\" % teResult)\n");
			// 3. Matlab
			matlabCode.append("fprintf('" + resultsPrefixString + "%.4f " + units + "\\n', teResult);\n");
			
			// Now set the code in the panel for the user
			javaCodeTextArea.setText(javaCode.toString());
			javaCodeTextArea.setCaretPosition(0); // Pull focus to the top
			pythonCodeTextArea.setText(pythonCode.toString());
			pythonCodeTextArea.setCaretPosition(0); // Pull focus to the top
			matlabCodeTextArea.setText(matlabCode.toString());
			matlabCodeTextArea.setCaretPosition(0); // Pull focus to the top
			
			// Now write the code to a file
			// 1. Java
			FileWriter codeFileWriter = new FileWriter("../java/infodynamics/demos/autoanalysis/GeneratedTECalculator.java");
			codeFileWriter.write(javaCode.toString());
			codeFileWriter.close();
			// 2. Python
			codeFileWriter = new FileWriter("GeneratedTECalculator.py");
			codeFileWriter.write(pythonCode.toString());
			codeFileWriter.close();
			// 3. Matlab
			codeFileWriter = new FileWriter("GeneratedTECalculator.m");
			codeFileWriter.write(matlabCode.toString());
			codeFileWriter.close();
			
			// Read the current property values back out (in case of 
			//  automated parameter assignment)
			for (String propName : propertyNames) {
				String propValue = null;
				try {
					propValue = teCalc.getProperty(propName);
					propertyValues.put(propName, propValue);
				} catch (Exception ex) {
					ex.printStackTrace(System.err);
					JOptionPane.showMessageDialog(this,
							ex.getMessage());
					resultsLabel.setText("Cannot find a value for property " + propName);
				}
			}
			propertiesTableModel.fireTableDataChanged(); // Alerts to refresh the table contents
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Calculation failed, please see console");
		}

	}
	
	protected class TableWithToolTip extends JTable {
		
		/**
		 * Default serialVersionUID
		 */
		private static final long serialVersionUID = 1L;

		public TableWithToolTip(AbstractTableModel atm) {
			super(atm);
		}
		
		public Component prepareRenderer(TableCellRenderer renderer,
				int rowIndex, int vColIndex) {
			Component c = super.prepareRenderer(renderer, rowIndex, vColIndex);
			if (c instanceof JComponent) {
				if (vColIndex == 0) {
					JComponent jc = (JComponent)c;
					try {
						jc.setToolTipText("<html>" + propertyFieldNames.get(rowIndex) + ": " + propertyDescriptions.get(rowIndex) + "</html>");
					} catch (ArrayIndexOutOfBoundsException aioobe) {
						// Catch if the row number was outside our array of descriptions (e.g. empty row)
					}
				}
			}
			return c;
		}
		/* Alternative method (this over-rides the above if both are in place)
		public String getToolTipText(MouseEvent event) {
			
			Point point = event.getPoint();
			int rowIndex = rowAtPoint(point);
			int columnIndex = columnAtPoint(point);
			
			if (columnIndex == 0) {
				try {
					return propertyFieldNames.get(rowIndex) + ": " + propertyDescriptions.get(rowIndex);
				} catch (ArrayIndexOutOfBoundsException aioobe) {
					// Catch if the row number was outside our array of descriptions (e.g. empty row)
					return null;
				}
			} else {
				// No tool tip for the property values
				return null;
			}
        }
        */
	}
	
	protected class PropertiesTableModel extends AbstractTableModel {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		
		public PropertiesTableModel() {
			super();
		}

		@Override
		public String getColumnName(int column) {
			if (column == 0) {
				return "Property name";
			} else {
				return "Property value";
			}
		}

		@Override
		public boolean isCellEditable(int rowIndex, int columnIndex) {
			if (columnIndex == 0) {
				return false;
			}
			return true;
		}

		@Override
		public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
			if (columnIndex == 0) {
				return;
			}
			// Else we've changed a property value -- which one?
			String propName = propertyNames.get(rowIndex);
			propertyValues.put(propName, (String) aValue);
			
			// And clear the result and code panels because of this change:
			resultsLabel.setText(" "); // Clear text for now if anything changes
			javaCodeTextArea.setText(codeDefaultText);
			pythonCodeTextArea.setText(codeDefaultText);
			matlabCodeTextArea.setText(codeDefaultText);
		}

		@Override
		public int getColumnCount() {
			return 2;
		}

		@Override
		public int getRowCount() {
			return propertyNames.size();
		}

		@Override
		public Object getValueAt(int rowIndex, int columnIndex) {
			String propName = propertyNames.get(rowIndex);
			if (columnIndex == 0) {
				return propName;
			} else {
				return propertyValues.get(propName);
			}
		}
		
	}
	
	public void putCalcPropertiesInTable() {
		System.out.println("Getting calc properties");
		String selectedCalcType = (String)
				calcTypeComboBox.getSelectedItem();
		@SuppressWarnings("rawtypes")
		Class teCalcClass = null;
		String[] classSpecificPropertyNames = null;
		String[] classSpecificPropertiesFieldNames = null;
		String[] classSpecificPropertyDescriptions = null;
		TransferEntropyCalculator teCalc = null;
		try {
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				throw new Exception("Not implemented yet");
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
				teCalcClass = TransferEntropyCalculatorGaussian.class;
				teCalc = new TransferEntropyCalculatorGaussian();
				classSpecificPropertyNames = gaussianProperties;
				classSpecificPropertiesFieldNames = gaussianPropertiesFieldNames;
				classSpecificPropertyDescriptions = gaussianPropertyDescriptions;
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
				teCalcClass = TransferEntropyCalculatorKraskov.class;
				teCalc = new TransferEntropyCalculatorKraskov();
				classSpecificPropertyNames = kraskovProperties;
				classSpecificPropertiesFieldNames = kraskovPropertiesFieldNames;
				classSpecificPropertyDescriptions = kraskovPropertyDescriptions;
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
				teCalcClass = TransferEntropyCalculatorKernel.class;
				teCalc = new TransferEntropyCalculatorKernel();
				classSpecificPropertyNames = kernelProperties;
				classSpecificPropertiesFieldNames = kernelPropertiesFieldNames;
				classSpecificPropertyDescriptions = kernelPropertyDescriptions;
			} else {
				throw new Exception("No recognised calculator selected: " +
						selectedCalcType);
			}
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Cannot load requested calculator");
		}
		
		// Now get all of the possible properties for this class:
		propertyNames = new Vector<String>();
		propertyFieldNames = new Vector<String>();
		propertyDescriptions = new Vector<String>();
		
		// First for the common properties
		int i = 0;
		for (String fieldName : commonContPropertiesFieldNames) {
			String propName = commonContPropertyNames[i];
			String propertyDescription = commonContPropertyDescriptions[i];
			i++;
			System.out.println("Adding property name TransferEntropyCalculator." + fieldName + " = \"" + propName + "\"");
			propertyFieldNames.add("TransferEntropyCalculator." + fieldName);
			propertyNames.add(propName);
			propertyDescriptions.add(propertyDescription);
		}
		// Then for the specific estimator types
		i = 0;
		for (String fieldName : classSpecificPropertiesFieldNames) {
			String propName = classSpecificPropertyNames[i];
			String propertyDescription = classSpecificPropertyDescriptions[i];
			i++;
			propertyNames.add(propName);
			propertyDescriptions.add(propertyDescription);
			if (fieldName.contains(".")) {
				System.out.println("Adding property name " + fieldName +
						" = \"" + propName + "\"");
				propertyFieldNames.add(fieldName);
			} else {
				System.out.println("Adding property name " + teCalcClass.getSimpleName() +
						"." + fieldName + " = \"" + propName + "\"");
				propertyFieldNames.add(teCalcClass.getSimpleName() + "." + fieldName);
			}
		}
		
		// Now extract the default values for all of these properties:
		propertyValues = new HashMap<String,String>();
		for (String propName : propertyNames) {
			String defaultPropValue = null;
			try {
				defaultPropValue = teCalc.getProperty(propName);
			} catch (Exception ex) {
				ex.printStackTrace(System.err);
				JOptionPane.showMessageDialog(this,
						ex.getMessage());
				propertyValues.put(propName, "Cannot find a value");
			}
			propertyValues.put(propName, defaultPropValue);
		}
	}

	@Override
	public void changedUpdate(DocumentEvent e) {
		// Source or dest col number updated
		resultsLabel.setText(" "); // Clear text for now if anything changes
		javaCodeTextArea.setText(codeDefaultText);
		pythonCodeTextArea.setText(codeDefaultText);
		matlabCodeTextArea.setText(codeDefaultText);
	}

	@Override
	public void insertUpdate(DocumentEvent e) {
		// Source or dest col number updated
		resultsLabel.setText(" "); // Clear text for now if anything changes
		javaCodeTextArea.setText(codeDefaultText);
		pythonCodeTextArea.setText(codeDefaultText);
		matlabCodeTextArea.setText(codeDefaultText);
	}

	@Override
	public void removeUpdate(DocumentEvent e) {
		// Source or dest col number updated
		resultsLabel.setText(" "); // Clear text for now if anything changes
		javaCodeTextArea.setText(codeDefaultText);
		pythonCodeTextArea.setText(codeDefaultText);
		matlabCodeTextArea.setText(codeDefaultText);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyser();
	}

	@Override
	public void mouseClicked(MouseEvent me) {
		// User clicked on the data file JLabel
		
		// Clear text fields for now if anything changes
		resultsLabel.setText(" ");
		javaCodeTextArea.setText(codeDefaultText);
		pythonCodeTextArea.setText(codeDefaultText);
		matlabCodeTextArea.setText(codeDefaultText);
		if (me.getSource() == dataFileTextField) {
			selectFileAction();
		}
	}

	@Override
	public void mouseEntered(MouseEvent me) {
		// Do nothing
	}

	@Override
	public void mouseExited(MouseEvent me) {
		// Do nothing
	}

	@Override
	public void mousePressed(MouseEvent me) {
		// Do nothing
	}

	@Override
	public void mouseReleased(MouseEvent me) {
		// Do nothing
	}
}
