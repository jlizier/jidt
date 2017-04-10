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

import infodynamics.measures.continuous.ChannelCalculatorCommon;
import infodynamics.measures.discrete.ChannelCalculatorDiscrete;
import infodynamics.utils.AnalyticMeasurementDistribution;
import infodynamics.utils.AnalyticNullDistributionComputer;
import infodynamics.utils.ArrayFileReader;
import infodynamics.utils.EmpiricalMeasurementDistribution;
import infodynamics.utils.MatrixUtils;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JComboBox;
import javax.swing.JCheckBox;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JFileChooser;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ToolTipManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
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

/**
 * This abstract class provides a GUI to build a simple calculation,
 *  and supply the code to execute it.
 * Child classes fix this to a TE or MI calculation.
 * 
 * 
 * @author Joseph Lizier
 *
 */
public abstract class AutoAnalyser extends JFrame
	implements ActionListener, DocumentListener, MouseListener, ChangeListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	// Set options for the Calculator type:
	public final static String CALC_TYPE_DISCRETE = "Discrete";
	public final static String CALC_TYPE_GAUSSIAN = "Gaussian";
	public final static String CALC_TYPE_KRASKOV  = "Kraskov (KSG)";
	public final static String CALC_TYPE_KERNEL   = "Kernel";
	protected String[] calcTypes = {
			CALC_TYPE_DISCRETE, CALC_TYPE_GAUSSIAN,
			CALC_TYPE_KRASKOV, CALC_TYPE_KERNEL};
	protected String[] unitsForEachCalc = {"bits", "nats", "nats", "bits"};
	
	// Properties for each calculator
	protected static final String DISCRETE_PROPNAME_BASE = "base";
	protected String[] discreteProperties; // Children to initialise
	protected String[] discretePropertyDefaultValues; // Children to initialise
	protected String[] discretePropertyDescriptions; // Children to initialise
	
	// Common property names for all continuous calculators:
	protected String[] commonContPropertyNames;
	protected String[] commonContPropertiesFieldNames;
	protected String[] commonContPropertyDescriptions;
	protected String[] gaussianProperties;
	protected String[] gaussianPropertiesFieldNames;
	protected String[] gaussianPropertyDescriptions;
	protected String[] kernelProperties;
	protected String[] kernelPropertiesFieldNames;
	protected String[] kernelPropertyDescriptions;
	protected String[] kraskovProperties;
	protected String[] kraskovPropertiesFieldNames;
	protected String[] kraskovPropertyDescriptions;

	// Store which calculator type we're using:
	@SuppressWarnings("rawtypes")
	protected Class calcClass = null;
	
	@SuppressWarnings("rawtypes")
	protected Class abstractContinuousClass = null;
	@SuppressWarnings("rawtypes")
	protected Class discreteClass = null;

	protected String measureAcronym = null;
	
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
	// The continuous data we're working with
	protected double[][] data = null;
	// The discrete data we're working with
	protected int[][] dataDiscrete = null;
	// Number of rows and columns of the data:
	protected int dataRows = 0;
	protected int dataColumns = 0;
	// Source column field
	protected JTextField sourceColTextField;
	// Destination column field
	protected JTextField destColTextField;
	// CheckBox for "all pairs?"
	protected JCheckBox allPairsCheckBox;
	// CheckBox for statistical significance
	protected JCheckBox statSigCheckBox;
	// CheckBox for statistical significance analytically
	protected JCheckBox statSigAnalyticCheckBox;
	// Number of permutations to use for statistical significance
	protected int numPermutationsToCheck = 100;
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
	// CheckBox for "Compute result?"
	protected JCheckBox computeResultCheckBox;
	// Results text
	protected JLabel resultsLabel;
	// Text area for the generated Java code
	protected JTextArea javaCodeTextArea;
	// Text area for the generated Python code
	protected JTextArea pythonCodeTextArea;
	// Text area for the generated Matlab code
	protected JTextArea matlabCodeTextArea;
	
	protected String codeDefaultText = "    ... Awaiting new parameter selection (press compute) ...";
	
	protected String buttonTextCodeAndCompute = "Generate code and Compute";
	protected String buttonTextCodeOnly = "Generate code";
	
	protected String appletTitle;
	
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
		
		makeSpecificInitialisations();
		
		// Build the swing applet
		
		ImageIcon icon = new ImageIcon("../../JIDT-logo.png"); // Location for distributions
		if (icon.getImageLoadStatus() != MediaTracker.COMPLETE) {
			// Try the alternative image location for SVN checkouts
			icon = new ImageIcon("../../web/JIDT-logo.png");
		}
		setIconImage(icon.getImage());
		
		Image watermarkImage = (new ImageIcon("JIDT-logo-watermark.png")).getImage();
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(1100,640);
		setTitle(appletTitle);
		// Centre in the middle of the screen
		setLocationRelativeTo(null);

		// Create the fields for the calc type:
		JLabel calcTypeLabel = new JLabel("Calculator Type:");
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
		sourceLabel.setToolTipText("First column is 0.");
		sourceColTextField = new JTextField(5);
		sourceColTextField.setEnabled(true);
		sourceColTextField.setText("0");
		// Must add document listener, not add action listener
		sourceColTextField.getDocument().addDocumentListener(this);
		// To column:
		JLabel destLabel = new JLabel("Destination column:");
		destLabel.setToolTipText("First column is 0.");
		destLabel.setBorder(BorderFactory.createEmptyBorder(0,0,10,0));
		destColTextField = new JTextField(5);
		destColTextField.setEnabled(true);
		destColTextField.setText("1");
		destColTextField.getDocument().addDocumentListener(this);

		// Checkbox for all pairs
		allPairsCheckBox = new JCheckBox("All pairs?");
		allPairsCheckBox.setToolTipText("Compute for all pairs of columns?");
		allPairsCheckBox.addChangeListener(this);
		
		// CheckBox for statistical signficance
		statSigCheckBox = new JCheckBox("Add stat. signif.?");
		statSigCheckBox.setToolTipText("Compute the null surrogate distribution under the null hypothesis that source and target are not related?");
		statSigCheckBox.addChangeListener(this);
		statSigAnalyticCheckBox = new JCheckBox("analytically?");
		statSigAnalyticCheckBox.setToolTipText("Compute the null surrogate distribution analytically?");
		statSigAnalyticCheckBox.addChangeListener(this);
		statSigAnalyticCheckBox.setEnabled(false);
		
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
		//  for Kraskov TE in without scrollbar)
		Dimension d = propertiesTable.getPreferredSize();
		int rowHeight = propertiesTable.getRowHeight();
		propsTableScrollPane.setPreferredSize(
		    new Dimension(d.width,rowHeight*17+6));
		propsTableScrollPane.setMinimumSize(
			new Dimension(d.width,rowHeight*17+6));
		System.out.println("Row height was " + rowHeight);
		
		// Checkbox for compute result
		computeResultCheckBox = new JCheckBox("Compute result?");
		computeResultCheckBox.setToolTipText("Compute result or only generate code?");
		computeResultCheckBox.setSelected(true);
		computeResultCheckBox.addChangeListener(this);

		// Button to compute
		computeButton = new JButton(buttonTextCodeAndCompute);
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
        int codeTextAreaHeight = 480;
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

        // Add the all pairs label and checkbox
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.gridx = 1;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        paramsPanel.add(allPairsCheckBox, c);
        c.gridx = -1; // reset to no indication
        
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

        // Add the CheckBoxes for statistical significance checks:
        c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
        c.fill = GridBagConstraints.NONE;      //reset to default
        c.weightx = 0.0;                       //reset to default
        paramsPanel.add(statSigCheckBox, c);
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        paramsPanel.add(statSigAnalyticCheckBox, c);

        // Add dummy label for spacing
        paramsPanel.add(dummyLabel1, c);
        // Add the properties table
        paramsPanel.add(propsTableScrollPane, c);
        // Add dummy label for spacing
        paramsPanel.add(dummyLabel2, c);
        // Add the "compute result" checkbox:
        c.gridwidth = GridBagConstraints.RELATIVE; //next-to-last
        c.fill = GridBagConstraints.NONE;      //reset to default
        c.weightx = 0.0;                       //reset to default
        paramsPanel.add(computeResultCheckBox);
        // Add the compute button
        c.gridwidth = GridBagConstraints.REMAINDER;     //end row
        c.fill = GridBagConstraints.HORIZONTAL;
        c.weightx = 1.0;
        paramsPanel.add(computeButton, c);
        
        // Status panel
		JPanel statusPanel = new JPanel();
		statusPanel.setBorder(BorderFactory.createCompoundBorder(
                                BorderFactory.createTitledBorder("Status"),
                                BorderFactory.createEmptyBorder(5,5,5,5)));
        statusPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
        // Add the results text
		statusPanel.add(resultsLabel);
        
		// Add all panels into the frame with Border layout
		add(paramsPanel, BorderLayout.WEST);
		add(codePanel, BorderLayout.EAST);
		add(statusPanel, BorderLayout.SOUTH);
		
		setVisible(true);
		
		// The default tool tip delay before dismissing was too short to read these, so 
		// I'm setting it to 30 sec.
		ToolTipManager.sharedInstance().setDismissDelay(30000);
		
		// Try to force the watermark images to come up
		javaCodeTextArea.repaint();
		pythonCodeTextArea.repaint();
		pythonCodeTextArea.repaint();
	}

	/**
	 * Child classes over-ride this to set strings etc specifically for their measure type
	 */
	protected abstract void makeSpecificInitialisations();
	
	@Override
	public void actionPerformed(ActionEvent e) {
		// Clear text fields for now if anything changes
		resultsLabel.setText(" ");
		javaCodeTextArea.setText(codeDefaultText);
		pythonCodeTextArea.setText(codeDefaultText);
		matlabCodeTextArea.setText(codeDefaultText);
		
		if (e.getSource() == computeButton) {
			compute();
		} else if (e.getSource() == openDataButton) {
			selectFileAction();
		} else if (e.getSource() == calcTypeComboBox) {
			putCalcPropertiesInTable();
			propertiesTableModel.fireTableDataChanged(); // Alerts to refresh the table contents
			System.out.println("Added properties for new calculator");
			if (statSigCheckBox.isSelected()) {
				// We might allow the analytic box to be checked if the selected
				//  calculator allows this:
				if (calcImplementsAnalyticStatSig()) {
					statSigAnalyticCheckBox.setEnabled(true);
				} else {
					statSigAnalyticCheckBox.setSelected(false);
					statSigAnalyticCheckBox.setEnabled(false);
				}
			} else {
				statSigAnalyticCheckBox.setSelected(false);
				statSigAnalyticCheckBox.setEnabled(false);
			}

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
			// Now load the file in to check it:
			loadData(false); // Assume is doubles for the moment
		}
		// Else do nothing
	}
	
	protected void loadData(boolean isInts) {
		ArrayFileReader afr = new ArrayFileReader(dataFile);
		try {
			if (isInts) {
				dataDiscrete = afr.getInt2DMatrix();
				dataRows = dataDiscrete.length;
				if (dataRows > 0) {
					dataColumns = dataDiscrete[0].length;
				} else {
					dataColumns = 0;
				}
			} else {
				data = afr.getDouble2DMatrix();
				dataRows = data.length;
				if (dataRows > 0) {
					dataColumns = data[0].length;
				} else {
					dataColumns = 0;
				}
			}
			dataFileDescriptorLabel.setText(
					String.format("Valid data file with %d rows and %d columns",
					dataRows, dataColumns));
			System.out.printf("Read in data with %d rows and %d columns\n",
					dataRows, dataColumns);
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
			JOptionPane.showMessageDialog(this, ex.getMessage() +
					(isInts ? ".\nFor Discrete calculator make sure you load a file with integer data!" : ""));
			dataFileDescriptorLabel.setText("Invalid data file, please load another");
			if (isInts) {
				dataDiscrete = null;
			} else {
				data = null;
			}
		}
	}
	
	protected void compute() {
		
		System.out.println("Compute button pressed ...");
		resultsLabel.setText("Computing ...");
		
		String selectedCalcType = (String)
				calcTypeComboBox.getSelectedItem();
		String units = unitsForEachCalc[calcTypeComboBox.getSelectedIndex()];
		
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
			// Try to read the data file as ints now:
			loadData(true);
			if (dataDiscrete == null) {
				// An error message will have been shown by loadData() 
				resultsLabel.setText(" ");
				return;
			}
		} else if (data == null) {
			// We need a file of continuous data but no file has been selected
			JOptionPane.showMessageDialog(this, "No valid data source selected");
			resultsLabel.setText(" ");
			return;
		}
		
		int singleSourceColumn = -1, singleDestColumn = -1;
		Vector<int[]> sourceDestPairs = new Vector<int[]>();
		if (allPairsCheckBox.isSelected()) {
			// We're doing all pairs
			for (int s = 0; s < dataColumns; s++) {
				for (int d = 0; d < dataColumns; d++) {
					sourceDestPairs.add(new int[] {s, d});
				}
			}
		} else {
			// we're doing a single source-destination pair
			singleSourceColumn = Integer.parseInt(sourceColTextField.getText());
			singleDestColumn = Integer.parseInt(destColTextField.getText());
			sourceDestPairs.add(new int[] {singleSourceColumn, singleDestColumn});
			if ((singleSourceColumn < 0) || (singleSourceColumn >= dataColumns)) {
				JOptionPane.showMessageDialog(this,
						String.format("Source column must be between 0 and %d for this data set",
								dataColumns-1));
				resultsLabel.setText(" ");				
				return;
			}
			if ((singleDestColumn < 0) || (singleDestColumn >= dataColumns)) {
				JOptionPane.showMessageDialog(this,
						String.format("Destination column must be between 0 and %d for this data set",
								dataColumns-1));
				resultsLabel.setText(" ");
				return;
			}
		}

		// Generate headers:
		// 1. Java
		StringBuffer javaCode = new StringBuffer();
		javaCode.append("package infodynamics.demos.autoanalysis;\n\n");
		javaCode.append("import infodynamics.utils.ArrayFileReader;\n");
		if (statSigCheckBox.isSelected()) {
			javaCode.append("import infodynamics.utils.EmpiricalMeasurementDistribution;\n");
		}
		javaCode.append("import infodynamics.utils.MatrixUtils;\n\n");
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
			// Cover all children:
			javaCode.append("import infodynamics.measures.discrete.*;\n");
		} else {
			// Cover the calculator and any common conditional MI classes
			//  used for property names
			javaCode.append("import infodynamics.measures.continuous.*;\n");
		}
		// 2. Python
		StringBuffer pythonCode = new StringBuffer();
		pythonCode.append("from jpype import *\n");
		pythonCode.append("import numpy\n");
		pythonCode.append("# I think this is a bit of a hack, python users will do better on this:\n");
		pythonCode.append("sys.path.append(\"../python\")\n");
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
			pythonCode.append("import readIntsFile\n\n");
		} else {
			pythonCode.append("import readFloatsFile\n\n");
		}
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
			// Create both discrete and continuous calculators to make
			//  following code simpler:
			ChannelCalculatorCommon calcContinuous = null;
			ChannelCalculatorDiscrete calcDiscrete = null;
			
			// Construct an instance of the selected calculator:
			String javaConstructorLine = null;
			String pythonPreConstructorLine = null;
			String matlabConstructorLine = null;
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				// Defer our processing for this to below ...
			} else {
				calcContinuous = assignCalcObjectContinuous(selectedCalcType);
				if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
					// Cover the calculator and any references to conditional MI calculator properties
					javaCode.append("import infodynamics.measures.continuous.gaussian.*;\n");
					javaConstructorLine = "    calc = new " + calcContinuous.getClass().getSimpleName() + "();\n";
					pythonPreConstructorLine = "calcClass = JPackage(\"infodynamics.measures.continuous.gaussian\")." +
							calcContinuous.getClass().getSimpleName() + "\n";
					matlabConstructorLine = "calc = javaObject('infodynamics.measures.continuous.gaussian." +
							calcContinuous.getClass().getSimpleName() + "');\n";
				} else if (selectedCalcType.startsWith(CALC_TYPE_KRASKOV)) {
					// The if statement will work for both MI Kraskov calculators
					// Cover the calculator and any references to conditional MI calculator properties
					javaCode.append("import infodynamics.measures.continuous.kraskov.*;\n");
					javaConstructorLine = "    calc = new " + calcContinuous.getClass().getSimpleName() + "();\n";
					pythonPreConstructorLine = "calcClass = JPackage(\"infodynamics.measures.continuous.kraskov\")." +
							calcContinuous.getClass().getSimpleName() + "\n";
					matlabConstructorLine = "calc = javaObject('infodynamics.measures.continuous.kraskov." +
							calcContinuous.getClass().getSimpleName() + "');\n";
				} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
					// Cover the calculator and any references to conditional MI calculator properties
					javaCode.append("import infodynamics.measures.continuous.kernel.*;\n");
					javaConstructorLine = "    calc = new " + calcContinuous.getClass().getSimpleName() + "();\n";
					pythonPreConstructorLine = "calcClass = JPackage(\"infodynamics.measures.continuous.kernel\")." +
							calcContinuous.getClass().getSimpleName() + "\n";
					matlabConstructorLine = "calc = javaObject('infodynamics.measures.continuous.kernel." +
							calcContinuous.getClass().getSimpleName() + "');\n";
				} else {
					throw new Exception("No recognised calculator selected: " +
							selectedCalcType);
				}
			}
			
			javaCode.append("\npublic class GeneratedCalculator {\n\n");
			javaCode.append("  public static void main(String[] args) throws Exception {\n\n");
			
			// Code to read in data:
			String loadDataComment = "0. Load/prepare the data:\n";
			// 1. Java
			javaCode.append("    // " + loadDataComment);
			javaCode.append("    String dataFile = \"" + dataFile.getAbsolutePath() + "\";\n");
			javaCode.append("    ArrayFileReader afr = new ArrayFileReader(dataFile);\n");
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				javaCode.append("    int[][] data = afr.getInt2DMatrix();\n");
				if (! allPairsCheckBox.isSelected()) {
					javaCode.append("    int[] source = MatrixUtils.selectColumn(data, " + singleSourceColumn + ");\n");
					javaCode.append("    int[] dest = MatrixUtils.selectColumn(data, " + singleDestColumn + ");\n\n");
				}
			} else {
				javaCode.append("    double[][] data = afr.getDouble2DMatrix();\n");
				if (! allPairsCheckBox.isSelected()) {
					javaCode.append("    double[] source = MatrixUtils.selectColumn(data, " + singleSourceColumn + ");\n");
					javaCode.append("    double[] dest = MatrixUtils.selectColumn(data, " + singleDestColumn + ");\n\n");
				}
			}
			// 2. Python
			pythonCode.append("# " + loadDataComment);
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				pythonCode.append("dataRaw = readIntsFile.readIntsFile(\"" + dataFile.getAbsolutePath() + "\")\n");
			} else {
				pythonCode.append("dataRaw = readFloatsFile.readFloatsFile(\"" + dataFile.getAbsolutePath() + "\")\n");
			}
			pythonCode.append("# As numpy array:\n");
			pythonCode.append("data = numpy.array(dataRaw)\n");
			if (! allPairsCheckBox.isSelected()) {
				pythonCode.append("source = JArray(JInt, 1)(data[:," + singleSourceColumn + "].tolist())\n");
				pythonCode.append("dest = JArray(JInt, 1)(data[:," + singleDestColumn + "].tolist())\n\n");
			}
			// 3. Matlab
			matlabCode.append("% " + loadDataComment);
			matlabCode.append("data = load('" + dataFile.getAbsolutePath() + "');\n");
			if (! allPairsCheckBox.isSelected()) {
				matlabCode.append("% Column indices start from 1 in Matlab:\n");
				if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
					matlabCode.append("source = octaveToJavaIntArray(data(:," + (singleSourceColumn+1) + "));\n");
					matlabCode.append("dest = octaveToJavaIntArray(data(:," + (singleDestColumn+1) + "));\n\n");
				} else {
					matlabCode.append("source = octaveToJavaDoubleArray(data(:," + (singleSourceColumn+1) + "));\n");
					matlabCode.append("dest = octaveToJavaDoubleArray(data(:," + (singleDestColumn+1) + "));\n\n");
				}
			}

			// Construct the calculator and set relevant properties:
			String constructComment = "1. Construct the calculator:\n";
			javaCode.append("    // " + constructComment);
			pythonCode.append("# " + constructComment);
			matlabCode.append("% " + constructComment);
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				DiscreteCalcAndArguments dcaa = assignCalcObjectDiscrete();
				if (dcaa == null) {
					return;
				}
				calcDiscrete = dcaa.calc;
				int base = dcaa.base;
				String args = dcaa.arguments;

				// Now check that the data is ok:
				int minInData = MatrixUtils.min(dataDiscrete);
				int maxInData = MatrixUtils.max(dataDiscrete);
				if ((minInData < 0) || (maxInData >= base)) {
					JOptionPane.showMessageDialog(this,
							"Values in data file (in range " + minInData +
							":" + maxInData + ") lie outside the expected range 0:" +
							(base-1) + " for base " + base);
					resultsLabel.setText(" ");
					return;
				}
				
				// 1. Java
				javaCode.append("    " + calcDiscrete.getClass().getSimpleName() + " calc\n");
				javaCode.append("        = new " + calcDiscrete.getClass().getSimpleName() +
								"(" + args + ");\n");
				// 2. Python
				pythonCode.append("calcClass = JPackage(\"infodynamics.measures.discrete\")." +
						calcDiscrete.getClass().getSimpleName() + "\n");
				pythonCode.append("calc = calcClass(" + args + ")\n");
				// 3. Matlab
				matlabCode.append("calc = javaObject('infodynamics.measures.discrete." +
						calcDiscrete.getClass().getSimpleName() + "', " +
						args + ");\n");
				String setPropertiesComment = "2. No other properties to set for discrete calculators.\n";
				javaCode.append("    // " + setPropertiesComment);
				pythonCode.append("# " + setPropertiesComment);
				matlabCode.append("% " + setPropertiesComment);
			} else {
				// Construct the calculator:
				// 1. Java
				javaCode.append("    " + abstractContinuousClass.getSimpleName() + " calc;\n");
				javaCode.append(javaConstructorLine);
				// 2. Python
				pythonCode.append(pythonPreConstructorLine);
				pythonCode.append("calc = calcClass()\n");
				// 3. Matlab
				matlabCode.append(matlabConstructorLine);
			
				// Set properties 
				String setPropertiesComment = "2. Set any properties to non-default values:\n";
				javaCode.append("    // " + setPropertiesComment);
				pythonCode.append("# " + setPropertiesComment);
				matlabCode.append("% " + setPropertiesComment);
				int i = 0;
				boolean anyPropValueChanged = false;
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
					if (!propValue.equalsIgnoreCase(calcContinuous.getProperty(propName))) {
						// We need to set this property:
						anyPropValueChanged = true;
						calcContinuous.setProperty(propName, propValue);
						// 1. Java Code -- use full field name here
						javaCode.append("    calc.setProperty(" + propFieldName +
									",\n        \"" +
									propValue + "\");\n");
						// 2. Python code
						pythonCode.append("calc.setProperty(\"" + propName + "\", \"" +
								propValue + "\")\n");
						// 3. Matlab code
						matlabCode.append("calc.setProperty('" + propName + "', '" +
								propValue + "');\n");
					}
				}
				if (!anyPropValueChanged) {
					String noPropValueChangedComment = "No properties were set to non-default values\n";
					javaCode.append("    // " + noPropValueChangedComment);
					pythonCode.append("# " + noPropValueChangedComment);
					matlabCode.append("% " + noPropValueChangedComment);
				}
			}
			
			String javaPrefix = "    ";
			String pythonPrefix = "";
			String matlabPrefix = "";
			if (sourceDestPairs.size() > 1) {
				// We're doing all pairs
				
				// Set up loops in the code:
				// 1. Java code
				javaCode.append("    \n");
				javaCode.append("    // Compute for all pairs:\n");
				javaCode.append("    for (int s = 0; s < " + dataColumns +
							"; s++) {\n");
				javaCode.append("        for (int d = 0; d < " + dataColumns +
						"; d++) {\n");
				javaPrefix = "            ";
				javaCode.append(javaPrefix + "// For each source-dest pair:\n");
				javaCode.append(javaPrefix + "if (s == d) {\n");
				javaCode.append(javaPrefix + "    continue;\n");
				javaCode.append(javaPrefix + "}\n");
				// 2. Python code
				pythonCode.append("\n");
				pythonCode.append("# Compute for all pairs:\n");
				pythonCode.append("for s in range(" + dataColumns + "):\n");
				pythonCode.append("\tfor d in range(" + dataColumns + "):\n");
				pythonPrefix = "\t\t";
				pythonCode.append(pythonPrefix+ "# For each source-dest pair:\n");
				pythonCode.append(pythonPrefix + "if (s == d):\n");
				pythonCode.append(pythonPrefix + "\tcontinue\n");
				// 3. Matlab code
				matlabCode.append("\n");
				matlabCode.append("% Compute for all pairs:\n");
				matlabCode.append("for s = 1:" + dataColumns + "\n");
				matlabCode.append("\tfor d = 1:" + dataColumns + "\n");
				matlabPrefix = "\t\t";
				matlabCode.append(matlabPrefix + "% For each source-dest pair:\n");
				matlabCode.append(matlabPrefix + "if (s == d)\n");
				matlabCode.append(matlabPrefix + "\tcontinue;\n");
				matlabCode.append(matlabPrefix + "end\n");

				// Read in the data for these columns for all pairs
				matlabCode.append(matlabPrefix + "% Column indices start from 1 in Matlab:\n");
				if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
					javaCode.append(javaPrefix + "int[] source = MatrixUtils.selectColumn(data, s);\n");
					javaCode.append(javaPrefix + "int[] dest = MatrixUtils.selectColumn(data, d);\n\n");
					matlabCode.append(matlabPrefix + "source = octaveToJavaIntArray(data(:, s));\n");
					matlabCode.append(matlabPrefix + "dest = octaveToJavaIntArray(data(:, d));\n\n");
					pythonCode.append(pythonPrefix + "source = JArray(JInt, 1)(data[:, s].tolist())\n");
					pythonCode.append(pythonPrefix + "dest = JArray(JInt, 1)(data[:, d].tolist())\n\n");
				} else {
					javaCode.append(javaPrefix + "double[] source = MatrixUtils.selectColumn(data, s);\n");
					javaCode.append(javaPrefix + "double[] dest = MatrixUtils.selectColumn(data, d);\n\n");
					matlabCode.append(matlabPrefix + "source = octaveToJavaDoubleArray(data(:, s));\n");
					matlabCode.append(matlabPrefix + "dest = octaveToJavaDoubleArray(data(:, d));\n\n");
					pythonCode.append(pythonPrefix + "source = data[:, s]\n");
					pythonCode.append(pythonPrefix + "dest = data[:, d]\n\n");
				}
			} else {
				// For a single pair, we don't need to set up loops, and
				// we've already read in the data to source and dest variables above
			}
			
			// Initialise
			String initialiseComment = "3. Initialise the calculator for (re-)use:\n";
			javaCode.append(javaPrefix + "// " + initialiseComment);
			javaCode.append(javaPrefix + "calc.initialise();\n");
			pythonCode.append(pythonPrefix + "# " + initialiseComment);
			pythonCode.append(pythonPrefix + "calc.initialise()\n");
			matlabCode.append(matlabPrefix + "% " + initialiseComment);
			matlabCode.append(matlabPrefix + "calc.initialise();\n");
				
			// Set observations
			String supplyDataComment = "4. Supply the sample data:\n";
			String setObservationsMethod;
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				setObservationsMethod = "addObservations";
			} else {
				setObservationsMethod = "setObservations";
			}
			// 1. Java
			javaCode.append(javaPrefix + "// " + supplyDataComment);
			javaCode.append(javaPrefix + "calc." + setObservationsMethod + "(source, dest);\n");
			// 2. Python
			pythonCode.append(pythonPrefix + "# " + supplyDataComment);
			pythonCode.append(pythonPrefix + "calc." + setObservationsMethod + "(source, dest)\n");
			// 3. Matlab
			matlabCode.append(matlabPrefix + "% " + supplyDataComment);
			matlabCode.append(matlabPrefix + "calc." + setObservationsMethod + "(source, dest);\n");
			
			// Compute the result
			String computeComment = "5. Compute the estimate:\n";
			javaCode.append(javaPrefix + "// " + computeComment);
			javaCode.append(javaPrefix + "double result = calc.computeAverageLocalOfObservations();\n");
			pythonCode.append(pythonPrefix + "# " + computeComment);
			pythonCode.append(pythonPrefix + "result = calc.computeAverageLocalOfObservations()\n");
			matlabCode.append(matlabPrefix + "% " + computeComment);
			matlabCode.append(matlabPrefix + "result = calc.computeAverageLocalOfObservations();\n");
			
			String resultsPrefixString, extraFormatTerms;
			if (sourceDestPairs.size() > 1) {
				resultsPrefixString = String.format(measureAcronym + "_%s(col_%%d -> col_%%d) = ",
						selectedCalcType);
				extraFormatTerms = "s, d, ";
			} else {
				resultsPrefixString = String.format(measureAcronym + "_%s(col_%d -> col_%d) = ",
					selectedCalcType, singleSourceColumn, singleDestColumn);
				extraFormatTerms = "";
			}
			
			String resultsSuffixString = "";
			String statSigFormatTerms = "";
			if (statSigCheckBox.isSelected()) {
				// Compute statistical significance
				String statSigComment = "6. Compute the (statistical significance via) null distribution (e.g. 100 permutations):\n";
				javaCode.append(javaPrefix + "// " + statSigComment);
				javaCode.append(javaPrefix + "EmpiricalMeasurementDistribution measDist = calc.computeSignificance(100);\n");
				pythonCode.append(pythonPrefix + "# " + statSigComment);
				pythonCode.append(pythonPrefix + "measDist = calc.computeSignificance(100)\n");
				matlabCode.append(matlabPrefix + "% " + statSigComment);
				matlabCode.append(matlabPrefix + "measDist = calc.computeSignificance(100);\n");
				if (statSigAnalyticCheckBox.isSelected()) {
					resultsSuffixString = " (analytic p(surrogate > measured)=%.5f)";
					statSigFormatTerms = ", measDist.pValue";
				} else {
					resultsSuffixString = " (null: %.4f +/- %.4f std dev.; p(surrogate > measured)=%.5f from %d surrogates)";
					statSigFormatTerms = String.format(", measDist.getMeanOfDistribution(), measDist.getStdOfDistribution(), measDist.pValue, %d",
							numPermutationsToCheck);
				}
			}
			
			// And generate code to write the results and finalise:
			// 1. Java
			javaCode.append("\n" + javaPrefix + "System.out.printf(\"" + resultsPrefixString +
					"%.4f " + units + resultsSuffixString + "\\n\",\n    " + javaPrefix +
					extraFormatTerms + "result" + statSigFormatTerms + ");\n");
			// 2. Python
			pythonCode.append("\n" + pythonPrefix + "print(\"" + resultsPrefixString +
					"%.4f " + units + resultsSuffixString + "\" %\n\t" + pythonPrefix + "(" +
					extraFormatTerms + "result" + statSigFormatTerms + "))\n");
			// 3. Matlab
			matlabCode.append("\n" + matlabPrefix + "fprintf('" + resultsPrefixString +
					"%.4f " + units + resultsSuffixString + "\\n', ...\n\t" + matlabPrefix +
					extraFormatTerms + "result" + statSigFormatTerms + ");\n");
			
			if (sourceDestPairs.size() > 1) {
				// We're doing all pairs
				
				// Finalise the loops in the code:
				// 1. Java code
				javaCode.append("        }\n");
				javaCode.append("    }\n");
				// 2. Python code
				// Nothing to do
				// 3. Matlab code
				matlabCode.append("\tend\n");
				matlabCode.append("end\n");
				
				// And tell the user to see results in the console:
				resultsLabel.setText("See console for all pairs calculation results");
				if (!selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
					System.out.println("Property values not read back into GUI for all pairs calculation");
				}
			}
			
			// Finalise the Java code:
			javaCode.append("  }\n");
			javaCode.append("}\n\n");
			
			// Now set the code in the panel for the user
			javaCodeTextArea.setText(javaCode.toString());
			javaCodeTextArea.setCaretPosition(0); // Pull focus to the top
			pythonCodeTextArea.setText(pythonCode.toString());
			pythonCodeTextArea.setCaretPosition(0); // Pull focus to the top
			matlabCodeTextArea.setText(matlabCode.toString());
			matlabCodeTextArea.setCaretPosition(0); // Pull focus to the top
			
			// Now write the code to a file
			// 1. Java
			FileWriter codeFileWriter = new FileWriter("../java/infodynamics/demos/autoanalysis/GeneratedCalculator.java");
			codeFileWriter.write(javaCode.toString());
			codeFileWriter.close();
			// 2. Python
			codeFileWriter = new FileWriter("GeneratedCalculator.py");
			codeFileWriter.write(pythonCode.toString());
			codeFileWriter.close();
			// 3. Matlab
			codeFileWriter = new FileWriter("GeneratedCalculator.m");
			codeFileWriter.write(matlabCode.toString());
			codeFileWriter.close();
			
			if (computeResultCheckBox.isSelected()) {
				// Now make our own calculations here:
				// TODO make the calculation in a separate thread, and have a 
				//  button to cancel the calculation whilst it's running
				for (int[] pair : sourceDestPairs) {
					int sourceColumn = pair[0];
					int destColumn = pair[1];
					
					if (sourceColumn == destColumn) {
						continue;
					}
					
					// Initialise
					if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
						calcDiscrete.initialise();
					} else {
						calcContinuous.initialise();
					}
					
					// Set observations
					if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
						calcDiscrete.addObservations(
								MatrixUtils.selectColumn(dataDiscrete, sourceColumn),
								MatrixUtils.selectColumn(dataDiscrete, destColumn));
					} else {
						// TODO We should be able to directly call setObservations(double[], double[])
						//  here but we can't right now because ChannelCalculator does not
						//  define this method. It would be complicated to fix this
						//  because it involves Conditional MI calculator it seems.
						// Revisit this later -- for now fix by deferring to child classes
						// calcContinuous.setObservations(
						//		MatrixUtils.selectColumn(data, sourceColumn),
						//		MatrixUtils.selectColumn(data, destColumn));
						setObservations(calcContinuous,
								MatrixUtils.selectColumn(data, sourceColumn),
								MatrixUtils.selectColumn(data, destColumn));
					}
					
					// Compute the result:
					double result;
					if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
						result = calcDiscrete.computeAverageLocalOfObservations();
					} else {
						result = calcContinuous.computeAverageLocalOfObservations();
					}
					
					String resultsSuffixStringFormatted = "";
					if (statSigCheckBox.isSelected()) {
						// Check statistical significance:
						if (statSigAnalyticCheckBox.isSelected()) {
							// analytic check of statistical significance:
							AnalyticNullDistributionComputer analyticCalc;
							if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
								analyticCalc =
										(AnalyticNullDistributionComputer) calcDiscrete;
							} else {
								analyticCalc =
										(AnalyticNullDistributionComputer) calcContinuous;
							}
							AnalyticMeasurementDistribution measDist =
									analyticCalc.computeSignificance();
							resultsSuffixStringFormatted = String.format(
									resultsSuffixString,
									measDist.pValue);
						} else {
							// permutation test of statistical significance:
							EmpiricalMeasurementDistribution measDist;
							if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
								measDist = calcDiscrete.computeSignificance(numPermutationsToCheck);
							} else {
								measDist = calcContinuous.computeSignificance(numPermutationsToCheck);
							}
							resultsSuffixStringFormatted = String.format(
									resultsSuffixString,
									measDist.getMeanOfDistribution(),
									measDist.getStdOfDistribution(),
									measDist.pValue,
									numPermutationsToCheck);
						}
					}
					
					String resultsText;
					if (sourceDestPairs.size() > 1) {
						// We're doing all pairs
						resultsText = String.format(
								resultsPrefixString + "%.4f %s" + resultsSuffixStringFormatted,
								sourceColumn, destColumn, result, units);
					} else {
						resultsText = String.format(
								resultsPrefixString + "%.4f %s" + resultsSuffixStringFormatted,
								result, units);					
						resultsLabel.setText(resultsText);
					}
					// Write the results to the console:
					System.out.println(resultsText);
				}
	
				if ((sourceDestPairs.size() == 1) && 
						!selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
					// Read the current property values back out (in case of 
					//  automated parameter assignment)
					for (String propName : propertyNames) {
						String propValue = null;
						try {
							propValue = calcContinuous.getProperty(propName);
							propertyValues.put(propName, propValue);
						} catch (Exception ex) {
							ex.printStackTrace(System.err);
							JOptionPane.showMessageDialog(this,
									ex.getMessage());
							resultsLabel.setText("Cannot find a value for property " + propName);
						}
					}
					propertiesTableModel.fireTableDataChanged(); // Alerts to refresh the table contents
				}
			} else {
				resultsLabel.setText("Code generated (no computation performed)");
			}
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Calculation failed, please see console");
		}

	}
	
	/**
	 * Method to assign and initialise our continuous calculator class
	 * @throws Exception 
	 */
	protected abstract ChannelCalculatorCommon assignCalcObjectContinuous(String selectedCalcType) throws Exception;
	
	/**
	 * Method to assign and initialise our discrete calculator class
	 * @throws Exception 
	 */
	protected abstract DiscreteCalcAndArguments assignCalcObjectDiscrete() throws Exception;
	
	/**
	 * Structure used for return calculator and arguments from constructing
	 *  the discrete calculator
	 *  
	 * @author Joseph Lizier
	 *
	 */
	protected class DiscreteCalcAndArguments {
		ChannelCalculatorDiscrete calc;
		int base;
		String arguments;
		
		DiscreteCalcAndArguments(ChannelCalculatorDiscrete calc, int base,
				String arguments) {
			this.calc = calc;
			this.base = base;
			this.arguments = arguments;
		}
	}
	
	/** 
	 * Method to set the observations on the underlying calculator
	 * 
	 * @param calc
	 * @param source
	 * @param dest
	 * @throws Exception
	 */
	protected abstract void setObservations(ChannelCalculatorCommon calc,
			double[] source, double[] dest) throws Exception;


	/**
	 * Extends JTable to add ToolTipText to the property names
	 * 
	 * @author joseph
	 *
	 */
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
		String[] classSpecificPropertyNames = null;
		String[] classSpecificPropertiesFieldNames = null;
		String[] classSpecificPropertyDescriptions = null;
		ChannelCalculatorCommon calc = null;
		try {
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				calcClass = discreteClass;
				calc = null; // Not used
				classSpecificPropertyNames = discreteProperties;
				classSpecificPropertiesFieldNames = null; // Not used
				classSpecificPropertyDescriptions = discretePropertyDescriptions;
			} else {
				calc = assignCalcObjectContinuous(selectedCalcType);
				calcClass = calc.getClass();
				if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
					classSpecificPropertyNames = gaussianProperties;
					classSpecificPropertiesFieldNames = gaussianPropertiesFieldNames;
					classSpecificPropertyDescriptions = gaussianPropertyDescriptions;
				} else if (selectedCalcType.startsWith(CALC_TYPE_KRASKOV)) {
					// The if statement will work for both MI Kraskov calculators
					classSpecificPropertyNames = kraskovProperties;
					classSpecificPropertiesFieldNames = kraskovPropertiesFieldNames;
					classSpecificPropertyDescriptions = kraskovPropertyDescriptions;
				} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
					classSpecificPropertyNames = kernelProperties;
					classSpecificPropertiesFieldNames = kernelPropertiesFieldNames;
					classSpecificPropertyDescriptions = kernelPropertyDescriptions;
				} else {
					calcClass = null;
					throw new Exception("No recognised calculator selected: " +
							selectedCalcType);
				}
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
		
		if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
			// Simply add the properties for this estimator, along
			//  with their default values
			int i = 0;
			propertyValues = new HashMap<String,String>();
			for (String propName : discreteProperties) {
				String propertyDescription = discretePropertyDescriptions[i];
				String defaultPropertyValue = discretePropertyDefaultValues[i];
				i++;
				propertyNames.add(propName);
				propertyDescriptions.add(propertyDescription);
				propertyValues.put(propName, defaultPropertyValue);
				System.out.println("Adding property name " + propName);
			}
		} else {
			// First for the common properties
			int i = 0;
			for (String fieldName : commonContPropertiesFieldNames) {
				String propName = commonContPropertyNames[i];
				String propertyDescription = commonContPropertyDescriptions[i];
				i++;
				System.out.println("Adding property name " + 
						abstractContinuousClass.getSimpleName() + "." + fieldName +
						" = \"" + propName + "\"");
				propertyFieldNames.add(abstractContinuousClass.getSimpleName() + "." + fieldName);
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
					System.out.println("Adding property name " + calcClass.getSimpleName() +
							"." + fieldName + " = \"" + propName + "\"");
					propertyFieldNames.add(calcClass.getSimpleName() + "." + fieldName);
				}
			}
			
			// Now extract the default values for all of these properties:
			propertyValues = new HashMap<String,String>();
			for (String propName : propertyNames) {
				String defaultPropValue = null;
				try {
					defaultPropValue = calc.getProperty(propName);
				} catch (Exception ex) {
					ex.printStackTrace(System.err);
					JOptionPane.showMessageDialog(this,
							ex.getMessage());
					propertyValues.put(propName, "Cannot find a value");
				}
				propertyValues.put(propName, defaultPropValue);
			}
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

	@Override
	public void stateChanged(ChangeEvent e) {
		// I'm not sure how to check which checkbox changed,
		//  so we'll just act on potential changes from both.
		// This won't cause a problem here

		if (e.getSource() == allPairsCheckBox) {
			// "All pairs" checkbox -- update in case changed:
			if (allPairsCheckBox.isSelected()) {
				sourceColTextField.setEnabled(false);
				destColTextField.setEnabled(false);
			} else {
				sourceColTextField.setEnabled(true);
				destColTextField.setEnabled(true);
			}
		} else if (e.getSource() == computeResultCheckBox) {
			// Compute checkbox -- update in case changed:
			if (computeResultCheckBox.isSelected()) {
				computeButton.setText(buttonTextCodeAndCompute);
			} else {
				computeButton.setText(buttonTextCodeOnly);
			}
		} else if (e.getSource() == statSigCheckBox) {
			System.out.println("StatSigCheckBox state changed to " + statSigCheckBox.isSelected());
			// statSigCheckBox -- update statSigAnalyticCheckBox in case changed:
			if (statSigCheckBox.isSelected()) {
				// We might allow the analytic box to be checked if the selected
				//  calculator allows this:
				if (calcImplementsAnalyticStatSig()) {
					statSigAnalyticCheckBox.setEnabled(true);
				} else {
					statSigAnalyticCheckBox.setSelected(false);
					statSigAnalyticCheckBox.setEnabled(false);
				}
			} else {
				statSigAnalyticCheckBox.setSelected(false);
				statSigAnalyticCheckBox.setEnabled(false);
			}
		}
	}
	
	/**
	 * Checks whether the current calculator class implements the 
	 * AnalyticNullDistributionComputer interface
	 * 
	 * @return
	 */
	protected boolean calcImplementsAnalyticStatSig() {
		Class<?> currentClass = calcClass;
		while (currentClass != null) {
			// This getInterfaces call doesn't seem to return the interfaces
			//  implemented by the superclass, so we need to go up the tree
			//  to check all of those. That seems odd, surely there is a simpler way?
			Class<?>[] interfacesOfThisClass = currentClass.getInterfaces();
			for (Class<?> theInterface : interfacesOfThisClass) {
				if (theInterface == AnalyticNullDistributionComputer.class) {
					return true;
				}
			}
			currentClass = currentClass.getSuperclass();
		}
		return false;
	}
}
