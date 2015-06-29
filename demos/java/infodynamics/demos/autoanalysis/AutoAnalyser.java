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
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JComboBox;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JFileChooser;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.EtchedBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.MediaTracker;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Vector;

public class AutoAnalyser extends JFrame implements ActionListener, DocumentListener {

	/**
	 * Need serialVersionUID to be serializable
	 */
	private static final long serialVersionUID = 1L;
	
	// Set options for the TE Calculator type:
	public final static String CALC_TYPE_DISCRETE = "Discrete";
	public final static String CALC_TYPE_GAUSSIAN = "Gaussian";
	public final static String CALC_TYPE_KRASKOV  = "Kraskov";
	public final static String CALC_TYPE_KERNEL   = "Kernel";
	protected String[] calcTypes = {
			CALC_TYPE_DISCRETE, CALC_TYPE_GAUSSIAN,
			CALC_TYPE_KRASKOV, CALC_TYPE_KERNEL};
	protected String[] unitsForEachCalc = {"bits", "nats", "nats", "bits"};
	
	// Store the properties for each calculator
	protected String[] commonContProperties = {
			TransferEntropyCalculator.K_PROP_NAME,
			
	};
	protected String[] gaussianProperties = {
			TransferEntropyCalculator.K_TAU_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.L_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.L_TAU_PROP_NAME, // Not common to Kernel
			TransferEntropyCalculator.DELAY_PROP_NAME, // Not common to Kernel
	};
	protected String[] kernelProperties = {
			TransferEntropyCalculatorKernel.KERNEL_WIDTH_PROP_NAME,
			TransferEntropyCalculatorKernel.DYN_CORR_EXCL_TIME_NAME,
			TransferEntropyCalculatorKernel.NORMALISE_PROP_NAME,			
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
	protected Vector<String> fieldNames;
	// Values of the properties
	protected HashMap<String,String> propertyValues;
	// Results text
	protected JLabel resultsLabel;
	// Text area for the generated code
	protected JTextArea codeTextArea;
	
	protected String codeDefaultText = "    ... Awaiting new parameter selection (press compute) ...";
	
	/**
	 * Constructor to generate the application windows
	 */
	public AutoAnalyser() {
		
		ImageIcon icon = new ImageIcon("../../JIDT-logo.png"); // Location for distributions
		if (icon.getImageLoadStatus() != MediaTracker.COMPLETE) {
			// Try the alternative image location for SVN checkouts
			icon = new ImageIcon("../../web/JIDT-logo.png");
		}
		setIconImage(icon.getImage());
			
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
		// Don't set border around this text field as it doesn't look right
		// TODO set action listener for on click of this field
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
		
		// TODO Need to set an appropriate width here
		putCalcPropertiesInTable();
		propertiesTableModel = new PropertiesTableModel();
		propertiesTable = new JTable(propertiesTableModel);
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
		codeTextArea = new JTextArea(codeDefaultText);
		codeTextArea.setOpaque(false);
		codeTextArea.setEditable(false);
		//codeTextArea.setLineWrap(true); // This makes it all unreadable
		//codeTextArea.setWrapStyleWord(true);
		JScrollPane areaScrollPane = new JScrollPane(codeTextArea);
        areaScrollPane.setVerticalScrollBarPolicy(
        		JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        areaScrollPane.setHorizontalScrollBarPolicy(
        		JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        int codeTextAreaWidth = 560;
        int codeTextAreaHeight = 460;
        Dimension codeTextAreaDimension = 
        		new Dimension(codeTextAreaWidth, codeTextAreaHeight);
        areaScrollPane.setPreferredSize(codeTextAreaDimension);
        areaScrollPane.setMinimumSize(codeTextAreaDimension);
        areaScrollPane.setMaximumSize(codeTextAreaDimension);
		codePanel.add(areaScrollPane);
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
        
		// TODO Add the JIDT logo somewhere
		
		// Add both panels into the frame with Border layout
		add(paramsPanel, BorderLayout.WEST);
		add(codePanel);
		
		setVisible(true);		
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		// Clear text fields for now if anything changes
		resultsLabel.setText(" ");
		codeTextArea.setText(codeDefaultText);
		
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
		
		StringBuffer javaCode = new StringBuffer();
		javaCode.append("package infodynamics.demos.autoanalysis;\n\n");
		javaCode.append("import infodynamics.utils.ArrayFileReader;\n");
		javaCode.append("import infodynamics.utils.MatrixUtils;\n\n");
		javaCode.append("import infodynamics.measures.continuous.TransferEntropyCalculator;\n");

		try{
			// Create a Kraskov TE calculator:
			TransferEntropyCalculator teCalc;
			
			// Construct an instance of the selected calculator:
			String javaConstructorLine = null;
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				throw new Exception("Not implemented yet");
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
				teCalc = new TransferEntropyCalculatorGaussian();
				javaCode.append("import infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian;\n");
				javaConstructorLine = "    teCalc = new TransferEntropyCalculatorGaussian();\n";
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
				teCalc = new TransferEntropyCalculatorKraskov();
				javaCode.append("import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov;\n");
				javaConstructorLine = "    teCalc = new TransferEntropyCalculatorKraskov();\n";
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
				teCalc = new TransferEntropyCalculatorKernel();
				javaCode.append("import infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel;\n");
				javaConstructorLine = "    teCalc = new TransferEntropyCalculatorKernel();\n";
			} else {
				throw new Exception("No recognised calculator selected: " +
						selectedCalcType);
			}
			
			javaCode.append("\npublic class GeneratedTECalculator {\n\n");
			javaCode.append("  public static void main(String[] args) throws Exception {\n");
			
			javaCode.append("    String dataFile = \"" + dataFile.getAbsolutePath() + "\";\n");
			javaCode.append("    ArrayFileReader afr = new ArrayFileReader(dataFile);\n");
			javaCode.append("    double[][] data = afr.getDouble2DMatrix();\n");
			
			javaCode.append("    TransferEntropyCalculator teCalc;\n");
			javaCode.append(javaConstructorLine);
			
			// Set properties 
			for (String propName : fieldNames) {
				String propValue = null;
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
					javaCode.append("    teCalc.setProperty(\"" + propName + "\", \"" +
								propValue + "\");\n");
				}
			}
			
			teCalc.initialise();
			javaCode.append("    teCalc.initialise();\n");
			// Compute TE 
			teCalc.setObservations(
					MatrixUtils.selectColumn(data, sourceColumn),
					MatrixUtils.selectColumn(data, destColumn));
			javaCode.append("    teCalc.setObservations(\n" +
						"        MatrixUtils.selectColumn(data, " + sourceColumn + "),\n" +
						"        MatrixUtils.selectColumn(data, " + destColumn + "));\n");
			double teResult = teCalc.computeAverageLocalOfObservations();
			javaCode.append("    double teResult = teCalc.computeAverageLocalOfObservations();\n");
			String resultsPrefixString = String.format("TE_%s(col_%d -> col_%d) = ",
					selectedCalcType, sourceColumn, destColumn);
			resultsLabel.setText(String.format(resultsPrefixString + "%.4f %s", teResult, units));
			javaCode.append("    System.out.printf(\"" + resultsPrefixString + "%.4f " + units + "\\n\", teResult);\n");
			javaCode.append("  }\n");
			javaCode.append("}\n\n");
			// System.out.println("Code to generate this result:");
			// System.out.print(code);
			codeTextArea.setText(javaCode.toString());
			
			// Now write the Java code to a file
			FileWriter javaFileWriter = new FileWriter("infodynamics/demos/autoanalysis/GeneratedTECalculator.java");
			javaFileWriter.write(javaCode.toString());
			javaFileWriter.close();
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
			JOptionPane.showMessageDialog(this,
					ex.getMessage());
			resultsLabel.setText("Calculation failed, please see console");
		}

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
			String propName = fieldNames.get(rowIndex);
			propertyValues.put(propName, (String) aValue);
		}

		@Override
		public int getColumnCount() {
			return 2;
		}

		@Override
		public int getRowCount() {
			return fieldNames.size();
		}

		@Override
		public Object getValueAt(int rowIndex, int columnIndex) {
			String propName = fieldNames.get(rowIndex);
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
		String[] classSpecificProperties = null;
		TransferEntropyCalculator teCalc = null;
		try {
			if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_DISCRETE)) {
				throw new Exception("Not implemented yet");
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_GAUSSIAN)) {
				teCalcClass = TransferEntropyCalculatorGaussian.class;
				teCalc = new TransferEntropyCalculatorGaussian();
				classSpecificProperties = gaussianProperties;
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KRASKOV)) {
				teCalcClass = TransferEntropyCalculatorKraskov.class;
				teCalc = new TransferEntropyCalculatorKraskov();
				classSpecificProperties = kraskovProperties;
			} else if (selectedCalcType.equalsIgnoreCase(CALC_TYPE_KERNEL)) {
				teCalcClass = TransferEntropyCalculatorKernel.class;
				teCalc = new TransferEntropyCalculatorKernel();
				classSpecificProperties = kernelProperties;
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
		fieldNames = new Vector<String>();
		Vector<String> fieldDescriptions = new Vector<String>();
		for (int i = 0; i < commonContProperties.length; i++) {
			fieldNames.add(commonContProperties[i]);
			System.out.println("Adding property name " + commonContProperties[i]);
		}
		for (int i = 0; i < classSpecificProperties.length; i++) {
			fieldNames.add(classSpecificProperties[i]);
			System.out.println("Adding property name " + classSpecificProperties[i]);
		}
		// Now extract the default values for all of these properties:
		propertyValues = new HashMap<String,String>();
		for (String propName : fieldNames) {
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
		codeTextArea.setText(codeDefaultText);
	}

	@Override
	public void insertUpdate(DocumentEvent e) {
		// Source or dest col number updated
		resultsLabel.setText(" "); // Clear text for now if anything changes
		codeTextArea.setText(codeDefaultText);
	}

	@Override
	public void removeUpdate(DocumentEvent e) {
		// Source or dest col number updated
		resultsLabel.setText(" "); // Clear text for now if anything changes
		codeTextArea.setText(codeDefaultText);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new AutoAnalyser();
	}
}
