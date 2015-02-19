package main;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Properties;

import enums.PARMAPropertiesEnum;

public final class PARMAProperties {

	public static Properties properties;
	public static String propertiesFileName = "parma.properties";

	public static void setProperty(PARMAPropertiesEnum propertyName,
			String newPropertyValue) {
		try {
			File newPropertiesFile = new File(propertiesFileName);
			if (!newPropertiesFile.exists()) {
				newPropertiesFile.createNewFile();
			}
			properties = new Properties();
			BufferedInputStream stream;
			stream = new BufferedInputStream(new FileInputStream(
					newPropertiesFile));
			properties.load(stream);

			OutputStream out = new FileOutputStream(newPropertiesFile);

			properties.setProperty(propertyName.toString(), newPropertyValue);
			properties.store(out, propertyName.toString());

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static String getProperty(PARMAPropertiesEnum propertyName) {
		try {
			File newPropertiesFile = new File(propertiesFileName);
			if (!newPropertiesFile.exists()) {
				newPropertiesFile.createNewFile();
				// FileWriter propertiesFileWriter = new FileWriter(
				// newPropertiesFile);
				// propertiesFileWriter
				// .write("# if the bwa parma mod is not in the class_path,"
				// + " please classify the path to parma here"
				// + System.getProperty("line.separator"));
				// propertiesFileWriter.write(PARMAPropertiesEnum.PARMA_LOCATION
				// .toString()
				// + "="
				// + System.getProperty("line.separator"));
				// propertiesFileWriter.write("# if bwa is not in the class_path,"
				// + " please classify the path to bwa here"
				// + System.getProperty("line.separator"));
				// propertiesFileWriter.write(PARMAPropertiesEnum.BWA_LOCATION
				// .toString()
				// + "="
				// + System.getProperty("line.separator"));
				// propertiesFileWriter
				// .write("# if bowtie2 is not in the class_path,"
				// + " please classify the path to bowtie2 here"
				// + System.getProperty("line.separator"));
				// propertiesFileWriter.write(PARMAPropertiesEnum.BT2_LOCATION
				// .toString()
				// + "="
				// + System.getProperty("line.separator"));
				// propertiesFileWriter.close();
				// // CREATE WARNING=??? OR CREATE STDIN SO THAT THE USER CAN
				// SPECIFY?!?!?!
			}
			properties = new Properties();
			BufferedInputStream stream;
			stream = new BufferedInputStream(new FileInputStream(
					newPropertiesFile));
			properties.load(stream);
			// properties.save(null, propertiesFileName);
			// properties.setProperty(propertiesFileName, propertiesFileName);
			// stream.close();
			return properties.getProperty(propertyName.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
		return "";
	}
}
