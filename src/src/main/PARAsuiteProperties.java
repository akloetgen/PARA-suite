package main;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Properties;

import enums.PARAsuitePropertiesEnum;

/**
 * Saves and loades properties for the PARA-suite toolkit from parasuite.properties, such
 * as the paths to the respective read aligners.
 * 
 * @author akloetgen
 * 
 */
public final class PARAsuiteProperties {

	public static Properties properties;
	public static String propertiesFileName = "parasuite.properties";

	public static void setProperty(PARAsuitePropertiesEnum propertyName,
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
			e.printStackTrace();
		}
	}

	public static String getProperty(PARAsuitePropertiesEnum propertyName) {
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
			return properties.getProperty(propertyName.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
		return "";
	}
}
