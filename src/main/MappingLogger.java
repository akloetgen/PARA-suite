package main;

import java.io.IOException;
import java.util.Map.Entry;
import java.util.Properties;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

public final class MappingLogger {

	private static Logger logger = Logger.getLogger("mappinglogger");
	private static String loggingFile = "mapping_errors.log";

	/**
	 * Returns logger for given properties-file.
	 * 
	 * @return logger with loaded properties.
	 */
	public static Logger getLogger() {
		try {
			Properties prop = new Properties();
			prop.load(MappingLogger.class.getClassLoader().getResourceAsStream(
					"mappinglogger.properties"));

			for (Entry<Object, Object> propEntry : prop.entrySet()) {
				if (propEntry.getKey().equals("log4j.appender.A6.File")) {
					propEntry.setValue(loggingFile);
				}
			}
			PropertyConfigurator.configure(prop);
		} catch (IOException e) {
			e.printStackTrace();
		}

		return logger;
	}

	public static void setLoggingFile(String loggingFile) {
		MappingLogger.loggingFile = loggingFile;
	}
}