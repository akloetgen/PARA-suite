package mapping;

import java.util.List;

public class MappingErrorException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = -720715720548525741L;

	private List<String> mappingCommand;

	public String getMappingCommand() {
		String mappingCommandString = "";
		for (int i = 0; i < mappingCommand.size(); i++) {
			mappingCommandString += mappingCommand.get(i) + " ";
		}
		return mappingCommandString;
	}

	public void setMappingCommand(List<String> mappingCommand) {
		this.mappingCommand = mappingCommand;
	}

}
