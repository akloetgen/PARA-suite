package mapping;

public class ExternalCallErrorException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = -720715720548525741L;

	private String mappingCommand;

	public ExternalCallErrorException(String mappingCommand) {
		this.mappingCommand = mappingCommand;
	}

	public String getMappingCommand() {
		return mappingCommand;
	}

	public void setMappingCommand(String mappingCommand) {
		this.mappingCommand = mappingCommand;
	}

}
