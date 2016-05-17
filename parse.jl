using ArgParse
using YAML

function parse_config(filename,verbose)
    config = YAML.load(open(filename))
    fields_config = config["fields-config"]
    trajsolver_config = config["trajsolver-config"]
    if verbose == true
        display(fields_config);println("")
        display(trajsolver_config);println("")
    end
    ##TODO: check config format
    return fields_config, trajsolver_config
end
