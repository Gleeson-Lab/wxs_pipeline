from itertools import chain

def get_options(rule_name, config, prefix="--", separator=" ", extra=[],
                ignore=set()):
    """
    Return extra parameters for the specified rule defined in the config file.
    :param rule_name: used to look in the appropriate section of the config file
    :param config: the config dictionary
    :param prefix: the string to prepend to parameters
    :param separator: the string to separate the parameter from its value
    :param extra: an iterable of (parameter, value) entries to additionally include
    :param ignore: a set of parameters to ignore
    """
    parameters = []
    for parameter, values in chain(config[rule_name].items(), extra):
        if parameter in ignore:
            continue
        if type(values) is not list:
            values = [values]
        for value in values:
            parameters.append(f"{prefix}{parameter}{separator}{value}")
    return " ".join(parameters)
