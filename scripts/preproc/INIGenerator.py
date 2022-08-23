#! /usr/bin/env python
'''#############################################################################
Converts a XML configuration file for CGS models into INI configuration file(s).

Brent Smith (February 2022), JHU APL
#############################################################################'''

# imports (std lib)
import os
import sys
import shutil
import pprint
import argparse
import logging
import configparser
import xml.etree.ElementTree as ET
# imports (3rd party)
# imports (local)

# module-based configurations
logging.basicConfig()
logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)
#===============================================================================
def parse_args(args=None):
    '''=========================================================================
    Command-Line Argument Parser
    ========================================================================='''
    parser = argparse.ArgumentParser()

    required = parser.add_argument_group('required positional arguments')
    required.add_argument('input', help='Input XML file.')
    required.add_argument('output', help='Output INI file.')
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Activate verbose execution mode. (default: %(default)s)'
    )

    try:
        if len(args) == 0:
            parser.print_help()
            sys.exit(2)
        args = parser.parse_args(args)
    except Exception as e:
        logging.getLogger(__name__).error(e)
        parser.print_help()
        sys.exit(2)

    return args

def get_xml_children(tree_element):
    '''=========================================================================
    Obtains recursively all children tags and attributes from the parent tree
    element
    ========================================================================='''
    if tree_element:
        children = []
        if len(list(tree_element)) > 1:
            children = {tree_element.tag:{}}
            for child in tree_element:
                child_dict = get_xml_children(child)
                for tag, attrib in child_dict.items():
                    children[tree_element.tag][tag] = attrib
            return children
        else:
            for child in tree_element:
                return {tree_element.tag: get_xml_children(child)}
    else:
        return {tree_element.tag: tree_element.attrib}

def set_ini_config(config, element, section=None):
    '''=========================================================================
    Uses configparser to assign section headings, section heading comments, and
    configuration parameters with values from an XML element tree.
    ========================================================================='''
    if element:
        # at least one child
        num_children = len(element)
        if num_children > 1:
            # more than 1 child
            for key, value in element.items():
                config.add_section(key)
                for _key, _value in value.items():
                    config.set(key, _key, _value)
        else:
            # one child {key:{key:value}} or {key:value}
            try:
                # values can be another dictionary or string
                for key, value in element.items():
                    set_ini_config(config, value, section=key)
            except:
                # it's a single dictionary {key:value}
                config[section] = element
    else:
        # no children of element
        print('...no children...')
        print(element)
    return config

def main(arguments=None):
    '''=========================================================================
    Main driver that converts an XML data file to an INI data file.
    ========================================================================='''
	# parse command-line arguments
    args = parse_args(arguments)

    if args.verbose:
        logging.getLogger(__name__).setLevel(logging.INFO)

	# we do not have a template INI since they are basically YAML (single-level)

	# Read XML file
    with open(args.input) as f:
        xml_data_str = f.read() # this contains all content
    root = ET.fromstring(xml_data_str) # this ignores top-level comments
    xml_data = get_xml_children(root)

    logging.getLogger(__name__).info('Retrieved XML data...')
    if args.verbose:
        pprint.pprint(xml_data)

    # remove kaiju level - case insensitive
    xml_data_keys = [key.lower() for key in xml_data.keys()]
    if 'kaiju' in xml_data_keys:
        index = xml_data_keys.index('kaiju')
        xml_data = xml_data[list(xml_data.keys())[index]]

    # now the keys of the dictionary are the models
    models = xml_data.keys()

    # Create INI configuration for each model
    for model in models:
        # write INI file (can't append to an existing INI file, so we replace anew)
        config = configparser.ConfigParser(allow_no_value=True)
        config.optionxform = str # preserve case of strings

        # create INI configuration
        logging.getLogger(__name__).info('Creating {}-based configuration...'.format(str(model)))
        config = set_ini_config(config, xml_data[model])

        # output to temp file named after model
        logging.getLogger(__name__).info('Writing temp file...')
        with open(model, 'w') as f:
            config.write(f)

        # append heading to configuration
        with open(model, 'r') as f:
            contents = f.read()

        contents = '## ' + model.upper() + ' ##\n' + contents

        logging.getLogger(__name__).info('Writing INI data...')
        with open(model, 'w') as f:
            f.write(contents)

    # concatenate all model-named temp files
    logging.getLogger(__name__).info('Removing all temp files...')
    with open(args.output, 'w') as f:
        for model in models:
            with open(model, 'r') as fi:
                shutil.copyfileobj(fi, f)
            os.remove(model)

    logging.getLogger(__name__).info('Complete')
    return 0


if __name__ == '__main__':
    # put main stuff in the MAIN function
    sys.exit(main(sys.argv[1:]))
