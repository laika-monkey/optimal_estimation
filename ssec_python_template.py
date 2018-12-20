#!/usr/bin/env python
# encoding: utf-8
"""
FIXME

Created by FIXME on 2009-FIXME-FIXME.
Copyright (c) 2009 University of Wisconsin SSEC. All rights reserved.
"""

import os, sys, logging

LOG = logging.getLogger(__name__)


def main():
    import optparse
    usage = """
%prog [options] 
run "%prog help" to list commands
"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--test', dest="self_test",
                    action="store_true", default=False, help="run unit tests")            
    parser.add_option('-q', '--quiet', dest="quiet",
                    action="store_true", default=False, help="only error output")
    parser.add_option('-v', '--verbose', dest="verbose",
                    action="store_true", default=False, help="enable more informational output")   
    parser.add_option('-w', '--debug', dest="debug",
                    action="store_true", default=False, help="enable debug output")   
    options, args = parser.parse_args()
    if options.self_test:
        import doctest
        doctest.testmod()
        sys.exit(2)

    lvl = logging.WARNING
    if options.debug: lvl = logging.DEBUG
    elif options.verbose: lvl = logging.INFO
    elif options.quiet: lvl = logging.ERROR
    logging.basicConfig(level = lvl)

    commands = {}
    prior = None
    prior = dict(locals())
    
    def create(*args):
        """create summary
        create extended info
        """
        LOG.info("creating database")
        
    def build(*args):
        """build summary
        build extended info
        """
        LOG.info("building database tables")
        
    def grant(*args):
        """grant summary
        grant extended info
        """
        LOG.info("granting permissions for tables")
        
    def index(*args):
        """index summary
        index extended info
        """
        LOG.info("creating indices for tables")
        
    def help(command=None):
        """print command help or list of commands
        e.g. help help
        """
        if command is None: 
            # print first line of docstring
            for cmd in commands:
                ds = commands[cmd].__doc__.split('\n')[0]
                print "%-16s %s" % (cmd,ds)
        else:
            print commands[command].__doc__
            
    def test():
        "run tests"
        test1()
        
    commands.update(dict(x for x in locals().items() if x[0] not in prior))    
    
    if (not args) or (args[0] not in commands): 
        parser.print_help()
        help()
        return 9
    else:
        locals()[args[0]](*args[1:])

    return 0


if __name__=='__main__':
    sys.exit(main())