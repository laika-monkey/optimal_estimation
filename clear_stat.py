#!/usr/bin/env python
#
# purpose	clear error state for jobs in queue because I can't figure out qmod -cq


def main():
	import os, re
	cmd = 'qstat'
	p = os.popen(cmd, 'r')
	stdio = p.readlines()

	re_qstat = re.compile('\s*(\d{6})\s+[\d.]+\s+.*?\s+wsessions\s+(\w+)\s+') #_1 proc, 2 state
	
	for line in stdio:
		reg = re_qstat.search(line)
		if not reg:
			continue

		#_pulling out state
		proc_id, state = reg.group(1), reg.group(2)

		if state == 'Eqw':
			print 'clearing process {0} has state {1}'.format(proc_id, state)
			os.popen('qmod -cj {0}'.format(proc_id))
		if state == 'r':
			print 'process {0} currently running'.format(proc_id)


if __name__ == '__main__':
	main()


'''
 541540 0.55500 run_oe.py  wsessions    r     07/09/2015 14:44:20 all.q@calzone.ssec.wisc.edu        1        
'''
