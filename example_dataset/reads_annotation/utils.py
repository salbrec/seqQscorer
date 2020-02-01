import subprocess

def getSystemCall(call):
	process = subprocess.Popen(call, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()
	out = out.decode(locale.getdefaultlocale()[1])
	err = err.decode(locale.getdefaultlocale()[1])
	if process.returncode:
		print("call failed, call was: %s" % ' '.join(call), file=sys.stderr)
		print("Message was: %s" % str(out), file=sys.stderr)
		print("Error code was %s, stderr: %s" % (process.returncode, err), file=sys.stderr, end='')
		raise Exception('runSystemCall Exception') 
	return out, err


