''' removing duplicated files from data folders '''

aatrj="pro.xtc"
gmxtools="$HOME/gmx/src/tools"

curdir = os.getcwd()
prjdir = os.path.join(curdir, "..")
native = os.path.join(prjdir, "naked.tpr")

def dif(id, f1 = None, f2 = None):
  p1 = os.path.join(curdir, "data%d" % id)
  p2 = os.path.join(curdir, "data%d" % (id+1))

  if f1: p1 = os.path.join(p1, f1)
  if f2: p2 = os.path.join(p2, f2)
  if not os.path.exists(p1):
    print "path %s doesn't exist id = %d" % (p1, id)
    return -1
  if not os.path.exists(p2):
    print "path %s doesn't exist, id = %d" % (p2, id)
    return -1

  if not f1: return 0 # just checking directories

  ret = os.system("cmp %s %s" % (p1, p2))
  if ret != 0:
    print "corruption: %s differs from %s" % (p1, p2)
    exit(1)
  else:
    print "removing %s" % p1
    os.remove(p1)

print curdir
id = 1
while 1:
  if dif(id) < 0:
    print "quit at %d" % id
    break
  dif(id, "mb.av", "mb.av0")
  dif(id, "MTSEED", "MTSEED0")
  dif(id, "hist.bin", "hist0.bin")
  dif(id, "state.cpt", "state0.cpt")
  id += 1

